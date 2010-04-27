
#include "sift.h"
#include "util.h"
#include "const.h"
#include "assert.h"
#include <stdlib.h>

SiftToolbox::SiftToolbox()
{
	m_baseImage = NULL;
}

SiftToolbox::~SiftToolbox()
{

}

bool SiftToolbox::InitImage(IplImage* image)
{
	IplImage* gray8 = NULL, *gray32 = NULL;

	gray8 = cvCreateImage(cvGetSize(image), IPL_DEPTH_8U, 1);
	gray32 = cvCreateImage(cvGetSize(image), IPL_DEPTH_32F, 1);
	if(image->nChannels == 1){
		gray8 = cvCloneImage(image);
	}else{
		cvCvtColor(image, gray8, CV_RGB2GRAY);
	}
	cvConvertScale(gray8, gray32, 1.0/255, 0);
	cvReleaseImage(&gray8);
	if(m_param.isPriorDouble){
		IplImage* doubleImage = cvCreateImage(cvSize(image->width*2, 
			image->height*2), IPL_DEPTH_32F, 1);
		cvResize(gray32, doubleImage, CV_INTER_LINEAR);
		cvReleaseImage(&gray32);
		m_baseImage = doubleImage;
	}else{
		m_baseImage = gray32;
	}
	cvSmooth(m_baseImage, m_baseImage, CV_GAUSSIAN, 0, 0, 
		m_param.sigmaOfInitGaussPyr, m_param.sigmaOfInitGaussPyr);
	return true;
}

IplImage* SiftToolbox::Process()
{
	return NULL;
}

bool SiftToolbox::BuildDogPyr()
{
	IplImage*** gaussPyr;
	int gaussScaleNum = m_param.numOfScalePerOctave + 3;
	int dogScaleNum = gaussScaleNum - 1;
	double* pSigmaArry = NULL;

	gaussPyr = (IplImage***)calloc(m_param.numOfOctave, sizeof(IplImage**));
	assert(gaussPyr != NULL);
	m_dogPyr = (IplImage***)calloc(m_param.numOfOctave, sizeof(IplImage**));
	assert(m_dogPyr != NULL);
	for(int i=0; i<m_param.numOfOctave; i++){
		gaussPyr[i] = (IplImage**)calloc(gaussScaleNum, sizeof(IplImage*));
		assert(gaussPyr[i] != NULL);
		m_dogPyr[i] = (IplImage**)calloc(dogScaleNum, sizeof(IplImage*));
		assert(m_dogPyr[i] != NULL);
	}

	pSigmaArry = (double*)calloc(gaussScaleNum, sizeof(double));
	pSigmaArry[0] = m_param.sigmaOfInitGaussPyr;
	double k = pow(2.0, 1.0/m_param.numOfScalePerOctave);
	double pre = pSigmaArry[0];
	for(int i=1; i<gaussScaleNum; i++){
		pSigmaArry[i] = pre*k;
		pre = pSigmaArry[i];
	}
	
	gaussPyr[0][0] = cvCloneImage(m_baseImage);
	for(int octave = 1; octave < m_param.numOfOctave; octave++){
		cvPyrDown(gaussPyr[octave-1][0], gaussPyr[octave][0], CV_GAUSSIAN_5x5);
	}
	
	for(int octave = 0; octave < m_param.numOfOctave; octave++){
		CvSize currOctaveSize = cvGetSize(gaussPyr[octave][0]);
		for(int scale = 1; scale < gaussScaleNum; scale++){
			gaussPyr[octave][scale] = cvCreateImage(currOctaveSize, IPL_DEPTH_32F, 1);
			cvSmooth(gaussPyr[octave][scale], gaussPyr[octave][scale], 
				CV_GAUSSIAN, 0, 0, pSigmaArry[scale], pSigmaArry[scale]);

			m_dogPyr[octave][scale-1] = cvCreateImage(currOctaveSize, IPL_DEPTH_32F, 1);
			cvSub(gaussPyr[octave][scale-1], gaussPyr[octave][scale], NULL);
		}
	}
	
	free(pSigmaArry);
	for(int i=0; i<m_param.numOfOctave; i++){
		free(&gaussPyr[i]);
	}
	free(&gaussPyr);
	return true;
}

void SiftToolbox::FindExtremePoint()
{
	for(int i=0; i<m_param.numOfOctave; i++){
		for(int j=1; j<m_param.numOfScalePerOctave-1; j++){
			CvSize scaleSize = cvGetSize(m_dogPyr[i][0]);
			for(int x = 1; x < scaleSize.width-1; x++){
				for(int y = 1; y< scaleSize.height-1; y++){
					if(IsExtrema(i,j,x,y)){
						SiftKeypoint_t* pKeyPoint = new SiftKeypoint_t;
						pKeyPoint->octave = i;
						pKeyPoint->scale = j;
						pKeyPoint->pos.x = x;
						pKeyPoint->pos.y = y;
						m_keyPointPool.push_back(*pKeyPoint);
					}
				}
			}
			
		}
	}
}
bool SiftToolbox::IsExtrema(int octave, int scale, int x, int y)
{
	float val = GetVal32f(m_dogPyr[octave][scale], x, y);
	float val1 = GetVal32f(m_dogPyr[octave][scale], x-1, y);
	if(val > val1){ //Check Max Value
		for(int i = -1; i<=1; i++){
			for(int j = -1; j<=1; j++){
				for(int k = -1; k<= 1; k++){
					if(val <= GetVal32f(m_dogPyr[octave][scale+i], x+j, y+k)){
						return false;
					}
				}
			}
		}
	}else if(val < val1){ //Check Min Value
		for(int i = -1; i<=1; i++){
			for(int j = -1; j<=1; j++){
				for(int k = -1; k<= 1; k++){
					if(val >= GetVal32f(m_dogPyr[octave][scale+i], x+j, y+k)){
						return false;
					}
				}
			}
		}
	}else{
		return false;
	}
	return true;
}