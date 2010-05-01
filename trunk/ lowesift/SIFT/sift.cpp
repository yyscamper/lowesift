
#include "sift.h"
#include "util.h"
#include "const.h"
#include "assert.h"
#include <stdlib.h>

SiftToolbox::SiftToolbox()
{
	m_param.isPriorDouble = true;
	m_param.nearMaxRatio = 0.8;
	m_param.numOfHistBins = 36;
	m_param.numOfScalePerOctave = 3;
	m_param.sigmaOfOriAssign = 1.5;
	m_param.sigmaOfInitGaussPyr = 1.6;
	m_param.sigmaOfPriorDouble = 1.6;
	m_param.ratioOfEdge = 10;
	m_param.radiusOfOriHistWnd = cvRound(m_param.sigmaOfOriAssign * 3);
	m_param.thresholdOfKeypointContrast = 0.03;
	
	m_dogPyr = NULL;
	m_pSigmaVal = NULL;
	m_pyrBaseImage = NULL;
	m_maxMag = 0.0;
}

SiftToolbox::~SiftToolbox()
{

}

IplImage* SiftToolbox::Process(IplImage* image)
{
	m_originImage = image;
	InitImage();
	BuildDogPyr();
	FindExtremePoint();
	OrientationAssignment();
	return NULL;
}

void SiftToolbox::InitImage()
{
	IplImage* gray8 = NULL, *gray32 = NULL;

	gray8 = cvCreateImage(cvGetSize(m_originImage), IPL_DEPTH_8U, 1);
	gray32 = cvCreateImage(cvGetSize(m_originImage), IPL_DEPTH_32F, 1);
	if(m_originImage->nChannels == 1){
		gray8 = cvCloneImage(m_originImage);
	}else{
		cvCvtColor(m_originImage, gray8, CV_RGB2GRAY);
	}
	cvConvertScale(gray8, gray32, 1.0/255, 0);
	cvReleaseImage(&gray8);
	if(m_param.isPriorDouble){
		IplImage* doubleImage = cvCreateImage(cvSize(m_originImage->width*2, 
			m_originImage->height*2), IPL_DEPTH_32F, 1);
		cvResize(gray32, doubleImage, CV_INTER_LINEAR);
		cvReleaseImage(&gray32);
		m_pyrBaseImage = doubleImage;
	}else{
		m_pyrBaseImage = gray32;
	}
	cvSmooth(m_pyrBaseImage, m_pyrBaseImage, CV_GAUSSIAN, 0, 0, 
		m_param.sigmaOfInitGaussPyr, m_param.sigmaOfInitGaussPyr);
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
	
	m_pSigmaVal = (double*)calloc(m_param.numOfScalePerOctave, sizeof(double));
	memcpy(m_pSigmaVal, pSigmaArry+1, m_param.numOfScalePerOctave);

	gaussPyr[0][0] = cvCloneImage(m_pyrBaseImage);
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
					if(IsExtrema(i,j,x,y) 
						|| IsRemovableForLowContrast(m_dogPyr[i][j], x, y)
						|| IsRemovableForEdge(m_dogPyr[i][j], x, y)){
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

bool SiftToolbox::IsExtrema(const int octave, const int scale, const int x, const int y) const 
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


bool SiftToolbox::IsRemovableForLowContrast(IplImage* img, int row, int col) const
{
	return false;
}

bool SiftToolbox::IsRemovableForEdge(IplImage* img, int row, int col) const
{
	double Dxx = GetVal32f(img, row, col+1) + GetVal32f(img, row, col-1) - 2*GetVal32f(img, row, col);
	double Dyy = GetVal32f(img, row-1, col) + GetVal32f(img, row+1, col) - 2*GetVal32f(img, row, col);
	double Dxy = (GetVal32f(img, row+1, col+1) + GetVal32f(img, row-1, col-1)
					- GetVal32f(img, row+1, col-1) - GetVal32f(img, row-1, col+1))/4.0;
	double trace = Dxx + Dyy;
	double det = Dxx*Dyy - Dxy*Dxy;
	if(det <= 0 || trace*trace/det >= pow(m_param.ratioOfEdge+1,2)/m_param.ratioOfEdge){
		return true;
	}
	return false;	
}


void SiftToolbox::OrientationAssignment()
{
	SiftKeypoint_t* pCurrKey = NULL;
	for(list<SiftKeypoint_t>::iterator pCurrKey = m_keyPointPool.begin(); pCurrKey != m_keyPointPool.end(); pCurrKey++){
		CalcDormOri(*pCurrKey);
	}
}

void SiftToolbox::CalcDormOri(SiftKeypoint_t& key)
{
	double* hist = CreateOriHist(m_dogPyr[key.octave][key.scale], key.pos.x, key.pos.y, m_pSigmaVal[key.scale]);	
	double max = hist[0];
	int maxCnt = 0;
	int maxIndex = 0;
	//Find the max one
	for(int i=1; i<m_param.numOfHistBins; i++){
		if(hist[i] > max){
			max = hist[0];
			maxIndex = i;
		}
	}
	
	key.mag = max;
	key.ori = maxIndex/(2*CV_PI);
	if(max > m_maxMag){
		m_maxMag = max;
	}

	//Find all of the value that near the max, add the assist keypoint
	double threshold = max*m_param.nearMaxRatio;
	for(int i=0; i < m_param.numOfHistBins; i++){
		if(hist[i] >= threshold){
			SiftKeypoint_t *pAssistKey = (SiftKeypoint_t*)calloc(1, sizeof(SiftKeypoint_t));
			memcpy(pAssistKey, &key, sizeof(SiftKeypoint_t));
			pAssistKey->assistFlag = true;
			pAssistKey->mag = hist[i];
			pAssistKey->ori = i/(2*CV_PI);	
			m_keyPointPool.push_back(*pAssistKey);
		}
	}
}

double* SiftToolbox::CreateOriHist(const IplImage* img, const int x, const int y, const double sigma)
{
	double den1, den2, weight, *hist, mag, ori;
	int radius = m_param.radiusOfOriHistWnd;
	int binIndex = 0;
	den1 = 2*sigma*sigma;
	den2 = sqrt(2*CV_PI)*sigma;
	hist = (double*)calloc(m_param.numOfHistBins, sizeof(double));
	for(int i=-radius; i<=radius; i++){
		for(int j=-radius; j<=radius; j++){
			GetGradMagOri(img, x+i, y+j,&mag, &ori);
			weight = exp((i*i + j*j)/den1)/den2;
			binIndex = cvRound(ori*m_param.numOfHistBins/(CV_PI*2));
			hist[binIndex] += weight*mag;
		}
	}
	return hist;
}

void SiftToolbox::GetGradMagOri(const IplImage* img, const int x, const int y, double* mag, double* ori) const
{
	double Dxx = GetVal32f(img, x+1,y) - GetVal32f(img, x-1, y);
	double Dyy = GetVal32f(img, x, y+1) - GetVal32f(img, x, y-1);
	*mag = sqrt(Dxx*Dxx + Dyy*Dyy);
	*ori = atan2(Dyy, Dxx) + CV_PI; //Change the radius range [0, 2*pi]
}

IplImage* SiftToolbox::PlotKeypoint()
{
	IplImage* img = cvCreateImage(cvGetSize(m_pyrBaseImage), m_pyrBaseImage->depth, m_pyrBaseImage->nChannels);
	img = cvCloneImage(m_pyrBaseImage);
	int	maxMagLen = min(img->width, img->height)/5;
	CvPoint	endPoint;
	for(list<SiftKeypoint_t>::iterator pIter = m_keyPointPool.begin(); pIter != m_keyPointPool.end(); pIter++){
		endPoint.x = pIter->pos.x + cvRound(cos(pIter->ori)*pIter->mag);
		endPoint.y = pIter->pos.y + cvRound(sin(pIter->ori)*pIter->mag);
		cvLine(img, pIter->pos, endPoint, CV_RGB(255,0,0), 1, 8, 0);
	}
	return img;
}