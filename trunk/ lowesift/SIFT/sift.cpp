
#include "sift.h"
//#include "util.h"
#include "const.h"
#include "assert.h"
#include <stdlib.h>
using namespace std;

SiftToolbox::SiftToolbox()
{
	m_param.numOfOctave = 0;
	m_param.topPyrSize = 16;
	m_param.isPriorDouble = false;
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
	double sigma = 0;

	m_param.numOfOctave = cvRound((log((double)min(m_pyrBaseImage->width, m_pyrBaseImage->height)) 
		- log((double)m_param.topPyrSize))/log(2.0));

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
	m_pSigmaVal = (double*)calloc(dogScaleNum, sizeof(double));
	assert(m_pSigmaVal != NULL);

	gaussPyr[0][0] = cvCloneImage(m_pyrBaseImage);
	for(int octave = 1; octave < m_param.numOfOctave; octave++){
		gaussPyr[octave][0] = cvCreateImage(cvSize(gaussPyr[octave-1][0]->width/2, gaussPyr[octave-1][0]->height/2),
			gaussPyr[octave-1][0]->depth, gaussPyr[octave-1][0]->nChannels);
		cvResize(gaussPyr[octave-1][0], gaussPyr[octave][0], CV_INTER_LINEAR);
	}
	
	double k = pow(2.0, 1.0/m_param.numOfScalePerOctave);
	m_pSigmaVal[0] = m_param.sigmaOfInitGaussPyr;
	for(int i=1; i<dogScaleNum; i++){
		m_pSigmaVal[i] = m_pSigmaVal[i-1]*k;
	}
	for(int octave = 0; octave < m_param.numOfOctave; octave++){
		CvSize currOctaveSize = cvGetSize(gaussPyr[octave][0]);
		for(int scale = 1; scale < gaussScaleNum; scale++){
			gaussPyr[octave][scale] = cvCreateImage(currOctaveSize, IPL_DEPTH_32F, 1);
			gaussPyr[octave][scale] = cvCloneImage(gaussPyr[octave][0]);
			cvSmooth(gaussPyr[octave][scale], gaussPyr[octave][scale], 
				CV_GAUSSIAN, 0, 0, m_pSigmaVal[scale-1], m_pSigmaVal[scale-1]);

			m_dogPyr[octave][scale-1] = cvCreateImage(currOctaveSize, IPL_DEPTH_32F, 1);
			cvSub(gaussPyr[octave][scale-1], gaussPyr[octave][scale], m_dogPyr[octave][scale-1], NULL);
		}
	}
	
	//for(int octave = 0; octave < m_param.numOfOctave; octave++){
	//	for(int scale = 0; scale < m_param.numOfScalePerOctave + 2; scale++){
	//		ShowImage(m_dogPyr[octave][scale]);
	//	}
	//}

	free(pSigmaArry);
	for(int i=0; i<m_param.numOfOctave; i++){
		for(int j=0; j<m_param.numOfScalePerOctave+3; j++){
			cvReleaseImage(&gaussPyr[i][j]);
		}
		free(gaussPyr[i]);
	}
	free(gaussPyr);
	return true;
}

void SiftToolbox::FindExtremePoint()
{
	for(int i=0; i<m_param.numOfOctave; i++){
		for(int j=1; j<m_param.numOfScalePerOctave+1; j++){
			CvSize scaleSize = cvGetSize(m_dogPyr[i][j]);
			for(int x = 4; x < (scaleSize.width-4); x++){
				for(int y = 4; y< (scaleSize.height-4); y++){
					//printf("\nOctave = %d, Scale = %d, x = %d, y = %d. ", i, j, x, y);
					if(IsExtrema(i,j,x,y) 
						&& !IsRemovableForLowContrast(i,j, x, y)
						&& !IsRemovableForEdge(m_dogPyr[i][j], x, y)){
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
	if(val > 0){ //Check Max Value
		for(int i = -1; i<=1; i++){
			for(int j = -1; j<=1; j++){
				for(int k = -1; k<= 1; k++){
					if(val < GetVal32f(m_dogPyr[octave][scale+i], x+j, y+k)){
						return false;
					}
				}
			}
		}
	}else 
	{ //Check Min Value
		for(int i = -1; i<=1; i++){
			for(int j = -1; j<=1; j++){
				for(int k = -1; k<= 1; k++){
					if(val > GetVal32f(m_dogPyr[octave][scale+i], x+j, y+k)){
						return false;
					}
				}
			}
		}
	}
	/*else{
		return false;
	}*/
	//printf("\nFind a extream point.");
	return true;
}

CvMat* SiftToolbox::GetDrivate_1(int octave, int scale, int x, int y) const
{
	double dx, dy, ds;
	dx = (GetVal32f(m_dogPyr[octave][scale], x+1, y) - GetVal32f(m_dogPyr[octave][scale], x-1, y))/2;
	dy = (GetVal32f(m_dogPyr[octave][scale], x, y+1) - GetVal32f(m_dogPyr[octave][scale], x, y-1))/2;
	ds = (GetVal32f(m_dogPyr[octave][scale+1], x, y) - GetVal32f(m_dogPyr[octave][scale-1], x, y))/2;
	CvMat* mat = cvCreateMat(3,1,CV_32FC1);
	cvmSet(mat, 0, 0, dx);
	cvmSet(mat, 1, 0, dy);
	cvmSet(mat, 2, 0, ds);
	return mat;
}

CvMat* SiftToolbox::GetDrivate_2(int octave, int scale, int x, int y) const
{
	CvMat* mat;
	double val, dxx, dyy, dss, dxy, dxs, dys;

	val = GetVal32f(m_dogPyr[octave][scale], x, y);
	dxx = (GetVal32f(m_dogPyr[octave][scale], x+1, y) + GetVal32f(m_dogPyr[octave][scale], x-1, y) - 2*val);
	dyy = (GetVal32f(m_dogPyr[octave][scale], x, y+1) + GetVal32f(m_dogPyr[octave][scale], x, y-1) - 2*val);
	dss = (GetVal32f(m_dogPyr[octave][scale+1], x, y) + GetVal32f(m_dogPyr[octave][scale-1], x, y) - 2*val);
	dxy = ( GetVal32f(m_dogPyr[octave][scale], x+1, y+1 ) + GetVal32f(m_dogPyr[octave][scale], x-1, y-1)
		  - GetVal32f(m_dogPyr[octave][scale], x+1, y-1 ) - GetVal32f(m_dogPyr[octave][scale], x-1, y+1) 
		  ) / 4.0;
	dxs = ( GetVal32f(m_dogPyr[octave][scale+1], x+1, y) + GetVal32f(m_dogPyr[octave][scale-1], x-1, y)
		- GetVal32f(m_dogPyr[octave][scale-1], x+1, y) - GetVal32f(m_dogPyr[octave][scale+1], x-1, y) 
		) / 4.0;
	dys = ( GetVal32f(m_dogPyr[octave][scale+1], x, y+1) + GetVal32f(m_dogPyr[octave][scale-1], x, y-1)
		- GetVal32f(m_dogPyr[octave][scale-1], x, y+1) - GetVal32f(m_dogPyr[octave][scale+1], x, y-1) 
		) / 4.0;

	mat = cvCreateMat( 3, 3, CV_32FC1 );
	cvmSet( mat, 0, 0, dxx );
	cvmSet( mat, 0, 1, dxy );
	cvmSet( mat, 0, 2, dxs );
	cvmSet( mat, 1, 0, dxy );
	cvmSet( mat, 1, 1, dyy );
	cvmSet( mat, 1, 2, dys );
	cvmSet( mat, 2, 0, dxs );
	cvmSet( mat, 2, 1, dys );
	cvmSet( mat, 2, 2, dss );
	return mat;
}

bool SiftToolbox::IsRemovableForLowContrast(int octave, int scale, int x, int y) const
{
	double Dx;
	CvMat* D1 = GetDrivate_1(octave, scale, x, y);
	CvMat* D2 = GetDrivate_2(octave, scale, x, y);
	CvMat* invD2 = cvCreateMat(3, 3, CV_32FC1);
	cvInvert(D2, invD2, CV_SVD);
	CvMat* X = cvCreateMat(3, 1, CV_32FC1);
	cvMul(invD2, D1, X, -1.0);
	int Xx = cvRound(cvmGet(X, 0, 0));
	int Xy = cvRound(cvmGet(X, 1, 0));
	int Xs = cvRound(cvmGet(X, 2, 0));

	CvMat* temp = cvCreateMat(1, 1, CV_32FC1);
	CvMat* tran = cvCreateMat(1, 3, CV_32FC1);
	cvTranspose(D1, tran);
	cvMul(tran, X, temp);
	Dx = GetVal32f(m_dogPyr[octave][Xs],Xx,Xy) + cvmGet(temp, 0, 0);

	cvTranspose(X, tran);
	cvMul(D2, X, X, 1);
	cvMul(tran, X, temp, 0.5);
	Dx += cvmGet(temp, 0, 0);

	cvReleaseMat(&D1);
	cvReleaseMat(&D2);
	cvReleaseMat(&invD2);
	cvReleaseMat(&X);
	cvReleaseMat(&temp);
	cvReleaseMat(&tran);

	if(abs(Dx) <= m_param.thresholdOfKeypointContrast){
		return true;
	}
	return false;
}

bool SiftToolbox::IsRemovableForEdge(IplImage* img, int x, int y) const
{
	double Dxx = GetVal32f(img, x, y+1) + GetVal32f(img, x, y-1) - 2*GetVal32f(img, x, y);
	double Dyy = GetVal32f(img, x-1, y) + GetVal32f(img, x+1, y) - 2*GetVal32f(img, x, y);
	double Dxy = (GetVal32f(img, x+1, y+1) + GetVal32f(img, x-1, y-1)
					- GetVal32f(img, x+1, y-1) - GetVal32f(img, x-1, y+1))/4.0;
	double trace = Dxx + Dyy;
	double det = Dxx*Dyy - Dxy*Dxy;
	if(det <= 0 || trace*trace/det >= pow(m_param.ratioOfEdge+1,2)/m_param.ratioOfEdge){
		return true;
	}
	return false;	
}


void SiftToolbox::OrientationAssignment()
{
	list<SiftKeypoint_t>::size_type total = m_keyPointPool.size(), i = 0;

	list<SiftKeypoint_t>::iterator pCurrKey = m_keyPointPool.begin();
	for(; i<total; pCurrKey++, i++){
		CalcDormOri(*pCurrKey);
		//printf("\n%d, Keypoint: mag=%.4f, ori=%.4f, x=%d", i, pCurrKey->mag, pCurrKey->ori, pCurrKey->pos.x);
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
			max = hist[i];
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
		if(i != maxIndex && hist[i] >= threshold){
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

	if(x+1 >= img->width || y+1 >= img->height || x-1 <0 || y-1 < 0){
		*mag = 0;
		*ori = 0;
		return;
	}

	double Dxx = GetVal32f(img, x+1,y) - GetVal32f(img, x-1, y);
	double Dyy = GetVal32f(img, x, y+1) - GetVal32f(img, x, y-1);
	*mag = sqrt(Dxx*Dxx + Dyy*Dyy);
	*ori = atan2(Dyy, Dxx) + CV_PI; //Change the radius range [0, 2*pi]
}

IplImage* SiftToolbox::PlotKeypoint()
{
	IplImage* img = cvCreateImage(cvGetSize(m_pyrBaseImage), m_pyrBaseImage->depth, m_pyrBaseImage->nChannels);
	img = cvCloneImage(m_originImage);
	int	maxMagLen = min(img->width, img->height)/2;
	CvPoint	startPoint;
	double adjustVal = 0.0;
	int	axesLength = 0;
	for(list<SiftKeypoint_t>::iterator pIter = m_keyPointPool.begin(); pIter != m_keyPointPool.end(); pIter++){
		if(m_param.isPriorDouble){
			adjustVal =  pow(2.0, (double)pIter->octave-1);
		}else{
			adjustVal =  pow(2.0, (double)pIter->octave);
		}
		startPoint.x = cvRound(pIter->pos.x * adjustVal);
		startPoint.y = cvRound(pIter->pos.y * adjustVal);
		axesLength = cvRound(pIter->mag * maxMagLen / m_maxMag);
		cvEllipse(img, startPoint, cvSize(axesLength, axesLength/4), pIter->ori*360/2*CV_PI, 0, 360, CV_RGB(255,0,0), 1, 8, 0);
		//endPoint.x = pIter->pos.x + cvRound(cos(pIter->ori)* pIter->mag * maxMagLen / m_maxMag);
		//endPoint.y = pIter->pos.y + cvRound(sin(pIter->ori)* pIter->mag * maxMagLen / m_maxMag);
		//cvLine(img, startPoint, endPoint, CV_RGB(255,0,0), 1, 8, 0);
	}
	return img;
}

inline float SiftToolbox::GetVal32f(const IplImage* img, int x, int y) const
{
	return ((float*)(img->imageData + img->widthStep*y))[x];
}

inline void SiftToolbox::Setval32f(IplImage* img, int x, int y, float val)
{
	((float*)(img->imageData + img->widthStep*y))[x] = val;
}

inline unsigned char SiftToolbox::GetVal8(const IplImage* img, int x, int y) const
{
	return ((unsigned char*)(img->imageData + img->widthStep*y))[x];
}

inline void SiftToolbox::Setval32f(IplImage* img, int x, int y, unsigned char val)
{
	((unsigned char*)(img->imageData + img->widthStep*y))[x] = val;
}

void SiftToolbox::ShowImage(IplImage* img)
{
	static bool winExistFlag = false;
	if(!winExistFlag){
		cvNamedWindow("Show Image");
		winExistFlag = true;
	}
	cvShowImage("Show Image", img);
	cvWaitKey(0);
}