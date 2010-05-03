
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
	m_param.numOfHistBinsOfDormOri = 36;
	m_param.numOfScalePerOctave = 3;
	m_param.sigmaOfOriAssign = 1.5;
	m_param.sigmaOfInitGaussPyr = 1.0;
	m_param.sigmaOfPriorDouble = 1.6;
	m_param.ratioOfEdge = 10;
	m_param.radiusOfOriHistWnd = cvRound(m_param.sigmaOfOriAssign * 3);
	m_param.thresholdOfKeypointContrast = 0.03;
	m_param.numOfMarginPixel = 5;
	m_param.maxNumOfInterp = 4;
	m_param.numOfDescriptorInRowOrCol = 4; //4*4 Descriptor
	m_param.sizeOfEachDescriptorWnd = 4;
	m_param.numOfMarginPixel = m_param.sizeOfEachDescriptorWnd*m_param.numOfDescriptorInRowOrCol/2 + 1;
	m_param.numOfDescirptorOri = 8;

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
	DescriptorRepresentation();
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
	double sigma = 0;

	m_param.numOfOctave = cvRound((log((double)min(m_pyrBaseImage->width, m_pyrBaseImage->height)) 
		- log((double)m_param.topPyrSize))/log(2.0));

	//allocate memory for the pyramid
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
	
	//calculate the sigma for each scale
	double k = pow(2.0, 1.0/m_param.numOfScalePerOctave);
	m_pSigmaVal[0] = m_param.sigmaOfInitGaussPyr;
	for(int i=1; i<dogScaleNum; i++){
		m_pSigmaVal[i] = m_pSigmaVal[i-1]*k;
	}

	//build the pyramid
	gaussPyr[0][0] = cvCreateImage(cvGetSize(m_pyrBaseImage), m_pyrBaseImage->depth, m_pyrBaseImage->nChannels);
	gaussPyr[0][0] = cvCloneImage(m_pyrBaseImage);
	for(int octave = 0; octave < m_param.numOfOctave; octave++){
		//the image size of current octave
		CvSize currOctaveSize;
		if(octave == 0){
			currOctaveSize = cvGetSize(gaussPyr[0][0]);
		}else{
			currOctaveSize = cvSize(gaussPyr[octave-1][0]->width/2, gaussPyr[octave-1][0]->height/2);
		}

		for(int scale = 0; scale < gaussScaleNum; scale++){
			if(scale ==0 && octave == 0){
				continue;
			}
			gaussPyr[octave][scale] = cvCreateImage(currOctaveSize, IPL_DEPTH_32F, 1);
			if(scale == 0 && octave != 0){ //down sample the previous octave image
				gaussPyr[octave][0] = cvCreateImage(currOctaveSize,
					gaussPyr[octave-1][0]->depth, gaussPyr[octave-1][0]->nChannels);
				cvResize(gaussPyr[octave-1][m_param.numOfScalePerOctave], gaussPyr[octave][0], CV_INTER_LINEAR);
				continue;
			}

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
			for(int x = 4; x < (scaleSize.width - m_param.numOfMarginPixel); x++){
				for(int y = 4; y< (scaleSize.height- m_param.numOfMarginPixel); y++){
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

inline bool SiftToolbox::IsRemovableForLowContrast(int octave, int scale, int x, int y) const
{
	//printf("\nOct: %d, Scale: %d, x: %d, y: %d", octave, scale, x, y);
	CvMat* D1 = GetDrivate_1(octave, scale, x, y);
	CvMat* D2 = GetDrivate_2(octave, scale, x, y);
	CvMat* invD2 = cvCreateMat(3, 3, CV_32FC1);
	CvMat* X = cvCreateMat(3, 1, CV_32FC1);
	int loopCnt = 0, Xx, Xy, Xs;
	for(; loopCnt < m_param.maxNumOfInterp; loopCnt++){
		//Find the extream position
		D1 = GetDrivate_1(octave, scale, x, y);
		D2 = GetDrivate_2(octave, scale, x, y);
		cvInvert(D2, invD2, CV_SVD);
		cvGEMM(invD2, D1, -1.0, NULL, 0, X, 0);

		Xx = cvRound(cvmGet(X, 0, 0));
		Xy = cvRound(cvmGet(X, 1, 0));
		Xs = cvRound(cvmGet(X, 2, 0));

		if(abs(Xx) < 0.5 && abs(Xy) < 0.5 && abs(Xs) < 0.5){
			break;
		}
		x += cvRound(Xx);
		y += cvRound(Xy);
		scale += cvRound(Xs);

		if(x < m_param.numOfMarginPixel || y < m_param.numOfMarginPixel || scale < 1 
			|| scale >= m_param.numOfScalePerOctave + 1
			|| x >= m_dogPyr[octave][scale]->width - m_param.numOfMarginPixel 
			|| y >= m_dogPyr[octave][scale]->height - m_param.numOfMarginPixel){
				return true;
		}
	}
	if(loopCnt >= m_param.maxNumOfInterp){
		return true;
	}

	CvMat* temp = cvCreateMat(1, 1, CV_32FC1);
	cvGEMM(D1, X, 0.5, NULL, 0, temp, CV_GEMM_A_T);
	double Dx = cvmGet(temp, 0, 0) + GetVal32f(m_dogPyr[octave][scale], x, y);
	//cvGEMM(D2, X, 0.5, NULL, 0, temp1, 0);
	//cvGEMM(D1, X, 1.0, NULL, 0, temp2, CV_GEMM_A_T);
	//cvGEMM(X, temp1, 1.0, temp2, 1.0, temp3, CV_GEMM_A_T);
	//Dx += cvmGet(temp3, 0, 0);

	cvReleaseMat(&D1);
	cvReleaseMat(&D2);
	cvReleaseMat(&invD2);
	cvReleaseMat(&X);
	cvReleaseMat(&temp);

	if(abs(Dx) <= m_param.thresholdOfKeypointContrast/m_param.numOfScalePerOctave){
		return true;
	}
	return false;
}

inline bool SiftToolbox::IsRemovableForEdge(IplImage* img, int x, int y) const
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

inline void SiftToolbox::CalcDormOri(SiftKeypoint_t& key)
{
	double* hist = CreateOriHist(m_dogPyr[key.octave][key.scale], key.pos.x, key.pos.y, m_pSigmaVal[key.scale]);	
	SmoothHist(hist, m_param.numOfHistBinsOfDormOri);

	double max = hist[0];
	int maxCnt = 0;
	int maxIndex = 0;
	//Find the max one
	for(int i=1; i<m_param.numOfHistBinsOfDormOri; i++){
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
	for(int i=0; i < m_param.numOfHistBinsOfDormOri; i++){
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

inline double* SiftToolbox::CreateOriHist(const IplImage* img, const int x, const int y, const double sigma)
{
	double den1, den2, weight, *hist, mag, ori;
	int radius = m_param.radiusOfOriHistWnd;
	int binIndex = 0;
	den1 = 2*sigma*sigma;
	den2 = sqrt(2*CV_PI)*sigma;
	hist = (double*)calloc(m_param.numOfHistBinsOfDormOri, sizeof(double));
	for(int i=-radius; i<=radius; i++){
		for(int j=-radius; j<=radius; j++){
			GetGradMagOri(img, x+i, y+j,&mag, &ori);
			weight = exp((i*i + j*j)/den1)/den2;
			binIndex = cvRound(ori*m_param.numOfHistBinsOfDormOri/(CV_PI*2));
			hist[binIndex] += weight*mag;
		}
	}
	return hist;
}

inline void SiftToolbox::SmoothHist(double* hist, int numOfBins)
{
	double pre = hist[0], temp = 0;
	/*printf("\nBefore smooth hist: \n");
	for(int i=0; i<numOfBins; i++){
		printf("%f, ", hist[i]);
	}*/
	hist[0] = 0.75*hist[0] + 0.25*hist[1];
	for(int i=1; i<numOfBins-1; i++){
		temp = hist[i];
		hist[i] = 0.6*hist[i] + 0.2*pre + 0.2*hist[i+1]; 
		pre = temp;
	}
	hist[numOfBins-1] = 0.75*hist[numOfBins-1] + 0.25*hist[numOfBins-2];
	/*printf("\nAfter smooth hist: \n");
	for(int i=0; i<numOfBins; i++){
		printf("%f, ", hist[i]);
	}*/
}

void SiftToolbox::DescriptorRepresentation(){
	list<SiftKeypoint_t>::iterator pCurrKey = m_keyPointPool.begin();
	for(; pCurrKey != m_keyPointPool.end(); pCurrKey++){
		GetDescriptor(*pCurrKey);
	}
}

void SiftToolbox::BuildGuassWndOfDscrp()
{
	int wndSize = m_param.numOfDescriptorInRowOrCol*m_param.sizeOfEachDescriptorWnd;
	double sigma, den1,den2;
	m_gaussWndOfDscp = (CvMat**)calloc(m_param.numOfScalePerOctave+2, sizeof(CvMat**));
	for(int i=0; i<m_param.numOfScalePerOctave+2; i++){
		m_gaussWndOfDscp[i] = cvCreateMat(wndSize, wndSize, CV_32FC1);
		sigma = wndSize/2;
		den1 = sqrt(2*CV_PI)*sigma;
		den2 = -2*sigma*sigma;
		for(int x= -wndSize; x<wndSize; x++){
			for(int y=-wndSize; y<wndSize; y++){
				cvmSet(m_gaussWndOfDscp[i], x, y, exp((x*x+y*y)/den2)/den1);
			}
		}
	}
}

inline void SiftToolbox::GetDescriptor(SiftKeypoint_t& key)
{
	int wndSize = m_param.sizeOfEachDescriptorWnd;
	CvPoint cord;
	key.pDescriptor = (double*)calloc(m_param.numOfDescriptorInRowOrCol*m_param.numOfDescriptorInRowOrCol
		*m_param.numOfDescirptorOri, sizeof(double));
	assert(key.pDescriptor != NULL);
	int offset = 0;
	for(int i=-m_param.numOfDescriptorInRowOrCol; i< m_param.numOfDescriptorInRowOrCol; i++){
		for(int j=-m_param.numOfDescriptorInRowOrCol; j< m_param.numOfDescriptorInRowOrCol; j++){
			cord.x = key.pos.x + i*m_param.sizeOfEachDescriptorWnd;
			cord.y = key.pos.y + j*m_param.sizeOfEachDescriptorWnd;
			GetDscrpWndVec(key.octave, key.scale, cord, m_param.sizeOfEachDescriptorWnd,
				key.pos, key.pDescriptor+offset);
			offset += m_param.numOfDescirptorOri;
		}
	}
}

inline void SiftToolbox::GetDscrpWndVec(int octave, int scale, CvPoint pt, int wndSize, CvPoint center, double* dst)
{
	double mag, ori, maxMag = 0;
	int numOfOriBins = m_param.numOfDescirptorOri;
	int maxMagPos = 0;
	IplImage* img = m_dogPyr[octave][scale];
	CvMat* gaussWnd = m_gaussWndOfDscp[scale];
	for(int x= pt.x; x < pt.x+wndSize; x++){
		for(int y=pt.y; y<pt.y+wndSize; y++){
			GetGradMagOri(img, x, y, &mag, &ori);
			mag *= cvmGet(gaussWnd, x-center.x+gaussWnd->width/2, y-center.y+gaussWnd->height/2);
			if(mag > maxMag){
				maxMagPos = x-pt.x;
			}
			dst[cvRound(ori*numOfOriBins/CV_PI)] += mag;
		}
	}
	ShiftOriHist(dst, numOfOriBins, maxMagPos);
	ThresholdAndNormalizeDsrp(dst, m_param.numOfDescirptorOri);
}

inline double* SiftToolbox::ShiftOriHist(double* hist, int num, int pos)
{
	if(pos == 0){
		return hist;
	}else{
		double* cpyHist = (double*)calloc(num, sizeof(double));
		memcpy(cpyHist, hist, num*sizeof(double));
		for(int i=pos; i<num; i++){
			hist[i-pos] = cpyHist[i];
		}
		for(int i=0; i<pos; i++){
			hist[i+num-pos] = cpyHist[i];
		}
		free(cpyHist);
	}
	
	return hist;
}

inline void SiftToolbox::ThresholdAndNormalizeDsrp(double* data, int num)
{
	double sum = 0.0;
	for(int i=0; i<num; i++){
		data[i] = data[i] > 0.2 ? 0.2: data[i];
		sum += data[i]*data[i];
	}
	double temp = 1/sqrt(sum);
	for(int i=0; i<num; i++){
		data[i] = data[i]*temp;
	}
}

inline void SiftToolbox::GetGradMagOri(const IplImage* img, const int x, const int y, double* mag, double* ori) const
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

IplImage* SiftToolbox::PlotKeypoint(int type)
{
	IplImage* img = cvCreateImage(cvGetSize(m_pyrBaseImage), m_pyrBaseImage->depth, m_pyrBaseImage->nChannels);
	img = cvCloneImage(m_originImage);
	int	maxMagLen = min(img->width, img->height)/4;
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
		switch (type)
		{
		case SIFT_PLOT_DOT:
			maxMagLen = 10;
			cvCircle(img, startPoint, cvRound(pIter->mag * maxMagLen / m_maxMag), CV_RGB(255, 0, 0), -1, 8, 0);
			break;
		case SIFT_PLOT_LINE:
			CvPoint endPoint;
			endPoint.x = pIter->pos.x + cvRound(cos(pIter->ori)* pIter->mag * maxMagLen / m_maxMag);
			endPoint.y = pIter->pos.y + cvRound(sin(pIter->ori)* pIter->mag * maxMagLen / m_maxMag);
			cvLine(img, startPoint, endPoint, CV_RGB(255,0,0), 1, 8, 0);
			break;
		case SIFT_PLOT_RECT:
			maxMagLen = 10;
			cvCircle(img, startPoint, cvRound(pIter->mag * maxMagLen / m_maxMag), CV_RGB(255, 0, 0), -1, 8, 0);
			break;
		default:
			cvEllipse(img, startPoint, cvSize(axesLength, axesLength/4), pIter->ori*360/(2*CV_PI), 0, 360, CV_RGB(255,0,0), 1, 8, 0);
		}
	}
	return img;
}

IplImage* SiftToolbox::GetMatchImage(SiftToolbox& another)
{
	if(this->GetDescriptorVectorDim() != another.GetDescriptorVectorDim()){
		return NULL;
	}

	IplImage* img1, *img2;
	img1 = this->m_originImage;
	img2 = another.m_originImage;
	list<SiftKeypoint_t>& keyPool1 = this->m_keyPointPool;
	list<SiftKeypoint_t>& keyPool2 = another.m_keyPointPool;

	int width = max(img1->width, img2->width);
	int height = max(img1->width, img2->width);
	int gap = 20;

	IplImage* matchImg = cvCreateImage(cvSize(width*2+gap, height), CV_32FC1, 1);
	for(int i=0; i<img1->width; i++){
		for(int j=0; j<img1->height; j++){
			Setval32f(matchImg, i, j, GetVal32f(img1, i, j));
		}
	}
	for(int i=0; i<img2->width; i++){
		for(int j=0; j<img2->height; j++){
			Setval32f(matchImg, i+width+gap, j, GetVal32f(img2, i, j));
		}
	}

	
	int numOfMatchFind = 100;
	if(numOfMatchFind > keyPool1.size()){
		numOfMatchFind = keyPool1.size();
	}
	SiftKeypoint_t *pMatchKey = NULL;		
	double adjust1, adjust2;
	CvPoint p1, p2;
	list<SiftKeypoint_t>::iterator pCurrKey = keyPool1.begin();
	for(int i=0; i<numOfMatchFind; i++, pCurrKey++){
		if(FindMatchPoint(*pCurrKey, keyPool2, GetDescriptorVectorDim(), pMatchKey)){
			if(m_param.isPriorDouble){
				adjust1 =  pow(2.0, (double)pCurrKey->octave-1);
			}else{
				adjust1 =  pow(2.0, (double)pCurrKey->octave);
			}

			if(another.m_param.isPriorDouble){
				adjust2 =  pow(2.0, (double)pCurrKey->octave-1);
			}else{
				adjust2 =  pow(2.0, (double)pCurrKey->octave);
			}
			
			p1.x = cvRound(pCurrKey->pos.x * adjust1);
			p1.y = cvRound(pCurrKey->pos.y * adjust1);
			p2.x = cvRound(pMatchKey->pos.x * adjust2) + width + gap;
			p2.y = cvRound(pMatchKey->pos.y * adjust2);
			cvLine(matchImg, p1, p2, CV_RGB(255,0,0), 1, 8, 0);
		}
	}

	return matchImg;
}
int SiftToolbox::GetDescriptorVectorDim()
{
	return m_param.numOfDescriptorInRowOrCol * m_param.numOfDescriptorInRowOrCol * m_param.numOfDescirptorOri;
}

inline bool SiftToolbox::FindMatchPoint(SiftKeypoint_t& srcKey, list<SiftKeypoint_t>& srcPool, int dim, SiftKeypoint_t* pMatchKey)
{
	double firstMin=1000000.0, secMin=1000000.0;
	double euDis = 0.0;
	double sum = 0.0;

	list<SiftKeypoint_t>::iterator pCurrKey = srcPool.begin();
	for(; pCurrKey != m_keyPointPool.end(); pCurrKey++){
		for(int i=0; i< dim; i++){
			sum += (srcKey.pDescriptor[i] - pCurrKey->pDescriptor[i]) 
				* (srcKey.pDescriptor[i] - pCurrKey->pDescriptor[i]);
		}
		euDis = sum;//sqrt(sum);
		if(euDis < firstMin){
			firstMin = euDis;
			pMatchKey = &(*pCurrKey);
		}else if(euDis < secMin){
			secMin = euDis;
		}
	}
	if(firstMin/secMin < 0.01){
		return true;
	}	
	return false;
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