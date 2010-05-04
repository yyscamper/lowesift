
/**
 * SIFT特征点的提取
 * 作者：袁宇
 * 时间：2010-05-01
 * 参考文献：
	[1] David G. Lowe, Distinctive Image Features form Scale-Invariant Keypoints, International Journal of Computer Vision, 2004
	[2] David G. Lowe, Object Recognition form Local Scale-Invariant Features, Proc. Of the International Conference on Computer Vision, Corfu, 1999
	[3]赵辉, SIFT特征特征匹配技术讲义
 */

#ifndef __SIFT_SIFT_H__
#define __SIFT_SIFT_H__

#include "const.h"
#include <cv.h>
#include <highgui.h>
#include <list>

typedef struct{
	int		octave;
	int		scale;
	double	ori;
	double	mag;
	CvPoint pos;
	bool	assistFlag;
	double* pDescriptor;
}SiftKeypoint_t;

typedef struct{
	int		topPyrSize;
	int		numOfOctave;
	int		numOfScalePerOctave;
	bool	isPriorDouble;
	double	sigmaOfPriorDouble;
	double	sigmaOfInitGaussPyr;
	double	thresholdOfKeypointContrast;
	double	ratioOfEdge; 
	double	sigmaOfOriAssign; //1.5
	int		radiusOfOriHistWnd; // 3*sigmaOfOriAssign
	int		numOfHistBinsOfDormOri; //36
	double	nearMaxRatio; //80%
	int		maxNumOfInterp;
	int		numOfMarginPixel;
	int		numOfDescriptorInRowOrCol; //4*4 or 2*2
	int		sizeOfEachDescriptorWnd; //
	int		numOfDescirptorOri;
}SiftParam_t;

enum Sift__Plot_Type{
	SIFT_PLOT_DOT = 0,
	SIFT_PLOT_LINE = 1,
	SIFT_PLOT_RECT = 2,
	SIFT_PLOT_ECLLIPSE = 3
};

class SiftToolbox
{
public:
	SiftToolbox();
	~SiftToolbox();
	IplImage* Process(IplImage* image);
	void SetParam(SiftParam_t& param);
	IplImage* PlotKeypoint(int type = SIFT_PLOT_ECLLIPSE);
	IplImage* GetMatchImage(SiftToolbox& another);
	int  GetDescriptorVectorDim();
private:
	void	InitImage();
	bool	BuildDogPyr();
	void	FindExtremePoint();
	bool	IsExtrema(const int octave, const int scale, const int x, const int y) const;
	inline bool	IsRemovableForLowContrast(int ocatve, int scale, int x, int y) const;
	inline bool	IsRemovableForEdge(IplImage* img, int x, int y) const;
	inline void	GetGradMagOri(const IplImage* img, const int x, const int y, double* mag, double* ori) const;
	void	OrientationAssignment();
	inline void	CalcDormOri(SiftKeypoint_t& key);
	inline double* CreateOriHist(const IplImage* img, const int x, const int y, const double sigma);
	inline void	SmoothHist(double* hist, int numOfBins);
	
	void BuildGuassWndOfDscrp();
	void	DescriptorRepresentation();
	inline void GetDescriptor(SiftKeypoint_t& key);
	inline void GetDscrpWndVec(int octave, int scale, CvPoint pt, int wndSize, CvPoint center, double oriOfKey, double* dst);
	inline double* ShiftOriHist(double* hist, int num, int pos);
	inline void ThresholdAndNormalizeDsrp(double* data, int num);

	bool FindMatchPoint(SiftKeypoint_t& srcKey, std::list<SiftKeypoint_t>& srcPool, int dim, SiftKeypoint_t* pMatchKey);

	inline CvMat*	GetDrivate_1(int octave, int scale, int x, int y) const;
	inline CvMat*	GetDrivate_2(int octave, int scale, int x, int y) const;

	inline float GetVal32f(const IplImage* img, int x, int y) const;
	inline void Setval32f(IplImage* img, int x, int y, float val);
	inline unsigned char GetVal8(const IplImage* img, int x, int y) const;
	inline void Setval32f(IplImage* img, int x, int y, unsigned char val);

	void	DebugInfo(char* format, ...);
	void	ShowImage(IplImage* img);
private:
	IplImage*				m_originImage;
	IplImage*				m_pyrBaseImage;
	IplImage***				m_dogPyr;
	double*					m_pSigmaVal;
	SiftParam_t				m_param;
	double					m_maxMag;
	CvMat**					m_gaussWndOfDscp;
	std::list<SiftKeypoint_t>	m_keyPointPool;
};

#endif