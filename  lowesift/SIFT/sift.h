
#ifndef __SIFT_SIFT_H__
#define __SIFT_SIFT_H__

#include "const.h"
#include <cv.h>
#include <highgui.h>
#include <list>
using namespace std;

typedef struct{

}SiftFeature_t;

typedef struct{
	int		octave;
	int		scale;
	double	ori;
	double	mag;
	CvPoint pos;
	bool	assistFlag;
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
	int		numOfHistBins; //36
	double	nearMaxRatio; //80%
	int		maxNumOfInterp;
	int		numOfMarginPixel;
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
	list<SiftKeypoint_t>	m_keyPointPool;
	double					m_maxMag;
};

#endif