
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
}SiftParam_t;

class SiftToolbox
{
public:
	SiftToolbox();
	~SiftToolbox();
	IplImage* Process(IplImage* image);
	void SetParam(SiftParam_t& param);
	IplImage* PlotKeypoint();
private:
	void	InitImage();
	bool	BuildDogPyr();
	void	FindExtremePoint();
	bool	IsExtrema(const int octave, const int scale, const int x, const int y) const;
	bool	IsRemovableForLowContrast(IplImage* img, int row, int col) const;
	bool	IsRemovableForEdge(IplImage* img, int row, int col) const;
	void	GetGradMagOri(const IplImage* img, const int x, const int y, double* mag, double* ori) const;
	void	OrientationAssignment();
	void	CalcDormOri(SiftKeypoint_t& key);
	double* CreateOriHist(const IplImage* img, const int x, const int y, const double sigma);

	void	DebugInfo();
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