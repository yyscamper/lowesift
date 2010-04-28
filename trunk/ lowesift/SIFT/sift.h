
#ifndef __SIFT_SIFT_H__
#define __SIFT_SIFT_H__

#include "const.h"
#include <cv.h>
#include <list>
using namespace std;

typedef struct{

}SiftFeature_t;

typedef struct{
	int octave;
	double scale;
	double orientation;
	CvPoint pos;
}SiftKeypoint_t;

typedef struct{
	int		numOfOctave;
	int		numOfScalePerOctave;
	bool	isPriorDouble;
	double	sigmaOfPriorDouble;
	double	sigmaOfInitGaussPyr;
	double	thresholdOfKeypointContrast;
	double	ratioOfEdge;
}SiftParam_t;

class SiftToolbox
{
public:
	SiftToolbox();
	~SiftToolbox();
	bool InitImage(IplImage* image);
	IplImage* Process();
	void SetParam(SiftParam_t& param);
private:
	void	InitImage();
	bool	BuildDogPyr();
	void	FindExtremePoint();
	bool	IsExtrema(const int octave, const int scale, const int x, const int y) const;
	bool	IsRemovableForLowContrast(IplImage* img, int row, int col) const;
	bool	IsRemovableForEdge(IplImage* img, int row, int col) const;
	void	GetGradMagOri(const IplImage* img, const int x, const int y, double* mag, double* ori) const;
	void	DebugInfo();
private:
	IplImage* m_baseImage;
	IplImage*** m_dogPyr;
	SiftParam_t m_param;
	list<SiftKeypoint_t> m_keyPointPool;
};

#endif