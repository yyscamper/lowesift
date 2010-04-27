
#include "util.h"

float GetVal32f(IplImage* img, int x, int y)
{
	return ((float*)(img->imageData + img->widthStep*x))[y];
}

inline void Setval32f(IplImage* img, int x, int y, float val)
{
	((float*)(img->imageData + img->widthStep*x))[y] = val;
}

inline unsigned char GetVal8(IplImage* img, int x, int y)
{
	return ((unsigned char*)(img->imageData + img->widthStep*x))[y];
}

inline void Setval32f(IplImage* img, int x, int y, unsigned char val)
{
	((unsigned char*)(img->imageData + img->widthStep*x))[y] = val;
}