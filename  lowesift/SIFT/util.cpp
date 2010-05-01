
#include "util.h"

 float GetVal32f(const IplImage* img, int x, int y)
{
	return ((float*)(img->imageData + img->widthStep*x))[y];
}

 void Setval32f(IplImage* img, int x, int y, float val)
{
	((float*)(img->imageData + img->widthStep*x))[y] = val;
}

 unsigned char GetVal8(const IplImage* img, int x, int y)
{
	return ((unsigned char*)(img->imageData + img->widthStep*x))[y];
}

 void Setval32f(IplImage* img, int x, int y, unsigned char val)
{
	((unsigned char*)(img->imageData + img->widthStep*x))[y] = val;
}