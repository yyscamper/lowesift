
#ifndef __SIFT_UTIL_H__
#define __SIFT_UTIL_H__

#include "cv.h"

 float GetVal32f(const IplImage* img, int x, int y);
 void Setval32f(IplImage* img, int x, int y, float val);
 unsigned char GetVal8(const IplImage* img, int x, int y);
 void Setval32f(IplImage* img, int x, int y, unsigned char val);
#endif