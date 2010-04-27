
#ifndef __SIFT_UTIL_H__
#define __SIFT_UTIL_H__

#include "cv.h"

inline float GetVal32f(IplImage* img, int x, int y);
inline void Setval32f(IplImage* img, int x, int y, float val);
inline unsigned char GetVal8(IplImage* img, int x, int y);
inline void Setval32f(IplImage* img, int x, int y, unsigned char val);
#endif