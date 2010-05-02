
#include "sift.h"
#include <stdio.h>

using namespace std;

int main()
{
	char* path = "book.jpg";
	IplImage* img;
	img = cvLoadImage(path, CV_LOAD_IMAGE_COLOR);
	if(img == NULL){
		printf("\nError: failed to load image from : \"%s\"", path);
		exit(0);
	}

	SiftToolbox sift;
	sift.Process(img);
	IplImage* siftImg = sift.PlotKeypoint();
	cvNamedWindow("SIFT Feature");
	cvShowImage("SIFT Feature", siftImg);
	cvWaitKey(0);
	cvReleaseImage(&img);
	cvReleaseImage(&siftImg);
	cvDestroyAllWindows();
	return 0;
}