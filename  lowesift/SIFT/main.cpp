
#include "sift.h"
#include <time.h>
#include <stdio.h>

using namespace std;
void GetBatchSiftFeature(char** imgNameList, int num);
void SaveImage(IplImage* img);

int main()
{
	char* imgNameList[3];
	imgNameList[0] = "food1.jpg";
	imgNameList[1] = "food2.jpg";
	imgNameList[2] = "food3.jpg";
	GetBatchSiftFeature(imgNameList, 3);
	return 0;
}

void GetBatchSiftFeature(char** imgNameList, int num)
{
	IplImage** siftImgArr = (IplImage**)calloc(num, sizeof(IplImage*));
	assert(siftImgArr != NULL);
	IplImage* img;
	SiftToolbox sift;
	for(int i=0; i<num; i++){
		img = cvLoadImage(imgNameList[i], CV_LOAD_IMAGE_COLOR);
		if(img == NULL){
			printf("\nError: failed to load image: \"%s\"", imgNameList[i]);
			return;
		}
		sift.Process(img);
		siftImgArr[i] = sift.PlotKeypoint(SIFT_PLOT_ECLLIPSE);
	}
	
	char wndName[] = "SIFT Feature 1";
	for(int i=0; i<num; i++){
		cvNamedWindow(wndName);
		cvShowImage(wndName, siftImgArr[i]);
		wndName[13]++;
	}
	cvWaitKey(0);
	for(int i=0; i<num; i++){
		cvReleaseImage(&siftImgArr[i]);
	}
	free(siftImgArr);
	cvReleaseImage(&img);
	cvDestroyAllWindows();
}

void SaveImage(IplImage* img)
{
	//Save the image
	char saveImageName[128];
	time_t now;
	now = time(NULL);
	struct tm *tmt = localtime(&now); 
	strftime(saveImageName, 128, "%m%d.%H%M%S.jpg", tmt);
	if(!cvSaveImage(saveImageName, img)){
		printf("Error: Failed to save image.");
		system("pause");
		exit(0);
	}
}