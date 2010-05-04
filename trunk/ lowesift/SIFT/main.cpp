
#include "sift.h"
#include <time.h>
#include <stdio.h>

using namespace std;
void GetBatchSiftFeature(char** imgNameList, int num);
void SaveImage(IplImage* img);
void MatchTest(char** imgNameList);
int main()
{
	char* imgNameList[10];

	/*imgNameList[0] = "haibao1.jpg";
	imgNameList[1] = "haibao2.jpg";
	imgNameList[2] = "haibao3.jpg";*/

	/*imgNameList[0] = "fuwa1.jpg";
	imgNameList[1] = "fuwa2.jpg";
	imgNameList[2] = "fuwa3.jpg";*/

	/*imgNameList[0] = "book1.jpg";
	imgNameList[1] = "book2.jpg";
	imgNameList[2] = "book3.jpg";
	imgNameList[3] = "book4.jpg";*/

	imgNameList[0] = "3611.jpg";
	imgNameList[1] = "3612.jpg";
	imgNameList[2] = "3613.jpg";
	imgNameList[3] = "3614.jpg";

	//GetBatchSiftFeature(imgNameList, 4);
	MatchTest(imgNameList);
	return 0;
}

void GetBatchSiftFeature(char** imgNameList, int num)
{
	IplImage** siftImgArr = (IplImage**)calloc(num, sizeof(IplImage*));
	assert(siftImgArr != NULL);
	IplImage* img;
	for(int i=0; i<num; i++){
		SiftToolbox sift;
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

void MatchTest(char** imgNameList)
{
	SiftToolbox sift1, sift2;
	IplImage *img1, *img2;
	//process the first image
	img1 = cvLoadImage(imgNameList[0], CV_LOAD_IMAGE_COLOR);
	if(img1 == NULL){
		printf("\nError: failed to load image: \"%s\"", imgNameList[0]);
		return;
	}
	sift1.Process(img1);

	//process the second image
	img2 = cvLoadImage(imgNameList[1], CV_LOAD_IMAGE_COLOR);
	if(img2 == NULL){
		printf("\nError: failed to load image: \"%s\"", imgNameList[1]);
		return;
	}
	sift2.Process(img2);

	IplImage* matchImage = sift1.GetMatchImage(sift2);
	
	char wndName[] = "SIFT Match Image";
	cvNamedWindow(wndName);
	cvShowImage(wndName, matchImage);
	cvWaitKey(0);
	cvDestroyAllWindows();
	cvReleaseImage(&img1);
	cvReleaseImage(&img2);
	cvReleaseImage(&matchImage);
	
}