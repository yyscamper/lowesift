
/**
 * SIFT特征点的提取
 * 作者：袁宇
 * 时间：2010-05-01
 * 参考文献：
	[1] David G. Lowe, Distinctive Image Features form Scale-Invariant Keypoints, International Journal of Computer Vision, 2004
	[2] David G. Lowe, Object Recognition form Local Scale-Invariant Features, Proc. Of the International Conference on Computer Vision, Corfu, 1999
	[3]赵辉, SIFT特征特征匹配技术讲义
 */

#include "sift.h"
#include <time.h>
#include <stdio.h>

using namespace std;
void GetBatchSiftFeature(char** imgNameList, int num);
void SaveImage(IplImage* img, char* imgName);
void MatchTest(char** imgNameList);
int main()
{
	char* imgNameList[10];

	/*imgNameList[0] = "haibao1.jpg";
	imgNameList[1] = "haibao2.jpg";
	imgNameList[2] = "haibao3.jpg";*/

	//imgNameList[0] = "fuwa1.jpg";
	//imgNameList[1] = "fuwa2.jpg";
	//imgNameList[2] = "fuwa3.jpg";
	//imgNameList[3] = "fuwa4.jpg";
	//imgNameList[4] = "fuwa9.jpg";

	/*imgNameList[0] = "book1.jpg";
	imgNameList[1] = "book2.jpg";
	imgNameList[2] = "book3.jpg";
	imgNameList[3] = "book4.jpg";
	imgNameList[4] = "book4.jpg";*/

	//imgNameList[0] = "2book1.jpg";
	//imgNameList[1] = "2book2.jpg";
	//imgNameList[2] = "2book3.jpg";
	//imgNameList[3] = "2book4.jpg";
	//imgNameList[4] = "2book5.jpg";
	//imgNameList[5] = "2book6.jpg";

	/*imgNameList[0] = "3611.jpg";
	imgNameList[1] = "3612.jpg";
	imgNameList[2] = "3613.jpg";
	imgNameList[3] = "3614.jpg";
	imgNameList[4] = "3615.jpg";*/

	//imgNameList[0] = "tattoo1.jpg";
	//imgNameList[1] = "tattoo2.jpg";
	//imgNameList[2] = "tattoo3.jpg";
	//imgNameList[3] = "tattoo4.jpg";
	//imgNameList[4] = "tattoo5.jpg";
	//imgNameList[5] = "tattoo9.jpg";

	
	//imgNameList[0] = "food1.jpg";
	//imgNameList[1] = "food2.jpg";
	//imgNameList[2] = "food3.jpg";
	//imgNameList[3] = "food4.jpg";

	imgNameList[0] = "deer.jpg";
	GetBatchSiftFeature(imgNameList,1);
	//MatchTest(imgNameList);
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
		siftImgArr[i] = sift.PlotKeypoint(SIFT_PLOT_LINE);
	}
	
	char wndName[] = "SIFT Feature 1";
	for(int i=0; i<num; i++){
		cvNamedWindow(wndName);
		cvShowImage(wndName, siftImgArr[i]);
		wndName[13]++;
	}
	cvWaitKey(0);
	
	char saveImgName[128];
	strcpy(saveImgName, "1_");
	for(int i=0; i<num; i++){
		strcpy(saveImgName+2, imgNameList[i]);
		SaveImage(siftImgArr[i], saveImgName);
		saveImgName[0]++;
	}

	for(int i=0; i<num; i++){
		cvReleaseImage(&siftImgArr[i]);
	}
	free(siftImgArr);
	cvReleaseImage(&img);
	cvDestroyAllWindows();
}

void SaveImage(IplImage* img, char* imgName)
{
	////Save the image
	//char saveImageName[128];
	//time_t now;
	//now = time(NULL);
	//struct tm *tmt = localtime(&now); 
	//strftime(saveImageName, 128, "%m%d.%H%M%S.jpg", tmt);
	if(!cvSaveImage(imgName, img)){
		printf("Error: Failed to save image: %s", imgName);
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