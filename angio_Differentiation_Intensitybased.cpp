/*
Author: Denny,	Date: August 2023

OCT Angiography - Intensity-based Differentiation 

the keys and their purposes are given below   
4 - increments the number of averages by 1 in Average mode
3 - decrements the number of averages by 1 in Average mode
] - increases threshold for display
[ - decreases threshold for display
. - increment n by 1, threshold stepsize being 10^n   
, - decrement n by 1, threshold stepsize being 10^n   
's' - saves live bscan, vibrational profile and stats
Esc or x - quits the program

*/

#include <iostream>
#include <sstream>
#include <opencv2/opencv.hpp>
#include "FDOCT.h"
#include "Spinnaker.h"
#include "SpinGenApi/SpinnakerGenApi.h"
#include <sys/time.h> // gettimeofday()
#include <sys/stat.h> // this is for mkdir
#include <array>
#include <cmath>

using namespace cv;
using namespace std;
using namespace Spinnaker;
using namespace Spinnaker::GenApi;

// Function declarations
void setCamera(CameraPtr pCam);
inline void makeonlypositive(Mat& src, Mat& dst);
void matwrite(const string& filename, const Mat& mat);
inline void savematasbin(char* p, char* d, char* f, Mat m);
inline void savematasimage(char* p, char* d, char* f, Mat m);

int main()
{
	int result = 0;
	bool bgframestaken = false, accummode = false;
	bool expchanged = false, skeypressed = false, doneflag = false, spacebarpressed = false;
	bool dir_created = false;
	double minval, maxval, minvalsum, maxvalsum, meanvalsum, minvaldiff, maxvaldiff, meanvaldiff, bscanthreshold = 0.0, threshold_N = 1.0;
	int dt, key;
	double fps = 0.0, framecount = 0.0;
	double Vmax = 0.0, Vmean = 0.0, Vmin, Vacc = 0.0, Vavg; // vibration amplitude
	double vibreadcount = 1.0;
	system("clear");
	// Print application build information
	cout << "Application build date: " << __DATE__ << " " << __TIME__ << endl << endl;

	SystemPtr system = System::GetInstance();
    
	// Retrieve list of cameras from the system
    CameraList camList = system->GetCameras();
    
	unsigned int numCameras = camList.GetSize();
    cout << "Number of cameras detected: " << numCameras << endl << endl;
	if (numCameras == 0)
    {
        // Clear camera list before releasing system
        camList.Clear();
        // Release system
        system->ReleaseInstance();
        cout << "Camera not detected. " << endl;
        cout << "Done! Press Enter to exit..." << endl;
        getchar();
        return -1;
    }
	CameraPtr pCam = nullptr;
	pCam = camList.GetByIndex(0);

	Mat statusimg = Mat::zeros(cv::Size(600, 300), CV_64F);
	Mat firstrowofstatusimg = statusimg(Rect(0, 0, 600, 50)); // x,y,width,height
	Mat secrowofstatusimg = statusimg(Rect(0, 50, 600, 50));
	Mat secrowofstatusimgRHS = statusimg(Rect(300, 50, 300, 50));
	Mat thirdrowofstatusimg = statusimg(Rect(0, 100, 600, 50));
	Mat fourthrowofstatusimg = statusimg(Rect(0, 150, 600, 50));
	char textbuffer[80];

	namedWindow("Interferogram", 0); // 0 = WINDOW_NORMAL
	moveWindow("Interferogram", 0, 80);

	namedWindow("(1) IbDiff", 0); // 0 = WINDOW_NORMAL
	moveWindow("(1) IbDiff", 0, 440);

	namedWindow("Live B-scan", 0); // 0 = WINDOW_NORMAL
	moveWindow("Live B-scan", 430, 80);

	namedWindow("Status", 0); // 0 = WINDOW_NORMAL
	moveWindow("Status", 900, 80);

	ifstream infile("spint.ini");
	string tempstring;

	int res;	

	char dirdescr[60];
	sprintf(dirdescr, "_");
	char dirname[80];
	char filename[20];
	char filenamec[20];
	unsigned int indexi = 0;
	char pathname[140];

	char lambdamaxstr[40];
	char lambdaminstr[40];
	double lambdamin, lambdamax;
	unsigned int averagescount, numdisplaypoints, fftpointsmultiplier;
	char thresholdstr[40];
	bool ROIflag = false;
	unsigned int ROIstartrow, ROIendrow, ROIstartcol, ROIendcol, heightROI, widthROI, ROIcentrey, ROIcentrex;
	double ROImeanVal, ROImaxVal;
	unsigned int binvaluex = 1, binvaluey = 1;

	// inputs from ini file
	if (infile.is_open())
	{
		// skip the first 22 lines - they contain camera settings
		for ( int i = 0; i < 22; i++)
		{
			infile >> tempstring;
		}

		infile >> tempstring;
		infile >> lambdaminstr;
		infile >> tempstring;
		infile >> lambdamaxstr;
		infile >> tempstring;
		infile >> averagescount;
		infile >> tempstring;
		infile >> numdisplaypoints;
		infile >> tempstring;
		infile >> fftpointsmultiplier;
		infile >> tempstring;
		infile >> thresholdstr;
		infile >> tempstring;
		infile >> dirdescr;
		infile >> tempstring;
		infile >> ROIstartrow;
		infile >> tempstring;
		infile >> ROIendrow;
		infile >> tempstring;
		infile >> ROIstartcol;
		infile >> tempstring;
		infile >> ROIendcol;
		infile >> tempstring;
		infile >> binvaluex;
		infile >> tempstring;
		infile >> binvaluey;
		
		infile.close();
		lambdamin = atof(lambdaminstr);
		lambdamax = atof(lambdamaxstr);
		bscanthreshold = atof(thresholdstr);
		cout << "lambdamin set to " << lambdamin << " ..." <<endl;
		cout << "lambdamax set to " << lambdamax << " ..." << endl;
	}

	try 
	{
		// Initialize camera
		pCam->Init();
		setCamera(pCam); // set the camera on non-trigerred mode
		
		unsigned int w, h, opw, oph, camtime;
		// double exposure method
		unsigned int camtime1, camtime2;
		camtime1 = 800; camtime2 = 400;
		w = pCam->Width.GetValue();
		h = pCam->Height.GetValue();
		opw = w / binvaluex;
		oph = h / binvaluey;
		cout << "binvaluex set to " << binvaluex << " ..." <<endl;
		cout << "binvaluey set to " << binvaluey << " ..." << endl;
		cout << "width set to " << opw << " ..." <<endl;
		cout << "height set to " << oph << " ..." << endl;
		camtime = pCam->ExposureTime.GetValue();

		////////// Initial computation for camera-based FD-OCT //////////////	
		const double pi = 3.141592653589793;
		Mat data_yb = Mat(oph, opw, CV_64F);// the Mat constructor Mat(rows,columns,type);
		Mat data_yI = Mat(oph, opw, CV_64F);// the Mat constructor Mat(rows,columns,type);
		Mat data_yO = Mat(oph, opw, CV_64F);// the Mat constructor Mat(rows,columns,type);
		double deltalambda = (lambdamax - lambdamin) / opw;
		// create modified Bartlett-Hann window
		Mat barthannwin	= Mat(1, opw, CV_64F);
		for (unsigned int p = 0; p<(opw); p++)
		{
			// https://in.mathworks.com/help/signal/ref/barthannwin.html
			float nn = p;
			float NN = opw - 1;
			barthannwin.at<double>(0, p) = 0.62 - 0.48*abs(nn / NN - 0.5) + 0.38*cos(2 * pi*(nn / NN - 0.5));
		}
		Mat klinear = Mat::zeros(cv::Size(1, opw), CV_64F);
		Mat fractionalk = Mat::zeros(cv::Size(1, opw), CV_64F);
		Mat nearestkindex = Mat::zeros(cv::Size(1, opw), CV_32S);
		Mat lambdas = Mat::zeros(cv::Size(1, opw), CV_64F);		//Size(cols,rows)
		Mat diffk = Mat::zeros(cv::Size(1, opw), CV_64F);
		unsigned int indextemp;
		// compute lambdas
		for (indextemp = 0; indextemp< (opw); indextemp++)
		{
			lambdas.at<double>(0, indextemp) = lambdamin + indextemp * deltalambda;
		}
		Mat k = 2 * pi / lambdas;
		double kmin = 2 * pi / (lambdamax - deltalambda);
		double kmax = 2 * pi / lambdamin;
		double deltak = (kmax - kmin) / opw;
		// compute klinear
		for (indextemp = 0; indextemp < opw; indextemp++)
		{
			klinear.at<double>(0, indextemp) = kmin + (indextemp + 1)*deltak;
		}
		// find the diff of the non-linear ks
		for (indextemp = 1; indextemp < opw; indextemp++)
		{
			// since this is a decreasing series, RHS is (i-1) - (i)
			diffk.at<double>(0, indextemp) = k.at<double>(0, indextemp - 1) - k.at<double>(0, indextemp);
		}
		// and initializing the first point separately
		diffk.at<double>(0, 0) = diffk.at<double>(0, 1);
		// find the index of the nearest k value, less than the linear k
		for (int f = 0; f < opw; f++)
		{
			for (indextemp = 0; indextemp < opw; indextemp++)
			{
				if (k.at<double>(0, indextemp) < klinear.at<double>(0, f))
				{
					nearestkindex.at<int>(0, f) = indextemp;
					break;
				}// end if
			}//end indextemp loop
		}// end f loop
	// now find the fractional amount by which the linearized k value is greater than the next lowest k
		for (int f = 0; f < opw; f++)
		{
			fractionalk.at<double>(0, f) = (klinear.at<double>(0, f) - k.at<double>(0, nearestkindex.at<int>(0, f))) / diffk.at<double>(0, nearestkindex.at<int>(0, f));
		}
		////////// End of initial computation for camera-based FD-OCT //////////////
	
		cout << "Acquiring images " << endl;
		pCam->BeginAcquisition();
		ImagePtr pResultImage;
		ImagePtr convertedImage;

		unsigned long time_start, time_end, tic, toc;	
		struct timeval tv;
		gettimeofday(&tv,NULL);	
		time_start = 1000000 * tv.tv_sec + tv.tv_usec;	

		Mat m, opm, Sk, bscantemp, mvector, bvector, bscanlive, bscanlivebg;
		Mat X1, X2, numerator, denominator, squareX1, squareX2, ADV;
		Mat num, denom, denom1, denom2, X1X2, sqX1, sqX2;
		Mat Fmreal, Fmimag, Fm1real, Fm1imag, Fsub, F, realFdiff, imagFdiff, magFdiff;
		Mat term1, term2, termreal, termReal, termimag, termImag;
		Mat jdiff, positivediff, bscanlog, bscandb, bscandisp;
		Mat Idiff, tempmat1, tempmat2;
		Mat Idiff1, Idiff2, Idiff5;
		Mat IDiff1, IDiff2, IDiff5;
		Mat cmagI, cmagI1;
		Mat ROI, tempmatdiff1, tempmatdiff2, tempmatdiff3, tempmatdiff4, tempmatdiff5;
		Mat bscandb1;
		Mat bscandisp1;	

		Idiff1 = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
		Idiff2 = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
		Idiff5 = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
		num = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
		denom1 = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
		denom2 = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
		termReal = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
		termImag = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);

		// ROI parameters
		Point ROItopleft(ROIstartcol,ROIstartrow), ROIbottright(ROIendcol,ROIendrow);
		enum ROIoption{position=1,size};
		int selectedROIoption;
		selectedROIoption = position;
		// Selected pixel
		Point selectedPixel((ROIstartcol+ROIendcol)/2,(ROIstartrow+ROIendrow)/2);

		enum displayStats4option{IbD=1, IPbD, IbDV, AD};
	
		resizeWindow("Live B-scan", oph, numdisplaypoints);		// (width,height)
		resizeWindow("(1) IbDiff", oph, numdisplaypoints);		// (width,height)

		int ret;
		unsigned int ii; // index
		bool stuck;

		while (1)	//camera frames acquisition loop, which is inside the try
		{
			if(bgframestaken == false)
			{
				res = 0;
				ret = 0;
				// save one image to Mat m
				while(ret == 0)
				{
					pResultImage = pCam->GetNextImage();
					if (pResultImage->IsIncomplete())
					{
						ret = 0;
					}
					else
					{
						ret = 1;
						convertedImage = pResultImage;
						m = Mat(h, w, CV_16UC1, convertedImage->GetData(), convertedImage->GetStride());
						// binning (averaging)
						resize(m, opm, Size(), 1.0 / binvaluex, 1.0 / binvaluey, INTER_AREA);
						opm.convertTo(data_yb, CV_64F);
					}
				}

				// pResultImage has to be released to avoid buffer filling up
				pResultImage->Release();
				imshow("Interferogram", opm);
				if(ret == 1)
				{
					////////////// compute Bscan //////////////////////	
					// DC removal and windowing
					for (int p = 0; p<(data_yb.rows); p++)
					{
						Scalar meanval = mean(data_yb.row(p));
						data_yb.row(p) = data_yb.row(p) - meanval(0);	// Only the first value of the scalar is useful for us
        				// uncomment the line below to use Bartlett-Hann window; commment for no windowing
						multiply(data_yb.row(p), barthannwin, data_yb.row(p));
					}
					Mat slopes = Mat::zeros(cv::Size(data_yb.rows, opw), CV_64F);
					Mat data_ylin(oph, opw, CV_64F);
					// interpolate to linear k space
					for (int p = 0; p<(data_yb.rows); p++)
					{
						for (int q = 1; q<(data_yb.cols); q++)
						{
							//find the slope of the data_y at each of the non-linear ks
							slopes.at<double>(p, q) = data_yb.at<double>(p, q) - data_yb.at<double>(p, q - 1);
							// in the .at notation, it is <double>(y,x)
						}
						// initialize the first slope separately
						slopes.at<double>(p, 0) = slopes.at<double>(p, 1);
						for (int q = 1; q<(data_ylin.cols - 1); q++)
						{
							//find the value of the data_ylin at each of the klinear points
							// data_ylin = data_y(nearestkindex) + fractionalk(nearestkindex)*slopes(nearestkindex)
							data_ylin.at<double>(p, q) = data_yb.at<double>(p, nearestkindex.at<int>(0, q)) + fractionalk.at<double>(nearestkindex.at<int>(0, q)) * slopes.at<double>(p, nearestkindex.at<int>(0, q));
						}
					}	
					
					Mat complexI, magI, bscansub;
					bscantemp = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
					// InvFFT - compute the bscan
					Mat planes[] = { Mat_<float>(data_ylin), Mat::zeros(data_ylin.size(), CV_32F) };
					merge(planes, 2, complexI);       // Add to the expanded another plane with zeros
					dft(complexI, complexI, DFT_ROWS | DFT_INVERSE);
					split(complexI, planes);          // planes[0] = Re(DFT(I)), planes[1] = Im(DFT(I))
					magnitude(planes[0], planes[1], magI);
					bscansub = magI.colRange(0, numdisplaypoints);
					bscansub.convertTo(bscantemp, CV_64F);
					accumulate(bscansub,bscantemp);
					
					////////////// end of compute Bscan //////////////////////	
					bscantemp.copyTo(bscanlivebg);
				} // end of if ret == 1 block 
				
				bgframestaken = true;				
			} // end of if(bgframestaken == false)

			for(ii = 1; ii <= averagescount; ii++)
			{
				if(ii % 2 == 1) // odd frame
				{
					if(IsReadable(pCam->ExposureTime) && IsWritable(pCam->ExposureTime))
					{
						pCam->ExposureTime.SetValue(camtime1);
					}
				}
				else // even frame
				{
					if(IsReadable(pCam->ExposureTime) && IsWritable(pCam->ExposureTime))
					{
						pCam->ExposureTime.SetValue(camtime2);
					}
				}
				stuck = false;
				ret = 0;
				usleep(3300);	
				// save one image to Mat m
				gettimeofday(&tv,NULL);
				tic = tv.tv_usec + tv.tv_sec * 1000000;	
				while(ret == 0)
				{
					pResultImage = pCam->GetNextImage();
					if (pResultImage->IsIncomplete())
					{   
						ret = 0;
					}
					else
					{
						gettimeofday(&tv,NULL);
						toc = tv.tv_usec + tv.tv_sec * 1000000;	
						if( (toc-tic) >= 10000 ) // 10 ms 
						{
							//cout << "Stuck1 ... breaking from loop " << endl;
							stuck = true;
							break;
						}
						else
						{	
							convertedImage = pResultImage;
							m = Mat(h, w, CV_16UC1, convertedImage->GetData(), convertedImage->GetStride());
							// binning (averaging)
							resize(m, opm, Size(), 1.0 / binvaluex, 1.0 / binvaluey, INTER_AREA);
							opm.convertTo(data_yI, CV_64F);
						}
						ret = 1;
					}
				}

				// pResultImage has to be released to avoid buffer filling up
				pResultImage->Release();
				if(stuck == true)
				{
					ii = ii -1;
					continue;
				}

				imshow("Interferogram", opm);
				if(ret == 1)
				{
					////////////// compute Bscan //////////////////////	
					// DC removal and windowing
					for (int p = 0; p<(data_yI.rows); p++)
					{
						Scalar meanval = mean(data_yI.row(p));
						data_yI.row(p) = data_yI.row(p) - meanval(0);	// Only the first value of the scalar is useful for us
        				// uncomment the line below to use Bartlett-Hann window; commment for no windowing
						multiply(data_yI.row(p), barthannwin, data_yI.row(p));
					}
					Mat slopes = Mat::zeros(cv::Size(data_yI.rows, opw), CV_64F);
					Mat data_ylin(oph, opw, CV_64F);
					// interpolate to linear k space
					for (int p = 0; p<(data_yI.rows); p++)
					{
						for (int q = 1; q<(data_yI.cols); q++)
						{
							//find the slope of the data_y at each of the non-linear ks
							slopes.at<double>(p, q) = data_yI.at<double>(p, q) - data_yI.at<double>(p, q - 1);
							// in the .at notation, it is <double>(y,x)
						}
						// initialize the first slope separately
						slopes.at<double>(p, 0) = slopes.at<double>(p, 1);
						for (int q = 1; q<(data_ylin.cols - 1); q++)
						{
							//find the value of the data_ylin at each of the klinear points
							// data_ylin = data_y(nearestkindex) + fractionalk(nearestkindex)*slopes(nearestkindex)
							data_ylin.at<double>(p, q) = data_yI.at<double>(p, nearestkindex.at<int>(0, q)) + fractionalk.at<double>(nearestkindex.at<int>(0, q)) * slopes.at<double>(p, nearestkindex.at<int>(0, q));
						}
					}
				
					Mat complexI, magI, bscansub;
					bscantemp = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
					Fmreal = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
					Fmimag = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
					// InvFFT - compute the bscan
					Mat planes[] = { Mat_<float>(data_ylin), Mat::zeros(data_ylin.size(), CV_32F) };
					merge(planes, 2, complexI);       // Add to the expanded another plane with zeros
					dft(complexI, complexI, DFT_ROWS | DFT_INVERSE);
					split(complexI, planes);          // planes[0] = Re(DFT(I)), planes[1] = Im(DFT(I))
					magnitude(planes[0], planes[1], magI);
					bscansub = magI.colRange(0, numdisplaypoints);
					bscansub.convertTo(bscansub, CV_64F);
					accumulate(bscansub,bscantemp);

					planes[0].copyTo(F);	
					Fsub = F.colRange(0, numdisplaypoints);
					Fsub.convertTo(Fsub, CV_64F);
					accumulate(Fsub,Fmreal);

					planes[1].copyTo(F);	
					Fsub = F.colRange(0, numdisplaypoints);
					Fsub.convertTo(Fsub, CV_64F);
					accumulate(Fsub,Fmimag);
					
					////////////// end of compute Bscan //////////////////////	
					bscantemp.copyTo(bscanlive);
					bscantemp.copyTo(X1);
					framecount ++;
					fps++;
				} // end of if ret == 1 block

				jdiff = bscanlive - bscanlivebg;
				jdiff.copyTo(positivediff);		// just to initialize the Mat
				makeonlypositive(jdiff, positivediff);
				positivediff += 0.000000001;			// to avoid log(0)
				log(positivediff, bscanlog);				// switch to logarithmic scale
				bscandb = 20.0 * bscanlog / 2.303;
				bscandb.row(4).copyTo(bscandb.row(1));	// masking out the DC in the display
				bscandb.row(4).copyTo(bscandb.row(0));
				tempmat1 = bscandb.colRange(0, numdisplaypoints);
			
				if(binvaluex > 1 || binvaluey > 1)
				{
					resize(tempmat1,tempmat2,Size(),binvaluex,binvaluey,INTER_AREA);
				}
				else
				{
					tempmat1.copyTo(tempmat2);
				}
				tempmat2.copyTo(bscandisp);
				bscandisp = max(bscandisp, bscanthreshold);
				normalize(bscandisp, bscandisp, 0, 1, NORM_MINMAX);	// normalize the log plot for display
				bscandisp.convertTo(bscandisp, CV_8UC1, 255.0);
				applyColorMap(bscandisp, cmagI, COLORMAP_JET);
				if (ROIflag == true)	
					rectangle(cmagI,ROItopleft,ROIbottright,Scalar(0,0,255),1, LINE_8);
				imshow("Live B-scan", cmagI);
				
				if(ii != 1)
				{
					// Intensity-based Differntiation method - Ref[22] in Li et.al. IbDiff
					absdiff(X1,X2,jdiff);
					// accumulate the absolute diff to Idiff
					accumulate(jdiff, Idiff1);
				}

				X1.copyTo(X2);

			} // end of for loop ii <= averagescount

			IDiff1 = Idiff1 / (1.0 * (averagescount-1));	
			IDiff1 += 0.000000001;			// to avoid log(0)
			log(IDiff1, bscanlog);				// switch to logarithmic scale
			bscandb1 = 20.0 * bscanlog / 2.303;
			bscandb1.row(4).copyTo(bscandb1.row(1));	// masking out the DC in the display
			bscandb1.row(4).copyTo(bscandb1.row(0));
			tempmat1 = bscandb1.colRange(0, numdisplaypoints);
			if(binvaluex > 1 || binvaluey > 1)
			{
				resize(tempmat1,tempmatdiff1,Size(),binvaluex,binvaluey,INTER_AREA);
			}
			else
			{
				tempmat1.copyTo(tempmatdiff1);
			}
			tempmatdiff1.copyTo(bscandisp1);
			bscandisp1 = max(bscandisp1, bscanthreshold);
			normalize(bscandisp1, bscandisp1, 0, 1, NORM_MINMAX);	// normalize the log plot for display
			bscandisp1.convertTo(bscandisp1, CV_8UC1, 255.0);
			applyColorMap(bscandisp1, cmagI1, COLORMAP_JET);
			if (ROIflag == true)	
				rectangle(cmagI1,ROItopleft,ROIbottright,Scalar(0,0,255),1, LINE_8);
			imshow("(1) IbDiff", cmagI1); // Intensity-based Differentiation

			gettimeofday(&tv,NULL);	
			time_end = 1000000 * tv.tv_sec + tv.tv_usec;	
			// update the image windows
			dt = time_end - time_start;
			if(dt > 500000) // 0.5 second in microseconds 
			{
				m.copyTo(mvector);
				mvector.reshape(0, 1);	//make it into a row array
				minMaxLoc(mvector, &minval, &maxval);
				sprintf(textbuffer, "fps = %d  Max val = %d", int(round(fps/dt*1e6)), int(floor(maxval)));
				firstrowofstatusimg = Mat::zeros(cv::Size(600, 50), CV_64F);
				putText(statusimg, textbuffer, Point(0, 30), FONT_HERSHEY_SIMPLEX, 1, Scalar(255, 255, 255), 3, 1);
				if(accummode == true)
					sprintf(textbuffer, "Accum mode");
				else if (averagescount == 1)	
					sprintf(textbuffer, "Live mode");
				else
					sprintf(textbuffer, "Avg mode N=%d",averagescount);

				secrowofstatusimgRHS = Mat::zeros(cv::Size(300, 50), CV_64F);
				putText(statusimg, textbuffer, Point(300, 80), FONT_HERSHEY_SIMPLEX, 1, Scalar(255, 255, 255), 3, 1);
				
				if(ROIflag == true)
				{
					// display max val of Bscan
					heightROI = ROIendrow - ROIstartrow;
					widthROI = ROIendcol - ROIstartcol;
					
					ROIcentrey = (ROIstartrow + ROIendrow) / 2;
					ROIcentrex = (ROIstartcol + ROIendcol) / 2;
					if(heightROI > 1 && widthROI > 1)
					{
						tempmatdiff1(Rect(ROIstartcol,ROIstartrow,widthROI,heightROI)).copyTo(ROI);
						ROI.reshape(0, 1);	//make it into a row array
						minMaxLoc(ROI, &minvaldiff, &maxvaldiff);
						meanvaldiff = mean(ROI)(0);
						Vmean = mean(ROI)(0);
					}
					sprintf(textbuffer,"6: mx=%.4e Av=%.4e",maxvaldiff,meanvaldiff);
				}
				else
				{
					sprintf(textbuffer, "Press i for ROI");
				}
				fourthrowofstatusimg = Mat::zeros(cv::Size(600, 50), CV_64F);
				putText(statusimg, textbuffer, Point(0,190), FONT_HERSHEY_SIMPLEX, 1, Scalar(255, 255, 255), 3, 1);
	

				resizeWindow("Status", 600, 300);
				imshow("Status", statusimg);
			
				fps = 0;
				gettimeofday(&tv,NULL);	
				time_start = 1000000 * tv.tv_sec + tv.tv_usec;	
			}		
			if(accummode == false)
			{
				Idiff1 = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
				Idiff2 = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
				Idiff5 = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
				num = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
				denom1 = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
				denom2 = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
				termReal = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
				termImag = Mat::zeros(Size(numdisplaypoints, oph), CV_64F);
				framecount = 0.0;
			}

			if (skeypressed == true)
			{
				if(dir_created == false)
				{
					// create a directory with time stamp
					struct tm *timenow;
					time_t now = time(NULL);
					timenow = localtime(&now);
					strftime(dirname, sizeof(dirname), "%Y-%m-%d_%H_%M_%S-", timenow);
					strcat(dirname, dirdescr);
					mkdir(dirname, 0755);
					strcpy(pathname, dirname);
					strcat(pathname, "/");
					dir_created = true;
				}
				indexi++;
				sprintf(filename, "bscan%03d", indexi);
				savematasbin(pathname, dirname, filename, bscandb);
				savematasimage(pathname, dirname, filename, bscandisp);
				sprintf(filenamec, "bscanc%03d", indexi);
				savematasimage(pathname, dirname, filenamec, cmagI);

				sprintf(filename, "IbDiff%03d", indexi);
				savematasbin(pathname, dirname, filename, bscandb1);
				savematasimage(pathname, dirname, filename, bscandisp1);
				sprintf(filenamec, "IbDiffc%03d", indexi);
				savematasimage(pathname, dirname, filenamec, cmagI1);

				sprintf(filename, "status%03d", indexi);
				savematasimage(pathname, dirname, filename, statusimg);
				skeypressed = false;
			}			

			key = waitKey(3); // wait for keypress
			switch (key)
			{

			case 27: //ESC key
			case 'x':
			case 'X':
				doneflag = true;
				break;
			
			case '+':
				camtime = camtime + 1;
				expchanged = true;
				break;

			case '=':
				camtime = camtime + 100;
				expchanged = true;
				break;

			case '-':
				if (camtime < 8)	// spinnaker has a min of 8 microsec
				{
					camtime = 8;
					break;
				}
				camtime = camtime - 100;
				expchanged = true;
				break;
		
			case '_':
				if (camtime < 8)	// spinnaker has a min of 8 microsec
				{
					camtime = 8;
					break;
				}
				camtime = camtime - 1;
				expchanged = true;
				break;
		
			case ']':
				bscanthreshold += pow(10,threshold_N);
				sprintf(textbuffer, "Threshold = %5.3f", bscanthreshold);
				thirdrowofstatusimg = Mat::zeros(cv::Size(600, 50), CV_64F);
				putText(statusimg, textbuffer, Point(0, 130), FONT_HERSHEY_SIMPLEX, 1, Scalar(255, 255, 255), 3, 1);
				imshow("Status", statusimg);
				break;

			case '[':
				bscanthreshold -= pow(10,threshold_N);   
				sprintf(textbuffer, "Threshold = %5.3f", bscanthreshold);
				thirdrowofstatusimg = Mat::zeros(cv::Size(600, 50), CV_64F);
				putText(statusimg, textbuffer, Point(0, 130), FONT_HERSHEY_SIMPLEX, 1, Scalar(255, 255, 255), 3, 1);
				imshow("Status", statusimg);
				break;
	
			case ',':
				threshold_N --;
				break;			
	
			case '.':
				threshold_N ++;
				break;			
	
			case '4':
				if (accummode == true)
					break;
				if (averagescount == 1)
					averagescount = 10; 
				else 
					averagescount += 10;
				bgframestaken = false;
				break;

			case '3':
				if (accummode == true)
					break;
				// decrement number of averages by 5
				if(averagescount >= 11)
					averagescount -= 10;
				else if(averagescount <= 10)
					averagescount = 1;
				bgframestaken = false;
				break;

			case 's':			
				// save vibrational profile of bscan			
				skeypressed = true;
				break;

			case 'b':
				// new bscan background			
				bgframestaken = false;
				break;

			case 'i':
				// toggle ROI display
				ROIflag = !ROIflag;
				break;
		
			case 'a':
				// accumulate mode
				accummode = true;
				//averagescount = 1;
				//bgframestaken = false;
				break;

			case 'l':
				// live mode
				accummode = false;
				averagescount = 1;
				//bgframestaken = false;
				break;
				
			case ' ': // spacebar pressed
				spacebarpressed = ! spacebarpressed;
				break;
	
				// keys to change ROI
			case 't':
				if (ROIflag == false)
					break;
				// select the ROI option to change position
				selectedROIoption = 1; // position
				break; 			 

			case 'r':
				if (ROIflag == false)
					break;
				// select the ROI option to change size
				selectedROIoption = 2; // size
				break; 
			 
			case 82: 							// up arrow key = R ?
				if (ROIflag == false)
					break;
				if(selectedROIoption == position)
				{
					// move ROI up
					if(ROIstartrow > 5)
					{
						ROIstartrow -= 1;
						ROIendrow -= 1;
					}
				}
				if(selectedROIoption == size)
				{
					// blow up the ROI vertically
					if(ROIstartrow > 5) 
						ROIstartrow -= 1;
					if(ROIendrow < (oph*binvaluey-5))
						ROIendrow += 1;
				}
				ROItopleft.y = ROIstartrow;
				ROIbottright.y = ROIendrow;
				selectedPixel.y = (ROIstartrow+ROIendrow)/2; 
				break;

			case 84: 						// down arrow key = T ?
				if (ROIflag == false)
					break;
				if(selectedROIoption == position)
				{
					// move ROI down
					if(ROIendrow < (oph*binvaluey-5))
					{
						ROIstartrow += 1;
						ROIendrow += 1;
					}
				}
				if(selectedROIoption == size)
				{
					// shrink the ROI vertically
					if(ROIendrow-ROIstartrow <= 1)
					{
						ROIstartrow -= 1;
						ROIendrow += 1;
					}
					else
					{
						ROIstartrow += 1;
						ROIendrow -= 1;
					}
				}
				ROItopleft.y = ROIstartrow;
				ROIbottright.y = ROIendrow;
				selectedPixel.y = (ROIstartrow+ROIendrow)/2; 
				break;

			case 81:						// left arrow key = Q ?
				if (ROIflag == false)
					break;
				if(selectedROIoption == position)
				{ 
					// move ROI left
					if(ROIstartcol > 5)
					{
						ROIstartcol -= 1;
						ROIendcol -= 1;
					}
				}
				if(selectedROIoption == size)
				{
					// blow up the ROI horizontally
					if(ROIstartcol > 5) 
						ROIstartcol -= 1;
					if(ROIendcol < (numdisplaypoints*binvaluey-5))
						ROIendcol += 1;
				}
				ROItopleft.x = ROIstartcol;
				ROIbottright.x = ROIendcol;
				selectedPixel.x = (ROIstartcol+ROIendcol)/2; 
				break;

			case 83:						// right arrow key = S ?
				if (ROIflag == false)
					break;
				if(selectedROIoption == position)
				{ 
					// move ROI right
					if(ROIendcol < (numdisplaypoints*binvaluex-5))
					{
						ROIstartcol += 1;
						ROIendcol += 1;
					}
				}
				if(selectedROIoption == size)
				{
					// shrink the ROI horizontally
					if(ROIendcol-ROIstartcol <= 1)
					{
						ROIstartcol -= 1;
						ROIendcol += 1;
					}
					else
					{
						ROIstartcol += 1;
						ROIendcol -= 1;
					}
				}
				ROItopleft.x = ROIstartcol;
				ROIbottright.x = ROIendcol;
				selectedPixel.x = (ROIstartcol+ROIendcol)/2; 
				break;
			

			default:
				break;
			}

			if (doneflag == 1)
			{
				break;
			}
			if (expchanged == true)
			{
				//Set exp with QuickSpin
				ret = 0;
				if (IsReadable(pCam->ExposureTime) && IsWritable(pCam->ExposureTime))
				{
					//pCam->ExposureTime.SetValue(camtime);
					ret = 1;
				}
				if (ret == 1)
				{
					sprintf(textbuffer, "Exp time hard-coded ");
					secrowofstatusimg = Mat::zeros(cv::Size(600, 50), CV_64F);
					putText(statusimg, textbuffer, Point(0, 80), FONT_HERSHEY_SIMPLEX, 1, Scalar(255, 255, 255), 3, 1);
					imshow("Status", statusimg);

				}
				else
				{
					sprintf(textbuffer, "CONTROL_EXPOSURE failed");
					secrowofstatusimg = Mat::zeros(cv::Size(600, 50), CV_64F);
					putText(statusimg, textbuffer, Point(0, 80), FONT_HERSHEY_SIMPLEX, 1, Scalar(255, 255, 255), 3, 1);
					imshow("Status", statusimg);
				}

			} // end of if expchanged

		}
		pCam->EndAcquisition();
		pCam->DeInit();
		pCam = nullptr;

		// Clear camera list before releasing system
    	camList.Clear();
    	
		// Release system
    	system->ReleaseInstance();
	}
	catch (Spinnaker::Exception &e)
	{
		cout << "Error: " << e.what() << endl;
		result = -1;
	}

    return result;
}

// Function definitions
void setCamera(CameraPtr pCam)
{
	int result = 0;    
	unsigned int w, h, camspeed, burstframecount,triggerdelay, camtime, camgain = 1, bpp;
	unsigned int offsetx = 0, offsety = 0;
	unsigned int cambinx, cambiny;
	
	ifstream infile("spint.ini");
	string tempstring;
	
	// inputs from ini file
	if (infile.is_open())
	{
		infile >> tempstring;
		infile >> tempstring;
		infile >> tempstring;
		// first three lines of ini file are comments
		infile >> camgain;
		infile >> tempstring;
		infile >> camtime;
		infile >> tempstring;
		infile >> bpp;
		infile >> tempstring;
		infile >> w;
		infile >> tempstring;
		infile >> h;
		infile >> tempstring;
		infile >> offsetx;
		infile >> tempstring;
		infile >> offsety;
		infile >> tempstring;
		infile >> camspeed;
		infile >> tempstring;
		infile >> cambinx;
		infile >> tempstring;
		infile >> cambiny;

		infile.close();
	}

	cout << "Initialising Camera settings ..." << endl;
	
	pCam->TLStream.StreamBufferHandlingMode.SetValue(StreamBufferHandlingMode_NewestOnly);
	pCam->AcquisitionMode.SetValue(AcquisitionMode_Continuous);
		
	// gain
	pCam->GainAuto.SetValue(GainAuto_Off);	
	pCam->Gain.SetValue(camgain);
	cout << "Gain set to " << pCam->Gain.GetValue() << " dB ..." << endl;

	// exposure time
	pCam->ExposureAuto.SetValue(ExposureAuto_Off);
	pCam->ExposureMode.SetValue(ExposureMode_Timed);
	pCam->ExposureTime.SetValue(camtime);
	cout << "Exp set to " << pCam->ExposureTime.GetValue() << " microsec ..." << endl;

	// bpp or cambitdepth 
	if (bpp == 16)
	{
		pCam->PixelFormat.SetValue(PixelFormat_Mono16);
		cout << "Pixel format set to " << pCam->PixelFormat.GetCurrentEntry()->GetSymbolic() << "..." << endl;
	}
	
	// cambinx
	pCam->BinningHorizontal.SetValue(cambinx);
	cout << "BinningHorizontal set to " << pCam->BinningHorizontal.GetValue() << "..." << endl;

	// cambiny
	pCam->BinningVertical.SetValue(cambiny);
	cout << "BinningVertical set to " << pCam->BinningVertical.GetValue() << "..." << endl;
	
	// width 
	if (IsReadable(pCam->Width) && IsWritable(pCam->Width))
	{
		pCam->Width.SetValue(w);
	}
	else
	{
		cout << "Width not available..." << endl;
	}
	
	// height 
	if (IsReadable(pCam->Height) && IsWritable(pCam->Height))
	{
		pCam->Height.SetValue(h);
	}
	else
	{
		cout << "Height not available..." << endl;
	}

	// offsetx
	if (IsReadable(pCam->OffsetX) && IsWritable(pCam->OffsetX))
	{
		pCam->OffsetX.SetValue(offsetx);
	}
	else
	{
		cout << "Offset X not available..." << endl;
	}
	
	// offsety
	if (IsReadable(pCam->OffsetY) && IsWritable(pCam->OffsetY))
	{
		pCam->OffsetY.SetValue(offsety);
	}
	else
	{
		cout << "Offset Y not available..." << endl;
	}

	// frame rate
	pCam->AcquisitionFrameRateEnable.SetValue(1);
	pCam->AcquisitionFrameRate.SetValue(camspeed);
	cout << "Frame rate set to " << camspeed << endl;

	// set the hardware trigger	     
	pCam->TriggerMode.SetValue(TriggerMode_Off);
	cout << "Camera set to trigger mode OFF "<< endl;
}

inline void makeonlypositive(Mat& src, Mat& dst)
{
	// from https://stackoverflow.com/questions/48313249/opencv-convert-all-negative-values-to-zero
	max(src, 0, dst);

}

inline void savematasbin(char* p, char* d, char* f, Mat m)
{
	// saves a Mat m by writing to a binary file  f appending .ocv, both windows and unix versions
	// p=pathname, d=dirname, f=filename

#ifdef __unix__
	strcpy(p, d);
	strcat(p, "/");
	strcat(p, f);
	strcat(p, ".ocv");
	matwrite(p, m);
#else

	strcpy(p, d);
	strcat(p, "\\");		// imwrite needs path with \\ separators, not /, on windows
	strcat(p, f);
	strcat(p, ".ocv");
	matwrite(p, m);
#endif	

}

inline void savematasimage(char* p, char* d, char* f, Mat m)
{
	// saves a Mat m using imwrite as filename f appending .png, both windows and unix versions
	// p=pathname, d=dirname, f=filename

#ifdef __unix__
	strcpy(p, d);
	strcat(p, "/");
	strcat(p, f);
	strcat(p, ".png");
	imwrite(p, m);
#else

	strcpy(p, d);
	strcat(p, "\\");		// imwrite needs path with \\ separators, not /, on windows
	strcat(p, f);
	strcat(p, ".png");
	imwrite(p, m);
#endif	

}

// from http://stackoverflow.com/a/32357875/5294258
void matwrite(const string& filename, const Mat& mat)
{
	ofstream fs(filename, fstream::binary);

	// Header
	int type = mat.type();
	int channels = mat.channels();
	fs.write((char*)&mat.rows, sizeof(int));    // rows
	fs.write((char*)&mat.cols, sizeof(int));    // cols
	fs.write((char*)&type, sizeof(int));        // type
	fs.write((char*)&channels, sizeof(int));    // channels

												// Data
	if (mat.isContinuous())
	{
		fs.write(mat.ptr<char>(0), (mat.dataend - mat.datastart));
	}
	else
	{
		int rowsz = CV_ELEM_SIZE(type) * mat.cols;
		for (int r = 0; r < mat.rows; ++r)
		{
			fs.write(mat.ptr<char>(r), rowsz);
		}
	}
}
