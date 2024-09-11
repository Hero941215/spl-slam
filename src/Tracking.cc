/**
* This file is part of ORB-SLAM2.
*
* Copyright (C) 2014-2016 Raúl Mur-Artal <raulmur at unizar dot es> (University of Zaragoza)
* For more information see <https://github.com/raulmur/ORB_SLAM2>
*
* ORB-SLAM2 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* ORB-SLAM2 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with ORB-SLAM2. If not, see <http://www.gnu.org/licenses/>.
*/


#include "Tracking.h"

#include<opencv2/core/core.hpp>
#include<opencv2/features2d/features2d.hpp>

#include"ORBmatcher.h"
#include"Linematcher.h" //(LineExpanding)
#include"FrameDrawer.h"
#include"Converter.h"
#include"Map.h"
#include"Initializer.h"

#include"Optimizer.h"
#include"PnPsolver.h"

#include<iostream>

#include<mutex>

#include"Timer.h"//直接从PL-SLAM中抄过来的一个时间计算器，内核是boost::chrono


using namespace std;

namespace PL_SLAM
{

Tracking::Tracking(System *pSys, ORBVocabulary* pVoc, FrameDrawer *pFrameDrawer, MapDrawer *pMapDrawer, Map *pMap, 
		   KeyFrameDatabase* pKFDB, const string &strSettingPath, const int sensor):
    mState(NO_IMAGES_YET), mSensor(sensor), mbOnlyTracking(false), mbVO(false), mpORBVocabulary(pVoc),
    mpKeyFrameDB(pKFDB), mpInitializer(static_cast<Initializer*>(NULL)), mpSystem(pSys), mpViewer(NULL),
    mpFrameDrawer(pFrameDrawer), mpMapDrawer(pMapDrawer), mpMap(pMap), mnLastRelocFrameId(0)
{
    // Load camera parameters from settings file

    cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);
    float fx = fSettings["Camera.fx"];
    float fy = fSettings["Camera.fy"];
    float cx = fSettings["Camera.cx"];
    float cy = fSettings["Camera.cy"];

    cv::Mat K = cv::Mat::eye(3,3,CV_32F);
    K.at<float>(0,0) = fx;
    K.at<float>(1,1) = fy;
    K.at<float>(0,2) = cx;
    K.at<float>(1,2) = cy;
    K.copyTo(mK);

    cv::Mat DistCoef(4,1,CV_32F);
    DistCoef.at<float>(0) = fSettings["Camera.k1"];
    DistCoef.at<float>(1) = fSettings["Camera.k2"];
    DistCoef.at<float>(2) = fSettings["Camera.p1"];
    DistCoef.at<float>(3) = fSettings["Camera.p2"];
    const float k3 = fSettings["Camera.k3"];
    if(k3!=0)
    {
        DistCoef.resize(5);
        DistCoef.at<float>(4) = k3;
    }
    DistCoef.copyTo(mDistCoef);

    mbf = fSettings["Camera.bf"];

    float fps = fSettings["Camera.fps"];
    if(fps==0)
        fps=30;

    // Max/Min Frames to insert keyframes and to check relocalisation
    mMinFrames = 0;
    mMaxFrames = fps;

    cout << endl << "Camera Parameters: " << endl;
    cout << "- fx: " << fx << endl;
    cout << "- fy: " << fy << endl;
    cout << "- cx: " << cx << endl;
    cout << "- cy: " << cy << endl;
    cout << "- k1: " << DistCoef.at<float>(0) << endl;
    cout << "- k2: " << DistCoef.at<float>(1) << endl;
    if(DistCoef.rows==5)
        cout << "- k3: " << DistCoef.at<float>(4) << endl;
    cout << "- p1: " << DistCoef.at<float>(2) << endl;
    cout << "- p2: " << DistCoef.at<float>(3) << endl;
    cout << "- fps: " << fps << endl;


    int nRGB = fSettings["Camera.RGB"];
    mbRGB = nRGB;

    if(mbRGB)
        cout << "- color order: RGB (ignored if grayscale)" << endl;
    else
        cout << "- color order: BGR (ignored if grayscale)" << endl;

    // Load ORB parameters
    int nFeatures = fSettings["ORBextractor.nFeatures"];
    float fScaleFactor = fSettings["ORBextractor.scaleFactor"];
    int nLevels = fSettings["ORBextractor.nLevels"];
    int fIniThFAST = fSettings["ORBextractor.iniThFAST"];
    int fMinThFAST = fSettings["ORBextractor.minThFAST"];

    mpORBextractorLeft = new ORBextractor(nFeatures,fScaleFactor,nLevels,fIniThFAST,fMinThFAST);

    if(sensor==System::STEREO)
        mpORBextractorRight = new ORBextractor(nFeatures,fScaleFactor,nLevels,fIniThFAST,fMinThFAST);

    if(sensor==System::MONOCULAR)
        mpIniORBextractor = new ORBextractor(2*nFeatures,fScaleFactor,nLevels,fIniThFAST,fMinThFAST);

    cout << endl  << "ORB Extractor Parameters: " << endl;
    cout << "- Number of Features: " << nFeatures << endl;
    cout << "- Scale Levels: " << nLevels << endl;
    cout << "- Scale Factor: " << fScaleFactor << endl;
    cout << "- Initial Fast Threshold: " << fIniThFAST << endl;
    cout << "- Minimum Fast Threshold: " << fMinThFAST << endl;
    
    mpORBmatcherLeftIni = new ORBmatcher(0.9,true);
    mpORBmatcherLeftmotion = new ORBmatcher(0.9,true);
    mpORBmatcherLeftBoW = new ORBmatcher(0.75,true);
    mpORBmatcherLeftRel = new ORBmatcher(0.9,true);
    
    //从配置文件中获取系统跟踪是否使用线特征
    int nUsingLine = fSettings["System.usingLine"];
    mbusingLine = nUsingLine;
    
    cout << endl << "mbusingLine: " << mbusingLine << endl;
    
    //若使用线特征，则判断使用LSD还是FLD作为关键线提取器
    if(mbusingLine)
    {
	int nUsingLsdFeature = fSettings["System.usingLsdFeature"];
	mbusingLsdFeature = nUsingLsdFeature;
	if(mbusingLsdFeature)
	    cout << "使用LSD线提取算法!" << endl;
	else
	    cout << "使用FLD线提取算法!" << endl;
    }
    
    if(sensor==System::MONOCULAR && mbusingLine == true)
    {
	if(mbusingLsdFeature)
	{
	    //从配置文件中获取图像的宽、高
	    int camerawidth = fSettings["Camera.width"];
	    int cameraheight = fSettings["Camera.height"];
	  
	    //从配置文件中获取线特征最小允许比例
	    float min_line_length_ratio = fSettings["Lineextractor.min_line_length_ratio"];
	    float min_line_length = min_line_length_ratio * std::min( camerawidth, cameraheight );
	  
	    //从配置文件中获取其他线特征提取所需参数
	    int nfeatures = fSettings["Lineextractor.nFeatures"];
	    int nlevels = fSettings["Lineextractor.nLevels"];
	    int refine = fSettings["Lineextractor.refine"];
	    double scale = fSettings["Lineextractor.scale"];
	    double sigma_scale = fSettings["Lineextractor.sigma_scale"]; 
	    double quant = fSettings["Lineextractor.quant"];
	    double ang_th = fSettings["Lineextractor.ang_th"];
	    double log_eps = fSettings["Lineextractor.log_eps"];
	    double density_th = fSettings["Lineextractor.density_th"];
	    int n_bins = fSettings["Lineextractor.n_bins"];
	    

	    mpIniLineextractor = new Lineextractor(2*nfeatures, nlevels, refine, scale, sigma_scale, 
						  quant, ang_th, log_eps, density_th, n_bins, min_line_length, mbusingLsdFeature);

	    mpLineextractorLeft = new Lineextractor(nfeatures, nlevels, refine, scale, sigma_scale, 
						  quant, ang_th, log_eps, density_th, n_bins, min_line_length, mbusingLsdFeature);
	
	    //输出线特征提取相关参数
	    cout << endl  << "LSD Extractor Parameters: " << endl;
	    cout << "- Lineextractor.nFeatures: " << nfeatures << endl;
	    cout << "- Lineextractor.nLevels: " << nlevels << endl;
	    cout << "- Lineextractor.refine: " << refine << endl;
	    cout << "- Lineextractor.scale: " << scale << endl;
	    cout << "- Lineextractor.sigma_scale: " << sigma_scale << endl;
	    cout << "- Lineextractor.quant: " << quant << endl;
	    cout << "- Lineextractor.ang_th: " << ang_th << endl;
	    cout << "- Lineextractor.log_eps: " << log_eps << endl;
	    cout << "- Lineextractor.density_th: " << density_th << endl;
	    cout << "- Lineextractor.n_bins: " << n_bins << endl;
	    cout << "- Lineextractor.min_line_length_ratio: " << min_line_length_ratio << endl << endl;
	    
	    //创建线特征匹配器指针（4组）
	    mpLinematcherLeftIni = new Linematcher(0.9,true,true,0.15);
	    mpLinematcherLeftmotion = new Linematcher(0.9,true,true,0.2);//第二个参数最好是true
	    mpLinematcherLeftKNN = new Linematcher(0.75,false,true,0.25);//同时在重定位和参考关键帧跟踪中使用,第二个参数必须是true
	    mpLinematcherLeftRel = new Linematcher(0.9,true,false,0.15);
	}
	else
	{
	  
	    //从配置文件中获取FLD线提取所需参数
	    int nfeatures = fSettings["Lineextractor.nFeatures"];
	    int nlevels = fSettings["Lineextractor.nLevels"];
	    double scale = fSettings["Lineextractor.scale"]; 
	    int threshold_length = fSettings["Lineextractor.threshold_length"];
	    float threshold_dist = fSettings["Lineextractor.threshold_dist"];
	    double canny_th1 = fSettings["Lineextractor.canny_th1"];
	    double canny_th2 = fSettings["Lineextractor.canny_th2"];
	    int canny_aperture_size = fSettings["Lineextractor.canny_aperture_size"];
	    int n_do_merge = fSettings["Lineextractor.do_merge"];
	    bool do_merge = n_do_merge;

	    mpIniLineextractor = new Lineextractor(2*nfeatures, nlevels, scale , threshold_length, threshold_dist, canny_th1, 
						  canny_th2, canny_aperture_size, do_merge);

	    mpLineextractorLeft = new Lineextractor(nfeatures, nlevels, scale, threshold_length, threshold_dist, canny_th1, 
						  canny_th2, canny_aperture_size, do_merge);
	
	    //输出线特征提取相关参数
	    cout << endl  << "FLD Extractor Parameters: " << endl;
	    cout << "- Lineextractor.nFeatures: " << nfeatures << endl;
	    cout << "- Lineextractor.nLevels: " << nlevels << endl;
	    cout << "- Lineextractor.threshold_length: " << threshold_length << endl;
	    cout << "- Lineextractor.threshold_dist: " << threshold_dist << endl;
	    cout << "- Lineextractor.canny_th1: " << canny_th1 << endl;
	    cout << "- Lineextractor.canny_th2: " << canny_th2 << endl;
	    cout << "- Lineextractor.canny_aperture_size: " << canny_aperture_size << endl;
	    cout << "- Lineextractor.do_merge: " << do_merge << endl << endl;
	    
	    
	    //创建线特征匹配器指针（4组）
	    mpLinematcherLeftIni = new Linematcher(0.9,true,true,0.15);
	    mpLinematcherLeftmotion = new Linematcher(0.9,true,true,0.2);
	    mpLinematcherLeftKNN = new Linematcher(0.75,false,true,0.25);
	    mpLinematcherLeftRel = new Linematcher(0.9,true,false,0.15);
	}     
      
    }
     
    if(sensor==System::STEREO || sensor==System::RGBD)
    {
        mThDepth = mbf*(float)fSettings["ThDepth"]/fx;
        cout << endl << "Depth Threshold (Close/Far Points): " << mThDepth << endl;
    }

    if(sensor==System::RGBD)
    {
        mDepthMapFactor = fSettings["DepthMapFactor"];
        if(fabs(mDepthMapFactor)<1e-5)
            mDepthMapFactor=1;
        else
            mDepthMapFactor = 1.0f/mDepthMapFactor;
    }
    
    TrackingAllTime = 0.0;
    TrackingAverageTime = 0.0;
    TrackingNumbers = 0;

}

void Tracking::SetLocalMapper(LocalMapping *pLocalMapper)
{
    mpLocalMapper=pLocalMapper;
}

void Tracking::SetLoopClosing(LoopClosing *pLoopClosing)
{
    mpLoopClosing=pLoopClosing;
}

void Tracking::SetViewer(Viewer *pViewer)
{
    mpViewer=pViewer;
}


cv::Mat Tracking::GrabImageStereo(const cv::Mat &imRectLeft, const cv::Mat &imRectRight, const double &timestamp)
{
    mImGray = imRectLeft;
    cv::Mat imGrayRight = imRectRight;

    if(mImGray.channels()==3)
    {
        if(mbRGB)
        {
            cvtColor(mImGray,mImGray,CV_RGB2GRAY);
            cvtColor(imGrayRight,imGrayRight,CV_RGB2GRAY);
        }
        else
        {
            cvtColor(mImGray,mImGray,CV_BGR2GRAY);
            cvtColor(imGrayRight,imGrayRight,CV_BGR2GRAY);
        }
    }
    else if(mImGray.channels()==4)
    {
        if(mbRGB)
        {
            cvtColor(mImGray,mImGray,CV_RGBA2GRAY);
            cvtColor(imGrayRight,imGrayRight,CV_RGBA2GRAY);
        }
        else
        {
            cvtColor(mImGray,mImGray,CV_BGRA2GRAY);
            cvtColor(imGrayRight,imGrayRight,CV_BGRA2GRAY);
        }
    }

    mCurrentFrame = Frame(mImGray,imGrayRight,timestamp,mpORBextractorLeft,mpORBextractorRight,mpORBVocabulary,mK,mDistCoef,mbf,mThDepth);

    Track();

    return mCurrentFrame.mTcw.clone();
}


cv::Mat Tracking::GrabImageRGBD(const cv::Mat &imRGB,const cv::Mat &imD, const double &timestamp)
{
    mImGray = imRGB;
    cv::Mat imDepth = imD;

    if(mImGray.channels()==3)
    {
        if(mbRGB)
            cvtColor(mImGray,mImGray,CV_RGB2GRAY);
        else
            cvtColor(mImGray,mImGray,CV_BGR2GRAY);
    }
    else if(mImGray.channels()==4)
    {
        if(mbRGB)
            cvtColor(mImGray,mImGray,CV_RGBA2GRAY);
        else
            cvtColor(mImGray,mImGray,CV_BGRA2GRAY);
    }

    if((fabs(mDepthMapFactor-1.0f)>1e-5) || imDepth.type()!=CV_32F)
        imDepth.convertTo(imDepth,CV_32F,mDepthMapFactor);

    mCurrentFrame = Frame(mImGray,imDepth,timestamp,mpORBextractorLeft,mpORBVocabulary,mK,mDistCoef,mbf,mThDepth);

    Track();

    return mCurrentFrame.mTcw.clone();
}

//修改为适合点和点线特征的单目跟踪器(根据数据成员)
cv::Mat Tracking::GrabImageMonocular(const cv::Mat &im, const double &timestamp)
{
    Timer timer;
    
    mImGray = im;
    
    if(mImGray.channels()==3)
    {
        if(mbRGB)
            cvtColor(mImGray,mImGray,CV_RGB2GRAY);
        else
            cvtColor(mImGray,mImGray,CV_BGR2GRAY);
    }
    else if(mImGray.channels()==4)
    {
        if(mbRGB)
            cvtColor(mImGray,mImGray,CV_RGBA2GRAY);
        else
            cvtColor(mImGray,mImGray,CV_BGRA2GRAY);
    }
    
    timer.start();
    
    if(mState==NOT_INITIALIZED || mState==NO_IMAGES_YET)
	if(mbusingLine)	    
	    mCurrentFrame = Frame(mImGray,timestamp,mpIniORBextractor,mpORBVocabulary,mpIniLineextractor,mK,mDistCoef,mbf,mThDepth);
	else
	    mCurrentFrame = Frame(mImGray,timestamp,mpIniORBextractor,mpORBVocabulary,mK,mDistCoef,mbf,mThDepth);
    else
	if(mbusingLine) 
	    mCurrentFrame = Frame(mImGray,timestamp,mpORBextractorLeft,mpORBVocabulary,mpLineextractorLeft,mK,mDistCoef,mbf,mThDepth);
	else
	    mCurrentFrame = Frame(mImGray,timestamp,mpORBextractorLeft,mpORBVocabulary,mK,mDistCoef,mbf,mThDepth);
	
    double FrameBuildTime = timer.stop();
    
    cout << "Frame创建时间: " << FrameBuildTime << " ms" <<endl;
    
    timer.start();
   
    if(mbusingLine)
	TrackBoth();
    else  
	Track();
    
    double TrackingTime = timer.stop();
    
    cout << "单步跟踪时间: " << TrackingTime << " ms" <<endl <<endl ;
    
    TrackingAllTime += TrackingTime;
    TrackingNumbers++;
    TrackingAverageTime = TrackingAllTime / TrackingNumbers;

    cout << "平均跟踪时间: " << TrackingAverageTime << " ms" <<endl <<endl ;
    
    return mCurrentFrame.mTcw.clone();
}

//仅使用点特征
void Tracking::Track()
{
    if(mState==NO_IMAGES_YET)
    {
        mState = NOT_INITIALIZED;
    }

    mLastProcessedState=mState;

    // Get Map Mutex -> Map cannot be changed
    unique_lock<mutex> lock(mpMap->mMutexMapUpdate);

    if(mState==NOT_INITIALIZED)
    {
        if(mSensor==System::STEREO || mSensor==System::RGBD)
            StereoInitialization();
        else
            MonocularInitialization();

        mpFrameDrawer->Update(this);

        if(mState!=OK)
            return;
    }
    else
    {
        // System is initialized. Track Frame.
        bool bOK;

        // Initial camera pose estimation using motion model or relocalization (if tracking is lost)
        if(!mbOnlyTracking)
        {
            // Local Mapping is activated. This is the normal behaviour, unless
            // you explicitly activate the "only tracking" mode.

            if(mState==OK)
            {
                // Local Mapping might have changed some MapPoints tracked in last frame
                CheckReplacedInLastFrame();

                if(mVelocity.empty() || mCurrentFrame.mnId<mnLastRelocFrameId+2)
                {
                    bOK = TrackReferenceKeyFrame();
		    cout << "初始化刚完成或者刚重定位成功,执行参考关键帧跟踪成功： " << true << endl <<endl;
                }
                else
                {
                    bOK = TrackWithMotionModel();
		    cout << "运动模型跟踪成功： " << true << endl <<endl;	    
                    if(!bOK)
		    {
			cout << "运动模型失败,执行参考关键帧跟踪： " << true << endl <<endl;
			bOK = TrackReferenceKeyFrame();
			if(OK)
			{
			    cout << "运动模型失败,执行参考关键帧跟踪成功： " << true << endl <<endl;
			}
			else
			{
			    cout << "运动模型失败,执行参考关键帧跟踪也失败,跟踪丢失： " << true << endl <<endl;
			}
			
		    }
                        
                }
            }
            else
            {
		cout << "执行重定位: " << true << endl << endl;
                bOK = Relocalization();
		if(bOK)
		{
		    cout << "重定位成功: " << true << endl << endl;
		}
		else
		{
		    cout << "重定位失败: " << true << endl << endl;
		}
		
            }
        }
        else
        {
            // Localization Mode: Local Mapping is deactivated

            if(mState==LOST)
            {
                bOK = Relocalization();
            }
            else
            {
                if(!mbVO)
                {
                    // In last frame we tracked enough MapPoints in the map

                    if(!mVelocity.empty())
                    {
                        bOK = TrackWithMotionModel();
                    }
                    else
                    {
                        bOK = TrackReferenceKeyFrame();
                    }
                }
                else
                {
                    // In last frame we tracked mainly "visual odometry" points.

                    // We compute two camera poses, one from motion model and one doing relocalization.
                    // If relocalization is sucessfull we choose that solution, otherwise we retain
                    // the "visual odometry" solution.

                    bool bOKMM = false;
                    bool bOKReloc = false;
                    vector<MapPoint*> vpMPsMM;
                    vector<bool> vbOutMM;
                    cv::Mat TcwMM;
                    if(!mVelocity.empty())
                    {
                        bOKMM = TrackWithMotionModel();
                        vpMPsMM = mCurrentFrame.mvpMapPoints;
                        vbOutMM = mCurrentFrame.mvbOutlier;
                        TcwMM = mCurrentFrame.mTcw.clone();
                    }
                    bOKReloc = Relocalization();

                    if(bOKMM && !bOKReloc)
                    {
                        mCurrentFrame.SetPose(TcwMM);
                        mCurrentFrame.mvpMapPoints = vpMPsMM;
                        mCurrentFrame.mvbOutlier = vbOutMM;

                        if(mbVO)
                        {
                            for(int i =0; i<mCurrentFrame.N; i++)
                            {
                                if(mCurrentFrame.mvpMapPoints[i] && !mCurrentFrame.mvbOutlier[i])
                                {
                                    mCurrentFrame.mvpMapPoints[i]->IncreaseFound();
                                }
                            }
                        }
                    }
                    else if(bOKReloc)
                    {
                        mbVO = false;
                    }

                    bOK = bOKReloc || bOKMM;
                }
            }
        }

        mCurrentFrame.mpReferenceKF = mpReferenceKF;

        // If we have an initial estimation of the camera pose and matching. Track the local map.
        if(!mbOnlyTracking)
        {
            if(bOK)
                bOK = TrackLocalMap();
        }
        else
        {
            // mbVO true means that there are few matches to MapPoints in the map. We cannot retrieve
            // a local map and therefore we do not perform TrackLocalMap(). Once the system relocalizes
            // the camera we will use the local map again.
            if(bOK && !mbVO)
                bOK = TrackLocalMap();
        }

        if(bOK)
            mState = OK;
        else
            mState=LOST;

        // Update drawer
        mpFrameDrawer->Update(this);

        // If tracking were good, check if we insert a keyframe
        if(bOK)
        {
            // Update motion model
            if(!mLastFrame.mTcw.empty())
            {
                cv::Mat LastTwc = cv::Mat::eye(4,4,CV_32F);
                mLastFrame.GetRotationInverse().copyTo(LastTwc.rowRange(0,3).colRange(0,3));
                mLastFrame.GetCameraCenter().copyTo(LastTwc.rowRange(0,3).col(3));
                mVelocity = mCurrentFrame.mTcw*LastTwc;
            }
            else
                mVelocity = cv::Mat();

            mpMapDrawer->SetCurrentCameraPose(mCurrentFrame.mTcw);

            // Clean VO matches
            for(int i=0; i<mCurrentFrame.N; i++)
            {
                MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];
                if(pMP)
                    if(pMP->Observations()<1)
                    {
                        mCurrentFrame.mvbOutlier[i] = false;
                        mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
                    }
            }

            // Delete temporal MapPoints
            for(list<MapPoint*>::iterator lit = mlpTemporalPoints.begin(), lend =  mlpTemporalPoints.end(); lit!=lend; lit++)
            {
                MapPoint* pMP = *lit;
                delete pMP;
            }
            mlpTemporalPoints.clear();

            // Check if we need to insert a new keyframe
            if(NeedNewKeyFrame())
                CreateNewKeyFrame();

            // We allow points with high innovation (considererd outliers by the Huber Function)
            // pass to the new keyframe, so that bundle adjustment will finally decide
            // if they are outliers or not. We don't want next frame to estimate its position
            // with those points so we discard them in the frame.
            for(int i=0; i<mCurrentFrame.N;i++)
            {
                if(mCurrentFrame.mvpMapPoints[i] && mCurrentFrame.mvbOutlier[i])
                    mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
            }
        }

        // Reset if the camera get lost soon after initialization
        if(mState==LOST)
        {
            if(mpMap->KeyFramesInMap()<=5)
            {
                cout << "Track lost soon after initialisation, reseting..." << endl;
                mpSystem->Reset();
                return;
            }
        }

        if(!mCurrentFrame.mpReferenceKF)
            mCurrentFrame.mpReferenceKF = mpReferenceKF;

        mLastFrame = Frame(mCurrentFrame);
    }

    // Store frame pose information to retrieve the complete camera trajectory afterwards.
    if(!mCurrentFrame.mTcw.empty())
    {
        cv::Mat Tcr = mCurrentFrame.mTcw*mCurrentFrame.mpReferenceKF->GetPoseInverse();
        mlRelativeFramePoses.push_back(Tcr);
        mlpReferences.push_back(mpReferenceKF);
        mlFrameTimes.push_back(mCurrentFrame.mTimeStamp);
        mlbLost.push_back(mState==LOST);
    }
    else
    {
        // This can happen if tracking is lost
        mlRelativeFramePoses.push_back(mlRelativeFramePoses.back());
        mlpReferences.push_back(mlpReferences.back());
        mlFrameTimes.push_back(mlFrameTimes.back());
        mlbLost.push_back(mState==LOST);
    }

}

//点线特征同时使用，完成相机位姿跟踪
void Tracking::TrackBoth()
{
    //判断当前状态是否为无图像，若是则设置为未初始化
    if(mState==NO_IMAGES_YET)
    {
        mState = NOT_INITIALIZED;
    }

    mLastProcessedState=mState;

    // Get Map Mutex -> Map cannot be changed
    unique_lock<mutex> lock(mpMap->mMutexMapUpdate);

    if(mState==NOT_INITIALIZED)
    {

        if(mSensor==System::MONOCULAR && mbusingLine == true)
	    MonocularInitializationBoth();

        mpFrameDrawer->UpdateBoth(this);
	
        if(mState!=OK)
            return;
    }
    else
    {
        // System is initialized. Track Frame.
        bool bOK;

        // Initial camera pose estimation using motion model or relocalization (if tracking is lost)
	// 使用地图对应的跟踪流程
        if(!mbOnlyTracking)
        {
            // Local Mapping is activated. This is the normal behaviour, unless
            // you explicitly activate the "only tracking" mode.

            if(mState==OK)
            {
                // Local Mapping might have changed some MapPoints tracked in last frame
                CheckReplacedInLastFrame();
		CheckReplacedInLastFrameLines();

		//若初始化刚完成或者刚重定位成功，这两种情况下均没有速度模型
                if(mVelocity.empty() || mCurrentFrame.mnId<mnLastRelocFrameId+2)
                {    
                    bOK = TrackReferenceKeyFrameBoth();
		    cout << "初始化刚完成或者刚重定位成功,执行参考关键帧跟踪成功： " << true << endl;
		}
                else
                {
		    //在有运动模型的情况下优先使用，因为速度快
                    bOK = TrackWithMotionModelBoth();
		    cout << "运动模型跟踪成功： " << true << endl;
		    
                    if(!bOK)
		    {
			cout << "运动模型失败,执行参考关键帧跟踪： " << true << endl;
                        bOK = TrackReferenceKeyFrameBoth();	
		    }
                }
            }
            else
            {
		//若当前的局部状态为LOST，则执行点线特征同时使用的重定位算法
		cout << "执行重定位: " << true << endl << endl;
                bOK = RelocalizationBoth();	
            }
        }
        else
        {
            // Localization Mode: Local Mapping is deactivated
            if(mState==LOST)
            {
                bOK = RelocalizationBoth();
            }
            else
            {	

                if(!mbVO)
                {
                    // In last frame we tracked enough MapPoints in the map

                    if(!mVelocity.empty())
                    {
                        bOK = TrackWithMotionModelBoth();
                    }
                    else
                    {
                        bOK = TrackReferenceKeyFrameBoth();
                    }
                }
                else
                {
                    // In last frame we tracked mainly "visual odometry" points.

                    // We compute two camera poses, one from motion model and one doing relocalization.
                    // If relocalization is sucessfull we choose that solution, otherwise we retain
                    // the "visual odometry" solution.

                    bool bOKMM = false;
                    bool bOKReloc = false;
                    vector<MapPoint*> vpMPsMM;
                    vector<MapLine*> vpMLsMM;
                    vector<bool> vbOutMM;
                    vector<bool> vbOutLineMM;
                    cv::Mat TcwMM;
		    
                    if(!mVelocity.empty())
                    {
                        bOKMM = TrackWithMotionModelBoth();
                        vpMPsMM = mCurrentFrame.mvpMapPoints;
			vpMLsMM =  mCurrentFrame.mvpMapLines;
                        vbOutMM = mCurrentFrame.mvbOutlier;
			vbOutLineMM =  mCurrentFrame.mvbOutlierLines;
                        TcwMM = mCurrentFrame.mTcw.clone();
                    }
                    bOKReloc = RelocalizationBoth();

                    if(bOKMM && !bOKReloc)
                    {
                        mCurrentFrame.SetPose(TcwMM);
                        mCurrentFrame.mvpMapPoints = vpMPsMM;
                        mCurrentFrame.mvbOutlier = vbOutMM;
			mCurrentFrame.mvpMapLines = vpMLsMM;
			mCurrentFrame.mvbOutlierLines = vbOutLineMM;

                        if(mbVO)
                        {

                            for(int i =0; i<mCurrentFrame.N; i++)
                            {
                                if(mCurrentFrame.mvpMapPoints[i] && !mCurrentFrame.mvbOutlier[i])
                                {
                                    mCurrentFrame.mvpMapPoints[i]->IncreaseFound();
                                }
                            }

                            for(int i =0; i<mCurrentFrame.NL; i++)
                            {
                                if(mCurrentFrame.mvpMapLines[i] && !mCurrentFrame.mvbOutlierLines[i])
                                {
                                    mCurrentFrame.mvpMapLines[i]->IncreaseFound();
                                }
                            }
                        }
                    }
                    else if(bOKReloc)
                    {
                        mbVO = false;
                    }

                    bOK = bOKReloc || bOKMM;
                }
            }
        }

        mCurrentFrame.mpReferenceKF = mpReferenceKF;
	mCurrentFrame.mpReferenceKFLines = mpReferenceKFLines;
	
        // If we have an initial estimation of the camera pose and matching. Track the local map.
        if(!mbOnlyTracking)
        {
            if(bOK)
		//跟踪当前图像帧的局部地图
                bOK = TrackLocalMapBoth();	    
        }
        else
        {
            // mbVO true means that there are few matches to MapPoints/MapLines in the map. We cannot retrieve
            // a local map and therefore we do not perform TrackLocalMap(). Once the system relocalizes
            // the camera we will use the local map again.
            if(bOK && !mbVO)
                bOK = TrackLocalMapBoth();
        }

        if(bOK)
            mState = OK;
        else
            mState=LOST;

        // Update drawer
        mpFrameDrawer->UpdateBoth(this);

        // If tracking were good, check if we insert a keyframe
        if(bOK)
        {
            // Update motion model
            if(!mLastFrame.mTcw.empty())
            {
                cv::Mat LastTwc = cv::Mat::eye(4,4,CV_32F);
                mLastFrame.GetRotationInverse().copyTo(LastTwc.rowRange(0,3).colRange(0,3));
                mLastFrame.GetCameraCenter().copyTo(LastTwc.rowRange(0,3).col(3));
                mVelocity = mCurrentFrame.mTcw*LastTwc;
            }
            else
                mVelocity = cv::Mat();

            mpMapDrawer->SetCurrentCameraPose(mCurrentFrame.mTcw);

            // Clean VO matches
            for(int i=0; i<mCurrentFrame.N; i++)
            {
                MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];
                if(pMP)
                    if(pMP->Observations()<1)
                    {
                        mCurrentFrame.mvbOutlier[i] = false;
                        mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
                    }
            }
            for(int i=0; i<mCurrentFrame.NL; i++)
            {
                MapLine* pML = mCurrentFrame.mvpMapLines[i];
                if(pML)
                    if(pML->Observations()<1)
                    {
                        mCurrentFrame.mvbOutlierLines[i] = false;
                        mCurrentFrame.mvpMapLines[i]=static_cast<MapLine*>(NULL);
                    }
            }

            // Check if we need to insert a new keyframe
            if(NeedNewKeyFrameBoth())
	    {
		CreateNewKeyFrameBoth();
	    }
                

            // We allow points with high innovation (considererd outliers by the Huber Function)
            // pass to the new keyframe, so that bundle adjustment will finally decide
            // if they are outliers or not. We don't want next frame to estimate its position
            // with those points so we discard them in the frame.
            for(int i=0; i<mCurrentFrame.N;i++)
            {
                if(mCurrentFrame.mvpMapPoints[i] && mCurrentFrame.mvbOutlier[i])
                    mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
            }
             for(int i=0; i<mCurrentFrame.NL;i++)
            {
                if(mCurrentFrame.mvpMapLines[i] && mCurrentFrame.mvbOutlierLines[i])
                    mCurrentFrame.mvpMapLines[i]=static_cast<MapLine*>(NULL);
            }
        }

        // Reset if the camera get lost soon after initialization
        if(mState==LOST)
        {
            if(mpMap->KeyFramesInMap()<=5)
            {
                cout << "Track lost soon after initialisation, reseting..." << endl;
                mpSystem->Reset();
                return;
            }
        }

        if(!mCurrentFrame.mpReferenceKF)
            mCurrentFrame.mpReferenceKF = mpReferenceKF;
	
	if(!mCurrentFrame.mpReferenceKFLines)
            mCurrentFrame.mpReferenceKFLines = mpReferenceKFLines;
	
        mLastFrame = Frame(mCurrentFrame);
    }

    // Store frame pose information to retrieve the complete camera trajectory afterwards.
    if(!mCurrentFrame.mTcw.empty())
    {
        cv::Mat Tcr = mCurrentFrame.mTcw*mCurrentFrame.mpReferenceKF->GetPoseInverse();
        mlRelativeFramePoses.push_back(Tcr);
        mlpReferences.push_back(mpReferenceKF);
        mlFrameTimes.push_back(mCurrentFrame.mTimeStamp);
        mlbLost.push_back(mState==LOST);
    }
    else
    {
        // This can happen if tracking is lost
        mlRelativeFramePoses.push_back(mlRelativeFramePoses.back());
        mlpReferences.push_back(mlpReferences.back());
        mlFrameTimes.push_back(mlFrameTimes.back());
        mlbLost.push_back(mState==LOST);
    }

}

void Tracking::StereoInitialization()
{
    if(mCurrentFrame.N>500)
    {
        // Set Frame pose to the origin
        mCurrentFrame.SetPose(cv::Mat::eye(4,4,CV_32F));

        // Create KeyFrame
        KeyFrame* pKFini = new KeyFrame(mCurrentFrame,mpMap,mpKeyFrameDB);

        // Insert KeyFrame in the map
        mpMap->AddKeyFrame(pKFini);

        // Create MapPoints and asscoiate to KeyFrame
        for(int i=0; i<mCurrentFrame.N;i++)
        {
            float z = mCurrentFrame.mvDepth[i];
            if(z>0)
            {
                cv::Mat x3D = mCurrentFrame.UnprojectStereo(i);
                MapPoint* pNewMP = new MapPoint(x3D,pKFini,mpMap);
                pNewMP->AddObservation(pKFini,i);
                pKFini->AddMapPoint(pNewMP,i);
                pNewMP->ComputeDistinctiveDescriptors();
                pNewMP->UpdateNormalAndDepth();
                mpMap->AddMapPoint(pNewMP);

                mCurrentFrame.mvpMapPoints[i]=pNewMP;
            }
        }

        cout << "New map created with " << mpMap->MapPointsInMap() << " points" << endl;

        mpLocalMapper->InsertKeyFrame(pKFini);

        mLastFrame = Frame(mCurrentFrame);
        mnLastKeyFrameId=mCurrentFrame.mnId;
        mpLastKeyFrame = pKFini;

        mvpLocalKeyFrames.push_back(pKFini);
        mvpLocalMapPoints=mpMap->GetAllMapPoints();
        mpReferenceKF = pKFini;
        mCurrentFrame.mpReferenceKF = pKFini;

        mpMap->SetReferenceMapPoints(mvpLocalMapPoints);

        mpMap->mvpKeyFrameOrigins.push_back(pKFini);

        mpMapDrawer->SetCurrentCameraPose(mCurrentFrame.mTcw);

        mState=OK;
    }
}

//Tracking线程中函数，间接调用匹配器中的匹配函数
//点匹配器：
void Tracking::SearchForInitialization(Frame &F1, Frame &F2, std::vector<cv::Point2f> &vbPrevMatched, std::vector<int> &vnMatches12,
				 int windowSize, int &nmatches)
{
    nmatches = mpORBmatcherLeftIni->SearchForInitialization(F1, F2, vbPrevMatched, vnMatches12, windowSize);
}

void Tracking::SearchByProjectionMotion(Frame &CurrentFrame, const Frame &LastFrame, const float th,
			    const bool bMono, int &nmatches)
{
    nmatches = mpORBmatcherLeftmotion->SearchByProjection( CurrentFrame, LastFrame, th, bMono);
}

void Tracking::SearchByBoW(KeyFrame *pKF, Frame &F, std::vector<MapPoint*> &vpMapPointMatches, int &nmatches)
{
    F.ComputeBoW();
    nmatches = mpORBmatcherLeftBoW->SearchByBoW(pKF, F, vpMapPointMatches);
}
    
void Tracking::SearchByProjection(Frame &CurrentFrame, KeyFrame* pKF, const std::set<MapPoint*> &sAlreadyFound,
			    const float th, const int ORBdist, int &nmatches)
{
    nmatches = mpORBmatcherLeftRel->SearchByProjection(CurrentFrame, pKF, sAlreadyFound, th, ORBdist);
}

//线匹配器：

void Tracking::SearchForInitializationLines(Frame &F1, Frame &F2, std::vector<cv::Point2f> &vbPrevMatched, std::vector<int> &vnMatches12,
				 int windowSize, int &nmatches)
{
    nmatches = mpLinematcherLeftIni->SearchForInitialization(F1, F2, vbPrevMatched, vnMatches12, windowSize);
}

void Tracking::SearchByProjectionLinesMotion(Frame &CurrentFrame, const Frame &LastFrame, const float th,
			    const bool bMono, int &nmatches)
{   
    nmatches = mpLinematcherLeftmotion->SearchByProjection(CurrentFrame, LastFrame, th, bMono);
}

void Tracking::SearchByKNNLines(KeyFrame *pKF, Frame &F, std::vector<MapLine*> &vpMapLineMatches, int &nmatches12)
{
    nmatches12 = mpLinematcherLeftKNN->SearchByKNN(pKF, F, vpMapLineMatches);
}

void Tracking::SearchByProjectionLines(Frame &CurrentFrame, KeyFrame* pKF, const std::set<MapLine*> &sAlreadyFound, 
			    const float th, const int Linedist, int &nmatches)
{
    nmatches = mpLinematcherLeftRel->SearchByProjection(CurrentFrame, pKF, sAlreadyFound, th, Linedist);
}

//点特征初始化器
void Tracking::MonocularInitialization()
{

    if(!mpInitializer)
    {
        // Set Reference Frame
        if(mCurrentFrame.mvKeys.size()>100)
        {
            mInitialFrame = Frame(mCurrentFrame);
            mLastFrame = Frame(mCurrentFrame);
            mvbPrevMatched.resize(mCurrentFrame.mvKeysUn.size());
            for(size_t i=0; i<mCurrentFrame.mvKeysUn.size(); i++)
                mvbPrevMatched[i]=mCurrentFrame.mvKeysUn[i].pt;

            if(mpInitializer)
                delete mpInitializer;

            mpInitializer =  new Initializer(mCurrentFrame,1.0,200,mbusingLine);

            fill(mvIniMatches.begin(),mvIniMatches.end(),-1);

            return;
        }
    }
    else
    {
        // Try to initialize
        if((int)mCurrentFrame.mvKeys.size()<=100)
        {
            delete mpInitializer;
            mpInitializer = static_cast<Initializer*>(NULL);
            fill(mvIniMatches.begin(),mvIniMatches.end(),-1);
            return;
        }

        // Find correspondences
        ORBmatcher matcher(0.9,true);
        int nmatches = matcher.SearchForInitialization(mInitialFrame,mCurrentFrame,mvbPrevMatched,mvIniMatches,100);

        // Check if there are enough correspondences
        if(nmatches<100)
        {
            delete mpInitializer;
            mpInitializer = static_cast<Initializer*>(NULL);
            return;
        }

        cv::Mat Rcw; // Current Camera Rotation
        cv::Mat tcw; // Current Camera Translation
        vector<bool> vbTriangulated; // Triangulated Correspondences (mvIniMatches)

        Timer timer;
	timer.start();
	
        if(mpInitializer->Initialize(mCurrentFrame, mvIniMatches, Rcw, tcw, mvIniP3D, vbTriangulated))
        {
	    
	  double InitializationTime = timer.stop();
    
	    cout << "初始化时间: " << InitializationTime << " ms" <<endl <<endl ;
	  
            for(size_t i=0, iend=mvIniMatches.size(); i<iend;i++)
            {
                if(mvIniMatches[i]>=0 && !vbTriangulated[i])
                {
                    mvIniMatches[i]=-1;
                    nmatches--;
                }
            }

            // Set Frame Poses
            mInitialFrame.SetPose(cv::Mat::eye(4,4,CV_32F));
            cv::Mat Tcw = cv::Mat::eye(4,4,CV_32F);
            Rcw.copyTo(Tcw.rowRange(0,3).colRange(0,3));
            tcw.copyTo(Tcw.rowRange(0,3).col(3));
            mCurrentFrame.SetPose(Tcw);

            CreateInitialMapMonocular();
        }
    }
}

//使用点线特征进行单目初始化
void Tracking::MonocularInitializationBoth()
{
    if(!mpInitializer)
    {
        // Set Reference Frame
        if(mCurrentFrame.mvKeys.size()>100 || mCurrentFrame.mvLines.size()>80)
	{
            mInitialFrame = Frame(mCurrentFrame);
            mLastFrame = Frame(mCurrentFrame);
	    
	    //预设第一图像帧的关键点和关键线中点向量
            mvbPrevMatched.resize(mCurrentFrame.mvKeysUn.size());    
            for(size_t i=0; i<mCurrentFrame.mvKeysUn.size(); i++)
                mvbPrevMatched[i]=mCurrentFrame.mvKeysUn[i].pt;
	    
	    mvbPrevMatchedMidPoints.resize(mCurrentFrame.mvMidPointsUn.size()); 
	    for(size_t i=0; i<mCurrentFrame.mvMidPointsUn.size(); i++)
                mvbPrevMatchedMidPoints[i]=mCurrentFrame.mvMidPointsUn[i].pt;

            if(mpInitializer)
                delete mpInitializer;

            mpInitializer =  new Initializer(mCurrentFrame,1.0,200,mbusingLine);

	    mvIniMatches.clear();
	    mvIniMatchesLines.clear();
	    
            return;
        }
    }
    else
    {
        // Try to initialize
        if((int)mCurrentFrame.mvKeys.size() <=100 && mCurrentFrame.mvLines.size()<=80)
        {
            delete mpInitializer;
            mpInitializer = static_cast<Initializer*>(NULL);
	    mvIniMatches.clear();
	    mvIniMatchesLines.clear();
            return;
        }
	
	int nmatches;
	int nmatchLines;
	
	thread threadSearchForInitialization(&Tracking::SearchForInitialization, this, ref(mInitialFrame), 
	 ref(mCurrentFrame), ref(mvbPrevMatched), ref(mvIniMatches), 100, ref(nmatches));
	
	thread threadSearchForInitializationLines(&Tracking::SearchForInitializationLines, this, ref(mInitialFrame), 
	 ref(mCurrentFrame), ref(mvbPrevMatchedMidPoints), ref(mvIniMatchesLines), 100, ref(nmatchLines));

	threadSearchForInitialization.join();
	threadSearchForInitializationLines.join();
	
        // Check if there are enough correspondences
        if(nmatches<=100 && nmatchLines<=80)
        {
            delete mpInitializer;
            mpInitializer = static_cast<Initializer*>(NULL);
            return;
        }

        cv::Mat Rcw; // Current Camera Rotation
        cv::Mat tcw; // Current Camera Translation
        vector<bool> vbTriangulated; // Triangulated Correspondences (mvIniMatches)
        vector<bool> vbTriangulatedMidP;

	Timer timer;
	timer.start();
	
        if(mpInitializer->InitializeBoth(mCurrentFrame, mvIniMatches,mvIniMatchesLines, Rcw, tcw, mvIniP3D,
	  mvIniMidP3D, mvIniFirP3D, mvIniEndP3D, vbTriangulated, vbTriangulatedMidP))
        {
	    
	    double InitializationTime = timer.stop();
    
	    cout << "初始化时间: " << InitializationTime << " ms" <<endl <<endl ;
	  
            for(size_t i=0, iend=mvIniMatches.size(); i<iend; i++)
            {
                if(mvIniMatches[i]>=0 && !vbTriangulated[i])
                {
                    mvIniMatches[i]=-1;
                    nmatches--;
                }
            }
            
            cout << "nmatches: " << nmatches << endl << endl;            	    
	    
            for(size_t iL=0, iLend=mvIniMatchesLines.size(); iL<iLend; iL++)
            {
		//输出每一个存在匹配位置的Bool值
		//cout << "vbTriangulatedMidP[iL]: " << vbTriangulatedMidP[iL] << endl;
		
                if(mvIniMatchesLines[iL]>=0 && !vbTriangulatedMidP[iL])
                {
                    mvIniMatchesLines[iL]=-1;
                    nmatchLines--;
                }
            }
            
	    cout << "nmatchLines: " << nmatchLines << endl << endl;
	    	               
            // Set Frame Poses
            mInitialFrame.SetPose(cv::Mat::eye(4,4,CV_32F));
            cv::Mat Tcw = cv::Mat::eye(4,4,CV_32F);
            Rcw.copyTo(Tcw.rowRange(0,3).colRange(0,3));
            tcw.copyTo(Tcw.rowRange(0,3).col(3));
            mCurrentFrame.SetPose(Tcw);
	    
            CreateInitialMapMonocularBoth();
        }
        
    }
    
}
//点特征
void Tracking::CreateInitialMapMonocular()
{
    // Create KeyFrames
    KeyFrame* pKFini = new KeyFrame(mInitialFrame,mpMap,mpKeyFrameDB);
    KeyFrame* pKFcur = new KeyFrame(mCurrentFrame,mpMap,mpKeyFrameDB);


    pKFini->ComputeBoW();
    pKFcur->ComputeBoW();

    // Insert KFs in the map
    mpMap->AddKeyFrame(pKFini);
    mpMap->AddKeyFrame(pKFcur);

    // Create MapPoints and asscoiate to keyframes
    for(size_t i=0; i<mvIniMatches.size();i++)
    {
        if(mvIniMatches[i]<0)
            continue;

        //Create MapPoint.
        cv::Mat worldPos(mvIniP3D[i]);

        MapPoint* pMP = new MapPoint(worldPos,pKFcur,mpMap);

        pKFini->AddMapPoint(pMP,i);
        pKFcur->AddMapPoint(pMP,mvIniMatches[i]);

        pMP->AddObservation(pKFini,i);
        pMP->AddObservation(pKFcur,mvIniMatches[i]);

        pMP->ComputeDistinctiveDescriptors();
        pMP->UpdateNormalAndDepth();

        //Fill Current Frame structure
        mCurrentFrame.mvpMapPoints[mvIniMatches[i]] = pMP;
        mCurrentFrame.mvbOutlier[mvIniMatches[i]] = false;

        //Add to Map
        mpMap->AddMapPoint(pMP);
    }

    // Update Connections
    pKFini->UpdateConnections();
    pKFcur->UpdateConnections();

    // Bundle Adjustment
    cout << "New Map created with " << mpMap->MapPointsInMap() << " points" << endl;

    Optimizer::GlobalBundleAdjustemnt(mpMap,20);

    // Set median depth to 1
    float medianDepth = pKFini->ComputeSceneMedianDepth(2);
    float invMedianDepth = 1.0f/medianDepth;

    if(medianDepth<0 || pKFcur->TrackedMapPoints(1)<100)
    {
        cout << "Wrong initialization, reseting..." << endl;
        Reset();
        return;
    }

    // Scale initial baseline
    cv::Mat Tc2w = pKFcur->GetPose();
    Tc2w.col(3).rowRange(0,3) = Tc2w.col(3).rowRange(0,3)*invMedianDepth;
    pKFcur->SetPose(Tc2w);

    // Scale points
    vector<MapPoint*> vpAllMapPoints = pKFini->GetMapPointMatches();
    for(size_t iMP=0; iMP<vpAllMapPoints.size(); iMP++)
    {
        if(vpAllMapPoints[iMP])
        {
            MapPoint* pMP = vpAllMapPoints[iMP];
            pMP->SetWorldPos(pMP->GetWorldPos()*invMedianDepth);
        }
    }

    mpLocalMapper->InsertKeyFrame(pKFini);
    mpLocalMapper->InsertKeyFrame(pKFcur);

    mCurrentFrame.SetPose(pKFcur->GetPose());
    mnLastKeyFrameId=mCurrentFrame.mnId;
    mpLastKeyFrame = pKFcur;

    mvpLocalKeyFrames.push_back(pKFcur);
    mvpLocalKeyFrames.push_back(pKFini);
    mvpLocalMapPoints=mpMap->GetAllMapPoints();
    mpReferenceKF = pKFcur;
    mCurrentFrame.mpReferenceKF = pKFcur;

    mLastFrame = Frame(mCurrentFrame);

    mpMap->SetReferenceMapPoints(mvpLocalMapPoints);

    mpMapDrawer->SetCurrentCameraPose(pKFcur->GetPose());

    mpMap->mvpKeyFrameOrigins.push_back(pKFini);

    mState=OK;
}

//点线特征同时使用
void Tracking::CreateInitialMapMonocularBoth()
{
    // Create KeyFrames
    KeyFrame* pKFini = new KeyFrame(mInitialFrame,mpMap,mpKeyFrameDB);
    KeyFrame* pKFcur = new KeyFrame(mCurrentFrame,mpMap,mpKeyFrameDB);


    pKFini->ComputeBoW();
    pKFcur->ComputeBoW();

    // Insert KFs in the map
    mpMap->AddKeyFrame(pKFini);
    mpMap->AddKeyFrame(pKFcur);

    // Create MapPoints and asscoiate to keyframes
    for(size_t i=0; i<mvIniMatches.size();i++)
    {
        if(mvIniMatches[i]<0)
            continue;

        //Create MapPoint.
        cv::Mat worldPos(mvIniP3D[i]);

        MapPoint* pMP = new MapPoint(worldPos,pKFcur,mpMap);

        pKFini->AddMapPoint(pMP,i);
        pKFcur->AddMapPoint(pMP,mvIniMatches[i]);

        pMP->AddObservation(pKFini,i);
        pMP->AddObservation(pKFcur,mvIniMatches[i]);

        pMP->ComputeDistinctiveDescriptors();
        pMP->UpdateNormalAndDepth();

        //Fill Current Frame structure
        mCurrentFrame.mvpMapPoints[mvIniMatches[i]] = pMP;
        mCurrentFrame.mvbOutlier[mvIniMatches[i]] = false;

        //Add to Map
        mpMap->AddMapPoint(pMP);
    }
   
    // Update Connections
    pKFini->UpdateConnections();
    pKFcur->UpdateConnections();
    
    for(size_t i=0; i<mvIniMatchesLines.size();i++)
    {
        if(mvIniMatchesLines[i]<0)
            continue;
	
        //Create MapPoint.
        cv::Mat worldFirPPos(mvIniFirP3D[i]);
	cv::Mat worldEndPPos(mvIniEndP3D[i]);
	cv::Mat worldMidPPos(mvIniMidP3D[i]);

        MapLine* pML = new MapLine(worldFirPPos, worldEndPPos,worldMidPPos,pKFcur,mpMap);
	
        pKFini->AddMapLine(pML,i);
        pKFcur->AddMapLine(pML,mvIniMatchesLines[i]);

        pML->AddObservation(pKFini,i);
        pML->AddObservation(pKFcur,mvIniMatchesLines[i]);
	
        pML->ComputeDistinctiveDescriptors();	
        pML->UpdateNormalAndDepth();	
	pML->Update2DLineLength();  
	
        //Fill Current Frame structure
        mCurrentFrame.mvpMapLines[mvIniMatchesLines[i]] = pML;
        mCurrentFrame.mvbOutlierLines[mvIniMatchesLines[i]] = false;

        //Add to Map
        mpMap->AddMapLine(pML);
    }

    // Update Connections
    pKFini->UpdateConnectionsLines();
    pKFcur->UpdateConnectionsLines();

    // Bundle Adjustment
    cout << "New Map created with " << mpMap->MapPointsInMap() << " points" << endl;
    cout << "New Map created with " <<  mpMap->MapLinesInMap() << " lines" << endl << endl;

    Optimizer::GlobalBundleAdjustemntIni(mpMap,10,true);

    float medianDepth = pKFini->ComputeSceneMedianDepthBoth(2);
    float invMedianDepth = 1.0f/medianDepth;

    if(medianDepth<0 || (pKFcur->TrackedMapPoints(1)<100 && pKFcur->TrackedMapLines(1)<80 && 
      (pKFcur->TrackedMapPoints(1)+pKFcur->TrackedMapLines(1))<150 ))
    {
        cout << "Wrong initialization, reseting..." << endl;
        ResetBoth();
        return;
    }

    // Scale initial baseline
    cv::Mat Tc2w = pKFcur->GetPose();
    Tc2w.col(3).rowRange(0,3) = Tc2w.col(3).rowRange(0,3)*invMedianDepth;
    pKFcur->SetPose(Tc2w);

    // Scale points
    vector<MapPoint*> vpAllMapPoints = pKFini->GetMapPointMatches();
    for(size_t iMP=0; iMP<vpAllMapPoints.size(); iMP++)
    {
        if(vpAllMapPoints[iMP])
        {
            MapPoint* pMP = vpAllMapPoints[iMP];
            pMP->SetWorldPos(pMP->GetWorldPos()*invMedianDepth);
        }
    }

    vector<MapLine*> vpAllMapLines = pKFini->GetMapLineMatches();
    for(size_t iML=0; iML<vpAllMapLines.size(); iML++)
    {
        if(vpAllMapLines[iML])
        {
            MapLine* pML = vpAllMapLines[iML];
            pML->SetFirPWorldPos(pML->GetFirPWorldPos()*invMedianDepth);
	    pML->SetEndPWorldPos(pML->GetEndPWorldPos()*invMedianDepth);
	    pML->SetMidPWorldPos(pML->GetMidPWorldPos()*invMedianDepth);
        }
    }
     
    mpLocalMapper->InsertKeyFrame(pKFini);
    mpLocalMapper->InsertKeyFrame(pKFcur);

    mCurrentFrame.SetPose(pKFcur->GetPose());
    mnLastKeyFrameId=mCurrentFrame.mnId;
    mpLastKeyFrame = pKFcur;

    mvpLocalKeyFrames.push_back(pKFcur);
    mvpLocalKeyFrames.push_back(pKFini);
    mvpLocalMapPoints=mpMap->GetAllMapPoints();
    mvpLocalMapLines=mpMap->GetAllMapLines();
    
    mpReferenceKF = pKFcur;
    mpReferenceKFLines = pKFcur;
    mCurrentFrame.mpReferenceKF = pKFcur;
    mCurrentFrame.mpReferenceKFLines = pKFcur;
    
    mLastFrame = Frame(mCurrentFrame);

    mpMap->SetReferenceMapPoints(mvpLocalMapPoints);
    mpMap->SetReferenceMapLines(mvpLocalMapLines);
    
    mpMapDrawer->SetCurrentCameraPose(pKFcur->GetPose());

    mpMap->mvpKeyFrameOrigins.push_back(pKFini);

    mState=OK;
    
    cout << "初始化成功:" << true << endl << endl;
    
}

void Tracking::CheckReplacedInLastFrame()
{
    for(int i =0; i<mLastFrame.N; i++)
    {
        MapPoint* pMP = mLastFrame.mvpMapPoints[i];

        if(pMP)
        {
            MapPoint* pRep = pMP->GetReplaced();
            if(pRep)
            {
                mLastFrame.mvpMapPoints[i] = pRep;
            }
        }
    }
}
//检测后端替换的地图线特征
void Tracking::CheckReplacedInLastFrameLines()
{
    for(int i =0; i<mLastFrame.NL; i++)
    {
        MapLine* pML = mLastFrame.mvpMapLines[i];

        if(pML)
        {
            MapLine* pReL = pML->GetReplaced();
            if(pReL)
            {
                mLastFrame.mvpMapLines[i] = pReL;
            }
        }
    }
}

bool Tracking::TrackReferenceKeyFrame()
{
    // Compute Bag of Words vector
    mCurrentFrame.ComputeBoW();

    // We perform first an ORB matching with the reference keyframe
    // If enough matches are found we setup a PnP solver
    ORBmatcher matcher(0.7,true);
    vector<MapPoint*> vpMapPointMatches;

    int nmatches = matcher.SearchByBoW(mpReferenceKF,mCurrentFrame,vpMapPointMatches);

    if(nmatches<15)
        return false;

    mCurrentFrame.mvpMapPoints = vpMapPointMatches;
    mCurrentFrame.SetPose(mLastFrame.mTcw);

    Optimizer::PoseOptimization(&mCurrentFrame);

    // Discard outliers
    int nmatchesMap = 0;
    for(int i =0; i<mCurrentFrame.N; i++)
    {
        if(mCurrentFrame.mvpMapPoints[i])
        {
            if(mCurrentFrame.mvbOutlier[i])
            {
                MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];

                mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
                mCurrentFrame.mvbOutlier[i]=false;
                pMP->mbTrackInView = false;
                pMP->mnLastFrameSeen = mCurrentFrame.mnId;
                nmatches--;
            }
            else if(mCurrentFrame.mvpMapPoints[i]->Observations()>0)
                nmatchesMap++;
        }
    }

    return nmatchesMap>=10;
}
//点线特征同时使用的参考关键帧跟踪函数
bool Tracking::TrackReferenceKeyFrameBoth()
{
     
    vector<MapPoint*> vpMapPointMatches;
    vector<MapLine*> vpMapLineMatches;

    int nmatches;
    int nmatchLines;

    thread threadSearchByBoW(&Tracking::SearchByBoW, this, mpReferenceKF, ref(mCurrentFrame), 
			     ref(vpMapPointMatches), ref(nmatches) );
    thread threadSearchByKNNLines(&Tracking::SearchByKNNLines, this, mpReferenceKFLines, 
				  ref(mCurrentFrame), ref(vpMapLineMatches), ref(nmatchLines) );
    
    threadSearchByBoW.join();
    threadSearchByKNNLines.join();
    
    if(nmatches<=12 && nmatchLines<=12)
        return false;

    mCurrentFrame.mvpMapPoints = vpMapPointMatches;
    mCurrentFrame.mvpMapLines = vpMapLineMatches;//当前图像帧保存地图线向来嗯
    mCurrentFrame.SetPose(mLastFrame.mTcw);

    int inlinersNum;
    int inlinerLinesNum;
    if(nmatches>12 && nmatchLines>12)
    {	
	Optimizer::PoseOptimizationmain(&mCurrentFrame, inlinersNum, inlinerLinesNum, 10);
    }
    else if(nmatches>15)
    {
	Optimizer::PoseOptimizationDoublePoints(&mCurrentFrame,inlinersNum);
	Optimizer::SetOutlierLinesForPose(&mCurrentFrame,inlinerLinesNum);
    }
    else if( nmatches>6 && nmatchLines>6 )
    {
	Optimizer::PoseOptimizationLowFeature(&mCurrentFrame, inlinersNum, inlinerLinesNum);
    }
    else
	return false;
            
    // Discard outliers
    mnMatchesInliers = 0;
    for(int i =0; i<mCurrentFrame.N; i++)
    {
        if(mCurrentFrame.mvpMapPoints[i])
        {
            if(mCurrentFrame.mvbOutlier[i])
            {
                MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];

                mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
                mCurrentFrame.mvbOutlier[i]=false;
                pMP->mbTrackInView = false;
                pMP->mnLastFrameSeen = mCurrentFrame.mnId;
                nmatches--;
            }
            else if(mCurrentFrame.mvpMapPoints[i]->Observations()>0)
                mnMatchesInliers++;
        }
    } 
    // Discard OutlierLines
    mnMatchesInlierLines = 0;
    for(int i =0; i<mCurrentFrame.NL; i++)
    {
        if(mCurrentFrame.mvpMapLines[i])
        {
            if(mCurrentFrame.mvbOutlierLines[i])
            {
                MapLine* pML = mCurrentFrame.mvpMapLines[i];

                mCurrentFrame.mvpMapLines[i]=static_cast<MapLine*>(NULL);
                mCurrentFrame.mvbOutlierLines[i]=false;
                pML->mbTrackInView = false;
                pML->mnLastFrameSeen = mCurrentFrame.mnId;
                nmatchLines--;
            }
            else if(mCurrentFrame.mvpMapLines[i]->Observations()>0)
                mnMatchesInlierLines++;
        }
    }
    
    if( (mnMatchesInliers>=8 && mnMatchesInlierLines>=8) )
	return true;
    else if(mnMatchesInliers>=10)
	return true; 
    else if(mnMatchesInlierLines+mnMatchesInliers>=10)
	return true; 
    else
	return false;
}

void Tracking::UpdateLastFrame()
{
    // Update pose according to reference keyframe
    KeyFrame* pRef = mLastFrame.mpReferenceKF;
    cv::Mat Tlr = mlRelativeFramePoses.back();

    mLastFrame.SetPose(Tlr*pRef->GetPose());

    if(mnLastKeyFrameId==mLastFrame.mnId || mSensor==System::MONOCULAR || !mbOnlyTracking)
        return;

    // Create "visual odometry" MapPoints
    // We sort points according to their measured depth by the stereo/RGB-D sensor
    vector<pair<float,int> > vDepthIdx;
    vDepthIdx.reserve(mLastFrame.N);
    for(int i=0; i<mLastFrame.N;i++)
    {
        float z = mLastFrame.mvDepth[i];
        if(z>0)
        {
            vDepthIdx.push_back(make_pair(z,i));
        }
    }

    if(vDepthIdx.empty())
        return;

    sort(vDepthIdx.begin(),vDepthIdx.end());

    // We insert all close points (depth<mThDepth)
    // If less than 100 close points, we insert the 100 closest ones.
    int nPoints = 0;
    for(size_t j=0; j<vDepthIdx.size();j++)
    {
        int i = vDepthIdx[j].second;

        bool bCreateNew = false;

        MapPoint* pMP = mLastFrame.mvpMapPoints[i];
        if(!pMP)
            bCreateNew = true;
        else if(pMP->Observations()<1)
        {
            bCreateNew = true;
        }

        if(bCreateNew)
        {
            cv::Mat x3D = mLastFrame.UnprojectStereo(i);
            MapPoint* pNewMP = new MapPoint(x3D,mpMap,&mLastFrame,i);

            mLastFrame.mvpMapPoints[i]=pNewMP;

            mlpTemporalPoints.push_back(pNewMP);
            nPoints++;
        }
        else
        {
            nPoints++;
        }

        if(vDepthIdx[j].first>mThDepth && nPoints>100)
            break;
    }
}

bool Tracking::TrackWithMotionModel()
{
    ORBmatcher matcher(0.9,true);

    // Update last frame pose according to its reference keyframe
    // Create "visual odometry" points if in Localization Mode
    UpdateLastFrame();

    mCurrentFrame.SetPose(mVelocity*mLastFrame.mTcw);

    fill(mCurrentFrame.mvpMapPoints.begin(),mCurrentFrame.mvpMapPoints.end(),static_cast<MapPoint*>(NULL));

    // Project points seen in previous frame
    int th;
    if(mSensor!=System::STEREO)
        th=15;
    else
        th=7;
    int nmatches = matcher.SearchByProjection(mCurrentFrame,mLastFrame,th,mSensor==System::MONOCULAR);

    // If few matches, uses a wider window search
    if(nmatches<20)
    {
        fill(mCurrentFrame.mvpMapPoints.begin(),mCurrentFrame.mvpMapPoints.end(),static_cast<MapPoint*>(NULL));
        nmatches = matcher.SearchByProjection(mCurrentFrame,mLastFrame,2*th,mSensor==System::MONOCULAR);
    }

    if(nmatches<20)
        return false;

    // Optimize frame pose with all matches
    Optimizer::PoseOptimization(&mCurrentFrame);

    // Discard outliers
    int nmatchesMap = 0;
    for(int i =0; i<mCurrentFrame.N; i++)
    {
        if(mCurrentFrame.mvpMapPoints[i])
        {
            if(mCurrentFrame.mvbOutlier[i])
            {
                MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];

                mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
                mCurrentFrame.mvbOutlier[i]=false;
                pMP->mbTrackInView = false;
                pMP->mnLastFrameSeen = mCurrentFrame.mnId;
                nmatches--;
            }
            else if(mCurrentFrame.mvpMapPoints[i]->Observations()>0)
                nmatchesMap++;
        }
    }    

    if(mbOnlyTracking)
    {
        mbVO = nmatchesMap<10;
        return nmatches>20;
    }

    return nmatchesMap>=10;
}
//点线特征同时使用的运动模型跟踪
bool Tracking::TrackWithMotionModelBoth()
{

    // Update last frame pose according to its reference keyframe
    UpdateLastFrame();

    mCurrentFrame.SetPose(mVelocity*mLastFrame.mTcw);
    
    fill(mCurrentFrame.mvpMapPoints.begin(),mCurrentFrame.mvpMapPoints.end(),static_cast<MapPoint*>(NULL));
    fill(mCurrentFrame.mvpMapLines.begin(),mCurrentFrame.mvpMapLines.end(),static_cast<MapLine*>(NULL));

    // Project points seen in previous frame
    int th=15;
    int thLines=30;

    int nmatches;
    int nmatchLines;
    
    //双线程分别进行点和线特征匹配
    thread threadSearchByProjection(&Tracking::SearchByProjectionMotion, this, ref(mCurrentFrame), ref(mLastFrame), 
				    th, mSensor==System::MONOCULAR, ref(nmatches) );
    thread threadSearchByProjectionLines(&Tracking::SearchByProjectionLinesMotion, this, ref(mCurrentFrame), ref(mLastFrame), 
					 thLines, mSensor==System::MONOCULAR, ref(nmatchLines) );               
    
    threadSearchByProjection.join();
    threadSearchByProjectionLines.join();
    
    // If few matches, uses a wider window search
    if(nmatches<=16 && nmatchLines<=12)
    {
        fill(mCurrentFrame.mvpMapPoints.begin(),mCurrentFrame.mvpMapPoints.end(),static_cast<MapPoint*>(NULL));
	fill(mCurrentFrame.mvpMapLines.begin(),mCurrentFrame.mvpMapLines.end(),static_cast<MapLine*>(NULL));
	
        thread threadSearchByProjection(&Tracking::SearchByProjectionMotion, this, ref(mCurrentFrame), ref(mLastFrame), 
				    2*th, mSensor==System::MONOCULAR, ref(nmatches) );
	thread threadSearchByProjectionLines(&Tracking::SearchByProjectionLinesMotion, this, ref(mCurrentFrame), ref(mLastFrame), 
					 2*thLines, mSensor==System::MONOCULAR, ref(nmatchLines) );               
    
	threadSearchByProjection.join();
	threadSearchByProjectionLines.join();
    }

    if(nmatches<=16 && nmatchLines<=12)
        return false;

    //定义返回的内点和内线参数引用
    int inlinersNum;
    int inlinerLinesNum;
    if(nmatches>16 && nmatchLines>12)
    {	
	Optimizer::PoseOptimizationmain(&mCurrentFrame, inlinersNum, inlinerLinesNum, 10);
	//Optimizer::PoseOptimizationLowFeature(&mCurrentFrame, inlinersNum, inlinerLinesNum);
    }
    else if(nmatches>20)
    {
	Optimizer::PoseOptimizationDoublePoints(&mCurrentFrame,inlinersNum);
	Optimizer::SetOutlierLinesForPose(&mCurrentFrame,inlinerLinesNum);
    }
    else if( nmatches>6 && nmatchLines>6 )
    {
	Optimizer::PoseOptimizationLowFeature(&mCurrentFrame, inlinersNum, inlinerLinesNum);
    }
    else
	return false;
    
    // Discard outliers
    mnMatchesInliers = 0;
    for(int i =0; i<mCurrentFrame.N; i++)
    {
        if(mCurrentFrame.mvpMapPoints[i])
        {
            if(mCurrentFrame.mvbOutlier[i])
            {
                MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];

                mCurrentFrame.mvpMapPoints[i]=static_cast<MapPoint*>(NULL);
                mCurrentFrame.mvbOutlier[i]=false;
                pMP->mbTrackInView = false;
                pMP->mnLastFrameSeen = mCurrentFrame.mnId;
                nmatches--;
            }
            else if(mCurrentFrame.mvpMapPoints[i]->Observations()>0)
                mnMatchesInliers++;
        }
    }    

    // Discard OutlierLines
    mnMatchesInlierLines = 0;
    for(int i =0; i<mCurrentFrame.NL; i++)
    {
        if(mCurrentFrame.mvpMapLines[i])
        {
            if(mCurrentFrame.mvbOutlierLines[i])
            {
                MapLine* pML = mCurrentFrame.mvpMapLines[i];

                mCurrentFrame.mvpMapLines[i]=static_cast<MapLine*>(NULL);
                mCurrentFrame.mvbOutlierLines[i]=false;
                pML->mbTrackInView = false;
                pML->mnLastFrameSeen = mCurrentFrame.mnId;
                nmatchLines--;
            }
            else if(mCurrentFrame.mvpMapLines[i]->Observations()>0)
                mnMatchesInlierLines++;
        }
    }

    if(mbOnlyTracking)
    {

        mbVO = (mnMatchesInliers<8 && mnMatchesInlierLines<8) ;
	
        return (nmatches>16 || nmatchLines>16);
    }
    
    //所有匹配有1/2为有效跟踪，则跟踪成功
    if( (mnMatchesInliers>=8 && mnMatchesInlierLines>=8) )
	return true;
    else if(mnMatchesInliers>=10)
	return true; 
    else if(mnMatchesInlierLines+mnMatchesInliers>=10)
	return true; 
    else
	return false;
}

bool Tracking::TrackLocalMap()
{
    // We have an estimation of the camera pose and some map points tracked in the frame.
    // We retrieve the local map and try to find matches to points in the local map.

    UpdateLocalMap();

    SearchLocalPoints();

    // Optimize Pose
    Optimizer::PoseOptimization(&mCurrentFrame);
    mnMatchesInliers = 0;

    // Update MapPoints Statistics
    for(int i=0; i<mCurrentFrame.N; i++)
    {
        if(mCurrentFrame.mvpMapPoints[i])
        {
            if(!mCurrentFrame.mvbOutlier[i])
            {
                mCurrentFrame.mvpMapPoints[i]->IncreaseFound();
                if(!mbOnlyTracking)
                {
                    if(mCurrentFrame.mvpMapPoints[i]->Observations()>0)
                        mnMatchesInliers++;
                }
                else
                    mnMatchesInliers++;
            }
            else if(mSensor==System::STEREO)
                mCurrentFrame.mvpMapPoints[i] = static_cast<MapPoint*>(NULL);

        }
    }

    // Decide if the tracking was succesful
    // More restrictive if there was a relocalization recently
    if(mCurrentFrame.mnId<mnLastRelocFrameId+mMaxFrames && mnMatchesInliers<50)
        return false;

    if(mnMatchesInliers<30)
        return false;
    else
        return true;
}

//点线特征同时使用的当前图像帧局部地图跟踪
bool Tracking::TrackLocalMapBoth()
{

    thread threadUpdateLocalMap(&Tracking::UpdateLocalMap, this );
    thread threadUpdateLocalMapLines(&Tracking::UpdateLocalMapLines, this );
    threadUpdateLocalMap.join();
    threadUpdateLocalMapLines.join();

    int addtionPointsNum;
    int addtionLinesNum;
    thread threadSearchLocalPoints(&Tracking::SearchLocalPoints2, this, ref(addtionPointsNum) );
    thread threadSearchLocalLines(&Tracking::SearchLocalLines, this, ref(addtionLinesNum) );
    threadSearchLocalPoints.join();
    threadSearchLocalLines.join();

    int addtionLinesNum2 = 0;
    //用于在地图线匹配数量较少或者没有时，通过点共视图对地图线信息进行补充
    if(mnMatchesInlierLines+addtionLinesNum <= 16)
    {
	addtionLinesNum2 = MapLineRenewing();
    }
    
    int inlinersNum;
    int inlinerLinesNum;
    if( (mnMatchesInliers+addtionPointsNum>27) &&
	(mnMatchesInlierLines+addtionLinesNum+addtionLinesNum2>16) )
    {	
	Optimizer::PoseOptimizationmain(&mCurrentFrame, inlinersNum, inlinerLinesNum);
	//Optimizer::PoseOptimizationLowFeature(&mCurrentFrame, inlinersNum, inlinerLinesNum);
	cout << "局部地图--正常跟踪成功！" << endl;
    }
    else if( mnMatchesInliers+addtionPointsNum>36 )
    {
	Optimizer::PoseOptimizationDoublePoints(&mCurrentFrame,inlinersNum);
	Optimizer::SetOutlierLinesForPose(&mCurrentFrame,inlinerLinesNum);
	cout << "局部地图--仅点跟踪成功！" << endl;
    }
    else if( (mnMatchesInliers+addtionPointsNum>6) &&
	(mnMatchesInlierLines+addtionLinesNum+addtionLinesNum2>6) )
    {
	Optimizer::PoseOptimizationLowFeature(&mCurrentFrame, inlinersNum, inlinerLinesNum);
	cout << "局部地图--低纹理跟踪成功！" << endl;
    }
    else
	return false;
       
    mnMatchesInliers = 0;
    // Update MapPoints Statistics
    for(int i=0; i<mCurrentFrame.N; i++)
    {
        if(mCurrentFrame.mvpMapPoints[i])
        {
            if(!mCurrentFrame.mvbOutlier[i])
            {
                mCurrentFrame.mvpMapPoints[i]->IncreaseFound();
                if(!mbOnlyTracking)
                {
                    if(mCurrentFrame.mvpMapPoints[i]->Observations()>0)
                        mnMatchesInliers++;
                }
                else
                    mnMatchesInliers++;
            }

        }
    }

    mnMatchesInlierLines = 0;
    // Update MapPoints Statistics
    for(int i=0; i<mCurrentFrame.NL; i++)
    {
        if(mCurrentFrame.mvpMapLines[i])
        {
            if(!mCurrentFrame.mvbOutlierLines[i])
            {
                mCurrentFrame.mvpMapLines[i]->IncreaseFound();		
                if(!mbOnlyTracking)
                {
                    if(mCurrentFrame.mvpMapLines[i]->Observations()>0)
                        mnMatchesInlierLines++;
                }
                else
                    mnMatchesInlierLines++;
            }

        }
    }  
    
    //修改重定位线的跟踪条件
    if(mCurrentFrame.mnId<mnLastRelocFrameId+mMaxFrames && mnMatchesInliers<30 && mnMatchesInlierLines<15 )
        return false;

    if( (mnMatchesInliers>=21 && mnMatchesInlierLines>=12) )
	return true;
    else if(mnMatchesInliers>=28)
	return true; 
    else if(mnMatchesInlierLines+mnMatchesInliers>=12)
	return true; 
    else
	return false;
}

//用于在地图线匹配数量较少或者没有时，通过点共视图对地图线信息进行补充
int Tracking::MapLineRenewing()
{
    int addtionLinesNum = 0;
  
    //获取上一关键帧的所有点共视图
    mvpLocalKeyFrames.clear();
    mvpLocalKeyFrames = mpLastKeyFrame->GetBestCovisibilityKeyFrames(15);
    
    //获取其中所有的地图线
    mvpLocalMapLines.clear();
    for(vector<KeyFrame*>::const_iterator itKF=mvpLocalKeyFrames.begin(), itEndKF=mvpLocalKeyFrames.end(); itKF!=itEndKF; itKF++)
    {
        KeyFrame* pKF = *itKF;
        const vector<MapLine*> vpMLs = pKF->GetMapLineMatches();

        for(vector<MapLine*>::const_iterator itML=vpMLs.begin(), itEndML=vpMLs.end(); itML!=itEndML; itML++)
        {
            MapLine* pML = *itML;
            if(!pML)
                continue;
	    //排除SearchLocalLines()中的地图线
            if(pML->mnTrackReferenceForFrame==mCurrentFrame.mnId)
                continue;
	    //排除当前图像帧中的地图线
	    if(pML->mnLastFrameSeen == mCurrentFrame.mnId)
            continue;
            if(!pML->isBad())
            {
                mvpLocalMapLines.push_back(pML);
                pML->mnTrackReferenceForFrame=mCurrentFrame.mnId;
            }
        }
    }
    
    //局部地图线中，对不为当前图像帧已经找到的地图线进行视野判断
    int nToMatch=0;
    // Project points in frame and check its visibility
    for(vector<MapLine*>::iterator vit=mvpLocalMapLines.begin(), vend=mvpLocalMapLines.end(); vit!=vend; vit++)
    {
        MapLine* pML = *vit;
        
        if(pML->isBad())
            continue;

        if(mCurrentFrame.isInFrustumLine(pML,0))
        {
	    pML->IncreaseVisible();
            nToMatch++;
        }
    }
    
    if(nToMatch>0)
    {
	//较运动模型，适当放宽长度容许误差
        Linematcher matcher(0.9,false,true,0.2);
        int th = 5;
        // If the camera has been relocalised recently, perform a coarser search
        if(mCurrentFrame.mnId<mnLastRelocFrameId+2)
            th = 8;
	
        addtionLinesNum = matcher.SearchByProjection(mCurrentFrame,mvpLocalMapLines,th);
	
	//cout <<"点共视图补充地图线数量:" << addtionLinesNum << endl <<endl <<endl;
    }
    
    return addtionLinesNum;
    
}

bool Tracking::NeedNewKeyFrame()
{
    if(mbOnlyTracking)
        return false;

    // If Local Mapping is freezed by a Loop Closure do not insert keyframes
    if(mpLocalMapper->isStopped() || mpLocalMapper->stopRequested())
        return false;

    const int nKFs = mpMap->KeyFramesInMap();

    // Do not insert keyframes if not enough frames have passed from last relocalisation
    if(mCurrentFrame.mnId<mnLastRelocFrameId+mMaxFrames && nKFs>mMaxFrames)
        return false;

    // Tracked MapPoints in the reference keyframe
    int nMinObs = 3;
    if(nKFs<=2)
        nMinObs=2;
    int nRefMatches = mpReferenceKF->TrackedMapPoints(nMinObs);

    // Local Mapping accept keyframes?
    bool bLocalMappingIdle = mpLocalMapper->AcceptKeyFrames();

    // Check how many "close" points are being tracked and how many could be potentially created.
    int nNonTrackedClose = 0;
    int nTrackedClose= 0;
    if(mSensor!=System::MONOCULAR)
    {
        for(int i =0; i<mCurrentFrame.N; i++)
        {
            if(mCurrentFrame.mvDepth[i]>0 && mCurrentFrame.mvDepth[i]<mThDepth)
            {
                if(mCurrentFrame.mvpMapPoints[i] && !mCurrentFrame.mvbOutlier[i])
                    nTrackedClose++;
                else
                    nNonTrackedClose++;
            }
        }
    }

    bool bNeedToInsertClose = (nTrackedClose<100) && (nNonTrackedClose>70);

    // Thresholds
    float thRefRatio = 0.75f;
    if(nKFs<2)
        thRefRatio = 0.4f;

    if(mSensor==System::MONOCULAR)
        thRefRatio = 0.9f;

    // Condition 1a: More than "MaxFrames" have passed from last keyframe insertion
    const bool c1a = mCurrentFrame.mnId>=mnLastKeyFrameId+mMaxFrames;
    // Condition 1b: More than "MinFrames" have passed and Local Mapping is idle
    const bool c1b = (mCurrentFrame.mnId>=mnLastKeyFrameId+mMinFrames && bLocalMappingIdle);
    //Condition 1c: tracking is weak
    const bool c1c =  mSensor!=System::MONOCULAR && (mnMatchesInliers<nRefMatches*0.25 || bNeedToInsertClose) ;
    // Condition 2: Few tracked points compared to reference keyframe. Lots of visual odometry compared to map matches.
    const bool c2 = ((mnMatchesInliers<nRefMatches*thRefRatio|| bNeedToInsertClose) && mnMatchesInliers>15);

    if((c1a||c1b||c1c)&&c2)
    {
        // If the mapping accepts keyframes, insert keyframe.
        // Otherwise send a signal to interrupt BA
        if(bLocalMappingIdle)
        {
            return true;
        }
        else
        {
            mpLocalMapper->InterruptBA();
            if(mSensor!=System::MONOCULAR)
            {
                if(mpLocalMapper->KeyframesInQueue()<3)
                    return true;
                else
                    return false;
            }
            else
                return false;
        }
    }
    else
        return false;
}

//单目基于点线特征判断是否创建新的关键帧
bool Tracking::NeedNewKeyFrameBoth()
{
    if(mbOnlyTracking)
        return false;

    // If Local Mapping is freezed by a Loop Closure do not insert keyframes
    if(mpLocalMapper->isStopped() || mpLocalMapper->stopRequested())
        return false;

    const int nKFs = mpMap->KeyFramesInMap();

    // Do not insert keyframes if not enough frames have passed from last relocalisation
    if(mCurrentFrame.mnId<mnLastRelocFrameId+mMaxFrames && nKFs>mMaxFrames)
        return false;

    // Tracked MapPoints in the reference keyframe
    int nMinObs = 3;
    if(nKFs<=2)
        nMinObs=2;
    int nRefMatches = mpReferenceKF->TrackedMapPoints(nMinObs);//获取所有当前图像帧可以观测到且不为坏的地图点
    int nRefMatchesLines = mpReferenceKFLines->TrackedMapLines(nMinObs);

    // Local Mapping accept keyframes?
    bool bLocalMappingIdle = mpLocalMapper->AcceptKeyFrames();

    float thRefRatio = 0.9f;
    float thRefRatioLines = 0.8f;
    
    // Condition 1a: More than "MaxFrames" have passed from last keyframe insertion
    const bool c1a = mCurrentFrame.mnId>=mnLastKeyFrameId+mMaxFrames;
    
    // Condition 1b: More than "MinFrames" have passed and Local Mapping is idle
    const bool c1b = (mCurrentFrame.mnId>=mnLastKeyFrameId+mMinFrames && bLocalMappingIdle);
    
    //Condition 1c: tracking is weak

    const bool c1c =  mSensor!=System::MONOCULAR && (mnMatchesInliers<nRefMatches*0.25 ) ;
    
    // Condition 2: Few tracked points compared to reference keyframe. Lots of visual odometry compared to map matches.
    const bool c2 = (( (mnMatchesInliers<nRefMatches*thRefRatio) || (mnMatchesInlierLines<nRefMatchesLines*thRefRatioLines)) && 
			(mnMatchesInliers>15 ||  mnMatchesInlierLines>10 || mnMatchesInliers+mnMatchesInlierLines>=12 ) );

    if((c1a||c1b||c1c)&&c2)
    {
        // If the mapping accepts keyframes, insert keyframe.
        // Otherwise send a signal to interrupt BA
        if(bLocalMappingIdle)
        {
            return true;
        }
        else
        {	
            mpLocalMapper->InterruptBA();
	    
            if(mSensor!=System::MONOCULAR)
            {
                if(mpLocalMapper->KeyframesInQueue()<3)
                    return true;
                else
                    return false;
            }
            else
                return false;
        }
    }
    else
        return false;
}

void Tracking::CreateNewKeyFrame()
{
    if(!mpLocalMapper->SetNotStop(true))
        return;

    KeyFrame* pKF = new KeyFrame(mCurrentFrame,mpMap,mpKeyFrameDB);

    mpReferenceKF = pKF;
    mCurrentFrame.mpReferenceKF = pKF;

    if(mSensor!=System::MONOCULAR)
    {
        mCurrentFrame.UpdatePoseMatrices();

        // We sort points by the measured depth by the stereo/RGBD sensor.
        // We create all those MapPoints whose depth < mThDepth.
        // If there are less than 100 close points we create the 100 closest.
        vector<pair<float,int> > vDepthIdx;
        vDepthIdx.reserve(mCurrentFrame.N);
        for(int i=0; i<mCurrentFrame.N; i++)
        {
            float z = mCurrentFrame.mvDepth[i];
            if(z>0)
            {
                vDepthIdx.push_back(make_pair(z,i));
            }
        }

        if(!vDepthIdx.empty())
        {
            sort(vDepthIdx.begin(),vDepthIdx.end());

            int nPoints = 0;
            for(size_t j=0; j<vDepthIdx.size();j++)
            {
                int i = vDepthIdx[j].second;

                bool bCreateNew = false;

                MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];
                if(!pMP)
                    bCreateNew = true;
                else if(pMP->Observations()<1)
                {
                    bCreateNew = true;
                    mCurrentFrame.mvpMapPoints[i] = static_cast<MapPoint*>(NULL);
                }

                if(bCreateNew)
                {
                    cv::Mat x3D = mCurrentFrame.UnprojectStereo(i);
                    MapPoint* pNewMP = new MapPoint(x3D,pKF,mpMap);
                    pNewMP->AddObservation(pKF,i);
                    pKF->AddMapPoint(pNewMP,i);
                    pNewMP->ComputeDistinctiveDescriptors();
                    pNewMP->UpdateNormalAndDepth();
                    mpMap->AddMapPoint(pNewMP);

                    mCurrentFrame.mvpMapPoints[i]=pNewMP;
                    nPoints++;
                }
                else
                {
                    nPoints++;
                }

                if(vDepthIdx[j].first>mThDepth && nPoints>100)
                    break;
            }
        }
    }

    mpLocalMapper->InsertKeyFrame(pKF);

    mpLocalMapper->SetNotStop(false);

    mnLastKeyFrameId = mCurrentFrame.mnId;
    mpLastKeyFrame = pKF;
}

//单目点线特征创建新的关键帧
void Tracking::CreateNewKeyFrameBoth()
{
    if(!mpLocalMapper->SetNotStop(true))
        return;

    KeyFrame* pKF = new KeyFrame(mCurrentFrame,mpMap,mpKeyFrameDB);

    mpReferenceKF = pKF;
    mpReferenceKFLines = pKF;//更新当前线程的线参考关键帧
    mCurrentFrame.mpReferenceKF = pKF;
    mCurrentFrame.mpReferenceKFLines = pKF;

    mpLocalMapper->InsertKeyFrame(pKF);

    mpLocalMapper->SetNotStop(false);

    mnLastKeyFrameId = mCurrentFrame.mnId;
    mpLastKeyFrame = pKF;
}

void Tracking::SearchLocalPoints()
{
    // Do not search map points already matched
    for(vector<MapPoint*>::iterator vit=mCurrentFrame.mvpMapPoints.begin(), vend=mCurrentFrame.mvpMapPoints.end(); vit!=vend; vit++)
    {
        MapPoint* pMP = *vit;
        if(pMP)
        {
            if(pMP->isBad())
            {
                *vit = static_cast<MapPoint*>(NULL);
            }
            else
            {
                pMP->IncreaseVisible();
                pMP->mnLastFrameSeen = mCurrentFrame.mnId;
                pMP->mbTrackInView = false;
            }
        }
    }

    int nToMatch=0;

    // Project points in frame and check its visibility
    for(vector<MapPoint*>::iterator vit=mvpLocalMapPoints.begin(), vend=mvpLocalMapPoints.end(); vit!=vend; vit++)
    {
        MapPoint* pMP = *vit;
        if(pMP->mnLastFrameSeen == mCurrentFrame.mnId)
            continue;
        if(pMP->isBad())
            continue;
        // Project (this fills MapPoint variables for matching)
        if(mCurrentFrame.isInFrustum(pMP,0.5))
        {
            pMP->IncreaseVisible();
            nToMatch++;
        }
    }

    if(nToMatch>0)
    {
        ORBmatcher matcher(0.8);
        int th = 1;
        if(mSensor==System::RGBD)
            th=3;
        // If the camera has been relocalised recently, perform a coarser search
        if(mCurrentFrame.mnId<mnLastRelocFrameId+2)
            th=5;
        matcher.SearchByProjection(mCurrentFrame,mvpLocalMapPoints,th);
    }
}

//点线特征框架需要返回额外找到的地图点数
void Tracking::SearchLocalPoints2(int &addtionPointsNum)
{
    // Do not search map points already matched
    for(vector<MapPoint*>::iterator vit=mCurrentFrame.mvpMapPoints.begin(), vend=mCurrentFrame.mvpMapPoints.end(); vit!=vend; vit++)
    {
        MapPoint* pMP = *vit;
        if(pMP)
        {
            if(pMP->isBad())
            {
                *vit = static_cast<MapPoint*>(NULL);
            }
            else
            {
                pMP->IncreaseVisible();
                pMP->mnLastFrameSeen = mCurrentFrame.mnId;
                pMP->mbTrackInView = false;
            }
        }
    }

    int nToMatch=0;

    // Project points in frame and check its visibility
    for(vector<MapPoint*>::iterator vit=mvpLocalMapPoints.begin(), vend=mvpLocalMapPoints.end(); vit!=vend; vit++)
    {
        MapPoint* pMP = *vit;
        if(pMP->mnLastFrameSeen == mCurrentFrame.mnId)
            continue;
        if(pMP->isBad())
            continue;
        // Project (this fills MapPoint variables for matching)
        if(mCurrentFrame.isInFrustum(pMP,0.5))
        {
            pMP->IncreaseVisible();
            nToMatch++;
        }
    }

    if(nToMatch>0)
    {
        ORBmatcher matcher(0.8);
        int th = 1;
        if(mSensor==System::RGBD)
            th=3;
        // If the camera has been relocalised recently, perform a coarser search
        if(mCurrentFrame.mnId<mnLastRelocFrameId+2)
            th=5;
        addtionPointsNum = matcher.SearchByProjection(mCurrentFrame,mvpLocalMapPoints,th);
    }
}

//将局部地图线在当前图像帧上进行投影，获得更多地图线匹配
void Tracking::SearchLocalLines(int &addtionLinesNum)
{
    // Do not search map Lines already matched
    for(vector<MapLine*>::iterator vit=mCurrentFrame.mvpMapLines.begin(), vend=mCurrentFrame.mvpMapLines.end(); vit!=vend; vit++)
    {
        MapLine* pML = *vit;
        if(pML)
        {
            if(pML->isBad())
            {
                *vit = static_cast<MapLine*>(NULL);
            }
            else
            {
                pML->IncreaseVisible();
                pML->mnLastFrameSeen = mCurrentFrame.mnId;
                pML->mbTrackInView = false;
            }
        }
    }

    int nToMatch=0;

    // Project points in frame and check its visibility
    for(vector<MapLine*>::iterator vit=mvpLocalMapLines.begin(), vend=mvpLocalMapLines.end(); vit!=vend; vit++)
    {
        MapLine* pML = *vit;
        if(pML->mnLastFrameSeen == mCurrentFrame.mnId)
            continue;
        if(pML->isBad())
            continue;

        if(mCurrentFrame.isInFrustumLine(pML,0))
        {
            pML->IncreaseVisible();
            nToMatch++;
        }
    }
    
    if(nToMatch>0)
    {
	//较运动模型，适当放宽长度容许误差
        Linematcher matcher(0.9,false,false,0.2);
        int th = 10;
        // If the camera has been relocalised recently, perform a coarser search
        if(mCurrentFrame.mnId<mnLastRelocFrameId+2)
            th = 15;
        addtionLinesNum = matcher.SearchByProjection(mCurrentFrame,mvpLocalMapLines,th);
    }
}

void Tracking::UpdateLocalMap()
{
    // This is for visualization
    mpMap->SetReferenceMapPoints(mvpLocalMapPoints);

    // Update
    UpdateLocalKeyFrames();
    UpdateLocalPoints();
}

void Tracking::UpdateLocalMapLines()
{
    // This is for visualization
    mpMap->SetReferenceMapLines(mvpLocalMapLines);

    // Update
    UpdateLocalKeyFramesLines();
    UpdateLocalLines();
}

void Tracking::UpdateLocalPoints()
{
    mvpLocalMapPoints.clear();

    for(vector<KeyFrame*>::const_iterator itKF=mvpLocalKeyFrames.begin(), itEndKF=mvpLocalKeyFrames.end(); itKF!=itEndKF; itKF++)
    {
        KeyFrame* pKF = *itKF;
        const vector<MapPoint*> vpMPs = pKF->GetMapPointMatches();

        for(vector<MapPoint*>::const_iterator itMP=vpMPs.begin(), itEndMP=vpMPs.end(); itMP!=itEndMP; itMP++)
        {
            MapPoint* pMP = *itMP;
            if(!pMP)
                continue;
            if(pMP->mnTrackReferenceForFrame==mCurrentFrame.mnId)
                continue;
            if(!pMP->isBad())
            {
                mvpLocalMapPoints.push_back(pMP);
                pMP->mnTrackReferenceForFrame=mCurrentFrame.mnId;
            }
        }
    }
}

//更新局部地图跟踪中的局部地图线集
void Tracking::UpdateLocalLines()
{
    mvpLocalMapLines.clear();

    for(vector<KeyFrame*>::const_iterator itKF=mvpLocalKeyFramesLines.begin(), itEndKF=mvpLocalKeyFramesLines.end(); itKF!=itEndKF; itKF++)
    {
        KeyFrame* pKF = *itKF;
        const vector<MapLine*> vpMLs = pKF->GetMapLineMatches();

        for(vector<MapLine*>::const_iterator itML=vpMLs.begin(), itEndML=vpMLs.end(); itML!=itEndML; itML++)
        {
            MapLine* pML = *itML;
            if(!pML)
                continue;
            if(pML->mnTrackReferenceForFrame==mCurrentFrame.mnId)
                continue;
            if(!pML->isBad())
            {
                mvpLocalMapLines.push_back(pML);
                pML->mnTrackReferenceForFrame=mCurrentFrame.mnId;
            }
        }
    }
}

void Tracking::UpdateLocalKeyFrames()
{
    // Each map point vote for the keyframes in which it has been observed
    map<KeyFrame*,int> keyframeCounter;
    for(int i=0; i<mCurrentFrame.N; i++)
    {
        if(mCurrentFrame.mvpMapPoints[i])
        {
            MapPoint* pMP = mCurrentFrame.mvpMapPoints[i];
            if(!pMP->isBad())
            {
                const map<KeyFrame*,size_t> observations = pMP->GetObservations();
                for(map<KeyFrame*,size_t>::const_iterator it=observations.begin(), itend=observations.end(); it!=itend; it++)
                    keyframeCounter[it->first]++;
            }
            else
            {
                mCurrentFrame.mvpMapPoints[i]=NULL;
            }
        }
    }

    if(keyframeCounter.empty())
        return;

    int max=0;
    KeyFrame* pKFmax= static_cast<KeyFrame*>(NULL);

    mvpLocalKeyFrames.clear();
    mvpLocalKeyFrames.reserve(3*keyframeCounter.size());

    // All keyframes that observe a map point are included in the local map. Also check which keyframe shares most points
    for(map<KeyFrame*,int>::const_iterator it=keyframeCounter.begin(), itEnd=keyframeCounter.end(); it!=itEnd; it++)
    {
        KeyFrame* pKF = it->first;

        if(pKF->isBad())
            continue;

        if(it->second>max)
        {
            max=it->second;
            pKFmax=pKF;
        }

        mvpLocalKeyFrames.push_back(it->first);
        pKF->mnTrackReferenceForFrame = mCurrentFrame.mnId;
    }


    // Include also some not-already-included keyframes that are neighbors to already-included keyframes
    for(vector<KeyFrame*>::const_iterator itKF=mvpLocalKeyFrames.begin(), itEndKF=mvpLocalKeyFrames.end(); itKF!=itEndKF; itKF++)
    {
        // Limit the number of keyframes
	if(mbusingLine)
	{
	    if(mvpLocalKeyFrames.size()>50)
	      break;
	}
	else
	{
	    if(mvpLocalKeyFrames.size()>80)
	      break;
	}

        KeyFrame* pKF = *itKF;

        const vector<KeyFrame*> vNeighs = pKF->GetBestCovisibilityKeyFrames(10);

        for(vector<KeyFrame*>::const_iterator itNeighKF=vNeighs.begin(), itEndNeighKF=vNeighs.end(); itNeighKF!=itEndNeighKF; itNeighKF++)
        {
            KeyFrame* pNeighKF = *itNeighKF;
            if(!pNeighKF->isBad())
            {
                if(pNeighKF->mnTrackReferenceForFrame!=mCurrentFrame.mnId)
                {
                    mvpLocalKeyFrames.push_back(pNeighKF);
                    pNeighKF->mnTrackReferenceForFrame=mCurrentFrame.mnId;
                    break;
                }
            }
        }

        const set<KeyFrame*> spChilds = pKF->GetChilds();
        for(set<KeyFrame*>::const_iterator sit=spChilds.begin(), send=spChilds.end(); sit!=send; sit++)
        {
            KeyFrame* pChildKF = *sit;
            if(!pChildKF->isBad())
            {
                if(pChildKF->mnTrackReferenceForFrame!=mCurrentFrame.mnId)
                {
                    mvpLocalKeyFrames.push_back(pChildKF);
                    pChildKF->mnTrackReferenceForFrame=mCurrentFrame.mnId;
                    break;
                }
            }
        }

        KeyFrame* pParent = pKF->GetParent();
        if(pParent)
        {
            if(pParent->mnTrackReferenceForFrame!=mCurrentFrame.mnId)
            {
                mvpLocalKeyFrames.push_back(pParent);
                pParent->mnTrackReferenceForFrame=mCurrentFrame.mnId;
                break;
            }
        }

    }

    if(pKFmax)
    {
        mpReferenceKF = pKFmax;
        mCurrentFrame.mpReferenceKF = mpReferenceKF;
    }
}

//当前图像帧的线特征局部关键帧更新
void Tracking::UpdateLocalKeyFramesLines()
{
    // Each map Line vote for the keyframes in which it has been observed
    map<KeyFrame*,int> keyframeCounter;
    for(int i=0; i<mCurrentFrame.NL; i++)
    {
        if(mCurrentFrame.mvpMapLines[i])
        {
            MapLine* pML = mCurrentFrame.mvpMapLines[i];
            if(!pML->isBad())
            {
                const map<KeyFrame*,size_t> observations = pML->GetObservations();
                for(map<KeyFrame*,size_t>::const_iterator it=observations.begin(), itend=observations.end(); it!=itend; it++)
                    keyframeCounter[it->first]++;
            }
            else
            {
                mCurrentFrame.mvpMapLines[i]=NULL;
            }
        }
    }
    
    if(keyframeCounter.empty())
        return;

    int max=0;
    KeyFrame* pKFmax= static_cast<KeyFrame*>(NULL);

    mvpLocalKeyFramesLines.clear();
    mvpLocalKeyFramesLines.reserve(3*keyframeCounter.size());

    // All keyframes that observe a map point are included in the local map. Also check which keyframe shares most points
    for(map<KeyFrame*,int>::const_iterator it=keyframeCounter.begin(), itEnd=keyframeCounter.end(); it!=itEnd; it++)
    {
        KeyFrame* pKF = it->first;

        if(pKF->isBadLines())
            continue;

        if(it->second>max)
        {
            max=it->second;
            pKFmax=pKF;
        }

        mvpLocalKeyFramesLines.push_back(it->first);
        pKF->mnTrackReferenceForFrameLines = mCurrentFrame.mnId;
    }


    // Include also some not-already-included keyframes that are neighbors to already-included keyframes
    for(vector<KeyFrame*>::const_iterator itKF=mvpLocalKeyFramesLines.begin(), itEndKF=mvpLocalKeyFramesLines.end(); itKF!=itEndKF; itKF++)
    {
        // Limit the number of keyframes
        if(mvpLocalKeyFramesLines.size()>80)
            break;

        KeyFrame* pKF = *itKF;

        const vector<KeyFrame*> vNeighs = pKF->GetBestCovisibilityKeyFramesLines(10);

        for(vector<KeyFrame*>::const_iterator itNeighKF=vNeighs.begin(), itEndNeighKF=vNeighs.end(); itNeighKF!=itEndNeighKF; itNeighKF++)
        {
            KeyFrame* pNeighKF = *itNeighKF;
            if(!pNeighKF->isBadLines())
            {
                if(pNeighKF->mnTrackReferenceForFrameLines!=mCurrentFrame.mnId)
                {
                    mvpLocalKeyFramesLines.push_back(pNeighKF);
                    pNeighKF->mnTrackReferenceForFrameLines=mCurrentFrame.mnId;
                    break;
                }
            }
        }

        const set<KeyFrame*> spChilds = pKF->GetChildsLines();
        for(set<KeyFrame*>::const_iterator sit=spChilds.begin(), send=spChilds.end(); sit!=send; sit++)
        {
            KeyFrame* pChildKF = *sit;
            if(!pChildKF->isBadLines())
            {
                if(pChildKF->mnTrackReferenceForFrameLines!=mCurrentFrame.mnId)
                {
                    mvpLocalKeyFramesLines.push_back(pChildKF);
                    pChildKF->mnTrackReferenceForFrameLines=mCurrentFrame.mnId;
                    break;
                }
            }
        }

        KeyFrame* pParent = pKF->GetParentLines();
        if(pParent)
        {
            if(pParent->mnTrackReferenceForFrameLines!=mCurrentFrame.mnId)
            {
                mvpLocalKeyFramesLines.push_back(pParent);
                pParent->mnTrackReferenceForFrameLines=mCurrentFrame.mnId;
                break;
            }
        }

    }

    if(pKFmax)
    {
        mpReferenceKFLines = pKFmax;
        mCurrentFrame.mpReferenceKFLines = mpReferenceKFLines;
    }
}

bool Tracking::Relocalization()
{
    // Compute Bag of Words Vector
    mCurrentFrame.ComputeBoW();

    // Relocalization is performed when tracking is lost
    // Track Lost: Query KeyFrame Database for keyframe candidates for relocalisation
    vector<KeyFrame*> vpCandidateKFs = mpKeyFrameDB->DetectRelocalizationCandidates(&mCurrentFrame);

    if(vpCandidateKFs.empty())
        return false;

    const int nKFs = vpCandidateKFs.size();

    // We perform first an ORB matching with each candidate
    // If enough matches are found we setup a PnP solver
    ORBmatcher matcher(0.75,true);

    vector<PnPsolver*> vpPnPsolvers;
    vpPnPsolvers.resize(nKFs);

    vector<vector<MapPoint*> > vvpMapPointMatches;
    vvpMapPointMatches.resize(nKFs);

    vector<bool> vbDiscarded;
    vbDiscarded.resize(nKFs);

    int nCandidates=0;

    for(int i=0; i<nKFs; i++)
    {
        KeyFrame* pKF = vpCandidateKFs[i];
        if(pKF->isBad())
            vbDiscarded[i] = true;
        else
        {
            int nmatches = matcher.SearchByBoW(pKF,mCurrentFrame,vvpMapPointMatches[i]);
            if(nmatches<15)
            {
                vbDiscarded[i] = true;
                continue;
            }
            else
            {
                PnPsolver* pSolver = new PnPsolver(mCurrentFrame,vvpMapPointMatches[i]);
                pSolver->SetRansacParameters(0.99,10,300,4,0.5,5.991);
                vpPnPsolvers[i] = pSolver;
                nCandidates++;
            }
        }
    }

    // Alternatively perform some iterations of P4P RANSAC
    // Until we found a camera pose supported by enough inliers
    bool bMatch = false;
    ORBmatcher matcher2(0.9,true);

    while(nCandidates>0 && !bMatch)
    {
        for(int i=0; i<nKFs; i++)
        {
            if(vbDiscarded[i])
                continue;

            // Perform 5 Ransac Iterations
            vector<bool> vbInliers;
            int nInliers;
            bool bNoMore;

            PnPsolver* pSolver = vpPnPsolvers[i];
	    
	    Timer timer;
	    timer.start();
            cv::Mat Tcw = pSolver->iterate(5,bNoMore,vbInliers,nInliers);
	    double pSolverTime = timer.stop();   
	    cout << "EPnP算法迭代计算初步位姿时间: " << pSolverTime << " ms" <<endl <<endl;
	    
            // If Ransac reachs max. iterations discard keyframe
            if(bNoMore)
            {
                vbDiscarded[i]=true;
                nCandidates--;
            }

            // If a Camera Pose is computed, optimize
            if(!Tcw.empty())
            {
                Tcw.copyTo(mCurrentFrame.mTcw);

                set<MapPoint*> sFound;

                const int np = vbInliers.size();

                for(int j=0; j<np; j++)
                {
                    if(vbInliers[j])
                    {
                        mCurrentFrame.mvpMapPoints[j]=vvpMapPointMatches[i][j];
                        sFound.insert(vvpMapPointMatches[i][j]);
                    }
                    else
                        mCurrentFrame.mvpMapPoints[j]=NULL;
                }

                int nGood = Optimizer::PoseOptimization(&mCurrentFrame);

                if(nGood<10)
                    continue;

                for(int io =0; io<mCurrentFrame.N; io++)
                    if(mCurrentFrame.mvbOutlier[io])
                        mCurrentFrame.mvpMapPoints[io]=static_cast<MapPoint*>(NULL);

                // If few inliers, search by projection in a coarse window and optimize again
                if(nGood<50)
                {
                    int nadditional =matcher2.SearchByProjection(mCurrentFrame,vpCandidateKFs[i],sFound,10,100);

                    if(nadditional+nGood>=50)
                    {
                        nGood = Optimizer::PoseOptimization(&mCurrentFrame);

                        // If many inliers but still not enough, search by projection again in a narrower window
                        // the camera has been already optimized with many points
                        if(nGood>30 && nGood<50)
                        {
                            sFound.clear();
                            for(int ip =0; ip<mCurrentFrame.N; ip++)
                                if(mCurrentFrame.mvpMapPoints[ip])
                                    sFound.insert(mCurrentFrame.mvpMapPoints[ip]);
                            nadditional =matcher2.SearchByProjection(mCurrentFrame,vpCandidateKFs[i],sFound,3,64);

                            // Final optimization
                            if(nGood+nadditional>=50)
                            {
                                nGood = Optimizer::PoseOptimization(&mCurrentFrame);

                                for(int io =0; io<mCurrentFrame.N; io++)
                                    if(mCurrentFrame.mvbOutlier[io])
                                        mCurrentFrame.mvpMapPoints[io]=NULL;
                            }
                        }
                    }
                }


                // If the pose is supported by enough inliers stop ransacs and continue
                if(nGood>=50)
                {
                    bMatch = true;
                    break;
                }
            }
        }
    }

    if(!bMatch)
    {
        return false;
    }
    else
    {
        mnLastRelocFrameId = mCurrentFrame.mnId;
        return true;
    }

}

//单目线特征的重定位主函数
bool Tracking::RelocalizationBoth()
{
    
    //重定位前，确认存在点线特征信息
    if(mCurrentFrame.mvKeys.empty()|| mCurrentFrame.mvLines.empty() )
	return false;
    
    // step1.1: Compute Bag of Words Vector
    mCurrentFrame.ComputeBoW();

    //  step1.2:Relocalization is performed when tracking is lost
    // step1.2:Track Lost: Query KeyFrame Database for keyframe candidates for relocalisation
    vector<KeyFrame*> vpCandidateKFs = mpKeyFrameDB->DetectRelocalizationCandidates(&mCurrentFrame);

    if(vpCandidateKFs.empty())
        return false;

    const int nKFs = vpCandidateKFs.size();
    
    // If enough matches are found we setup a PnP solver
    //候选关键帧对应的EPnL求解器
    vector<PnPsolver*> vpPnPsolvers;
    vpPnPsolvers.resize(nKFs);

    vector<vector<MapPoint*> > vvpMapPointMatches;
    vvpMapPointMatches.resize(nKFs);
    vector<vector<MapLine*> > vvpMapLineMatches;
    vvpMapLineMatches.resize(nKFs);

    vector<bool> vbDiscarded;
    vbDiscarded.resize(nKFs);
    
    cout << "重定位准备完成！" << endl<< endl<< endl;

    int nCandidates=0;

    for(int i=0; i<nKFs; i++)
    {
        KeyFrame* pKF = vpCandidateKFs[i];

        if(pKF->isBad() || pKF->isBadLines())
            vbDiscarded[i] = true;
        else
        {
	    //通过引用返回的匹配数
	    int nmatches;
	    int nmatchLines;
	    
	    cout << "进入参考关键帧分线程匹配！" << endl<< endl<< endl;

	    //注意：两个线程使用的参考关键帧是相同的
	    thread threadSearchByBoW(&Tracking::SearchByBoW, this, pKF, ref(mCurrentFrame), 
				  ref(vvpMapPointMatches[i]), ref(nmatches) );
	    thread threadSearchByKNNLines(&Tracking::SearchByKNNLines, this, pKF, 
				  ref(mCurrentFrame), ref(vvpMapLineMatches[i]), ref(nmatchLines) );
    
	    threadSearchByBoW.join();
	    threadSearchByKNNLines.join();
	    
	    cout << "参考关键帧分线程匹配成功！" << endl<< endl<< endl;
	    
	    //若最初的两种匹配数均很少，则直接放弃该候选关键帧
            if(nmatches<15 || nmatchLines<8)
            {
                vbDiscarded[i] = true;
                continue;
            }
            else
            {	
		PnPsolver* pSolver = new PnPsolver(mCurrentFrame,vvpMapLineMatches[i]);
		//设置参数
		pSolver->SetRansacParametersLines(0.99,8,300,4,0.5,3.841);
		vpPnPsolvers[i] = pSolver;
		nCandidates++;
            }
        }
    }
    
    // Alternatively perform some iterations of P4P RANSAC
    // Until we found a camera pose supported by enough inliers
    bool bMatch = false;

    while(nCandidates>0 && !bMatch)
    {
        for(int i=0; i<nKFs; i++)
        {
            if(vbDiscarded[i])
                continue;

            // Perform 5 Ransac Iterations
	    vector<bool> vbInlierLines;
            int nInlierLines;
            bool bNoMore;

            PnPsolver* pSolver = vpPnPsolvers[i];

	    Timer timer;
	    timer.start();	    
            cv::Mat Tcw = pSolver->iterateLines(5,bNoMore,vbInlierLines,nInlierLines);

	    double pSolverTime = timer.stop();   
	    cout << "EPnL算法迭代计算初步位姿时间: " << pSolverTime << " ms" <<endl <<endl ;
    
	    if(bNoMore)
	    {
		cout << "EPnL算法迭代计算初步位姿失败！" << endl;
	    }
	    else
	    {
		cout << "EPnL算法迭代计算初步位姿成功！" << endl << endl;
	    }
	    
            // If Ransac reachs max. iterations discard keyframe
            if(bNoMore)
            {
		//当前候选关键帧失败
                vbDiscarded[i]=true;
                nCandidates--;
            }

            // If a Camera Pose is computed, optimize
            if(!Tcw.empty())
            {
		//找到初步的位姿估计，并拷贝到当前图像帧
                Tcw.copyTo(mCurrentFrame.mTcw);

                set<MapPoint*> sFound;
		set<MapLine*>  sFoundLines;

		const int nl = vbInlierLines.size();
		
		//补充函数：基于当前帧位姿，通过重投影误差剔除无效的地图点匹配，
		int nInliers = SetCurrentFrameMappointsAndInliers(mCurrentFrame, vvpMapPointMatches[i],sFound);
		
		cout << "nInliers: " << nInliers << endl;
		
		//说明当前计算位姿过于偏向线特征
		if( nInliers <15 )
		    return false;
		
                //设置当前图像帧的地图线向量和找到线向量
                for(int j=0; j<nl; j++)
                {
                    if(vbInlierLines[j])
                    {
                        mCurrentFrame.mvpMapLines[j]=vvpMapLineMatches[i][j];
                        sFoundLines.insert(vvpMapLineMatches[i][j]);
                    }
                    else
                        mCurrentFrame.mvpMapLines[j]=NULL;
                }
                  
                  
                //return false;
		
		int inlinersNum;
		int inlinerLinesNum;
		    Optimizer::PoseOptimizationmain(&mCurrentFrame, inlinersNum, inlinerLinesNum);
					
		//优化后，发现内点和内线均太少，即当前计算位姿并不可靠
                if( inlinersNum<10 || inlinerLinesNum<8)
		    return false;
		
		//将外点和外线置为空
                for(int io =0; io<mCurrentFrame.N; io++)
                    if(mCurrentFrame.mvbOutlier[io])
                        mCurrentFrame.mvpMapPoints[io]=static_cast<MapPoint*>(NULL);
		for(int io =0; io<mCurrentFrame.NL; io++)
                    if(mCurrentFrame.mvbOutlierLines[io])
                        mCurrentFrame.mvpMapLines[io]=static_cast<MapLine*>(NULL);

                // If few inliers, search by projection in a coarse window and optimize again
                if( inlinersNum<40 && inlinerLinesNum<20)
                { 
	
		    //通过引用返回的匹配数
		    int nadditional;
		    int nadditionalLines;

		    //注意：两个线程使用的参考关键帧是相同的
		    thread threadSearchByProjection(&Tracking::SearchByProjection, this, ref(mCurrentFrame), vpCandidateKFs[i], 
				  ref(sFound),10,100, ref(nadditional) );
		    thread threadSearchByProjectionLines(&Tracking::SearchByProjectionLines, this,  ref(mCurrentFrame),
				  vpCandidateKFs[i], ref(sFoundLines),30,100, ref(nadditionalLines) );
    
		    threadSearchByProjection.join();
		    threadSearchByProjectionLines.join();	    
		    
                    if( (nadditionalLines+inlinerLinesNum>=20) && (nadditional+inlinersNum>=40) )
                    {
			Optimizer::PoseOptimizationmain(&mCurrentFrame, inlinersNum, inlinerLinesNum);
			RelocalizationBothTwiceSearch(vpCandidateKFs[i],inlinersNum, inlinerLinesNum);
                    }
                    else
		      return false;
                    	     
                }
                else if(inlinersNum<40 && inlinerLinesNum>=20)
		{
		    ORBmatcher matcher(0.9,true);
		    int nadditional = matcher.SearchByProjection(mCurrentFrame,vpCandidateKFs[i],sFound,10,100);

		    if( nadditional+inlinersNum>=40 )
		    {
			Optimizer::PoseOptimizationmain(&mCurrentFrame, inlinersNum, inlinerLinesNum);
			RelocalizationBothTwiceSearch(vpCandidateKFs[i],inlinersNum, inlinerLinesNum);		
                    }
                    else
		      return false;
		    
		}
		else if(inlinerLinesNum<20 && inlinersNum>=40)
		{
		    
		    Linematcher matcher(0.9,true,false,0.15);
		    int nadditionalLines = matcher.SearchByProjection(mCurrentFrame,vpCandidateKFs[i],sFoundLines,30,100);

		    if( nadditionalLines+inlinerLinesNum>=30)
		    {
			Optimizer::PoseOptimizationmain(&mCurrentFrame, inlinersNum, inlinerLinesNum);
			RelocalizationBothTwiceSearch(vpCandidateKFs[i],inlinersNum, inlinerLinesNum);
                    }
                    else
		      return false;
		    
		}

                if( inlinerLinesNum>=20 && inlinersNum>=40 )
                {
                    bMatch = true;
                    break;
                }
                else
		  return false;
            }
        }
    }

    if(!bMatch)
    {
        return false;
    }
    else
    {
        mnLastRelocFrameId = mCurrentFrame.mnId;
        return true;
    }

}

void Tracking::RelocalizationBothTwiceSearch(KeyFrame* vpCandidateKFs, int &inlinersNum, int &inlinerLinesNum )
{
    int nadditional;
    int nadditionalLines;
    
    std::set<MapPoint*> sFound;
    std::set<MapLine*> sFoundLines;
  
    if( (inlinersNum>30 && inlinersNum<40) && (inlinerLinesNum>10 && inlinerLinesNum<20) )
    {
	//获得已经找到的地图点和地图线
	for(int ip =0; ip<mCurrentFrame.N; ip++)
	    if(mCurrentFrame.mvpMapPoints[ip])
		sFound.insert(mCurrentFrame.mvpMapPoints[ip]);
	    
	for(int ip =0; ip<mCurrentFrame.NL; ip++)
	    if(mCurrentFrame.mvpMapLines[ip])
		sFoundLines.insert(mCurrentFrame.mvpMapLines[ip]);

	//注意：两个线程使用的参考关键帧是相同的
	thread threadSearchByProjection(&Tracking::SearchByProjection, this, ref(mCurrentFrame), vpCandidateKFs, 
		      ref(sFound),3,64, ref(nadditional) );
	thread threadSearchByProjectionLines(&Tracking::SearchByProjectionLines, this,  ref(mCurrentFrame),
		      vpCandidateKFs, ref(sFoundLines),30,64, ref(nadditionalLines) );

	threadSearchByProjection.join();
	threadSearchByProjectionLines.join();
	
	// Final optimization
	if( (nadditionalLines+inlinerLinesNum>=20) && (nadditional+inlinersNum>=40) )
	{
	    Optimizer::PoseOptimizationmain(&mCurrentFrame, inlinersNum, inlinerLinesNum);

	    for(int io =0; io<mCurrentFrame.N; io++)
		if(mCurrentFrame.mvbOutlier[io])
		    mCurrentFrame.mvpMapPoints[io]=NULL;
		
	    for(int io =0; io<mCurrentFrame.NL; io++)
		if(mCurrentFrame.mvbOutlierLines[io])
		    mCurrentFrame.mvpMapLines[io]=NULL;
	}
	  
    }
    else if( (inlinersNum>30 && inlinersNum<40) && inlinerLinesNum>=20 )
    {
	//生成点匹配器
	ORBmatcher matcher(0.9,true);
	
	//获得已经找到的地图点
	for(int ip =0; ip<mCurrentFrame.N; ip++)
	    if(mCurrentFrame.mvpMapPoints[ip])
		sFound.insert(mCurrentFrame.mvpMapPoints[ip]);
	
	//最后一次投影匹配，希望找到足够的新匹配地图点
	nadditional = matcher.SearchByProjection(mCurrentFrame,vpCandidateKFs,sFound,3,64);
	
	// Final optimization
	if( inlinerLinesNum>=20 && (nadditional+inlinersNum>=40) )
	{
	    Optimizer::PoseOptimizationmain(&mCurrentFrame, inlinersNum, inlinerLinesNum);

	    for(int io =0; io<mCurrentFrame.N; io++)
		if(mCurrentFrame.mvbOutlier[io])
		    mCurrentFrame.mvpMapPoints[io]=NULL;
		
	    for(int io =0; io<mCurrentFrame.NL; io++)
		if(mCurrentFrame.mvbOutlierLines[io])
		    mCurrentFrame.mvpMapLines[io]=NULL;
	}
	
    }
    else if( (inlinerLinesNum>10 && inlinerLinesNum<20) && inlinersNum>=40)
    {
	//生成线匹配器
	Linematcher matcher(0.9,true,false,0.25);
	
	//获得已经找到的地图线
	for(int ip =0; ip<mCurrentFrame.NL; ip++)
	    if(mCurrentFrame.mvpMapLines[ip])
		sFoundLines.insert(mCurrentFrame.mvpMapLines[ip]);
	    
	//最后一次投影匹配，希望找到足够的新匹配地图点
	nadditionalLines = matcher.SearchByProjection(mCurrentFrame,vpCandidateKFs,sFoundLines,30,64);

	// Final optimization
	if( (nadditionalLines+inlinerLinesNum>=20) && inlinersNum>=40 )
	{
	    Optimizer::PoseOptimizationmain(&mCurrentFrame, inlinersNum, inlinerLinesNum);

	    for(int io =0; io<mCurrentFrame.N; io++)
		if(mCurrentFrame.mvbOutlier[io])
		    mCurrentFrame.mvpMapPoints[io]=NULL;
		
	    for(int io =0; io<mCurrentFrame.NL; io++)
		if(mCurrentFrame.mvbOutlierLines[io])
		    mCurrentFrame.mvpMapLines[io]=NULL;
	}
		
    }
}

//初步位姿估计结束，基于位姿计算当前图像帧的有效地图点匹配和内点数量
int Tracking::SetCurrentFrameMappointsAndInliers(Frame &CurrentFrame, std::vector<MapPoint*> vpMapPointMatches, std::set<MapPoint*> &sFound)
{
    int inliners = 0;
    
    //获取标定矩阵
    cv::Mat K = CurrentFrame.mK;
    const float fx = K.at<float>(0,0);
    const float fy = K.at<float>(1,1);
    const float cx = K.at<float>(0,2);
    const float cy = K.at<float>(1,2);
    
    //获取当前图像帧的旋转和平移矩阵
    cv::Mat Rcw = CurrentFrame.mTcw.rowRange(0,3).colRange(0,3);
    cv::Mat tcw = CurrentFrame.mTcw.rowRange(0,3).col(3);
    
    float thPoints = 5.991;
    
    for(size_t i=0, iend=vpMapPointMatches.size(); i<iend; i++)
    {
	//获取通过参考关键帧匹配到的地图点
	MapPoint* pMP = vpMapPointMatches[i];

	if(pMP)
	{
	    if(!pMP->isBad())
	    {
		//获取匹配地图点对应的点特征位置
		const cv::KeyPoint &kp = CurrentFrame.mvKeysUn[i];
		
		//获取地图点位置，计算图像坐标系坐标
		cv::Mat PosWorld = pMP->GetWorldPos();
		cv::Mat PosCamera = Rcw*PosWorld+tcw;
		float im1x, im1y;
		float invZ1 = 1.0/PosCamera.at<float>(2);
		im1x = fx*PosCamera.at<float>(0)*invZ1+cx;
		im1y = fy*PosCamera.at<float>(1)*invZ1+cy;
		 
		//计算重投影误差
		float squareError = (im1x-kp.pt.x)*(im1x-kp.pt.x)+(im1y-kp.pt.y)*(im1y-kp.pt.y);
		
		if(squareError < CurrentFrame.mvLevelSigma2[kp.octave]*thPoints)
		{
		    inliners++;
		    mCurrentFrame.mvpMapPoints[i]=pMP;
		    sFound.insert(pMP);
		}
		else
		{
		    mCurrentFrame.mvpMapPoints[i]=NULL;
		}
		 
	    }
	}
	
    }
    //cout << "inliners: " << inliners << endl;
    
    return inliners;
    
}

//原ORB-SLAM2中的重置函数
void Tracking::Reset()
{

    cout << "System Reseting" << endl;
    if(mpViewer)
    {
        mpViewer->RequestStop();
        while(!mpViewer->isStopped())
            usleep(3000);
    }

    // Reset Local Mapping
    cout << "Reseting Local Mapper...";
    mpLocalMapper->RequestReset();
    cout << " done" << endl;

    // Reset Loop Closing
    cout << "Reseting Loop Closing...";
    mpLoopClosing->RequestReset();
    cout << " done" << endl;

    // Clear BoW Database
    cout << "Reseting Database...";
    mpKeyFrameDB->clear();
    cout << " done" << endl;

    // Clear Map (this erase MapPoints and KeyFrames)
    mpMap->clear();

    KeyFrame::nNextId = 0;
    Frame::nNextId = 0;
    mState = NO_IMAGES_YET;

    if(mpInitializer)
    {
        delete mpInitializer;
        mpInitializer = static_cast<Initializer*>(NULL);
    }

    mlRelativeFramePoses.clear();
    mlpReferences.clear();
    mlFrameTimes.clear();
    mlbLost.clear();

    if(mpViewer)
        mpViewer->Release();
}


//使用线特征的删除函数
void Tracking::ResetBoth()
{

    cout << "System Reseting" << endl;
    if(mpViewer)
    {
        mpViewer->RequestStop();
        while(!mpViewer->isStopped())
            usleep(3000);
    }

    // Reset Local Mapping
    cout << "Reseting Local Mapper...";
    mpLocalMapper->RequestReset();
    cout << " done" << endl;

    // Reset Loop Closing
    cout << "Reseting Loop Closing...";
    mpLoopClosing->RequestReset();
    cout << " done" << endl;

    // Clear BoW Database
    cout << "Reseting Database...";
    mpKeyFrameDB->clear();
    cout << " done" << endl;

    // Clear Map (this erase MapPoints and KeyFrames 若使用线特征，还会删除地图中的地图线)
    cout << "Reseting Map...";
    mpMap->clearBoth();
    cout << " done" << endl;
    
    KeyFrame::nNextId = 0;
    Frame::nNextId = 0;
    mState = NO_IMAGES_YET;

    cout << "Reseting Initializer...";
    if(mpInitializer)
    {
        delete mpInitializer;
        mpInitializer = static_cast<Initializer*>(NULL);
    }
    cout << " done" << endl;
    
    mlRelativeFramePoses.clear();
    mlpReferences.clear();
    mlFrameTimes.clear();
    mlbLost.clear();

    cout << "Reseting Viewer..."; 
    if(mpViewer)
        mpViewer->Release(); 
    cout << " done" << endl;
    
}

void Tracking::ChangeCalibration(const string &strSettingPath)
{
    cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);
    float fx = fSettings["Camera.fx"];
    float fy = fSettings["Camera.fy"];
    float cx = fSettings["Camera.cx"];
    float cy = fSettings["Camera.cy"];

    cv::Mat K = cv::Mat::eye(3,3,CV_32F);
    K.at<float>(0,0) = fx;
    K.at<float>(1,1) = fy;
    K.at<float>(0,2) = cx;
    K.at<float>(1,2) = cy;
    K.copyTo(mK);

    cv::Mat DistCoef(4,1,CV_32F);
    DistCoef.at<float>(0) = fSettings["Camera.k1"];
    DistCoef.at<float>(1) = fSettings["Camera.k2"];
    DistCoef.at<float>(2) = fSettings["Camera.p1"];
    DistCoef.at<float>(3) = fSettings["Camera.p2"];
    const float k3 = fSettings["Camera.k3"];
    if(k3!=0)
    {
        DistCoef.resize(5);
        DistCoef.at<float>(4) = k3;
    }
    DistCoef.copyTo(mDistCoef);

    mbf = fSettings["Camera.bf"];

    Frame::mbInitialComputations = true;
}

void Tracking::InformOnlyTracking(const bool &flag)
{
    mbOnlyTracking = flag;
}

} //namespace ORB_SLAM
