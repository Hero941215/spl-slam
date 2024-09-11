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

#include "LocalMapping.h"
#include "LoopClosing.h"
#include "ORBmatcher.h"
#include "Linematcher.h"
#include "Optimizer.h"

#include <thread>
#include <mutex>

#include"Timer.h"

namespace PL_SLAM
{

LocalMapping::LocalMapping(Map *pMap, const float bMonocular):
    mbMonocular(bMonocular), mbResetRequested(false), mbFinishRequested(false), mbFinished(true), mpMap(pMap),
    mbAbortBA(false), mbStopped(false), mbStopRequested(false), mbNotStop(false), mbAcceptKeyFrames(true)
{
    //局部BA时间记录变量初始化
    localBAAllTime = 0.0;
    localBAAverageTime = 0.0;
    localBANumbers = 0;
}

void LocalMapping::SetLoopCloser(LoopClosing* pLoopCloser)
{
    mpLoopCloser = pLoopCloser;
}

void LocalMapping::SetTracker(Tracking *pTracker)
{
    mpTracker=pTracker;
}

void LocalMapping::Run()
{

    mbFinished = false;

    while(1)
    {
        // Tracking will see that Local Mapping is busy
        SetAcceptKeyFrames(false);

        // 等待处理的关键帧列表不为空
        if(CheckNewKeyFrames())
        {
            // BoW conversion and insertion in Map
            ProcessNewKeyFrame();
	    
            MapPointCulling();
            
            // Triangulate new MapPoints
            CreateNewMapPoints();
	    
	    // 等待处理的关键帧列表不为空
            if(!CheckNewKeyFrames())
            {
                SearchInNeighbors();
            }

            mbAbortBA = false;

            if(!CheckNewKeyFrames() && !stopRequested())
            {
	      
		Timer timer;
		timer.start();
                if(mpMap->KeyFramesInMap()>2)
                    Optimizer::LocalBundleAdjustment(mpCurrentKeyFrame,&mbAbortBA, mpMap);		
		double LocalBATime = timer.stop();
		cout << "局部BA时间: " << LocalBATime << " ms" <<endl;

                // Check redundant local Keyframes
                KeyFrameCulling();
            }

            mpLoopCloser->InsertKeyFrame(mpCurrentKeyFrame);
        }
        else if(Stop())
        {
            // Safe area to stop
            while(isStopped() && !CheckFinish())
            {
                usleep(3000);
            }
            if(CheckFinish())
                break;
        }

        ResetIfRequested();

        // Tracking will see that Local Mapping is not busy
	//后端不工作3毫秒
        SetAcceptKeyFrames(true);

        if(CheckFinish())
            break;

        usleep(3000);
    }

    SetFinish();
}

//点线特征后端工作主函数
void LocalMapping::RunBoth()
{

    mbFinished = false;

    while(1)
    {
        // Tracking will see that Local Mapping is busy
        SetAcceptKeyFrames(false);

        if(CheckNewKeyFrames())
        {
	    Timer timer;
	    timer.start();
            ProcessNewKeyFrameBoth();
	    double ProcessNewKeyFrameTime = timer.stop();
	    cout << "分线程处理当前关键帧时间: " << ProcessNewKeyFrameTime << " ms" <<endl;
	    
	    
	    timer.start();
	    thread threadMapPointCulling(&LocalMapping::MapPointCulling, this );
	    thread threadMapLineCulling(&LocalMapping::MapLineCulling, this );
	    threadMapPointCulling.join();
	    threadMapLineCulling.join();  
	    double MapCullingTime = timer.stop();
	    cout << "剔除临时地图时间: " << MapCullingTime << " ms" <<endl;

	    
	    //线特征发生彻底跟踪丢失时，使用点的共视图作为完成局部线地图的替代
	    vector<KeyFrame*> vpLocalKeyFramesPoints = mpCurrentKeyFrame->GetVectorCovisibleKeyFrames();
	    vector<KeyFrame*> vpLocalKeyFramesLines = mpCurrentKeyFrame->GetVectorCovisibleKeyFramesLines();
	    
	    if(!vpLocalKeyFramesPoints.empty() && !vpLocalKeyFramesLines.empty())
	    {
		timer.start();
		thread threadCreateNewMapPoints(&LocalMapping::CreateNewMapPoints, this );
		thread threadCreateNewMapLines(&LocalMapping::CreateNewMapLines, this );
		threadCreateNewMapPoints.join();
		threadCreateNewMapLines.join(); 
		double CreateNewMapTime = timer.stop();
		cout << "创建新地图时间: " << CreateNewMapTime << " ms" <<endl;
	    }
	    else if(!vpLocalKeyFramesPoints.empty())
	    {
		timer.start();
		thread threadCreateNewMapPoints(&LocalMapping::CreateNewMapPoints, this );
		thread threadCreateNewMapLines(&LocalMapping::CreateNewMapLinesWhenPointsOnly, this );
		threadCreateNewMapPoints.join();
		threadCreateNewMapLines.join(); 
		double CreateNewMapTime = timer.stop();
		cout << "线退化情况下，创建新地图时间: " << CreateNewMapTime << " ms" <<endl<<endl<<endl;
	    }
	    

            if(!CheckNewKeyFrames())
            {
                // 分线程：Find more matches in neighbor keyframes and fuse point and line duplications
		timer.start();
                thread threadSearchInNeighbors(&LocalMapping::SearchInNeighbors, this );
		thread threadSearchInNeighborsLines(&LocalMapping::SearchInNeighborsLines, this );
		threadSearchInNeighbors.join();
		threadSearchInNeighborsLines.join();	
		double SearchInNeighborsTime = timer.stop();
		cout << "冗余信息融合时间: " << SearchInNeighborsTime << " ms" <<endl;
		
            }

            mbAbortBA = false;
            if(!CheckNewKeyFrames() && !stopRequested())
            {
                // Local BA
                if(mpMap->KeyFramesInMap()>2)
		{
		    vector<KeyFrame*> vpLocalKeyFramesPoints = mpCurrentKeyFrame->GetVectorCovisibleKeyFrames();
		    vector<KeyFrame*> vpLocalKeyFramesLines = mpCurrentKeyFrame->GetVectorCovisibleKeyFramesLines();
	    
		    timer.start();   
		    if(!vpLocalKeyFramesPoints.empty() && !vpLocalKeyFramesLines.empty()) 
		    {
			Optimizer::LocalBundleAdjustmentmainOld(mpCurrentKeyFrame,&mbAbortBA, mpMap);
			cout << "执行点线同时使用的局部BA成功: 1" << endl;
		    }
		    else if(vpLocalKeyFramesLines.empty())
		    {
			Optimizer::LocalBundleAdjustment(mpCurrentKeyFrame,&mbAbortBA, mpMap);
			cout << "执行仅点同时使用的局部BA成功: 2" << endl;
			
		    }
		    else if(vpLocalKeyFramesPoints.empty())
		    {
			//理论上这一条不会执行，因为本系统以点特征为主特征，线特征为辅特征，若主特征全部丢失，则系统直接崩溃
			Optimizer::LocalBundleAdjustmentDoubleLines(mpCurrentKeyFrame,&mbAbortBA, mpMap);
			cout << "执行仅线同时使用的局部BA成功: 3" << endl;
		    }
		    
		    double LocalBATime = timer.stop();
		    cout << "当步局部BA时间: " << LocalBATime << " ms" <<endl<<endl<<endl;
		    
		    localBAAllTime += LocalBATime;
		    localBANumbers++;
		    localBAAverageTime = localBAAllTime / localBANumbers;
		    
		    cout << "平均局部BA时间: " << localBAAverageTime << " ms" <<endl<<endl<<endl;
		}
		
                timer.start();    
                KeyFrameCullingBoth();
		double  KeyFrameCullingTime = timer.stop();
		cout << "冗余关键帧剔除时间: " <<  KeyFrameCullingTime << " ms" <<endl<<endl;
            }

            mpLoopCloser->InsertKeyFrame(mpCurrentKeyFrame);
        }
        else if(Stop())
        {
            // Safe area to stop
            while(isStopped() && !CheckFinish())
            {
                usleep(3000);
            }
            if(CheckFinish())
                break;
        }

        //使用点线的重置函数
        ResetIfRequestedBoth();

        // Tracking will see that Local Mapping is not busy
	//后端不工作3毫秒
        SetAcceptKeyFrames(true);

        if(CheckFinish())
            break;

        usleep(3000);
	
	//cout << "后端执行结束: " << true << endl << endl;
    }

    SetFinish();
}

void LocalMapping::InsertKeyFrame(KeyFrame *pKF)
{
    unique_lock<mutex> lock(mMutexNewKFs);
    mlNewKeyFrames.push_back(pKF);
    mbAbortBA=true;
}


bool LocalMapping::CheckNewKeyFrames()
{
    unique_lock<mutex> lock(mMutexNewKFs);
    return(!mlNewKeyFrames.empty());
}

void LocalMapping::ProcessNewKeyFrame()
{
    {
        unique_lock<mutex> lock(mMutexNewKFs);
        mpCurrentKeyFrame = mlNewKeyFrames.front();
        mlNewKeyFrames.pop_front();
    }

    // Compute Bags of Words structures
    mpCurrentKeyFrame->ComputeBoW();

    // Associate MapPoints to the new keyframe and update normal and descriptor
    const vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();

    for(size_t i=0; i<vpMapPointMatches.size(); i++)
    {
        MapPoint* pMP = vpMapPointMatches[i];
        if(pMP)
        {
            if(!pMP->isBad())
            {
                if(!pMP->IsInKeyFrame(mpCurrentKeyFrame))
                {
                    pMP->AddObservation(mpCurrentKeyFrame, i);
                    pMP->UpdateNormalAndDepth();
                    pMP->ComputeDistinctiveDescriptors();
                }
                else // this can only happen for new stereo points inserted by the Tracking
                {
                    mlpRecentAddedMapPoints.push_back(pMP);
                }
            }
        }
    }    

    // Update links in the Covisibility Graph
    mpCurrentKeyFrame->UpdateConnections();

    // Insert Keyframe in Map
    mpMap->AddKeyFrame(mpCurrentKeyFrame);
}
//处理当前帧，更新点线特征观测信息和对应共视图信息，在地图中插入关键帧：
void LocalMapping::ProcessNewKeyFrameBoth()
{
  
    {
        unique_lock<mutex> lock(mMutexNewKFs);
        mpCurrentKeyFrame = mlNewKeyFrames.front();
        mlNewKeyFrames.pop_front();
    }
	
    thread threadProcessNewKeyFramePoints(&LocalMapping::ProcessNewKeyFramePoints, this );
    thread threadProcessNewKeyFrameLines(&LocalMapping::ProcessNewKeyFrameLines, this );
    threadProcessNewKeyFramePoints.join();
    threadProcessNewKeyFrameLines.join();
    
    // Insert Keyframe in Map
    mpMap->AddKeyFrame(mpCurrentKeyFrame);
}

void LocalMapping::ProcessNewKeyFramePoints()
{
    
    // Compute Bags of Words structures
    mpCurrentKeyFrame->ComputeBoW();

    // Associate MapPoints to the new keyframe and update normal and descriptor
    const vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();

    for(size_t i=0; i<vpMapPointMatches.size(); i++)
    {
        MapPoint* pMP = vpMapPointMatches[i];
        if(pMP)
        {
            if(!pMP->isBad())
            {
                if(!pMP->IsInKeyFrame(mpCurrentKeyFrame))
                {
                    pMP->AddObservation(mpCurrentKeyFrame, i);
                    pMP->UpdateNormalAndDepth();
                    pMP->ComputeDistinctiveDescriptors();
                }
                else //单目不运行这里
                {
                    mlpRecentAddedMapPoints.push_back(pMP);
                }
            }
        }
    }    
    
    mpCurrentKeyFrame->UpdateConnections();

}
//对当前关键帧的地图线向量进行处理，更新线共视图
void LocalMapping::ProcessNewKeyFrameLines()
{

    const vector<MapLine*> vpMapLineMatches = mpCurrentKeyFrame->GetMapLineMatches();

    for(size_t i=0; i<vpMapLineMatches.size(); i++)
    {
        MapLine* pML = vpMapLineMatches[i];
        if(pML)
        {
            if(!pML->isBad())
            {
                if(!pML->IsInKeyFrame(mpCurrentKeyFrame))
                {
                    pML->AddObservation(mpCurrentKeyFrame, i);
                    pML->UpdateNormalAndDepth();
                    pML->ComputeDistinctiveDescriptors();
		    pML->Update2DLineLength();
                }
                else // this can only happen for new stereo points inserted by the Tracking
                {
                    mlpRecentAddedMapLines.push_back(pML);
                }
            }
        }
    }    

    mpCurrentKeyFrame->UpdateConnectionsLines();
    
}

void LocalMapping::MapPointCulling()
{
    // Check Recent Added MapPoints
    list<MapPoint*>::iterator lit = mlpRecentAddedMapPoints.begin();
    const unsigned long int nCurrentKFid = mpCurrentKeyFrame->mnId;

    int nThObs;
    if(mbMonocular)
        nThObs = 2;
    else
        nThObs = 3;
    const int cnThObs = nThObs;

    while(lit!=mlpRecentAddedMapPoints.end())
    {
        MapPoint* pMP = *lit;
        if(pMP->isBad())
        {
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(pMP->GetFoundRatio()<0.25f )
        {
            pMP->SetBadFlag();
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(((int)nCurrentKFid-(int)pMP->mnFirstKFid)>=2 && pMP->Observations()<=cnThObs)
        {
            pMP->SetBadFlag();
            lit = mlpRecentAddedMapPoints.erase(lit);
        }
        else if(((int)nCurrentKFid-(int)pMP->mnFirstKFid)>=3)
            lit = mlpRecentAddedMapPoints.erase(lit);
        else
            lit++;
    }
}

//单双目通用
void LocalMapping::MapLineCulling()
{
    // Check Recent Added MapLines
    list<MapLine*>::iterator lit = mlpRecentAddedMapLines.begin();
    const unsigned long int nCurrentKFid = mpCurrentKeyFrame->mnId;

    int nThObs;
    if(mbMonocular)
        nThObs = 2;
    else
        nThObs = 3;
    const int cnThObs = nThObs;

    //
    while(lit!=mlpRecentAddedMapLines.end())
    {
        MapLine* pML = *lit;
        if(pML->isBad()) 
        {
            lit = mlpRecentAddedMapLines.erase(lit);
        }
        else if(pML->GetFoundRatio()<0.25f ) 
        {
            pML->SetBadFlag();
            lit = mlpRecentAddedMapLines.erase(lit);
        }
        else if(((int)nCurrentKFid-(int)pML->mnFirstKFid)>=2 && pML->Observations()<=cnThObs)      
        {
            pML->SetBadFlag();
            lit = mlpRecentAddedMapLines.erase(lit);
        }
        else if(((int)nCurrentKFid-(int)pML->mnFirstKFid)>=3)
            lit = mlpRecentAddedMapLines.erase(lit);
        else
            lit++;
    }
}

void LocalMapping::CreateNewMapPoints()
{
    // Retrieve neighbor keyframes in covisibility graph
    int nn = 10;
    if(mbMonocular)
        nn=20;

    const vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);

    ORBmatcher matcher(0.6,false);

    cv::Mat Rcw1 = mpCurrentKeyFrame->GetRotation();
    cv::Mat Rwc1 = Rcw1.t();
    cv::Mat tcw1 = mpCurrentKeyFrame->GetTranslation();
    cv::Mat Tcw1(3,4,CV_32F);
    Rcw1.copyTo(Tcw1.colRange(0,3));
    tcw1.copyTo(Tcw1.col(3));

    cv::Mat Ow1 = mpCurrentKeyFrame->GetCameraCenter();

    const float &fx1 = mpCurrentKeyFrame->fx;
    const float &fy1 = mpCurrentKeyFrame->fy;
    const float &cx1 = mpCurrentKeyFrame->cx;
    const float &cy1 = mpCurrentKeyFrame->cy;
    const float &invfx1 = mpCurrentKeyFrame->invfx;
    const float &invfy1 = mpCurrentKeyFrame->invfy;

    const float ratioFactor = 1.5f*mpCurrentKeyFrame->mfScaleFactor;

    int nnew=0;

    // Search matches with epipolar restriction and triangulate

    for(size_t i=0; i<vpNeighKFs.size(); i++)
    {
        if(i>0 && CheckNewKeyFrames())
            return;

        KeyFrame* pKF2 = vpNeighKFs[i];

        // Check first that baseline is not too short
        cv::Mat Ow2 = pKF2->GetCameraCenter();
        cv::Mat vBaseline = Ow2-Ow1;
        const float baseline = cv::norm(vBaseline);

        if(!mbMonocular)
        {
            if(baseline<pKF2->mb)
            continue;
        }
        else
        {
            const float medianDepthKF2 = pKF2->ComputeSceneMedianDepth(2);
            const float ratioBaselineDepth = baseline/medianDepthKF2;
            if(ratioBaselineDepth<0.01)
                continue;
        }

        // Compute Fundamental Matrix
        cv::Mat F12 = ComputeF12(mpCurrentKeyFrame,pKF2);

        // Search matches that fullfil epipolar constraint
        vector<pair<size_t,size_t> > vMatchedIndices;
        matcher.SearchForTriangulation(mpCurrentKeyFrame,pKF2,F12,vMatchedIndices,false);

        cv::Mat Rcw2 = pKF2->GetRotation();
        cv::Mat Rwc2 = Rcw2.t();
        cv::Mat tcw2 = pKF2->GetTranslation();
        cv::Mat Tcw2(3,4,CV_32F);
        Rcw2.copyTo(Tcw2.colRange(0,3));
        tcw2.copyTo(Tcw2.col(3));

        const float &fx2 = pKF2->fx;
        const float &fy2 = pKF2->fy;
        const float &cx2 = pKF2->cx;
        const float &cy2 = pKF2->cy;
        const float &invfx2 = pKF2->invfx;
        const float &invfy2 = pKF2->invfy;

        // Triangulate each match
        const int nmatches = vMatchedIndices.size();
        for(int ikp=0; ikp<nmatches; ikp++)
        {

            const int &idx1 = vMatchedIndices[ikp].first;
            const int &idx2 = vMatchedIndices[ikp].second;

            const cv::KeyPoint &kp1 = mpCurrentKeyFrame->mvKeysUn[idx1];
            const float kp1_ur=mpCurrentKeyFrame->mvuRight[idx1];
            bool bStereo1 = kp1_ur>=0;

            const cv::KeyPoint &kp2 = pKF2->mvKeysUn[idx2];
            const float kp2_ur = pKF2->mvuRight[idx2];
            bool bStereo2 = kp2_ur>=0;

            // Check parallax between rays
            cv::Mat xn1 = (cv::Mat_<float>(3,1) << (kp1.pt.x-cx1)*invfx1, (kp1.pt.y-cy1)*invfy1, 1.0);
            cv::Mat xn2 = (cv::Mat_<float>(3,1) << (kp2.pt.x-cx2)*invfx2, (kp2.pt.y-cy2)*invfy2, 1.0);

            cv::Mat ray1 = Rwc1*xn1;
            cv::Mat ray2 = Rwc2*xn2;
            const float cosParallaxRays = ray1.dot(ray2)/(cv::norm(ray1)*cv::norm(ray2));

            float cosParallaxStereo = cosParallaxRays+1;
            float cosParallaxStereo1 = cosParallaxStereo;
            float cosParallaxStereo2 = cosParallaxStereo;

            if(bStereo1)
                cosParallaxStereo1 = cos(2*atan2(mpCurrentKeyFrame->mb/2,mpCurrentKeyFrame->mvDepth[idx1]));
            else if(bStereo2)
                cosParallaxStereo2 = cos(2*atan2(pKF2->mb/2,pKF2->mvDepth[idx2]));

            cosParallaxStereo = min(cosParallaxStereo1,cosParallaxStereo2);

            cv::Mat x3D;
            if(cosParallaxRays<cosParallaxStereo && cosParallaxRays>0 && (bStereo1 || bStereo2 || cosParallaxRays<0.9998))
            {
                // Linear Triangulation Method
                // 见Initializer.cpp的Triangulate函数
                cv::Mat A(4,4,CV_32F);
                A.row(0) = xn1.at<float>(0)*Tcw1.row(2)-Tcw1.row(0);
                A.row(1) = xn1.at<float>(1)*Tcw1.row(2)-Tcw1.row(1);
                A.row(2) = xn2.at<float>(0)*Tcw2.row(2)-Tcw2.row(0);
                A.row(3) = xn2.at<float>(1)*Tcw2.row(2)-Tcw2.row(1);

                cv::Mat w,u,vt;
                cv::SVD::compute(A,w,u,vt,cv::SVD::MODIFY_A| cv::SVD::FULL_UV);

                x3D = vt.row(3).t();

                if(x3D.at<float>(3)==0)
                    continue;

                // Euclidean coordinates
                x3D = x3D.rowRange(0,3)/x3D.at<float>(3);
            }
            else if(bStereo1 && cosParallaxStereo1<cosParallaxStereo2)
            {
                x3D = mpCurrentKeyFrame->UnprojectStereo(idx1);                
            }
            else if(bStereo2 && cosParallaxStereo2<cosParallaxStereo1)
            {
                x3D = pKF2->UnprojectStereo(idx2);
            }
            else
                continue; //No stereo and very low parallax

            cv::Mat x3Dt = x3D.t();

            //Check triangulation in front of cameras
            float z1 = Rcw1.row(2).dot(x3Dt)+tcw1.at<float>(2);
            if(z1<=0)
                continue;

            float z2 = Rcw2.row(2).dot(x3Dt)+tcw2.at<float>(2);
            if(z2<=0)
                continue;

            //Check reprojection error in first keyframe
            const float &sigmaSquare1 = mpCurrentKeyFrame->mvLevelSigma2[kp1.octave];
            const float x1 = Rcw1.row(0).dot(x3Dt)+tcw1.at<float>(0);
            const float y1 = Rcw1.row(1).dot(x3Dt)+tcw1.at<float>(1);
            const float invz1 = 1.0/z1;

            if(!bStereo1)
            {
                float u1 = fx1*x1*invz1+cx1;
                float v1 = fy1*y1*invz1+cy1;
                float errX1 = u1 - kp1.pt.x;
                float errY1 = v1 - kp1.pt.y;
                if((errX1*errX1+errY1*errY1)>5.991*sigmaSquare1)
                    continue;
            }
            else
            {
                float u1 = fx1*x1*invz1+cx1;
                float u1_r = u1 - mpCurrentKeyFrame->mbf*invz1;
                float v1 = fy1*y1*invz1+cy1;
                float errX1 = u1 - kp1.pt.x;
                float errY1 = v1 - kp1.pt.y;
                float errX1_r = u1_r - kp1_ur;
                if((errX1*errX1+errY1*errY1+errX1_r*errX1_r)>7.8*sigmaSquare1)
                    continue;
            }

            //Check reprojection error in second keyframe
            const float sigmaSquare2 = pKF2->mvLevelSigma2[kp2.octave];
            const float x2 = Rcw2.row(0).dot(x3Dt)+tcw2.at<float>(0);
            const float y2 = Rcw2.row(1).dot(x3Dt)+tcw2.at<float>(1);
            const float invz2 = 1.0/z2;
            if(!bStereo2)
            {
                float u2 = fx2*x2*invz2+cx2;
                float v2 = fy2*y2*invz2+cy2;
                float errX2 = u2 - kp2.pt.x;
                float errY2 = v2 - kp2.pt.y;
                if((errX2*errX2+errY2*errY2)>5.991*sigmaSquare2)
                    continue;
            }
            else
            {
                float u2 = fx2*x2*invz2+cx2;
                float u2_r = u2 - mpCurrentKeyFrame->mbf*invz2;
                float v2 = fy2*y2*invz2+cy2;
                float errX2 = u2 - kp2.pt.x;
                float errY2 = v2 - kp2.pt.y;
                float errX2_r = u2_r - kp2_ur;
                if((errX2*errX2+errY2*errY2+errX2_r*errX2_r)>7.8*sigmaSquare2)
                    continue;
            }

            //Check scale consistency
            cv::Mat normal1 = x3D-Ow1;
            float dist1 = cv::norm(normal1);

            cv::Mat normal2 = x3D-Ow2;
            float dist2 = cv::norm(normal2);

            if(dist1==0 || dist2==0)
                continue;

            const float ratioDist = dist2/dist1;
            const float ratioOctave = mpCurrentKeyFrame->mvScaleFactors[kp1.octave]/pKF2->mvScaleFactors[kp2.octave];

            if(ratioDist*ratioFactor<ratioOctave || ratioDist>ratioOctave*ratioFactor)
                continue;

            // Triangulation is succesfull
            MapPoint* pMP = new MapPoint(x3D,mpCurrentKeyFrame,mpMap);

            pMP->AddObservation(mpCurrentKeyFrame,idx1);            
            pMP->AddObservation(pKF2,idx2);

            mpCurrentKeyFrame->AddMapPoint(pMP,idx1);
            pKF2->AddMapPoint(pMP,idx2);
            pMP->ComputeDistinctiveDescriptors();
            pMP->UpdateNormalAndDepth();

            mpMap->AddMapPoint(pMP);
            mlpRecentAddedMapPoints.push_back(pMP);

            nnew++;
        }
    }
}

//仅单目
void LocalMapping::CreateNewMapLines()
{
    // Retrieve neighbor keyframes in covisibility graph
    int nn = 10;
    if(mbMonocular)
        nn=20; 
    const vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFramesLines(nn);

    Linematcher matcher(0.7,false,false,0.15);

    cv::Mat Rcw1 = mpCurrentKeyFrame->GetRotation();
    cv::Mat Rwc1 = Rcw1.t();
    cv::Mat tcw1 = mpCurrentKeyFrame->GetTranslation();
    cv::Mat Tcw1(3,4,CV_32F);
    Rcw1.copyTo(Tcw1.colRange(0,3));
    tcw1.copyTo(Tcw1.col(3));
    cv::Mat Ow1 = mpCurrentKeyFrame->GetCameraCenter();

    //相机标定矩阵
    const float &fx1 = mpCurrentKeyFrame->fx;
    const float &fy1 = mpCurrentKeyFrame->fy;
    const float &cx1 = mpCurrentKeyFrame->cx;
    const float &cy1 = mpCurrentKeyFrame->cy;
    const float &invfx1 = mpCurrentKeyFrame->invfx;
    const float &invfy1 = mpCurrentKeyFrame->invfy;

    const float ratioFactor = 1.5f*mpCurrentKeyFrame->mfScaleFactorLines;
    int nnew=0;

    // Search matches with epipolar restriction and triangulate
    for(size_t i=0; i<vpNeighKFs.size(); i++)
    {

        if(i>0 && CheckNewKeyFrames())
            return;

        KeyFrame* pKF2 = vpNeighKFs[i];

        // Check first that baseline is not too short
        cv::Mat Ow2 = pKF2->GetCameraCenter();
        cv::Mat vBaseline = Ow2-Ow1;
        const float baseline = cv::norm(vBaseline);
        const float medianDepthKF2 = pKF2->ComputeSceneMedianDepthBoth(2);
        const float ratioBaselineDepth = baseline/medianDepthKF2;
        if(ratioBaselineDepth<0.01)
            continue;
        
        // Compute Fundamental Matrix
        cv::Mat F12 = ComputeF12(mpCurrentKeyFrame,pKF2);

        // Search matches that fullfil epipolar constraint
        vector<pair<size_t,size_t> > vMatchedIndices;
        matcher.SearchForTriangulation(mpCurrentKeyFrame,pKF2,F12,vMatchedIndices);

        cv::Mat Rcw2 = pKF2->GetRotation();
        cv::Mat Rwc2 = Rcw2.t();
        cv::Mat tcw2 = pKF2->GetTranslation();
        cv::Mat Tcw2(3,4,CV_32F);
        Rcw2.copyTo(Tcw2.colRange(0,3));
        tcw2.copyTo(Tcw2.col(3));

        const float &fx2 = pKF2->fx;
        const float &fy2 = pKF2->fy;
        const float &cx2 = pKF2->cx;
        const float &cy2 = pKF2->cy;
        const float &invfx2 = pKF2->invfx;
        const float &invfy2 = pKF2->invfy;

        // Triangulate each match
        const int nmatches = vMatchedIndices.size();
        for(int ikp=0; ikp<nmatches; ikp++)
        {
            const int &idx1 = vMatchedIndices[ikp].first;
            const int &idx2 = vMatchedIndices[ikp].second;

            const cv::KeyPoint &kp1 = mpCurrentKeyFrame->mvMidPointsUn[idx1];	    
            const cv::KeyPoint &kp2 = pKF2->mvMidPointsUn[idx2];
	    
            // Check parallax between rays
            cv::Mat xn1 = (cv::Mat_<float>(3,1) << (kp1.pt.x-cx1)*invfx1, (kp1.pt.y-cy1)*invfy1, 1.0);
            cv::Mat xn2 = (cv::Mat_<float>(3,1) << (kp2.pt.x-cx2)*invfx2, (kp2.pt.y-cy2)*invfy2, 1.0);
	    
            cv::Mat ray1 = Rwc1*xn1;
            cv::Mat ray2 = Rwc2*xn2;

            const float cosParallaxRays = ray1.dot(ray2)/(cv::norm(ray1)*cv::norm(ray2));

            cv::Mat x3D;
            if(cosParallaxRays>0 && cosParallaxRays<0.9998)
            {
                // Linear Triangulation Method
                cv::Mat A(4,4,CV_32F);
                A.row(0) = xn1.at<float>(0)*Tcw1.row(2)-Tcw1.row(0);
                A.row(1) = xn1.at<float>(1)*Tcw1.row(2)-Tcw1.row(1);
                A.row(2) = xn2.at<float>(0)*Tcw2.row(2)-Tcw2.row(0);
                A.row(3) = xn2.at<float>(1)*Tcw2.row(2)-Tcw2.row(1);

                cv::Mat w,u,vt;
                cv::SVD::compute(A,w,u,vt,cv::SVD::MODIFY_A| cv::SVD::FULL_UV);

                x3D = vt.row(3).t();

                if(x3D.at<float>(3)==0)
                    continue;

                // Euclidean coordinates
                x3D = x3D.rowRange(0,3)/x3D.at<float>(3);

            }
            else
                continue; //No stereo and very low parallax

            cv::Mat x3Dt = x3D.t();

            //Check triangulation in front of cameras
            float z1 = Rcw1.row(2).dot(x3Dt)+tcw1.at<float>(2);
            if(z1<=0)
                continue;

            float z2 = Rcw2.row(2).dot(x3Dt)+tcw2.at<float>(2);
            if(z2<=0)
                continue;

            //Check reprojection error in first keyframe
	    const KeyLine &kL1 = mpCurrentKeyFrame->mvLinesUn[idx1];
	    Vector3f sp_l1; sp_l1 << kL1.startPointX, kL1.startPointY, 1.0;
	    Vector3f ep_l1; ep_l1 << kL1.endPointX,   kL1.endPointY,   1.0;
	    Vector3f le_l1; le_l1 << sp_l1.cross(ep_l1);
	    //le_l1 = le_l1 / std::sqrt( le_l1(0)*le_l1(0) + le_l1(1)*le_l1(1) + le_l1(2)*le_l1(2));
	    le_l1 = le_l1 / std::sqrt( le_l1(0)*le_l1(0) + le_l1(1)*le_l1(1));
	    
            const float &sigmaSquare1 = mpCurrentKeyFrame->mvLevelSigma2Lines[kp1.octave];
            const float x1 = Rcw1.row(0).dot(x3Dt)+tcw1.at<float>(0);
            const float y1 = Rcw1.row(1).dot(x3Dt)+tcw1.at<float>(1);
            const float invz1 = 1.0/z1;

            float u1 = fx1*x1*invz1+cx1;
            float v1 = fy1*y1*invz1+cy1;       
	    float error1 = fabs(le_l1(0)*u1 + le_l1(1)*v1 + le_l1(2));
            if( error1>3.841*sigmaSquare1 )
                continue;
           
            //Check reprojection error in second keyframe
	    const KeyLine &kL2 = pKF2->mvLinesUn[idx2];
	    Vector3f sp_l2; sp_l2 << kL2.startPointX, kL2.startPointY, 1.0;
	    Vector3f ep_l2; ep_l2 << kL2.endPointX,   kL2.endPointY,   1.0;
	    Vector3f le_l2; le_l2 << sp_l2.cross(ep_l2);
	    //le_l2 = le_l2 / std::sqrt( le_l2(0)*le_l2(0) + le_l2(1)*le_l2(1) + le_l2(2)*le_l2(2));
	    le_l2 = le_l2 / std::sqrt( le_l2(0)*le_l2(0) + le_l2(1)*le_l2(1));
	    
            const float sigmaSquare2 = pKF2->mvLevelSigma2Lines[kp2.octave];
            const float x2 = Rcw2.row(0).dot(x3Dt)+tcw2.at<float>(0);
            const float y2 = Rcw2.row(1).dot(x3Dt)+tcw2.at<float>(1);
            const float invz2 = 1.0/z2;

	    float u2 = fx2*x2*invz2+cx2;
            float v2 = fy2*y2*invz2+cy2;
	    float error2 = fabs(le_l2(0)*u2 + le_l2(1)*v2 + le_l2(2));
            if(error2>3.841*sigmaSquare2)
                continue;
            
            //Check scale consistency
            cv::Mat normal1 = x3D-Ow1;
            float dist1 = cv::norm(normal1);

            cv::Mat normal2 = x3D-Ow2;
            float dist2 = cv::norm(normal2);

            if(dist1==0 || dist2==0)
                continue;

            const float ratioDist = dist2/dist1;
            const float ratioOctave = mpCurrentKeyFrame->mvScaleFactorsLines[kp1.octave]/pKF2->mvScaleFactorsLines[kp2.octave];

            if(ratioDist*ratioFactor<ratioOctave || ratioDist>ratioOctave*ratioFactor)
                continue;
	    
	    cv::Mat FirPx3D;
	    cv::Mat EndPx3D;

	    xn1 = (cv::Mat_<float>(3,1) << (kL1.startPointX-cx1)*invfx1, (kL1.startPointY-cy1)*invfy1, 1.0);
            xn2 = (cv::Mat_<float>(3,1) << (kL2.startPointX-cx2)*invfx2, (kL2.startPointY-cy2)*invfy2, 1.0);
	    
	    cv::Mat AFir(4,4,CV_32F);
            AFir.row(0) = xn1.at<float>(0)*Tcw1.row(2)-Tcw1.row(0);
            AFir.row(1) = xn1.at<float>(1)*Tcw1.row(2)-Tcw1.row(1);
            AFir.row(2) = xn2.at<float>(0)*Tcw2.row(2)-Tcw2.row(0);
            AFir.row(3) = xn2.at<float>(1)*Tcw2.row(2)-Tcw2.row(1);
	    cv::Mat wFir,uFir,vtFir;
	    cv::SVD::compute(AFir,wFir,uFir,vtFir,cv::SVD::MODIFY_A| cv::SVD::FULL_UV);

	    FirPx3D = vtFir.row(3).t();
	    if(FirPx3D.at<float>(3)==0)
                continue; 
            // Euclidean coordinates
            FirPx3D = FirPx3D.rowRange(0,3)/FirPx3D.at<float>(3);
	    
	    cv::Mat FirPx3Dt = FirPx3D.t();
            //Check triangulation in front of cameras
            float zFir1 = Rcw1.row(2).dot(FirPx3Dt)+tcw1.at<float>(2);
            if(zFir1<=0)
                continue;

            float zFir2 = Rcw2.row(2).dot(FirPx3Dt)+tcw2.at<float>(2);
            if(zFir2<=0)
                continue;

	    xn1 = (cv::Mat_<float>(3,1) << (kL1.endPointX-cx1)*invfx1, (kL1.endPointY-cy1)*invfy1, 1.0);
            xn2 = (cv::Mat_<float>(3,1) << (kL2.endPointX-cx2)*invfx2, (kL2.endPointY-cy2)*invfy2, 1.0);

	    cv::Mat AEnd(4,4,CV_32F);
            AEnd.row(0) = xn1.at<float>(0)*Tcw1.row(2)-Tcw1.row(0);
            AEnd.row(1) = xn1.at<float>(1)*Tcw1.row(2)-Tcw1.row(1);
            AEnd.row(2) = xn2.at<float>(0)*Tcw2.row(2)-Tcw2.row(0);
            AEnd.row(3) = xn2.at<float>(1)*Tcw2.row(2)-Tcw2.row(1);
	    cv::Mat wEnd,uEnd,vtEnd;
	    cv::SVD::compute(AEnd,wEnd,uEnd,vtEnd,cv::SVD::MODIFY_A| cv::SVD::FULL_UV);

	    EndPx3D = vtEnd.row(3).t();
	    if(EndPx3D.at<float>(3)==0)
                continue;
            // Euclidean coordinates
            EndPx3D = EndPx3D.rowRange(0,3)/EndPx3D.at<float>(3);
	    
	    cv::Mat EndPx3Dt = EndPx3D.t();
            //Check triangulation in front of cameras
            float zEnd1 = Rcw1.row(2).dot(EndPx3Dt)+tcw1.at<float>(2);
            if(zEnd1<=0)
                continue;

            float zEnd2 = Rcw2.row(2).dot(EndPx3Dt)+tcw2.at<float>(2);
            if(zEnd2<=0)
                continue;
	    	      
            // Triangulation is succesfull
            MapLine* pML = new MapLine(FirPx3D, EndPx3D, x3D, mpCurrentKeyFrame, mpMap);

            pML->AddObservation(mpCurrentKeyFrame,idx1);            
            pML->AddObservation(pKF2,idx2);

            mpCurrentKeyFrame->AddMapLine(pML,idx1);
            pKF2->AddMapLine(pML,idx2);

            pML->ComputeDistinctiveDescriptors();
            pML->UpdateNormalAndDepth();    
	    pML->Update2DLineLength();

            mpMap->AddMapLine(pML);
            mlpRecentAddedMapLines.push_back(pML);

            nnew++;
        }
    }
}

//当前端线特征彻底丢失的情况下，使用点共视图完成地图创建
//这里创建的地图线无需经过地图线融合的步骤，因为这是最新创建的肯定不会发生冲突
void LocalMapping::CreateNewMapLinesWhenPointsOnly()
{
    // Retrieve neighbor keyframes in covisibility graph
    int nn = 10;
    if(mbMonocular)
        nn=20; 
    
    //获取足够数量的点共视图
    const vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);

    //使用线特征匹配器
    Linematcher matcher(0.6,false,false,0.15);

    cv::Mat Rcw1 = mpCurrentKeyFrame->GetRotation();
    cv::Mat Rwc1 = Rcw1.t();
    cv::Mat tcw1 = mpCurrentKeyFrame->GetTranslation();
    cv::Mat Tcw1(3,4,CV_32F);
    Rcw1.copyTo(Tcw1.colRange(0,3));
    tcw1.copyTo(Tcw1.col(3));
    cv::Mat Ow1 = mpCurrentKeyFrame->GetCameraCenter();

    //相机标定矩阵
    const float &fx1 = mpCurrentKeyFrame->fx;
    const float &fy1 = mpCurrentKeyFrame->fy;
    const float &cx1 = mpCurrentKeyFrame->cx;
    const float &cy1 = mpCurrentKeyFrame->cy;
    const float &invfx1 = mpCurrentKeyFrame->invfx;
    const float &invfy1 = mpCurrentKeyFrame->invfy;

    const float ratioFactor = 1.5f*mpCurrentKeyFrame->mfScaleFactorLines;
    int nnew=0;

    // Search matches with epipolar restriction and triangulate
    for(size_t i=0; i<vpNeighKFs.size(); i++)
    {

        if(i>0 && CheckNewKeyFrames())
            return;

        KeyFrame* pKF2 = vpNeighKFs[i];

        // Check first that baseline is not too short
        cv::Mat Ow2 = pKF2->GetCameraCenter();
        cv::Mat vBaseline = Ow2-Ow1;
        const float baseline = cv::norm(vBaseline);
	//这里此时实际上只有点特征，因为线特征已经彻底丢失了
        const float medianDepthKF2 = pKF2->ComputeSceneMedianDepthBoth(2);
        const float ratioBaselineDepth = baseline/medianDepthKF2;
        if(ratioBaselineDepth<0.01)
            continue;
        
        // Compute Fundamental Matrix
        cv::Mat F12 = ComputeF12(mpCurrentKeyFrame,pKF2);

        // Search matches that fullfil epipolar constraint
        vector<pair<size_t,size_t> > vMatchedIndices;
        matcher.SearchForTriangulation(mpCurrentKeyFrame,pKF2,F12,vMatchedIndices);

        cv::Mat Rcw2 = pKF2->GetRotation();
        cv::Mat Rwc2 = Rcw2.t();
        cv::Mat tcw2 = pKF2->GetTranslation();
        cv::Mat Tcw2(3,4,CV_32F);
        Rcw2.copyTo(Tcw2.colRange(0,3));
        tcw2.copyTo(Tcw2.col(3));

        const float &fx2 = pKF2->fx;
        const float &fy2 = pKF2->fy;
        const float &cx2 = pKF2->cx;
        const float &cy2 = pKF2->cy;
        const float &invfx2 = pKF2->invfx;
        const float &invfy2 = pKF2->invfy;

        // Triangulate each match
        const int nmatches = vMatchedIndices.size();
        for(int ikp=0; ikp<nmatches; ikp++)
        {
            const int &idx1 = vMatchedIndices[ikp].first;
            const int &idx2 = vMatchedIndices[ikp].second;

            const cv::KeyPoint &kp1 = mpCurrentKeyFrame->mvMidPointsUn[idx1];	    
            const cv::KeyPoint &kp2 = pKF2->mvMidPointsUn[idx2];
	    
            // Check parallax between rays
            cv::Mat xn1 = (cv::Mat_<float>(3,1) << (kp1.pt.x-cx1)*invfx1, (kp1.pt.y-cy1)*invfy1, 1.0);
            cv::Mat xn2 = (cv::Mat_<float>(3,1) << (kp2.pt.x-cx2)*invfx2, (kp2.pt.y-cy2)*invfy2, 1.0);
	    
            cv::Mat ray1 = Rwc1*xn1;
            cv::Mat ray2 = Rwc2*xn2;

            const float cosParallaxRays = ray1.dot(ray2)/(cv::norm(ray1)*cv::norm(ray2));

            cv::Mat x3D;
            if(cosParallaxRays>0 && cosParallaxRays<0.9998)
            {
                // Linear Triangulation Method
                cv::Mat A(4,4,CV_32F);
                A.row(0) = xn1.at<float>(0)*Tcw1.row(2)-Tcw1.row(0);
                A.row(1) = xn1.at<float>(1)*Tcw1.row(2)-Tcw1.row(1);
                A.row(2) = xn2.at<float>(0)*Tcw2.row(2)-Tcw2.row(0);
                A.row(3) = xn2.at<float>(1)*Tcw2.row(2)-Tcw2.row(1);

                cv::Mat w,u,vt;
                cv::SVD::compute(A,w,u,vt,cv::SVD::MODIFY_A| cv::SVD::FULL_UV);

                x3D = vt.row(3).t();

                if(x3D.at<float>(3)==0)
                    continue;

                // Euclidean coordinates
                x3D = x3D.rowRange(0,3)/x3D.at<float>(3);

            }
            else
                continue; //No stereo and very low parallax

            cv::Mat x3Dt = x3D.t();

            //Check triangulation in front of cameras
            float z1 = Rcw1.row(2).dot(x3Dt)+tcw1.at<float>(2);
            if(z1<=0)
                continue;

            float z2 = Rcw2.row(2).dot(x3Dt)+tcw2.at<float>(2);
            if(z2<=0)
                continue;

            //Check reprojection error in first keyframe
	    const KeyLine &kL1 = mpCurrentKeyFrame->mvLinesUn[idx1];
	    Vector3f sp_l1; sp_l1 << kL1.startPointX, kL1.startPointY, 1.0;
	    Vector3f ep_l1; ep_l1 << kL1.endPointX,   kL1.endPointY,   1.0;
	    Vector3f le_l1; le_l1 << sp_l1.cross(ep_l1);
	    //le_l1 = le_l1 / std::sqrt( le_l1(0)*le_l1(0) + le_l1(1)*le_l1(1) + le_l1(2)*le_l1(2));
	    le_l1 = le_l1 / std::sqrt( le_l1(0)*le_l1(0) + le_l1(1)*le_l1(1));
	    
            const float &sigmaSquare1 = mpCurrentKeyFrame->mvLevelSigma2Lines[kp1.octave];
            const float x1 = Rcw1.row(0).dot(x3Dt)+tcw1.at<float>(0);
            const float y1 = Rcw1.row(1).dot(x3Dt)+tcw1.at<float>(1);
            const float invz1 = 1.0/z1;

            float u1 = fx1*x1*invz1+cx1;
            float v1 = fy1*y1*invz1+cy1;       
	    float error1 = fabs(le_l1(0)*u1 + le_l1(1)*v1 + le_l1(2));
            if( error1>3.841*sigmaSquare1 )
                continue;
           
            //Check reprojection error in second keyframe
	    const KeyLine &kL2 = pKF2->mvLinesUn[idx2];
	    Vector3f sp_l2; sp_l2 << kL2.startPointX, kL2.startPointY, 1.0;
	    Vector3f ep_l2; ep_l2 << kL2.endPointX,   kL2.endPointY,   1.0;
	    Vector3f le_l2; le_l2 << sp_l2.cross(ep_l2);
	    //le_l2 = le_l2 / std::sqrt( le_l2(0)*le_l2(0) + le_l2(1)*le_l2(1) + le_l2(2)*le_l2(2));
	    le_l2 = le_l2 / std::sqrt( le_l2(0)*le_l2(0) + le_l2(1)*le_l2(1));
	    
            const float sigmaSquare2 = pKF2->mvLevelSigma2Lines[kp2.octave];
            const float x2 = Rcw2.row(0).dot(x3Dt)+tcw2.at<float>(0);
            const float y2 = Rcw2.row(1).dot(x3Dt)+tcw2.at<float>(1);
            const float invz2 = 1.0/z2;

	    float u2 = fx2*x2*invz2+cx2;
            float v2 = fy2*y2*invz2+cy2;
	    float error2 = fabs(le_l2(0)*u2 + le_l2(1)*v2 + le_l2(2));
            if(error2>3.841*sigmaSquare2)
                continue;
            
            //Check scale consistency
            cv::Mat normal1 = x3D-Ow1;
            float dist1 = cv::norm(normal1);

            cv::Mat normal2 = x3D-Ow2;
            float dist2 = cv::norm(normal2);

            if(dist1==0 || dist2==0)
                continue;

            const float ratioDist = dist2/dist1;
            const float ratioOctave = mpCurrentKeyFrame->mvScaleFactorsLines[kp1.octave]/pKF2->mvScaleFactorsLines[kp2.octave];

            if(ratioDist*ratioFactor<ratioOctave || ratioDist>ratioOctave*ratioFactor)
                continue;
	    
	    cv::Mat FirPx3D;
	    cv::Mat EndPx3D;

	    xn1 = (cv::Mat_<float>(3,1) << (kL1.startPointX-cx1)*invfx1, (kL1.startPointY-cy1)*invfy1, 1.0);
            xn2 = (cv::Mat_<float>(3,1) << (kL2.startPointX-cx2)*invfx2, (kL2.startPointY-cy2)*invfy2, 1.0);
	    
	    cv::Mat AFir(4,4,CV_32F);
            AFir.row(0) = xn1.at<float>(0)*Tcw1.row(2)-Tcw1.row(0);
            AFir.row(1) = xn1.at<float>(1)*Tcw1.row(2)-Tcw1.row(1);
            AFir.row(2) = xn2.at<float>(0)*Tcw2.row(2)-Tcw2.row(0);
            AFir.row(3) = xn2.at<float>(1)*Tcw2.row(2)-Tcw2.row(1);
	    cv::Mat wFir,uFir,vtFir;
	    cv::SVD::compute(AFir,wFir,uFir,vtFir,cv::SVD::MODIFY_A| cv::SVD::FULL_UV);

	    FirPx3D = vtFir.row(3).t();
	    if(FirPx3D.at<float>(3)==0)
                continue;
	    
            // Euclidean coordinates
            FirPx3D = FirPx3D.rowRange(0,3)/FirPx3D.at<float>(3);
	    
	    cv::Mat FirPx3Dt = FirPx3D.t();
            //Check triangulation in front of cameras
            float zFir1 = Rcw1.row(2).dot(FirPx3Dt)+tcw1.at<float>(2);
            if(zFir1<=0)
                continue;

            float zFir2 = Rcw2.row(2).dot(FirPx3Dt)+tcw2.at<float>(2);
            if(zFir2<=0)
                continue;

	    xn1 = (cv::Mat_<float>(3,1) << (kL1.endPointX-cx1)*invfx1, (kL1.endPointY-cy1)*invfy1, 1.0);
            xn2 = (cv::Mat_<float>(3,1) << (kL2.endPointX-cx2)*invfx2, (kL2.endPointY-cy2)*invfy2, 1.0);

	    cv::Mat AEnd(4,4,CV_32F);
            AEnd.row(0) = xn1.at<float>(0)*Tcw1.row(2)-Tcw1.row(0);
            AEnd.row(1) = xn1.at<float>(1)*Tcw1.row(2)-Tcw1.row(1);
            AEnd.row(2) = xn2.at<float>(0)*Tcw2.row(2)-Tcw2.row(0);
            AEnd.row(3) = xn2.at<float>(1)*Tcw2.row(2)-Tcw2.row(1);
	    cv::Mat wEnd,uEnd,vtEnd;
	    cv::SVD::compute(AEnd,wEnd,uEnd,vtEnd,cv::SVD::MODIFY_A| cv::SVD::FULL_UV);

	    EndPx3D = vtEnd.row(3).t();
	    if(EndPx3D.at<float>(3)==0)
                continue;
            // Euclidean coordinates
            EndPx3D = EndPx3D.rowRange(0,3)/EndPx3D.at<float>(3);
	    	
	    cv::Mat EndPx3Dt = EndPx3D.t();
            //Check triangulation in front of cameras
            float zEnd1 = Rcw1.row(2).dot(EndPx3Dt)+tcw1.at<float>(2);
            if(zEnd1<=0)
                continue;

            float zEnd2 = Rcw2.row(2).dot(EndPx3Dt)+tcw2.at<float>(2);
            if(zEnd2<=0)
                continue;
	    
            // Triangulation is succesfull
            MapLine* pML = new MapLine(FirPx3D, EndPx3D, x3D, mpCurrentKeyFrame, mpMap);

            pML->AddObservation(mpCurrentKeyFrame,idx1);            
            pML->AddObservation(pKF2,idx2);

            mpCurrentKeyFrame->AddMapLine(pML,idx1);
            pKF2->AddMapLine(pML,idx2);

            pML->ComputeDistinctiveDescriptors();
            pML->UpdateNormalAndDepth();    
	    pML->Update2DLineLength();

            mpMap->AddMapLine(pML);
            mlpRecentAddedMapLines.push_back(pML);

            nnew++;
        }
    }
}


void LocalMapping::SearchInNeighbors()
{
    // Retrieve neighbor keyframes
    int nn = 10;
    if(mbMonocular)
        nn=20;
    const vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);
    vector<KeyFrame*> vpTargetKFs;
    for(vector<KeyFrame*>::const_iterator vit=vpNeighKFs.begin(), vend=vpNeighKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;
        if(pKFi->isBad() || pKFi->mnFuseTargetForKF == mpCurrentKeyFrame->mnId)
            continue;
        vpTargetKFs.push_back(pKFi);
        pKFi->mnFuseTargetForKF = mpCurrentKeyFrame->mnId;

        // Extend to some second neighbors
        const vector<KeyFrame*> vpSecondNeighKFs = pKFi->GetBestCovisibilityKeyFrames(5);
        for(vector<KeyFrame*>::const_iterator vit2=vpSecondNeighKFs.begin(), vend2=vpSecondNeighKFs.end(); vit2!=vend2; vit2++)
        {
            KeyFrame* pKFi2 = *vit2;
            if(pKFi2->isBad() || pKFi2->mnFuseTargetForKF==mpCurrentKeyFrame->mnId || pKFi2->mnId==mpCurrentKeyFrame->mnId)
                continue;
            vpTargetKFs.push_back(pKFi2);
        }
    }


    // Search matches by projection from current KF in target KFs
    ORBmatcher matcher;
    vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    for(vector<KeyFrame*>::iterator vit=vpTargetKFs.begin(), vend=vpTargetKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;

        matcher.Fuse(pKFi,vpMapPointMatches);
    }

    // Search matches by projection from target KFs in current KF
    vector<MapPoint*> vpFuseCandidates;
    vpFuseCandidates.reserve(vpTargetKFs.size()*vpMapPointMatches.size());

    for(vector<KeyFrame*>::iterator vitKF=vpTargetKFs.begin(), vendKF=vpTargetKFs.end(); vitKF!=vendKF; vitKF++)
    {
        KeyFrame* pKFi = *vitKF;

        vector<MapPoint*> vpMapPointsKFi = pKFi->GetMapPointMatches();

        for(vector<MapPoint*>::iterator vitMP=vpMapPointsKFi.begin(), vendMP=vpMapPointsKFi.end(); vitMP!=vendMP; vitMP++)
        {
            MapPoint* pMP = *vitMP;
            if(!pMP)
                continue;
            if(pMP->isBad() || pMP->mnFuseCandidateForKF == mpCurrentKeyFrame->mnId)
                continue;
            pMP->mnFuseCandidateForKF = mpCurrentKeyFrame->mnId;
            vpFuseCandidates.push_back(pMP);
        }
    }

    matcher.Fuse(mpCurrentKeyFrame,vpFuseCandidates);


    // Update points
    vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    for(size_t i=0, iend=vpMapPointMatches.size(); i<iend; i++)
    {
        MapPoint* pMP=vpMapPointMatches[i];
        if(pMP)
        {
            if(!pMP->isBad())
            {
                pMP->ComputeDistinctiveDescriptors();
                pMP->UpdateNormalAndDepth();
            }
        }
    }

    // Update connections in covisibility graph
    mpCurrentKeyFrame->UpdateConnections();
}

void LocalMapping::SearchInNeighborsLines()
{
    // Retrieve neighbor keyframes
    int nn = 10;
    if(mbMonocular)
        nn=20;
    const vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFramesLines(nn);

    vector<KeyFrame*> vpTargetKFs;
    for(vector<KeyFrame*>::const_iterator vit=vpNeighKFs.begin(), vend=vpNeighKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;
        if(pKFi->isBadLines() || pKFi->mnFuseTargetForKFLines == mpCurrentKeyFrame->mnId)
            continue;
        vpTargetKFs.push_back(pKFi);
        pKFi->mnFuseTargetForKFLines = mpCurrentKeyFrame->mnId;

        // Extend to some second neighbors
        const vector<KeyFrame*> vpSecondNeighKFs = pKFi->GetBestCovisibilityKeyFramesLines(5);
        for(vector<KeyFrame*>::const_iterator vit2=vpSecondNeighKFs.begin(), vend2=vpSecondNeighKFs.end(); vit2!=vend2; vit2++)
        {
            KeyFrame* pKFi2 = *vit2;
            if(pKFi2->isBadLines() || pKFi2->mnFuseTargetForKFLines==mpCurrentKeyFrame->mnId || pKFi2->mnId==mpCurrentKeyFrame->mnId)
                continue;
            vpTargetKFs.push_back(pKFi2);
        }
    }


    // Search matches by projection from current KF in target KFs
    Linematcher matcher(0.9, false, false, 0.25);
    vector<MapLine*> vpMapLineMatches = mpCurrentKeyFrame->GetMapLineMatches();
    for(vector<KeyFrame*>::iterator vit=vpTargetKFs.begin(), vend=vpTargetKFs.end(); vit!=vend; vit++)
    {
        KeyFrame* pKFi = *vit;

        matcher.Fuse(pKFi,vpMapLineMatches);
    }

    // Search matches by projection from target KFs in current KF
    vector<MapLine*> vpFuseCandidates;
    vpFuseCandidates.reserve(vpTargetKFs.size()*vpMapLineMatches.size());

    for(vector<KeyFrame*>::iterator vitKF=vpTargetKFs.begin(), vendKF=vpTargetKFs.end(); vitKF!=vendKF; vitKF++)
    {
        KeyFrame* pKFi = *vitKF;

        vector<MapLine*> vpMapLinesKFi = pKFi->GetMapLineMatches();

        for(vector<MapLine*>::iterator vitML=vpMapLinesKFi.begin(), vendML=vpMapLinesKFi.end(); vitML!=vendML; vitML++)
        {
            MapLine* pML = *vitML;
            if(!pML)
                continue;
            if(pML->isBad() || pML->mnFuseCandidateForKF == mpCurrentKeyFrame->mnId)
                continue;
            pML->mnFuseCandidateForKF = mpCurrentKeyFrame->mnId;
            vpFuseCandidates.push_back(pML);
        }
    }

    matcher.Fuse(mpCurrentKeyFrame,vpFuseCandidates);


    // Update Lines
    vpMapLineMatches = mpCurrentKeyFrame->GetMapLineMatches();
    for(size_t i=0, iend=vpMapLineMatches.size(); i<iend; i++)
    {
        MapLine* pML = vpMapLineMatches[i];
        if(pML)
        {
            if(!pML->isBad())
            {
                pML->ComputeDistinctiveDescriptors();
                pML->UpdateNormalAndDepth();
		pML->Update2DLineLength();
            }
        }
    }

    mpCurrentKeyFrame->UpdateConnectionsLines();
}

cv::Mat LocalMapping::ComputeF12(KeyFrame *&pKF1, KeyFrame *&pKF2)
{
    cv::Mat R1w = pKF1->GetRotation();
    cv::Mat t1w = pKF1->GetTranslation();
    cv::Mat R2w = pKF2->GetRotation();
    cv::Mat t2w = pKF2->GetTranslation();

    cv::Mat R12 = R1w*R2w.t();
    cv::Mat t12 = -R1w*R2w.t()*t2w+t1w;

    cv::Mat t12x = SkewSymmetricMatrix(t12);

    const cv::Mat &K1 = pKF1->mK;
    const cv::Mat &K2 = pKF2->mK;


    return K1.t().inv()*t12x*R12*K2.inv();
}

void LocalMapping::RequestStop()
{
    unique_lock<mutex> lock(mMutexStop);
    mbStopRequested = true;
    unique_lock<mutex> lock2(mMutexNewKFs);
    mbAbortBA = true;
}

bool LocalMapping::Stop()
{
    unique_lock<mutex> lock(mMutexStop);
    if(mbStopRequested && !mbNotStop)
    {
        mbStopped = true;
        cout << "Local Mapping STOP" << endl;
        return true;
    }

    return false;
}

bool LocalMapping::isStopped()
{
    unique_lock<mutex> lock(mMutexStop);
    return mbStopped;
}

bool LocalMapping::stopRequested()
{
    unique_lock<mutex> lock(mMutexStop);
    return mbStopRequested;
}

void LocalMapping::Release()
{
    unique_lock<mutex> lock(mMutexStop);
    unique_lock<mutex> lock2(mMutexFinish);
    if(mbFinished)
        return;
    mbStopped = false;
    mbStopRequested = false;
    for(list<KeyFrame*>::iterator lit = mlNewKeyFrames.begin(), lend=mlNewKeyFrames.end(); lit!=lend; lit++)
        delete *lit;
    mlNewKeyFrames.clear();

    cout << "Local Mapping RELEASE" << endl;
}

bool LocalMapping::AcceptKeyFrames()
{
    unique_lock<mutex> lock(mMutexAccept);
    return mbAcceptKeyFrames;
}

void LocalMapping::SetAcceptKeyFrames(bool flag)
{
    unique_lock<mutex> lock(mMutexAccept);
    mbAcceptKeyFrames=flag;
}

bool LocalMapping::SetNotStop(bool flag)
{
    unique_lock<mutex> lock(mMutexStop);

    if(flag && mbStopped)
        return false;

    mbNotStop = flag;

    return true;
}

void LocalMapping::InterruptBA()
{
    mbAbortBA = true;
}

void LocalMapping::KeyFrameCulling()
{
    // Check redundant keyframes (only local keyframes)
    // A keyframe is considered redundant if the 90% of the MapPoints it sees, are seen
    // in at least other 3 keyframes (in the same or finer scale)
    // We only consider close stereo points
    vector<KeyFrame*> vpLocalKeyFrames = mpCurrentKeyFrame->GetVectorCovisibleKeyFrames();

    for(vector<KeyFrame*>::iterator vit=vpLocalKeyFrames.begin(), vend=vpLocalKeyFrames.end(); vit!=vend; vit++)
    {
        KeyFrame* pKF = *vit;
        if(pKF->mnId==0)
            continue;
        const vector<MapPoint*> vpMapPoints = pKF->GetMapPointMatches();

        int nObs = 3;
        const int thObs=nObs;
        int nRedundantObservations=0;
        int nMPs=0;
        for(size_t i=0, iend=vpMapPoints.size(); i<iend; i++)
        {
            MapPoint* pMP = vpMapPoints[i];
            if(pMP)
            {
                if(!pMP->isBad())
                {
                    if(!mbMonocular)
                    {
                        if(pKF->mvDepth[i]>pKF->mThDepth || pKF->mvDepth[i]<0)
                            continue;
                    }

                    nMPs++;
                    if(pMP->Observations()>thObs)
                    {
                        const int &scaleLevel = pKF->mvKeysUn[i].octave;
                        const map<KeyFrame*, size_t> observations = pMP->GetObservations();
                        int nObs=0;
                        for(map<KeyFrame*, size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
                        {
                            KeyFrame* pKFi = mit->first;
                            if(pKFi==pKF)
                                continue;
                            const int &scaleLeveli = pKFi->mvKeysUn[mit->second].octave;

                            if(scaleLeveli<=scaleLevel+1)
                            {
                                nObs++;
                                if(nObs>=thObs)
                                    break;
                            }
                        }
                        if(nObs>=thObs)
                        {
                            nRedundantObservations++;
                        }
                    }
                }
            }
        }  

        if(nRedundantObservations>0.9*nMPs)
            pKF->SetBadFlag();
    }
}

//剔除策略的主函数，执行完下面两个子函数后，分三种情况：点线均有共视关键帧、仅点有关键帧和仅线有关键帧
void LocalMapping::KeyFrameCullingBoth()
{	
    thread threadKeyFrameCullingPoints(&LocalMapping::KeyFrameCullingPoints, this );
    thread threadKeyFrameCullingLines(&LocalMapping::KeyFrameCullingLines, this );
    threadKeyFrameCullingPoints.join();
    threadKeyFrameCullingLines.join();
    
    set<KeyFrame*> spLocalKeyFramesPoints = mpCurrentKeyFrame->GetConnectedKeyFrames();
    set<KeyFrame*> spLocalKeyFramesLines = mpCurrentKeyFrame->GetConnectedKeyFramesLines();
    if(!spLocalKeyFramesPoints.empty() && !spLocalKeyFramesLines.empty())
    {
      
	for(set<KeyFrame*>::iterator vit=spLocalKeyFramesPoints.begin(), vend=spLocalKeyFramesPoints.end(); vit!=vend; vit++)
	{
	    KeyFrame* pKF = *vit;
	    if(pKF->isBad())
	    {
		if(spLocalKeyFramesLines.count(pKF) && pKF->isBadLines())
		{
		    mpMap->EraseKeyFrame(pKF);
		}
	    }
	}
    }
    else if(spLocalKeyFramesLines.empty())
    {
      
	 for(set<KeyFrame*>::iterator vit=spLocalKeyFramesPoints.begin(), vend=spLocalKeyFramesPoints.end(); vit!=vend; vit++)
	 {
	     KeyFrame* pKF = *vit;
	     if(pKF->isBad())
	     {
		 mpMap->EraseKeyFrame(pKF);
	     }
	 }
    }
    else if(spLocalKeyFramesPoints.empty())
    {
	for(set<KeyFrame*>::iterator vit=spLocalKeyFramesLines.begin(), vend=spLocalKeyFramesLines.end(); vit!=vend; vit++)
	{
	    KeyFrame* pKF = *vit;
	    if(pKF->isBadLines())
	    {
		 mpMap->EraseKeyFrame(pKF);
	    }
	}	
    }
    
}

//（分线程，达到90%共视，同原点特征的关键帧剔除函数，这里可以直接删除用于重定位和回环检测的关键帧集中的关键帧）
void LocalMapping::KeyFrameCullingPoints()
{
    // Check redundant keyframes (only local keyframes)
    // A keyframe is considered redundant if the 90% of the MapPoints it sees, are seen
    // in at least other 3 keyframes (in the same or finer scale)
    vector<KeyFrame*> vpLocalKeyFrames = mpCurrentKeyFrame->GetVectorCovisibleKeyFrames();

    for(vector<KeyFrame*>::iterator vit=vpLocalKeyFrames.begin(), vend=vpLocalKeyFrames.end(); vit!=vend; vit++)
    {
        KeyFrame* pKF = *vit;
        if(pKF->mnId==0)
            continue;
        const vector<MapPoint*> vpMapPoints = pKF->GetMapPointMatches();

        int nObs = 3;
        const int thObs=nObs;
        int nRedundantObservations=0;
        int nMPs=0;
        for(size_t i=0, iend=vpMapPoints.size(); i<iend; i++)
        {
            MapPoint* pMP = vpMapPoints[i];
            if(pMP)
            {
                if(!pMP->isBad())
                {
                    if(!mbMonocular)
                    {
                        if(pKF->mvDepth[i]>pKF->mThDepth || pKF->mvDepth[i]<0)
                            continue;
                    }

                    nMPs++;
                    if(pMP->Observations()>thObs)
                    {
                        const int &scaleLevel = pKF->mvKeysUn[i].octave;
                        const map<KeyFrame*, size_t> observations = pMP->GetObservations();
                        int nObs=0;
                        for(map<KeyFrame*, size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
                        {
                            KeyFrame* pKFi = mit->first;
                            if(pKFi==pKF)
                                continue;
                            const int &scaleLeveli = pKFi->mvKeysUn[mit->second].octave;

                            if(scaleLeveli<=scaleLevel+1)
                            {
                                nObs++;
                                if(nObs>=thObs)
                                    break;
                            }
                        }
                        if(nObs>=thObs)
                        {
                            nRedundantObservations++;
                        }
                    }
                }
            }
        }  

        if(nRedundantObservations>0.9*nMPs)
            pKF->SetBadFlagPoints();
    }
}

//（分线程，达到90%共视）
void LocalMapping::KeyFrameCullingLines()
{
    // Check redundant keyframes (only local keyframes)
    // A keyframe is considered redundant if the 90% of the MapLines it sees, are seen
    // in at least other 3 keyframes (in the same or finer scale)
    vector<KeyFrame*> vpLocalKeyFrames = mpCurrentKeyFrame->GetVectorCovisibleKeyFramesLines();

    for(vector<KeyFrame*>::iterator vit=vpLocalKeyFrames.begin(), vend=vpLocalKeyFrames.end(); vit!=vend; vit++)
    {
        KeyFrame* pKF = *vit;
        if(pKF->mnId==0)
            continue;
        const vector<MapLine*> vpMapLines = pKF->GetMapLineMatches();

        int nObs = 3;
        const int thObs=nObs;
        int nRedundantObservations=0;
        int nMLs=0;
        for(size_t i=0, iend=vpMapLines.size(); i<iend; i++)
        {
            MapLine* pML = vpMapLines[i];
            if(pML)
            {
                if(!pML->isBad())
                {
                    nMLs++;
                    if(pML->Observations()>thObs)
                    {
                        const int &scaleLevel = pKF->mvLinesUn[i].octave;
                        const map<KeyFrame*, size_t> observations = pML->GetObservations();
                        int nObs=0;
                        for(map<KeyFrame*, size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
                        {
                            KeyFrame* pKFi = mit->first;
                            if(pKFi==pKF)
                                continue;
                            const int &scaleLeveli = pKFi->mvLinesUn[mit->second].octave;

                            if(scaleLeveli<=scaleLevel+1)
                            {
                                nObs++;
                                if(nObs>=thObs)
                                    break;
                            }
                        }
                        if(nObs>=thObs)
                        {
                            nRedundantObservations++;
                        }
                    }
                }
            }
        }  

        if(nRedundantObservations>0.9*nMLs)
            pKF->SetBadFlagLines();
    }
}

cv::Mat LocalMapping::SkewSymmetricMatrix(const cv::Mat &v)
{
    return (cv::Mat_<float>(3,3) <<             0, -v.at<float>(2), v.at<float>(1),
            v.at<float>(2),               0,-v.at<float>(0),
            -v.at<float>(1),  v.at<float>(0),              0);
}

void LocalMapping::RequestReset()
{
    {
        unique_lock<mutex> lock(mMutexReset);
        mbResetRequested = true;
    }

    while(1)
    {
        {
            unique_lock<mutex> lock2(mMutexReset);
            if(!mbResetRequested)
                break;
        }
        usleep(3000);
    }
}

void LocalMapping::ResetIfRequested()
{
    unique_lock<mutex> lock(mMutexReset);
    if(mbResetRequested)
    {
        mlNewKeyFrames.clear();
        mlpRecentAddedMapPoints.clear();
        mbResetRequested=false;
    }
}

//使用点线特征的重置函数
void LocalMapping::ResetIfRequestedBoth()
{
    unique_lock<mutex> lock(mMutexReset);
    if(mbResetRequested)
    {
        mlNewKeyFrames.clear();
        mlpRecentAddedMapPoints.clear();
	mlpRecentAddedMapLines.clear();
        mbResetRequested=false;
    }
}

void LocalMapping::RequestFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    mbFinishRequested = true;
}

bool LocalMapping::CheckFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinishRequested;
}

void LocalMapping::SetFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    mbFinished = true;    
    unique_lock<mutex> lock2(mMutexStop);
    mbStopped = true;
}

bool LocalMapping::isFinished()
{
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinished;
}

} //namespace ORB_SLAM
