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

#ifndef LOCALMAPPING_H
#define LOCALMAPPING_H

#include "KeyFrame.h"
#include "Map.h"
#include "LoopClosing.h"
#include "Tracking.h"
#include "KeyFrameDatabase.h"

#include <list>

#include <mutex>


namespace PL_SLAM
{

class Tracking;
class LoopClosing;
class Map;

class LocalMapping
{
public:
    LocalMapping(Map* pMap, const float bMonocular);

    void SetLoopCloser(LoopClosing* pLoopCloser);

    void SetTracker(Tracking* pTracker);

    // Main function
    void Run();
    void RunBoth();

    void InsertKeyFrame(KeyFrame* pKF);

    // Thread Synch
    void RequestStop();
    void RequestReset();
    bool Stop();
    void Release();
    bool isStopped();
    bool stopRequested();
    bool AcceptKeyFrames();
    void SetAcceptKeyFrames(bool flag);
    bool SetNotStop(bool flag);

    void InterruptBA();

    void RequestFinish();
    bool isFinished();

    int KeyframesInQueue(){
        unique_lock<std::mutex> lock(mMutexNewKFs);
        return mlNewKeyFrames.size();
    }

protected:

    bool CheckNewKeyFrames();
    
    //当前关键帧的3D动态点线特征检测：(以后的工作)
    void DynamicPointsDetecting();
    void DynamicLinesDetecting();
    
    //处理当前帧，更新点线特征观测信息和对应共视图信息，在地图中插入关键帧：
    void ProcessNewKeyFrame();
    void ProcessNewKeyFrameBoth();
    void ProcessNewKeyFramePoints();
    void ProcessNewKeyFrameLines();
    
    
    //创建新的地图点线：
    void CreateNewMapPoints();
    void CreateNewMapLines();
    void CreateNewMapLinesWhenPointsOnly();
    
    //临时3D点线特征信息剔除：
    void MapPointCulling();
    void MapLineCulling();
    
    //当前关键帧与相邻帧（两级相邻）3D点线特征融合：
    void SearchInNeighbors();
    void SearchInNeighborsLines();

    //冗余关键帧剔除策略：
    void KeyFrameCulling();
    void KeyFrameCullingBoth();
    void KeyFrameCullingPoints();
    void KeyFrameCullingLines();
    
    cv::Mat ComputeF12(KeyFrame* &pKF1, KeyFrame* &pKF2);

    cv::Mat SkewSymmetricMatrix(const cv::Mat &v);

    bool mbMonocular;

    void ResetIfRequested();
    void ResetIfRequestedBoth();
    
    bool mbResetRequested;
    std::mutex mMutexReset;

    bool CheckFinish();
    void SetFinish();
    bool mbFinishRequested;
    bool mbFinished;
    std::mutex mMutexFinish;

    Map* mpMap;

    LoopClosing* mpLoopCloser;
    Tracking* mpTracker;

    std::list<KeyFrame*> mlNewKeyFrames;

    KeyFrame* mpCurrentKeyFrame;

    std::list<MapPoint*> mlpRecentAddedMapPoints;
    std::list<MapLine*> mlpRecentAddedMapLines;
    
    std::mutex mMutexNewKFs;

    bool mbAbortBA;

    bool mbStopped;
    bool mbStopRequested;
    bool mbNotStop;
    std::mutex mMutexStop;

    bool mbAcceptKeyFrames;
    
    //时间记录变量
    double localBAAllTime;
    double localBAAverageTime;
    int localBANumbers;
    
    
    std::mutex mMutexAccept;
};

} //namespace ORB_SLAM

#endif // LOCALMAPPING_H
