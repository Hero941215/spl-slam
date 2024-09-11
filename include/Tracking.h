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


#ifndef TRACKING_H
#define TRACKING_H

#include<opencv2/core/core.hpp>
#include<opencv2/features2d/features2d.hpp>

#include"Viewer.h"
#include"FrameDrawer.h"
#include"Map.h"
#include"LocalMapping.h"
#include"LoopClosing.h"
#include"Frame.h"
#include "ORBVocabulary.h"
#include"KeyFrameDatabase.h"
#include"ORBextractor.h"
#include"Lineextractor.h"//(LineExpanding)
#include"ORBmatcher.h"
#include"Linematcher.h" //(LineExpanding)
#include "Initializer.h"
#include "MapDrawer.h"
#include "System.h"

#include <mutex>

namespace PL_SLAM
{

class Viewer;
class FrameDrawer;
class Map;
class LocalMapping;
class LoopClosing;
class System;

class Tracking
{  

public:
    //添加输入参数：是否使用线特征
    Tracking(System* pSys, ORBVocabulary* pVoc, FrameDrawer* pFrameDrawer, MapDrawer* pMapDrawer, Map* pMap,
             KeyFrameDatabase* pKFDB, const string &strSettingPath, const int sensor);

    // Preprocess the input and call Track(). Extract features and performs stereo matching.
    cv::Mat GrabImageStereo(const cv::Mat &imRectLeft,const cv::Mat &imRectRight, const double &timestamp);
    cv::Mat GrabImageRGBD(const cv::Mat &imRGB,const cv::Mat &imD, const double &timestamp);
    cv::Mat GrabImageMonocular(const cv::Mat &im, const double &timestamp);

    void SetLocalMapper(LocalMapping* pLocalMapper);
    void SetLoopClosing(LoopClosing* pLoopClosing);
    void SetViewer(Viewer* pViewer);

    // Load new settings
    // The focal lenght should be similar or scale prediction will fail when projecting points
    // TODO: Modify MapPoint::PredictScale to take into account focal lenght
    void ChangeCalibration(const string &strSettingPath);

    // Use this function if you have deactivated local mapping and you only want to localize the camera.
    void InformOnlyTracking(const bool &flag);


public:

    // Tracking states
    enum eTrackingState{
        SYSTEM_NOT_READY=-1,
        NO_IMAGES_YET=0,
        NOT_INITIALIZED=1,
        OK=2,
        LOST=3
    };

    eTrackingState mState;
    eTrackingState mLastProcessedState;

    // Input sensor
    int mSensor;

    bool mbusingLine;
    bool mbusingLsdFeature;

    // Current Frame
    Frame mCurrentFrame;
    cv::Mat mImGray;

    // Initialization Variables (Monocular)
    std::vector<int> mvIniLastMatches;//未使用
    
    std::vector<int> mvIniMatches;
    std::vector<int> mvIniMatchesLines;
    
    std::vector<cv::Point2f> mvbPrevMatched;
    std::vector<cv::Point2f> mvbPrevMatchedMidPoints;
    
    std::vector<cv::Point3f> mvIniP3D;
    std::vector<cv::Point3f> mvIniFirP3D;
    std::vector<cv::Point3f> mvIniEndP3D;
    std::vector<cv::Point3f> mvIniMidP3D;//这是最先恢复的线特征3D点
    
    Frame mInitialFrame;

    // Lists used to recover the full camera trajectory at the end of the execution.
    // Basically we store the reference keyframe for each frame and its relative transformation
    list<cv::Mat> mlRelativeFramePoses;
    list<KeyFrame*> mlpReferences;
    list<double> mlFrameTimes;
    list<bool> mlbLost;

    // True if local mapping is deactivated and we are performing only localization
    bool mbOnlyTracking;

    void Reset();
    void ResetBoth();

protected:
    

    // Main tracking function. It is independent of the input sensor.
    void Track();
    void TrackBoth (); 
    
    // Map initialization for stereo and RGB-D（目前还未扩展）
    void StereoInitialization();

    //单目初始化：
    // Map initialization for monocular
    void MonocularInitialization();
    void MonocularInitializationBoth();
    void CreateInitialMapMonocular();
    void CreateInitialMapMonocularBoth ();
    
    //检查上一帧图像中点线特征是否被后端替换：
    void CheckReplacedInLastFrame();
    void CheckReplacedInLastFrameLines ();
    
    //当前图像帧的初步位姿估计
    bool TrackReferenceKeyFrame();
    bool TrackReferenceKeyFrameBoth ();
    
    void UpdateLastFrame();//点特征与点线特征通用
    bool TrackWithMotionModel();
    bool TrackWithMotionModelBoth ();

    //重定位：
    bool Relocalization();
    bool RelocalizationBoth();
    //第二次点线同时优化结束后，找到的有效信息又出现少于50的情况（三种），分情况在投影匹配查找和优化
    void RelocalizationBothTwiceSearch(KeyFrame* vpCandidateKFs, int &inlinersNum, int &inlinerLinesNum );
    //初步位姿估计结束，基于位姿计算当前图像帧的有效地图点匹配和内点数量
    int SetCurrentFrameMappointsAndInliers(Frame &CurrentFrame, std::vector<MapPoint*> vpMapPointMatches, std::set<MapPoint*> &sFound);

    //初步位姿估计后，为了提高估计精度，进一步使用当前图像帧的局部地图进行跟踪
    void UpdateLocalMap();
    void UpdateLocalPoints();
    void UpdateLocalKeyFrames();
    void UpdateLocalMapLines();
    void UpdateLocalKeyFramesLines();
    void UpdateLocalLines();
    
    //搜索当前帧与局部地图点和线中点的匹配
    void SearchLocalPoints();
    void SearchLocalPoints2(int &addtionPointsNum);
    void SearchLocalLines(int &addtionLinesNum);
    
    bool TrackLocalMap();
    bool TrackLocalMapBoth();
    //用于在地图线匹配数量较少或者没有时，通过点共视图对地图线信息进行补充
    int MapLineRenewing();

    //判断是否需要关键帧：
    bool NeedNewKeyFrame();
    bool NeedNewKeyFrameBoth ();
    
    //创建新的关键帧：
    void CreateNewKeyFrame();
    void CreateNewKeyFrameBoth();
    
    //调用前端匹配器：
    //点匹配器调用：
    void SearchForInitialization(Frame &F1, Frame &F2, std::vector<cv::Point2f> &vbPrevMatched, std::vector<int> &vnMatches12,
				 int windowSize, int &nmatches);
    
    void SearchByProjectionMotion(Frame &CurrentFrame, const Frame &LastFrame, const float th,
			    const bool bMono, int &nmatches);
    
    void SearchByBoW(KeyFrame *pKF, Frame &F, std::vector<MapPoint*> &vpMapPointMatches, int &nmatches);

    void SearchByProjection(Frame &CurrentFrame, KeyFrame* pKF, const std::set<MapPoint*> &sAlreadyFound,
			    const float th, const int ORBdist, int &nmatches);
    
    //线匹配器调用
    void SearchForInitializationLines(Frame &F1, Frame &F2, std::vector<cv::Point2f> &vbPrevMatched, std::vector<int> &vnMatches12,
				 int windowSize, int &nmatches);
  
    void SearchByProjectionLinesMotion(Frame &CurrentFrame, const Frame &LastFrame, const float th,
			    const bool bMono, int &nmatches);
    
    void SearchByKNNLines(KeyFrame *pKF, Frame &F, std::vector<MapLine*> &vpMapLineMatches, int &nmatches12);

    void SearchByProjectionLines(Frame &CurrentFrame, KeyFrame* pKF, const std::set<MapLine*> &sAlreadyFound, 
			    const float th, const int Linedist, int &nmatches);

    // In case of performing only localization, this flag is true when there are no matches to
    // points in the map. Still tracking will continue if there are enough matches with temporal points.
    // In that case we are doing visual odometry. The system will try to do relocalization to recover
    // "zero-drift" localization to the map.
    bool mbVO;

    //Other Thread Pointers
    LocalMapping* mpLocalMapper;
    LoopClosing* mpLoopClosing;

    //ORB
    ORBextractor* mpORBextractorLeft, *mpORBextractorRight;
    ORBextractor* mpIniORBextractor;
    //线特征提取器：
    Lineextractor* mpLineextractorLeft;
    Lineextractor* mpIniLineextractor;
    
    //匹配器：
    ORBmatcher* mpORBmatcherLeftIni;
    ORBmatcher* mpORBmatcherLeftmotion;
    ORBmatcher* mpORBmatcherLeftBoW;
    ORBmatcher* mpORBmatcherLeftRel;
    //线匹配器指针：
    Linematcher* mpLinematcherLeftIni;
    Linematcher* mpLinematcherLeftmotion;
    Linematcher* mpLinematcherLeftKNN;
    Linematcher* mpLinematcherLeftRel;
    
    //BoW
    ORBVocabulary* mpORBVocabulary;
    KeyFrameDatabase* mpKeyFrameDB;

    // Initalization (only for monocular)
    Initializer* mpInitializer;

    //Local Map
    KeyFrame* mpReferenceKF;
    KeyFrame* mpReferenceKFLines;//该关键帧为与当前图像帧有最多匹配线的关键帧
    
    std::vector<KeyFrame*> mvpLocalKeyFrames;
    std::vector<MapPoint*> mvpLocalMapPoints;
    std::vector<KeyFrame*> mvpLocalKeyFramesLines;
    std::vector<MapLine*> mvpLocalMapLines;
    
    // System
    System* mpSystem;
    
    //Drawers
    Viewer* mpViewer;
    FrameDrawer* mpFrameDrawer;
    MapDrawer* mpMapDrawer;

    //Map
    Map* mpMap;

    //Calibration matrix
    cv::Mat mK;
    cv::Mat mDistCoef;
    float mbf;

    //New KeyFrame rules (according to fps)
    int mMinFrames;
    int mMaxFrames;

    // Threshold close/far points
    // Points seen as close by the stereo/RGBD sensor are considered reliable
    // and inserted from just one frame. Far points requiere a match in two keyframes.
    float mThDepth;

    // For RGB-D inputs only. For some datasets (e.g. TUM) the depthmap values are scaled.
    float mDepthMapFactor;

    //时间记录变量
    double TrackingAllTime;
    double TrackingAverageTime;
    int TrackingNumbers;
    
    //Current matches in frame
    int mnMatchesInliers;
    int mnMatchesInlierLines;

    //Last Frame, KeyFrame and Relocalisation Info
    KeyFrame* mpLastKeyFrame;
    Frame mLastFrame;
    unsigned int mnLastKeyFrameId;
    unsigned int mnLastRelocFrameId;

    //Motion Model
    cv::Mat mVelocity;

    //Color order (true RGB, false BGR, ignored if grayscale)
    bool mbRGB;

    list<MapPoint*> mlpTemporalPoints;
};

} //namespace ORB_SLAM

#endif // TRACKING_H
