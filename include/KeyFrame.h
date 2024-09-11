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

#ifndef KEYFRAME_H
#define KEYFRAME_H

#include "MapPoint.h"
#include "MapLine.h"
#include "Thirdparty/DBoW2/DBoW2/BowVector.h"
#include "Thirdparty/DBoW2/DBoW2/FeatureVector.h"
#include "ORBVocabulary.h"
#include "ORBextractor.h"
#include "Lineextractor.h"
#include "Frame.h"
#include "KeyFrameDatabase.h"

#include <mutex>


namespace PL_SLAM
{

class Map;
class MapPoint;
class MapLine;
class Frame;
class KeyFrameDatabase;

class KeyFrame
{
public:
    KeyFrame(Frame &F, Map* pMap, KeyFrameDatabase* pKFDB);

    // Pose functions
    void SetPose(const cv::Mat &Tcw);
    void SetPosePoints(const cv::Mat &TcwPoints);
    void SetPoseLines(const cv::Mat &TcwLines);
    
    cv::Mat GetPose();
    cv::Mat GetPoseInverse();
    cv::Mat GetCameraCenter();
    cv::Mat GetStereoCenter();
    cv::Mat GetRotation();
    cv::Mat GetTranslation();

    // Bag of Words Representation
    void ComputeBoW();

    // Covisibility graph functions
    void AddConnection(KeyFrame* pKF, const int &weight);
    void EraseConnection(KeyFrame* pKF);
    void UpdateConnections();
    void UpdateBestCovisibles();
    std::set<KeyFrame *> GetConnectedKeyFrames();
    std::vector<KeyFrame* > GetVectorCovisibleKeyFrames();
    std::vector<KeyFrame*> GetBestCovisibilityKeyFrames(const int &N);
    std::vector<KeyFrame*> GetCovisiblesByWeight(const int &w);
    int GetWeight(KeyFrame* pKF);
    //线特征的共视图函数：
    void AddConnectionLines(KeyFrame* pKF, const int &weight);
    void EraseConnectionLines(KeyFrame* pKF);
    void UpdateConnectionsLines();
    void UpdateBestCovisiblesLines();
    std::set<KeyFrame *> GetConnectedKeyFramesLines();
    std::vector<KeyFrame* > GetVectorCovisibleKeyFramesLines();
    std::vector<KeyFrame*> GetBestCovisibilityKeyFramesLines(const int &N);
    std::vector<KeyFrame*> GetCovisiblesByWeightLines(const int &w);
    int GetWeightLines(KeyFrame* pKF);
    
    // Spanning tree functions
    void AddChild(KeyFrame* pKF);
    void EraseChild(KeyFrame* pKF);
    void ChangeParent(KeyFrame* pKF);
    std::set<KeyFrame*> GetChilds();
    KeyFrame* GetParent();
    bool hasChild(KeyFrame* pKF);
    //线特征的最小生成树：
    void AddChildLines(KeyFrame* pKF);
    void EraseChildLines(KeyFrame* pKF);
    void ChangeParentLine(KeyFrame* pKF);
    std::set<KeyFrame*> GetChildsLines();
    KeyFrame* GetParentLines();
    bool hasChildLines(KeyFrame* pKF);
    
    // Loop Edges
    void AddLoopEdge(KeyFrame* pKF);
    std::set<KeyFrame*> GetLoopEdges();

    // MapPoint observation functions
    void AddMapPoint(MapPoint* pMP, const size_t &idx);
    void EraseMapPointMatch(const size_t &idx);
    void EraseMapPointMatch(MapPoint* pMP);
    void ReplaceMapPointMatch(const size_t &idx, MapPoint* pMP);
    std::set<MapPoint*> GetMapPoints();
    std::vector<MapPoint*> GetMapPointMatches();
    int TrackedMapPoints(const int &minObs);
    MapPoint* GetMapPoint(const size_t &idx);
    //地图线观测函数：//(LineExpanding)
    void AddMapLine(MapLine* pML, const size_t &idx);
    void EraseMapLineMatch(const size_t &idx);
    void EraseMapLineMatch(MapLine* pML);
    void ReplaceMapLineMatch(const size_t &idx, MapLine* pML);
    std::set<MapLine*> GetMapLines();
    std::vector<MapLine*> GetMapLineMatches();
    int TrackedMapLines(const int &minObs);
    MapLine* GetMapLine(const size_t &idx);
    
    // KeyPoint functions
    std::vector<size_t> GetFeaturesInArea(const float &x, const float  &y, const float  &r) const;
    std::vector<size_t> GetFeaturesInAreaLines(const float &x, const float  &y, const float  &r) const;
    cv::Mat UnprojectStereo(int i);

    // Image
    bool IsInImage(const float &x, const float &y) const;

    // Enable/Disable bad flag changes
    void SetNotErase();
    void SetErase();

    // Set/check bad flag
    void SetBadFlag();
    bool isBad();
    bool isBadLines();
    
    void SetBadFlagPoints();
    void SetBadFlagLines();
    
    // Compute Scene Depth (q=2 median). Used in monocular.
    float ComputeSceneMedianDepth(const int q);
    float ComputeSceneMedianDepthBoth(const int q);

    static bool weightComp( int a, int b){
        return a>b;
    }

    static bool lId(KeyFrame* pKF1, KeyFrame* pKF2){
        return pKF1->mnId<pKF2->mnId;
    }


    // The following variables are accesed from only 1 thread or never change (no mutex needed).
public:

    static long unsigned int nNextId;
    long unsigned int mnId;
    const long unsigned int mnFrameId;

    const double mTimeStamp;

    // Grid (to speed up feature matching)
    const int mnGridCols;
    const int mnGridRows;
    const float mfGridElementWidthInv;
    const float mfGridElementHeightInv;
    //图像线特征分块记录：
    int mnGridColsLines;
    int mnGridRowsLines;
    float mfGridElementWidthInvLines;
    float mfGridElementHeightInvLines;
    
    // Variables used by the tracking
    long unsigned int mnTrackReferenceForFrame;
    long unsigned int mnFuseTargetForKF;

    // Variables used by the local mapping
    long unsigned int mnBALocalForKF;
    long unsigned int mnBAFixedForKF;
    
    long unsigned int mnTrackReferenceForFrameLines;
    long unsigned int mnFuseTargetForKFLines;
    long unsigned int mnBALocalForKFLines;
    long unsigned int mnBAFixedForKFLines;
    
    long unsigned int mBAPoseLinesandPointsForKF;
    long unsigned int mnBALocalForKFBoth;
    long unsigned int mnBAFixedForKFBoth;
    
    // Variables used by the keyframe database
    long unsigned int mnLoopQuery;
    int mnLoopWords;
    float mLoopScore;
    long unsigned int mnRelocQuery;
    int mnRelocWords;
    float mRelocScore;

    // Variables used by loop closing
    cv::Mat mTcwGBA;
    cv::Mat mTcwBefGBA;
    long unsigned int mnBAGlobalForKF;

    // Calibration parameters
    const float fx, fy, cx, cy, invfx, invfy, mbf, mb, mThDepth;

    const int N;
    int NL;
    
    // KeyPoints, stereo coordinate and descriptors (all associated by an index)
    const std::vector<cv::KeyPoint> mvKeys;
    const std::vector<cv::KeyPoint> mvKeysUn;
    const std::vector<float> mvuRight; // negative value for monocular points
    const std::vector<float> mvDepth; // negative value for monocular points
    const cv::Mat mDescriptors;
    

    std::vector<KeyLine> mvLines; 
    std::vector<cv::KeyPoint> mvMidPoints;
    std::vector<KeyLine> mvLinesUn; 
    std::vector<cv::KeyPoint>  mvMidPointsUn;
    cv::Mat mDescriptorLines;
    
    //BoW
    DBoW2::BowVector mBowVec;
    DBoW2::FeatureVector mFeatVec;

    // Pose relative to parent (this is computed when bad flag is activated)
    cv::Mat mTcp;
    cv::Mat mTcpLines;

    // Scale
    const int mnScaleLevels;
    const float mfScaleFactor;
    const float mfLogScaleFactor;
    const std::vector<float> mvScaleFactors;
    const std::vector<float> mvLevelSigma2;
    const std::vector<float> mvInvLevelSigma2;
    //线特征提取的尺度金字塔信息：(LineExpanding)
    int mnScaleLevelLines;
    float mfScaleFactorLines;
    float mfLogScaleFactorsLines;
    vector<float> mvScaleFactorsLines;
    vector<float> mvInvScaleFactorsLines;
    vector<float> mvLevelSigma2Lines;
    vector<float> mvInvLevelSigma2Lines;

    // Image bounds and calibration
    const int mnMinX;
    const int mnMinY;
    const int mnMaxX;
    const int mnMaxY;
    const cv::Mat mK;
    
    cv::Mat TcwPoints;
    cv::Mat TcwLines;

    // The following variables need to be accessed trough a mutex to be thread safe.
protected:

    // SE3 Pose and camera center
    cv::Mat Tcw;
    cv::Mat Twc;
    cv::Mat Ow;

    cv::Mat Cw; // Stereo middel point. Only for visualization

    // MapPoints associated to keypoints
    std::vector<MapPoint*> mvpMapPoints;
    //地图线：
    std::vector<MapLine*> mvpMapLines;
    
    // BoW
    KeyFrameDatabase* mpKeyFrameDB;
    ORBVocabulary* mpORBvocabulary;

    // Grid over the image to speed up feature matching
    std::vector< std::vector <std::vector<size_t> > > mGrid;
    std::vector< std::vector <std::vector<size_t> > > mGridLines;

    std::map<KeyFrame*,int> mConnectedKeyFrameWeights;
    std::vector<KeyFrame*> mvpOrderedConnectedKeyFrames;
    std::vector<int> mvOrderedWeights;
    //线特征共视图：
    std::map<KeyFrame*,int> mConnectedKeyFrameWeightsLines; 
    std::vector<KeyFrame*> mvpOrderedConnectedKeyFramesLines; 
    std::vector<int> mvOrderedWeightsLines; 

    // Spanning Tree and Loop Edges
    bool mbFirstConnection;
    KeyFrame* mpParent;
    std::set<KeyFrame*> mspChildrens;
    std::set<KeyFrame*> mspLoopEdges;

    bool mbFirstConnectionLines;
    KeyFrame* mpParentLines;
    std::set<KeyFrame*> mspChildrensLines;
    
    // Bad flags
    bool mbNotErase;
    bool mbToBeErased;
    bool mbBad;    
    bool mbBadLines;
    
    float mHalfBaseline; // Only for visualization

    Map* mpMap;

    std::mutex mMutexPose;
    std::mutex mMutexConnections;
    std::mutex mMutexFeatures;
    std::mutex mMutexConnectionsLines;
    std::mutex mMutexFeaturesLines;
};

} //namespace ORB_SLAM

#endif // KEYFRAME_H
