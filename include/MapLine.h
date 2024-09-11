/**
*  This file is part of SPL-SLAM.
* 
* 
* 
* SPL-SLAM is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SPL-SLAM is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details. 
*
* 
*/

#ifndef MAPLINE_H
#define MAPLINE_H

#include"KeyFrame.h"
#include"Frame.h"
#include"Map.h"

#include<opencv2/core/core.hpp>
#include<mutex>

#include<map>
namespace PL_SLAM
{
  class KeyFrame;
  class Map;
  class Frame;
  
class MapLine
{
public:
    //前2个参数是3D线段端点，在构造函数体中建立3D的线特征中点。
    MapLine(const cv::Mat &FirPPos, const cv::Mat &EndPPos, 
	    const cv::Mat &MidPPos, KeyFrame* pRefKF, Map* pMap);
  
    void SetFirPWorldPos(const cv::Mat &FirPPos);
    void SetEndPWorldPos(const cv::Mat &EndPPos);
    void SetMidPWorldPos(const cv::Mat &MidPPos);
    cv::Mat GetFirPWorldPos();
    cv::Mat GetEndPWorldPos();
    cv::Mat GetMidPWorldPos();

    cv::Mat GetNormal();

    KeyFrame* GetReferenceKeyFrame();
    
    std::map<KeyFrame*,size_t> GetObservations();

    int Observations();
    void AddObservation(KeyFrame* pKF,size_t idx);
    void EraseObservation(KeyFrame* pKF);

    int GetIndexInKeyFrame(KeyFrame* pKF);
    bool IsInKeyFrame(KeyFrame* pKF);

    void SetBadFlag();
    bool isBad();

    void Replace(MapLine* pML);    

    MapLine* GetReplaced();

    void IncreaseVisible(int n=1);
    void IncreaseFound(int n=1);
    float GetFoundRatio();
    
    inline int GetFound(){
	return mnFound;
    }

    void ComputeDistinctiveDescriptors();
    cv::Mat GetDescriptor();
    
    float Get2DLineLengthAverage();
    void UpdateNormalAndDepth();
    void Update2DLineLength();

    //1.2倍最大深度 / 0.8倍最小深度：
    float GetMinDistanceInvariance();
    float GetMaxDistanceInvariance();
    
    //尺度预测：
    int PredictScale(const float &currentDist, KeyFrame*pKF);
    int PredictScale(const float &currentDist, Frame* pF);
public:   

    long unsigned int mnId;
    static long unsigned int nNextId;
    long int mnFirstKFid;

    long int mnFirstFrame;
    int nObs;

    //前端跟踪使用的变量
    float mTrackProjX;
    float mTrackProjY;
    bool mbTrackInView;
    int mnTrackScaleLevel;
    float mTrackViewCos;
    long unsigned int mnTrackReferenceForFrame;
    long unsigned int mnLastFrameSeen;

    //后端使用的变量：
    long unsigned int mnBALocalForKF;

    long unsigned int mnBALocalForKFBoth;

    long unsigned int mnFuseCandidateForKF;
    
    static std::mutex mGlobalMutex;
    
protected: 

    cv::Mat mWorldFirPPos;
    cv::Mat mWorldEndPPos;
    cv::Mat mWorldMidPPos;

    std::map<KeyFrame*,size_t> mObservations;
    cv::Mat mNormalVector;
    float m2DLineLengthAverage;

    //描述子：
    cv::Mat mDescriptor;
    KeyFrame* mpRefKF;

    int mnVisible;
    int mnFound;

    bool mbBad;

    MapLine* mpReplaced;

    float mfMinDistance;
    float mfMaxDistance;

    Map* mpMap;

    std::mutex mMutexPos;
    std::mutex mMutexFeatures;
  
}; 

}//namespace ORB_SLAM

#endif // MAPLINE_H
