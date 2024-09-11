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

#include "MapLine.h"
#include "Linematcher.h"

#include<mutex>

namespace PL_SLAM
{

long unsigned int MapLine::nNextId=0;
mutex MapLine::mGlobalMutex;

MapLine::MapLine(const cv::Mat &FirPPos, const cv::Mat &EndPPos, 
	    const cv::Mat &MidPPos, KeyFrame* pRefKF, Map* pMap):
    mnFirstKFid(pRefKF->mnId), mnFirstFrame(pRefKF->mnFrameId), nObs(0),mnTrackReferenceForFrame(0), mnLastFrameSeen(0), mnBALocalForKF(0),
    mnBALocalForKFBoth(0), mnFuseCandidateForKF(0),mpRefKF(pRefKF), mnVisible(1), mnFound(1), mbBad(false),mpReplaced(static_cast<MapLine*>(NULL)),
    mfMinDistance(0), mfMaxDistance(0), mpMap(pMap)
	    
{
  FirPPos.copyTo(mWorldFirPPos);
  EndPPos.copyTo(mWorldEndPPos);
  MidPPos.copyTo(mWorldMidPPos);

  mNormalVector = cv::Mat::zeros(3,1,CV_32F);

  unique_lock<mutex> lock(mpMap->mMutexLineCreation);
  mnId=nNextId++;
}

void MapLine::SetFirPWorldPos(const cv::Mat &FirPPos)
{
    unique_lock<mutex> lock2(mGlobalMutex);
    unique_lock<mutex> lock(mMutexPos);
    FirPPos.copyTo(mWorldFirPPos);
}

void MapLine::SetEndPWorldPos(const cv::Mat &EndPPos)
{
    unique_lock<mutex> lock2(mGlobalMutex);
    unique_lock<mutex> lock(mMutexPos);
    EndPPos.copyTo(mWorldEndPPos);
}

void MapLine::SetMidPWorldPos(const cv::Mat &MidPPos)
{
    unique_lock<mutex> lock2(mGlobalMutex);
    unique_lock<mutex> lock(mMutexPos);
    MidPPos.copyTo(mWorldMidPPos);
}
 
cv::Mat MapLine::GetFirPWorldPos()
{
    unique_lock<mutex> lock(mMutexPos);
    return mWorldFirPPos.clone();
}

cv::Mat MapLine::GetEndPWorldPos()
{
    unique_lock<mutex> lock(mMutexPos);
    return mWorldEndPPos.clone();
}

cv::Mat MapLine::GetMidPWorldPos()
{
    unique_lock<mutex> lock(mMutexPos);
    return mWorldMidPPos.clone();
}

cv::Mat MapLine::GetNormal()
{
    unique_lock<mutex> lock(mMutexPos);
    return mNormalVector.clone();
}

KeyFrame* MapLine::GetReferenceKeyFrame()
{
    unique_lock<mutex> lock(mMutexFeatures);
    return mpRefKF;
}

map<KeyFrame*,size_t> MapLine::GetObservations()
{
    unique_lock<mutex> lock(mMutexFeatures);
    return mObservations;
}

int MapLine::Observations()
{
    unique_lock<mutex> lock(mMutexFeatures);
    return nObs;
}

void MapLine::AddObservation(KeyFrame* pKF,size_t idx)
{
    unique_lock<mutex> lock(mMutexFeatures);
    if(mObservations.count(pKF))
        return;
    mObservations[pKF]=idx;
    
    if(pKF->mvuRight[idx]>=0)
        nObs+=2;
    else
        nObs++;
}
void MapLine::EraseObservation(KeyFrame* pKF)
{
    bool bBad=false;
    {
        unique_lock<mutex> lock(mMutexFeatures);
        if(mObservations.count(pKF))
        {
            int idx = mObservations[pKF];
            if(pKF->mvuRight[idx]>=0)
                nObs-=2;
            else
                nObs--;

            mObservations.erase(pKF);

            if(mpRefKF==pKF)
                mpRefKF=mObservations.begin()->first;

            // If only 2 observations or less, discard point
            if(nObs<=2)
                bBad=true;
        }
    }

    if(bBad)
        SetBadFlag();
}

int MapLine::GetIndexInKeyFrame(KeyFrame* pKF)
{
    unique_lock<mutex> lock(mMutexFeatures);
    if(mObservations.count(pKF))
        return mObservations[pKF];
    else
        return -1;
}

bool MapLine::IsInKeyFrame(KeyFrame* pKF)
{
    unique_lock<mutex> lock(mMutexFeatures);
    return (mObservations.count(pKF));
}

void MapLine::SetBadFlag()
{

    map<KeyFrame*,size_t> obs;
    {
        unique_lock<mutex> lock1(mMutexFeatures);
        unique_lock<mutex> lock2(mMutexPos);
        mbBad=true;
        obs = mObservations;
        mObservations.clear();
    }

    for(map<KeyFrame*,size_t>::iterator mit=obs.begin(), mend=obs.end(); mit!=mend; mit++)
    {
        KeyFrame* pKF = mit->first;
        pKF->EraseMapLineMatch(mit->second);
    }

    mpMap->EraseMapLine(this);
}

bool MapLine::isBad()
{
    unique_lock<mutex> lock(mMutexFeatures);
    unique_lock<mutex> lock2(mMutexPos);
    return mbBad;
}

void MapLine::Replace(MapLine* pML)
{
    if(pML->mnId==this->mnId)
        return;

    int nvisible, nfound;
    map<KeyFrame*,size_t> obs;
    {
        unique_lock<mutex> lock1(mMutexFeatures);
        unique_lock<mutex> lock2(mMutexPos);
        obs=mObservations;
        mObservations.clear();
        mbBad=true;
        nvisible = mnVisible;
        nfound = mnFound;
        mpReplaced = pML;
    }

    for(map<KeyFrame*,size_t>::iterator mit=obs.begin(), mend=obs.end(); mit!=mend; mit++)
    {
        // Replace measurement in keyframe
        KeyFrame* pKF = mit->first;

        if(!pML->IsInKeyFrame(pKF))
        {	
	    // 无冲突
            pKF->ReplaceMapLineMatch(mit->second, pML);
            pML->AddObservation(pKF,mit->second);
        }
        else
        {
            // 产生冲突
            pKF->EraseMapLineMatch(mit->second);
        }
    }
    pML->IncreaseFound(nfound);
    pML->IncreaseVisible(nvisible);
    pML->ComputeDistinctiveDescriptors();

    mpMap->EraseMapLine(this);
}

MapLine* MapLine::GetReplaced()
{
    unique_lock<mutex> lock1(mMutexFeatures);
    unique_lock<mutex> lock2(mMutexPos);
    return mpReplaced;
}

void MapLine::IncreaseVisible(int n)
{
    unique_lock<mutex> lock(mMutexFeatures);
    mnVisible+=n;
}

void MapLine::IncreaseFound(int n)
{
    unique_lock<mutex> lock(mMutexFeatures);
    mnFound+=n;
}

float MapLine::GetFoundRatio()
{
    unique_lock<mutex> lock(mMutexFeatures);
    return static_cast<float>(mnFound)/mnVisible;
}

//计算最佳点特征描述子
void MapLine::ComputeDistinctiveDescriptors()
{
  // Retrieve all observed descriptors
    vector<cv::Mat> vDescriptors;

    map<KeyFrame*,size_t> observations;

    {
        unique_lock<mutex> lock1(mMutexFeatures);
        if(mbBad)
            return;
        observations=mObservations;
    }

    if(observations.empty())
        return;

    vDescriptors.reserve(observations.size());

    for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
    {
        KeyFrame* pKF = mit->first;

        if(!pKF->isBadLines())
            vDescriptors.push_back(pKF->mDescriptorLines.row(mit->second));
    }

    if(vDescriptors.empty())
        return;

    // Compute distances between them
    const size_t N = vDescriptors.size();
	
    //float Distances[N][N];
	std::vector<std::vector<float> > Distances;
	Distances.resize(N, vector<float>(N, 0));
	for (size_t i = 0; i<N; i++)
    {
        Distances[i][i]=0;
        for(size_t j=i+1;j<N;j++)
        {
            int distij = Linematcher::DescriptorDistance(vDescriptors[i],vDescriptors[j]);
            Distances[i][j]=distij;
            Distances[j][i]=distij;
        }
    }

    // Take the descriptor with least median distance to the rest
    int BestMedian = INT_MAX;
    int BestIdx = 0;
    for(size_t i=0;i<N;i++)
    {
	vector<int> vDists(Distances[i].begin(), Distances[i].end());
	sort(vDists.begin(), vDists.end());

        int median = vDists[0.5*(N-1)];

        if(median<BestMedian)
        {
            BestMedian = median;
            BestIdx = i;
        }
    }

    {
        unique_lock<mutex> lock(mMutexFeatures);
        mDescriptor = vDescriptors[BestIdx].clone();       
    }
}

cv::Mat MapLine::GetDescriptor()
{
    unique_lock<mutex> lock(mMutexFeatures);
    return mDescriptor.clone();
}

float MapLine::Get2DLineLengthAverage()
{
    unique_lock<mutex> lock(mMutexFeatures);
    return m2DLineLengthAverage;
}

void MapLine::UpdateNormalAndDepth()
{
    map<KeyFrame*,size_t> observations;
    KeyFrame* pRefKF;
    cv::Mat Pos;
    {
        unique_lock<mutex> lock1(mMutexFeatures);
        unique_lock<mutex> lock2(mMutexPos);
        if(mbBad)
            return;
        observations=mObservations;
        pRefKF=mpRefKF;
        Pos = mWorldMidPPos.clone();
    }

    if(observations.empty())
        return;
    
    cv::Mat normal = cv::Mat::zeros(3,1,CV_32F);
    int n=0;
    for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
    {
        KeyFrame* pKF = mit->first;
        cv::Mat Owi = pKF->GetCameraCenter();
        cv::Mat normali = Pos - Owi;
        normal = normal + normali/cv::norm(normali);
        n++;
    }
    
    cv::Mat PC = Pos - pRefKF->GetCameraCenter();
    const float dist = cv::norm(PC);
    
    const int level = pRefKF->mvMidPointsUn[observations[pRefKF]].octave;
    const float levelScaleFactor =  pRefKF->mvScaleFactorsLines[level];
    const int nLevels = pRefKF->mnScaleLevelLines;

    {
        unique_lock<mutex> lock3(mMutexPos);
        mfMaxDistance = dist*levelScaleFactor;
        mfMinDistance = mfMaxDistance/pRefKF->mvScaleFactorsLines[nLevels-1];
        mNormalVector = normal/n;
    }
    
}

//1.2倍最大深度 / 0.8倍最小深度：
float MapLine::GetMinDistanceInvariance()
{
    unique_lock<mutex> lock(mMutexPos);
    return 0.8f*mfMinDistance;
}

float MapLine::GetMaxDistanceInvariance()
{
    unique_lock<mutex> lock(mMutexPos);
    return 1.2f*mfMaxDistance;
}

void MapLine::Update2DLineLength()
{
    map<KeyFrame*,size_t> observations;

    {
        unique_lock<mutex> lock1(mMutexFeatures);
        if(mbBad)
            return;
        observations=mObservations;
    }

    if(observations.empty())
        return;
    
    float m2DLineLengthSum = 0.0;
    int KFNumber = 0;
    
    for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
    {
        KeyFrame* pKF = mit->first;
    
        if(!pKF->isBadLines())
	{
	  m2DLineLengthSum = m2DLineLengthSum + pKF->mvLinesUn[mit->second].lineLength;	
	  KFNumber++;
	}
	
    }
    
    {
	unique_lock<mutex> lock(mMutexFeatures);
	m2DLineLengthAverage = m2DLineLengthSum/KFNumber;
    }  
}

int MapLine::PredictScale(const float &currentDist, KeyFrame*pKF)
{
    float ratio;
    {
        unique_lock<mutex> lock(mMutexPos);
        ratio = mfMaxDistance/currentDist;
    }

    int nScale = ceil(log(ratio)/pKF->mfLogScaleFactorsLines);
    if(nScale<0)
        nScale = 0;
    else if(nScale>=pKF->mnScaleLevelLines)
        nScale = pKF->mnScaleLevelLines-1;

    return nScale;
}

//用于前端
int MapLine::PredictScale(const float &currentDist, Frame* pF)
{
    float ratio;
    {
        unique_lock<mutex> lock(mMutexPos);
        ratio = mfMaxDistance/currentDist;
    }

    int nScale = ceil(log(ratio)/pF->mfLogScaleFactorsLines);
    if(nScale<0)
        nScale = 0;
    else if(nScale>=pF->mnScaleLevelLines)
        nScale = pF->mnScaleLevelLines-1;

    return nScale;
}


}//namespace ORB_SLAM