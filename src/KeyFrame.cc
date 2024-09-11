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

#include "KeyFrame.h"
#include "Converter.h"
#include "ORBmatcher.h"
#include<mutex>

namespace PL_SLAM
{

long unsigned int KeyFrame::nNextId=0;

KeyFrame::KeyFrame(Frame &F, Map *pMap, KeyFrameDatabase *pKFDB):
    mnFrameId(F.mnId),  mTimeStamp(F.mTimeStamp), mnGridCols(FRAME_GRID_COLS), mnGridRows(FRAME_GRID_ROWS),
    mfGridElementWidthInv(F.mfGridElementWidthInv), mfGridElementHeightInv(F.mfGridElementHeightInv),
    mnTrackReferenceForFrame(0), mnFuseTargetForKF(0), mnBALocalForKF(0), mnBAFixedForKF(0),
    mnLoopQuery(0), mnLoopWords(0), mnRelocQuery(0), mnRelocWords(0), mnBAGlobalForKF(0),
    fx(F.fx), fy(F.fy), cx(F.cx), cy(F.cy), invfx(F.invfx), invfy(F.invfy),
    mbf(F.mbf), mb(F.mb), mThDepth(F.mThDepth), N(F.N), mvKeys(F.mvKeys), mvKeysUn(F.mvKeysUn),
    mvuRight(F.mvuRight), mvDepth(F.mvDepth), mDescriptors(F.mDescriptors.clone()),
    mBowVec(F.mBowVec), mFeatVec(F.mFeatVec), mnScaleLevels(F.mnScaleLevels), mfScaleFactor(F.mfScaleFactor),
    mfLogScaleFactor(F.mfLogScaleFactor), mvScaleFactors(F.mvScaleFactors), mvLevelSigma2(F.mvLevelSigma2),
    mvInvLevelSigma2(F.mvInvLevelSigma2), mnMinX(F.mnMinX), mnMinY(F.mnMinY), mnMaxX(F.mnMaxX),
    mnMaxY(F.mnMaxY), mK(F.mK), mvpMapPoints(F.mvpMapPoints), mpKeyFrameDB(pKFDB),
    mpORBvocabulary(F.mpORBvocabulary), mbFirstConnection(true), mpParent(NULL), mbNotErase(false),
    mbToBeErased(false), mbBad(false), mHalfBaseline(F.mb/2), mpMap(pMap)
{
    mnId=nNextId++;

    mGrid.resize(mnGridCols);
    for(int i=0; i<mnGridCols;i++)
    {
        mGrid[i].resize(mnGridRows);
        for(int j=0; j<mnGridRows; j++)
            mGrid[i][j] = F.mGrid[i][j];
    }

    SetPose(F.mTcw);    
    
    if(F.mpLineextractor!=NULL)
    {
      //图像线特征分块记录
      mnGridColsLines = FRAME_GRID_COLS_Lines;
      mnGridRowsLines = FRAME_GRID_ROWS_Lines;
      mfGridElementWidthInvLines = F.mfGridElementWidthInvLines;
      mfGridElementHeightInvLines = F.mfGridElementHeightInvLines;
      
      mnTrackReferenceForFrameLines = 0;
      mnFuseTargetForKFLines = 0;
      mnBALocalForKFLines = 0;
      mnBAFixedForKFLines = 0;
      mBAPoseLinesandPointsForKF = 0;
      mnBALocalForKFBoth = 0;
      mnBAFixedForKFBoth = 0;
      
      NL = F.NL;
      mvLines = F.mvLines;
      mvMidPoints= F.mvMidPoints;
      mvLinesUn = F.mvLinesUn;
      mvMidPointsUn = F.mvMidPointsUn;
      mDescriptorLines = F.mDescriptorLines.clone();
      
      mnScaleLevelLines = F.mnScaleLevelLines;
      mfScaleFactorLines = F.mfScaleFactorLines;
      mfLogScaleFactorsLines = F.mfLogScaleFactorsLines;
      mvScaleFactorsLines = F.mvScaleFactorsLines;
      mvInvScaleFactorsLines = F.mvInvScaleFactorsLines;
      mvLevelSigma2Lines = F.mvLevelSigma2Lines;
      mvInvLevelSigma2Lines = F.mvInvLevelSigma2Lines;
	
      //地图线：
      mvpMapLines = F.mvpMapLines;
      mbFirstConnectionLines = true;
      mpParentLines = NULL;
      mbBadLines = false;
      
      //线特征网格分布拷贝
      mGridLines.resize(mnGridColsLines);
      for(int i=0; i<mnGridColsLines;i++)
      {
        mGridLines[i].resize(mnGridRowsLines);
        for(int j=0; j<mnGridRowsLines; j++)
            mGridLines[i][j] = F.mGridLines[i][j];
      }
    }
}

void KeyFrame::ComputeBoW()
{
    if(mBowVec.empty() || mFeatVec.empty())
    {
        vector<cv::Mat> vCurrentDesc = Converter::toDescriptorVector(mDescriptors);
        // Feature vector associate features with nodes in the 4th level (from leaves up)
        // We assume the vocabulary tree has 6 levels, change the 4 otherwise
        mpORBvocabulary->transform(vCurrentDesc,mBowVec,mFeatVec,4);
    }
}

void KeyFrame::SetPose(const cv::Mat &Tcw_)
{
    unique_lock<mutex> lock(mMutexPose);
    Tcw_.copyTo(Tcw);
    cv::Mat Rcw = Tcw.rowRange(0,3).colRange(0,3);
    cv::Mat tcw = Tcw.rowRange(0,3).col(3);
    cv::Mat Rwc = Rcw.t();
    Ow = -Rwc*tcw;

    Twc = cv::Mat::eye(4,4,Tcw.type());
    Rwc.copyTo(Twc.rowRange(0,3).colRange(0,3));
    Ow.copyTo(Twc.rowRange(0,3).col(3));
    cv::Mat center = (cv::Mat_<float>(4,1) << mHalfBaseline, 0 , 0, 1);
    Cw = Twc*center;
}

//设置关键帧的临时位姿T矩阵,局部BA优化时使用：(Point+LineExpanding)
void KeyFrame::SetPosePoints(const cv::Mat &TcwPoints_)
{
    unique_lock<mutex> lock(mMutexPose);
    TcwPoints_.copyTo(TcwPoints);
}
//设置关键帧的临时位姿T矩阵,局部BA优化时使用：(Point+LineExpanding)
void KeyFrame::SetPoseLines(const cv::Mat &TcwLines_)
{
    unique_lock<mutex> lock(mMutexPose);
    TcwLines_.copyTo(TcwLines);
}

cv::Mat KeyFrame::GetPose()
{
    unique_lock<mutex> lock(mMutexPose);
    return Tcw.clone();
}

cv::Mat KeyFrame::GetPoseInverse()
{
    unique_lock<mutex> lock(mMutexPose);
    return Twc.clone();
}

cv::Mat KeyFrame::GetCameraCenter()
{
    unique_lock<mutex> lock(mMutexPose);
    return Ow.clone();
}

cv::Mat KeyFrame::GetStereoCenter()
{
    unique_lock<mutex> lock(mMutexPose);
    return Cw.clone();
}


cv::Mat KeyFrame::GetRotation()
{
    unique_lock<mutex> lock(mMutexPose);
    return Tcw.rowRange(0,3).colRange(0,3).clone();
}

cv::Mat KeyFrame::GetTranslation()
{
    unique_lock<mutex> lock(mMutexPose);
    return Tcw.rowRange(0,3).col(3).clone();
}

void KeyFrame::AddConnection(KeyFrame *pKF, const int &weight)
{
    {
        unique_lock<mutex> lock(mMutexConnections);//注意这里互斥锁是线共视图互斥
        if(!mConnectedKeyFrameWeights.count(pKF))
            mConnectedKeyFrameWeights[pKF]=weight;
        else if(mConnectedKeyFrameWeights[pKF]!=weight)
            mConnectedKeyFrameWeights[pKF]=weight;
        else
            return;
    }

    UpdateBestCovisibles();
}

//添加线共视图及其权重：
void KeyFrame::AddConnectionLines(KeyFrame* pKF, const int &weight)
{
    {
	  unique_lock<mutex> lock(mMutexConnectionsLines);
	  if(!mConnectedKeyFrameWeightsLines.count(pKF))
	      mConnectedKeyFrameWeightsLines[pKF]=weight;
	  else if(mConnectedKeyFrameWeightsLines[pKF]!=weight)
	      mConnectedKeyFrameWeightsLines[pKF]=weight;
	  else
	      return;
      }

      UpdateBestCovisiblesLines();
}

void KeyFrame::UpdateBestCovisibles()
{
    unique_lock<mutex> lock(mMutexConnections);
    vector<pair<int,KeyFrame*> > vPairs;
    vPairs.reserve(mConnectedKeyFrameWeights.size());
    for(map<KeyFrame*,int>::iterator mit=mConnectedKeyFrameWeights.begin(), mend=mConnectedKeyFrameWeights.end(); mit!=mend; mit++)
       vPairs.push_back(make_pair(mit->second,mit->first));

    sort(vPairs.begin(),vPairs.end());
    list<KeyFrame*> lKFs;
    list<int> lWs;
    for(size_t i=0, iend=vPairs.size(); i<iend;i++)
    {
        lKFs.push_front(vPairs[i].second);
        lWs.push_front(vPairs[i].first);
    }

    mvpOrderedConnectedKeyFrames = vector<KeyFrame*>(lKFs.begin(),lKFs.end());
    mvOrderedWeights = vector<int>(lWs.begin(), lWs.end());    
}

//更新线最佳共视图
void KeyFrame::UpdateBestCovisiblesLines()
{
    unique_lock<mutex> lock(mMutexConnectionsLines);
    vector<pair<int,KeyFrame*> > vPairs;
    vPairs.reserve(mConnectedKeyFrameWeightsLines.size());
    for(map<KeyFrame*,int>::iterator mit=mConnectedKeyFrameWeightsLines.begin(), mend=mConnectedKeyFrameWeightsLines.end(); mit!=mend; mit++)
       vPairs.push_back(make_pair(mit->second,mit->first));

    sort(vPairs.begin(),vPairs.end());
    list<KeyFrame*> lKFs;
    list<int> lWs;
    for(size_t i=0, iend=vPairs.size(); i<iend;i++)
    {
        lKFs.push_front(vPairs[i].second);
        lWs.push_front(vPairs[i].first);
    }

    mvpOrderedConnectedKeyFramesLines = vector<KeyFrame*>(lKFs.begin(),lKFs.end());
    mvOrderedWeightsLines = vector<int>(lWs.begin(), lWs.end());  
}

set<KeyFrame*> KeyFrame::GetConnectedKeyFrames()
{
    unique_lock<mutex> lock(mMutexConnections);
    set<KeyFrame*> s;
    for(map<KeyFrame*,int>::iterator mit=mConnectedKeyFrameWeights.begin();mit!=mConnectedKeyFrameWeights.end();mit++)
        s.insert(mit->first);
    return s;
}

//获取原始无序的线共视图中的关键帧集
set<KeyFrame*> KeyFrame::GetConnectedKeyFramesLines()
{
    unique_lock<mutex> lock(mMutexConnectionsLines);
    set<KeyFrame*> s;
    for(map<KeyFrame*,int>::iterator mit=mConnectedKeyFrameWeightsLines.begin();mit!=mConnectedKeyFrameWeightsLines.end();mit++)
        s.insert(mit->first);
    return s;
}

vector<KeyFrame*> KeyFrame::GetVectorCovisibleKeyFrames()
{
    unique_lock<mutex> lock(mMutexConnections);
    return mvpOrderedConnectedKeyFrames;
}

//获取线共视图的有序关键帧集
vector<KeyFrame*> KeyFrame::GetVectorCovisibleKeyFramesLines()
{
    unique_lock<mutex> lock(mMutexConnectionsLines);
    return mvpOrderedConnectedKeyFramesLines;
}

vector<KeyFrame*> KeyFrame::GetBestCovisibilityKeyFrames(const int &N)
{
    unique_lock<mutex> lock(mMutexConnections);
    if((int)mvpOrderedConnectedKeyFrames.size()<N)
        return mvpOrderedConnectedKeyFrames;
    else
        return vector<KeyFrame*>(mvpOrderedConnectedKeyFrames.begin(),mvpOrderedConnectedKeyFrames.begin()+N);

}

//获取线共视图中的最佳前N个关键帧
vector<KeyFrame*> KeyFrame::GetBestCovisibilityKeyFramesLines(const int &N)
{
    unique_lock<mutex> lock(mMutexConnectionsLines);
    if((int)mvpOrderedConnectedKeyFramesLines.size()<N)
        return mvpOrderedConnectedKeyFramesLines;
    else
        return vector<KeyFrame*>(mvpOrderedConnectedKeyFramesLines.begin(),mvpOrderedConnectedKeyFramesLines.begin()+N);

}

vector<KeyFrame*> KeyFrame::GetCovisiblesByWeight(const int &w)
{
    unique_lock<mutex> lock(mMutexConnections);

    if(mvpOrderedConnectedKeyFrames.empty())
        return vector<KeyFrame*>();

    vector<int>::iterator it = upper_bound(mvOrderedWeights.begin(),mvOrderedWeights.end(),w,KeyFrame::weightComp);
    if(it==mvOrderedWeights.end())
        return vector<KeyFrame*>();
    else
    {
        int n = it-mvOrderedWeights.begin();
        return vector<KeyFrame*>(mvpOrderedConnectedKeyFrames.begin(), mvpOrderedConnectedKeyFrames.begin()+n);
    }
}

//通过权重获取有序线共视图关键帧集中的关键帧（获取大于某权重的所有关键帧）
vector<KeyFrame*> KeyFrame::GetCovisiblesByWeightLines(const int &w)
{
    unique_lock<mutex> lock(mMutexConnectionsLines);

    if(mvpOrderedConnectedKeyFramesLines.empty())
        return vector<KeyFrame*>();

    vector<int>::iterator it = upper_bound(mvOrderedWeightsLines.begin(),mvOrderedWeightsLines.end(),w,KeyFrame::weightComp);
    if(it==mvOrderedWeightsLines.end())
        return vector<KeyFrame*>();
    else
    {
        int n = it-mvOrderedWeightsLines.begin();
        return vector<KeyFrame*>(mvpOrderedConnectedKeyFramesLines.begin(), mvpOrderedConnectedKeyFramesLines.begin()+n);
    }
}

int KeyFrame::GetWeight(KeyFrame *pKF)
{
    unique_lock<mutex> lock(mMutexConnections);
    if(mConnectedKeyFrameWeights.count(pKF))
        return mConnectedKeyFrameWeights[pKF];
    else
        return 0;
}

//通过线共视图关键帧，获取对应的权重值
int KeyFrame::GetWeightLines(KeyFrame *pKF)
{
    unique_lock<mutex> lock(mMutexConnectionsLines);
    if(mConnectedKeyFrameWeightsLines.count(pKF))
        return mConnectedKeyFrameWeightsLines[pKF];
    else
        return 0;
}

void KeyFrame::AddMapPoint(MapPoint *pMP, const size_t &idx)
{
    unique_lock<mutex> lock(mMutexFeatures);
    mvpMapPoints[idx]=pMP;
}

//地图下线向量更新
//添加地图线
void KeyFrame::AddMapLine(MapLine *pML, const size_t &idx)
{
    unique_lock<mutex> lock(mMutexFeaturesLines);
    mvpMapLines[idx]=pML;
}

void KeyFrame::EraseMapPointMatch(const size_t &idx)
{
    unique_lock<mutex> lock(mMutexFeatures);
    mvpMapPoints[idx]=static_cast<MapPoint*>(NULL);
}

//删除地图线：通过索引删除
void KeyFrame::EraseMapLineMatch(const size_t &idx)
{
    unique_lock<mutex> lock(mMutexFeaturesLines);
    mvpMapLines[idx]=static_cast<MapLine*>(NULL);
}

void KeyFrame::EraseMapPointMatch(MapPoint* pMP)
{
    int idx = pMP->GetIndexInKeyFrame(this);
    if(idx>=0)
        mvpMapPoints[idx]=static_cast<MapPoint*>(NULL);
}

//删除地图线：直接给出地图线
void KeyFrame::EraseMapLineMatch(MapLine* pML)
{
    int idx = pML->GetIndexInKeyFrame(this);
    if(idx>=0)
        mvpMapLines[idx]=static_cast<MapLine*>(NULL);
}

void KeyFrame::ReplaceMapPointMatch(const size_t &idx, MapPoint* pMP)
{
    mvpMapPoints[idx]=pMP;
}

//替换指定位置的地图线向量
void KeyFrame::ReplaceMapLineMatch(const size_t &idx, MapLine* pML)
{
    mvpMapLines[idx]=pML;
}

set<MapPoint*> KeyFrame::GetMapPoints()
{
    unique_lock<mutex> lock(mMutexFeatures);
    set<MapPoint*> s;
    for(size_t i=0, iend=mvpMapPoints.size(); i<iend; i++)
    {
        if(!mvpMapPoints[i])
            continue;
        MapPoint* pMP = mvpMapPoints[i];
        if(!pMP->isBad())
            s.insert(pMP);
    }
    return s;
}

//获取基于红黑树存储的(非坏)地图线
set<MapLine*> KeyFrame::GetMapLines()
{
    unique_lock<mutex> lock(mMutexFeaturesLines);
    set<MapLine*> s;
    for(size_t i=0, iend=mvpMapLines.size(); i<iend; i++)
    {
        if(!mvpMapLines[i])
            continue;
        MapLine* pML = mvpMapLines[i];
        if(!pML->isBad())
            s.insert(pML);
    }
    return s;
}

int KeyFrame::TrackedMapPoints(const int &minObs)
{
    unique_lock<mutex> lock(mMutexFeatures);

    int nPoints=0;
    const bool bCheckObs = minObs>0;
    for(int i=0; i<N; i++)
    {
        MapPoint* pMP = mvpMapPoints[i];
        if(pMP)
        {
            if(!pMP->isBad())
            {
                if(bCheckObs)
                {
                    if(mvpMapPoints[i]->Observations()>=minObs)
                        nPoints++;
                }
                else
                    nPoints++;
            }
        }
    }

    return nPoints;
}

//获取大于指定观测数的当前关键帧中的所有地图线数量
int KeyFrame::TrackedMapLines(const int &minObs)
{
    unique_lock<mutex> lock(mMutexFeaturesLines);

    int nLines=0;
    const bool bCheckObs = minObs>0;
    for(int i=0; i<NL; i++)
    {
        MapLine* pML = mvpMapLines[i];
        if(pML)
        {
            if(!pML->isBad())
            {
                if(bCheckObs)
                {
                    if(mvpMapLines[i]->Observations()>=minObs)
                        nLines++;
                }
                else
                    nLines++;
            }
        }
    }

    return nLines;
}

vector<MapPoint*> KeyFrame::GetMapPointMatches()
{
    unique_lock<mutex> lock(mMutexFeatures);
    return mvpMapPoints;
}

//获取关键帧中存储的地图线向量
vector<MapLine*> KeyFrame::GetMapLineMatches()
{
    unique_lock<mutex> lock(mMutexFeaturesLines);
    return mvpMapLines;
}

MapPoint* KeyFrame::GetMapPoint(const size_t &idx)
{
    unique_lock<mutex> lock(mMutexFeatures);
    return mvpMapPoints[idx];
}

//获取地图线向量指定索引出的地图线
MapLine* KeyFrame::GetMapLine(const size_t &idx)
{
    unique_lock<mutex> lock(mMutexFeaturesLines);
    return mvpMapLines[idx];
}

void KeyFrame::UpdateConnections()
{
    map<KeyFrame*,int> KFcounter;

    vector<MapPoint*> vpMP;

    {
        unique_lock<mutex> lockMPs(mMutexFeatures);
        vpMP = mvpMapPoints;
    }

    //For all map points in keyframe check in which other keyframes are they seen
    //Increase counter for those keyframes
    for(vector<MapPoint*>::iterator vit=vpMP.begin(), vend=vpMP.end(); vit!=vend; vit++)
    {
        MapPoint* pMP = *vit;

        if(!pMP)
            continue;

        if(pMP->isBad())
            continue;

        map<KeyFrame*,size_t> observations = pMP->GetObservations();

        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            if(mit->first->mnId==mnId)
                continue;
            KFcounter[mit->first]++;
        }
    }

    // This should not happen
    if(KFcounter.empty())
        return;

    //If the counter is greater than threshold add connection
    //In case no keyframe counter is over threshold add the one with maximum counter
    int nmax=0;
    KeyFrame* pKFmax=NULL;
    int th = 15;

    vector<pair<int,KeyFrame*> > vPairs;
    vPairs.reserve(KFcounter.size());
    for(map<KeyFrame*,int>::iterator mit=KFcounter.begin(), mend=KFcounter.end(); mit!=mend; mit++)
    {
        if(mit->second>nmax)
        {
            nmax=mit->second;
            pKFmax=mit->first;
        }
        if(mit->second>=th)
        {
            vPairs.push_back(make_pair(mit->second,mit->first));
            (mit->first)->AddConnection(this,mit->second);
        }
    }

    if(vPairs.empty())
    {
        vPairs.push_back(make_pair(nmax,pKFmax));
        pKFmax->AddConnection(this,nmax);
    }

    sort(vPairs.begin(),vPairs.end());
    list<KeyFrame*> lKFs;
    list<int> lWs;
    for(size_t i=0; i<vPairs.size();i++)
    {
        lKFs.push_front(vPairs[i].second);
        lWs.push_front(vPairs[i].first);
    }

    {
        unique_lock<mutex> lockCon(mMutexConnections);

        // mspConnectedKeyFrames = spConnectedKeyFrames;
        mConnectedKeyFrameWeights = KFcounter;
        mvpOrderedConnectedKeyFrames = vector<KeyFrame*>(lKFs.begin(),lKFs.end());
        mvOrderedWeights = vector<int>(lWs.begin(), lWs.end());

        if(mbFirstConnection && mnId!=0)
        {
            mpParent = mvpOrderedConnectedKeyFrames.front();
            mpParent->AddChild(this);
            mbFirstConnection = false;
        }

    }
}

//更新线特征共视图连接
void KeyFrame::UpdateConnectionsLines()
{
    map<KeyFrame*,int> KFcounter;

    vector<MapLine*> vpML;

    {
        unique_lock<mutex> lockMPs(mMutexFeaturesLines);
        vpML = mvpMapLines;
    }

    //For all map Lines in keyframe check in which other keyframes are they seen
    //Increase counter for those keyframes
    for(vector<MapLine*>::iterator vit=vpML.begin(), vend=vpML.end(); vit!=vend; vit++)
    {
        MapLine* pML = *vit;

        if(!pML)
            continue;

        if(pML->isBad())
            continue;

        map<KeyFrame*,size_t> observations = pML->GetObservations();

        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            if(mit->first->mnId==mnId)
                continue;
            KFcounter[mit->first]++;
        }
    }

    // This should not happen
    if(KFcounter.empty())
        return;

    //If the counter is greater than threshold add connection
    //In case no keyframe counter is over threshold add the one with maximum counter
    int nmax=0;
    KeyFrame* pKFmax=NULL;
    int th = 15; 

    vector<pair<int,KeyFrame*> > vPairs;
    vPairs.reserve(KFcounter.size());
    for(map<KeyFrame*,int>::iterator mit=KFcounter.begin(), mend=KFcounter.end(); mit!=mend; mit++)
    {
        if(mit->second>nmax)
        {
            nmax=mit->second;
            pKFmax=mit->first;
        }
        if(mit->second>=th)
        {
            vPairs.push_back(make_pair(mit->second,mit->first));
            (mit->first)->AddConnectionLines(this,mit->second);
        }
    }

    if(vPairs.empty())
    {
        vPairs.push_back(make_pair(nmax,pKFmax));
        pKFmax->AddConnectionLines(this,nmax);
    }

    sort(vPairs.begin(),vPairs.end());
    list<KeyFrame*> lKFs;
    list<int> lWs;
    for(size_t i=0; i<vPairs.size();i++)
    {
        lKFs.push_front(vPairs[i].second);
        lWs.push_front(vPairs[i].first);
    }

    {
        unique_lock<mutex> lockCon(mMutexConnectionsLines);

        // mspConnectedKeyFrames = spConnectedKeyFrames;
        mConnectedKeyFrameWeightsLines = KFcounter;
        mvpOrderedConnectedKeyFramesLines = vector<KeyFrame*>(lKFs.begin(),lKFs.end());
        mvOrderedWeightsLines = vector<int>(lWs.begin(), lWs.end());

        if(mbFirstConnectionLines && mnId!=0)
        {
            mpParentLines = mvpOrderedConnectedKeyFramesLines.front();
            mpParentLines->AddChildLines(this);
            mbFirstConnectionLines = false;
        }

    }
}

void KeyFrame::AddChild(KeyFrame *pKF)
{
    unique_lock<mutex> lockCon(mMutexConnections);
    mspChildrens.insert(pKF);
}

//更新前端局部地图跟踪使用的父子生成树
void KeyFrame::AddChildLines(KeyFrame *pKF)
{
    unique_lock<mutex> lockCon(mMutexConnectionsLines);
    mspChildrensLines.insert(pKF);
}

void KeyFrame::EraseChild(KeyFrame *pKF)
{
    unique_lock<mutex> lockCon(mMutexConnections);
    mspChildrens.erase(pKF);
}

//添加线孩子关键帧
void KeyFrame::EraseChildLines(KeyFrame *pKF)
{
    unique_lock<mutex> lockCon(mMutexConnectionsLines);
    mspChildrensLines.erase(pKF);
}

void KeyFrame::ChangeParent(KeyFrame *pKF)
{
    unique_lock<mutex> lockCon(mMutexConnections);
    mpParent = pKF;
    pKF->AddChild(this);
}

//改变当前关键帧的线父亲关键帧
void KeyFrame::ChangeParentLine(KeyFrame *pKF)
{
    unique_lock<mutex> lockCon(mMutexConnectionsLines);
    mpParentLines = pKF;
    pKF->AddChildLines(this);
}

set<KeyFrame*> KeyFrame::GetChilds()
{
    unique_lock<mutex> lockCon(mMutexConnections);
    return mspChildrens;
}

//获取所有的线孩子关键帧
set<KeyFrame*> KeyFrame::GetChildsLines()
{
    unique_lock<mutex> lockCon(mMutexConnectionsLines);
    return mspChildrensLines;
}

KeyFrame* KeyFrame::GetParent()
{
    unique_lock<mutex> lockCon(mMutexConnections);
    return mpParent;
}

//获取线父亲关键帧
KeyFrame* KeyFrame::GetParentLines()
{
    unique_lock<mutex> lockCon(mMutexConnectionsLines);
    return mpParentLines;
}

bool KeyFrame::hasChild(KeyFrame *pKF)
{
    unique_lock<mutex> lockCon(mMutexConnections);
    return mspChildrens.count(pKF);
}

//判断线孩子关键帧集是否输入关键帧
bool KeyFrame::hasChildLines(KeyFrame *pKF)
{
    unique_lock<mutex> lockCon(mMutexConnectionsLines);
    return mspChildrensLines.count(pKF);
}

void KeyFrame::AddLoopEdge(KeyFrame *pKF)
{
    unique_lock<mutex> lockCon(mMutexConnections);
    mbNotErase = true;
    mspLoopEdges.insert(pKF);
}

set<KeyFrame*> KeyFrame::GetLoopEdges()
{
    unique_lock<mutex> lockCon(mMutexConnections);
    return mspLoopEdges;
}

void KeyFrame::SetNotErase()
{
    unique_lock<mutex> lock(mMutexConnections);
    mbNotErase = true;
}

void KeyFrame::SetErase()
{
    {
        unique_lock<mutex> lock(mMutexConnections);
        if(mspLoopEdges.empty())
        {
            mbNotErase = false;
        }
    }

    if(mbToBeErased)
    {
        SetBadFlag();
    }
}

void KeyFrame::SetBadFlag()
{   
    {
        unique_lock<mutex> lock(mMutexConnections);
        if(mnId==0)
            return;
        else if(mbNotErase)
        {
            mbToBeErased = true;
            return;
        }
    }
    
    for(map<KeyFrame*,int>::iterator mit = mConnectedKeyFrameWeights.begin(), mend=mConnectedKeyFrameWeights.end(); mit!=mend; mit++)
        mit->first->EraseConnection(this);

    for(size_t i=0; i<mvpMapPoints.size(); i++)
        if(mvpMapPoints[i])
            mvpMapPoints[i]->EraseObservation(this);
	
    {
        unique_lock<mutex> lock(mMutexConnections);
        unique_lock<mutex> lock1(mMutexFeatures);
	
	//清除有序的
        mConnectedKeyFrameWeights.clear();
        mvpOrderedConnectedKeyFrames.clear();

        // Update Spanning Tree
        set<KeyFrame*> sParentCandidates;
        sParentCandidates.insert(mpParent);

        // Assign at each iteration one children with a parent (the pair with highest covisibility weight)
        // Include that children as new parent candidate for the rest
        while(!mspChildrens.empty())
        {
            bool bContinue = false;

            int max = -1;
            KeyFrame* pC;
            KeyFrame* pP;

            for(set<KeyFrame*>::iterator sit=mspChildrens.begin(), send=mspChildrens.end(); sit!=send; sit++)
            {
                KeyFrame* pKF = *sit;
                if(pKF->isBad())
                    continue;

                // Check if a parent candidate is connected to the keyframe
                vector<KeyFrame*> vpConnected = pKF->GetVectorCovisibleKeyFrames();
                for(size_t i=0, iend=vpConnected.size(); i<iend; i++)
                {
                    for(set<KeyFrame*>::iterator spcit=sParentCandidates.begin(), spcend=sParentCandidates.end(); spcit!=spcend; spcit++)
                    {
                        if(vpConnected[i]->mnId == (*spcit)->mnId)
                        {
                            int w = pKF->GetWeight(vpConnected[i]);
                            if(w>max)
                            {
                                pC = pKF;
                                pP = vpConnected[i];
                                max = w;
                                bContinue = true;
                            }
                        }
                    }
                }
            }

            if(bContinue)
            {
                pC->ChangeParent(pP);
                sParentCandidates.insert(pC);
                mspChildrens.erase(pC);
            }
            else
                break;
        }

        // If a children has no covisibility links with any parent candidate, assign to the original parent of this KF
        if(!mspChildrens.empty())
            for(set<KeyFrame*>::iterator sit=mspChildrens.begin(); sit!=mspChildrens.end(); sit++)
            {
                (*sit)->ChangeParent(mpParent);
            }

        mpParent->EraseChild(this);
        mTcp = Tcw*mpParent->GetPoseInverse();
        mbBad = true;
    }


    mpMap->EraseKeyFrame(this);
    mpKeyFrameDB->erase(this);
}

//点线特征使用的关键帧设置：点关键帧删除
void KeyFrame::SetBadFlagPoints()
{   
    //判断是否为不允许删除的关键帧
    {
        unique_lock<mutex> lock(mMutexConnections);
        if(mnId==0)
            return;
        else if(mbNotErase)
        {
            mbToBeErased = true;//假意删除，回环检测结束后才真正删除
            return;
        }
    }
    
    //删除无序的共视图集中所有关键帧与当前关键帧的连接
    for(map<KeyFrame*,int>::iterator mit = mConnectedKeyFrameWeights.begin(), mend=mConnectedKeyFrameWeights.end(); mit!=mend; mit++)
        mit->first->EraseConnection(this);
    
    //对所有当前关键帧中的地图点，删除对当前关键帧的观测关系
    for(size_t i=0; i<mvpMapPoints.size(); i++)
        if(mvpMapPoints[i])
            mvpMapPoints[i]->EraseObservation(this);
	
    {
        unique_lock<mutex> lock(mMutexConnections);
        unique_lock<mutex> lock1(mMutexFeatures);//这个锁有意义么？

	//清除有序的共视图关键帧集、权重集
        mConnectedKeyFrameWeights.clear();
        mvpOrderedConnectedKeyFrames.clear();

        // Update Spanning Tree
        set<KeyFrame*> sParentCandidates;
        sParentCandidates.insert(mpParent);

        // Assign at each iteration one children with a parent (the pair with highest covisibility weight)
        // Include that children as new parent candidate for the rest
        while(!mspChildrens.empty())
        {
            bool bContinue = false;

            int max = -1;
            KeyFrame* pC;
            KeyFrame* pP;
	    
	    //当步大循环:为某一个孩子找到新父亲
            for(set<KeyFrame*>::iterator sit=mspChildrens.begin(), send=mspChildrens.end(); sit!=send; sit++)
            {
                KeyFrame* pKF = *sit;
                if(pKF->isBad())
                    continue;

                // Check if a parent candidate is connected to the keyframe
                vector<KeyFrame*> vpConnected = pKF->GetVectorCovisibleKeyFrames();
                for(size_t i=0, iend=vpConnected.size(); i<iend; i++)
                {
                    for(set<KeyFrame*>::iterator spcit=sParentCandidates.begin(), spcend=sParentCandidates.end(); spcit!=spcend; spcit++)
                    {
                        if(vpConnected[i]->mnId == (*spcit)->mnId)
                        {
                            int w = pKF->GetWeight(vpConnected[i]);
                            if(w>max)
                            {
                                pC = pKF;
                                pP = vpConnected[i];
                                max = w;
                                bContinue = true;
                            }
                        }
                    }
                }
                
                //当前孩子关键帧计算结束，判断该孩子是否找到新父亲，若找到则提前退出当步大循环
                if(bContinue)
		    break;
                
            }
	    
	    //若找到，则更新父子关系
            if(bContinue)
            {
                pC->ChangeParent(pP);
                sParentCandidates.insert(pC);
                mspChildrens.erase(pC);
            }
            else
                break;
        }

        // If a children has no covisibility links with any parent candidate, assign to the original parent of this KF
        if(!mspChildrens.empty())
            for(set<KeyFrame*>::iterator sit=mspChildrens.begin(); sit!=mspChildrens.end(); sit++)
            {
                (*sit)->ChangeParent(mpParent);
            }

        mpParent->EraseChild(this);
        mTcp = Tcw*mpParent->GetPoseInverse();
        mbBad = true;
    }
    
    mpKeyFrameDB->erase(this);
}


void KeyFrame::SetBadFlagLines()
{   

    for(map<KeyFrame*,int>::iterator mit = mConnectedKeyFrameWeightsLines.begin(), mend=mConnectedKeyFrameWeightsLines.end(); mit!=mend; mit++)
        mit->first->EraseConnectionLines(this);

    for(size_t i=0; i<mvpMapLines.size(); i++)
        if(mvpMapLines[i])
            mvpMapLines[i]->EraseObservation(this);
	
    {
        unique_lock<mutex> lock(mMutexConnectionsLines);
        unique_lock<mutex> lock1(mMutexFeaturesLines);

        mConnectedKeyFrameWeightsLines.clear();
        mvpOrderedConnectedKeyFramesLines.clear();

        // Update Spanning Tree
        set<KeyFrame*> sParentCandidates;
        sParentCandidates.insert(mpParentLines);

        // Assign at each iteration one children with a parent (the pair with highest covisibility weight)
        // Include that children as new parent candidate for the rest
        while(!mspChildrensLines.empty())
        {
            bool bContinue = false;

            int max = -1;
            KeyFrame* pC;
            KeyFrame* pP;

            for(set<KeyFrame*>::iterator sit=mspChildrensLines.begin(), send=mspChildrensLines.end(); sit!=send; sit++)
            {
                KeyFrame* pKF = *sit;
                if(pKF->isBadLines())
                    continue;

                // Check if a parent candidate is connected to the keyframe
                vector<KeyFrame*> vpConnected = pKF->GetVectorCovisibleKeyFramesLines();
                for(size_t i=0, iend=vpConnected.size(); i<iend; i++)
                {
                    for(set<KeyFrame*>::iterator spcit=sParentCandidates.begin(), spcend=sParentCandidates.end(); spcit!=spcend; spcit++)
                    {
                        if(vpConnected[i]->mnId == (*spcit)->mnId)
                        {
                            int w = pKF->GetWeightLines(vpConnected[i]);
                            if(w>max)
                            {
                                pC = pKF;
                                pP = vpConnected[i];
                                max = w;
                                bContinue = true;
                            }
                        }
                    }
                }
                
                //当前孩子关键帧计算结束，判断该孩子是否找到新父亲，若找到则提前退出当步大循环
                if(bContinue)
		    break;
                
            }
	    
	    //若找到，则更新父子关系
            if(bContinue)
            {
                pC->ChangeParentLine(pP);
                sParentCandidates.insert(pC);
                mspChildrensLines.erase(pC);
            }
            else
                break;
        }

        // If a children has no covisibility links with any parent candidate, assign to the original parent of this KF
        if(!mspChildrensLines.empty())
            for(set<KeyFrame*>::iterator sit=mspChildrensLines.begin(); sit!=mspChildrensLines.end(); sit++)
            {
                (*sit)->ChangeParentLine(mpParent);
            }

        mpParentLines->EraseChildLines(this);
        mTcpLines = Tcw*mpParentLines->GetPoseInverse();
        mbBadLines = true;
    }

}

bool KeyFrame::isBad()
{
    unique_lock<mutex> lock(mMutexConnections);
    return mbBad;
}

//判断当前关键帧是否为坏线关键帧
bool KeyFrame::isBadLines()
{
    unique_lock<mutex> lock(mMutexConnectionsLines);
    return mbBadLines;
}

void KeyFrame::EraseConnection(KeyFrame* pKF)
{
    bool bUpdate = false;
    {
        unique_lock<mutex> lock(mMutexConnections);
        if(mConnectedKeyFrameWeights.count(pKF))
        {
            mConnectedKeyFrameWeights.erase(pKF);
            bUpdate=true;
        }
    }

    if(bUpdate)
        UpdateBestCovisibles();
}

void KeyFrame::EraseConnectionLines(KeyFrame* pKF)
{
    bool bUpdate = false;
    {
        unique_lock<mutex> lock(mMutexConnectionsLines);
        if(mConnectedKeyFrameWeightsLines.count(pKF))
        {
            mConnectedKeyFrameWeightsLines.erase(pKF);
            bUpdate=true;
        }
    }

    if(bUpdate)
        UpdateBestCovisiblesLines();
}

vector<size_t> KeyFrame::GetFeaturesInArea(const float &x, const float &y, const float &r) const
{
    vector<size_t> vIndices;
    vIndices.reserve(N);

    const int nMinCellX = max(0,(int)floor((x-mnMinX-r)*mfGridElementWidthInv));
    if(nMinCellX>=mnGridCols)
        return vIndices;

    const int nMaxCellX = min((int)mnGridCols-1,(int)ceil((x-mnMinX+r)*mfGridElementWidthInv));
    if(nMaxCellX<0)
        return vIndices;

    const int nMinCellY = max(0,(int)floor((y-mnMinY-r)*mfGridElementHeightInv));
    if(nMinCellY>=mnGridRows)
        return vIndices;

    const int nMaxCellY = min((int)mnGridRows-1,(int)ceil((y-mnMinY+r)*mfGridElementHeightInv));
    if(nMaxCellY<0)
        return vIndices;

    for(int ix = nMinCellX; ix<=nMaxCellX; ix++)
    {
        for(int iy = nMinCellY; iy<=nMaxCellY; iy++)
        {
            const vector<size_t> vCell = mGrid[ix][iy];
            for(size_t j=0, jend=vCell.size(); j<jend; j++)
            {
                const cv::KeyPoint &kpUn = mvKeysUn[vCell[j]];
                const float distx = kpUn.pt.x-x;
                const float disty = kpUn.pt.y-y;

                if(fabs(distx)<r && fabs(disty)<r)
                    vIndices.push_back(vCell[j]);
            }
        }
    }

    return vIndices;
}

//获取关键帧指定大小半径内所有线特征索引
vector<size_t> KeyFrame::GetFeaturesInAreaLines(const float &x, const float &y, const float &r) const
{
    vector<size_t> vIndices;
    vIndices.reserve(NL);

    const int nMinCellX = max(0,(int)floor((x-mnMinX-r)*mfGridElementWidthInvLines));
    if(nMinCellX>=mnGridColsLines)
        return vIndices;

    const int nMaxCellX = min((int)mnGridColsLines-1,(int)ceil((x-mnMinX+r)*mfGridElementWidthInvLines));
    if(nMaxCellX<0)
        return vIndices;

    const int nMinCellY = max(0,(int)floor((y-mnMinY-r)*mfGridElementHeightInvLines));
    if(nMinCellY>=mnGridRowsLines)
        return vIndices;

    const int nMaxCellY = min((int)mnGridRowsLines-1,(int)ceil((y-mnMinY+r)*mfGridElementHeightInvLines));
    if(nMaxCellY<0)
        return vIndices;

    for(int ix = nMinCellX; ix<=nMaxCellX; ix++)
    {
        for(int iy = nMinCellY; iy<=nMaxCellY; iy++)
        {
            const vector<size_t> vCell = mGridLines[ix][iy];
            for(size_t j=0, jend=vCell.size(); j<jend; j++)
            {
                const cv::KeyPoint &kpUn = mvMidPointsUn[vCell[j]];
                const float distx = kpUn.pt.x-x;
                const float disty = kpUn.pt.y-y;

                if(fabs(distx)<r && fabs(disty)<r)
                    vIndices.push_back(vCell[j]);
            }
        }
    }

    return vIndices;
}

bool KeyFrame::IsInImage(const float &x, const float &y) const
{
    return (x>=mnMinX && x<mnMaxX && y>=mnMinY && y<mnMaxY);
}

cv::Mat KeyFrame::UnprojectStereo(int i)
{
    const float z = mvDepth[i];
    if(z>0)
    {
        const float u = mvKeys[i].pt.x;
        const float v = mvKeys[i].pt.y;
        const float x = (u-cx)*z*invfx;
        const float y = (v-cy)*z*invfy;
        cv::Mat x3Dc = (cv::Mat_<float>(3,1) << x, y, z);

        unique_lock<mutex> lock(mMutexPose);
        return Twc.rowRange(0,3).colRange(0,3)*x3Dc+Twc.rowRange(0,3).col(3);
    }
    else
        return cv::Mat();
}

float KeyFrame::ComputeSceneMedianDepth(const int q)
{
    vector<MapPoint*> vpMapPoints;
    cv::Mat Tcw_;
    {
        unique_lock<mutex> lock(mMutexFeatures);
        unique_lock<mutex> lock2(mMutexPose);
        vpMapPoints = mvpMapPoints;
        Tcw_ = Tcw.clone();
    }

    vector<float> vDepths;
    vDepths.reserve(N);
    cv::Mat Rcw2 = Tcw_.row(2).colRange(0,3);
    Rcw2 = Rcw2.t();
    float zcw = Tcw_.at<float>(2,3);
    for(int i=0; i<N; i++)
    {
        if(mvpMapPoints[i])
        {
            MapPoint* pMP = mvpMapPoints[i];
            cv::Mat x3Dw = pMP->GetWorldPos();
            float z = Rcw2.dot(x3Dw)+zcw;
            vDepths.push_back(z);
        }
    }

    sort(vDepths.begin(),vDepths.end());

    return vDepths[(vDepths.size()-1)/q];
}

//计算点特征+线特征中点的平均景深：
float KeyFrame::ComputeSceneMedianDepthBoth(const int q)
{
    vector<MapPoint*> vpMapPoints;
    cv::Mat Tcw_;
    {
        unique_lock<mutex> lock(mMutexFeatures);
        unique_lock<mutex> lock2(mMutexPose);
        vpMapPoints = mvpMapPoints;
        Tcw_ = Tcw.clone();
    }
    
    vector<MapLine*> vpMapLines;
    {
	unique_lock<mutex> lock(mMutexFeaturesLines);
	vpMapLines = mvpMapLines;
    }
    
    vector<float> vDepths;
    vDepths.reserve(N+NL);
    
    cv::Mat Rcw2 = Tcw_.row(2).colRange(0,3);
    Rcw2 = Rcw2.t();
    float zcw = Tcw_.at<float>(2,3);
    for(int i=0; i<N; i++)
    {
        if(vpMapPoints[i])
        {
            MapPoint* pMP = vpMapPoints[i];
            cv::Mat x3Dw = pMP->GetWorldPos();
            float z = Rcw2.dot(x3Dw)+zcw;
            vDepths.push_back(z);
        }
    }
    
    for(int i=0; i<NL; i++)
    {
        if(vpMapLines[i])
        {
            MapLine* pML = vpMapLines[i];
            cv::Mat x3Dw = pML->GetMidPWorldPos();
            float z = Rcw2.dot(x3Dw)+zcw;
            vDepths.push_back(z);
        }
    }

    sort(vDepths.begin(),vDepths.end());

    return vDepths[(vDepths.size()-1)/q];
}

} //namespace ORB_SLAM
