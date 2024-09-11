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

#include "Frame.h"
#include "Converter.h"
#include "ORBmatcher.h"
#include "Linematcher.h"
#include <thread>

namespace PL_SLAM
{

long unsigned int Frame::nNextId=0;
bool Frame::mbInitialComputations=true;
float Frame::cx, Frame::cy, Frame::fx, Frame::fy, Frame::invfx, Frame::invfy;
float Frame::mnMinX, Frame::mnMinY, Frame::mnMaxX, Frame::mnMaxY;
float Frame::mfGridElementWidthInv, Frame::mfGridElementHeightInv;
float Frame::mfGridElementWidthInvLines;
float Frame::mfGridElementHeightInvLines;

Frame::Frame()
{}

//Copy Constructor (LineExpanding)
Frame::Frame(const Frame &frame)
    :mpORBvocabulary(frame.mpORBvocabulary), mpORBextractorLeft(frame.mpORBextractorLeft), mpORBextractorRight(frame.mpORBextractorRight),mpLineextractor(frame.mpLineextractor),
     mTimeStamp(frame.mTimeStamp), mK(frame.mK.clone()), mDistCoef(frame.mDistCoef.clone()),
     mbf(frame.mbf), mb(frame.mb), mThDepth(frame.mThDepth), N(frame.N), mvKeys(frame.mvKeys),
     mvKeysRight(frame.mvKeysRight), mvKeysUn(frame.mvKeysUn),  mvuRight(frame.mvuRight),
     mvDepth(frame.mvDepth), mBowVec(frame.mBowVec), mFeatVec(frame.mFeatVec),
     mDescriptors(frame.mDescriptors.clone()), mDescriptorsRight(frame.mDescriptorsRight.clone()),
     mvpMapPoints(frame.mvpMapPoints), mvbOutlier(frame.mvbOutlier), mnId(frame.mnId),
     mpReferenceKF(frame.mpReferenceKF), mnScaleLevels(frame.mnScaleLevels),
     mfScaleFactor(frame.mfScaleFactor), mfLogScaleFactor(frame.mfLogScaleFactor),
     mvScaleFactors(frame.mvScaleFactors), mvInvScaleFactors(frame.mvInvScaleFactors),
     mvLevelSigma2(frame.mvLevelSigma2), mvInvLevelSigma2(frame.mvInvLevelSigma2)
{
    for(int i=0;i<FRAME_GRID_COLS;i++)
        for(int j=0; j<FRAME_GRID_ROWS; j++)
            mGrid[i][j]=frame.mGrid[i][j];

    if(!frame.mTcw.empty())
        SetPose(frame.mTcw);
    
    //线特征相关信息拷贝
    if(frame.mpLineextractor!=NULL)
    {	
      //线特征提取图像金字塔参数：
      mnScaleLevelLines = frame.mnScaleLevelLines;
      mfScaleFactorLines = frame.mfScaleFactorLines;
      mfLogScaleFactorsLines = frame.mfLogScaleFactorsLines;
      mvScaleFactorsLines = frame.mvScaleFactorsLines;
      mvInvScaleFactorsLines = frame.mvInvScaleFactorsLines;
      mvLevelSigma2Lines = frame.mvLevelSigma2Lines;
      mvInvLevelSigma2Lines = frame.mvInvLevelSigma2Lines;
      
      //线特征对应关键线、中点和描述子：
      mvLines = frame.mvLines; 
      mvMidPoints = frame.mvMidPoints;
      mDescriptorLines = frame.mDescriptorLines.clone();
      mvLinesUn = frame.mvLinesUn;
      mvMidPointsUn = frame.mvMidPointsUn;
      
      //线特征总数：
      NL = frame.NL;
      
      //当前图像帧的线参考关键帧拷贝(对应于与当前图像帧共视地图线最多的关键帧)
      mpReferenceKFLines = frame.mpReferenceKFLines;
      
      //地图线、外线信息初始化为点和线特征数量的向量：
      mvpMapLines = frame.mvpMapLines;//内含地图线中点
      mvbOutlierLines = frame.mvbOutlierLines;
      
      //线特征网格信息拷贝
      for(int i=0;i<FRAME_GRID_COLS_Lines;i++)
        for(int j=0; j<FRAME_GRID_ROWS_Lines; j++)
            mGridLines[i][j]=frame.mGridLines[i][j];
    }
}


Frame::Frame(const cv::Mat &imLeft, const cv::Mat &imRight, const double &timeStamp, ORBextractor* extractorLeft, ORBextractor* extractorRight, ORBVocabulary* voc, cv::Mat &K, cv::Mat &distCoef, const float &bf, const float &thDepth)
    :mpORBvocabulary(voc),mpORBextractorLeft(extractorLeft),mpORBextractorRight(extractorRight), mTimeStamp(timeStamp), mK(K.clone()),mDistCoef(distCoef.clone()), mbf(bf), mThDepth(thDepth),
     mpReferenceKF(static_cast<KeyFrame*>(NULL))
{
    // Frame ID
    mnId=nNextId++;

    // Scale Level Info
    mnScaleLevels = mpORBextractorLeft->GetLevels();
    mfScaleFactor = mpORBextractorLeft->GetScaleFactor();
    mfLogScaleFactor = log(mfScaleFactor);
    mvScaleFactors = mpORBextractorLeft->GetScaleFactors();
    mvInvScaleFactors = mpORBextractorLeft->GetInverseScaleFactors();
    mvLevelSigma2 = mpORBextractorLeft->GetScaleSigmaSquares();
    mvInvLevelSigma2 = mpORBextractorLeft->GetInverseScaleSigmaSquares();

    // ORB extraction
    thread threadLeft(&Frame::ExtractORB,this,0,imLeft);
    thread threadRight(&Frame::ExtractORB,this,1,imRight);
    threadLeft.join();
    threadRight.join();

    N = mvKeys.size();

    if(mvKeys.empty())
        return;

    UndistortKeyPoints();

    ComputeStereoMatches();

    mvpMapPoints = vector<MapPoint*>(N,static_cast<MapPoint*>(NULL));    
    mvbOutlier = vector<bool>(N,false);


    // This is done only for the first Frame (or after a change in the calibration)
    if(mbInitialComputations)
    {
        ComputeImageBounds(imLeft);

        mfGridElementWidthInv=static_cast<float>(FRAME_GRID_COLS)/(mnMaxX-mnMinX);
        mfGridElementHeightInv=static_cast<float>(FRAME_GRID_ROWS)/(mnMaxY-mnMinY);

        fx = K.at<float>(0,0);
        fy = K.at<float>(1,1);
        cx = K.at<float>(0,2);
        cy = K.at<float>(1,2);
        invfx = 1.0f/fx;
        invfy = 1.0f/fy;

        mbInitialComputations=false;
    }

    mb = mbf/fx;

    AssignFeaturesToGrid();
}

Frame::Frame(const cv::Mat &imGray, const cv::Mat &imDepth, const double &timeStamp, ORBextractor* extractor,ORBVocabulary* voc, cv::Mat &K, cv::Mat &distCoef, const float &bf, const float &thDepth)
    :mpORBvocabulary(voc),mpORBextractorLeft(extractor),mpORBextractorRight(static_cast<ORBextractor*>(NULL)),
     mTimeStamp(timeStamp), mK(K.clone()),mDistCoef(distCoef.clone()), mbf(bf), mThDepth(thDepth)
{
    // Frame ID
    mnId=nNextId++;

    // Scale Level Info
    mnScaleLevels = mpORBextractorLeft->GetLevels();
    mfScaleFactor = mpORBextractorLeft->GetScaleFactor();    
    mfLogScaleFactor = log(mfScaleFactor);
    mvScaleFactors = mpORBextractorLeft->GetScaleFactors();
    mvInvScaleFactors = mpORBextractorLeft->GetInverseScaleFactors();
    mvLevelSigma2 = mpORBextractorLeft->GetScaleSigmaSquares();
    mvInvLevelSigma2 = mpORBextractorLeft->GetInverseScaleSigmaSquares();

    // ORB extraction
    ExtractORB(0,imGray);

    N = mvKeys.size();

    if(mvKeys.empty())
        return;

    UndistortKeyPoints();

    ComputeStereoFromRGBD(imDepth);

    mvpMapPoints = vector<MapPoint*>(N,static_cast<MapPoint*>(NULL));
    mvbOutlier = vector<bool>(N,false);

    // This is done only for the first Frame (or after a change in the calibration)
    if(mbInitialComputations)
    {
        ComputeImageBounds(imGray);

        mfGridElementWidthInv=static_cast<float>(FRAME_GRID_COLS)/static_cast<float>(mnMaxX-mnMinX);
        mfGridElementHeightInv=static_cast<float>(FRAME_GRID_ROWS)/static_cast<float>(mnMaxY-mnMinY);

        fx = K.at<float>(0,0);
        fy = K.at<float>(1,1);
        cx = K.at<float>(0,2);
        cy = K.at<float>(1,2);
        invfx = 1.0f/fx;
        invfy = 1.0f/fy;

        mbInitialComputations=false;
    }

    mb = mbf/fx;

    AssignFeaturesToGrid();
}

//修改后的单目仅点的构造函数，添加空的线特征提取器
Frame::Frame(const cv::Mat &imGray, const double &timeStamp, ORBextractor* extractor,ORBVocabulary* voc, cv::Mat &K, cv::Mat &distCoef, const float &bf, const float &thDepth)
    :mpORBvocabulary(voc),mpORBextractorLeft(extractor),mpORBextractorRight(static_cast<ORBextractor*>(NULL)),mpLineextractor(static_cast<Lineextractor*>(NULL)),
     mTimeStamp(timeStamp), mK(K.clone()),mDistCoef(distCoef.clone()), mbf(bf), mThDepth(thDepth)
{
    // Frame ID
    mnId=nNextId++;

    // Scale Level Info
    mnScaleLevels = mpORBextractorLeft->GetLevels();
    mfScaleFactor = mpORBextractorLeft->GetScaleFactor();
    mfLogScaleFactor = log(mfScaleFactor);
    mvScaleFactors = mpORBextractorLeft->GetScaleFactors();
    mvInvScaleFactors = mpORBextractorLeft->GetInverseScaleFactors();
    mvLevelSigma2 = mpORBextractorLeft->GetScaleSigmaSquares();
    mvInvLevelSigma2 = mpORBextractorLeft->GetInverseScaleSigmaSquares();

    // ORB extraction
    ExtractORB(0,imGray);

    N = mvKeys.size();

    if(mvKeys.empty())
        return;

    UndistortKeyPoints();

    // Set no stereo information
    mvuRight = vector<float>(N,-1);
    mvDepth = vector<float>(N,-1);

    mvpMapPoints = vector<MapPoint*>(N,static_cast<MapPoint*>(NULL));
    mvbOutlier = vector<bool>(N,false);

    // This is done only for the first Frame (or after a change in the calibration)
    if(mbInitialComputations)
    {
        ComputeImageBounds(imGray);

        mfGridElementWidthInv=static_cast<float>(FRAME_GRID_COLS)/static_cast<float>(mnMaxX-mnMinX);
        mfGridElementHeightInv=static_cast<float>(FRAME_GRID_ROWS)/static_cast<float>(mnMaxY-mnMinY);

        fx = K.at<float>(0,0);
        fy = K.at<float>(1,1);
        cx = K.at<float>(0,2);
        cy = K.at<float>(1,2);
        invfx = 1.0f/fx;
        invfy = 1.0f/fy;

        mbInitialComputations=false;
    }

    mb = mbf/fx;

    AssignFeaturesToGrid();
}

//点线特征同时使用单目构造函数
Frame::Frame(const cv::Mat &imGray, const double &timeStamp, ORBextractor* extractor, ORBVocabulary* voc, Lineextractor* lineextractor ,cv::Mat &K, cv::Mat &distCoef, 
	     const float &bf, const float &thDepth):mpORBvocabulary(voc),mpORBextractorLeft(extractor),mpORBextractorRight(static_cast<ORBextractor*>(NULL)),mpLineextractor(lineextractor),
             mTimeStamp(timeStamp), mK(K.clone()),mDistCoef(distCoef.clone()), mbf(bf), mThDepth(thDepth)     
{
      // Frame ID
      mnId=nNextId++;
      
      // 点特征 Scale Level Info
      mnScaleLevels = mpORBextractorLeft->GetLevels();
      mfScaleFactor = mpORBextractorLeft->GetScaleFactor();
      mfLogScaleFactor = log(mfScaleFactor);
      mvScaleFactors = mpORBextractorLeft->GetScaleFactors();
      mvInvScaleFactors = mpORBextractorLeft->GetInverseScaleFactors();
      mvLevelSigma2 = mpORBextractorLeft->GetScaleSigmaSquares();
      mvInvLevelSigma2 = mpORBextractorLeft->GetInverseScaleSigmaSquares();

      // 线特征 Scale Level Info
      mnScaleLevelLines = mpLineextractor->GetLevels();
      mfScaleFactorLines = mpLineextractor->GetScaleFactor();
      mfLogScaleFactorsLines = log(mfScaleFactorLines);
      mvScaleFactorsLines = mpLineextractor->GetScaleFactors();
      mvInvScaleFactorsLines = mpLineextractor->GetInverseScaleFactors();
      mvLevelSigma2Lines = mpLineextractor->GetScaleSigmaSquares();
      mvInvLevelSigma2Lines = mpLineextractor->GetInverseScaleSigmaSquares();
      
      //分线程提取点线特征
      // ORB extraction
      Mat imGray1 = imGray.clone();
      Mat imGray2 = imGray.clone();
          
      if(mpLineextractor->busingLSD)
      {
	  thread threadPoints(&Frame::ExtractORB,this,0, ref(imGray1));
	  thread threadLines(&Frame::ExtractLsdWithLBD,this, ref(imGray2));
	  threadPoints.join();
	  threadLines.join();
      }  
      else
      {
	  thread threadPoints(&Frame::ExtractORB,this,0, ref(imGray1));
	  thread threadLines(&Frame::ExtractFldWithLBD,this, ref(imGray2));
	  threadPoints.join();
	  threadLines.join();  
      }

      //点线提取特征数拷贝
      N = mvKeys.size();
      NL = mvLines.size();
      
      //点线特征均存在
      if( mvKeys.empty() || mvLines.empty() )
	  return;
      
      //基于输入的相机畸变系数对关键点和关键线、关键线中点矫正
      UndistortKeyPoints();
      UndistortKeyLines();
               
      // Set no stereo information
      mvuRight = vector<float>(N,-1);
      mvDepth = vector<float>(N,-1);

      //设置地图点、外点、地图线、外线存储容器
      mvpMapPoints = vector<MapPoint*>(N,static_cast<MapPoint*>(NULL));
      mvbOutlier = vector<bool>(N,false);
      mvpMapLines = vector<MapLine*>(NL,static_cast<MapLine*>(NULL));
      mvbOutlierLines = vector<bool>(NL,false);

      // This is done only for the first Frame (or after a change in the calibration)
      if(mbInitialComputations)
      {
	  ComputeImageBounds(imGray);

	  mfGridElementWidthInv=static_cast<float>(FRAME_GRID_COLS)/static_cast<float>(mnMaxX-mnMinX);
	  mfGridElementHeightInv=static_cast<float>(FRAME_GRID_ROWS)/static_cast<float>(mnMaxY-mnMinY);
	  //计算线特征对应单位网格像素数的逆
	  mfGridElementWidthInvLines=static_cast<float>(FRAME_GRID_COLS_Lines)/static_cast<float>(mnMaxX-mnMinX);
	  mfGridElementHeightInvLines=static_cast<float>(FRAME_GRID_ROWS_Lines)/static_cast<float>(mnMaxY-mnMinY);
	  
	  fx = K.at<float>(0,0);
	  fy = K.at<float>(1,1);
	  cx = K.at<float>(0,2);
	  cy = K.at<float>(1,2);
	  invfx = 1.0f/fx;
	  invfy = 1.0f/fy;

	  mbInitialComputations=false;
      }

      mb = mbf/fx;
      
      //分别基于矫正后的关键点和关键线中点将特征分配到网格中
      AssignFeaturesToGrid();
      AssignFeaturesToGridLines();
      
}

void Frame::AssignFeaturesToGrid()
{
    int nReserve = 0.5f*N/(FRAME_GRID_COLS*FRAME_GRID_ROWS);
    for(unsigned int i=0; i<FRAME_GRID_COLS;i++)
        for (unsigned int j=0; j<FRAME_GRID_ROWS;j++)
            mGrid[i][j].reserve(nReserve);

    for(int i=0;i<N;i++)
    {
        const cv::KeyPoint &kp = mvKeysUn[i];

        int nGridPosX, nGridPosY;
        if(PosInGrid(kp,nGridPosX,nGridPosY))
            mGrid[nGridPosX][nGridPosY].push_back(i);
    }
}

//根据线特征中点位置，将关键线索引存储到对应网格中
void Frame::AssignFeaturesToGridLines()
{
    int nReserve = 0.5f*NL/(FRAME_GRID_COLS_Lines*FRAME_GRID_ROWS_Lines);
    for(unsigned int i=0; i<FRAME_GRID_COLS_Lines;i++)
        for (unsigned int j=0; j<FRAME_GRID_ROWS_Lines;j++)
            mGridLines[i][j].reserve(nReserve);

    for(int i=0;i<NL;i++)
    {
        const KeyLine &kl = mvLinesUn[i];
	const cv::KeyPoint &kp = mvMidPointsUn[i];

        int nGridPosX, nGridPosY;
        if(PosInGridLines(kl,kp,nGridPosX,nGridPosY))
            mGridLines[nGridPosX][nGridPosY].push_back(i);
    }
}


void Frame::ExtractORB(int flag, const cv::Mat &im)
{
    if(flag==0)
        (*mpORBextractorLeft)(im,cv::Mat(),mvKeys,mDescriptors);
    else
        (*mpORBextractorRight)(im,cv::Mat(),mvKeysRight,mDescriptorsRight);
}

//线特征提取器，使用ComputeLsdWithLbd()函数
void Frame::ExtractLsdWithLBD(const cv::Mat &im)
{	
    mpLineextractor->ComputeLsdWithLbd(im,mvLines,mvMidPoints,mDescriptorLines);
}

//线特征提取器，使用ComputeLsdWithLbd()函数
void Frame::ExtractFldWithLBD(cv::Mat &im)
{	
    mpLineextractor->ComputeFldWithLbd(im,mvLines,mvMidPoints,mDescriptorLines);
}

//图像帧位姿设置：最终，仅点，仅线。其中，仅点和仅线为临时存储，有效数据在mTcw中
void Frame::SetPose(cv::Mat Tcw)
{
    mTcw = Tcw.clone();
    UpdatePoseMatrices();
}
void Frame::SetPosePoints(cv::Mat TcwPoints)
{
    mTcwPoints = TcwPoints.clone();
}
void Frame::SetPoseLines(cv::Mat TcwLines)
{
    mTcwLines = TcwLines.clone();
}


void Frame::UpdatePoseMatrices()
{ 
    mRcw = mTcw.rowRange(0,3).colRange(0,3);
    mRwc = mRcw.t();
    mtcw = mTcw.rowRange(0,3).col(3);
    mOw = -mRcw.t()*mtcw;
}

bool Frame::isInFrustum(MapPoint *pMP, float viewingCosLimit)
{
    pMP->mbTrackInView = false;

    // 3D in absolute coordinates
    cv::Mat P = pMP->GetWorldPos(); 

    // 3D in camera coordinates
    const cv::Mat Pc = mRcw*P+mtcw;
    const float &PcX = Pc.at<float>(0);
    const float &PcY = Pc.at<float>(1);
    const float &PcZ = Pc.at<float>(2);

    // Check positive depth
    if(PcZ<0.0f)
        return false;

    // Project in image and check it is not outside
    const float invz = 1.0f/PcZ;
    const float u=fx*PcX*invz+cx;
    const float v=fy*PcY*invz+cy;

    if(u<mnMinX || u>mnMaxX)
        return false;
    if(v<mnMinY || v>mnMaxY)
        return false;

    // Check distance is in the scale invariance region of the MapPoint
    const float maxDistance = pMP->GetMaxDistanceInvariance();
    const float minDistance = pMP->GetMinDistanceInvariance();
    const cv::Mat PO = P-mOw;
    const float dist = cv::norm(PO);

    if(dist<minDistance || dist>maxDistance)
        return false;

   // Check viewing angle
    cv::Mat Pn = pMP->GetNormal();

    const float viewCos = PO.dot(Pn)/dist;

    if(viewCos<viewingCosLimit)
        return false;

    // Predict scale in the image
    const int nPredictedLevel = pMP->PredictScale(dist,this);

    // Data used by the tracking
    pMP->mbTrackInView = true;
    pMP->mTrackProjX = u;
    pMP->mTrackProjXR = u - mbf*invz;
    pMP->mTrackProjY = v;
    pMP->mnTrackScaleLevel= nPredictedLevel;
    pMP->mTrackViewCos = viewCos;

    return true;
}
//判断地图线特征中点是否在视野内：
bool Frame::isInFrustumLine(MapLine *pML, float viewingCosLimit)
{
    //将地图线恢复默认设置，即不进行地图跟踪
    pML->mbTrackInView = false;

    // 获取地图线中点的世界坐标
    cv::Mat P = pML->GetMidPWorldPos(); 

    // 基于初步估计位姿（跟踪运动模型或跟踪参考关键帧）计算相机坐标系下的地图线中点
    const cv::Mat Pc = mRcw*P+mtcw;
    const float &PcX = Pc.at<float>(0);
    const float &PcY = Pc.at<float>(1);
    const float &PcZ = Pc.at<float>(2);

    // Check positive depth
    if(PcZ<0.0f)
        return false;

    // Project in image and check it is not outside
    const float invz = 1.0f/PcZ;
    const float u=fx*PcX*invz+cx;
    const float v=fy*PcY*invz+cy;

    if(u<mnMinX || u>mnMaxX)
        return false;
    if(v<mnMinY || v>mnMaxY)
        return false;

    // 计算地图线中点到相机中心的距离，是否在尺度范围内
    const float maxDistance = pML->GetMaxDistanceInvariance();
    const float minDistance = pML->GetMinDistanceInvariance();
    const cv::Mat PO = P-mOw;
    const float dist = cv::norm(PO);

    if(dist<minDistance || dist>maxDistance)
        return false;

   // 计算并判断地图线中点视角和平均视角是否小于60度。
    cv::Mat Pn = pML->GetNormal();

    const float viewCos = PO.dot(Pn)/dist;

    if(viewCos<viewingCosLimit)
        return false;

    // 根据地图线中点到相机中心的距离，预测地图线中点对应尺度（对应提取层数）
    const int nPredictedLevel = pML->PredictScale(dist,this);

    // Data used by the tracking
    pML->mbTrackInView = true;
    pML->mTrackProjX = u;
    pML->mTrackProjY = v;
    pML->mnTrackScaleLevel= nPredictedLevel;
    pML->mTrackViewCos = viewCos;

    return true;
}

vector<size_t> Frame::GetFeaturesInArea(const float &x, const float  &y, const float  &r, const int minLevel, const int maxLevel) const
{
    vector<size_t> vIndices;
    vIndices.reserve(N);

    const int nMinCellX = max(0,(int)floor((x-mnMinX-r)*mfGridElementWidthInv));
    if(nMinCellX>=FRAME_GRID_COLS)
        return vIndices;

    const int nMaxCellX = min((int)FRAME_GRID_COLS-1,(int)ceil((x-mnMinX+r)*mfGridElementWidthInv));
    if(nMaxCellX<0)
        return vIndices;

    const int nMinCellY = max(0,(int)floor((y-mnMinY-r)*mfGridElementHeightInv));
    if(nMinCellY>=FRAME_GRID_ROWS)
        return vIndices;

    const int nMaxCellY = min((int)FRAME_GRID_ROWS-1,(int)ceil((y-mnMinY+r)*mfGridElementHeightInv));
    if(nMaxCellY<0)
        return vIndices;

    const bool bCheckLevels = (minLevel>0) || (maxLevel>=0);

    for(int ix = nMinCellX; ix<=nMaxCellX; ix++)
    {
        for(int iy = nMinCellY; iy<=nMaxCellY; iy++)
        {
            const vector<size_t> vCell = mGrid[ix][iy];
            if(vCell.empty())
                continue;

            for(size_t j=0, jend=vCell.size(); j<jend; j++)
            {
                const cv::KeyPoint &kpUn = mvKeysUn[vCell[j]];
                if(bCheckLevels)
                {
                    if(kpUn.octave<minLevel)
                        continue;
                    if(maxLevel>=0)
                        if(kpUn.octave>maxLevel)
                            continue;
                }

                const float distx = kpUn.pt.x-x;
                const float disty = kpUn.pt.y-y;

                if(fabs(distx)<r && fabs(disty)<r)
                    vIndices.push_back(vCell[j]);
            }
        }
    }

    return vIndices;
}
//获取当前线特征中点半径搜索网格区域的所有线特征索引：
vector<size_t> Frame::GetFeaturesInAreaLines(const float &x, const float  &y, const float  &r, const int minLevel, const int maxLevel) const
{
    vector<size_t> vIndices;
    vIndices.reserve(NL);
    //获得以线特征中点为中心，r为输入半径的最大网格宽高边界，并检测是否超出
    //向下取整
    const int nMinCellX = max(0,(int)floor((x-mnMinX-r)*mfGridElementWidthInvLines));
    if(nMinCellX>=FRAME_GRID_COLS_Lines)
        return vIndices;

    const int nMaxCellX = min((int)FRAME_GRID_COLS_Lines-1,(int)ceil((x-mnMinX+r)*mfGridElementWidthInvLines));
    if(nMaxCellX<0)
        return vIndices;

    const int nMinCellY = max(0,(int)floor((y-mnMinY-r)*mfGridElementHeightInvLines));
    if(nMinCellY>=FRAME_GRID_ROWS_Lines)
        return vIndices;

    const int nMaxCellY = min((int)FRAME_GRID_ROWS_Lines-1,(int)ceil((y-mnMinY+r)*mfGridElementHeightInvLines));
    if(nMaxCellY<0)
        return vIndices;
    
    //是否检测尺度标志
    const bool bCheckLevels = (minLevel>0) || (maxLevel>=0);
    
    
    for(int ix = nMinCellX; ix<=nMaxCellX; ix++)
    {
        for(int iy = nMinCellY; iy<=nMaxCellY; iy++)
        {
	    //修改为使用线特征网格分布
            const vector<size_t> vCell = mGridLines[ix][iy];
            if(vCell.empty())
                continue;

            for(size_t j=0, jend=vCell.size(); j<jend; j++)
            {
		//修改为使用矫正后的关键线中点位置（使用关键点容器存储）
                const cv::KeyPoint &kpUn = mvMidPointsUn[vCell[j]];
                if(bCheckLevels)
                {
                    if(kpUn.octave<minLevel)
                        continue;
                    if(maxLevel>=0)
                        if(kpUn.octave>maxLevel)
                            continue;
                }

                const float distx = kpUn.pt.x-x;
                const float disty = kpUn.pt.y-y;
		
		//这里必须使用fabs()函数，而不是abs()函数
                if(fabs(distx)<r && fabs(disty)<r)
                    vIndices.push_back(vCell[j]);
            }
        }
    }

    return vIndices;
}

bool Frame::PosInGrid(const cv::KeyPoint &kp, int &posX, int &posY)
{
    posX = round((kp.pt.x-mnMinX)*mfGridElementWidthInv);
    posY = round((kp.pt.y-mnMinY)*mfGridElementHeightInv);

    //Keypoint's coordinates are undistorted, which could cause to go out of the image
    if(posX<0 || posX>=FRAME_GRID_COLS || posY<0 || posY>=FRAME_GRID_ROWS)
        return false;

    return true;
}

//判断关键线是否在图像中，对提取关键线的两端点是否在网格中进行判断，计算线特征中点在网格中的索引并通过引用返回
bool Frame::PosInGridLines(const KeyLine &kL, const cv::KeyPoint &kp, int &posMidX1, int &posMidY1)
{
    //获取线特征首端点网格索引
    int posStartPointX = round((kL.startPointX-mnMinX)*mfGridElementWidthInvLines);
    int posStartPointY = round((kL.startPointY-mnMinY)*mfGridElementHeightInvLines);
    
    
    //KeyLine's coordinates are undistorted, which could cause to go out of the image
    if(posStartPointX<0 || posStartPointX>=FRAME_GRID_COLS_Lines || posStartPointY<0 || posStartPointY>=FRAME_GRID_ROWS_Lines)
        return false;

    //获取线特征尾端点网格索引
    int posEndPointX = round((kL.endPointX-mnMinX)*mfGridElementWidthInvLines);
    int posEndPointY = round((kL.endPointY-mnMinY)*mfGridElementHeightInvLines);
    
    
    //KeyLine's coordinates are undistorted, which could cause to go out of the image
    if(posEndPointX<0 || posEndPointX>=FRAME_GRID_COLS_Lines || posEndPointY<0 || posEndPointY>=FRAME_GRID_ROWS_Lines)
        return false;
    
    //若线特征的首尾端点均在图像内，进一步判断线特征中点是否在图像内（需要判断边界，如果存在畸变系数）
    posMidX1 = round((kp.pt.x-mnMinX)*mfGridElementWidthInvLines);
    posMidY1 = round((kp.pt.y-mnMinY)*mfGridElementHeightInvLines);
    
    
    //KeyLine's coordinates are undistorted, which could cause to go out of the image
    if(posMidX1<0 || posMidX1>=FRAME_GRID_COLS_Lines || posMidY1<0 || posMidY1>=FRAME_GRID_ROWS_Lines)
        return false;
    
    
    return true;
}

void Frame::ComputeBoW()
{
    if(mBowVec.empty())
    {
        vector<cv::Mat> vCurrentDesc = Converter::toDescriptorVector(mDescriptors);
        mpORBvocabulary->transform(vCurrentDesc,mBowVec,mFeatVec,4);
    }
}

void Frame::UndistortKeyPoints()
{
    if(mDistCoef.at<float>(0)==0.0)
    {
        mvKeysUn=mvKeys;
        return;
    }

    // Fill matrix with points
    cv::Mat mat(N,2,CV_32F);
    for(int i=0; i<N; i++)
    {
        mat.at<float>(i,0)=mvKeys[i].pt.x;
        mat.at<float>(i,1)=mvKeys[i].pt.y;
    }

    // Undistort points
    mat=mat.reshape(2);
    cv::undistortPoints(mat,mat,mK,mDistCoef,cv::Mat(),mK);
    mat=mat.reshape(1);

    // Fill undistorted keypoint vector
    mvKeysUn.resize(N);
    for(int i=0; i<N; i++)
    {
        cv::KeyPoint kp = mvKeys[i];
        kp.pt.x=mat.at<float>(i,0);
        kp.pt.y=mat.at<float>(i,1);
        mvKeysUn[i]=kp;
    }
}

//通过配置文件中的畸变参数，对线特征提取算法获得的关键线和关键线中点位置进行矫正
void Frame::UndistortKeyLines()
{
    //基于输入畸变参数判断是否需要矫正
    if(mDistCoef.at<float>(0)==0.0)
    {
        mvLinesUn=mvLines;
	mvMidPointsUn=mvMidPoints;
        return;
    }
    // step1:矫正关键线中点
    cv::Mat mat0(NL,2,CV_32F);//NL*2的单通道矩阵
    for(int i=0; i<NL; i++)
    {
        mat0.at<float>(i,0)=mvMidPoints[i].pt.x;
        mat0.at<float>(i,1)=mvMidPoints[i].pt.y;
    }
    
    // 矫正关键线中点
    mat0=mat0.reshape(2);//NL*1的双通道矩阵
    cv::undistortPoints(mat0,mat0,mK,mDistCoef,cv::Mat(),mK);
    mat0=mat0.reshape(1);//NL*2的单通道矩阵
    
    //存储矫正后的关键线中点
    mvMidPointsUn.resize(NL);
    for(int i=0; i<NL; i++)
    {
        cv::KeyPoint kp = mvMidPoints[i];//很关键，因为还要拷贝提取层数等信息
        kp.pt.x=mat0.at<float>(i,0);
        kp.pt.y=mat0.at<float>(i,1);
        mvMidPointsUn[i]=kp;
    }
    
    //矫正关键线
    cv::Mat mat1(NL,2,CV_32F);//线特征首端点NL*2的单通道矩阵
    cv::Mat mat2(NL,2,CV_32F);//线特征为尾端点NL*2的单通道矩阵
    for(int i=0; i<NL; i++)
    {
        mat1.at<float>(i,0)=mvLines[i].startPointX;
        mat1.at<float>(i,1)=mvLines[i].startPointY;	
	mat2.at<float>(i,0)=mvLines[i].endPointX;
        mat2.at<float>(i,1)=mvLines[i].endPointY;
    }
    
    mat1=mat1.reshape(2);//NL*1的双通道矩阵
    mat2=mat2.reshape(2);//NL*1的双通道矩阵
    cv::undistortPoints(mat1,mat1,mK,mDistCoef,cv::Mat(),mK);
    cv::undistortPoints(mat2,mat2,mK,mDistCoef,cv::Mat(),mK);
    mat1=mat1.reshape(1);//NL*2的单通道矩阵
    mat2=mat2.reshape(1);//NL*2的单通道矩阵
    
    mvLinesUn.resize(NL);
    for(int i=0; i<NL; i++)
    {
        KeyLine kl = mvLines[i];
	kl.startPointX = mat1.at<float>(i,0);
	kl.startPointY = mat1.at<float>(i,1);
	kl.endPointX = mat2.at<float>(i,0);
	kl.endPointY = mat2.at<float>(i,1);
        mvLinesUn[i]=kl;
    }
}

//没有使用这个新补充的函数哦～
void Frame::ComputeLines2DLengthUn()
{

    
    for(int i=0; i<NL; i++)
    {
      //与层数相关的尺度因子
      float octaveScale = pow( (float)mfScaleFactorLines, mvLinesUn[i].octave );
      
      //重新输入尺度意义下首尾端点及2D线长度
      mvLinesUn[i].sPointInOctaveX = mvLinesUn[i].startPointX / octaveScale;
      mvLinesUn[i].sPointInOctaveY = mvLinesUn[i].startPointY / octaveScale;
      mvLinesUn[i].ePointInOctaveX = mvLinesUn[i].endPointX / octaveScale;
      mvLinesUn[i].ePointInOctaveY = mvLinesUn[i].endPointY / octaveScale;
      
      mvLinesUn[i].lineLength = (float) sqrt( pow( mvLinesUn[i].sPointInOctaveX - mvLinesUn[i].ePointInOctaveX, 2 ) +
      pow( mvLinesUn[i].sPointInOctaveY - mvLinesUn[i].ePointInOctaveY, 2 ) );

    }
}


void Frame::ComputeImageBounds(const cv::Mat &imLeft)
{
    if(mDistCoef.at<float>(0)!=0.0)
    {
        cv::Mat mat(4,2,CV_32F);
        mat.at<float>(0,0)=0.0; mat.at<float>(0,1)=0.0;
        mat.at<float>(1,0)=imLeft.cols; mat.at<float>(1,1)=0.0;
        mat.at<float>(2,0)=0.0; mat.at<float>(2,1)=imLeft.rows;
        mat.at<float>(3,0)=imLeft.cols; mat.at<float>(3,1)=imLeft.rows;

        // Undistort corners
        mat=mat.reshape(2);
        cv::undistortPoints(mat,mat,mK,mDistCoef,cv::Mat(),mK);
        mat=mat.reshape(1);

        mnMinX = min(mat.at<float>(0,0),mat.at<float>(2,0));
        mnMaxX = max(mat.at<float>(1,0),mat.at<float>(3,0));
        mnMinY = min(mat.at<float>(0,1),mat.at<float>(1,1));
        mnMaxY = max(mat.at<float>(2,1),mat.at<float>(3,1));

    }
    else
    {
        mnMinX = 0.0f;
        mnMaxX = imLeft.cols;
        mnMinY = 0.0f;
        mnMaxY = imLeft.rows;
    }
}

void Frame::ComputeStereoMatches()
{
    mvuRight = vector<float>(N,-1.0f);
    mvDepth = vector<float>(N,-1.0f);

    const int thOrbDist = (ORBmatcher::TH_HIGH+ORBmatcher::TH_LOW)/2;

    const int nRows = mpORBextractorLeft->mvImagePyramid[0].rows;

    //Assign keypoints to row table
    vector<vector<size_t> > vRowIndices(nRows,vector<size_t>());

    for(int i=0; i<nRows; i++)
        vRowIndices[i].reserve(200);

    const int Nr = mvKeysRight.size();

    for(int iR=0; iR<Nr; iR++)
    {
        const cv::KeyPoint &kp = mvKeysRight[iR];
        const float &kpY = kp.pt.y;
        const float r = 2.0f*mvScaleFactors[mvKeysRight[iR].octave];
        const int maxr = ceil(kpY+r);
        const int minr = floor(kpY-r);

        for(int yi=minr;yi<=maxr;yi++)
            vRowIndices[yi].push_back(iR);
    }

    // Set limits for search
    const float minZ = mb;
    const float minD = 0;
    const float maxD = mbf/minZ;

    // For each left keypoint search a match in the right image
    vector<pair<int, int> > vDistIdx;
    vDistIdx.reserve(N);

    for(int iL=0; iL<N; iL++)
    {
        const cv::KeyPoint &kpL = mvKeys[iL];
        const int &levelL = kpL.octave;
        const float &vL = kpL.pt.y;
        const float &uL = kpL.pt.x;

        const vector<size_t> &vCandidates = vRowIndices[vL];

        if(vCandidates.empty())
            continue;

        const float minU = uL-maxD;
        const float maxU = uL-minD;

        if(maxU<0)
            continue;

        int bestDist = ORBmatcher::TH_HIGH;
        size_t bestIdxR = 0;

        const cv::Mat &dL = mDescriptors.row(iL);

        // Compare descriptor to right keypoints
        for(size_t iC=0; iC<vCandidates.size(); iC++)
        {
            const size_t iR = vCandidates[iC];
            const cv::KeyPoint &kpR = mvKeysRight[iR];

            if(kpR.octave<levelL-1 || kpR.octave>levelL+1)
                continue;

            const float &uR = kpR.pt.x;

            if(uR>=minU && uR<=maxU)
            {
                const cv::Mat &dR = mDescriptorsRight.row(iR);
                const int dist = ORBmatcher::DescriptorDistance(dL,dR);

                if(dist<bestDist)
                {
                    bestDist = dist;
                    bestIdxR = iR;
                }
            }
        }

        // Subpixel match by correlation
        if(bestDist<thOrbDist)
        {
            // coordinates in image pyramid at keypoint scale
            const float uR0 = mvKeysRight[bestIdxR].pt.x;
            const float scaleFactor = mvInvScaleFactors[kpL.octave];
            const float scaleduL = round(kpL.pt.x*scaleFactor);
            const float scaledvL = round(kpL.pt.y*scaleFactor);
            const float scaleduR0 = round(uR0*scaleFactor);

            // sliding window search
            const int w = 5;
            cv::Mat IL = mpORBextractorLeft->mvImagePyramid[kpL.octave].rowRange(scaledvL-w,scaledvL+w+1).colRange(scaleduL-w,scaleduL+w+1);
            IL.convertTo(IL,CV_32F);
            IL = IL - IL.at<float>(w,w) *cv::Mat::ones(IL.rows,IL.cols,CV_32F);

            int bestDist = INT_MAX;
            int bestincR = 0;
            const int L = 5;
            vector<float> vDists;
            vDists.resize(2*L+1);

            const float iniu = scaleduR0+L-w;
            const float endu = scaleduR0+L+w+1;
            if(iniu<0 || endu >= mpORBextractorRight->mvImagePyramid[kpL.octave].cols)
                continue;

            for(int incR=-L; incR<=+L; incR++)
            {
                cv::Mat IR = mpORBextractorRight->mvImagePyramid[kpL.octave].rowRange(scaledvL-w,scaledvL+w+1).colRange(scaleduR0+incR-w,scaleduR0+incR+w+1);
                IR.convertTo(IR,CV_32F);
                IR = IR - IR.at<float>(w,w) *cv::Mat::ones(IR.rows,IR.cols,CV_32F);

                float dist = cv::norm(IL,IR,cv::NORM_L1);
                if(dist<bestDist)
                {
                    bestDist =  dist;
                    bestincR = incR;
                }

                vDists[L+incR] = dist;
            }

            if(bestincR==-L || bestincR==L)
                continue;

            // Sub-pixel match (Parabola fitting)
            const float dist1 = vDists[L+bestincR-1];
            const float dist2 = vDists[L+bestincR];
            const float dist3 = vDists[L+bestincR+1];

            const float deltaR = (dist1-dist3)/(2.0f*(dist1+dist3-2.0f*dist2));

            if(deltaR<-1 || deltaR>1)
                continue;

            // Re-scaled coordinate
            float bestuR = mvScaleFactors[kpL.octave]*((float)scaleduR0+(float)bestincR+deltaR);

            float disparity = (uL-bestuR);

            if(disparity>=minD && disparity<maxD)
            {
                if(disparity<=0)
                {
                    disparity=0.01;
                    bestuR = uL-0.01;
                }
                mvDepth[iL]=mbf/disparity;
                mvuRight[iL] = bestuR;
                vDistIdx.push_back(pair<int,int>(bestDist,iL));
            }
        }
    }

    sort(vDistIdx.begin(),vDistIdx.end());
    const float median = vDistIdx[vDistIdx.size()/2].first;
    const float thDist = 1.5f*1.4f*median;

    for(int i=vDistIdx.size()-1;i>=0;i--)
    {
        if(vDistIdx[i].first<thDist)
            break;
        else
        {
            mvuRight[vDistIdx[i].second]=-1;
            mvDepth[vDistIdx[i].second]=-1;
        }
    }
}


void Frame::ComputeStereoFromRGBD(const cv::Mat &imDepth)
{
    mvuRight = vector<float>(N,-1);
    mvDepth = vector<float>(N,-1);

    for(int i=0; i<N; i++)
    {
        const cv::KeyPoint &kp = mvKeys[i];
        const cv::KeyPoint &kpU = mvKeysUn[i];

        const float &v = kp.pt.y;
        const float &u = kp.pt.x;

        const float d = imDepth.at<float>(v,u);

        if(d>0)
        {
            mvDepth[i] = d;
            mvuRight[i] = kpU.pt.x-mbf/d;
        }
    }
}

cv::Mat Frame::UnprojectStereo(const int &i)
{
    const float z = mvDepth[i];
    if(z>0)
    {
        const float u = mvKeysUn[i].pt.x;
        const float v = mvKeysUn[i].pt.y;
        const float x = (u-cx)*z*invfx;
        const float y = (v-cy)*z*invfy;
        cv::Mat x3Dc = (cv::Mat_<float>(3,1) << x, y, z);
        return mRwc*x3Dc+mOw;
    }
    else
        return cv::Mat();
}




} //namespace ORB_SLAM
