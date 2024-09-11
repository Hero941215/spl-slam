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

#include "Linematcher.h"

#include<limits.h>

#include<opencv2/core/core.hpp>
#include<opencv2/features2d/features2d.hpp>

#include<stdint.h>

#include <thread> 

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense> 
using namespace Eigen;

using namespace std;

namespace PL_SLAM
{

const int Linematcher::TH_HIGH = 100;
const int Linematcher::TH_LOW = 50;
const int Linematcher::HISTO_LENGTH = 30;

//线特征匹配器构造函数：
Linematcher::Linematcher(float nnratio, bool checkOri, bool checklen, float lengtherr):
	    mfNNratio(nnratio),mbCheckOrientation(checkOri),mbchecklen(checklen),mflengtherr(lengtherr)
{
}

//线特征的二进制描述子比较函数：
int Linematcher::DescriptorDistance(const cv::Mat &a, const cv::Mat &b)
{
    const int *pa = a.ptr<int32_t>();
    const int *pb = b.ptr<int32_t>();

    int dist=0;

    for(int i=0; i<8; i++, pa++, pb++)
    {
        unsigned  int v = *pa ^ *pb;
        v = v - ((v >> 1) & 0x55555555);
        v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
        dist += (((v + (v >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
    }

    return dist;
}

//最佳旋转方向直方图计算
void Linematcher::ComputeThreeMaxima(vector<int>* histo, const int L, int &ind1, int &ind2, int &ind3)
{
    int max1=0;
    int max2=0;
    int max3=0;

    for(int i=0; i<L; i++)
    {
        const int s = histo[i].size();
        if(s>max1)
        {
            max3=max2;
            max2=max1;
            max1=s;
            ind3=ind2;
            ind2=ind1;
            ind1=i;
        }
        else if(s>max2)
        {
            max3=max2;
            max2=s;
            ind3=ind2;
            ind2=i;
        }
        else if(s>max3)
        {
            max3=s;
            ind3=i;
        }
    }

    if(max2<0.1f*(float)max1)
    {
        ind2=-1;
        ind3=-1;
    }
    else if(max3<0.1f*(float)max1)
    {
        ind3=-1;
    }
}

float Linematcher::RadiusByViewingCos(const float &viewCos)
{
    if(viewCos>0.998)
        return 2.5;
    else
        return 4.0;
}

//相较点特征，线特征中点和核线约束取，单自由度的显著性水平为0.1
bool Linematcher::CheckDistEpipolarLine(const cv::KeyPoint &kp1,const cv::KeyPoint &kp2,const cv::Mat &F12,const KeyFrame* pKF2)
{
    // Epipolar line in second image l = x1'F12 = [a b c]
    const float a = kp1.pt.x*F12.at<float>(0,0)+kp1.pt.y*F12.at<float>(1,0)+F12.at<float>(2,0);
    const float b = kp1.pt.x*F12.at<float>(0,1)+kp1.pt.y*F12.at<float>(1,1)+F12.at<float>(2,1);
    const float c = kp1.pt.x*F12.at<float>(0,2)+kp1.pt.y*F12.at<float>(1,2)+F12.at<float>(2,2);

    const float num = a*kp2.pt.x+b*kp2.pt.y+c;

    const float den = a*a+b*b;

    if(den==0)
        return false;

    const float dsqr = num*num/den;
    
//     float xxx = 0.0;
//     xxx = dsqr / pKF2->mvLevelSigma2Lines[kp2.octave];
//     cout << "xxx" << xxx <<endl<<endl<<endl;

    //return dsqr<3.841*pKF2->mvLevelSigma2Lines[kp2.octave];
    return dsqr<3.841*pKF2->mvLevelSigma2Lines[kp2.octave];
}

//单目初始化所需的线特征匹配函数
int Linematcher::SearchForInitialization(Frame &F1, Frame &F2, std::vector<cv::Point2f> &vbPrevMatched, 
			    std::vector<int> &vnMatches12, int windowSize)
{
    int nmatches = 0;
    
    vnMatches12 = vector<int>(F1.mvMidPointsUn.size(),-1);

    vector<int> rotHist[HISTO_LENGTH];
    for(int i=0;i<HISTO_LENGTH;i++)
        rotHist[i].reserve(500);
    const float factor = 1.0f/HISTO_LENGTH;

    vector<int> vMatchedDistance(F2.mvMidPointsUn.size(),INT_MAX);
    vector<int> vnMatches21(F2.mvMidPointsUn.size(),-1);

    
    for(size_t i1=0, iend1=F1.mvMidPointsUn.size(); i1<iend1; i1++)
    {
	
        cv::KeyPoint kp1 = F1.mvMidPointsUn[i1];
        int level1 = kp1.octave;
        if(level1>0)
            continue;

        vector<size_t> vIndices2 = F2.GetFeaturesInAreaLines(vbPrevMatched[i1].x,vbPrevMatched[i1].y, windowSize,level1,level1);

        if(vIndices2.empty())
            continue;

        cv::Mat d1 = F1.mDescriptorLines.row(i1);

        int bestDist = INT_MAX;
        int bestDist2 = INT_MAX;
        int bestIdx2 = -1;
	
	//获得候选线特征进行遍历
        for(vector<size_t>::iterator vit=vIndices2.begin(); vit!=vIndices2.end(); vit++)
        {
            size_t i2 = *vit;

            cv::Mat d2 = F2.mDescriptorLines.row(i2);

            int dist = DescriptorDistance(d1,d2);
	    
            if(vMatchedDistance[i2]<=dist)
                continue;

            if(dist<bestDist)
            {
                bestDist2=bestDist;
                bestDist=dist;
                bestIdx2=i2;
            }
            else if(dist<bestDist2)
            {
                bestDist2=dist;
            }
        }

        // 判断是否描述子距离小于低阈值（50）
        if(bestDist<=TH_LOW)
        {

            if(bestDist<(float)bestDist2*mfNNratio)
            {

                if(vnMatches21[bestIdx2]>=0)
                {
                    vnMatches12[vnMatches21[bestIdx2]]=-1;
                    nmatches--;
                }

                vnMatches12[i1]=bestIdx2;
                vnMatches21[bestIdx2]=i1;
                vMatchedDistance[bestIdx2]=bestDist;
                nmatches++;
		

                if(mbCheckOrientation)
                {

                    float rot = F1.mvLinesUn[i1].angle-F2.mvLinesUn[bestIdx2].angle;
                    if(rot<0.0)
                        rot+=360.0f;
                    int bin = round(rot*factor);
                    if(bin==HISTO_LENGTH)
                        bin=0;
                    assert(bin>=0 && bin<HISTO_LENGTH);
                    rotHist[bin].push_back(i1);
                }
            }
        }

    }

    if(mbCheckOrientation)
    {
        int ind1=-1;
        int ind2=-1;
        int ind3=-1;

        ComputeThreeMaxima(rotHist,HISTO_LENGTH,ind1,ind2,ind3);

        for(int i=0; i<HISTO_LENGTH; i++)
        {
            if(i==ind1 || i==ind2 || i==ind3)
                continue;
            for(size_t j=0, jend=rotHist[i].size(); j<jend; j++)
            {
                int idx1 = rotHist[i][j];
                if(vnMatches12[idx1]>=0)
                {
                    vnMatches12[idx1]=-1;
                    nmatches--;
                }
            }
        }

    }

    if(mbchecklen)
    {
       for(size_t i1=0, iend1=vnMatches12.size(); i1<iend1; i1++)

	  if(vnMatches12[i1]>=0)

	    if( ((1-mflengtherr)*F1.mvLinesUn[i1].lineLength>F2.mvLinesUn[vnMatches12[i1]].lineLength )
	       || ( (1+mflengtherr)*F1.mvLinesUn[i1].lineLength<F2.mvLinesUn[vnMatches12[i1]].lineLength) )
	    {
		vnMatches12[i1]=-1;
                nmatches--;
            }
	   
    }
    
    for(size_t i1=0, iend1=vnMatches12.size(); i1<iend1; i1++)
        if(vnMatches12[i1]>=0)
            vbPrevMatched[i1]=F2.mvMidPointsUn[vnMatches12[i1]].pt;

    return nmatches;
}

//跟踪运动模型的线特征匹配函数，仅适用于单目，未考虑双目情况
int Linematcher::SearchByProjection(Frame &CurrentFrame, const Frame &LastFrame, const float th, const bool bMono)
{
    int nmatches = 0;
  
    vector<int> rotHist[HISTO_LENGTH];
    for(int i=0;i<HISTO_LENGTH;i++)
        rotHist[i].reserve(500);
    const float factor = 1.0f/HISTO_LENGTH;
    
    const cv::Mat Rcw = CurrentFrame.mTcw.rowRange(0,3).colRange(0,3);
    const cv::Mat tcw = CurrentFrame.mTcw.rowRange(0,3).col(3);
    
    for(int i=0; i<LastFrame.NL; i++)
    {
        MapLine* pML = LastFrame.mvpMapLines[i];
	
        if(pML)
        {
            if(!LastFrame.mvbOutlierLines[i])
            {

                cv::Mat x3Dw = pML->GetMidPWorldPos();
                cv::Mat x3Dc = Rcw*x3Dw+tcw;

                const float xc = x3Dc.at<float>(0);
                const float yc = x3Dc.at<float>(1);
                const float invzc = 1.0/x3Dc.at<float>(2);

                if(invzc<0)
                    continue;
		
                float u = CurrentFrame.fx*xc*invzc+CurrentFrame.cx;
                float v = CurrentFrame.fy*yc*invzc+CurrentFrame.cy;
		
                if(u<CurrentFrame.mnMinX || u>CurrentFrame.mnMaxX)
                    continue;
                if(v<CurrentFrame.mnMinY || v>CurrentFrame.mnMaxY)
                    continue;
		
                int nLastOctave = LastFrame.mvMidPointsUn[i].octave;

                // Search in a window. Size depends on scale
                float radius = th*CurrentFrame.mvScaleFactorsLines[nLastOctave]; 

		vector<size_t> vIndices2;
                vIndices2 = CurrentFrame.GetFeaturesInAreaLines(u,v, radius, nLastOctave-1, nLastOctave+1);

                if(vIndices2.empty())
                    continue;
		
                const cv::Mat dML = pML->GetDescriptor();

                int bestDist = 256;
                int bestIdx2 = -1;

                for(vector<size_t>::const_iterator vit=vIndices2.begin(), vend=vIndices2.end(); vit!=vend; vit++)
                {

                    const size_t i2 = *vit;
                    if(CurrentFrame.mvpMapLines[i2])
                        if(CurrentFrame.mvpMapLines[i2]->Observations()>0)
                            continue;

                    const cv::Mat &d = CurrentFrame.mDescriptorLines.row(i2);
		    
		    int dist = DescriptorDistance(dML,d);		    
		    
                    if(dist<bestDist)
                    {
                        bestDist=dist;
                        bestIdx2=i2;
                    }
                }

                // 大于高阈值100
                if(bestDist<=TH_HIGH)
                {
                    CurrentFrame.mvpMapLines[bestIdx2]=pML;
                    nmatches++;

                    if(mbCheckOrientation)
                    {
                        float rot = LastFrame.mvLinesUn[i].angle-CurrentFrame.mvLinesUn[bestIdx2].angle;
                        if(rot<0.0)
                            rot+=360.0f;
                        int bin = round(rot*factor);
                        if(bin==HISTO_LENGTH)
                            bin=0;
                        assert(bin>=0 && bin<HISTO_LENGTH);
                        rotHist[bin].push_back(bestIdx2);
                    }
                }
            }
        }
    }

    //Apply rotation consistency
    if(mbCheckOrientation)
    {
        int ind1=-1;
        int ind2=-1;
        int ind3=-1;

        ComputeThreeMaxima(rotHist,HISTO_LENGTH,ind1,ind2,ind3);

        for(int i=0; i<HISTO_LENGTH; i++)
        {
            if(i!=ind1 && i!=ind2 && i!=ind3)
            {
                for(size_t j=0, jend=rotHist[i].size(); j<jend; j++)
                {
                    CurrentFrame.mvpMapLines[rotHist[i][j]]=static_cast<MapLine*>(NULL);
                    nmatches--;
                }
            }
        }
    }   
    
    //检测地图线中存储的平均2D长度与当前帧线特征长度的误差
    if(mbchecklen)
    {

	for(int i=0; i<CurrentFrame.NL; i++)
	{
	  
	    MapLine* pML = CurrentFrame.mvpMapLines[i];

	    if(pML)
	    {	
		//获取地图线中存储的平均2D长度    
		float LineAverageLength = pML->Get2DLineLengthAverage();;  

		if( ( (1-mflengtherr)*LineAverageLength>CurrentFrame.mvLinesUn[i].lineLength )
		     || ( (1+mflengtherr)*LineAverageLength<CurrentFrame.mvLinesUn[i].lineLength) )
		{
		      CurrentFrame.mvpMapLines[i]=static_cast<MapLine*>(NULL);
		      nmatches--;
		}
		
	    }
	    
	}
    }
    
    return nmatches;

}

int Linematcher::SearchByKNN(KeyFrame *pKF, Frame &F, std::vector<MapLine*> &vpMapLineMatches)
{

    int nmatches12 = 0;
    int nmatches21 = 0;
    
    cv::Mat desc1 = pKF->mDescriptorLines;
    cv::Mat desc2 = F.mDescriptorLines;
    
    //正向匹配记录（pKF --> F）
    std::vector<int> matches_12;
    std::vector<int> matches_21;
    
    const vector<MapLine*> vpMapLinesKF = pKF->GetMapLineMatches();
    vpMapLineMatches = vector<MapLine*>(F.NL,static_cast<MapLine*>(NULL));
    
    //调用OPENCV中的KNN算法，分线程执行双向匹配
    thread threadKFToF(&Linematcher::matchNNR, this, ref(desc1), ref(desc2), mfNNratio, ref(matches_12), ref(nmatches12) );
    thread threadFToKF(&Linematcher::matchNNR, this, ref(desc2), ref(desc1), mfNNratio, ref(matches_21), ref(nmatches21));
    threadKFToF.join();
    threadFToKF.join();
    
    //通过反向匹配剔除误匹配
    for (size_t i1 = 0; i1 < matches_12.size(); ++i1)
    {   
         int &i2 = matches_12[i1];
	 
	 size_t Matches_21Now = matches_21[i2];
	 
         if (i2 >= 0 && Matches_21Now != i1)
	 {    
	   i2 = -1;
           nmatches12--;
         }  
    }

    //使用C++独有的引用
    if(mbchecklen)
    {
	
	for (size_t i1 = 0; i1 < matches_12.size(); ++i1)
	{
	    int &i2 = matches_12[i1];

	    if (i2 >= 0)
	    {	
		
		MapLine* pML = vpMapLinesKF[i1];
	    
		if(!pML)
		{
		    i2 = -1;
		    nmatches12--;
		    continue;
		}	    
		if(pML->isBad())
		{
		    i2 = -1;
		    nmatches12--;
		    continue;
		}
   
		float LineAverageLength = pML->Get2DLineLengthAverage();

		if( ( (1-mflengtherr)*LineAverageLength<=F.mvLinesUn[i2].lineLength )
		      || ( (1+mflengtherr)*LineAverageLength>=F.mvLinesUn[i2].lineLength) )
		{
		    vpMapLineMatches[i2] = pML;
		}
		else
		{
		    i2 = -1;
		    nmatches12--;
		}
	    }
	}
    }
    
    return nmatches12;
    
}

//SearchForTriangulation()和 SearchByKNN()调用的KNN匹配子函数，内部除了匹配，还完成了最优和次优描述子距离比较：
void Linematcher::matchNNR(const cv::Mat &desc1, const cv::Mat &desc2, float nnr, std::vector<int> &matches_12,int &nmatches )              
{
  
    matches_12.resize(desc1.rows, -1);

    std::vector<std::vector<cv::DMatch>> matches_;
    cv::Ptr<cv::BFMatcher> bfm = cv::BFMatcher::create(cv::NORM_HAMMING, false); // cross-check
    bfm->knnMatch(desc1, desc2, matches_, 2);

    size_t Desc1RowsNow = desc1.rows;
    
    if(Desc1RowsNow != matches_.size())
        throw std::runtime_error("[matchNNR] Different size for matches and descriptors!");

    for(int idx = 0; idx < desc1.rows; ++idx) {
        if(matches_[idx][0].distance < matches_[idx][1].distance * nnr) {
            matches_12[idx] = matches_[idx][0].trainIdx;
            nmatches++;
        }
    }
    
}

//局部地图跟踪函数
int Linematcher::SearchByProjection(Frame &F, const std::vector<MapLine*> &vpMapLines, const float th)
{
  
    int nmatches=0;
    
    const bool bFactor = th!=1.0;

    for(size_t iML=0; iML<vpMapLines.size(); iML++)
    {
        MapLine* pML = vpMapLines[iML];

        if(!pML->mbTrackInView)
            continue;

        if(pML->isBad())
            continue;

        const int &nPredictedLevel = pML->mnTrackScaleLevel;

        // The size of the window will depend on the viewing direction
        float r = RadiusByViewingCos(pML->mTrackViewCos);
        
        if(bFactor)
            r*=th;

        const vector<size_t> vIndices =
                F.GetFeaturesInAreaLines(pML->mTrackProjX,pML->mTrackProjY,r*F.mvScaleFactorsLines[nPredictedLevel],nPredictedLevel-1,nPredictedLevel);

        if(vIndices.empty())
            continue;

        const cv::Mat MPdescriptor = pML->GetDescriptor();

        int bestDist=256;
        int bestLevel= -1; 
        int bestDist2=256;
        int bestLevel2 = -1;
        int bestIdx =-1 ;

        // Get best and second matches with near keypoints
        for(vector<size_t>::const_iterator vit=vIndices.begin(), vend=vIndices.end(); vit!=vend; vit++)
        {
            const size_t idx = *vit;

            if(F.mvpMapLines[idx])
                if(F.mvpMapLines[idx]->Observations()>0)
                    continue;

            const cv::Mat &d = F.mDescriptorLines.row(idx);

            const int dist = DescriptorDistance(MPdescriptor,d);

            if(dist<bestDist)
            {
                bestDist2=bestDist;
                bestDist=dist;
                bestLevel2 = bestLevel;
                bestLevel = F.mvMidPointsUn[idx].octave;
                bestIdx=idx;
            }
            else if(dist<bestDist2)
            {
                bestLevel2 = F.mvMidPointsUn[idx].octave;
                bestDist2=dist;
            }
        }

        // Apply ratio to second match (only if best and second are in the same scale level)
        if(bestDist<=TH_HIGH)
        {
            if(bestLevel==bestLevel2 && bestDist>mfNNratio*bestDist2)
                continue;

            F.mvpMapLines[bestIdx]=pML; 
            nmatches++;
        }
         
    }

    if(mbchecklen)
    {
	for(int i=0; i<F.NL; i++)
	{  
	    MapLine* pML = F.mvpMapLines[i];

	    if(pML)
	    {
		float LineAverageLength = pML->Get2DLineLengthAverage();

		if( ( (1-mflengtherr)*LineAverageLength>F.mvLinesUn[i].lineLength )
		  || ( (1+mflengtherr)*LineAverageLength<F.mvLinesUn[i].lineLength) )
		{
		    F.mvpMapLines[i]=static_cast<MapLine*>(NULL);
		    nmatches--;
		}
	    }
	    
	}
    } 
    
    return nmatches;
    
}

int Linematcher::SearchByProjection(Frame &CurrentFrame, KeyFrame* pKF, const std::set<MapLine*> &sAlreadyFound, const float th, const int Linedist)
{
    int nmatches = 0;

    const cv::Mat Rcw = CurrentFrame.mTcw.rowRange(0,3).colRange(0,3);
    const cv::Mat tcw = CurrentFrame.mTcw.rowRange(0,3).col(3);
    const cv::Mat Ow = -Rcw.t()*tcw;

    // Rotation Histogram (to check rotation consistency)
    vector<int> rotHist[HISTO_LENGTH];
    for(int i=0;i<HISTO_LENGTH;i++)
        rotHist[i].reserve(500);
    const float factor = 1.0f/HISTO_LENGTH;
    
    const vector<MapLine*> vpMLs = pKF->GetMapLineMatches();

    for(size_t i=0, iend=vpMLs.size(); i<iend; i++)
    {
        MapLine* pML = vpMLs[i];

        if(pML)
        {	

            if(!pML->isBad() && !sAlreadyFound.count(pML))
            {

                cv::Mat x3Dw = pML->GetMidPWorldPos();
                cv::Mat x3Dc = Rcw*x3Dw+tcw;

                const float xc = x3Dc.at<float>(0);
                const float yc = x3Dc.at<float>(1);
                const float invzc = 1.0/x3Dc.at<float>(2);

                const float u = CurrentFrame.fx*xc*invzc+CurrentFrame.cx;
                const float v = CurrentFrame.fy*yc*invzc+CurrentFrame.cy;

                if(u<CurrentFrame.mnMinX || u>CurrentFrame.mnMaxX)
                    continue;
                if(v<CurrentFrame.mnMinY || v>CurrentFrame.mnMaxY)
                    continue;

                // Compute predicted scale level
                cv::Mat PO = x3Dw-Ow;
                float dist3D = cv::norm(PO);

                const float maxDistance = pML->GetMaxDistanceInvariance();
                const float minDistance = pML->GetMinDistanceInvariance();

                // Depth must be inside the scale pyramid of the image
                if(dist3D<minDistance || dist3D>maxDistance)
                    continue;

                int nPredictedLevel = pML->PredictScale(dist3D,&CurrentFrame);

                // Search in a window
                const float radius = th*CurrentFrame.mvScaleFactorsLines[nPredictedLevel];

                const vector<size_t> vIndices2 = CurrentFrame.GetFeaturesInAreaLines(u, v, radius, nPredictedLevel-1, nPredictedLevel+1);

                if(vIndices2.empty())
                    continue;

                const cv::Mat dML = pML->GetDescriptor();

                int bestDist = 256;
                int bestIdx2 = -1;

                for(vector<size_t>::const_iterator vit=vIndices2.begin(); vit!=vIndices2.end(); vit++)
                {
                    const size_t i2 = *vit;
                    if(CurrentFrame.mvpMapLines[i2])
                        continue;

                    const cv::Mat &d = CurrentFrame.mDescriptorLines.row(i2);

                    const int dist = DescriptorDistance(dML,d);

                    if(dist<bestDist)
                    {
                        bestDist=dist;
                        bestIdx2=i2;
                    }
                }

                if(bestDist<=Linedist)
                {
                    CurrentFrame.mvpMapLines[bestIdx2]=pML;
                    nmatches++;

                    if(mbCheckOrientation)
                    {
                        float rot = pKF->mvLinesUn[i].angle-CurrentFrame.mvLinesUn[bestIdx2].angle;
                        if(rot<0.0)
                            rot+=360.0f;
                        int bin = round(rot*factor);
                        if(bin==HISTO_LENGTH)
                            bin=0;
                        assert(bin>=0 && bin<HISTO_LENGTH);
                        rotHist[bin].push_back(bestIdx2);
                    }
                }

            }
        }
    }

    if(mbCheckOrientation)
    {
        int ind1=-1;
        int ind2=-1;
        int ind3=-1;

        ComputeThreeMaxima(rotHist,HISTO_LENGTH,ind1,ind2,ind3);

        for(int i=0; i<HISTO_LENGTH; i++)
        {
            if(i!=ind1 && i!=ind2 && i!=ind3)
            {
                for(size_t j=0, jend=rotHist[i].size(); j<jend; j++)
                {
                    CurrentFrame.mvpMapLines[rotHist[i][j]]=static_cast<MapLine*>(NULL);
                    nmatches--;
                }
            }
        }
    }

    if(mbchecklen)
    {

	for(int i=0; i<CurrentFrame.NL; i++)
	{
	    MapLine* pML = CurrentFrame.mvpMapLines[i];

	    if(pML)
	    {
  
		float LineAverageLength = pML->Get2DLineLengthAverage();

		if( ( (1-mflengtherr)*LineAverageLength>CurrentFrame.mvLinesUn[i].lineLength )
		    || ( (1+mflengtherr)*LineAverageLength<CurrentFrame.mvLinesUn[i].lineLength) )
		{
		     CurrentFrame.mvpMapLines[i]=static_cast<MapLine*>(NULL);
		     nmatches--;
		}
		
	    }
	    
	}
	
    }
    
    return nmatches;
  
}

int Linematcher::SearchForTriangulation(KeyFrame *pKF1, KeyFrame* pKF2, cv::Mat F12, 
			   std::vector<pair<size_t, size_t> > &vMatchedPairs)
{
    //Compute epipole in second image
    cv::Mat Cw = pKF1->GetCameraCenter();
    cv::Mat R2w = pKF2->GetRotation();
    cv::Mat t2w = pKF2->GetTranslation();
    cv::Mat C2 = R2w*Cw+t2w;
    const float invz = 1.0f/C2.at<float>(2);
    const float ex =pKF2->fx*C2.at<float>(0)*invz+pKF2->cx;
    const float ey =pKF2->fy*C2.at<float>(1)*invz+pKF2->cy;

    cv::Mat desc1 = pKF1->mDescriptorLines;
    cv::Mat desc2 = pKF2->mDescriptorLines;
    
    std::vector<int> matches_12;
    std::vector<int> matches_21;
    
    int nmatches12 = 0;
    int nmatches21 = 0;
    //调用OPENCV中的KNN算法，分线程执行双向匹配
    thread threadKFToF(&Linematcher::matchNNR, this, ref(desc1), ref(desc2), mfNNratio, ref(matches_12), ref(nmatches12) );
    thread threadFToKF(&Linematcher::matchNNR, this, ref(desc2), ref(desc1), mfNNratio, ref(matches_21), ref(nmatches21));
    threadKFToF.join();
    threadFToKF.join();

    const vector<MapLine*> vpMapLinesKF1 = pKF1->GetMapLineMatches();
    const vector<MapLine*> vpMapLinesKF2 = pKF2->GetMapLineMatches();
    
    for (size_t i1 = 0; i1 < matches_12.size(); ++i1)
    {   
         int &i2 = matches_12[i1];
	    
	 size_t Matches_12Now = matches_21[i2];
	 
         if (i2 >= 0 &&  Matches_12Now != i1 && vpMapLinesKF1[i1] && vpMapLinesKF2[i2] )
	 {    
	   i2 = -1;
           nmatches12--;
         }  
    }

     for (size_t i1 = 0; i1 < matches_12.size(); ++i1)
    {

	int &i2 = matches_12[i1];

	const cv::KeyPoint &kp1 = pKF1->mvMidPointsUn[i1];
	const cv::KeyPoint &kp2 = pKF2->mvMidPointsUn[i2];
	
	const float distex = ex-kp2.pt.x;
        const float distey = ey-kp2.pt.y;

        if(distex*distex+distey*distey<100*pKF2->mvScaleFactorsLines[kp2.octave])
	{
	    i2 = -1;
            nmatches12--;
	}

        if(!CheckDistEpipolarLine(kp1,kp2,F12,pKF2))
        {
            i2 = -1;
            nmatches12--;
        }                   
    }

    for (size_t i1 = 0; i1 < matches_12.size(); ++i1)
    {
	if(matches_12[i1]<0)
	    continue;
	
	vMatchedPairs.push_back(make_pair(i1, matches_12[i1]));
    }
    
    return nmatches12;
}

int Linematcher::Fuse(KeyFrame* pKF, const vector<MapLine*> &vpMapLines, const float th)
{

    cv::Mat Rcw = pKF->GetRotation();
    cv::Mat tcw = pKF->GetTranslation();

    const float &fx = pKF->fx;
    const float &fy = pKF->fy;
    const float &cx = pKF->cx;
    const float &cy = pKF->cy;

    cv::Mat Ow = pKF->GetCameraCenter();

    int nFused=0;

    const int nMLs = vpMapLines.size();

    for(int i=0; i<nMLs; i++)
    {
        MapLine* pML = vpMapLines[i];

        if(!pML)
            continue;

        if(pML->isBad() || pML->IsInKeyFrame(pKF))
            continue;

        cv::Mat p3Dw = pML->GetMidPWorldPos();
        cv::Mat p3Dc = Rcw*p3Dw + tcw;

        // Depth must be positive
        if(p3Dc.at<float>(2)<0.0f)
            continue;

        const float invz = 1/p3Dc.at<float>(2);
        const float x = p3Dc.at<float>(0)*invz;
        const float y = p3Dc.at<float>(1)*invz;

        const float u = fx*x+cx;
        const float v = fy*y+cy;

        // Point must be inside the image
        if(!pKF->IsInImage(u,v))
            continue;

        const float maxDistance = pML->GetMaxDistanceInvariance();
        const float minDistance = pML->GetMinDistanceInvariance();
        cv::Mat PO = p3Dw-Ow;
        const float dist3D = cv::norm(PO);

        // Depth must be inside the scale pyramid of the image
        if(dist3D<minDistance || dist3D>maxDistance )
            continue;

        // Viewing angle must be less than 60 deg
        cv::Mat Pn = pML->GetNormal();

	//将线的改为90度角
        if(PO.dot(Pn)<0*dist3D)
            continue;

        int nPredictedLevel = pML->PredictScale(dist3D,pKF);

        // Search in a radius
        const float radius = th*pKF->mvScaleFactorsLines[nPredictedLevel];

        const vector<size_t> vIndices = pKF->GetFeaturesInAreaLines(u,v,radius);

        if(vIndices.empty())
            continue;

        // Match to the most similar keypoint in the radius

        const cv::Mat dML = pML->GetDescriptor();

        int bestDist = 256;
        int bestIdx = -1;
        for(vector<size_t>::const_iterator vit=vIndices.begin(), vend=vIndices.end(); vit!=vend; vit++)
        {
	  
            const size_t idx = *vit;
            const KeyLine &kL = pKF->mvLinesUn[idx];
            const int &kLLevel= kL.octave;

            if(kLLevel<nPredictedLevel-1 || kLLevel>nPredictedLevel)
                continue;

	    Vector3f sp_l; sp_l << kL.startPointX, kL.startPointY, 1.0;
	    Vector3f ep_l; ep_l << kL.endPointX,   kL.endPointY,   1.0;
	    Vector3f le_l; le_l << sp_l.cross(ep_l);
	    //le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) + le_l(2)*le_l(2));
	    le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1));
            const float error = le_l(0)*u + le_l(1)*v +le_l(2);
            if(fabs(error)*pKF->mvInvLevelSigma2Lines[kLLevel]>3.841)
                continue;

            const cv::Mat &dKF = pKF->mDescriptorLines.row(idx);

            const int dist = DescriptorDistance(dML,dKF);

            if(dist<bestDist)
            {
                bestDist = dist;
                bestIdx = idx;
            }
        }

        // If there is already a MapPoint replace otherwise add new measurement
        if(bestDist<=TH_LOW)
        {
            MapLine* pMLinKF = pKF->GetMapLine(bestIdx);
            if(pMLinKF)
            {
                if(!pMLinKF->isBad())
                {
                    if(pMLinKF->Observations()>pML->Observations())
                        pML->Replace(pMLinKF);
                    else
                        pMLinKF->Replace(pML);
                }
            }
            else
            {
                pML->AddObservation(pKF,bestIdx);
                pKF->AddMapLine(pML,bestIdx);
            }
            nFused++;
        }
    }

    return nFused;
}

} //namespace ORB_SLAM

