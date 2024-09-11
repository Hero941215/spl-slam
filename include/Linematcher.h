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

#ifndef LINEMATCHER_H
#define LINEMATCHER_H

#include<vector>
#include<opencv2/core/core.hpp>
#include<opencv2/features2d/features2d.hpp>

#include"MapLine.h"
#include"KeyFrame.h"
#include"Frame.h"

namespace PL_SLAM
{
  
 class Linematcher 
{
public:  
  
    //构造函数：
    Linematcher(float nnratio=0.6, bool checkOri=true, bool checklen= true, float lengtherr= 0.1);

    //线特征描述子距离：
    static int DescriptorDistance(const cv::Mat &a, const cv::Mat &b); 
  
    //a.前端
    int SearchForInitialization(Frame &F1, Frame &F2, std::vector<cv::Point2f> &vbPrevMatched, std::vector<int> &vnMatches12, int windowSize=10);
    int SearchByProjection(Frame &CurrentFrame, const Frame &LastFrame, const float th, const bool bMono );
    int SearchByKNN(KeyFrame *pKF, Frame &F, std::vector<MapLine*> &vpMapLineMatches );
    int SearchByProjection(Frame &F, const std::vector<MapLine*> &vpMapLines, const float th=3);
    int SearchByProjection(Frame &CurrentFrame, KeyFrame* pKF, const std::set<MapLine*> &sAlreadyFound, const float th, const int Linedist);

    //b.后端
    int SearchForTriangulation(KeyFrame *pKF1, KeyFrame* pKF2, cv::Mat F12, std::vector<pair<size_t, size_t> > &vMatchedPairs);
    int Fuse(KeyFrame* pKF, const vector<MapLine *> &vpMapLines, const float th=8.0);
    
public:  
  
    //线特征描述子的距离阈值：
    static const int TH_HIGH;
    static const int TH_LOW;
    static const int HISTO_LENGTH;
    
protected:
  
    //检查核线约束：
    bool CheckDistEpipolarLine(const cv::KeyPoint &kp1, const cv::KeyPoint &kp2, const cv::Mat &F12, const KeyFrame *pKF);
    void ComputeThreeMaxima(std::vector<int>* histo, const int L, int &ind1, int &ind2, int &ind3);

    //SearchForTriangulation()和 SearchByKNN()调用的KNN匹配子函数，内部除了匹配，还完成了最优和次优描述子距离比较：
    void matchNNR(const cv::Mat &desc1, const cv::Mat &desc2, float nnr, std::vector<int> &matches_12,int &nmatches);

    //通过视角计算初步的投影匹配搜索半径：
    float RadiusByViewingCos(const float &viewCos);
  
    float mfNNratio;
    bool mbCheckOrientation;
    bool mbchecklen;
    float mflengtherr;
    
};

}// namespace PL_SLAM


#endif //LINEMATCHER_H

