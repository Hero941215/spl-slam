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
#ifndef INITIALIZER_H
#define INITIALIZER_H

#include<opencv2/opencv.hpp>
#include "Frame.h"


namespace PL_SLAM
{

// THIS IS THE INITIALIZER FOR MONOCULAR SLAM. NOT USED IN THE STEREO OR RGBD CASE.
class Initializer
{
    typedef pair<int,int> Match;

public:

    // Fix the reference frame
    //Initializer(const Frame &ReferenceFrame, float sigma = 1.0, int iterations = 200);
    Initializer(const Frame &ReferenceFrame, float sigma = 1.0, int iterations = 200, const bool usingLine = true);

    // Computes in parallel a fundamental matrix and a homography
    // Selects a model and tries to recover the motion and the structure from motion
    bool Initialize(const Frame &CurrentFrame, const vector<int> &vMatches12,
                    cv::Mat &R21, cv::Mat &t21, vector<cv::Point3f> &vP3D, vector<bool> &vbTriangulated);
    //点线特征初始化主函数：
    bool InitializeBoth(const Frame &CurrentFrame, const vector<int> &vMatches12, const vector<int> &vMatchesLines12, cv::Mat &R21, cv::Mat &t21,
		    vector<cv::Point3f> &vP3D, vector<cv::Point3f> &vMidP3D ,vector<cv::Point3f> &vFirP3D ,vector<cv::Point3f> &vEndP3D ,
		    vector<bool> &vbTriangulated, vector<bool> &vbTriangulatedMidP);

    //是否使用线特征
    bool mbusingLine;

private:

    void FindHomography(vector<bool> &vbMatchesInliers, float &score, cv::Mat &H21);
    void FindFundamental(vector<bool> &vbInliers, float &score, cv::Mat &F21);
    //基于匹配的点特征和线特征中点计算H和F矩阵，计算最高得分的主函数：
    void FindHomographyBoth(vector<bool> &vbMatchesInliers, vector<bool> &vbMatchesInlierLines , float &score, cv::Mat &H21);
    void FindFundamentalBoth(vector<bool> &vbInliers, vector<bool> &vbInlierLines,float &score, cv::Mat &F21);
    
    
    cv::Mat ComputeH21(const vector<cv::Point2f> &vP1, const vector<cv::Point2f> &vP2);
    cv::Mat ComputeF21(const vector<cv::Point2f> &vP1, const vector<cv::Point2f> &vP2);

    
    float CheckHomography(const cv::Mat &H21, const cv::Mat &H12, vector<bool> &vbMatchesInliers, float sigma);
    float CheckFundamental(const cv::Mat &F21, vector<bool> &vbMatchesInliers, float sigma);
    //对匹配的线特征中点，计算内线以及综合得分：
    float CheckHomographyLines(const cv::Mat &H21, const cv::Mat &H12, vector<bool> &vbMatchesInlierLines, float sigma);
    float CheckFundamentalLines(const cv::Mat &F21 , vector<bool> &vbMatchesInlierLines, float sigma);

    
    bool ReconstructF(vector<bool> &vbMatchesInliers, cv::Mat &F21, cv::Mat &K,
                      cv::Mat &R21, cv::Mat &t21, vector<cv::Point3f> &vP3D, vector<bool> &vbTriangulated, float minParallax, int minTriangulated);  
    bool ReconstructH(vector<bool> &vbMatchesInliers, cv::Mat &H21, cv::Mat &K,
                      cv::Mat &R21, cv::Mat &t21, vector<cv::Point3f> &vP3D, vector<bool> &vbTriangulated, float minParallax, int minTriangulated);
    //传入线特征的内线集，恢复出初步的3D线特征：
    void ReconstructFLines(vector<bool> &vbMatchesInlierLines,  cv::Mat &F21, cv::Mat &K, cv::Mat &R21, cv::Mat &t21, vector<cv::Point3f> &vMidP3D, 
			   vector<cv::Point3f> &vFirP3D, vector<cv::Point3f> &vEndP3D , vector<bool> &vbTriangulatedMidP ,float minParallax,
			   int minTriangulated, bool &RFLines);
    void ReconstructHLines(vector<bool> &vbMatchesInlierLines, cv::Mat &H21, cv::Mat &K, cv::Mat &R21, cv::Mat &t21, vector<cv::Point3f> &vMidP3D,
			   vector<cv::Point3f> &vFirP3D, vector<cv::Point3f> &vEndP3D, vector<bool> &vbTriangulatedMidP, float minParallax,
			   int minTriangulated, bool &RHLines);
    //传入点特征的内点集，恢复出初步的3D点特征：
    void ReconstructFPoints(vector<bool> &vbMatchesInliers, cv::Mat &F21, cv::Mat &K,cv::Mat &R21, cv::Mat &t21, vector<cv::Point3f> &vP3D, 
			    vector<bool> &vbTriangulated, float minParallax, int minTriangulated, bool &RFPoints);  
    void ReconstructHPoints(vector<bool> &vbMatchesInliers, cv::Mat &H21, cv::Mat &K,cv::Mat &R21, cv::Mat &t21, vector<cv::Point3f> &vP3D, 
			    vector<bool> &vbTriangulated, float minParallax, int minTriangulated, bool &RHPoints);
    
    
    void Triangulate(const cv::KeyPoint &kp1, const cv::KeyPoint &kp2, const cv::Mat &P1, const cv::Mat &P2, cv::Mat &x3D);
    void TriangulateLine(const KeyLine &kL1, const KeyLine &kL2, const cv::Mat &P1, const cv::Mat &P2, cv::Mat &xFirP3D, cv::Mat &xEndP3D);

    
    void Normalize(const vector<cv::KeyPoint> &vKeys, vector<cv::Point2f> &vNormalizedPoints, cv::Mat &T);
    void NormalizeBoth(const vector<cv::KeyPoint> &vKeys, const vector<cv::KeyPoint> &vKeysMidP, vector<cv::Point2f> &vNormalizedPoints,
		   vector<cv::Point2f> &vNormalizedMidPoints,  cv::Mat &T);
    
    //th2 = 5.991
    int CheckRT(const cv::Mat &R, const cv::Mat &t, const vector<cv::KeyPoint> &vKeys1, const vector<cv::KeyPoint> &vKeys2,
                       const vector<Match> &vMatches12, vector<bool> &vbInliers,
                       const cv::Mat &K, vector<cv::Point3f> &vP3D, float th2, vector<bool> &vbGood, float &parallax);
    
    //th2 = 3.841
    int CheckRTLines(const cv::Mat &R, const cv::Mat &t, const vector<cv::KeyPoint> &vKeys1, const vector<cv::KeyPoint> &vKeys2,
		 const vector<KeyLine> &vKeyLines1, const vector<KeyLine> &vKeyLines2, const vector<Match> &vMatches12, 
		 vector<bool> &vbInliers, const cv::Mat &K, vector<cv::Point3f> &vMidP3D, vector<cv::Point3f> &vFirP3D,
		 vector<cv::Point3f> &vEndP3D ,float th2, vector<bool> &vbGood, float &parallax);

    void DecomposeE(const cv::Mat &E, cv::Mat &R1, cv::Mat &R2, cv::Mat &t);


    // Keypoints from Reference Frame (Frame 1)
    vector<cv::KeyPoint> mvKeys1;
    // Keypoints from Current Frame (Frame 2)
    vector<cv::KeyPoint> mvKeys2;
    //线特征中点存储：
    vector<cv::KeyPoint> mvKeyMidP1; 
    vector<cv::KeyPoint> mvKeyMidP2; 
    //线特征关键线存储：
    std::vector<KeyLine> mvLinesUn1; 
    std::vector<KeyLine> mvLinesUn2; 

    // Current Matches from Reference to Current
    vector<Match> mvMatches12;
    vector<bool> mvbMatched1;
    //线特征的匹配记录：
    vector<Match> mvMatchesLines12; 
    vector<bool> mvbMatchedLines1; 

    // Calibration
    cv::Mat mK;

    // Standard Deviation and Variance
    float mSigma, mSigma2;

    // Ransac max iterations
    int mMaxIterations;

    // Ransac sets
    vector<vector<size_t> > mvSets;   

};

} //namespace ORB_SLAM

#endif // INITIALIZER_H
