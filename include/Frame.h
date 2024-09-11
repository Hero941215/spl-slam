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

#ifndef FRAME_H
#define FRAME_H

#include<vector>

#include "MapPoint.h"
#include "MapLine.h"
#include "Thirdparty/DBoW2/DBoW2/BowVector.h"
#include "Thirdparty/DBoW2/DBoW2/FeatureVector.h"
#include "ORBVocabulary.h"
#include "KeyFrame.h"
#include "ORBextractor.h"
#include "Lineextractor.h"

#include <opencv2/opencv.hpp>

namespace PL_SLAM
{
#define FRAME_GRID_ROWS 48
#define FRAME_GRID_COLS 64

//640*480 -- 24/32
//Kitti --12/16
  
#define FRAME_GRID_ROWS_Lines 12
#define FRAME_GRID_COLS_Lines 16
  
  
class MapPoint;
class MapLine; 
class KeyFrame;

class Frame
{
public:
  
    Frame();

    // Copy constructor.
    Frame(const Frame &frame);

    // Constructor for stereo cameras.
    Frame(const cv::Mat &imLeft, const cv::Mat &imRight, const double &timeStamp, ORBextractor* extractorLeft, ORBextractor* extractorRight, ORBVocabulary* voc, cv::Mat &K, cv::Mat &distCoef, const float &bf, const float &thDepth);

    // Constructor for RGB-D cameras.
    Frame(const cv::Mat &imGray, const cv::Mat &imDepth, const double &timeStamp, ORBextractor* extractor,ORBVocabulary* voc, cv::Mat &K, cv::Mat &distCoef, const float &bf, const float &thDepth);

    // Constructor for Monocular cameras.
    Frame(const cv::Mat &imGray, const double &timeStamp, ORBextractor* extractor,ORBVocabulary* voc, cv::Mat &K, cv::Mat &distCoef, const float &bf, const float &thDepth);
 
    //点线特征同时使用的单目图像帧(LineExpanding)
    Frame(const cv::Mat &imGray, const double &timeStamp, ORBextractor* extractor, ORBVocabulary* voc, Lineextractor* lineextractor ,cv::Mat &K, cv::Mat &distCoef, const float &bf, const float &thDepth);
    
    // Extract ORB on the image. 0 for left image and 1 for right image.
    void ExtractORB(int flag, const cv::Mat &im);

    void ExtractLsdWithLBD(const cv::Mat &im);
    
    void ExtractFldWithLBD(cv::Mat &im);
    
    // Compute Bag of Words representation.
    void ComputeBoW();

    // Set the camera pose.(Point+LineExpanding)
    void SetPose(cv::Mat Tcw);
    void SetPosePoints(cv::Mat TcwPoints);
    void SetPoseLines(cv::Mat TcwLines);
    
    // Computes rotation, translation and camera center matrices from the camera pose.
    void UpdatePoseMatrices();

    // Returns the camera center.
    inline cv::Mat GetCameraCenter(){
        return mOw.clone();
    }

    // Returns inverse of rotation
    inline cv::Mat GetRotationInverse(){
        return mRwc.clone();
    }

    // Check if a MapPoint is in the frustum of the camera
    // and fill variables of the MapPoint to be used by the tracking
    bool isInFrustum(MapPoint* pMP, float viewingCosLimit);   
    bool isInFrustumLine(MapLine* pML, float viewingCosLimit);
    
    // Compute the cell of a keypoint (return false if outside the grid)
    bool PosInGrid(const cv::KeyPoint &kp, int &posX, int &posY);
    bool PosInGridLines(const KeyLine &kL, const cv::KeyPoint &kp, int &posMidX1, int &posMidY1);

    vector<size_t> GetFeaturesInArea(const float &x, const float  &y, const float  &r, const int minLevel=-1, const int maxLevel=-1) const;
    vector<size_t> GetFeaturesInAreaLines(const float &x, const float  &y, const float  &r, const int minLevel=-1, const int maxLevel=-1) const;
  
    // Search a match for each keypoint in the left image to a keypoint in the right image.
    // If there is a match, depth is computed and the right coordinate associated to the left keypoint is stored.
    void ComputeStereoMatches();

    // Associate a "right" coordinate to a keypoint if there is valid depth in the depthmap.
    void ComputeStereoFromRGBD(const cv::Mat &imDepth);

    // Backprojects a keypoint (if stereo/depth info available) into 3D world coordinates.
    cv::Mat UnprojectStereo(const int &i);

public:
    // Vocabulary used for relocalization.
    ORBVocabulary* mpORBvocabulary;

    // Feature extractor. The right is used only in the stereo case.
    ORBextractor* mpORBextractorLeft, *mpORBextractorRight;
    Lineextractor* mpLineextractor;
    
    // Frame timestamp.
    double mTimeStamp;

    // Calibration matrix and OpenCV distortion parameters.
    cv::Mat mK;
    static float fx;
    static float fy;
    static float cx;
    static float cy;
    static float invfx;
    static float invfy;
    cv::Mat mDistCoef;

    // Stereo baseline multiplied by fx.
    float mbf;

    // Stereo baseline in meters.
    float mb;

    // Threshold close/far points. Close points are inserted from 1 view.
    // Far points are inserted as in the monocular case from 2 views.
    float mThDepth;

    // Number of KeyPoints.
    int N;
    int NL;
    
    
    // Vector of keypoints (original for visualization) and undistorted (actually used by the system).
    // In the stereo case, mvKeysUn is redundant as images must be rectified.
    // In the RGB-D case, RGB images can be distorted.
    std::vector<cv::KeyPoint> mvKeys, mvKeysRight;
    std::vector<cv::KeyPoint> mvKeysUn;
    // Corresponding stereo coordinate and depth for each keypoint.
    // "Monocular" keypoints have a negative value.
    std::vector<float> mvuRight;
    std::vector<float> mvDepth;

    
    std::vector<KeyLine> mvLines; 
    std::vector<cv::KeyPoint> mvMidPoints;
    std::vector<KeyLine> mvLinesUn; 
    std::vector<cv::KeyPoint>  mvMidPointsUn;
   
    
    // Bag of Words Vector structures.
    DBoW2::BowVector mBowVec;
    DBoW2::FeatureVector mFeatVec;

    
    // ORB descriptor, each row associated to a keypoint.
    cv::Mat mDescriptors, mDescriptorsRight;   
    //线特征描述子：(LineExpanding)
    cv::Mat mDescriptorLines;

    // MapPoints associated to keypoints, NULL pointer if no association.
    std::vector<MapPoint*> mvpMapPoints;
    std::vector<MapLine*> mvpMapLines;//内含地图线中点
    
    
    // Flag to identify outlier associations.
    std::vector<bool> mvbOutlier;
    //被检测出的外线：(LineExpanding)
    std::vector<bool> mvbOutlierLines; 
    
    // Keypoints are assigned to cells in a grid to reduce matching complexity when projecting MapPoints.
    static float mfGridElementWidthInv;
    static float mfGridElementHeightInv;
    std::vector<std::size_t> mGrid[FRAME_GRID_COLS][FRAME_GRID_ROWS];

    static float mfGridElementWidthInvLines;
    static float mfGridElementHeightInvLines;
    std::vector<std::size_t> mGridLines[FRAME_GRID_COLS_Lines][FRAME_GRID_ROWS_Lines];
    
    // Camera pose.
    cv::Mat mTcw;
    cv::Mat mTcwPoints;
    cv::Mat mTcwLines;
    
    // Current and Next Frame id.
    static long unsigned int nNextId;
    long unsigned int mnId;

    // Reference Keyframe.
    KeyFrame* mpReferenceKF;//该关键帧为与当前图像帧有最多匹配点的关键帧
    KeyFrame* mpReferenceKFLines;

    // Scale pyramid info.
    int mnScaleLevels;
    float mfScaleFactor;
    float mfLogScaleFactor;
    vector<float> mvScaleFactors;
    vector<float> mvInvScaleFactors;
    vector<float> mvLevelSigma2;
    vector<float> mvInvLevelSigma2;
    //线特征提取的尺度金字塔信息：(LineExpanding)
    int mnScaleLevelLines;
    float mfScaleFactorLines;
    float mfLogScaleFactorsLines;
    vector<float> mvScaleFactorsLines;
    vector<float> mvInvScaleFactorsLines;
    vector<float> mvLevelSigma2Lines;
    vector<float> mvInvLevelSigma2Lines;
    
    
    // Undistorted Image Bounds (computed once).
    static float mnMinX;
    static float mnMaxX;
    static float mnMinY;
    static float mnMaxY;

    static bool mbInitialComputations;


private:

    // Undistort keypoints given OpenCV distortion parameters.
    // Only for the RGB-D case. Stereo must be already rectified!
    // (called in the constructor).
    void UndistortKeyPoints();
    void UndistortKeyLines();
    
    //计算校正后关键线的2D长度向量，重新存入KeyLine.lineLength中：(LineExpanding)
    void ComputeLines2DLengthUn();
    
    // Computes image bounds for the undistorted image (called in the constructor).
    void ComputeImageBounds(const cv::Mat &imLeft);

    // Assign keypoints to the grid for speed up feature matching (called in the constructor).
    void AssignFeaturesToGrid();
    void AssignFeaturesToGridLines();

    // Rotation, translation and camera center
    cv::Mat mRcw;
    cv::Mat mtcw;
    cv::Mat mRwc;
    cv::Mat mOw; //==mtwc
};

}// namespace ORB_SLAM

#endif // FRAME_H
