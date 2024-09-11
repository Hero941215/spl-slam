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

#ifndef LINEEXTRACTOR_H
#define LINEEXTRACTOR_H

#include <vector>
#include <opencv/cv.h>

#include <opencv2/core.hpp>

#include "Thirdparty/line_descriptor/include/line_descriptor_custom.hpp"
#include "Thirdparty/line_descriptor/include/line_descriptor/descriptor_custom.hpp"

#include "Precomp.h" 

using namespace cv;
using namespace line_descriptor;

namespace PL_SLAM
{
  
//FLD关键线提取结构体
struct SEGMENT
{
    float x1, y1, x2, y2, angle;
};
  
class Lineextractor
{	
public:
    
    //LSD-LBD构造函数
    Lineextractor(int nfeatures=240, int nlevels=3, int refine=0, double scale=1.05, double sigma_scale=0.6, 
		  double quant=2.0,double ang_th=22.5, double log_eps=1.0, double density_th=0.7, int n_bins=1024,
		  double min_line_length = 32.0, bool busingLSD = true); 
    
    
    //FLD-LBD构造函数
    Lineextractor(int _nfeatures = 240, int _nlevels = 1, double _scale = 1.05, int _length_threshold = 10, 
		  float _distance_threshold = 1.414213562f,double _canny_th1 = 50.0, double _canny_th2 = 100.0,
		  int _canny_aperture_size = 3, bool _do_merge = false, bool _busingLSD = false);
    
    ~Lineextractor(){}

    void ComputeLsdWithLbd(const cv::Mat &image ,std::vector<KeyLine>& keyLines, std::vector<cv::KeyPoint> &keypoints, cv::Mat &descriptors);

    
    //基于关键线响应对关键线进行比较
    struct sort_lines_by_response
    {
	inline bool operator()(const KeyLine& a, const KeyLine& b){
	    return ( a.response > b.response );
	}
    };
    
    //用于计算图像上的FLD关键线、关键线中点位置和LBD描述子的主函数
    void ComputeFldWithLbdOld(cv::Mat &image ,std::vector<KeyLine>& keyLines, 
			   std::vector<cv::KeyPoint> &keypoints, cv::Mat &descriptors);

    void ComputeFldWithLbd(cv::Mat &image ,std::vector<KeyLine>& keyLines, 
			   std::vector<cv::KeyPoint> &keypoints, cv::Mat &descriptors);
    
    //单层FLD关键线提取
    void detect(InputArray _image, OutputArray _lines);
    
    //扩展为基于图像金字塔的FLD关键线提取
    void detectFldWithPyramid(std::vector<std::vector<cv::Vec4f> > &lines_fld);
    
    struct sort_flines_by_length
    {
	inline bool operator()(const Vec4f& a, const Vec4f& b){
	    return ( sqrt(pow(a(0)-a(2),2.0)+pow(a(1)-a(3),2.0)) > sqrt(pow(b(0)-b(2),2.0)+pow(b(1)-b(3),2.0)) );
	}
    };

    int inline GetLevels(){
        return nlevels;
    }

    float inline GetScaleFactor(){
        return scale;      
    }

    std::vector<float> inline GetScaleFactors(){
        return mvScaleFactor;
    }

    std::vector<float> inline GetInverseScaleFactors(){
	return mvInvScaleFactor;
    } 

    std::vector<float> inline GetScaleSigmaSquares(){
        return mvLevelSigma2;
    }

    std::vector<float> inline GetInverseScaleSigmaSquares(){
        return mvInvLevelSigma2;
    }
    
protected:
  
    int nfeatures;
    int nlevels;

    //LSD线特征检测参数：
    int refine;//默认0
    double scale; //尺度金字塔参数，尺度因子，默认1.05
    double sigma_scale; //下采样，默认0.6
    double quant;//像素误差，默认2
    double ang_th;//容许线梯度误差，默认22.5度
    double log_eps;//最终NFA容错检测阈值，默认1.0
    double density_th;//线特征矩阵有效像素比例，默认0.7
    int n_bins;//所有像素根据梯度分级排序，默认1024级别
    double min_line_length;//最小线段长度
    
    //尺度金字塔具体参数：
    std::vector<float> mvScaleFactor;
    std::vector<float> mvInvScaleFactor;    
    std::vector<float> mvLevelSigma2;
    std::vector<float> mvInvLevelSigma2;

    //预设各层图像保留线特征数
    std::vector<int> mnFeaturesPerLevel;
     
    //FLD线特征检测参数：
    int imagewidth, imageheight, threshold_length;
    float threshold_dist;
    double canny_th1, canny_th2;
    int canny_aperture_size;
    bool do_merge;
    
    //Lineextractor& operator= (const Lineextractor&); // to quiet MSVC
    
    template<class T>
    void incidentPoint(const Mat& l, T& pt);

    void mergeLines(const SEGMENT& seg1, const SEGMENT& seg2, SEGMENT& seg_merged);

    bool mergeSegments(const SEGMENT& seg1, const SEGMENT& seg2, SEGMENT& seg_merged);

    bool getPointChain(const Mat& img, Point pt, Point& chained_pt, float& direction, int step);

    double distPointLine(const Mat& p, Mat& l);

    void extractSegments(const std::vector<Point2i>& points, std::vector<SEGMENT>& segments );

    void lineDetection(const Mat& src, std::vector<SEGMENT>& segments_all);

    void pointInboardTest(const Mat& src, Point2i& pt);

    inline void getAngle(SEGMENT& seg);

    void additionalOperationsOnSegment(const Mat& src, SEGMENT& seg);

    void drawSegment(Mat& mat, const SEGMENT& seg, Scalar bgr = Scalar(0,255,0),
	    int thickness = 1, bool directed = true);
    
    void ComputePyramid(cv::Mat image);
    
public:
    //LSD-FLD辨识Bool量
    bool busingLSD;
    
    //FLD尺度金字塔图像存储向量
    std::vector<cv::Mat> mvImagePyramid;
  
};

  
}//namespace PL_SLAM

#endif