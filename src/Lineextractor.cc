/**
*  This file is part of SPL-SLAM.
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
* 使用了OPENCV中的LSD关键线提取和LBD线描述子计算
*
* 使用了OPENCV扩展包中的FLD源码
*/

#include "Lineextractor.h"
#include <vector>

using namespace std;
using namespace cv;

namespace PL_SLAM
{

const int EDGE_THRESHOLD = 19; 
    
//LSD-LBD线特征提取构造函数
Lineextractor::Lineextractor(int _nfeatures, int _nlevels, int _refine, double _scale, double _sigma_scale, double _quant,double _ang_th, double _log_eps, double _density_th, int _n_bins,double _min_line_length, bool _busingLSD ):          
	nfeatures(_nfeatures), nlevels(_nlevels), refine(_refine), scale(_scale), sigma_scale(_sigma_scale), quant(_quant),
	ang_th(_ang_th), log_eps(_log_eps), density_th(_density_th), n_bins(_n_bins), min_line_length(_min_line_length), busingLSD(_busingLSD)
{
	mvScaleFactor.resize(nlevels);
	mvLevelSigma2.resize(nlevels);
	mvScaleFactor[0]=1.0f;
	mvLevelSigma2[0]=1.0f;
	for(int i=1; i<nlevels; i++)
	{
	    mvScaleFactor[i]=mvScaleFactor[i-1]*scale;
	    mvLevelSigma2[i]=mvScaleFactor[i]*mvScaleFactor[i];
	}

	mvInvScaleFactor.resize(nlevels);
	mvInvLevelSigma2.resize(nlevels);
	for(int i=0; i<nlevels; i++)
	{
	    mvInvScaleFactor[i]=1.0f/mvScaleFactor[i];
	    mvInvLevelSigma2[i]=1.0f/mvLevelSigma2[i];
	}

	mnFeaturesPerLevel.resize(nlevels);

	float factor = 1.0f / scale;
	float nDesiredFeaturesPerScale = nfeatures*(1 - factor)/(1 - (float)pow((double)factor, (double)nlevels));

	int sumFeatures = 0;
	for(int level = 0; level < nlevels-1; level++ )
	{
	    mnFeaturesPerLevel[level] = cvRound(nDesiredFeaturesPerScale);
	    sumFeatures += mnFeaturesPerLevel[level];
	    nDesiredFeaturesPerScale *= factor;
	}
	mnFeaturesPerLevel[nlevels-1] = std::max(nfeatures - sumFeatures, 0);
}

//FLD-LBD线特征提取构造函数
Lineextractor::Lineextractor(int _nfeatures, int _nlevels, double _scale, int _length_threshold, float _distance_threshold,
        double _canny_th1, double _canny_th2, int _canny_aperture_size, bool _do_merge, bool _busingLSD)
    :nfeatures(_nfeatures), nlevels(_nlevels), scale(_scale), threshold_length(_length_threshold), threshold_dist(_distance_threshold),
    canny_th1(_canny_th1), canny_th2(_canny_th2), canny_aperture_size(_canny_aperture_size), do_merge(_do_merge), busingLSD(_busingLSD)
{
	CV_Assert(_length_threshold > 0 && _distance_threshold > 0 &&
		_canny_th1 > 0 && _canny_th2 > 0 && _canny_aperture_size > 0);
	
	mvScaleFactor.resize(nlevels);
	mvLevelSigma2.resize(nlevels);
	mvScaleFactor[0]=1.0f;
	mvLevelSigma2[0]=1.0f;
	for(int i=1; i<nlevels; i++)
	{
	    mvScaleFactor[i]=mvScaleFactor[i-1]*scale;
	    mvLevelSigma2[i]=mvScaleFactor[i]*mvScaleFactor[i];
	}

	mvInvScaleFactor.resize(nlevels);
	mvInvLevelSigma2.resize(nlevels);
	for(int i=0; i<nlevels; i++)
	{
	    mvInvScaleFactor[i]=1.0f/mvScaleFactor[i];
	    mvInvLevelSigma2[i]=1.0f/mvLevelSigma2[i];
	}

	mnFeaturesPerLevel.resize(nlevels);

	float factor = 1.0f / scale;
	float nDesiredFeaturesPerScale = nfeatures*(1 - factor)/(1 - (float)pow((double)factor, (double)nlevels));

	int sumFeatures = 0;
	for(int level = 0; level < nlevels-1; level++ )
	{
	    mnFeaturesPerLevel[level] = cvRound(nDesiredFeaturesPerScale);
	    sumFeatures += mnFeaturesPerLevel[level];
	    nDesiredFeaturesPerScale *= factor;
	}
	mnFeaturesPerLevel[nlevels-1] = std::max(nfeatures - sumFeatures, 0);
      
}
  
void Lineextractor::ComputeLsdWithLbd(const cv::Mat &image ,std::vector<KeyLine>& keyLines, 
			     std::vector<cv::KeyPoint>& keypoints, cv::Mat &descriptors) 
{
	if(image.empty())
        return;	
	
	//创建LBD线特征描述子动态指针
	Ptr<BinaryDescriptor>  lbd = BinaryDescriptor::createBinaryDescriptor();
	
	//创建LSD线特征提取器动态指针
	Ptr<line_descriptor::LSDDetectorC> lsd = line_descriptor::LSDDetectorC::createLSDDetectorC();
	
	//根据线特征类中参数，定义LSD线特征提取器提取参数
	line_descriptor::LSDDetectorC::LSDOptions opts;
	opts.refine       = refine;
        opts.scale        = scale;
        opts.sigma_scale  = sigma_scale;
        opts.quant        = quant;
        opts.ang_th       = ang_th;
        opts.log_eps      = log_eps;
        opts.density_th   = density_th;
        opts.n_bins       = n_bins;
        opts.min_length   = min_line_length;
	
	lsd->detect( image, keyLines, 2, nlevels, opts);
		
	//step1: 重新获取每一层关键线
	int octaveIdx = 0;
	vector<KeyLine> vDetectKeyLineByOctave;
	vector<vector<KeyLine> > vvDetectKeyLineByOctaveAndResponse;
	//将当前层数的所有关键线存入对应的vDetectKeyLineByOctaveAndResponse中
	for(vector<KeyLine>::iterator keyline=keyLines.begin(),keylineEnd=keyLines.end(); keyline!=keylineEnd; keyline++)
	{
	    
	    if(keyline->octave==octaveIdx)
	    {
		vDetectKeyLineByOctave.push_back(*keyline);
	    }
	    else
	    {
		vvDetectKeyLineByOctaveAndResponse.push_back(vDetectKeyLineByOctave);
		vDetectKeyLineByOctave.clear();
		vDetectKeyLineByOctave.push_back(*keyline);
		octaveIdx++;
	    }
	}
	//存入最后一层所有关键线
	vvDetectKeyLineByOctaveAndResponse.push_back(vDetectKeyLineByOctave);
	
	//step2: 将每一层关键线根据Response进行保存，最多为mnFeaturesPerLevel中对应层的数目
	for(size_t octaveIdy=0; octaveIdy<vvDetectKeyLineByOctaveAndResponse.size(); octaveIdy++)
	{
	    //获取当前层的预设关键线提取数
	    size_t mnFeaturesPerLevelNow = mnFeaturesPerLevel[octaveIdy];
	    
	    //判断当前层的关键线数是否小于预设需要提取的数目
	    if(vvDetectKeyLineByOctaveAndResponse[octaveIdy].size()<=mnFeaturesPerLevelNow)
	    {
		continue;//无需处理
	    }
	    else
	    {
		//基于关键线响应对关键线进行比较
		sort( vvDetectKeyLineByOctaveAndResponse[octaveIdy].begin(), 
		    vvDetectKeyLineByOctaveAndResponse[octaveIdy].end(), sort_lines_by_response() );
		
		//重新设置关键线数量
		vvDetectKeyLineByOctaveAndResponse[octaveIdy].resize(mnFeaturesPerLevel[octaveIdy]);
	    }    
	}
	
	//清除原关键线向量
	keyLines.clear();	
	
	//step3: 将每一层关键线重新放入输入参数引用keyLines，同时放入输入参数引用keypoints
	for(size_t octaveIdz=0; octaveIdz<vvDetectKeyLineByOctaveAndResponse.size(); octaveIdz++)
	{
	    for(size_t k=0; k<vvDetectKeyLineByOctaveAndResponse[octaveIdz].size(); k++)
	    {
		//存储KeyLines 
		KeyLine k1 = vvDetectKeyLineByOctaveAndResponse[octaveIdz][k];
		keyLines.push_back( k1 );
		//存储keypoints
		KeyPoint kp;
		kp.pt.x = ( k1.startPointX + k1.endPointX ) / 2;
		kp.pt.y = ( k1.startPointY + k1.endPointY ) / 2;
		kp.octave = k1.octave;
		keypoints.push_back( kp );
	    }
	}	
	
	//重新设置对应关键线的ID
	for(size_t k=0; k<keyLines.size(); k++)
	{
	  keyLines[k].class_id = k;
	}
	  
	//step4: 计算关键线的LBD描述子
	lbd->compute( image, keyLines, descriptors);
	
}

//后验检测FLD关键线是否超出图像边界
inline void checkLineExtremes( cv::Vec4f& extremes, cv::Size imageSize )
{
      if( extremes[0] < 0 )
	extremes[0] = 0;

      if( extremes[0] >= imageSize.width )
	extremes[0] = (float)imageSize.width - 1.0f;

      if( extremes[2] < 0 )
	extremes[2] = 0;

      if( extremes[2] >= imageSize.width )
	extremes[2] = (float)imageSize.width - 1.0f;

      if( extremes[1] < 0 )
	extremes[1] = 0;

      if( extremes[1] >= imageSize.height )
	extremes[1] = (float)imageSize.height - 1.0f;

      if( extremes[3] < 0 )
	extremes[3] = 0;

      if( extremes[3] >= imageSize.height )
	extremes[3] = (float)imageSize.height - 1.0f;
}

void Lineextractor::ComputeFldWithLbd(cv::Mat &image ,std::vector<KeyLine>& keyLines, 
			   std::vector<cv::KeyPoint> &keypoints, cv::Mat &descriptors)
{
      if(image.empty())
	  return;
      
      //创建LBD线特征描述子动态指针
      Ptr<BinaryDescriptor>  lbd = BinaryDescriptor::createBinaryDescriptor();
	  
      //确保满足FLD提取器的图像输入格式要求
      cv::Mat fld_img;    
      image.convertTo( fld_img, CV_8UC1);
      
      //计算尺度金字塔
      ComputePyramid(fld_img);
      
      //存储不同层的FLD关键线
      std::vector<std::vector<cv::Vec4f> > vvlines_fld;
      
      //提取不同层的FLD关键线
      detectFldWithPyramid(vvlines_fld);
      
      //step2: 将每一层关键线根据Length进行保存，最多为mnFeaturesPerLevel中对应层的数目
      for(size_t octaveIdx=0; octaveIdx<vvlines_fld.size(); octaveIdx++)
      {
	  //获取当前层的预设关键线提取数
	  size_t mnFeaturesPerLevelNow = mnFeaturesPerLevel[octaveIdx];
	  
	  //判断当前层的关键线数是否小于预设需要提取的数目
	  if(vvlines_fld[octaveIdx].size()<=mnFeaturesPerLevelNow)
	  {
	      continue;
	  }
	  else
	  {
	      //基于关键线长度对关键线进行比较
	      sort( vvlines_fld[octaveIdx].begin(), 
		  vvlines_fld[octaveIdx].end(), sort_flines_by_length());
	      
	      //重新设置关键线数量
	      vvlines_fld[octaveIdx].resize(mnFeaturesPerLevel[octaveIdx]);
	  }    
      }
      
      //清除原关键线向量
      keyLines.clear();
      int class_counter = -1;
      for(size_t octaveIdy=0; octaveIdy<vvlines_fld.size(); octaveIdy++)
      {
	  //获取当层尺度
	  float octaveScale = pow( (float)scale, octaveIdy );
	  
	  for(size_t k=0; k<vvlines_fld[octaveIdy].size(); k++)
	  {
	      //存储KeyLines 
	      KeyLine kl;	    
	      cv::Vec4f extremes = vvlines_fld[octaveIdy][k];
	      
	      //保证提取关键线不超过图像边界
	      checkLineExtremes( extremes, mvImagePyramid[octaveIdy].size() );
	      
	      kl.startPointX = extremes[0] * octaveScale;
	      kl.startPointY = extremes[1] * octaveScale;
	      kl.endPointX = extremes[2] * octaveScale;
	      kl.endPointY = extremes[3] * octaveScale;
	      kl.sPointInOctaveX = extremes[0];
	      kl.sPointInOctaveY = extremes[1];
	      kl.ePointInOctaveX = extremes[2];
	      kl.ePointInOctaveY = extremes[3];
	      
	      kl.lineLength = (float) sqrt( pow( extremes[0] - extremes[2], 2 ) + pow( extremes[1] - extremes[3], 2 ) );
	      
	      /* compute number of pixels covered by line */
	      LineIterator li( mvImagePyramid[octaveIdy], Point2f( extremes[0], extremes[1] ), Point2f( extremes[2], extremes[3] ) );
	      kl.numOfPixels = li.count;

	      kl.angle = atan2( ( kl.endPointY - kl.startPointY ), ( kl.endPointX - kl.startPointX ) );
	      kl.class_id = ++class_counter;
	      kl.octave = octaveIdy;
	      kl.size = ( kl.endPointX - kl.startPointX ) * ( kl.endPointY - kl.startPointY );
	      kl.response = kl.lineLength / max( mvImagePyramid[octaveIdy].cols, mvImagePyramid[octaveIdy].rows );
	      
	      keyLines.push_back( kl );
	      //存储keypoints
	      KeyPoint kp;
	      kp.pt.x = ( kl.startPointX + kl.endPointX ) / 2;
	      kp.pt.y = ( kl.startPointY + kl.endPointY ) / 2;
	      kp.octave = kl.octave;
	      keypoints.push_back( kp );
	  }
      }         
      
      lbd->compute( image, keyLines, descriptors);
    
}

void Lineextractor::ComputeFldWithLbdOld(cv::Mat &image ,std::vector<KeyLine>& keyLines, 
			   std::vector<cv::KeyPoint> &keypoints, cv::Mat &descriptors)
{
    if(image.empty())
        return;
    
    //创建LBD线特征描述子动态指针
    Ptr<BinaryDescriptor>  lbd = BinaryDescriptor::createBinaryDescriptor();
	
    //确保满足FLD提取器的图像输入格式要求
    cv::Mat fld_img;    
    image.convertTo( fld_img, CV_8UC1);
    
    vector<Vec4f> fld_lines;
       
    detect( fld_img, fld_lines );
    
    size_t featurenum = nfeatures;
    
     // filter lines
    if( fld_lines.size()>featurenum && featurenum!=0  )
    {
	// sort lines by their length
	sort( fld_lines.begin(), fld_lines.end(), sort_flines_by_length() );
	fld_lines.resize(featurenum);
    }
    
    keyLines.clear();
    keypoints.clear();
    
    keyLines.reserve(fld_lines.size());
    keypoints.reserve(fld_lines.size());
    
    for( size_t i = 0; i < fld_lines.size(); i++ )
    {
	KeyLine kl;
	float octaveScale = 1.f;
	int    octaveIdx   = 0;

	kl.startPointX     = fld_lines[i][0] * octaveScale;
	kl.startPointY     = fld_lines[i][1] * octaveScale;
	kl.endPointX       = fld_lines[i][2] * octaveScale;
	kl.endPointY       = fld_lines[i][3] * octaveScale;
	
	kl.sPointInOctaveX = fld_lines[i][0];
	kl.sPointInOctaveY = fld_lines[i][1];
	kl.ePointInOctaveX = fld_lines[i][2];
	kl.ePointInOctaveY = fld_lines[i][3];

	kl.lineLength = (float) sqrt( pow( fld_lines[i][0] - fld_lines[i][2], 2 ) + pow( fld_lines[i][1] - fld_lines[i][3], 2 ) );

	kl.angle    = atan2( ( kl.endPointY - kl.startPointY ), ( kl.endPointX - kl.startPointX ) );
	kl.class_id = i;
	kl.octave   = octaveIdx;
	kl.size     = ( kl.endPointX - kl.startPointX ) * ( kl.endPointY - kl.startPointY );
	kl.pt       = Point2f( ( kl.endPointX + kl.startPointX ) / 2, ( kl.endPointY + kl.startPointY ) / 2 );

	kl.response = kl.lineLength / max( fld_img.cols, fld_img.rows );
	cv::LineIterator li( fld_img, Point2f( fld_lines[i][0], fld_lines[i][1] ), Point2f( fld_lines[i][2], fld_lines[i][3] ) );
	kl.numOfPixels = li.count;

	keyLines.push_back( kl );
	
	cv::KeyPoint kp;
	kp.pt.x = ( kl.startPointX + kl.endPointX ) / 2;
	kp.pt.y = ( kl.startPointY + kl.endPointY ) / 2;	
	kp.octave = kl.octave;
	keypoints.push_back( kp );
    }
      
    // compute lbd descriptor
    lbd->compute( fld_img, keyLines, descriptors);
    
}

void Lineextractor::ComputePyramid(cv::Mat image)
{
    /* clear class fields */
  mvImagePyramid.clear();

  /* insert input image into pyramid */
  cv::Mat currentMat = image.clone();
  //cv::GaussianBlur( currentMat, currentMat, cv::Size( 5, 5 ), 1 );
  mvImagePyramid.push_back( currentMat );

  /* fill Gaussian pyramid */
  for ( int pyrCounter = 1; pyrCounter < nlevels; pyrCounter++ )
  {
    /* compute and store next image in pyramid and its size */
    pyrDown( currentMat, currentMat, Size( currentMat.cols / 2, currentMat.rows / 2 ) );
    mvImagePyramid.push_back( currentMat );
  }   
}

//基于尺度金字塔的线特征提取函数
void Lineextractor::detectFldWithPyramid(std::vector<std::vector<cv::Vec4f> > &lines_fld)
{
    for ( int i = 0; i < nlevels; i++ )
    {
      std::vector<Vec4f> octave_lines;
      detect( mvImagePyramid[i], octave_lines );
      lines_fld.push_back( octave_lines );
    }
}

void Lineextractor::detect(InputArray _image, OutputArray _lines)
{
    //CV_INSTRUMENT_REGION();

    Mat image = _image.getMat();
    CV_Assert(!image.empty() && image.type() == CV_8UC1);

    std::vector<Vec4f> lines;
    std::vector<SEGMENT> segments;
    lineDetection(image, segments);
    for(size_t i = 0; i < segments.size(); ++i)
    {
        const SEGMENT seg = segments[i];
        Vec4f line(seg.x1, seg.y1, seg.x2, seg.y2);
        lines.push_back(line);
    }
    Mat(lines).copyTo(_lines);
}

void Lineextractor::mergeLines(const SEGMENT& seg1, const SEGMENT& seg2, SEGMENT& seg_merged)
{
    double xg = 0.0, yg = 0.0;
    double delta1x = 0.0, delta1y = 0.0, delta2x = 0.0, delta2y = 0.0;
    float ax = 0, bx = 0, cx = 0, dx = 0;
    float ay = 0, by = 0, cy = 0, dy = 0;
    double li = 0.0, lj = 0.0;
    double thi = 0.0, thj = 0.0, thr = 0.0;
    double axg = 0.0, bxg = 0.0, cxg = 0.0, dxg = 0.0, delta1xg = 0.0, delta2xg = 0.0;

    ax = seg1.x1;
    ay = seg1.y1;

    bx = seg1.x2;
    by = seg1.y2;
    cx = seg2.x1;
    cy = seg2.y1;

    dx = seg2.x2;
    dy = seg2.y2;

    float dlix = (bx - ax);
    float dliy = (by - ay);
    float dljx = (dx - cx);
    float dljy = (dy - cy);

    li = sqrt((double) (dlix * dlix) + (double) (dliy * dliy));
    lj = sqrt((double) (dljx * dljx) + (double) (dljy * dljy));

    xg = (li * (double) (ax + bx) + lj * (double) (cx + dx))
        / (double) (2.0 * (li + lj));
    yg = (li * (double) (ay + by) + lj * (double) (cy + dy))
        / (double) (2.0 * (li + lj));

    if(dlix == 0.0f) thi = CV_PI / 2.0;
    else thi = atan(dliy / dlix);

    if(dljx == 0.0f) thj = CV_PI / 2.0;
    else thj = atan(dljy / dljx);

    if (fabs(thi - thj) <= CV_PI / 2.0)
    {
        thr = (li * thi + lj * thj) / (li + lj);
    }
    else
    {
        double tmp = thj - CV_PI * (thj / fabs(thj));
        thr = li * thi + lj * tmp;
        thr /= (li + lj);
    }

    axg = ((double) ay - yg) * sin(thr) + ((double) ax - xg) * cos(thr);
    bxg = ((double) by - yg) * sin(thr) + ((double) bx - xg) * cos(thr);
    cxg = ((double) cy - yg) * sin(thr) + ((double) cx - xg) * cos(thr);
    dxg = ((double) dy - yg) * sin(thr) + ((double) dx - xg) * cos(thr);

    delta1xg = min(axg,min(bxg,min(cxg,dxg)));
    delta2xg = max(axg,max(bxg,max(cxg,dxg)));

    delta1x = delta1xg * cos(thr) + xg;
    delta1y = delta1xg * sin(thr) + yg;
    delta2x = delta2xg * cos(thr) + xg;
    delta2y = delta2xg * sin(thr) + yg;

    seg_merged.x1 = (float)delta1x;
    seg_merged.y1 = (float)delta1y;
    seg_merged.x2 = (float)delta2x;
    seg_merged.y2 = (float)delta2y;
}

double Lineextractor::distPointLine(const Mat& p, Mat& l)
{
    double x = l.at<double>(0,0);
    double y = l.at<double>(1,0);
    double w = sqrt(x*x+y*y);

    l.at<double>(0,0) = x / w;
    l.at<double>(1,0) = y / w;
    l.at<double>(2,0) = l.at<double>(2,0) / w;

    return l.dot(p);
}

bool Lineextractor::mergeSegments(const SEGMENT& seg1, const SEGMENT& seg2, SEGMENT& seg_merged)
{
    double o[] = { 0.0, 0.0, 1.0 };
    double a[] = { 0.0, 0.0, 1.0 };
    double b[] = { 0.0, 0.0, 1.0 };
    double c[3];

    o[0] = ( seg2.x1 + seg2.x2 ) / 2.0;
    o[1] = ( seg2.y1 + seg2.y2 ) / 2.0;

    a[0] = seg1.x1;
    a[1] = seg1.y1;
    b[0] = seg1.x2;
    b[1] = seg1.y2;

    Mat ori = Mat(3, 1, CV_64FC1, o).clone();
    Mat p1 = Mat(3, 1, CV_64FC1, a).clone();
    Mat p2 = Mat(3, 1, CV_64FC1, b).clone();
    Mat l1 = Mat(3, 1, CV_64FC1, c).clone();

    l1 = p1.cross(p2);

    Point2f seg1mid, seg2mid;
    seg1mid.x = (seg1.x1 + seg1.x2) /2.0f;
    seg1mid.y = (seg1.y1 + seg1.y2) /2.0f;
    seg2mid.x = (seg2.x1 + seg2.x2) /2.0f;
    seg2mid.y = (seg2.y1 + seg2.y2) /2.0f;

    float seg1len = sqrt((seg1.x1 - seg1.x2)*(seg1.x1 - seg1.x2)+(seg1.y1 - seg1.y2)*(seg1.y1 - seg1.y2));
    float seg2len = sqrt((seg2.x1 - seg2.x2)*(seg2.x1 - seg2.x2)+(seg2.y1 - seg2.y2)*(seg2.y1 - seg2.y2));
    float middist = sqrt((seg1mid.x - seg2mid.x)*(seg1mid.x - seg2mid.x) + (seg1mid.y - seg2mid.y)*(seg1mid.y - seg2mid.y));
    float angdiff = fabs(seg1.angle - seg2.angle);

    float dist = (float)distPointLine(ori, l1);

    if ( fabs( dist ) <= threshold_dist * 2.0f && middist <= seg1len / 2.0f + seg2len / 2.0f + 20.0f
            && angdiff <= CV_PI / 180.0f * 5.0f)
    {
        mergeLines(seg1, seg2, seg_merged);
        return true;
    }
    else
    {
        return false;
    }
}

template<class T>
    void Lineextractor::incidentPoint(const Mat& l, T& pt)
    {
        double a[] = { (double)pt.x, (double)pt.y, 1.0 };
        double b[] = { l.at<double>(0,0), l.at<double>(1,0), 0.0 };
        double c[3];

        Mat xk = Mat(3, 1, CV_64FC1, a).clone();
        Mat lh = Mat(3, 1, CV_64FC1, b).clone();
        Mat lk = Mat(3, 1, CV_64FC1, c).clone();

        lk = xk.cross(lh);
        xk = lk.cross(l);

        xk.convertTo(xk, -1, 1.0 / xk.at<double>(2,0));

        Point2f pt_tmp;
        pt_tmp.x = (float)xk.at<double>(0,0) < 0.0f ? 0.0f : (float)xk.at<double>(0,0)
            >= (imagewidth - 1.0f) ? (imagewidth - 1.0f) : (float)xk.at<double>(0,0);
        pt_tmp.y = (float)xk.at<double>(1,0) < 0.0f ? 0.0f : (float)xk.at<double>(1,0)
            >= (imageheight - 1.0f) ? (imageheight - 1.0f) : (float)xk.at<double>(1,0);
        pt = T(pt_tmp);
    }

void Lineextractor::extractSegments(const std::vector<Point2i>& points, std::vector<SEGMENT>& segments )
{
    bool is_line;

    int i, j;
    SEGMENT seg;
    Point2i ps, pe, pt;

    std::vector<Point2i> l_points;

    int total = (int)points.size();

    for ( i = 0; i + threshold_length < total; i++ )
    {
        ps = points[i];
        pe = points[i + threshold_length];

        double a[] = { (double)ps.x, (double)ps.y, 1 };
        double b[] = { (double)pe.x, (double)pe.y, 1 };
        double c[3], d[3];

        Mat p1 = Mat(3, 1, CV_64FC1, a).clone();
        Mat p2 = Mat(3, 1, CV_64FC1, b).clone();
        Mat p = Mat(3, 1, CV_64FC1, c).clone();
        Mat l = Mat(3, 1, CV_64FC1, d).clone();
        l = p1.cross(p2);

        is_line = true;

        l_points.clear();
        l_points.push_back(ps);

        for ( j = 1; j < threshold_length; j++ )
        {
            pt.x = points[i+j].x;
            pt.y = points[i+j].y;

            p.at<double>(0,0) = (double)pt.x;
            p.at<double>(1,0) = (double)pt.y;
            p.at<double>(2,0) = 1.0;

            double dist = distPointLine(p, l);

            if ( fabs( dist ) > threshold_dist )
            {
                is_line = false;
                break;
            }
            l_points.push_back(pt);
        }

        // Line check fail, test next point
        if ( is_line == false )
            continue;

        l_points.push_back(pe);

        Vec4f line;
        fitLine( Mat(l_points), line, DIST_L2, 0, 0.01, 0.01);
        a[0] = line[2];
        a[1] = line[3];
        b[0] = line[2] + line[0];
        b[1] = line[3] + line[1];

        p1 = Mat(3, 1, CV_64FC1, a).clone();
        p2 = Mat(3, 1, CV_64FC1, b).clone();

        l = p1.cross(p2);

        incidentPoint(l, ps);

        // Extending line
        for ( j = threshold_length + 1; i + j < total; j++ )
        {
            pt.x = points[i+j].x;
            pt.y = points[i+j].y;

            p.at<double>(0,0) = (double)pt.x;
            p.at<double>(1,0) = (double)pt.y;
            p.at<double>(2,0) = 1.0;

            double dist = distPointLine(p, l);
            if ( fabs( dist ) > threshold_dist )
            {
                fitLine( Mat(l_points), line, DIST_L2, 0, 0.01, 0.01);
                a[0] = line[2];
                a[1] = line[3];
                b[0] = line[2] + line[0];
                b[1] = line[3] + line[1];

                p1 = Mat(3, 1, CV_64FC1, a).clone();
                p2 = Mat(3, 1, CV_64FC1, b).clone();

                l = p1.cross(p2);
                dist = distPointLine(p, l);
                if ( fabs( dist ) > threshold_dist ) {
                    j--;
                    break;
                }
            }
            pe = pt;
            l_points.push_back(pt);
        }
        fitLine( Mat(l_points), line, DIST_L2, 0, 0.01, 0.01);
        a[0] = line[2];
        a[1] = line[3];
        b[0] = line[2] + line[0];
        b[1] = line[3] + line[1];

        p1 = Mat(3, 1, CV_64FC1, a).clone();
        p2 = Mat(3, 1, CV_64FC1, b).clone();

        l = p1.cross(p2);

        Point2f e1, e2;
        e1.x = (float)ps.x;
        e1.y = (float)ps.y;
        e2.x = (float)pe.x;
        e2.y = (float)pe.y;

        incidentPoint(l, e1);
        incidentPoint(l, e2);
        seg.x1 = e1.x;
        seg.y1 = e1.y;
        seg.x2 = e2.x;
        seg.y2 = e2.y;

        segments.push_back(seg);
        i = i + j;
    }
}

void Lineextractor::pointInboardTest(const Mat& src, Point2i& pt)
{
    pt.x = pt.x <= 5 ? 5 : pt.x >= src.cols - 5 ? src.cols - 5 : pt.x;
    pt.y = pt.y <= 5 ? 5 : pt.y >= src.rows - 5 ? src.rows - 5 : pt.y;
}

bool Lineextractor::getPointChain(const Mat& img, Point pt,
        Point& chained_pt, float& direction, int step)
{
    int ri, ci;
    int indices[8][2] = { {1,1}, {1,0}, {1,-1}, {0,-1},
        {-1,-1},{-1,0}, {-1,1}, {0,1} };

    float min_dir_diff = 7.0f;
    Point consistent_pt;
    int consistent_direction = 0;
    for ( int i = 0; i < 8; i++ )
    {
        ci = pt.x + indices[i][1];
        ri = pt.y + indices[i][0];

        if ( ri < 0 || ri == img.rows || ci < 0 || ci == img.cols )
            continue;

        if ( img.at<unsigned char>(ri, ci) == 0 )
            continue;

        if(step == 0)
        {
            chained_pt.x = ci;
            chained_pt.y = ri;
            // direction = (float)i;
            direction = i > 4 ? (float)(i - 8) : (float)i;
            return true;
        }
        else
        {
            float curr_dir = i > 4 ? (float)(i - 8) : (float)i;
            float dir_diff = abs(curr_dir - direction);
            dir_diff = dir_diff > 4.0f ? 8.0f - dir_diff : dir_diff;
            if(dir_diff <= min_dir_diff)
            {
                min_dir_diff = dir_diff;
                consistent_pt.x = ci;
                consistent_pt.y = ri;
                consistent_direction = i > 4 ? i - 8 : i;
            }
        }
    }
    if(min_dir_diff < 2.0f)
    {
        chained_pt.x = consistent_pt.x;
        chained_pt.y = consistent_pt.y;
        direction = (direction * (float)step + (float)consistent_direction)
            / (float)(step + 1);
        return true;
    }
    return false;
}

void Lineextractor::lineDetection(const Mat& src, std::vector<SEGMENT>& segments_all)
{
    int r, c;
    imageheight=src.rows; imagewidth=src.cols;

    std::vector<Point2i> points;
    std::vector<SEGMENT> segments, segments_tmp;
    Mat canny;
    Canny(src, canny, canny_th1, canny_th2, canny_aperture_size);

    canny.colRange(0,6).rowRange(0,6) = 0;
    canny.colRange(src.cols-5,src.cols).rowRange(src.rows-5,src.rows) = 0;

    SEGMENT seg, seg1, seg2;

    for ( r = 0; r < imageheight; r++ )
    {
        for ( c = 0; c < imagewidth; c++ )
        {
            // Find seeds - skip for non-seeds
            if ( canny.at<unsigned char>(r,c) == 0 )
                continue;

            // Found seeds
            Point2i pt = Point2i(c,r);

            points.push_back(pt);
            canny.at<unsigned char>(pt.y, pt.x) = 0;

            float direction = 0.0f;
            int step = 0;
            while(getPointChain(canny, pt, pt, direction, step))
            {
                points.push_back(pt);
                step++;
                canny.at<unsigned char>(pt.y, pt.x) = 0;
            }

            if ( points.size() < (unsigned int)threshold_length + 1 )
            {
                points.clear();
                continue;
            }

            extractSegments(points, segments);

            if ( segments.size() == 0 )
            {
                points.clear();
                continue;
            }
            for ( int i = 0; i < (int)segments.size(); i++ )
            {
                seg = segments[i];
                float length = sqrt((seg.x1 - seg.x2)*(seg.x1 - seg.x2) +
                        (seg.y1 - seg.y2)*(seg.y1 - seg.y2));
                if(length < threshold_length)
                    continue;
                if( (seg.x1 <= 5.0f && seg.x2 <= 5.0f) ||
                    (seg.y1 <= 5.0f && seg.y2 <= 5.0f) ||
                    (seg.x1 >= imagewidth - 5.0f && seg.x2 >= imagewidth - 5.0f) ||
                    (seg.y1 >= imageheight - 5.0f && seg.y2 >= imageheight - 5.0f) )
                    continue;
                additionalOperationsOnSegment(src, seg);
                if(!do_merge)
                    segments_all.push_back(seg);
                segments_tmp.push_back(seg);
            }
            points.clear();
            segments.clear();
        }
    }
    if(!do_merge)
        return;

    bool is_merged = false;
    int ith = (int)segments_tmp.size() - 1;
    int jth = ith - 1;
    while(ith > 1 || jth > 0)
    {
        seg1 = segments_tmp[ith];
        seg2 = segments_tmp[jth];
        SEGMENT seg_merged;
        is_merged = mergeSegments(seg1, seg2, seg_merged);
        if(is_merged == true)
        {
            seg2 = seg_merged;
            additionalOperationsOnSegment(src, seg2);
            std::vector<SEGMENT>::iterator it = segments_tmp.begin() + ith;
            *it = seg2;
            segments_tmp.erase(segments_tmp.begin()+jth);
            ith--;
            jth = ith - 1;
        }
        else
        {
            jth--;
        }
        if(jth < 0) {
            ith--;
            jth = ith - 1;
        }
    }
    segments_all = segments_tmp;
}

inline void Lineextractor::getAngle(SEGMENT& seg)
{
    seg.angle = (float)(fastAtan2(seg.y2 - seg.y1, seg.x2 - seg.x1) / 180.0f * CV_PI);
}

void Lineextractor::additionalOperationsOnSegment(const Mat& src, SEGMENT& seg)
{
    if(seg.x1 == 0.0f && seg.x2 == 0.0f && seg.y1 == 0.0f && seg.y2 == 0.0f)
        return;

    getAngle(seg);
    double ang = (double)seg.angle;

    Point2f start = Point2f(seg.x1, seg.y1);
    Point2f end = Point2f(seg.x2, seg.y2);

    double dx = 0.0, dy = 0.0;
    dx = (double) end.x - (double) start.x;
    dy = (double) end.y - (double) start.y;

    int num_points = 10;
    Point2f *points = new Point2f[num_points];

    points[0] = start;
    points[num_points - 1] = end;
    for (int i = 0; i < num_points; i++)
    {
        if (i == 0 || i == num_points - 1)
            continue;
        points[i].x = points[0].x + ((float)dx / float(num_points - 1) * (float) i);
        points[i].y = points[0].y + ((float)dy / float(num_points - 1) * (float) i);
    }

    Point2i *points_right = new Point2i[num_points];
    Point2i *points_left = new Point2i[num_points];
    double gap = 1.0;

    for(int i = 0; i < num_points; i++)
    {
        points_right[i].x = cvRound(points[i].x + gap*cos(90.0 * CV_PI / 180.0 + ang));
        points_right[i].y = cvRound(points[i].y + gap*sin(90.0 * CV_PI / 180.0 + ang));
        points_left[i].x = cvRound(points[i].x - gap*cos(90.0 * CV_PI / 180.0 + ang));
        points_left[i].y = cvRound(points[i].y - gap*sin(90.0 * CV_PI / 180.0 + ang));
        pointInboardTest(src, points_right[i]);
        pointInboardTest(src, points_left[i]);
    }

    int iR = 0, iL = 0;
    for(int i = 0; i < num_points; i++)
    {
        iR += src.at<unsigned char>(points_right[i].y, points_right[i].x);
        iL += src.at<unsigned char>(points_left[i].y, points_left[i].x);
    }

    if(iR > iL)
    {
        std::swap(seg.x1, seg.x2);
        std::swap(seg.y1, seg.y2);
        getAngle(seg);
    }

    delete[] points;
    delete[] points_right;
    delete[] points_left;

    return;
}

void Lineextractor::drawSegment(Mat& mat, const SEGMENT& seg, Scalar bgr, int thickness, bool directed)
{
    double gap = 10.0;
    double ang = (double)seg.angle;
    double arrow_angle = 30.0;

    Point2i p1;
    p1.x = cvRound(seg.x2 - gap*cos(arrow_angle * CV_PI / 180.0 + ang));
    p1.y = cvRound(seg.y2 - gap*sin(arrow_angle * CV_PI / 180.0 + ang));
    pointInboardTest(mat, p1);

    line(mat, Point(cvRound(seg.x1), cvRound(seg.y1)),
            Point(cvRound(seg.x2), cvRound(seg.y2)), bgr, thickness, 1);
    if(directed)
        line(mat, Point(cvRound(seg.x2), cvRound(seg.y2)), p1, bgr, thickness, 1);
}

}//namespace ORB_SLAM2
