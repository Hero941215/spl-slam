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

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "Map.h"
#include "MapPoint.h"
#include "KeyFrame.h"
#include "LoopClosing.h"
#include "Frame.h"

#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"
#include "Thirdparty/g2o/g2o/core/base_vertex.h"
#include "Thirdparty/g2o/g2o/core/base_unary_edge.h"

#include <Eigen/Geometry>

using namespace Eigen;

typedef Matrix<double,6,1> Vector6d;
typedef Matrix<double,6,6> Matrix6d;

namespace PL_SLAM
{

class LoopClosing;

//设置前端使用的仅线优化的边模板参数
class  EdgeSE3ProjectXYZOnlyPoseLines: public  g2o::BaseUnaryEdge<3, Vector3d, g2o::VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeSE3ProjectXYZOnlyPoseLines(){}

  bool read(std::istream& is);
  bool write(std::ostream& os) const;

  void computeError()  {
    const g2o::VertexSE3Expmap* v1 = static_cast<const g2o::VertexSE3Expmap*>(_vertices[0]);
    Vector3d obs(_measurement);
    Vector2d uv = cam_project(v1->estimate().map(Xw));
    _error(0,0) = obs(0) * uv(0) + obs(1) * uv(1) + obs(2);
    _error(1,0) = 0;
    _error(2,0) = 0;
  }

  bool isDepthPositive() {
    const g2o::VertexSE3Expmap* v1 = static_cast<const g2o::VertexSE3Expmap*>(_vertices[0]);
    return (v1->estimate().map(Xw))(2)>0.0;
  }

  virtual void linearizeOplus();

  Vector2d cam_project(const Vector3d & trans_xyz) const;

  Vector3d Xw;
  double fx, fy, cx, cy;
};

//地图线端点构成的边
class  EdgeSE3ProjectXYZLines: public  g2o::BaseBinaryEdge<3, Vector3d, g2o::VertexSBAPointXYZ, g2o::VertexSE3Expmap>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgeSE3ProjectXYZLines();

  bool read(std::istream& is);
  bool write(std::ostream& os) const;

  void computeError()  {
    const g2o::VertexSE3Expmap* v1 = static_cast<const g2o::VertexSE3Expmap*>(_vertices[1]);
    const g2o::VertexSBAPointXYZ* v2 = static_cast<const g2o::VertexSBAPointXYZ*>(_vertices[0]);
    Vector3d obs(_measurement);
    Vector2d uv = cam_project(v1->estimate().map(v2->estimate()));
    _error(0,0) = obs(0) * uv(0) + obs(1) * uv(1) + obs(2);
    _error(1,0) = 0;
    _error(2,0) = 0;
  }
  
  bool isDepthPositive() {
    const g2o::VertexSE3Expmap* v1 = static_cast<const g2o::VertexSE3Expmap*>(_vertices[1]);
    const g2o::VertexSBAPointXYZ* v2 = static_cast<const g2o::VertexSBAPointXYZ*>(_vertices[0]);
    return (v1->estimate().map(v2->estimate()))(2)>0.0;
  }
    
  virtual void linearizeOplus();

  Vector2d cam_project(const Vector3d & trans_xyz) const;

  double fx, fy, cx, cy;
};

class Optimizer
{
public:
    void static BundleAdjustment(const std::vector<KeyFrame*> &vpKF, const std::vector<MapPoint*> &vpMP,
                                 int nIterations = 5, bool *pbStopFlag=NULL, const unsigned long nLoopKF=0,
                                 const bool bRobust = true);
    void static GlobalBundleAdjustemnt(Map* pMap, int nIterations=5, bool *pbStopFlag=NULL,
                                       const unsigned long nLoopKF=0, const bool bRobust = true);
    void static LocalBundleAdjustment(KeyFrame* pKF, bool *pbStopFlag, Map *pMap);
    
    int static PoseOptimization(Frame* pFrame);

    // if bFixScale is true, 6DoF optimization (stereo,rgbd), 7DoF otherwise (mono)
    void static OptimizeEssentialGraph(Map* pMap, KeyFrame* pLoopKF, KeyFrame* pCurKF,
                                       const LoopClosing::KeyFrameAndPose &NonCorrectedSim3,
                                       const LoopClosing::KeyFrameAndPose &CorrectedSim3,
                                       const map<KeyFrame *, set<KeyFrame *> > &LoopConnections,
                                       const bool &bFixScale);

    // if bFixScale is true, optimize SE3 (stereo,rgbd), Sim3 otherwise (mono)
    static int OptimizeSim3(KeyFrame* pKF1, KeyFrame* pKF2, std::vector<MapPoint *> &vpMatches1,
                            g2o::Sim3 &g2oS12, const float th2, const bool bFixScale);
    

    // 单目初始化：（系统运行只使用一次）
    void static GlobalBundleAdjustemntIni(Map* pMap, int nIterations=5, const bool bRobust = true);
    void static BundleAdjustmentPointsIni(std::vector<KeyFrame*> &vpKFs, std::vector<MapPoint*> &vpMP, int nIterations, const bool bRobust, float &PoseErrorForPoints);
    void static BundleAdjustmentLinesIni(std::vector<KeyFrame*> &vpKFs, std::vector<MapLine*> &vpML, int nIterations, const bool bRobust, float &PoseErrorForLines);
    void static BundleAdjustmentBothIni(std::vector<KeyFrame*> &vpKFs, std::vector<MapPoint*> &vpMP, std::vector<MapLine*> &vpML, int nIterations = 5, const bool bRobust = true);
    
    //复合点线特征位姿优化函数：（前端）
    void static PoseOptimizationmain(Frame* pFrame, int &inlinersNum, int &inlinerLinesNum, int nIterations = 10 );
    void static PoseOptimizationDoublePoints(Frame* pFrame, int &inlinersNum ); //内部子函数决定迭代次数
    void static PoseOptimizationLowFeature(Frame* pFrame, int &inlinersNum, int &inlinerLinesNum, int nIterations = 10 );
    
    //复合优化子函数：
    void static PoseOptimizationPoints(Frame* pFrame, float &PoseErrorForPoints, int &nmatches, int nIterations = 5 );
    void static PoseOptimizationLines(Frame* pFrame, float &PoseErrorForLines, int &nmatchLines, int nIterations = 5 );
    void static PoseOptimizationBothOld(Frame* pFrame, int &inlinersNum, int &inlinerLinesNum, int nIterations = 5 );
    void static PoseOptimizationBoth(Frame* pFrame, int &inlinersNum, int &inlinerLinesNum, int nIterations = 10 );
  
    void static SetOutliersForPose(Frame* pFrame, int &inlinersNum);
    void static SetOutlierLinesForPose(Frame* pFrame, int &inlinerLinesNum);  
    
    //复合局部BA优化函数：（后端）
    void static LocalBundleAdjustmentmain(KeyFrame* pKF, bool *pbStopFlag, Map *pMap);  
    //两次均同时使用点线特征进行优化,与LocalBundleAdjustment()思想类似
    void static LocalBundleAdjustmentmainOld(KeyFrame* pKF, bool *pbStopFlag, Map *pMap);
    //先执行点优化，再执行点线同时的局部BA优化
    void static LocalBundleAdjustmentmainNew(KeyFrame* pKF, bool *pbStopFlag, Map *pMap);
    void static LocalBundleAdjustmentDoubleLines(KeyFrame* pKF, bool *pbStopFlag, Map *pMap);
    void static LocalBundleAdjustmentPoints(KeyFrame* pKF, bool *pbStopFlag, Map *pMap, std::map<KeyFrame*,float> &PoseUnitErrorForPointsLocalBA);
    void static LocalBundleAdjustmentLines(KeyFrame* pKF, bool *pbStopFlag, Map *pMap, std::map<KeyFrame*,float> &PoseUnitErrorForLinesLocalBA);
    void static LocalBundleAdjustmentBoth(KeyFrame* pKF, bool *pbStopFlag, Map *pMap);
    
protected:
  
    //a.前端：
    void static GaussNewtonOptimizationForPose(Frame* pFrame, int nIterations = 5);
    void static OptimizeFunctionForPose(Frame* pFrame, g2o::SE3Quat &TcwBothUpdate, Matrix6d &H, Vector6d &g, double &error);
    float static RobustWeightHuber( float PointReProjectionErrorNorm2, float deta);
    void static RecheckBothOutliersForPose(Frame* pFrame, int &inlinersNum, int &inlinerLinesNum);
    
    //b.后端：
    void static LocalBAPoseDecidingBetweenLinesAndPoints(KeyFrame* pKF, Map *pMap, std::map<KeyFrame*,float> &PoseUnitErrorForPointsLocalBA, 
							 std::map<KeyFrame*,float > &PoseUnitErrorForLinesLocalBA );
    
};

} //namespace ORB_SLAM

#endif // OPTIMIZER_H
