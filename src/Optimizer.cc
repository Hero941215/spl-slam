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

#include "Optimizer.h"

#include "Thirdparty/g2o/g2o/core/block_solver.h"
#include "Thirdparty/g2o/g2o/core/optimization_algorithm_levenberg.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_eigen.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"
#include "Thirdparty/g2o/g2o/core/robust_kernel_impl.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_dense.h"
#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"

#include<Eigen/StdVector>

#include "Converter.h"

#include <map>

#include<mutex>
#include<thread> 

namespace PL_SLAM
{

Vector2d project2d(const Vector3d& v)  {
  Vector2d res;
  res(0) = v(0)/v(2);
  res(1) = v(1)/v(2);
  return res;
} 
//前端边定义：
bool EdgeSE3ProjectXYZOnlyPoseLines::read(std::istream& is){
  for (int i=0; i<3; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<3; i++)
    for (int j=i; j<3; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeSE3ProjectXYZOnlyPoseLines::write(std::ostream& os) const {

  for (int i=0; i<3; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<3; i++)
    for (int j=i; j<3; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}

Vector2d EdgeSE3ProjectXYZOnlyPoseLines::cam_project(const Vector3d & trans_xyz) const{
  Vector2d proj = project2d(trans_xyz);
  Vector2d res;
  res[0] = proj[0]*fx + cx;
  res[1] = proj[1]*fy + cy;
  return res;
}

void EdgeSE3ProjectXYZOnlyPoseLines::linearizeOplus() {
  g2o::VertexSE3Expmap * vi = static_cast<g2o::VertexSE3Expmap *>(_vertices[0]);
  Vector3d xyz_trans = vi->estimate().map(Xw);

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double invz = 1.0/xyz_trans[2];
  double invz_2 = invz*invz;

  Vector3d obs(_measurement);

  _jacobianOplusXi(0,0) =  -x*y*invz_2*obs(0)*fx - (1+y*y*invz_2)*obs(1)*fy;
  _jacobianOplusXi(0,1) = (1+(x*x*invz_2)) *obs(0)*fx + x*y*invz_2 *obs(1)*fy;
  _jacobianOplusXi(0,2) = -y*invz *obs(0)*fx + x*invz *obs(1)*fy;
  _jacobianOplusXi(0,3) = invz *obs(0)*fx;
  _jacobianOplusXi(0,4) = invz *obs(1)*fy;
  _jacobianOplusXi(0,5) = -x*invz_2 *obs(0)*fx -y*invz_2 *obs(1)*fy;

  _jacobianOplusXi(1,0) = 0;
  _jacobianOplusXi(1,1) = 0;
  _jacobianOplusXi(1,2) = 0;
  _jacobianOplusXi(1,3) = 0;
  _jacobianOplusXi(1,4) = 0;
  _jacobianOplusXi(1,5) = 0;
  
  _jacobianOplusXi(2,0) = 0;
  _jacobianOplusXi(2,1) = 0;
  _jacobianOplusXi(2,2) = 0;
  _jacobianOplusXi(2,3) = 0;
  _jacobianOplusXi(2,4) = 0;
  _jacobianOplusXi(2,5) = 0;
  
}

//初始化、后端边定义
EdgeSE3ProjectXYZLines::EdgeSE3ProjectXYZLines() : BaseBinaryEdge<3, Vector3d, g2o::VertexSBAPointXYZ, g2o::VertexSE3Expmap>() {
}

bool EdgeSE3ProjectXYZLines::read(std::istream& is){
  for (int i=0; i<3; i++){
    is >> _measurement[i];
  }
  for (int i=0; i<3; i++)
    for (int j=i; j<3; j++) {
      is >> information()(i,j);
      if (i!=j)
        information()(j,i)=information()(i,j);
    }
  return true;
}

bool EdgeSE3ProjectXYZLines::write(std::ostream& os) const {

  for (int i=0; i<3; i++){
    os << measurement()[i] << " ";
  }

  for (int i=0; i<3; i++)
    for (int j=i; j<3; j++){
      os << " " <<  information()(i,j);
    }
  return os.good();
}


Vector2d EdgeSE3ProjectXYZLines::cam_project(const Vector3d & trans_xyz) const{
  Vector2d proj = project2d(trans_xyz);
  Vector2d res;
  res[0] = proj[0]*fx + cx;
  res[1] = proj[1]*fy + cy;
  return res;
}

void EdgeSE3ProjectXYZLines::linearizeOplus() {

  g2o::VertexSE3Expmap * vj = static_cast<g2o::VertexSE3Expmap *>(_vertices[1]);
  g2o::SE3Quat T(vj->estimate());

  g2o::VertexSBAPointXYZ* vi = static_cast<g2o::VertexSBAPointXYZ*>(_vertices[0]);
  Vector3d xyz = vi->estimate();
  Vector3d xyz_trans = T.map(xyz);

  double x = xyz_trans[0];
  double y = xyz_trans[1];
  double z = xyz_trans[2];
  double invz = 1. / z;
  double invz_2 = invz*invz;

  Vector3d obs(_measurement);
  
  Matrix<double,3,3> tmp;
  tmp(0,0) = obs(0)*fx;
  tmp(0,1) = obs(1)*fy;
  tmp(0,2) = -x/z*obs(0)*fx - y/z*obs(1)*fy;

  tmp(1,0) = 0;
  tmp(1,1) = 0;
  tmp(1,2) = 0;
  
  tmp(2,0) = 0;
  tmp(2,1) = 0;
  tmp(2,2) = 0;

  _jacobianOplusXi =  1./z * tmp * T.rotation().toRotationMatrix();

  //图像帧位姿的雅克比矩阵
  _jacobianOplusXj(0,0) =  -x*y*invz_2*obs(0)*fx - (1+y*y*invz_2)*obs(1)*fy;
  _jacobianOplusXj(0,1) = (1+(x*x*invz_2)) *obs(0)*fx + x*y*invz_2 *obs(1)*fy;
  _jacobianOplusXj(0,2) = -y*invz *obs(0)*fx + x*invz *obs(1)*fy;
  _jacobianOplusXj(0,3) = invz *obs(0)*fx;
  _jacobianOplusXj(0,4) = invz *obs(1)*fy;
  _jacobianOplusXj(0,5) = -x*invz_2 *obs(0)*fx -y*invz_2 *obs(1)*fy;

  _jacobianOplusXj(1,0) = 0;
  _jacobianOplusXj(1,1) = 0;
  _jacobianOplusXj(1,2) = 0;
  _jacobianOplusXj(1,3) = 0;
  _jacobianOplusXj(1,4) = 0;
  _jacobianOplusXj(1,5) = 0;
  
  _jacobianOplusXj(2,0) = 0;
  _jacobianOplusXj(2,1) = 0;
  _jacobianOplusXj(2,2) = 0;
  _jacobianOplusXj(2,3) = 0;
  _jacobianOplusXj(2,4) = 0;
  _jacobianOplusXj(2,5) = 0;
}

void Optimizer::GlobalBundleAdjustemnt(Map* pMap, int nIterations, bool* pbStopFlag, const unsigned long nLoopKF, const bool bRobust)
{
    vector<KeyFrame*> vpKFs = pMap->GetAllKeyFrames();
    vector<MapPoint*> vpMP = pMap->GetAllMapPoints();
    BundleAdjustment(vpKFs,vpMP,nIterations,pbStopFlag, nLoopKF, bRobust);
}

void Optimizer::BundleAdjustment(const vector<KeyFrame *> &vpKFs, const vector<MapPoint *> &vpMP,
                                 int nIterations, bool* pbStopFlag, const unsigned long nLoopKF, const bool bRobust)
{
    vector<bool> vbNotIncludedMP;
    vbNotIncludedMP.resize(vpMP.size());

    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    if(pbStopFlag)
        optimizer.setForceStopFlag(pbStopFlag);

    long unsigned int maxKFid = 0;

    // Set KeyFrame vertices
    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];
        if(pKF->isBad())
            continue;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKF->GetPose()));
        vSE3->setId(pKF->mnId);
        vSE3->setFixed(pKF->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKF->mnId>maxKFid)
            maxKFid=pKF->mnId;
    }

    const float thHuber2D = sqrt(5.99);
    const float thHuber3D = sqrt(7.815);

    // Set MapPoint vertices
    for(size_t i=0; i<vpMP.size(); i++)
    {
        MapPoint* pMP = vpMP[i];
        if(pMP->isBad())
            continue;
        g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
        vPoint->setEstimate(Converter::toVector3d(pMP->GetWorldPos()));
        const int id = pMP->mnId+maxKFid+1;
        vPoint->setId(id);
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

       const map<KeyFrame*,size_t> observations = pMP->GetObservations();

        int nEdges = 0;
        //SET EDGES
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(); mit!=observations.end(); mit++)
        {

            KeyFrame* pKF = mit->first;
            if(pKF->isBad() || pKF->mnId>maxKFid)
                continue;

            nEdges++;

            const cv::KeyPoint &kpUn = pKF->mvKeysUn[mit->second];

            if(pKF->mvuRight[mit->second]<0)
            {
                Eigen::Matrix<double,2,1> obs;
                obs << kpUn.pt.x, kpUn.pt.y;

                g2o::EdgeSE3ProjectXYZ* e = new g2o::EdgeSE3ProjectXYZ();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKF->mnId)));
                e->setMeasurement(obs);
                const float &invSigma2 = pKF->mvInvLevelSigma2[kpUn.octave];
                e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                if(bRobust)
                {
                    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuber2D);
                }

                e->fx = pKF->fx;
                e->fy = pKF->fy;
                e->cx = pKF->cx;
                e->cy = pKF->cy;

                optimizer.addEdge(e);
            }
            else
            {
                Eigen::Matrix<double,3,1> obs;
                const float kp_ur = pKF->mvuRight[mit->second];
                obs << kpUn.pt.x, kpUn.pt.y, kp_ur;

                g2o::EdgeStereoSE3ProjectXYZ* e = new g2o::EdgeStereoSE3ProjectXYZ();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKF->mnId)));
                e->setMeasurement(obs);
                const float &invSigma2 = pKF->mvInvLevelSigma2[kpUn.octave];
                Eigen::Matrix3d Info = Eigen::Matrix3d::Identity()*invSigma2;
                e->setInformation(Info);

                if(bRobust)
                {
                    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuber3D);
                }

                e->fx = pKF->fx;
                e->fy = pKF->fy;
                e->cx = pKF->cx;
                e->cy = pKF->cy;
                e->bf = pKF->mbf;

                optimizer.addEdge(e);
            }
        }

        if(nEdges==0)
        {
            optimizer.removeVertex(vPoint);
            vbNotIncludedMP[i]=true;
        }
        else
        {
            vbNotIncludedMP[i]=false;
        }
    }

    // Optimize!
    optimizer.initializeOptimization();
    optimizer.optimize(nIterations);

    // Recover optimized data

    //Keyframes
    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];
        if(pKF->isBad())
            continue;
        g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
        if(nLoopKF==0)
        {
            pKF->SetPose(Converter::toCvMat(SE3quat));
        }
        else
        {
            pKF->mTcwGBA.create(4,4,CV_32F);
            Converter::toCvMat(SE3quat).copyTo(pKF->mTcwGBA);
            pKF->mnBAGlobalForKF = nLoopKF;
        }
    }

    //Points
    for(size_t i=0; i<vpMP.size(); i++)
    {
        if(vbNotIncludedMP[i])
            continue;

        MapPoint* pMP = vpMP[i];

        if(pMP->isBad())
            continue;
        g2o::VertexSBAPointXYZ* vPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(pMP->mnId+maxKFid+1));

        if(nLoopKF==0)
        {
            pMP->SetWorldPos(Converter::toCvMat(vPoint->estimate()));
            pMP->UpdateNormalAndDepth();
        }
        else
        {
            pMP->mPosGBA.create(3,1,CV_32F);
            Converter::toCvMat(vPoint->estimate()).copyTo(pMP->mPosGBA);
            pMP->mnBAGlobalForKF = nLoopKF;
        }
    }

}

int Optimizer::PoseOptimization(Frame *pFrame)
{
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    int nInitialCorrespondences=0;

    // Set Frame vertex
    g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
    vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
    vSE3->setId(0);
    vSE3->setFixed(false);
    optimizer.addVertex(vSE3);

    // Set MapPoint vertices
    const int N = pFrame->N;

    vector<g2o::EdgeSE3ProjectXYZOnlyPose*> vpEdgesMono;
    vector<size_t> vnIndexEdgeMono;
    vpEdgesMono.reserve(N);
    vnIndexEdgeMono.reserve(N);

    vector<g2o::EdgeStereoSE3ProjectXYZOnlyPose*> vpEdgesStereo;
    vector<size_t> vnIndexEdgeStereo;
    vpEdgesStereo.reserve(N);
    vnIndexEdgeStereo.reserve(N);

    const float deltaMono = sqrt(5.991);
    const float deltaStereo = sqrt(7.815);


    {
    unique_lock<mutex> lock(MapPoint::mGlobalMutex);

    for(int i=0; i<N; i++)
    {
        MapPoint* pMP = pFrame->mvpMapPoints[i];
        if(pMP)
        {
            // Monocular observation
            if(pFrame->mvuRight[i]<0)
            {
                nInitialCorrespondences++;
                pFrame->mvbOutlier[i] = false;

                Eigen::Matrix<double,2,1> obs;
                const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
                obs << kpUn.pt.x, kpUn.pt.y;

                g2o::EdgeSE3ProjectXYZOnlyPose* e = new g2o::EdgeSE3ProjectXYZOnlyPose();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
                e->setMeasurement(obs);
                const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
                e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                e->setRobustKernel(rk);
                rk->setDelta(deltaMono);

                e->fx = pFrame->fx;
                e->fy = pFrame->fy;
                e->cx = pFrame->cx;
                e->cy = pFrame->cy;
                cv::Mat Xw = pMP->GetWorldPos();
                e->Xw[0] = Xw.at<float>(0);
                e->Xw[1] = Xw.at<float>(1);
                e->Xw[2] = Xw.at<float>(2);

                optimizer.addEdge(e);

                vpEdgesMono.push_back(e);
                vnIndexEdgeMono.push_back(i);
            }
            else  // Stereo observation
            {
                nInitialCorrespondences++;
                pFrame->mvbOutlier[i] = false;

                //SET EDGE
                Eigen::Matrix<double,3,1> obs;
                const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
                const float &kp_ur = pFrame->mvuRight[i];
                obs << kpUn.pt.x, kpUn.pt.y, kp_ur;

                g2o::EdgeStereoSE3ProjectXYZOnlyPose* e = new g2o::EdgeStereoSE3ProjectXYZOnlyPose();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
                e->setMeasurement(obs);
                const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
                Eigen::Matrix3d Info = Eigen::Matrix3d::Identity()*invSigma2;
                e->setInformation(Info);

                g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                e->setRobustKernel(rk);
                rk->setDelta(deltaStereo);

                e->fx = pFrame->fx;
                e->fy = pFrame->fy;
                e->cx = pFrame->cx;
                e->cy = pFrame->cy;
                e->bf = pFrame->mbf;
                cv::Mat Xw = pMP->GetWorldPos();
                e->Xw[0] = Xw.at<float>(0);
                e->Xw[1] = Xw.at<float>(1);
                e->Xw[2] = Xw.at<float>(2);

                optimizer.addEdge(e);

                vpEdgesStereo.push_back(e);
                vnIndexEdgeStereo.push_back(i);
            }
        }

    }
    }


    if(nInitialCorrespondences<3)
        return 0;

    // We perform 4 optimizations, after each optimization we classify observation as inlier/outlier
    // At the next optimization, outliers are not included, but at the end they can be classified as inliers again.
    const float chi2Mono[4]={5.991,5.991,5.991,5.991};
    const float chi2Stereo[4]={7.815,7.815,7.815, 7.815};
    const int its[4]={10,10,10,10};    

    int nBad=0;
    for(size_t it=0; it<4; it++)
    {

        vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
        optimizer.initializeOptimization(0);
        optimizer.optimize(its[it]);

        nBad=0;
        for(size_t i=0, iend=vpEdgesMono.size(); i<iend; i++)
        {
            g2o::EdgeSE3ProjectXYZOnlyPose* e = vpEdgesMono[i];

            const size_t idx = vnIndexEdgeMono[i];

            if(pFrame->mvbOutlier[idx])
            {
                e->computeError();
            }

            const float chi2 = e->chi2();

            if(chi2>chi2Mono[it])
            {                
                pFrame->mvbOutlier[idx]=true;
                e->setLevel(1);
                nBad++;
            }
            else
            {
                pFrame->mvbOutlier[idx]=false;
                e->setLevel(0);
            }

            if(it==2)
                e->setRobustKernel(0);
        }

        for(size_t i=0, iend=vpEdgesStereo.size(); i<iend; i++)
        {
            g2o::EdgeStereoSE3ProjectXYZOnlyPose* e = vpEdgesStereo[i];

            const size_t idx = vnIndexEdgeStereo[i];

            if(pFrame->mvbOutlier[idx])
            {
                e->computeError();
            }

            const float chi2 = e->chi2();

            if(chi2>chi2Stereo[it])
            {
                pFrame->mvbOutlier[idx]=true;
                e->setLevel(1);
                nBad++;
            }
            else
            {                
                e->setLevel(0);
                pFrame->mvbOutlier[idx]=false;
            }

            if(it==2)
                e->setRobustKernel(0);
        }

        if(optimizer.edges().size()<10)
            break;
    }    

    // Recover optimized pose and return number of inliers
    g2o::VertexSE3Expmap* vSE3_recov = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
    g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
    cv::Mat pose = Converter::toCvMat(SE3quat_recov);
    pFrame->SetPose(pose);

    return nInitialCorrespondences-nBad;
}

void Optimizer::LocalBundleAdjustment(KeyFrame *pKF, bool* pbStopFlag, Map* pMap)
{    
    // Local KeyFrames: First Breath Search from Current Keyframe
    list<KeyFrame*> lLocalKeyFrames;

    lLocalKeyFrames.push_back(pKF);
    pKF->mnBALocalForKF = pKF->mnId;

    const vector<KeyFrame*> vNeighKFs = pKF->GetVectorCovisibleKeyFrames();
    for(int i=0, iend=vNeighKFs.size(); i<iend; i++)
    {
        KeyFrame* pKFi = vNeighKFs[i];
        pKFi->mnBALocalForKF = pKF->mnId;
        if(!pKFi->isBad())
            lLocalKeyFrames.push_back(pKFi);
    }

    // Local MapPoints seen in Local KeyFrames
    list<MapPoint*> lLocalMapPoints;
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        vector<MapPoint*> vpMPs = (*lit)->GetMapPointMatches();
        for(vector<MapPoint*>::iterator vit=vpMPs.begin(), vend=vpMPs.end(); vit!=vend; vit++)
        {
            MapPoint* pMP = *vit;
            if(pMP)
                if(!pMP->isBad())
                    if(pMP->mnBALocalForKF!=pKF->mnId)
                    {
                        lLocalMapPoints.push_back(pMP);
                        pMP->mnBALocalForKF=pKF->mnId;
                    }
        }
    }

    // Fixed Keyframes. Keyframes that see Local MapPoints but that are not Local Keyframes
    list<KeyFrame*> lFixedCameras;
    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        map<KeyFrame*,size_t> observations = (*lit)->GetObservations();
        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(pKFi->mnBALocalForKF!=pKF->mnId && pKFi->mnBAFixedForKF!=pKF->mnId)
            {                
                pKFi->mnBAFixedForKF=pKF->mnId;
                if(!pKFi->isBad())
                    lFixedCameras.push_back(pKFi);
            }
        }
    }

    // Setup optimizer
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    if(pbStopFlag)
        optimizer.setForceStopFlag(pbStopFlag);

    unsigned long maxKFid = 0;

    // Set Local KeyFrame vertices
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(pKFi->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }

    // Set Fixed KeyFrame vertices
    for(list<KeyFrame*>::iterator lit=lFixedCameras.begin(), lend=lFixedCameras.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(true);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }

    // Set MapPoint vertices
    const int nExpectedSize = (lLocalKeyFrames.size()+lFixedCameras.size())*lLocalMapPoints.size();

    vector<g2o::EdgeSE3ProjectXYZ*> vpEdgesMono;
    vpEdgesMono.reserve(nExpectedSize);

    vector<KeyFrame*> vpEdgeKFMono;
    vpEdgeKFMono.reserve(nExpectedSize);

    vector<MapPoint*> vpMapPointEdgeMono;
    vpMapPointEdgeMono.reserve(nExpectedSize);

    vector<g2o::EdgeStereoSE3ProjectXYZ*> vpEdgesStereo;
    vpEdgesStereo.reserve(nExpectedSize);

    vector<KeyFrame*> vpEdgeKFStereo;
    vpEdgeKFStereo.reserve(nExpectedSize);

    vector<MapPoint*> vpMapPointEdgeStereo;
    vpMapPointEdgeStereo.reserve(nExpectedSize);

    const float thHuberMono = sqrt(5.991);
    const float thHuberStereo = sqrt(7.815);

    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
        vPoint->setEstimate(Converter::toVector3d(pMP->GetWorldPos()));
        int id = pMP->mnId+maxKFid+1;
        vPoint->setId(id);
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

        const map<KeyFrame*,size_t> observations = pMP->GetObservations();

        //Set edges
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(!pKFi->isBad())
            {                
                const cv::KeyPoint &kpUn = pKFi->mvKeysUn[mit->second];

                // Monocular observation
                if(pKFi->mvuRight[mit->second]<0)
                {
                    Eigen::Matrix<double,2,1> obs;
                    obs << kpUn.pt.x, kpUn.pt.y;

                    g2o::EdgeSE3ProjectXYZ* e = new g2o::EdgeSE3ProjectXYZ();

                    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                    e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
                    e->setMeasurement(obs);
                    const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                    e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuberMono);

                    e->fx = pKFi->fx;
                    e->fy = pKFi->fy;
                    e->cx = pKFi->cx;
                    e->cy = pKFi->cy;

                    optimizer.addEdge(e);
                    vpEdgesMono.push_back(e);
                    vpEdgeKFMono.push_back(pKFi);
                    vpMapPointEdgeMono.push_back(pMP);
                }
                else // Stereo observation
                {
                    Eigen::Matrix<double,3,1> obs;
                    const float kp_ur = pKFi->mvuRight[mit->second];
                    obs << kpUn.pt.x, kpUn.pt.y, kp_ur;

                    g2o::EdgeStereoSE3ProjectXYZ* e = new g2o::EdgeStereoSE3ProjectXYZ();

                    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                    e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
                    e->setMeasurement(obs);
                    const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                    Eigen::Matrix3d Info = Eigen::Matrix3d::Identity()*invSigma2;
                    e->setInformation(Info);

                    g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuberStereo);

                    e->fx = pKFi->fx;
                    e->fy = pKFi->fy;
                    e->cx = pKFi->cx;
                    e->cy = pKFi->cy;
                    e->bf = pKFi->mbf;

                    optimizer.addEdge(e);
                    vpEdgesStereo.push_back(e);
                    vpEdgeKFStereo.push_back(pKFi);
                    vpMapPointEdgeStereo.push_back(pMP);
                }
            }
        }
    }

    if(pbStopFlag)
        if(*pbStopFlag)
            return;

    optimizer.initializeOptimization();
    optimizer.optimize(5);

    bool bDoMore= true;

    if(pbStopFlag)
        if(*pbStopFlag)
            bDoMore = false;

    if(bDoMore)
    {

    // Check inlier observations
    for(size_t i=0, iend=vpEdgesMono.size(); i<iend;i++)
    {
        g2o::EdgeSE3ProjectXYZ* e = vpEdgesMono[i];
        MapPoint* pMP = vpMapPointEdgeMono[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>5.991 || !e->isDepthPositive())
        {
            e->setLevel(1);
        }

        e->setRobustKernel(0);
    }

    for(size_t i=0, iend=vpEdgesStereo.size(); i<iend;i++)
    {
        g2o::EdgeStereoSE3ProjectXYZ* e = vpEdgesStereo[i];
        MapPoint* pMP = vpMapPointEdgeStereo[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>7.815 || !e->isDepthPositive())
        {
            e->setLevel(1);
        }

        e->setRobustKernel(0);
    }

    // Optimize again without the outliers

    optimizer.initializeOptimization(0);
    optimizer.optimize(10);

    }

    vector<pair<KeyFrame*,MapPoint*> > vToErase;
    vToErase.reserve(vpEdgesMono.size()+vpEdgesStereo.size());

    // Check inlier observations       
    for(size_t i=0, iend=vpEdgesMono.size(); i<iend;i++)
    {
        g2o::EdgeSE3ProjectXYZ* e = vpEdgesMono[i];
        MapPoint* pMP = vpMapPointEdgeMono[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>5.991 || !e->isDepthPositive())
        {
            KeyFrame* pKFi = vpEdgeKFMono[i];
            vToErase.push_back(make_pair(pKFi,pMP));
        }
    }

    for(size_t i=0, iend=vpEdgesStereo.size(); i<iend;i++)
    {
        g2o::EdgeStereoSE3ProjectXYZ* e = vpEdgesStereo[i];
        MapPoint* pMP = vpMapPointEdgeStereo[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>7.815 || !e->isDepthPositive())
        {
            KeyFrame* pKFi = vpEdgeKFStereo[i];
            vToErase.push_back(make_pair(pKFi,pMP));
        }
    }

    // Get Map Mutex
    unique_lock<mutex> lock(pMap->mMutexMapUpdate);

    if(!vToErase.empty())
    {
        for(size_t i=0;i<vToErase.size();i++)
        {
            KeyFrame* pKFi = vToErase[i].first;
            MapPoint* pMPi = vToErase[i].second;
            pKFi->EraseMapPointMatch(pMPi);
            pMPi->EraseObservation(pKFi);
        }
    }

    // Recover optimized data

    //Keyframes
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKF = *lit;
        g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
        pKF->SetPose(Converter::toCvMat(SE3quat));
    }

    //Points
    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexSBAPointXYZ* vPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(pMP->mnId+maxKFid+1));
        pMP->SetWorldPos(Converter::toCvMat(vPoint->estimate()));
        pMP->UpdateNormalAndDepth();
    }
}


void Optimizer::OptimizeEssentialGraph(Map* pMap, KeyFrame* pLoopKF, KeyFrame* pCurKF,
                                       const LoopClosing::KeyFrameAndPose &NonCorrectedSim3,
                                       const LoopClosing::KeyFrameAndPose &CorrectedSim3,
                                       const map<KeyFrame *, set<KeyFrame *> > &LoopConnections, const bool &bFixScale)
{
    // Setup optimizer
    g2o::SparseOptimizer optimizer;
    optimizer.setVerbose(false);
    g2o::BlockSolver_7_3::LinearSolverType * linearSolver =
           new g2o::LinearSolverEigen<g2o::BlockSolver_7_3::PoseMatrixType>();
    g2o::BlockSolver_7_3 * solver_ptr= new g2o::BlockSolver_7_3(linearSolver);
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);

    solver->setUserLambdaInit(1e-16);
    optimizer.setAlgorithm(solver);

    const vector<KeyFrame*> vpKFs = pMap->GetAllKeyFrames();
    const vector<MapPoint*> vpMPs = pMap->GetAllMapPoints();

    const unsigned int nMaxKFid = pMap->GetMaxKFid();

    vector<g2o::Sim3,Eigen::aligned_allocator<g2o::Sim3> > vScw(nMaxKFid+1);
    vector<g2o::Sim3,Eigen::aligned_allocator<g2o::Sim3> > vCorrectedSwc(nMaxKFid+1);
    vector<g2o::VertexSim3Expmap*> vpVertices(nMaxKFid+1);

    const int minFeat = 100;

    // Set KeyFrame vertices
    for(size_t i=0, iend=vpKFs.size(); i<iend;i++)
    {
        KeyFrame* pKF = vpKFs[i];
        if(pKF->isBad())
            continue;
        g2o::VertexSim3Expmap* VSim3 = new g2o::VertexSim3Expmap();

        const int nIDi = pKF->mnId;

        LoopClosing::KeyFrameAndPose::const_iterator it = CorrectedSim3.find(pKF);

        if(it!=CorrectedSim3.end())
        {
            vScw[nIDi] = it->second;
            VSim3->setEstimate(it->second);
        }
        else
        {
            Eigen::Matrix<double,3,3> Rcw = Converter::toMatrix3d(pKF->GetRotation());
            Eigen::Matrix<double,3,1> tcw = Converter::toVector3d(pKF->GetTranslation());
            g2o::Sim3 Siw(Rcw,tcw,1.0);
            vScw[nIDi] = Siw;
            VSim3->setEstimate(Siw);
        }

        if(pKF==pLoopKF)
            VSim3->setFixed(true);

        VSim3->setId(nIDi);
        VSim3->setMarginalized(false);
        VSim3->_fix_scale = bFixScale;

        optimizer.addVertex(VSim3);

        vpVertices[nIDi]=VSim3;
    }


    set<pair<long unsigned int,long unsigned int> > sInsertedEdges;

    const Eigen::Matrix<double,7,7> matLambda = Eigen::Matrix<double,7,7>::Identity();

    // Set Loop edges
    for(map<KeyFrame *, set<KeyFrame *> >::const_iterator mit = LoopConnections.begin(), mend=LoopConnections.end(); mit!=mend; mit++)
    {
        KeyFrame* pKF = mit->first;
        const long unsigned int nIDi = pKF->mnId;
        const set<KeyFrame*> &spConnections = mit->second;
        const g2o::Sim3 Siw = vScw[nIDi];
        const g2o::Sim3 Swi = Siw.inverse();

        for(set<KeyFrame*>::const_iterator sit=spConnections.begin(), send=spConnections.end(); sit!=send; sit++)
        {
            const long unsigned int nIDj = (*sit)->mnId;
            if((nIDi!=pCurKF->mnId || nIDj!=pLoopKF->mnId) && pKF->GetWeight(*sit)<minFeat)
                continue;

            const g2o::Sim3 Sjw = vScw[nIDj];
            const g2o::Sim3 Sji = Sjw * Swi;

            g2o::EdgeSim3* e = new g2o::EdgeSim3();
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDj)));
            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDi)));
            e->setMeasurement(Sji);

            e->information() = matLambda;

            optimizer.addEdge(e);

            sInsertedEdges.insert(make_pair(min(nIDi,nIDj),max(nIDi,nIDj)));
        }
    }

    // Set normal edges
    for(size_t i=0, iend=vpKFs.size(); i<iend; i++)
    {
        KeyFrame* pKF = vpKFs[i];

        const int nIDi = pKF->mnId;

        g2o::Sim3 Swi;

        LoopClosing::KeyFrameAndPose::const_iterator iti = NonCorrectedSim3.find(pKF);

        if(iti!=NonCorrectedSim3.end())
            Swi = (iti->second).inverse();
        else
            Swi = vScw[nIDi].inverse();

        KeyFrame* pParentKF = pKF->GetParent();

        // Spanning tree edge
        if(pParentKF)
        {
            int nIDj = pParentKF->mnId;

            g2o::Sim3 Sjw;

            LoopClosing::KeyFrameAndPose::const_iterator itj = NonCorrectedSim3.find(pParentKF);

            if(itj!=NonCorrectedSim3.end())
                Sjw = itj->second;
            else
                Sjw = vScw[nIDj];

            g2o::Sim3 Sji = Sjw * Swi;

            g2o::EdgeSim3* e = new g2o::EdgeSim3();
            e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDj)));
            e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDi)));
            e->setMeasurement(Sji);

            e->information() = matLambda;
            optimizer.addEdge(e);
        }

        // Loop edges
        const set<KeyFrame*> sLoopEdges = pKF->GetLoopEdges();
        for(set<KeyFrame*>::const_iterator sit=sLoopEdges.begin(), send=sLoopEdges.end(); sit!=send; sit++)
        {
            KeyFrame* pLKF = *sit;
            if(pLKF->mnId<pKF->mnId)
            {
                g2o::Sim3 Slw;

                LoopClosing::KeyFrameAndPose::const_iterator itl = NonCorrectedSim3.find(pLKF);

                if(itl!=NonCorrectedSim3.end())
                    Slw = itl->second;
                else
                    Slw = vScw[pLKF->mnId];

                g2o::Sim3 Sli = Slw * Swi;
                g2o::EdgeSim3* el = new g2o::EdgeSim3();
                el->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pLKF->mnId)));
                el->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDi)));
                el->setMeasurement(Sli);
                el->information() = matLambda;
                optimizer.addEdge(el);
            }
        }

        // Covisibility graph edges
        const vector<KeyFrame*> vpConnectedKFs = pKF->GetCovisiblesByWeight(minFeat);
        for(vector<KeyFrame*>::const_iterator vit=vpConnectedKFs.begin(); vit!=vpConnectedKFs.end(); vit++)
        {
            KeyFrame* pKFn = *vit;
            if(pKFn && pKFn!=pParentKF && !pKF->hasChild(pKFn) && !sLoopEdges.count(pKFn))
            {
                if(!pKFn->isBad() && pKFn->mnId<pKF->mnId)
                {
                    if(sInsertedEdges.count(make_pair(min(pKF->mnId,pKFn->mnId),max(pKF->mnId,pKFn->mnId))))
                        continue;

                    g2o::Sim3 Snw;

                    LoopClosing::KeyFrameAndPose::const_iterator itn = NonCorrectedSim3.find(pKFn);

                    if(itn!=NonCorrectedSim3.end())
                        Snw = itn->second;
                    else
                        Snw = vScw[pKFn->mnId];

                    g2o::Sim3 Sni = Snw * Swi;

                    g2o::EdgeSim3* en = new g2o::EdgeSim3();
                    en->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFn->mnId)));
                    en->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(nIDi)));
                    en->setMeasurement(Sni);
                    en->information() = matLambda;
                    optimizer.addEdge(en);
                }
            }
        }
    }

    // Optimize!
    optimizer.initializeOptimization();
    optimizer.optimize(20);

    unique_lock<mutex> lock(pMap->mMutexMapUpdate);

    // SE3 Pose Recovering. Sim3:[sR t;0 1] -> SE3:[R t/s;0 1]
    for(size_t i=0;i<vpKFs.size();i++)
    {
        KeyFrame* pKFi = vpKFs[i];

        const int nIDi = pKFi->mnId;

        g2o::VertexSim3Expmap* VSim3 = static_cast<g2o::VertexSim3Expmap*>(optimizer.vertex(nIDi));
        g2o::Sim3 CorrectedSiw =  VSim3->estimate();
        vCorrectedSwc[nIDi]=CorrectedSiw.inverse();
        Eigen::Matrix3d eigR = CorrectedSiw.rotation().toRotationMatrix();
        Eigen::Vector3d eigt = CorrectedSiw.translation();
        double s = CorrectedSiw.scale();

        eigt *=(1./s); //[R t/s;0 1]

        cv::Mat Tiw = Converter::toCvSE3(eigR,eigt);

        pKFi->SetPose(Tiw);
    }

    // Correct points. Transform to "non-optimized" reference keyframe pose and transform back with optimized pose
    for(size_t i=0, iend=vpMPs.size(); i<iend; i++)
    {
        MapPoint* pMP = vpMPs[i];

        if(pMP->isBad())
            continue;

        int nIDr;
        if(pMP->mnCorrectedByKF==pCurKF->mnId)
        {
            nIDr = pMP->mnCorrectedReference;
        }
        else
        {
            KeyFrame* pRefKF = pMP->GetReferenceKeyFrame();
            nIDr = pRefKF->mnId;
        }


        g2o::Sim3 Srw = vScw[nIDr];
        g2o::Sim3 correctedSwr = vCorrectedSwc[nIDr];

        cv::Mat P3Dw = pMP->GetWorldPos();
        Eigen::Matrix<double,3,1> eigP3Dw = Converter::toVector3d(P3Dw);
        Eigen::Matrix<double,3,1> eigCorrectedP3Dw = correctedSwr.map(Srw.map(eigP3Dw));

        cv::Mat cvCorrectedP3Dw = Converter::toCvMat(eigCorrectedP3Dw);
        pMP->SetWorldPos(cvCorrectedP3Dw);

        pMP->UpdateNormalAndDepth();
    }
}

int Optimizer::OptimizeSim3(KeyFrame *pKF1, KeyFrame *pKF2, vector<MapPoint *> &vpMatches1, g2o::Sim3 &g2oS12, const float th2, const bool bFixScale)
{
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolverX::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>();

    g2o::BlockSolverX * solver_ptr = new g2o::BlockSolverX(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    // Calibration
    const cv::Mat &K1 = pKF1->mK;
    const cv::Mat &K2 = pKF2->mK;

    // Camera poses
    const cv::Mat R1w = pKF1->GetRotation();
    const cv::Mat t1w = pKF1->GetTranslation();
    const cv::Mat R2w = pKF2->GetRotation();
    const cv::Mat t2w = pKF2->GetTranslation();

    // Set Sim3 vertex
    g2o::VertexSim3Expmap * vSim3 = new g2o::VertexSim3Expmap();    
    vSim3->_fix_scale=bFixScale;
    vSim3->setEstimate(g2oS12);
    vSim3->setId(0);
    vSim3->setFixed(false);
    vSim3->_principle_point1[0] = K1.at<float>(0,2);
    vSim3->_principle_point1[1] = K1.at<float>(1,2);
    vSim3->_focal_length1[0] = K1.at<float>(0,0);
    vSim3->_focal_length1[1] = K1.at<float>(1,1);
    vSim3->_principle_point2[0] = K2.at<float>(0,2);
    vSim3->_principle_point2[1] = K2.at<float>(1,2);
    vSim3->_focal_length2[0] = K2.at<float>(0,0);
    vSim3->_focal_length2[1] = K2.at<float>(1,1);
    optimizer.addVertex(vSim3);

    // Set MapPoint vertices
    const int N = vpMatches1.size();
    const vector<MapPoint*> vpMapPoints1 = pKF1->GetMapPointMatches();
    vector<g2o::EdgeSim3ProjectXYZ*> vpEdges12;
    vector<g2o::EdgeInverseSim3ProjectXYZ*> vpEdges21;
    vector<size_t> vnIndexEdge;

    vnIndexEdge.reserve(2*N);
    vpEdges12.reserve(2*N);
    vpEdges21.reserve(2*N);

    const float deltaHuber = sqrt(th2);

    int nCorrespondences = 0;

    for(int i=0; i<N; i++)
    {
        if(!vpMatches1[i])
            continue;

        MapPoint* pMP1 = vpMapPoints1[i];
        MapPoint* pMP2 = vpMatches1[i];

        const int id1 = 2*i+1;
        const int id2 = 2*(i+1);

        const int i2 = pMP2->GetIndexInKeyFrame(pKF2);

        if(pMP1 && pMP2)
        {
            if(!pMP1->isBad() && !pMP2->isBad() && i2>=0)
            {
                g2o::VertexSBAPointXYZ* vPoint1 = new g2o::VertexSBAPointXYZ();
                cv::Mat P3D1w = pMP1->GetWorldPos();
                cv::Mat P3D1c = R1w*P3D1w + t1w;
                vPoint1->setEstimate(Converter::toVector3d(P3D1c));
                vPoint1->setId(id1);
                vPoint1->setFixed(true);
                optimizer.addVertex(vPoint1);

                g2o::VertexSBAPointXYZ* vPoint2 = new g2o::VertexSBAPointXYZ();
                cv::Mat P3D2w = pMP2->GetWorldPos();
                cv::Mat P3D2c = R2w*P3D2w + t2w;
                vPoint2->setEstimate(Converter::toVector3d(P3D2c));
                vPoint2->setId(id2);
                vPoint2->setFixed(true);
                optimizer.addVertex(vPoint2);
            }
            else
                continue;
        }
        else
            continue;

        nCorrespondences++;

        // Set edge x1 = S12*X2
        Eigen::Matrix<double,2,1> obs1;
        const cv::KeyPoint &kpUn1 = pKF1->mvKeysUn[i];
        obs1 << kpUn1.pt.x, kpUn1.pt.y;

        g2o::EdgeSim3ProjectXYZ* e12 = new g2o::EdgeSim3ProjectXYZ();
        e12->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id2)));
        e12->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
        e12->setMeasurement(obs1);
        const float &invSigmaSquare1 = pKF1->mvInvLevelSigma2[kpUn1.octave];
        e12->setInformation(Eigen::Matrix2d::Identity()*invSigmaSquare1);

        g2o::RobustKernelHuber* rk1 = new g2o::RobustKernelHuber;
        e12->setRobustKernel(rk1);
        rk1->setDelta(deltaHuber);
        optimizer.addEdge(e12);

        // Set edge x2 = S21*X1
        Eigen::Matrix<double,2,1> obs2;
        const cv::KeyPoint &kpUn2 = pKF2->mvKeysUn[i2];
        obs2 << kpUn2.pt.x, kpUn2.pt.y;

        g2o::EdgeInverseSim3ProjectXYZ* e21 = new g2o::EdgeInverseSim3ProjectXYZ();

        e21->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id1)));
        e21->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
        e21->setMeasurement(obs2);
        float invSigmaSquare2 = pKF2->mvInvLevelSigma2[kpUn2.octave];
        e21->setInformation(Eigen::Matrix2d::Identity()*invSigmaSquare2);

        g2o::RobustKernelHuber* rk2 = new g2o::RobustKernelHuber;
        e21->setRobustKernel(rk2);
        rk2->setDelta(deltaHuber);
        optimizer.addEdge(e21);

        vpEdges12.push_back(e12);
        vpEdges21.push_back(e21);
        vnIndexEdge.push_back(i);
    }

    // Optimize!
    optimizer.initializeOptimization();
    optimizer.optimize(5);

    // Check inliers
    int nBad=0;
    for(size_t i=0; i<vpEdges12.size();i++)
    {
        g2o::EdgeSim3ProjectXYZ* e12 = vpEdges12[i];
        g2o::EdgeInverseSim3ProjectXYZ* e21 = vpEdges21[i];
        if(!e12 || !e21)
            continue;

        if(e12->chi2()>th2 || e21->chi2()>th2)
        {
            size_t idx = vnIndexEdge[i];
            vpMatches1[idx]=static_cast<MapPoint*>(NULL);
            optimizer.removeEdge(e12);
            optimizer.removeEdge(e21);
            vpEdges12[i]=static_cast<g2o::EdgeSim3ProjectXYZ*>(NULL);
            vpEdges21[i]=static_cast<g2o::EdgeInverseSim3ProjectXYZ*>(NULL);
            nBad++;
        }
    }

    int nMoreIterations;
    if(nBad>0)
        nMoreIterations=10;
    else
        nMoreIterations=5;

    if(nCorrespondences-nBad<10)
        return 0;

    // Optimize again only with inliers

    optimizer.initializeOptimization();
    optimizer.optimize(nMoreIterations);

    int nIn = 0;
    for(size_t i=0; i<vpEdges12.size();i++)
    {
        g2o::EdgeSim3ProjectXYZ* e12 = vpEdges12[i];
        g2o::EdgeInverseSim3ProjectXYZ* e21 = vpEdges21[i];
        if(!e12 || !e21)
            continue;

        if(e12->chi2()>th2 || e21->chi2()>th2)
        {
            size_t idx = vnIndexEdge[i];
            vpMatches1[idx]=static_cast<MapPoint*>(NULL);
        }
        else
            nIn++;
    }

    // Recover optimized Sim3
    g2o::VertexSim3Expmap* vSim3_recov = static_cast<g2o::VertexSim3Expmap*>(optimizer.vertex(0));
    g2oS12= vSim3_recov->estimate();

    return nIn;
}

//复合点线特征位姿优化函数：（前端）
void Optimizer::PoseOptimizationmain(Frame* pFrame, int &inlinersNum, int &inlinerLinesNum, int nIterations )
{	
    float PoseErrorForPoints;
    int nmatches;
    
    PoseOptimizationPoints(pFrame, PoseErrorForPoints, nmatches, 10);

    pFrame->SetPose(pFrame->mTcwPoints);

    PoseOptimizationBoth(pFrame,inlinersNum,inlinerLinesNum,nIterations);
    
}

//b.两次点优化：
void Optimizer::PoseOptimizationDoublePoints(Frame* pFrame, int &inlinersNum )
{
    float PoseErrorForPoints;
    int nmatches;
    
    PoseOptimizationPoints(pFrame,PoseErrorForPoints,nmatches,10);
    
    pFrame->SetPose(pFrame->mTcwPoints);
    
    PoseOptimizationPoints(pFrame,PoseErrorForPoints,inlinersNum,10);
    
    pFrame->SetPose(pFrame->mTcwPoints);
}

//c.两次线优化：
void Optimizer::PoseOptimizationLowFeature(Frame* pFrame, int &inlinersNum, int &inlinerLinesNum, int nIterations )
{
    PoseOptimizationBoth(pFrame,inlinersNum,inlinerLinesNum,nIterations);
    
    PoseOptimizationBoth(pFrame,inlinersNum,inlinerLinesNum,nIterations);
    
}
//复合优化子函数：
void Optimizer::PoseOptimizationPoints(Frame* pFrame, float &PoseErrorForPoints, int &nmatches, int nIterations )
{
    //初始化G2O优化器
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    int nInitialCorrespondences=0;

    // Set Frame vertex
    g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
    vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
    vSE3->setId(0);
    vSE3->setFixed(false);
    optimizer.addVertex(vSE3);

    // Set MapPoint vertices
    const int N = pFrame->N;

    vector<g2o::EdgeSE3ProjectXYZOnlyPose*> vpEdgesMono; 
    vector<size_t> vnIndexEdgeMono;
    vpEdgesMono.reserve(N);
    vnIndexEdgeMono.reserve(N);

    const float deltaMono = sqrt(5.991);

    {
    unique_lock<mutex> lock(MapPoint::mGlobalMutex);

    for(int i=0; i<N; i++)
    {
        MapPoint* pMP = pFrame->mvpMapPoints[i];
        if(pMP)
        {
            // Monocular observation
            if(pFrame->mvuRight[i]<0)
            {
                nInitialCorrespondences++;

                Eigen::Matrix<double,2,1> obs;
                const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
                obs << kpUn.pt.x, kpUn.pt.y;

                g2o::EdgeSE3ProjectXYZOnlyPose* e = new g2o::EdgeSE3ProjectXYZOnlyPose();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
                e->setMeasurement(obs);
                const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
                e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                g2o::RobustKernelCauchy* rk = new g2o::RobustKernelCauchy;
                e->setRobustKernel(rk);
                rk->setDelta(deltaMono);

                e->fx = pFrame->fx;
                e->fy = pFrame->fy;
                e->cx = pFrame->cx;
                e->cy = pFrame->cy;
                cv::Mat Xw = pMP->GetWorldPos();
                e->Xw[0] = Xw.at<float>(0);
                e->Xw[1] = Xw.at<float>(1);
                e->Xw[2] = Xw.at<float>(2);
		
		//新添加
		if(pFrame->mvbOutlier[i])
		{
		    e->setLevel(1);
		}

                optimizer.addEdge(e);
		
                vpEdgesMono.push_back(e);
                vnIndexEdgeMono.push_back(i);
            }
	}
    }    
    }
    
    optimizer.initializeOptimization(0);
    optimizer.optimize(nIterations);
    
    int nBad=0;
    int ninliersNum = 0;
    PoseErrorForPoints = 0.0;
    
    for(size_t i=0, iend=vpEdgesMono.size(); i<iend; i++)
    {
         g2o::EdgeSE3ProjectXYZOnlyPose* e = vpEdgesMono[i];

         const size_t idx = vnIndexEdgeMono[i];

         if(pFrame->mvbOutlier[idx])
         {
              e->computeError();
         }

         const float chi2 = e->chi2();

         if(chi2>5.991)
         {                
             pFrame->mvbOutlier[idx]=true;
             nBad++;
         }
         else
         {
             pFrame->mvbOutlier[idx]=false;  
	     PoseErrorForPoints += chi2;
	     ninliersNum ++;
         }
    }
    PoseErrorForPoints = PoseErrorForPoints / ninliersNum;
    
    // Recover optimized pose and return number of inliers
    g2o::VertexSE3Expmap* vSE3_recov = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
    g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
    cv::Mat pose = Converter::toCvMat(SE3quat_recov);
    pFrame->SetPosePoints(pose);
    
    nmatches = nInitialCorrespondences - nBad;
}

//b.线特征优化（g2o）
void Optimizer::PoseOptimizationLines(Frame* pFrame, float &PoseErrorForLines, int &nmatchLines, int nIterations )
{
    //初始化G2O优化器
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    int nInitialCorrespondences=0;

    // Set Frame vertex
    g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
    vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
    vSE3->setId(0);
    vSE3->setFixed(false);
    optimizer.addVertex(vSE3);

    // Set MapLine vertices
    const int NL = pFrame->NL;

    vector<EdgeSE3ProjectXYZOnlyPoseLines*> vpEdgesMono; 
    vector<size_t> vnIndexEdgeMono;
    vpEdgesMono.reserve(NL);
    vnIndexEdgeMono.reserve(NL);

    const float deltaMono = sqrt(3.841);

    {
    unique_lock<mutex> lock(MapLine::mGlobalMutex);

    for(int i=0; i<NL; i++)
    {
        MapLine* pML = pFrame->mvpMapLines[i];
        if(pML)
        {
	    nInitialCorrespondences++;

	    Eigen::Matrix<double,3,1> obs;
	    const KeyLine &kL = pFrame->mvLinesUn[i];
	    Vector3f sp_l; sp_l << kL.startPointX, kL.startPointY, 1.0;
	    Vector3f ep_l; ep_l << kL.endPointX,   kL.endPointY,   1.0;
	    Vector3f le_l; le_l << sp_l.cross(ep_l);
	    //le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) + le_l(2)*le_l(2));
	    le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1));
	    obs << le_l(0), le_l(1), le_l(2);
	    
	    EdgeSE3ProjectXYZOnlyPoseLines* e = new EdgeSE3ProjectXYZOnlyPoseLines();

	    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
	    e->setMeasurement(obs);
	    const float invSigma2 = pFrame->mvInvLevelSigma2Lines[kL.octave];
	    e->setInformation(Eigen::Matrix3d::Identity()*invSigma2);

	    g2o::RobustKernelCauchy* rk = new g2o::RobustKernelCauchy;
	    e->setRobustKernel(rk);
	    rk->setDelta(deltaMono);

	    e->fx = pFrame->fx;
	    e->fy = pFrame->fy;
	    e->cx = pFrame->cx;
	    e->cy = pFrame->cy;
	    cv::Mat Xw = pML->GetMidPWorldPos();
	    e->Xw[0] = Xw.at<float>(0);
	    e->Xw[1] = Xw.at<float>(1);
	    e->Xw[2] = Xw.at<float>(2);
	    
	    //新添加
	    if(pFrame->mvbOutlierLines[i])
	    {
		e->setLevel(1);
	    }

	    optimizer.addEdge(e);    
	    vpEdgesMono.push_back(e);
	    vnIndexEdgeMono.push_back(i);      
	}
    }   
    }

    
    optimizer.initializeOptimization(0);
    optimizer.optimize(nIterations);
    
    int nBad=0;
    int ninliersNum = 0;
    PoseErrorForLines = 0.0;

    for(size_t i=0, iend=vpEdgesMono.size(); i<iend; i++)
    {
         EdgeSE3ProjectXYZOnlyPoseLines* e = vpEdgesMono[i];

         const size_t idx = vnIndexEdgeMono[i];

         if(pFrame->mvbOutlierLines[idx])
         {
              e->computeError();
         }

         const float chi2 = e->chi2();

         if(chi2>3.841)
         {                
             pFrame->mvbOutlierLines[idx]=true;
             nBad++;
         }
         else
         {
             pFrame->mvbOutlierLines[idx]=false;  
	     PoseErrorForLines += chi2;
	     ninliersNum ++;
         }
    }
    PoseErrorForLines = PoseErrorForLines / ninliersNum;
    
    // Recover optimized pose and return number of inliers
    g2o::VertexSE3Expmap* vSE3_recov = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
    g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
    cv::Mat pose = Converter::toCvMat(SE3quat_recov);
    pFrame->SetPoseLines(pose);
    
    nmatchLines = nInitialCorrespondences - nBad;
}

void Optimizer::PoseOptimizationBothOld(Frame* pFrame, int &inlinersNum, int &inlinerLinesNum, int nIterations )
{
  
    GaussNewtonOptimizationForPose(pFrame,nIterations);
    
    RecheckBothOutliersForPose(pFrame,inlinersNum,inlinerLinesNum);
    
}

// 调用此函数（g2o），进行鲁棒的LM优化，分别返回内点和内线的数量。
void Optimizer::PoseOptimizationBoth(Frame* pFrame, int &inlinersNum, int &inlinerLinesNum, int nIterations )
{
    //初始化G2O优化器
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    int nInitialCorrespondences=0;
    int nInitialCorrespondenceLines=0;

    // Set Frame vertex
    g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
    vSE3->setEstimate(Converter::toSE3Quat(pFrame->mTcw));
    vSE3->setId(0);
    vSE3->setFixed(false);
    optimizer.addVertex(vSE3);
    
    // Set MapPoint vertices
    const int N = pFrame->N;

    vector<g2o::EdgeSE3ProjectXYZOnlyPose*> vpEdgesMono; 
    vector<size_t> vnIndexEdgeMono;
    vpEdgesMono.reserve(N);
    vnIndexEdgeMono.reserve(N);

    const float deltaMono = sqrt(5.991);
      
    // Set MapLine vertices
    const int NL = pFrame->NL;

    vector<EdgeSE3ProjectXYZOnlyPoseLines*> vpEdgesMonoLines;
    vector<size_t> vnIndexEdgeMonoLines;
    vpEdgesMonoLines.reserve(NL);
    vnIndexEdgeMonoLines.reserve(NL);

    const float deltaMonoLines = sqrt(3.841);

    {
    unique_lock<mutex> lock(MapPoint::mGlobalMutex);

    for(int i=0; i<N; i++)
    {
        MapPoint* pMP = pFrame->mvpMapPoints[i];
        if(pMP)
        {
            // Monocular observation
            if(pFrame->mvuRight[i]<0)
            {
                nInitialCorrespondences++;

                Eigen::Matrix<double,2,1> obs;
                const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
                obs << kpUn.pt.x, kpUn.pt.y;

                g2o::EdgeSE3ProjectXYZOnlyPose* e = new g2o::EdgeSE3ProjectXYZOnlyPose();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
                e->setMeasurement(obs);
                const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
                e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                g2o::RobustKernelCauchy* rk = new g2o::RobustKernelCauchy;
                e->setRobustKernel(rk);
                rk->setDelta(deltaMono);

                e->fx = pFrame->fx;
                e->fy = pFrame->fy;
                e->cx = pFrame->cx;
                e->cy = pFrame->cy;
                cv::Mat Xw = pMP->GetWorldPos();
                e->Xw[0] = Xw.at<float>(0);
                e->Xw[1] = Xw.at<float>(1);
                e->Xw[2] = Xw.at<float>(2);

		if(pFrame->mvbOutlier[i])
		{
		    e->setLevel(1);
		}

                optimizer.addEdge(e);
		
                vpEdgesMono.push_back(e);
                vnIndexEdgeMono.push_back(i);
            }
	}
    }    
    }
    
    {
    unique_lock<mutex> lock(MapLine::mGlobalMutex);

    for(int i=0; i<NL; i++)
    {
        MapLine* pML = pFrame->mvpMapLines[i];
        if(pML)
        {
	    nInitialCorrespondenceLines++;

	    Eigen::Matrix<double,3,1> obs;
	    const KeyLine &kL = pFrame->mvLinesUn[i];
	    Vector3f sp_l; sp_l << kL.startPointX, kL.startPointY, 1.0;
	    Vector3f ep_l; ep_l << kL.endPointX,   kL.endPointY,   1.0;
	    Vector3f le_l; le_l << sp_l.cross(ep_l);
	    //le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) + le_l(2)*le_l(2));
	    le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1));
	    obs << le_l(0), le_l(1), le_l(2);
	    
	    EdgeSE3ProjectXYZOnlyPoseLines* e = new EdgeSE3ProjectXYZOnlyPoseLines();

	    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
	    e->setMeasurement(obs);
	    const float invSigma2 = pFrame->mvInvLevelSigma2Lines[kL.octave];
	    e->setInformation(Eigen::Matrix3d::Identity()*invSigma2);

	    g2o::RobustKernelCauchy* rk = new g2o::RobustKernelCauchy;
	    e->setRobustKernel(rk);
	    rk->setDelta(deltaMonoLines);

	    e->fx = pFrame->fx;
	    e->fy = pFrame->fy;
	    e->cx = pFrame->cx;
	    e->cy = pFrame->cy;
	    cv::Mat Xw = pML->GetMidPWorldPos();
	    e->Xw[0] = Xw.at<float>(0);
	    e->Xw[1] = Xw.at<float>(1);
	    e->Xw[2] = Xw.at<float>(2);
	    
	    //新添加
	    if(pFrame->mvbOutlierLines[i])
	    {
		e->setLevel(1);
	    }

	    optimizer.addEdge(e);
	    
	    vpEdgesMonoLines.push_back(e);
	    vnIndexEdgeMonoLines.push_back(i);      
	}
    }
    }
    
    optimizer.initializeOptimization(0);
    optimizer.optimize(nIterations);
    
    int nBad = 0;
    int nBadLines = 0;

    for(size_t i=0, iend=vpEdgesMono.size(); i<iend; i++)
    {
         g2o::EdgeSE3ProjectXYZOnlyPose* e = vpEdgesMono[i];

         const size_t idx = vnIndexEdgeMono[i];

         if(pFrame->mvbOutlier[idx])
         {
              e->computeError();
         }

         const float chi2 = e->chi2();

         if(chi2>5.991)
         {                
             pFrame->mvbOutlier[idx]=true;
             nBad++;
         }
         else
         {
             pFrame->mvbOutlier[idx]=false;  
         }
    }

    for(size_t i=0, iend=vpEdgesMonoLines.size(); i<iend; i++)
    {
         EdgeSE3ProjectXYZOnlyPoseLines* e = vpEdgesMonoLines[i];

         const size_t idx = vnIndexEdgeMonoLines[i];

         if(pFrame->mvbOutlierLines[idx])
         {
              e->computeError();
         }

         const float chi2 = e->chi2();

         if(chi2>3.841)
         {                
             pFrame->mvbOutlierLines[idx]=true;
             nBadLines++;
         }
         else
         {
             pFrame->mvbOutlierLines[idx]=false;  
         }
    }
    
    // Recover optimized pose and return number of inliers
    g2o::VertexSE3Expmap* vSE3_recov = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
    g2o::SE3Quat SE3quat_recov = vSE3_recov->estimate();
    cv::Mat pose = Converter::toCvMat(SE3quat_recov);
    pFrame->SetPose(pose);
    
    inlinersNum = nInitialCorrespondences - nBad;
    inlinerLinesNum = nInitialCorrespondenceLines - nBadLines;
    
}

//a.前端：
void Optimizer::GaussNewtonOptimizationForPose(Frame* pFrame, int nIterations)
{
      Matrix6d H;
      Vector6d g, Tcw_inc;
      double Err, Err_prev = 999999999.9;
      bool Solution_is_good = true;
      g2o::SE3Quat TcwBothUpdateG2o = Converter::toSE3Quat(pFrame->mTcw);
      
      for(int i=0; i<nIterations; i++)
      {
	  //基于当前的位姿矩阵，计算当步海森矩阵和g矩阵
	  //内部使用SE3Quat.map()，计算相机坐标系下的3D点
	  OptimizeFunctionForPose(pFrame, TcwBothUpdateG2o, H, g, Err);
	  
	  //如果点线特征的单位误差很小或者误差变化量很小
	  if(Err<1e-7 || fabs(Err-Err_prev)<1e-7  )
	    break;
	  
	  //更新计算位姿增量,判断位姿变化量结果的有效性
	  ColPivHouseholderQR<Matrix6d> solver(H);
	  Tcw_inc = solver.solve(g);
	  if( solver.logAbsDeterminant() < 0.0 || solver.info() != Success )
	  {
	      Solution_is_good = false;
	      break;
	  }
	  
	  //位姿矩阵更新
	  TcwBothUpdateG2o = g2o::SE3Quat::exp(Tcw_inc)*TcwBothUpdateG2o;
	  
	  //更新上一步的单位误差
	  Err_prev = Err;
      }
      
      //根据结果是否正确决定存储到当前图像帧中mat格式位姿矩阵的信息
      if(Solution_is_good)
      {
	  //存储优化后的位姿信息
	  cv::Mat pose = Converter::toCvMat(TcwBothUpdateG2o);
	  pFrame->SetPose(pose);
      }
      else
      {
	  //当前图像帧位姿为空，必须置空哦，表示跟踪失败～
	  pFrame->SetPose(cv::Mat());
      }
}

//计算H、g和当步总误差error：
void Optimizer::OptimizeFunctionForPose(Frame* pFrame, g2o::SE3Quat &TcwBothUpdate, Matrix6d &H, Vector6d &g, double &error)
{
    //标定矩阵
    cv::Mat K = pFrame->mK;
    const float fx = K.at<float>(0,0);
    const float fy = K.at<float>(1,1);
    const float cx = K.at<float>(0,2);
    const float cy = K.at<float>(1,2);
    
    //设置点和线特征各自的局部变量：(1)H矩阵：H_p和H_l。(2)g矩阵：g_p和g_l。(3)总误差：e_p和e_l
    Matrix6d H_l, H_p;
    Vector6d g_l, g_p;
    double   e_l = 0.0, e_p = 0.0;
    
    //设置点和线特征整合的局部变量：
    H = Matrix6d::Zero(); H_l = H; H_p = H;
    g = Vector6d::Zero(); g_l = g; g_p = g;
    error = 0.0;
    
    //设置鲁棒核函数中的单像素误差边界
    const float deltaMonoPoint = sqrt(5.991);
    const float deltaMonoLine = sqrt(3.841);
    
    int N_p = 0;
    // Set MapPoint vertices
    const int N = pFrame->N;
    {
	unique_lock<mutex> lock(MapPoint::mGlobalMutex);
	for(int i=0; i<N; i++)
	{
	    MapPoint* pMP = pFrame->mvpMapPoints[i];
	    if(pMP && !pFrame->mvbOutlier[i])
	    {
		//获取观测值
		const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
		Vector2d Point2D; Point2D << kpUn.pt.x, kpUn.pt.y;
		
		//设置标准信息矩阵
		const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
		Matrix2d Information = Eigen::Matrix2d::Identity()*invSigma2;
				
		//获取3D预测值和2D预测值
		cv::Mat Xw = pMP->GetWorldPos();
		//点特征雅克比矩阵参数
		float PX = Xw.at<float>(0);
		float PY = Xw.at<float>(1);
		float PZ = Xw.at<float>(2);
		float PZinv = 1.0f / PZ;
		float PZ2inv = 1.0f / ( PZ * PZ );
		
		//debug:这里不能用const
		Vector3d XYZ; XYZ <<  PX, PY, PZ;//地图点世界坐标
		Vector3d Trans_XYZ = TcwBothUpdate.map(XYZ);//地图点相机坐标
		Vector2d PointProjection2D;//地图点的图像坐标系坐标		
		PointProjection2D(0) = fx*Trans_XYZ(0)/Trans_XYZ(2) + cx;
		PointProjection2D(1) = fy*Trans_XYZ(1)/Trans_XYZ(2) + cy;
		Vector2d PointError2D; PointError2D << Point2D(0)-PointProjection2D(0), Point2D(1)-PointProjection2D(1);
		
		//标准信息矩阵下重投影误差的模平方
		float PointReProjectionErrorNorm2 = PointError2D.dot(Information*PointError2D);
		//基于Huber核函数计算权重
		float rho = RobustWeightHuber(PointReProjectionErrorNorm2, deltaMonoPoint);
		Matrix2d RobWInformation = rho*Information;
		
		
		//float Eu = PointError2D(0);
		//float Ev = PointError2D(1);
		
		//点特征雅克比矩阵(PL-SLAM(stereo)的这里有问题)
		Matrix<double,2,6> _jacobianOplusXi;
		_jacobianOplusXi(0,0) =  PX*PY*PZ2inv *fx;
		_jacobianOplusXi(0,1) = -(1+(PX*PX*PZ2inv)) *fx;
		_jacobianOplusXi(0,2) = PY*PZinv *fx;
		_jacobianOplusXi(0,3) = -PZinv *fx;
		_jacobianOplusXi(0,4) = 0;
		_jacobianOplusXi(0,5) = PX*PZ2inv *fx;

		_jacobianOplusXi(1,0) = (1+PY*PY*PZ2inv) *fy;
		_jacobianOplusXi(1,1) = -PX*PY*PZ2inv *fy;
		_jacobianOplusXi(1,2) = -PX*PZinv *fy;
		_jacobianOplusXi(1,3) = 0;
		_jacobianOplusXi(1,4) = -PZinv *fy;
		_jacobianOplusXi(1,5) = PY*PZ2inv *fy;
		
		//点特征的H、g、参与优化的点误差累加和、参与优化的点特征数量
		H_p += _jacobianOplusXi.transpose()* RobWInformation* _jacobianOplusXi;
		g_p -= _jacobianOplusXi.transpose()* RobWInformation* PointError2D;
		e_p += PointError2D.transpose()* RobWInformation* PointError2D;
		N_p++;
	    }	    
	}  
    }
    
    int N_l = 0;
    // Set MapPoint vertices
    const int NL = pFrame->NL;
    {
	unique_lock<mutex> lock(MapLine::mGlobalMutex);
	for(int i=0; i<NL; i++)
	{
	    MapLine* pML = pFrame->mvpMapLines[i];
	    if(pML && !pFrame->mvbOutlierLines[i])
	    {
		//获取观测值:线方向
                const KeyLine &kLUn = pFrame->mvLinesUn[i];//线特征
		Vector3d sp_l; sp_l << kLUn.startPointX, kLUn.startPointY, 1.0;
		Vector3d ep_l; ep_l << kLUn.endPointX,   kLUn.endPointY,   1.0;
		Vector3d le_l; le_l << sp_l.cross(ep_l);
		//le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) + le_l(2)*le_l(2));
		le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1));
		
		//设置标准信息矩阵(变量）
		const float invSigma2 = pFrame->mvInvLevelSigma2Lines[kLUn.octave];
		const float Information = invSigma2;
		float Informationsqr = sqrt(Information);
				
		//获取3D预测值和2D预测值
		cv::Mat MidXw = pML->GetMidPWorldPos();
		Vector3d MidXYZ; MidXYZ <<  MidXw.at<float>(0), MidXw.at<float>(1), MidXw.at<float>(2);//地图线中点世界坐标
		Vector3d Trans_MidXYZ = TcwBothUpdate.map(MidXYZ);//地图线中点相机坐标
		Vector2d MidPointProjection2D;//地图线中点的图像坐标系坐标		
		MidPointProjection2D(0) = fx*Trans_MidXYZ(0)/Trans_MidXYZ(2) + cx;
		MidPointProjection2D(1) = fy*Trans_MidXYZ(1)/Trans_MidXYZ(2) + cy;
		float MidPointReProjectionError = le_l(0)*MidPointProjection2D(0) + le_l(1)*MidPointProjection2D(1) + le_l(2);
		
		//标准信息矩阵下重投影误差的模平方
		float MidPointReProjectionErrorNorm = Informationsqr*MidPointReProjectionError;         
		float MidPointReProjectionErrorNorm2 = MidPointReProjectionErrorNorm*MidPointReProjectionErrorNorm;
		//基于Huber核函数计算权重
		float rho = RobustWeightHuber(MidPointReProjectionErrorNorm2, deltaMonoLine);
		float RobWInformation = rho*Information;
		
		//线特征雅克比矩阵参数
		float PX = MidXw.at<float>(0);
		float PY = MidXw.at<float>(1);
		float PZinv = 1 / MidXw.at<float>(2);
		float PZ2inv = 1 / ( MidXw.at<float>(2) * MidXw.at<float>(2) );
		
		//线特征雅克比矩阵(PL-SLAM(stereo)的这里有问题)
		//先旋转后平移
		Matrix<double,1,6> _jacobianOplusXi;
		_jacobianOplusXi(0,0) =  -PX*PY*PZ2inv *le_l(0)*fx - (1+PY*PY*PZ2inv) *le_l(1)*fy;
		_jacobianOplusXi(0,1) = (1+(PX*PX*PZ2inv)) *le_l(0)*fx + PX*PY*PZ2inv *le_l(1)*fy;
		_jacobianOplusXi(0,2) = -PY*PZinv *le_l(0)*fx + PX*PZinv *le_l(1)*fy;
		_jacobianOplusXi(0,3) = PZinv *le_l(0)*fx;
		_jacobianOplusXi(0,4) = PZinv *le_l(1)*fy;
		_jacobianOplusXi(0,5) = -PX*PZ2inv *le_l(0)*fx - PY*PZ2inv *le_l(1)*fy;
		
		//点特征的H、g、参与优化的点误差累加和、参与优化的点特征数量
		H_l += _jacobianOplusXi.transpose()* RobWInformation* _jacobianOplusXi;
		g_l -= _jacobianOplusXi.transpose()* RobWInformation* MidPointReProjectionError;
		e_l += MidPointReProjectionError* RobWInformation* MidPointReProjectionError;
		N_l++;
	    }   
	}
    }
    // 整合：H,g,总体error
    H = H_p + H_l;
    g = g_p + g_l;
    error = e_p + e_l;

    // 单位误差计算
    error /= (N_l+N_p);
    
}

//通过输入的重投影误差的模平方以及地图点的双自由度单像素误差
//或者地图线的单自由度单像素误差，通过返回值返回权重值rho（float）
float Optimizer::RobustWeightHuber( float PointReProjectionErrorNorm2, float deta)
{
   float DetaSqr = deta*deta;
   
   if(PointReProjectionErrorNorm2<=DetaSqr)
   {
      return 1.0;
   }
   else
   {
      double sqrte = sqrt(PointReProjectionErrorNorm2); // absolut value of the error    
      return deta / sqrte;
   }
}

//基于当前帧位姿设置外点和外线，内点数和内线数：
void Optimizer::RecheckBothOutliersForPose(Frame* pFrame, int &inlinersNum, int &inlinerLinesNum)
{
  
    if(pFrame->mTcw.empty())
      return;
    
    //基于当前图像帧位姿，分别进行所有（无论之前是否设置为外点）地图点和地图线中点重投影误差计算，并设置外点
    g2o::SE3Quat TcwBothUpdate = Converter::toSE3Quat(pFrame->mTcw);
       
    //标定矩阵
    cv::Mat K = pFrame->mK;
    const float fx = K.at<float>(0,0);
    const float fy = K.at<float>(1,1);
    const float cx = K.at<float>(0,2);
    const float cy = K.at<float>(1,2);
    
    inlinersNum = 0;
    inlinerLinesNum = 0;
    
    const int N = pFrame->N;
    for(int i = 0; i<N; i++)
    {
	MapPoint* pMP = pFrame->mvpMapPoints[i];
	if(pMP)
	{	    
	    //获取观测值
	    const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
	    Vector2d Point2D; Point2D << kpUn.pt.x, kpUn.pt.y;
		
	    //设置标准信息矩阵
	    const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
	    Matrix2d Information = Matrix2d::Identity()*invSigma2;
				
	    //获取3D预测值和2D预测值
	    cv::Mat Xw = pMP->GetWorldPos();
	    Vector3d XYZ; XYZ <<  Xw.at<float>(0), Xw.at<float>(1), Xw.at<float>(2);//地图点世界坐标
	    Vector3d Trans_XYZ = TcwBothUpdate.map(XYZ);//地图点相机坐标
	    Vector2d PointProjection2D;//地图点的图像坐标系坐标		
	    PointProjection2D(0) = fx*Trans_XYZ(0)/Trans_XYZ(2) + cx;
	    PointProjection2D(1) = fy*Trans_XYZ(1)/Trans_XYZ(2) + cy;
	    Vector2d PointError2D; PointError2D << Point2D(0)-PointProjection2D(0), Point2D(1)-PointProjection2D(1);
		
	    //标准信息矩阵下重投影误差的模平方
	    float PointReProjectionErrorNorm2 = PointError2D.dot(Information*PointError2D);
		
	    //判断是否外点重置
	    if(PointReProjectionErrorNorm2<=5.991)
	    {
		pFrame->mvbOutlier[i]=false;
		inlinersNum++;
	    }
	    else
		pFrame->mvbOutlier[i]=true;	    
	}
    }
    
    const int NL = pFrame->NL;
    for(int i = 0; i<NL; i++)
    {
	MapLine* pML = pFrame->mvpMapLines[i];
	if(pML)
	{
	    //获取观测值:线方向
	    const KeyLine &kLUn = pFrame->mvLinesUn[i];//线特征
	    Vector3d sp_l; sp_l << kLUn.startPointX, kLUn.startPointY, 1.0;
	    Vector3d ep_l; ep_l << kLUn.endPointX,   kLUn.endPointY,   1.0;
	    Vector3d le_l; le_l << sp_l.cross(ep_l);
	    //le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) + le_l(2)*le_l(2));
	    le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1));
		
	    //设置标准信息矩阵(变量）
	    const float invSigma2 = pFrame->mvInvLevelSigma2Lines[kLUn.octave];
	    float Informationsqr = sqrt(invSigma2);
				
	    //获取3D预测值和2D预测值
	    cv::Mat MidXw = pML->GetMidPWorldPos();
	    Vector3d MidXYZ; MidXYZ <<  MidXw.at<float>(0), MidXw.at<float>(1), MidXw.at<float>(2);//地图线中点世界坐标
	    Vector3d Trans_MidXYZ = TcwBothUpdate.map(MidXYZ);//地图线中点相机坐标
	    Vector2d MidPointProjection2D;//地图线中点的图像坐标系坐标		
	    MidPointProjection2D(0) = fx*Trans_MidXYZ(0)/Trans_MidXYZ(2) + cx;
	    MidPointProjection2D(1) = fy*Trans_MidXYZ(1)/Trans_MidXYZ(2) + cy;
	    float MidPointReProjectionError = le_l(0)*MidPointProjection2D(0) + le_l(1)*MidPointProjection2D(1) + le_l(2);
		
	    //标准信息矩阵下重投影误差的模平方
	    float MidPointReProjectionErrorNorm = Informationsqr*MidPointReProjectionError;         
	    float MidPointReProjectionErrorNorm2 = MidPointReProjectionErrorNorm*MidPointReProjectionErrorNorm;
		
	    //内外线重置
	    if(MidPointReProjectionErrorNorm2<=3.841)
	    {
		pFrame->mvbOutlierLines[i]=false;
		inlinerLinesNum++;
	    }
	    else
		pFrame->mvbOutlierLines[i]=true;
	}
    }
}

//基于当前帧位姿设置外点，内点数：
//需要直接调用，因此为公有
void Optimizer::SetOutliersForPose(Frame* pFrame, int &inlinersNum)
{
    if(pFrame->mTcw.empty())
      return;
    
    //基于当前图像帧位姿，分别进行所有（无论之前是否设置为外点）地图点和地图线中点重投影误差计算，并设置外点
    g2o::SE3Quat TcwBothUpdate = Converter::toSE3Quat(pFrame->mTcw);
       
    //标定矩阵
    cv::Mat K = pFrame->mK;
    const float fx = K.at<float>(0,0);
    const float fy = K.at<float>(1,1);
    const float cx = K.at<float>(0,2);
    const float cy = K.at<float>(1,2);
    
    inlinersNum = 0;
    
    const int N = pFrame->N;
    for(int i = 0; i<N; i++)
    {
	MapPoint* pMP = pFrame->mvpMapPoints[i];
	if(pMP)
	{	    
	    //获取观测值
	    const cv::KeyPoint &kpUn = pFrame->mvKeysUn[i];
	    Vector2d Point2D; Point2D << kpUn.pt.x, kpUn.pt.y;
		
	    //设置标准信息矩阵
	    const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
	    Matrix2d Information = Matrix2d::Identity()*invSigma2;
				
	    //获取3D预测值和2D预测值
	    cv::Mat Xw = pMP->GetWorldPos();
	    Vector3d XYZ; XYZ <<  Xw.at<float>(0), Xw.at<float>(1), Xw.at<float>(2);//地图点世界坐标
	    Vector3d Trans_XYZ = TcwBothUpdate.map(XYZ);//地图点相机坐标
	    Vector2d PointProjection2D;//地图点的图像坐标系坐标		
	    PointProjection2D(0) = fx*Trans_XYZ(0)/Trans_XYZ(2) + cx;
	    PointProjection2D(1) = fy*Trans_XYZ(1)/Trans_XYZ(2) + cy;
	    Vector2d PointError2D; PointError2D << Point2D(0)-PointProjection2D(0), Point2D(1)-PointProjection2D(1);
		
	    //标准信息矩阵下重投影误差的模平方
	    float PointReProjectionErrorNorm2 = PointError2D.dot(Information*PointError2D);
		
	    //判断是否外点重置
	    if(PointReProjectionErrorNorm2<=5.991)
	    {
		pFrame->mvbOutlier[i]=false;
		inlinersNum++;
	    }
	    else
		pFrame->mvbOutlier[i]=true;	    
	}
    }
}

//基于当前帧位姿设置外线，内线数：
//需要直接调用，因此为公有
void Optimizer::SetOutlierLinesForPose(Frame* pFrame, int &inlinerLinesNum)
{
    if(pFrame->mTcw.empty())
      return;
    
    //基于当前图像帧位姿，分别进行所有（无论之前是否设置为外点）地图点和地图线中点重投影误差计算，并设置外点
    g2o::SE3Quat TcwBothUpdate = Converter::toSE3Quat(pFrame->mTcw);
       
    //标定矩阵
    cv::Mat K = pFrame->mK;
    const float fx = K.at<float>(0,0);
    const float fy = K.at<float>(1,1);
    const float cx = K.at<float>(0,2);
    const float cy = K.at<float>(1,2);
    
    inlinerLinesNum = 0;
    
    const int NL = pFrame->NL;
    for(int i = 0; i<NL; i++)
    {
	MapLine* pML = pFrame->mvpMapLines[i];
	if(pML)
	{

	    //获取观测值:线方向
	    const KeyLine &kLUn = pFrame->mvLinesUn[i];//线特征
	    Vector3d sp_l; sp_l << kLUn.startPointX, kLUn.startPointY, 1.0;
	    Vector3d ep_l; ep_l << kLUn.endPointX,   kLUn.endPointY,   1.0;
	    Vector3d le_l; le_l << sp_l.cross(ep_l);
	    //le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) + le_l(2)*le_l(2));
	    le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1));
		
	    //设置标准信息矩阵(变量）
	    const float invSigma2 = pFrame->mvInvLevelSigma2Lines[kLUn.octave];
	    float Informationsqr = sqrt(invSigma2);
				
	    //获取3D预测值和2D预测值
	    cv::Mat MidXw = pML->GetMidPWorldPos();
	    Vector3d MidXYZ; MidXYZ <<  MidXw.at<float>(0), MidXw.at<float>(1), MidXw.at<float>(2);//地图线中点世界坐标
	    Vector3d Trans_MidXYZ = TcwBothUpdate.map(MidXYZ);//地图线中点相机坐标
	    Vector2d MidPointProjection2D;//地图线中点的图像坐标系坐标		
	    MidPointProjection2D(0) = fx*Trans_MidXYZ(0)/Trans_MidXYZ(2) + cx;
	    MidPointProjection2D(1) = fy*Trans_MidXYZ(1)/Trans_MidXYZ(2) + cy;
	    float MidPointReProjectionError = le_l(0)*MidPointProjection2D(0) + le_l(1)*MidPointProjection2D(1) + le_l(2);
		
	    //标准信息矩阵下重投影误差的模平方
	    float MidPointReProjectionErrorNorm = Informationsqr*MidPointReProjectionError;         
	    float MidPointReProjectionErrorNorm2 = MidPointReProjectionErrorNorm*MidPointReProjectionErrorNorm;
		
	    //内外线重置
	    if(MidPointReProjectionErrorNorm2<=3.841)
	    {
		pFrame->mvbOutlierLines[i]=false;
		inlinerLinesNum++;
	    }
	    else
		pFrame->mvbOutlierLines[i]=true;
	}
    }
}

//执行两次点线同时优化
void Optimizer::LocalBundleAdjustmentmainOld(KeyFrame* pKF, bool *pbStopFlag, Map *pMap)
{
    
  // Local KeyFrames: First Breath Search from Current Keyframe
    list<KeyFrame*> lLocalKeyFrames;

    lLocalKeyFrames.push_back(pKF);
    pKF->mnBALocalForKFBoth = pKF->mnId;
    
    const vector<KeyFrame*> vNeighKFs = pKF->GetVectorCovisibleKeyFrames();
    for(int i=0, iend=vNeighKFs.size(); i<iend; i++)
    {
        KeyFrame* pKFi = vNeighKFs[i];
	pKFi->mnBALocalForKFBoth = pKF->mnId; 
        if(!pKFi->isBad() && !pKFi->isBadLines() )
            lLocalKeyFrames.push_back(pKFi);
    }

    list<MapPoint*> lLocalMapPoints;
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        vector<MapPoint*> vpMPs = (*lit)->GetMapPointMatches();
        for(vector<MapPoint*>::iterator vit=vpMPs.begin(), vend=vpMPs.end(); vit!=vend; vit++)
        {
            MapPoint* pMP = *vit;
            if(pMP)
                if(!pMP->isBad())
                    if(pMP->mnBALocalForKFBoth!=pKF->mnId)
                    {
                        lLocalMapPoints.push_back(pMP);
                        pMP->mnBALocalForKFBoth=pKF->mnId;
                    }
        }
    }
    
    list<KeyFrame*> lLocalKeyFramesLines;
    lLocalKeyFramesLines.push_back(pKF);
    
    const vector<KeyFrame*> vNeighKFsLines = pKF->GetVectorCovisibleKeyFramesLines();
    for(int i=0, iend=vNeighKFsLines.size(); i<iend; i++)
    {
        KeyFrame* pKFi = vNeighKFsLines[i];
	
	if(pKFi->isBadLines())
	    continue;
	
	if(pKFi->isBad())
	    continue;
	
	if(pKFi->mnBALocalForKFBoth!=pKF->mnId)
	{
	    pKFi->mnBALocalForKFBoth = pKF->mnId;
	    lLocalKeyFrames.push_back(pKFi);
	}	  
	
	lLocalKeyFramesLines.push_back(pKFi);
    }
    
    list<MapLine*> lLocalMapLines;
    for(list<KeyFrame*>::iterator lit=lLocalKeyFramesLines.begin() , lend=lLocalKeyFramesLines.end(); lit!=lend; lit++)
    {
        vector<MapLine*> vpMLs = (*lit)->GetMapLineMatches();
        for(vector<MapLine*>::iterator vit=vpMLs.begin(), vend=vpMLs.end(); vit!=vend; vit++)
        {
            MapLine* pML = *vit;
            if(pML)
                if(!pML->isBad())
                    if(pML->mnBALocalForKFBoth!=pKF->mnId)
                    {
                        lLocalMapLines.push_back(pML);
                        pML->mnBALocalForKFBoth=pKF->mnId;
                    }
        }
    }
    
    list<KeyFrame*> lFixedCameras;
    //a.基于局部地图点的固定关键帧
    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        map<KeyFrame*,size_t> observations = (*lit)->GetObservations();
        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(pKFi->mnBALocalForKFBoth!=pKF->mnId && pKFi->mnBAFixedForKFBoth!=pKF->mnId)
            {                
                pKFi->mnBAFixedForKFBoth=pKF->mnId;
                if(!pKFi->isBad() && !pKFi->isBadLines())
                    lFixedCameras.push_back(pKFi);
            }
        }
    }
    //b.基于局部地图线的固定关键帧
    for(list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        map<KeyFrame*,size_t> observations = (*lit)->GetObservations();
        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(pKFi->mnBALocalForKFBoth!=pKF->mnId && pKFi->mnBAFixedForKFBoth!=pKF->mnId)
            {                
                pKFi->mnBAFixedForKFBoth=pKF->mnId;
                if(!pKFi->isBadLines() && !pKFi->isBad())
                    lFixedCameras.push_back(pKFi);
            }
        }
    }
    

    // 地图线顶点的首尾端点必须各自独立成边，否则不满足BlockSolver_6_3形式
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType *linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);
    
    if(pbStopFlag)
        optimizer.setForceStopFlag(pbStopFlag);
    
    unsigned long maxKFid = 0;

    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(pKFi->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }
    
    for(list<KeyFrame*>::iterator lit=lFixedCameras.begin(), lend=lFixedCameras.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(true);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }
    
    // Set MapPoint vertices
    const int nExpectedSizePoints = (lLocalKeyFrames.size()+lFixedCameras.size())*lLocalMapPoints.size();

    vector<g2o::EdgeSE3ProjectXYZ*> vpEdgesMonoPoints;
    vpEdgesMonoPoints.reserve(nExpectedSizePoints);

    vector<KeyFrame*> vpEdgeKFMonoPoints;
    vpEdgeKFMonoPoints.reserve(nExpectedSizePoints);

    vector<MapPoint*> vpMapPointEdgeMono;
    vpMapPointEdgeMono.reserve(nExpectedSizePoints);
    
    const float thHuberMonoPoints = sqrt(5.991);
    
    unsigned long number = 0;
    //基于局部地图点设置边
    unsigned long maxMapPointAndKFid = 0;
    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
        vPoint->setEstimate(Converter::toVector3d(pMP->GetWorldPos()));
	
        unsigned long id = pMP->mnId+maxKFid+1;
        vPoint->setId(id);
	if(id>maxMapPointAndKFid)
            maxMapPointAndKFid=id;
	
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

        const map<KeyFrame*,size_t> observations = pMP->GetObservations();

        //Set edges
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(!pKFi->isBad() && !pKFi->isBadLines())
            {                
                const cv::KeyPoint &kpUn = pKFi->mvKeysUn[mit->second];

                // Monocular observation
                if(pKFi->mvuRight[mit->second]<0)
                {
                    Eigen::Matrix<double,2,1> obs;
                    obs << kpUn.pt.x, kpUn.pt.y;

                    g2o::EdgeSE3ProjectXYZ* e = new g2o::EdgeSE3ProjectXYZ();

                    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                    e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
                    e->setMeasurement(obs);
                    const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                    e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                    g2o::RobustKernelCauchy* rk = new g2o::RobustKernelCauchy;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuberMonoPoints);

                    e->fx = pKFi->fx;
                    e->fy = pKFi->fy;
                    e->cx = pKFi->cx;
                    e->cy = pKFi->cy;

                    optimizer.addEdge(e);
                    vpEdgesMonoPoints.push_back(e);
                    vpEdgeKFMonoPoints.push_back(pKFi);
                    vpMapPointEdgeMono.push_back(pMP);
		    number++;
		}
	    }
	}
    }
  
    if(pbStopFlag)
        if(*pbStopFlag)
            return;
	
     // Set MapPoint vertices
    const int nExpectedSizeLines = (lLocalKeyFrames.size()+lFixedCameras.size())*lLocalMapLines.size();

    vector<EdgeSE3ProjectXYZLines*> vpEdgesMonoLines;
    vpEdgesMonoLines.reserve(2*nExpectedSizeLines);

    vector<KeyFrame*> vpEdgeKFMonoLines;
    vpEdgeKFMonoLines.reserve(nExpectedSizeLines);

    vector<MapLine*> vpMapLineEdgeMono;
    vpMapLineEdgeMono.reserve(nExpectedSizeLines);	
    
    const float thHuberMonoLines = sqrt(3.841);
	
    for(list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        MapLine* pML = *lit;

        g2o::VertexSBAPointXYZ* vFirPoint = new g2o::VertexSBAPointXYZ();
        vFirPoint->setEstimate(Converter::toVector3d(pML->GetFirPWorldPos()));
        int idFir = 2*pML->mnId+maxMapPointAndKFid+1;
        vFirPoint->setId(idFir);
        vFirPoint->setMarginalized(true);
        optimizer.addVertex(vFirPoint);

        g2o::VertexSBAPointXYZ* vEndPoint = new g2o::VertexSBAPointXYZ();
        vEndPoint->setEstimate(Converter::toVector3d(pML->GetEndPWorldPos()));
        int idEnd = 2*pML->mnId+1+maxMapPointAndKFid+1;
        vEndPoint->setId(idEnd);
        vEndPoint->setMarginalized(true);
        optimizer.addVertex(vEndPoint);

        const map<KeyFrame*,size_t> observations = pML->GetObservations();

        //Set edges
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
	    
            KeyFrame* pKFi = mit->first;   

            if(!pKFi->isBadLines() && !pKFi->isBad() )
            {                
                const KeyLine &kLUn = pKFi->mvLinesUn[mit->second];
		Vector3d sp_l; sp_l << kLUn.startPointX, kLUn.startPointY, 1.0;
		Vector3d ep_l; ep_l << kLUn.endPointX,   kLUn.endPointY,   1.0;
		Vector3d le_l; le_l << sp_l.cross(ep_l);
		//le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) + le_l(2)*le_l(2));
		le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) );
         
		//输入地图线首尾端点的观测值
		Eigen::Matrix<double,3,1> obs;
		obs << le_l(0) , le_l(1), le_l(2);

		EdgeSE3ProjectXYZLines* eFir = new EdgeSE3ProjectXYZLines();

		eFir->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(idFir)));
		eFir->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
		eFir->setMeasurement(obs);
		const float &invSigma2 = pKFi->mvInvLevelSigma2Lines[kLUn.octave];
		eFir->setInformation(Eigen::Matrix3d::Identity()*invSigma2);
		
		g2o::RobustKernelCauchy* rkFir = new g2o::RobustKernelCauchy;
		eFir->setRobustKernel(rkFir);
		rkFir->setDelta(thHuberMonoLines);
		
		eFir->fx = pKFi->fx;
		eFir->fy = pKFi->fy;
		eFir->cx = pKFi->cx;
		eFir->cy = pKFi->cy;
		number++;
		optimizer.addEdge(eFir);
		vpEdgesMonoLines.push_back(eFir);//存入首端点对应边
		
		
		EdgeSE3ProjectXYZLines* eEnd = new EdgeSE3ProjectXYZLines();

		eEnd->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(idEnd)));
		eEnd->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
		eEnd->setMeasurement(obs);
		eEnd->setInformation(Eigen::Matrix3d::Identity()*invSigma2);

		g2o::RobustKernelCauchy* rkEnd = new g2o::RobustKernelCauchy;
		eEnd->setRobustKernel(rkEnd);
		rkEnd->setDelta(thHuberMonoLines);

		eEnd->fx = pKFi->fx;
		eEnd->fy = pKFi->fy;
		eEnd->cx = pKFi->cx;
		eEnd->cy = pKFi->cy;

		number++;		
		optimizer.addEdge(eEnd);
		vpEdgesMonoLines.push_back(eEnd);
			
		vpEdgeKFMonoLines.push_back(pKFi);		    		 	    		    
		vpMapLineEdgeMono.push_back(pML);
		    		    	    		    
	    }
	}
    }   
	
    optimizer.initializeOptimization();
    optimizer.optimize(5);
    
    
    bool bDoMore= true;

    if(pbStopFlag)
        if(*pbStopFlag)
            bDoMore = false;
	
    if(bDoMore)
    {

    // Check inlier observations
    for(size_t i=0, iend=vpEdgesMonoPoints.size(); i<iend;i++)
    {
        g2o::EdgeSE3ProjectXYZ* e = vpEdgesMonoPoints[i];
        MapPoint* pMP = vpMapPointEdgeMono[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>5.991 || !e->isDepthPositive())
        {
            e->setLevel(1);
        }

        e->setRobustKernel(0);
    }
    
    //检测内线
    for(size_t i=0, iend=vpEdgesMonoLines.size(); i<iend;i=i+2)
    {
        EdgeSE3ProjectXYZLines* eFir = vpEdgesMonoLines[i];
	EdgeSE3ProjectXYZLines* eEnd = vpEdgesMonoLines[i+1];
	
        MapLine* pML = vpMapLineEdgeMono[i / 2];
	
	if(pML->isBad())
            continue;
	 if(eFir->chi2()+eEnd->chi2()>5.991 || !eFir->isDepthPositive() || !eEnd->isDepthPositive())
        {
	    eFir->setLevel(1);
	    eEnd->setLevel(1);
	}
	
	eFir->setRobustKernel(0);
	eEnd->setRobustKernel(0);
    }
    
    optimizer.initializeOptimization(0);
    optimizer.optimize(5);
    
    }
    
    vector<pair<KeyFrame*,MapPoint*> > vToErasePoints;
    vToErasePoints.reserve(vpEdgesMonoPoints.size());
    for(size_t i=0, iend=vpEdgesMonoPoints.size(); i<iend;i++)
    {
        g2o::EdgeSE3ProjectXYZ* e = vpEdgesMonoPoints[i];
        MapPoint* pMP = vpMapPointEdgeMono[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>5.991 || !e->isDepthPositive())
        {
            KeyFrame* pKFi = vpEdgeKFMonoPoints[i];
            vToErasePoints.push_back(make_pair(pKFi,pMP));
        }
    }
    
    vector<pair<KeyFrame*,MapLine*> > vToEraseLines;
    vToEraseLines.reserve(vpEdgesMonoLines.size());
    
     for(size_t i=0, iend=vpEdgesMonoLines.size(); i<iend;i=i+2)
    {
        EdgeSE3ProjectXYZLines* eFir = vpEdgesMonoLines[i];
	EdgeSE3ProjectXYZLines* eEnd = vpEdgesMonoLines[i+1];
	
        MapLine* pML = vpMapLineEdgeMono[i / 2];
	
	if(pML->isBad())
            continue;
	 if(eFir->chi2()+eEnd->chi2()>5.991 || !eFir->isDepthPositive() || !eEnd->isDepthPositive())
        {
	    KeyFrame* pKFi = vpEdgeKFMonoLines[i / 2];
	    vToEraseLines.push_back(make_pair(pKFi,pML));
	}
    }
    
    // Get Map Mutex
    unique_lock<mutex> lock(pMap->mMutexMapUpdate);
    
    //删除坏点的观测
    if(!vToErasePoints.empty())
    {
        for(size_t i=0;i<vToErasePoints.size();i++)
        {
            KeyFrame* pKFi = vToErasePoints[i].first;
            MapPoint* pMPi = vToErasePoints[i].second;
            pKFi->EraseMapPointMatch(pMPi);
            pMPi->EraseObservation(pKFi);
        }
    }
    
    //删除坏的线观测
    if(!vToEraseLines.empty())
    {
        for(size_t i=0;i<vToEraseLines.size();i++)
        {
            KeyFrame* pKFi = vToEraseLines[i].first;
            MapLine* pMLi = vToEraseLines[i].second;
            pKFi->EraseMapLineMatch(pMLi);
            pMLi->EraseObservation(pKFi);
        }
    }
    
    // Recover optimized data
    
    //Keyframes
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKF = *lit;
        g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
        pKF->SetPose(Converter::toCvMat(SE3quat));
    }
    
    //Points
    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexSBAPointXYZ* vPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(pMP->mnId+maxKFid+1));
        pMP->SetWorldPos(Converter::toCvMat(vPoint->estimate()));
        pMP->UpdateNormalAndDepth();
    }
    
    //Lines
    for(list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        MapLine* pML = *lit;

        g2o::VertexSBAPointXYZ* vFirPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(2*pML->mnId+maxMapPointAndKFid+1));
	cv::Mat FirPWorldPos = Converter::toCvMat(vFirPoint->estimate());
        pML->SetFirPWorldPos(FirPWorldPos);

	g2o::VertexSBAPointXYZ* vEndPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(2*pML->mnId+1+maxMapPointAndKFid+1));
	cv::Mat EndPWorldPos = Converter::toCvMat(vEndPoint->estimate());
        pML->SetEndPWorldPos(EndPWorldPos);

	cv::Mat MidPWorldPos(3,1,CV_32F);
	for(int i=0;i<3;i++)
            MidPWorldPos.at<float>(i)=
            0.5*( FirPWorldPos.at<float>(i)+EndPWorldPos.at<float>(i) );	
	pML->SetMidPWorldPos(MidPWorldPos);

        pML->UpdateNormalAndDepth();
    }
    
}
//复合局部BA优化函数：（后端）
void Optimizer::LocalBundleAdjustmentmain(KeyFrame* pKF, bool *pbStopFlag, Map *pMap)
{
  
    //设置点特征及线特征的关键帧单位误差返回数据结构
    std::map<KeyFrame*,float> PoseUnitErrorForPointsLocalBA;
    std::map<KeyFrame*,float> PoseUnitErrorForLinesLocalBA;
    
    
    //调用调用g2o完成仅点和仅线的局部BA优化
    thread threadLocalBundleAdjustmentPoints(&Optimizer::LocalBundleAdjustmentPoints, pKF, pbStopFlag, 
				    pMap, ref(PoseUnitErrorForPointsLocalBA) );
    thread threadLocalBundleAdjustmentLines(&Optimizer::LocalBundleAdjustmentLines, pKF, pbStopFlag, 
				    pMap, ref(PoseUnitErrorForLinesLocalBA) );               
    
    threadLocalBundleAdjustmentPoints.join();
    threadLocalBundleAdjustmentLines.join();
    
    //更新当前关键帧的共视图连接关系（点和线共视图）
    LocalBAPoseDecidingBetweenLinesAndPoints(pKF, pMap, PoseUnitErrorForPointsLocalBA, PoseUnitErrorForLinesLocalBA);

    if(pbStopFlag)
        if(*pbStopFlag)
            return;
    
    //基于当前关键帧的点和线共视图进行局部BA优化
    LocalBundleAdjustmentBoth(pKF, pbStopFlag, pMap);
    
}

void Optimizer::LocalBundleAdjustmentmainNew(KeyFrame* pKF, bool *pbStopFlag, Map *pMap)
{
    //先执行仅点的局部BA优化
    std::map<KeyFrame*,float> PoseUnitErrorForPointsLocalBA;
    LocalBundleAdjustmentPoints(pKF,pbStopFlag,pMap,PoseUnitErrorForPointsLocalBA);
    
    // 点特征的局部共视图
    list<KeyFrame*> lLocalKeyFrames;
    lLocalKeyFrames.push_back(pKF);
    const vector<KeyFrame*> vNeighKFs = pKF->GetVectorCovisibleKeyFrames();
    for(int i=0, iend=vNeighKFs.size(); i<iend; i++)
    {
	KeyFrame* pKFi = vNeighKFs[i];
	if(!pKFi->isBad())
	    lLocalKeyFrames.push_back(pKFi);
    }
    
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
	KeyFrame* pKFi = *lit;
	pKFi->SetPose(pKFi->TcwPoints);
    }
    
    if(pbStopFlag)
        if(*pbStopFlag)
            return;
	
    pKF->UpdateConnections();
    
    //执行点线同时优化
    LocalBundleAdjustmentBoth(pKF, pbStopFlag, pMap);
    
}

//仅点的局部BA优化，与ORB-SLAM2的局部BA非常类似
void Optimizer::LocalBundleAdjustmentPoints(KeyFrame* pKF, bool *pbStopFlag, Map *pMap, std::map<KeyFrame*,float> &PoseUnitErrorForPointsLocalBA)
{
    // Local KeyFrames: First Breath Search from Current Keyframe
    list<KeyFrame*> lLocalKeyFrames;

    lLocalKeyFrames.push_back(pKF);
    pKF->mnBALocalForKF = pKF->mnId;
    
    const vector<KeyFrame*> vNeighKFs = pKF->GetVectorCovisibleKeyFrames();
    for(int i=0, iend=vNeighKFs.size(); i<iend; i++)
    {
        KeyFrame* pKFi = vNeighKFs[i];
        pKFi->mnBALocalForKF = pKF->mnId;
        if(!pKFi->isBad())
            lLocalKeyFrames.push_back(pKFi);
    }

    // Local MapPoints seen in Local KeyFrames
    list<MapPoint*> lLocalMapPoints;
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        vector<MapPoint*> vpMPs = (*lit)->GetMapPointMatches();
        for(vector<MapPoint*>::iterator vit=vpMPs.begin(), vend=vpMPs.end(); vit!=vend; vit++)
        {
            MapPoint* pMP = *vit;
            if(pMP)
                if(!pMP->isBad())
                    if(pMP->mnBALocalForKF!=pKF->mnId)
                    {
                        lLocalMapPoints.push_back(pMP);
                        pMP->mnBALocalForKF=pKF->mnId;
                    }
        }
    }

    // Fixed Keyframes. Keyframes that see Local MapPoints but that are not Local Keyframes
    list<KeyFrame*> lFixedCameras;
    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        map<KeyFrame*,size_t> observations = (*lit)->GetObservations();
        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(pKFi->mnBALocalForKF!=pKF->mnId && pKFi->mnBAFixedForKF!=pKF->mnId)
            {                
                pKFi->mnBAFixedForKF=pKF->mnId;
                if(!pKFi->isBad())
                    lFixedCameras.push_back(pKFi);
            }
        }
    }

    // Setup optimizer
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    if(pbStopFlag)
        optimizer.setForceStopFlag(pbStopFlag);

    unsigned long maxKFid = 0;

    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(pKFi->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }

    for(list<KeyFrame*>::iterator lit=lFixedCameras.begin(), lend=lFixedCameras.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(true);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }

    // Set MapPoint vertices
    const int nExpectedSize = (lLocalKeyFrames.size()+lFixedCameras.size())*lLocalMapPoints.size();

    vector<g2o::EdgeSE3ProjectXYZ*> vpEdgesMono;
    vpEdgesMono.reserve(nExpectedSize);

    vector<KeyFrame*> vpEdgeKFMono;
    vpEdgeKFMono.reserve(nExpectedSize);

    vector<MapPoint*> vpMapPointEdgeMono;
    vpMapPointEdgeMono.reserve(nExpectedSize);

    map<KeyFrame*, int> PoseUnitErrorNumForPointsLocalBA;
    const float thHuberMono = sqrt(5.991);
     
     for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
        vPoint->setEstimate(Converter::toVector3d(pMP->GetWorldPos()));
        int id = pMP->mnId+maxKFid+1;
        vPoint->setId(id);
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

        const map<KeyFrame*,size_t> observations = pMP->GetObservations();

        //Set edges
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(!pKFi->isBad())
            {                
                const cv::KeyPoint &kpUn = pKFi->mvKeysUn[mit->second];

                // Monocular observation
                if(pKFi->mvuRight[mit->second]<0)
                {
                    Eigen::Matrix<double,2,1> obs;
                    obs << kpUn.pt.x, kpUn.pt.y;

                    g2o::EdgeSE3ProjectXYZ* e = new g2o::EdgeSE3ProjectXYZ();

                    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                    e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
                    e->setMeasurement(obs);
                    const float invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                    e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                    g2o::RobustKernelCauchy* rk = new g2o::RobustKernelCauchy;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuberMono);

                    e->fx = pKFi->fx;
                    e->fy = pKFi->fy;
                    e->cx = pKFi->cx;
                    e->cy = pKFi->cy;

                    optimizer.addEdge(e);
                    vpEdgesMono.push_back(e);
                    vpEdgeKFMono.push_back(pKFi);
                    vpMapPointEdgeMono.push_back(pMP);
		}
	    }
	}
    }

    if(pbStopFlag)
        if(*pbStopFlag)
            return;

    optimizer.initializeOptimization();
    optimizer.optimize(5);
    
    vector<pair<KeyFrame*,MapPoint*> > vToErase;
    vToErase.reserve(vpEdgesMono.size());

    // Check inlier observations       
    for(size_t i=0, iend=vpEdgesMono.size(); i<iend;i++)
    {
        g2o::EdgeSE3ProjectXYZ* e = vpEdgesMono[i];
        MapPoint* pMP = vpMapPointEdgeMono[i];
	KeyFrame* pKFi = vpEdgeKFMono[i];
	
        if(pMP->isBad())
            continue;

        if(e->chi2()>5.991 || !e->isDepthPositive())
        {

            vToErase.push_back(make_pair(pKFi,pMP));
        }
        else
	{
	    if(!PoseUnitErrorForPointsLocalBA.count(pKFi))
	    {
		PoseUnitErrorForPointsLocalBA[pKFi] = e->chi2();
		PoseUnitErrorNumForPointsLocalBA[pKFi] = 1;
	    }
	    else
	    {
		PoseUnitErrorForPointsLocalBA[pKFi] += e->chi2();
		PoseUnitErrorNumForPointsLocalBA[pKFi] += 1;
	    }
	}
    }
    
    for(map<KeyFrame*, float>::iterator mit=PoseUnitErrorForPointsLocalBA.begin(), 
      mend=PoseUnitErrorForPointsLocalBA.end();mit!=mend;mit++)
    {
	mit->second = mit->second / PoseUnitErrorNumForPointsLocalBA[mit->first];                                                       
    }
      
    // Get Map Mutex
    unique_lock<mutex> lock(pMap->mMutexMapUpdate);

    if(!vToErase.empty())
    {
        for(size_t i=0;i<vToErase.size();i++)
        {
            KeyFrame* pKFi = vToErase[i].first;
            MapPoint* pMPi = vToErase[i].second;
            pKFi->EraseMapPointMatch(pMPi);
            pMPi->EraseObservation(pKFi);
        }
    }

    // Recover optimized data

    //Keyframes
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKF = *lit;
        g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
        pKF->SetPosePoints(Converter::toCvMat(SE3quat));
    }

    //Points
    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexSBAPointXYZ* vPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(pMP->mnId+maxKFid+1));
        pMP->SetWorldPos(Converter::toCvMat(vPoint->estimate()));
        pMP->UpdateNormalAndDepth();
    }
    
}

//仅线的局部BA优化，与ORB-SLAM2的局部BA非常类似
void Optimizer::LocalBundleAdjustmentLines(KeyFrame* pKF, bool *pbStopFlag, Map *pMap, std::map<KeyFrame*,float> &PoseUnitErrorForLinesLocalBA)
{
    // Local KeyFrames: First Breath Search from Current Keyframe
    list<KeyFrame*> lLocalKeyFrames;

    lLocalKeyFrames.push_back(pKF);
    pKF->mnBALocalForKFLines = pKF->mnId;

    const vector<KeyFrame*> vNeighKFs = pKF->GetVectorCovisibleKeyFramesLines();
    for(int i=0, iend=vNeighKFs.size(); i<iend; i++)
    {
        KeyFrame* pKFi = vNeighKFs[i];
        pKFi->mnBALocalForKFLines = pKF->mnId;
        if(!pKFi->isBadLines())
            lLocalKeyFrames.push_back(pKFi);
    }

    // Local MapLines seen in Local KeyFrames
    list<MapLine*> lLocalMapLines;
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        vector<MapLine*> vpMLs = (*lit)->GetMapLineMatches();
        for(vector<MapLine*>::iterator vit=vpMLs.begin(), vend=vpMLs.end(); vit!=vend; vit++)
        {
            MapLine* pML = *vit;
            if(pML)
                if(!pML->isBad())
                    if(pML->mnBALocalForKF!=pKF->mnId)
                    {
                        lLocalMapLines.push_back(pML);
                        pML->mnBALocalForKF=pKF->mnId;
                    }
        }
    }

    // Fixed Keyframes. Keyframes that see Local MapPoints but that are not Local Keyframes
    list<KeyFrame*> lFixedCameras;
    for(list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        map<KeyFrame*,size_t> observations = (*lit)->GetObservations();
        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(pKFi->mnBALocalForKFLines!=pKF->mnId && pKFi->mnBAFixedForKFLines!=pKF->mnId)
            {                
                pKFi->mnBAFixedForKFLines=pKF->mnId;
                if(!pKFi->isBadLines())
                    lFixedCameras.push_back(pKFi);
            }
        }
    }

    // Setup optimizer
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    if(pbStopFlag)
        optimizer.setForceStopFlag(pbStopFlag);

    unsigned long maxKFid = 0;

    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(pKFi->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }

    for(list<KeyFrame*>::iterator lit=lFixedCameras.begin(), lend=lFixedCameras.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(true);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }

    // Set MapPoint vertices
    const int nExpectedSize = (lLocalKeyFrames.size()+lFixedCameras.size())*lLocalMapLines.size();

    vector<EdgeSE3ProjectXYZLines*> vpEdgesMono;
    vpEdgesMono.reserve(2*nExpectedSize);

    vector<KeyFrame*> vpEdgeKFMono;
    vpEdgeKFMono.reserve(nExpectedSize);

    vector<MapLine*> vpMapLineEdgeMono;
    vpMapLineEdgeMono.reserve(nExpectedSize);
    
    map<KeyFrame*, int> PoseUnitErrorNumForLinesLocalBA;
    const float thHuberMono = sqrt(3.841);
     
     for(list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        MapLine* pML = *lit;

        g2o::VertexSBAPointXYZ* vFirPoint = new g2o::VertexSBAPointXYZ();
        vFirPoint->setEstimate(Converter::toVector3d(pML->GetFirPWorldPos()));
        int idFir = 2*pML->mnId+maxKFid+1;
        vFirPoint->setId(idFir);
        vFirPoint->setMarginalized(true);
        optimizer.addVertex(vFirPoint);

        g2o::VertexSBAPointXYZ* vEndPoint = new g2o::VertexSBAPointXYZ();
        vEndPoint->setEstimate(Converter::toVector3d(pML->GetEndPWorldPos()));
        int idEnd = 2*pML->mnId+1+maxKFid+1;
        vEndPoint->setId(idEnd);
        vEndPoint->setMarginalized(true);
        optimizer.addVertex(vEndPoint);

        const map<KeyFrame*,size_t> observations = pML->GetObservations();

        //Set edges
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(!pKFi->isBadLines())
            {                
                const KeyLine &kLUn = pKFi->mvLinesUn[mit->second];
		Vector3d sp_l; sp_l << kLUn.startPointX, kLUn.startPointY, 1.0;
		Vector3d ep_l; ep_l << kLUn.endPointX,   kLUn.endPointY,   1.0;
		Vector3d le_l; le_l << sp_l.cross(ep_l);
		//le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) + le_l(2)*le_l(2));
		le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1));
		
		Eigen::Matrix<double,3,1> obs;
		obs << le_l(0) , le_l(1), le_l(2);
		
		EdgeSE3ProjectXYZLines* eFir = new EdgeSE3ProjectXYZLines();

		eFir->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(idFir)));
		eFir->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
		eFir->setMeasurement(obs);
		const float invSigma2 = pKFi->mvInvLevelSigma2Lines[kLUn.octave];
		eFir->setInformation(Eigen::Matrix3d::Identity()*invSigma2);

		g2o::RobustKernelCauchy* rkFir = new g2o::RobustKernelCauchy;
		eFir->setRobustKernel(rkFir);
		rkFir->setDelta(thHuberMono);

		eFir->fx = pKFi->fx;
		eFir->fy = pKFi->fy;
		eFir->cx = pKFi->cx;
		eFir->cy = pKFi->cy;

		optimizer.addEdge(eFir);
		vpEdgesMono.push_back(eFir);
		EdgeSE3ProjectXYZLines* eEnd = new EdgeSE3ProjectXYZLines();

		eEnd->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(idEnd)));
		eEnd->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
		eEnd->setMeasurement(obs);
		eEnd->setInformation(Eigen::Matrix3d::Identity()*invSigma2);

		g2o::RobustKernelCauchy* rkEnd = new g2o::RobustKernelCauchy;
		eEnd->setRobustKernel(rkEnd);
		rkEnd->setDelta(thHuberMono);

		eEnd->fx = pKFi->fx;
		eEnd->fy = pKFi->fy;
		eEnd->cx = pKFi->cx;
		eEnd->cy = pKFi->cy;

		optimizer.addEdge(eEnd);
		vpEdgesMono.push_back(eEnd);
		vpEdgeKFMono.push_back(pKFi);		    		 	    		    
		vpMapLineEdgeMono.push_back(pML);
		    		    	    		    
	    }
	}
    }
    
    if(pbStopFlag)
        if(*pbStopFlag)
            return;

    optimizer.initializeOptimization();
    optimizer.optimize(5);
    
    vector<pair<KeyFrame*,MapLine*> > vToErase;
    vToErase.reserve(vpEdgesMono.size());

    // Check inlier observations       
    for(size_t i=0, iend=vpEdgesMono.size(); i<iend;i=i+2)
    {
        EdgeSE3ProjectXYZLines* eFir = vpEdgesMono[i];
	EdgeSE3ProjectXYZLines* eEnd = vpEdgesMono[i+1];
	
        MapLine* pML = vpMapLineEdgeMono[i / 2];
	KeyFrame* pKFi = vpEdgeKFMono[i / 2];
	
        if(pML->isBad())
            continue;

        if(eFir->chi2()+eEnd->chi2()>5.991 || !eFir->isDepthPositive() || !eEnd->isDepthPositive())
        {
	  
            vToErase.push_back(make_pair(pKFi,pML));
        }
        else
	{
	    if(!PoseUnitErrorForLinesLocalBA.count(pKFi))
	    {
		double LineError = eFir->chi2()+eEnd->chi2();
		PoseUnitErrorForLinesLocalBA[pKFi] = LineError;
		PoseUnitErrorNumForLinesLocalBA[pKFi] = 1;
	    }
	    else
	    {
		double LineError = eFir->chi2()+eEnd->chi2();
		PoseUnitErrorForLinesLocalBA[pKFi] += LineError;
		PoseUnitErrorNumForLinesLocalBA[pKFi] += 1;
	    }
	}
    }
    
    for(map<KeyFrame*, float>::iterator mit=PoseUnitErrorForLinesLocalBA.begin(), 
      mend=PoseUnitErrorForLinesLocalBA.end();mit!=mend;mit++)
    {
	mit->second = mit->second / PoseUnitErrorNumForLinesLocalBA[mit->first];                                                       
    }
    
    // Get Map Mutex
    unique_lock<mutex> lock(pMap->mMutexMapUpdate);

    if(!vToErase.empty())
    {
        for(size_t i=0;i<vToErase.size();i++)
        {
            KeyFrame* pKFi = vToErase[i].first;
            MapLine* pMLi = vToErase[i].second;
            pKFi->EraseMapLineMatch(pMLi);
            pMLi->EraseObservation(pKFi);
        }
    }
  
    // Recover optimized data

    //Keyframes
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKF = *lit;
        g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
        pKF->SetPoseLines(Converter::toCvMat(SE3quat));
    }

    //Lines
    for(list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        MapLine* pML = *lit;
	
        g2o::VertexSBAPointXYZ* vFirPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(2*pML->mnId+maxKFid+1));
	cv::Mat FirPWorldPos = Converter::toCvMat(vFirPoint->estimate());
        pML->SetFirPWorldPos(FirPWorldPos);
	
	g2o::VertexSBAPointXYZ* vEndPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(2*pML->mnId+1+maxKFid+1));
	cv::Mat EndPWorldPos = Converter::toCvMat(vEndPoint->estimate());
        pML->SetEndPWorldPos(EndPWorldPos);

	cv::Mat MidPWorldPos(3,1,CV_32F);
	for(int i=0;i<3;i++)
            MidPWorldPos.at<float>(i)=
            0.5*( FirPWorldPos.at<float>(i)+EndPWorldPos.at<float>(i) );	
	pML->SetMidPWorldPos(MidPWorldPos);

        pML->UpdateNormalAndDepth();
    }   
}

//b.后端的局部BA在执行点线特征同时优化前，需要对第一步分级优化的局部位姿进行全局设置
void Optimizer::LocalBAPoseDecidingBetweenLinesAndPoints(KeyFrame* pKF, Map *pMap, std::map<KeyFrame*,float> &PoseUnitErrorForPointsLocalBA, 
							 std::map<KeyFrame*,float > &PoseUnitErrorForLinesLocalBA )
{
    unique_lock<mutex> lock(pMap->mMutexMapUpdate);
    
    //上一步中，点线均执行优化
    if(!pKF->TcwLines.empty() && !pKF->TcwPoints.empty())
    {
      
	//（1）获取当前关键帧的点共视图和线的局部共视图（实际上就是中发生位姿调整的所有局部关键帧）。
	// 点特征的局部共视图
	list<KeyFrame*> lLocalKeyFrames;
	lLocalKeyFrames.push_back(pKF);
	const vector<KeyFrame*> vNeighKFs = pKF->GetVectorCovisibleKeyFrames();
	for(int i=0, iend=vNeighKFs.size(); i<iend; i++)
	{
	    KeyFrame* pKFi = vNeighKFs[i];
	    if(!pKFi->isBad())
		lLocalKeyFrames.push_back(pKFi);
	}
	// 线特征的局部共视图
	list<KeyFrame*> lLocalKeyFramesLines;
	lLocalKeyFramesLines.push_back(pKF);
	const vector<KeyFrame*> vNeighKFsLines = pKF->GetVectorCovisibleKeyFramesLines();
	for(int i=0, iend=vNeighKFsLines.size(); i<iend; i++)
	{
	    KeyFrame* pKFi = vNeighKFsLines[i];
	    if(!pKFi->isBadLines())
		lLocalKeyFramesLines.push_back(pKFi);
	}
	
	//（2）遍历点共视图所有关键帧，判断当前关键帧是否在步骤（1）返回的线关键帧map<KeyFrame*, int UnitError>集合中，
	// 若在，则进行单位误差比较，决定更优位姿，设置调整关键帧的mBAPoseLinesandPointsForKF为当前关键帧ID，
	// 若不在，则直接设置该关键帧位姿为点位姿。
	for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
	{
	    KeyFrame* pKFi = *lit ;
	    
	    if(PoseUnitErrorForLinesLocalBA.count(pKFi))
	    {

		//相同关键帧的单位误差比较,因为这里存在两组当前关键帧的位姿信息，择优录取
		if(PoseUnitErrorForPointsLocalBA[pKFi]<=PoseUnitErrorForLinesLocalBA[pKFi])
		{
		    pKFi->SetPose(pKFi->TcwPoints);
		    pKFi->mBAPoseLinesandPointsForKF = pKFi->mnId;
		}
		else
		{
		    pKFi->SetPose(pKFi->TcwLines);
		    pKFi->mBAPoseLinesandPointsForKF = pKFi->mnId;
		}
		    
	    }
	    else
	    {
		pKFi->SetPose(pKFi->TcwPoints);//表示当前关键帧只是点关键帧,不在线共视局部关键帧集中
		pKFi->mBAPoseLinesandPointsForKF = pKFi->mnId;
	    }
		
	}
	
	//（3）遍历线共视图所有关键帧。对于不为当前关键帧ID的关键帧，将其关键帧位姿设置为线位姿
	for(list<KeyFrame*>::iterator lit=lLocalKeyFramesLines.begin() , lend=lLocalKeyFramesLines.end(); lit!=lend; lit++)
	{
	    KeyFrame* pKFi = *lit;
	    
	    if(pKFi->mBAPoseLinesandPointsForKF != pKFi->mnId)
	    {
		pKFi->SetPose(pKFi->TcwLines);
	    }
	    else
		continue;
	}
	
	//（4）更新当前关键帧点和线的共视图连接关系。
	pKF->UpdateConnections();
	pKF->UpdateConnectionsLines();
    }
    else if(pKF->TcwLines.empty())
    {
	// 点特征的局部共视图
	list<KeyFrame*> lLocalKeyFrames;
	lLocalKeyFrames.push_back(pKF);
	const vector<KeyFrame*> vNeighKFs = pKF->GetVectorCovisibleKeyFrames();
	for(int i=0, iend=vNeighKFs.size(); i<iend; i++)
	{
	    KeyFrame* pKFi = vNeighKFs[i];
	    if(!pKFi->isBad())
		lLocalKeyFrames.push_back(pKFi);
	}
	
	for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
	{
	    KeyFrame* pKFi = *lit;
	    pKFi->SetPose(pKFi->TcwPoints);
	}
	
    }
    else if(pKF->TcwPoints.empty())
    {
	// 线特征的局部共视图
	list<KeyFrame*> lLocalKeyFramesLines;
	lLocalKeyFramesLines.push_back(pKF);
	const vector<KeyFrame*> vNeighKFsLines = pKF->GetVectorCovisibleKeyFramesLines();
	for(int i=0, iend=vNeighKFsLines.size(); i<iend; i++)
	{
	    KeyFrame* pKFi = vNeighKFsLines[i];
	    if(!pKFi->isBadLines())
		lLocalKeyFramesLines.push_back(pKFi);
	}
	
	for(list<KeyFrame*>::iterator lit=lLocalKeyFramesLines.begin() , lend=lLocalKeyFramesLines.end(); lit!=lend; lit++)
	{
	    KeyFrame* pKFi = *lit;
	    pKFi->SetPose(pKFi->TcwLines);
	}
	
    }
    else
	return;
    
}

//后端局部BA分级优化的第二步
void Optimizer::LocalBundleAdjustmentBoth(KeyFrame* pKF, bool *pbStopFlag, Map *pMap)
{
    // Local KeyFrames: First Breath Search from Current Keyframe
    list<KeyFrame*> lLocalKeyFrames;

    lLocalKeyFrames.push_back(pKF);
    pKF->mnBALocalForKFBoth = pKF->mnId;
    
    const vector<KeyFrame*> vNeighKFs = pKF->GetVectorCovisibleKeyFrames();
    for(int i=0, iend=vNeighKFs.size(); i<iend; i++)
    {
        KeyFrame* pKFi = vNeighKFs[i];
	pKFi->mnBALocalForKFBoth = pKF->mnId; 
        if(!pKFi->isBad() && !pKFi->isBadLines() )
            lLocalKeyFrames.push_back(pKFi);
    }

    list<MapPoint*> lLocalMapPoints;
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        vector<MapPoint*> vpMPs = (*lit)->GetMapPointMatches();
        for(vector<MapPoint*>::iterator vit=vpMPs.begin(), vend=vpMPs.end(); vit!=vend; vit++)
        {
            MapPoint* pMP = *vit;
            if(pMP)
                if(!pMP->isBad())
                    if(pMP->mnBALocalForKFBoth!=pKF->mnId)
                    {
                        lLocalMapPoints.push_back(pMP);
                        pMP->mnBALocalForKFBoth=pKF->mnId;
                    }
        }
    }
    
    list<KeyFrame*> lLocalKeyFramesLines;
    lLocalKeyFramesLines.push_back(pKF);
    
    const vector<KeyFrame*> vNeighKFsLines = pKF->GetVectorCovisibleKeyFramesLines();
    for(int i=0, iend=vNeighKFsLines.size(); i<iend; i++)
    {
        KeyFrame* pKFi = vNeighKFsLines[i];
	
	if(pKFi->isBadLines())
	    continue;
	
	if(pKFi->isBad())
	    continue;
	
	if(pKFi->mnBALocalForKFBoth!=pKF->mnId)
	{
	    pKFi->mnBALocalForKFBoth = pKF->mnId;
	    lLocalKeyFrames.push_back(pKFi);
	}	  
	
	lLocalKeyFramesLines.push_back(pKFi);
    }
    
    list<MapLine*> lLocalMapLines;
    for(list<KeyFrame*>::iterator lit=lLocalKeyFramesLines.begin() , lend=lLocalKeyFramesLines.end(); lit!=lend; lit++)
    {
        vector<MapLine*> vpMLs = (*lit)->GetMapLineMatches();
        for(vector<MapLine*>::iterator vit=vpMLs.begin(), vend=vpMLs.end(); vit!=vend; vit++)
        {
            MapLine* pML = *vit;
            if(pML)
                if(!pML->isBad())
                    if(pML->mnBALocalForKFBoth!=pKF->mnId)
                    {
                        lLocalMapLines.push_back(pML);
                        pML->mnBALocalForKFBoth=pKF->mnId;
                    }
        }
    }
    
    list<KeyFrame*> lFixedCameras;
    //a.基于局部地图点的固定关键帧
    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        map<KeyFrame*,size_t> observations = (*lit)->GetObservations();
        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(pKFi->mnBALocalForKFBoth!=pKF->mnId && pKFi->mnBAFixedForKFBoth!=pKF->mnId)
            {                
                pKFi->mnBAFixedForKFBoth=pKF->mnId;
                if(!pKFi->isBad() && !pKFi->isBadLines())
                    lFixedCameras.push_back(pKFi);
            }
        }
    }
    //b.基于局部地图线的固定关键帧
    for(list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        map<KeyFrame*,size_t> observations = (*lit)->GetObservations();
        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(pKFi->mnBALocalForKFBoth!=pKF->mnId && pKFi->mnBAFixedForKFBoth!=pKF->mnId)
            {                
                pKFi->mnBAFixedForKFBoth=pKF->mnId;
                if(!pKFi->isBadLines() && !pKFi->isBad())
                    lFixedCameras.push_back(pKFi);
            }
        }
    }
    

    // 地图线顶点的首尾端点必须各自独立成边，否则不满足BlockSolver_6_3形式
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType *linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);
    
    if(pbStopFlag)
        optimizer.setForceStopFlag(pbStopFlag);
    
    unsigned long maxKFid = 0;

    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(pKFi->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }
    
    for(list<KeyFrame*>::iterator lit=lFixedCameras.begin(), lend=lFixedCameras.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(true);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }
    
    // Set MapPoint vertices
    const int nExpectedSizePoints = (lLocalKeyFrames.size()+lFixedCameras.size())*lLocalMapPoints.size();

    vector<g2o::EdgeSE3ProjectXYZ*> vpEdgesMonoPoints;
    vpEdgesMonoPoints.reserve(nExpectedSizePoints);

    vector<KeyFrame*> vpEdgeKFMonoPoints;
    vpEdgeKFMonoPoints.reserve(nExpectedSizePoints);

    vector<MapPoint*> vpMapPointEdgeMono;
    vpMapPointEdgeMono.reserve(nExpectedSizePoints);
    
    const float thHuberMonoPoints = sqrt(5.991);
    
    unsigned long number = 0;
    //基于局部地图点设置边
    unsigned long maxMapPointAndKFid = 0;
    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
        vPoint->setEstimate(Converter::toVector3d(pMP->GetWorldPos()));
	
        unsigned long id = pMP->mnId+maxKFid+1;
        vPoint->setId(id);
	if(id>maxMapPointAndKFid)
            maxMapPointAndKFid=id;
	
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

        const map<KeyFrame*,size_t> observations = pMP->GetObservations();

        //Set edges
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(!pKFi->isBad() && !pKFi->isBadLines())
            {                
                const cv::KeyPoint &kpUn = pKFi->mvKeysUn[mit->second];

                // Monocular observation
                if(pKFi->mvuRight[mit->second]<0)
                {
                    Eigen::Matrix<double,2,1> obs;
                    obs << kpUn.pt.x, kpUn.pt.y;

                    g2o::EdgeSE3ProjectXYZ* e = new g2o::EdgeSE3ProjectXYZ();

                    e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                    e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
                    e->setMeasurement(obs);
                    const float &invSigma2 = pKFi->mvInvLevelSigma2[kpUn.octave];
                    e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                    g2o::RobustKernelCauchy* rk = new g2o::RobustKernelCauchy;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuberMonoPoints);

                    e->fx = pKFi->fx;
                    e->fy = pKFi->fy;
                    e->cx = pKFi->cx;
                    e->cy = pKFi->cy;

                    optimizer.addEdge(e);
                    vpEdgesMonoPoints.push_back(e);
                    vpEdgeKFMonoPoints.push_back(pKFi);
                    vpMapPointEdgeMono.push_back(pMP);
		    number++;
		}
	    }
	}
    }
  
    if(pbStopFlag)
        if(*pbStopFlag)
            return;
	
     // Set MapPoint vertices
    const int nExpectedSizeLines = (lLocalKeyFrames.size()+lFixedCameras.size())*lLocalMapLines.size();

    vector<EdgeSE3ProjectXYZLines*> vpEdgesMonoLines;
    vpEdgesMonoLines.reserve(2*nExpectedSizeLines);

    vector<KeyFrame*> vpEdgeKFMonoLines;
    vpEdgeKFMonoLines.reserve(nExpectedSizeLines);

    vector<MapLine*> vpMapLineEdgeMono;
    vpMapLineEdgeMono.reserve(nExpectedSizeLines);	
    
    const float thHuberMonoLines = sqrt(3.841);
	
    for(list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        MapLine* pML = *lit;

        g2o::VertexSBAPointXYZ* vFirPoint = new g2o::VertexSBAPointXYZ();
        vFirPoint->setEstimate(Converter::toVector3d(pML->GetFirPWorldPos()));
        int idFir = 2*pML->mnId+maxMapPointAndKFid+1;
        vFirPoint->setId(idFir);
        vFirPoint->setMarginalized(true);
        optimizer.addVertex(vFirPoint);

        g2o::VertexSBAPointXYZ* vEndPoint = new g2o::VertexSBAPointXYZ();
        vEndPoint->setEstimate(Converter::toVector3d(pML->GetEndPWorldPos()));
        int idEnd = 2*pML->mnId+1+maxMapPointAndKFid+1;
        vEndPoint->setId(idEnd);
        vEndPoint->setMarginalized(true);
        optimizer.addVertex(vEndPoint);

        const map<KeyFrame*,size_t> observations = pML->GetObservations();

        //Set edges
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
	    
            KeyFrame* pKFi = mit->first;   

            if(!pKFi->isBadLines() && !pKFi->isBad() )
            {                
                const KeyLine &kLUn = pKFi->mvLinesUn[mit->second];
		Vector3d sp_l; sp_l << kLUn.startPointX, kLUn.startPointY, 1.0;
		Vector3d ep_l; ep_l << kLUn.endPointX,   kLUn.endPointY,   1.0;
		Vector3d le_l; le_l << sp_l.cross(ep_l);
		//le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) + le_l(2)*le_l(2));
		le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1));
		
		//输入地图线首尾端点的观测值
		Eigen::Matrix<double,3,1> obs;
		obs << le_l(0) , le_l(1), le_l(2);

		EdgeSE3ProjectXYZLines* eFir = new EdgeSE3ProjectXYZLines();

		eFir->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(idFir)));
		eFir->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
		eFir->setMeasurement(obs);
		const float &invSigma2 = pKFi->mvInvLevelSigma2Lines[kLUn.octave];
		eFir->setInformation(Eigen::Matrix3d::Identity()*invSigma2);
		
		g2o::RobustKernelCauchy* rkFir = new g2o::RobustKernelCauchy;
		eFir->setRobustKernel(rkFir);
		rkFir->setDelta(thHuberMonoLines);
		
		eFir->fx = pKFi->fx;
		eFir->fy = pKFi->fy;
		eFir->cx = pKFi->cx;
		eFir->cy = pKFi->cy;
		number++;
		optimizer.addEdge(eFir);
		vpEdgesMonoLines.push_back(eFir);//存入首端点对应边
		
		
		EdgeSE3ProjectXYZLines* eEnd = new EdgeSE3ProjectXYZLines();

		eEnd->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(idEnd)));
		eEnd->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
		eEnd->setMeasurement(obs);
		eEnd->setInformation(Eigen::Matrix3d::Identity()*invSigma2);

		g2o::RobustKernelCauchy* rkEnd = new g2o::RobustKernelCauchy;
		eEnd->setRobustKernel(rkEnd);
		rkEnd->setDelta(thHuberMonoLines);

		eEnd->fx = pKFi->fx;
		eEnd->fy = pKFi->fy;
		eEnd->cx = pKFi->cx;
		eEnd->cy = pKFi->cy;

		number++;		
		optimizer.addEdge(eEnd);
		vpEdgesMonoLines.push_back(eEnd);
			
		vpEdgeKFMonoLines.push_back(pKFi);		    		 	    		    
		vpMapLineEdgeMono.push_back(pML);
		    		    	    		    
	    }
	}
    }   
	
    optimizer.initializeOptimization();
    optimizer.optimize(5);
    
    vector<pair<KeyFrame*,MapPoint*> > vToErasePoints;
    vToErasePoints.reserve(vpEdgesMonoPoints.size());
    for(size_t i=0, iend=vpEdgesMonoPoints.size(); i<iend;i++)
    {
        g2o::EdgeSE3ProjectXYZ* e = vpEdgesMonoPoints[i];
        MapPoint* pMP = vpMapPointEdgeMono[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>5.991 || !e->isDepthPositive())
        {
            KeyFrame* pKFi = vpEdgeKFMonoPoints[i];
            vToErasePoints.push_back(make_pair(pKFi,pMP));
        }
    }
    
    vector<pair<KeyFrame*,MapLine*> > vToEraseLines;
    vToEraseLines.reserve(vpEdgesMonoLines.size());
    
     for(size_t i=0, iend=vpEdgesMonoLines.size(); i<iend;i=i+2)
    {
        EdgeSE3ProjectXYZLines* eFir = vpEdgesMonoLines[i];
	EdgeSE3ProjectXYZLines* eEnd = vpEdgesMonoLines[i+1];
	
        MapLine* pML = vpMapLineEdgeMono[i / 2];
	
	if(pML->isBad())
            continue;
	 if(eFir->chi2()+eEnd->chi2()>5.991 || !eFir->isDepthPositive() || !eEnd->isDepthPositive())
        {
	    KeyFrame* pKFi = vpEdgeKFMonoLines[i / 2];
	    vToEraseLines.push_back(make_pair(pKFi,pML));
	}
    }
    
    // Get Map Mutex
    unique_lock<mutex> lock(pMap->mMutexMapUpdate);
    
    //删除坏点的观测
    if(!vToErasePoints.empty())
    {
        for(size_t i=0;i<vToErasePoints.size();i++)
        {
            KeyFrame* pKFi = vToErasePoints[i].first;
            MapPoint* pMPi = vToErasePoints[i].second;
            pKFi->EraseMapPointMatch(pMPi);
            pMPi->EraseObservation(pKFi);
        }
    }
    
    //删除坏的线观测
    if(!vToEraseLines.empty())
    {
        for(size_t i=0;i<vToEraseLines.size();i++)
        {
            KeyFrame* pKFi = vToEraseLines[i].first;
            MapLine* pMLi = vToEraseLines[i].second;
            pKFi->EraseMapLineMatch(pMLi);
            pMLi->EraseObservation(pKFi);
        }
    }

    // Recover optimized data

    //Keyframes
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKF = *lit;
        g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
        pKF->SetPose(Converter::toCvMat(SE3quat));
    }
    
    //Points
    for(list<MapPoint*>::iterator lit=lLocalMapPoints.begin(), lend=lLocalMapPoints.end(); lit!=lend; lit++)
    {
        MapPoint* pMP = *lit;
        g2o::VertexSBAPointXYZ* vPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(pMP->mnId+maxKFid+1));
        pMP->SetWorldPos(Converter::toCvMat(vPoint->estimate()));
        pMP->UpdateNormalAndDepth();
    }
    
    //Lines
    for(list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        MapLine* pML = *lit;

        g2o::VertexSBAPointXYZ* vFirPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(2*pML->mnId+maxMapPointAndKFid+1));
	cv::Mat FirPWorldPos = Converter::toCvMat(vFirPoint->estimate());
        pML->SetFirPWorldPos(FirPWorldPos);

	g2o::VertexSBAPointXYZ* vEndPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(2*pML->mnId+1+maxMapPointAndKFid+1));
	cv::Mat EndPWorldPos = Converter::toCvMat(vEndPoint->estimate());
        pML->SetEndPWorldPos(EndPWorldPos);

	cv::Mat MidPWorldPos(3,1,CV_32F);
	for(int i=0;i<3;i++)
            MidPWorldPos.at<float>(i)=
            0.5*( FirPWorldPos.at<float>(i)+EndPWorldPos.at<float>(i) );	
	pML->SetMidPWorldPos(MidPWorldPos);

        pML->UpdateNormalAndDepth();
    }  

}

//这里使用，应该是对应于点特征已经彻底丢失，线特征仍能勉强维持的情况
void Optimizer::LocalBundleAdjustmentDoubleLines(KeyFrame* pKF, bool *pbStopFlag, Map *pMap)
{
  // Local KeyFrames: First Breath Search from Current Keyframe
    list<KeyFrame*> lLocalKeyFrames;

    lLocalKeyFrames.push_back(pKF);
    pKF->mnBALocalForKFLines = pKF->mnId;

    const vector<KeyFrame*> vNeighKFs = pKF->GetVectorCovisibleKeyFramesLines();
    for(int i=0, iend=vNeighKFs.size(); i<iend; i++)
    {
        KeyFrame* pKFi = vNeighKFs[i];
        pKFi->mnBALocalForKFLines = pKF->mnId;
        if(!pKFi->isBadLines())
            lLocalKeyFrames.push_back(pKFi);
    }

    // Local MapPoints seen in Local KeyFrames
    list<MapLine*> lLocalMapLines;
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin() , lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        vector<MapLine*> vpMLs = (*lit)->GetMapLineMatches();
        for(vector<MapLine*>::iterator vit=vpMLs.begin(), vend=vpMLs.end(); vit!=vend; vit++)
        {
            MapLine* pML = *vit;
            if(pML)
                if(!pML->isBad())
                    if(pML->mnBALocalForKF!=pKF->mnId)
                    {
                        lLocalMapLines.push_back(pML);
                        pML->mnBALocalForKF=pKF->mnId;
                    }
        }
    }

    // Fixed Keyframes. Keyframes that see Local MapPoints but that are not Local Keyframes
    list<KeyFrame*> lFixedCameras;
    for(list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        map<KeyFrame*,size_t> observations = (*lit)->GetObservations();
        for(map<KeyFrame*,size_t>::iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(pKFi->mnBALocalForKFLines!=pKF->mnId && pKFi->mnBAFixedForKFLines!=pKF->mnId)
            {                
                pKFi->mnBAFixedForKFLines=pKF->mnId;
                if(!pKFi->isBadLines())
                    lFixedCameras.push_back(pKFi);
            }
        }
    }

    // Setup optimizer
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    if(pbStopFlag)
        optimizer.setForceStopFlag(pbStopFlag);

    unsigned long maxKFid = 0;

    // Set Local KeyFrame vertices
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(pKFi->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }

    // Set Fixed KeyFrame vertices
    for(list<KeyFrame*>::iterator lit=lFixedCameras.begin(), lend=lFixedCameras.end(); lit!=lend; lit++)
    {
        KeyFrame* pKFi = *lit;
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKFi->GetPose()));
        vSE3->setId(pKFi->mnId);
        vSE3->setFixed(true);
        optimizer.addVertex(vSE3);
        if(pKFi->mnId>maxKFid)
            maxKFid=pKFi->mnId;
    }

    // Set MapPoint vertices
    const int nExpectedSize = (lLocalKeyFrames.size()+lFixedCameras.size())*lLocalMapLines.size();

    vector<EdgeSE3ProjectXYZLines*> vpEdgesMono;
    vpEdgesMono.reserve(2*nExpectedSize);

    vector<KeyFrame*> vpEdgeKFMono;
    vpEdgeKFMono.reserve(nExpectedSize);

    vector<MapLine*> vpMapLineEdgeMono;
    vpMapLineEdgeMono.reserve(nExpectedSize);

    const float thHuberMono = sqrt(3.841);

    //Set edges
    for(list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        MapLine* pML = *lit;
	//设置地图线首端点
        g2o::VertexSBAPointXYZ* vFirPoint = new g2o::VertexSBAPointXYZ();
        vFirPoint->setEstimate(Converter::toVector3d(pML->GetFirPWorldPos()));
        int idFir = 2*pML->mnId+maxKFid+1;
        vFirPoint->setId(idFir);
        vFirPoint->setMarginalized(true);
        optimizer.addVertex(vFirPoint);
	//设置地图线尾端点
        g2o::VertexSBAPointXYZ* vEndPoint = new g2o::VertexSBAPointXYZ();
        vEndPoint->setEstimate(Converter::toVector3d(pML->GetEndPWorldPos()));
        int idEnd = 2*pML->mnId+1+maxKFid+1;
        vEndPoint->setId(idEnd);
        vEndPoint->setMarginalized(true);
        optimizer.addVertex(vEndPoint);

        const map<KeyFrame*,size_t> observations = pML->GetObservations();

        //Set edges
        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(), mend=observations.end(); mit!=mend; mit++)
        {
            KeyFrame* pKFi = mit->first;

            if(!pKFi->isBadLines())
            {                
                const KeyLine &kLUn = pKFi->mvLinesUn[mit->second];
		Vector3d sp_l; sp_l << kLUn.startPointX, kLUn.startPointY, 1.0;
		Vector3d ep_l; ep_l << kLUn.endPointX,   kLUn.endPointY,   1.0;
		Vector3d le_l; le_l << sp_l.cross(ep_l);
		//le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) + le_l(2)*le_l(2));
		le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1));

		Eigen::Matrix<double,3,1> obs;
		obs << le_l(0) , le_l(1), le_l(2);

		EdgeSE3ProjectXYZLines* eFir = new EdgeSE3ProjectXYZLines();

		eFir->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(idFir)));
		eFir->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
		eFir->setMeasurement(obs);
		const float &invSigma2 = pKFi->mvInvLevelSigma2Lines[kLUn.octave];
		eFir->setInformation(Eigen::Matrix3d::Identity()*invSigma2);

		g2o::RobustKernelCauchy* rkFir = new g2o::RobustKernelCauchy;
		eFir->setRobustKernel(rkFir);
		rkFir->setDelta(thHuberMono);

		eFir->fx = pKFi->fx;
		eFir->fy = pKFi->fy;
		eFir->cx = pKFi->cx;
		eFir->cy = pKFi->cy;

		optimizer.addEdge(eFir);
		vpEdgesMono.push_back(eFir);
			
		EdgeSE3ProjectXYZLines* eEnd = new EdgeSE3ProjectXYZLines();

		eEnd->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(idEnd)));
		eEnd->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
		eEnd->setMeasurement(obs);
		eEnd->setInformation(Eigen::Matrix3d::Identity()*invSigma2);

		g2o::RobustKernelCauchy* rkEnd = new g2o::RobustKernelCauchy;
		eEnd->setRobustKernel(rkEnd);
		rkEnd->setDelta(thHuberMono);

		eEnd->fx = pKFi->fx;
		eEnd->fy = pKFi->fy;
		eEnd->cx = pKFi->cx;
		eEnd->cy = pKFi->cy;

		optimizer.addEdge(eEnd);
		vpEdgesMono.push_back(eEnd);
			
		vpEdgeKFMono.push_back(pKFi);		    		 	    		    
		vpMapLineEdgeMono.push_back(pML);		    		    	    		    
	    }
	}
    }   

    if(pbStopFlag)
        if(*pbStopFlag)
            return;

    optimizer.initializeOptimization();
    optimizer.optimize(5);

    bool bDoMore= true;

    if(pbStopFlag)
        if(*pbStopFlag)
            bDoMore = false;

    if(bDoMore)
    {
	// Check inlier observations
	for(size_t i=0, iend=vpEdgesMono.size(); i<iend; i = i+2)
	{
	    EdgeSE3ProjectXYZLines* eFir = vpEdgesMono[i];
	    EdgeSE3ProjectXYZLines* eEnd = vpEdgesMono[i+1];
	    
	    MapLine* pML = vpMapLineEdgeMono[i];

	    if(pML->isBad())
		continue;

	    if(eFir->chi2()+eEnd->chi2()>5.991 || !eFir->isDepthPositive() || !eEnd->isDepthPositive())
	    {
		eFir->setLevel(1);
		eEnd->setLevel(1);
	    }

	    eFir->setRobustKernel(0);
	    eEnd->setRobustKernel(0);
	}
	// Optimize again without the outliers

	optimizer.initializeOptimization(0);
	optimizer.optimize(10);
    }

    vector<pair<KeyFrame*,MapLine*> > vToErase;
    vToErase.reserve(vpEdgesMono.size());

    // Check inlier observations       
    for(size_t i=0, iend=vpEdgesMono.size(); i<iend;i=i+2)
    {
        EdgeSE3ProjectXYZLines* eFir = vpEdgesMono[i];
	EdgeSE3ProjectXYZLines* eEnd = vpEdgesMono[i+1];
	
        MapLine* pML = vpMapLineEdgeMono[i / 2];
	
        if(pML->isBad())
            continue;

        if(eFir->chi2()+eEnd->chi2()>5.991 || !eFir->isDepthPositive() || !eEnd->isDepthPositive())
        {
	    KeyFrame* pKFi = vpEdgeKFMono[i / 2];
            vToErase.push_back(make_pair(pKFi,pML));
        }

    }
    // Get Map Mutex
    unique_lock<mutex> lock(pMap->mMutexMapUpdate);

    if(!vToErase.empty())
    {
        for(size_t i=0;i<vToErase.size();i++)
        {
            KeyFrame* pKFi = vToErase[i].first;
            MapLine* pMLi = vToErase[i].second;
            pKFi->EraseMapLineMatch(pMLi);
            pMLi->EraseObservation(pKFi);
        }
    }

    // Recover optimized data

    //Keyframes
    for(list<KeyFrame*>::iterator lit=lLocalKeyFrames.begin(), lend=lLocalKeyFrames.end(); lit!=lend; lit++)
    {
        KeyFrame* pKF = *lit;
        g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
        pKF->SetPose(Converter::toCvMat(SE3quat));
    }

    //Lines
    for(list<MapLine*>::iterator lit=lLocalMapLines.begin(), lend=lLocalMapLines.end(); lit!=lend; lit++)
    {
        MapLine* pML = *lit;

        g2o::VertexSBAPointXYZ* vFirPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(2*pML->mnId+maxKFid+1));
	cv::Mat FirPWorldPos = Converter::toCvMat(vFirPoint->estimate());
        pML->SetFirPWorldPos(FirPWorldPos);

	g2o::VertexSBAPointXYZ* vEndPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(2*pML->mnId+1+maxKFid+1));
	cv::Mat EndPWorldPos = Converter::toCvMat(vEndPoint->estimate());
        pML->SetEndPWorldPos(EndPWorldPos);
	
	cv::Mat MidPWorldPos(3,1,CV_32F);
	for(int i=0;i<3;i++)
            MidPWorldPos.at<float>(i)=
            0.5*( FirPWorldPos.at<float>(i)+EndPWorldPos.at<float>(i) );	
	pML->SetMidPWorldPos(MidPWorldPos);

        pML->UpdateNormalAndDepth();
    }
}

// 单目初始化：（系统运行只使用一次）
void Optimizer::GlobalBundleAdjustemntIni(Map* pMap, int nIterations, const bool bRobust)
{

    vector<KeyFrame*> vpKFs = pMap->GetAllKeyFrames();
    vector<MapPoint*> vpMP = pMap->GetAllMapPoints();
    vector<MapLine*> vpML = pMap->GetAllMapLines();

    float PoseErrorForPoints;
    float PoseErrorForLines;
    
    //cout << "debug 分线程优化数据准备成功: " << true << endl << endl; 
    
    //（2）分线程：调用仅点和仅线的BA子函数对前两帧关键帧进行优化。
    thread threadBundleAdjustmentPointsIni(&Optimizer::BundleAdjustmentPointsIni, ref(vpKFs), ref(vpMP), 5, true, ref(PoseErrorForPoints) );                      
    thread threadBundleAdjustmentLinesIni(&Optimizer::BundleAdjustmentLinesIni, ref(vpKFs), ref(vpML), 5, true, ref(PoseErrorForLines) );
    threadBundleAdjustmentPointsIni.join();
    threadBundleAdjustmentLinesIni.join();
    
    vector<KeyFrame*> vpKFs2 = pMap->GetAllKeyFrames(); //实际这里应该没有发生改变
    vector<MapPoint*> vpMP2 = pMap->GetAllMapPoints();
    vector<MapLine*> vpML2 = pMap->GetAllMapLines();
    
    //根据单位误差设置当前帧的更优位姿
    for(vector<KeyFrame*>::iterator vit = vpKFs2.begin(), vend=vpKFs2.end(); vit!=vend; vit++)
    {
	KeyFrame* pKFi = *vit;
	
	if(pKFi->mnId !=0)
	{
	    if(PoseErrorForPoints<=PoseErrorForLines)
	    {
		pKFi->SetPose(pKFi->TcwPoints);
	    }
	    else
	    {
		pKFi->SetPose(pKFi->TcwLines);
	    }	
	}	  
    }
    
    BundleAdjustmentBothIni(vpKFs2,vpMP2,vpML2,5,true);
    
}

void Optimizer::BundleAdjustmentPointsIni(std::vector<KeyFrame*> &vpKFs, std::vector<MapPoint*> &vpMP, 
					  int nIterations, const bool bRobust, float &PoseErrorForPoints)
{
    PoseErrorForPoints = 0.0;
  
    const int nExpectedSize = vpKFs.size()*vpMP.size();	
    
    //用于处理将要删除的坏地图点
    vector<g2o::EdgeSE3ProjectXYZ*> vpEdgesMono;
    vpEdgesMono.reserve(nExpectedSize);

    vector<KeyFrame*> vpEdgeKFMono;
    vpEdgeKFMono.reserve(nExpectedSize);

    vector<MapPoint*> vpMapPointEdgeMono;
    vpMapPointEdgeMono.reserve(nExpectedSize);

    //初始化g2o优化器
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);
    
    //记录最大关键帧
    long unsigned int maxKFid = 0;

    // Set KeyFrame vertices
    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];

        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKF->GetPose()));
        vSE3->setId(pKF->mnId);
        vSE3->setFixed(pKF->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKF->mnId>maxKFid)
            maxKFid=pKF->mnId;
    }
    
    const float thHuber2D = sqrt(5.991);
    // Set MapPoint vertices
    for(size_t i=0; i<vpMP.size(); i++)
    {
        MapPoint* pMP = vpMP[i];
	
        g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
        vPoint->setEstimate(Converter::toVector3d(pMP->GetWorldPos()));
        const int id = pMP->mnId+maxKFid+1;
        vPoint->setId(id);
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

       const map<KeyFrame*,size_t> observations = pMP->GetObservations();

        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(); mit!=observations.end(); mit++)
        {

            KeyFrame* pKF = mit->first;

            const cv::KeyPoint &kpUn = pKF->mvKeysUn[mit->second];

            if(pKF->mvuRight[mit->second]<0)
            {
                Eigen::Matrix<double,2,1> obs;
                obs << kpUn.pt.x, kpUn.pt.y;

                g2o::EdgeSE3ProjectXYZ* e = new g2o::EdgeSE3ProjectXYZ();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKF->mnId)));
                e->setMeasurement(obs);
                const float &invSigma2 = pKF->mvInvLevelSigma2[kpUn.octave];
                e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                if(bRobust)
                {
                    g2o::RobustKernelCauchy* rk = new g2o::RobustKernelCauchy;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuber2D);
                }

                e->fx = pKF->fx;
                e->fy = pKF->fy;
                e->cx = pKF->cx;
                e->cy = pKF->cy;

                optimizer.addEdge(e);
		vpEdgesMono.push_back(e);
		vpEdgeKFMono.push_back(pKF);
                vpMapPointEdgeMono.push_back(pMP);
            }
	}
	
    }
       
    //优化
    optimizer.initializeOptimization();
    optimizer.optimize(nIterations);
    
    
    vector<pair<KeyFrame*,MapPoint*> > vToErase;
    vToErase.reserve(vpEdgesMono.size());
    
    int nInlierKFNum = 0;
    // Check inlier observations       
    for(size_t i=0, iend=vpEdgesMono.size(); i<iend;i++)
    {
        g2o::EdgeSE3ProjectXYZ* e = vpEdgesMono[i];
        MapPoint* pMP = vpMapPointEdgeMono[i];
	
	KeyFrame* pKFi = vpEdgeKFMono[i];
        if(e->chi2()>5.991 || !e->isDepthPositive())
        {
            
            vToErase.push_back(make_pair(pKFi,pMP));
        }
        else
	{
	    if(pKFi->mnId !=0)
	    {
		//计算当前关键帧的总有效误差
		PoseErrorForPoints += e->chi2();
		nInlierKFNum++;
	    }
	}	  
    }
    //计算当前关键帧的单位有效误差
    PoseErrorForPoints = PoseErrorForPoints / nInlierKFNum;

    if(!vToErase.empty())
    {
        for(size_t i=0;i<vToErase.size();i++)
        {
            KeyFrame* pKFi = vToErase[i].first;
            MapPoint* pMPi = vToErase[i].second;
            pKFi->EraseMapPointMatch(pMPi);
            pMPi->EraseObservation(pKFi);
        }
    }
    
    cout << "PoseErrorForPoints:" << PoseErrorForPoints << endl<< endl<< endl;
    
    // Recover optimized data

    //Keyframes
    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];
	
	g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
	
	if(pKF->mnId == 0)
	{
	    pKF->SetPose(Converter::toCvMat(SE3quat));
	}
	else
	{
	    pKF->SetPosePoints(Converter::toCvMat(SE3quat));
	}
    }
    
    //Points
    for(size_t i=0; i<vpMP.size(); i++)
    {

        MapPoint* pMP = vpMP[i];

         if(pMP->isBad())
             continue;
        g2o::VertexSBAPointXYZ* vPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(pMP->mnId+maxKFid+1));
       
        pMP->SetWorldPos(Converter::toCvMat(vPoint->estimate()));  
    }
}

void Optimizer::BundleAdjustmentLinesIni(std::vector<KeyFrame*> &vpKFs,std::vector<MapLine*> &vpML,
				     int nIterations, const bool bRobust, float &PoseErrorForLines)
{
    PoseErrorForLines = 0.0;
    
    const int nExpectedSize = vpKFs.size()*vpML.size();	
    
    //用于处理将要删除的坏地图点
    vector<EdgeSE3ProjectXYZLines*> vpEdgesMono;
    vpEdgesMono.reserve(2*nExpectedSize);

    vector<KeyFrame*> vpEdgeKFMono;
    vpEdgeKFMono.reserve(nExpectedSize);

    vector<MapLine*> vpMapLineEdgeMono;
    vpMapLineEdgeMono.reserve(nExpectedSize);

    //初始化g2o优化器
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);
    
    //记录最大关键帧
    long unsigned int maxKFid = 0;

    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];

        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKF->GetPose()));
        vSE3->setId(pKF->mnId);
        vSE3->setFixed(pKF->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKF->mnId>maxKFid)
            maxKFid=pKF->mnId;
    }
    
//     float PoseErrorForLinesIni = 0.0;
//     int LinesNum = 0;
    
    const float thHuber2D = sqrt(3.841);
    // Set MapPoint vertices
    for(size_t i=0; i<vpML.size(); i++)
    {
        MapLine* pML = vpML[i];

        g2o::VertexSBAPointXYZ* vFirPoint = new g2o::VertexSBAPointXYZ();
        vFirPoint->setEstimate(Converter::toVector3d(pML->GetFirPWorldPos()));
        int idFir = 2*pML->mnId+maxKFid+1;
        vFirPoint->setId(idFir);
        vFirPoint->setMarginalized(true);
        optimizer.addVertex(vFirPoint);

        g2o::VertexSBAPointXYZ* vEndPoint = new g2o::VertexSBAPointXYZ();
        vEndPoint->setEstimate(Converter::toVector3d(pML->GetEndPWorldPos()));
        int idEnd = 2*pML->mnId+1+maxKFid+1;
        vEndPoint->setId(idEnd);
        vEndPoint->setMarginalized(true);
        optimizer.addVertex(vEndPoint);

        const map<KeyFrame*,size_t> observations = pML->GetObservations();

        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(); mit!=observations.end(); mit++)
        {

            KeyFrame* pKFi = mit->first;
	    const KeyLine &kLUn = pKFi->mvLinesUn[mit->second];

	    Vector3d sp_l; sp_l << kLUn.startPointX, kLUn.startPointY, 1.0;
	    Vector3d ep_l; ep_l << kLUn.endPointX,   kLUn.endPointY,   1.0;
	    Vector3d le_l; le_l << sp_l.cross(ep_l);
	    //le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) + le_l(2)*le_l(2));
	    le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1));

	    Eigen::Matrix<double,3,1> obs;
	    obs << le_l(0) , le_l(1), le_l(2);
	    
	    EdgeSE3ProjectXYZLines* eFir = new EdgeSE3ProjectXYZLines();

	    eFir->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(idFir)));
	    eFir->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
	    eFir->setMeasurement(obs);
	    const float invSigma2 = pKFi->mvInvLevelSigma2Lines[kLUn.octave];
	    eFir->setInformation(Eigen::Matrix3d::Identity()*invSigma2);
	    
	    if(bRobust)
	    {
		g2o::RobustKernelCauchy* rkFir = new g2o::RobustKernelCauchy;
		eFir->setRobustKernel(rkFir);
		rkFir->setDelta(thHuber2D);
	    }
	    
	    eFir->fx = pKFi->fx;
	    eFir->fy = pKFi->fy;
	    eFir->cx = pKFi->cx;
	    eFir->cy = pKFi->cy;

	    optimizer.addEdge(eFir);
	    vpEdgesMono.push_back(eFir);
	    	    
	    EdgeSE3ProjectXYZLines* eEnd = new EdgeSE3ProjectXYZLines();

	    eEnd->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(idEnd)));
	    eEnd->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
	    eEnd->setMeasurement(obs);
	    eEnd->setInformation(Eigen::Matrix3d::Identity()*invSigma2);

	    if(bRobust)
	    {
		g2o::RobustKernelCauchy* rkEnd = new g2o::RobustKernelCauchy;
		eEnd->setRobustKernel(rkEnd);
		rkEnd->setDelta(thHuber2D);
	    }

	    eEnd->fx = pKFi->fx;
	    eEnd->fy = pKFi->fy;
	    eEnd->cx = pKFi->cx;
	    eEnd->cy = pKFi->cy;

	    optimizer.addEdge(eEnd);
	    vpEdgesMono.push_back(eEnd);
		    
	    vpEdgeKFMono.push_back(pKFi);		    		 	    		    
	    vpMapLineEdgeMono.push_back(pML);  
	    
	    
// 	    PoseErrorForLinesIni += eFir->chi2();
// 	    PoseErrorForLinesIni += eEnd->chi2();
// 	    LinesNum ++;
	}
	
    }
    
//     PoseErrorForLinesIni /= LinesNum;
//     cout << "PoseErrorForLinesIni:" << PoseErrorForLinesIni << endl<< endl<< endl;
//     cout << "LinesNum:" << LinesNum << endl<< endl<< endl;
    
    //优化
    optimizer.initializeOptimization();
    optimizer.optimize(nIterations); 
    
    vector<pair<KeyFrame*,MapLine*> > vToErase;
    vToErase.reserve(vpEdgesMono.size());
    
    int nInlierKFNum = 0;
    // Check inlier observations       
    for(size_t i=0, iend=vpEdgesMono.size(); i<iend; i=i+2 )
    {
        EdgeSE3ProjectXYZLines* eFir = vpEdgesMono[i];
	EdgeSE3ProjectXYZLines* eEnd = vpEdgesMono[i+1];
	
        MapLine* pML = vpMapLineEdgeMono[i / 2];
	
	KeyFrame* pKFi = vpEdgeKFMono[i / 2];
        if(eFir->chi2()+eEnd->chi2()>5.991 || !eFir->isDepthPositive() || !eEnd->isDepthPositive())
        {
            
            vToErase.push_back(make_pair(pKFi,pML));
        }
        else
	{
	    if(pKFi->mnId !=0)
	    {
		//计算当前关键帧的总有效误差
		float LineError = eFir->chi2()+eEnd->chi2();
		PoseErrorForLines += LineError;
		nInlierKFNum++;
	    }
	}	  
    }

    PoseErrorForLines = PoseErrorForLines / nInlierKFNum;
    
    cout << "PoseErrorForLines:" << PoseErrorForLines << endl<< endl<< endl;
    //cout << "nInlierKFNum:" << nInlierKFNum << endl<< endl<< endl;

    if(!vToErase.empty())
    {
        for(size_t i=0;i<vToErase.size();i++)
        {
            KeyFrame* pKFi = vToErase[i].first;
            MapLine* pMLi = vToErase[i].second;
            pKFi->EraseMapLineMatch(pMLi);
            pMLi->EraseObservation(pKFi);
        }
    }
    
    // Recover optimized data

    //Keyframes
    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];
	
	g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
	
	if(pKF->mnId == 0)
	{
	    pKF->SetPose(Converter::toCvMat(SE3quat));
	}
	else
	{
	    pKF->SetPoseLines(Converter::toCvMat(SE3quat));
	}
    }
    
    //Lines
    for(size_t i=0; i<vpML.size(); i++)
    {

        MapLine* pML = vpML[i];

         if(pML->isBad())
             continue;
	
        g2o::VertexSBAPointXYZ* vFirPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(2*pML->mnId+maxKFid+1));
	cv::Mat FirPWorldPos = Converter::toCvMat(vFirPoint->estimate());
        pML->SetFirPWorldPos(FirPWorldPos);
	
	g2o::VertexSBAPointXYZ* vEndPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(2*pML->mnId+1+maxKFid+1));
	cv::Mat EndPWorldPos = Converter::toCvMat(vEndPoint->estimate());
        pML->SetEndPWorldPos(EndPWorldPos);

	cv::Mat MidPWorldPos(3,1,CV_32F);
	for(int i=0;i<3;i++)
            MidPWorldPos.at<float>(i)=
            0.5*( FirPWorldPos.at<float>(i)+EndPWorldPos.at<float>(i) );	
	pML->SetMidPWorldPos(MidPWorldPos);
        
    }
}

//初始化分级优化策略的第二步，点线特征同时优化
void Optimizer::BundleAdjustmentBothIni(std::vector<KeyFrame*> &vpKFs,std::vector<MapPoint*> &vpMP ,
					std::vector<MapLine*> &vpML, int nIterations, const bool bRobust)
{
    const int nExpectedSizePoints = vpKFs.size()*vpMP.size();	
    
    //用于处理将要删除的坏地图点
    vector<g2o::EdgeSE3ProjectXYZ*> vpEdgesMonoPoints;
    vpEdgesMonoPoints.reserve(nExpectedSizePoints);

    vector<KeyFrame*> vpEdgeKFMonoPoints;
    vpEdgeKFMonoPoints.reserve(nExpectedSizePoints);

    vector<MapPoint*> vpMapPointEdgeMono;
    vpMapPointEdgeMono.reserve(nExpectedSizePoints);
    
    const int nExpectedSizeLines = vpKFs.size()*vpML.size();	
    
    //用于处理将要删除的坏地图点
    vector<EdgeSE3ProjectXYZLines*> vpEdgesMonoLines;
    vpEdgesMonoLines.reserve(nExpectedSizeLines);

    vector<KeyFrame*> vpEdgeKFMonoLines;
    vpEdgeKFMonoLines.reserve(nExpectedSizeLines);

    vector<MapLine*> vpMapLineEdgeMono;
    vpMapLineEdgeMono.reserve(nExpectedSizeLines);
    
    //初始化g2o优化器
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);

    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);
    

    long unsigned int maxKFid = 0;
    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];
	
        g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
        vSE3->setEstimate(Converter::toSE3Quat(pKF->GetPose()));
        vSE3->setId(pKF->mnId);
        vSE3->setFixed(pKF->mnId==0);
        optimizer.addVertex(vSE3);
        if(pKF->mnId>maxKFid)
            maxKFid=pKF->mnId;
    }
    const float thHuber2DPoints = sqrt(5.991);
    const float thHuber2DLines = sqrt(3.841);
    
    unsigned long maxMapPointAndKFid = 0;
    
    for(size_t i=0; i<vpMP.size(); i++)
    {
        MapPoint* pMP = vpMP[i];
	
        g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
        vPoint->setEstimate(Converter::toVector3d(pMP->GetWorldPos()));
        
	unsigned long id = pMP->mnId+maxKFid+1;
        vPoint->setId(id);
	if(id>maxMapPointAndKFid)
            maxMapPointAndKFid=id;
	
        vPoint->setMarginalized(true);
        optimizer.addVertex(vPoint);

       const map<KeyFrame*,size_t> observations = pMP->GetObservations();

        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(); mit!=observations.end(); mit++)
        {

            KeyFrame* pKF = mit->first;

            const cv::KeyPoint &kpUn = pKF->mvKeysUn[mit->second];

            if(pKF->mvuRight[mit->second]<0)
            {
                Eigen::Matrix<double,2,1> obs;
                obs << kpUn.pt.x, kpUn.pt.y;

                g2o::EdgeSE3ProjectXYZ* e = new g2o::EdgeSE3ProjectXYZ();

                e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
                e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKF->mnId)));
                e->setMeasurement(obs);
                const float &invSigma2 = pKF->mvInvLevelSigma2[kpUn.octave];
                e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

                if(bRobust)
                {
                    g2o::RobustKernelCauchy* rk = new g2o::RobustKernelCauchy;
                    e->setRobustKernel(rk);
                    rk->setDelta(thHuber2DPoints);
                }

                e->fx = pKF->fx;
                e->fy = pKF->fy;
                e->cx = pKF->cx;
                e->cy = pKF->cy;

                optimizer.addEdge(e);
		vpEdgesMonoPoints.push_back(e);
		vpEdgeKFMonoPoints.push_back(pKF);
                vpMapPointEdgeMono.push_back(pMP);
            }
	}
    }
    
    for(size_t i=0; i<vpML.size(); i++)
    {
        MapLine* pML = vpML[i];

        g2o::VertexSBAPointXYZ* vFirPoint = new g2o::VertexSBAPointXYZ();
        vFirPoint->setEstimate(Converter::toVector3d(pML->GetFirPWorldPos()));
        int idFir = 2*pML->mnId+maxMapPointAndKFid+1;
        vFirPoint->setId(idFir);
        vFirPoint->setMarginalized(true);
        optimizer.addVertex(vFirPoint);
	
        g2o::VertexSBAPointXYZ* vEndPoint = new g2o::VertexSBAPointXYZ();
        vEndPoint->setEstimate(Converter::toVector3d(pML->GetEndPWorldPos()));
        int idEnd = 2*pML->mnId+1+maxMapPointAndKFid+1;
        vEndPoint->setId(idEnd);
        vEndPoint->setMarginalized(true);
        optimizer.addVertex(vEndPoint);

        const map<KeyFrame*,size_t> observations = pML->GetObservations();

        for(map<KeyFrame*,size_t>::const_iterator mit=observations.begin(); mit!=observations.end(); mit++)
        {

            KeyFrame* pKFi = mit->first;
	    const KeyLine &kLUn = pKFi->mvLinesUn[mit->second];
	    Vector3d sp_l; sp_l << kLUn.startPointX, kLUn.startPointY, 1.0;
	    Vector3d ep_l; ep_l << kLUn.endPointX,   kLUn.endPointY,   1.0;
	    Vector3d le_l; le_l << sp_l.cross(ep_l);
	    //le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1) + le_l(2)*le_l(2));
	    le_l = le_l / std::sqrt( le_l(0)*le_l(0) + le_l(1)*le_l(1));

	    Eigen::Matrix<double,3,1> obs;
	    obs << le_l(0) , le_l(1), le_l(2);

	    EdgeSE3ProjectXYZLines* eFir = new EdgeSE3ProjectXYZLines();

	    eFir->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(idFir)));
	    eFir->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
	    eFir->setMeasurement(obs);
	    const float invSigma2 = pKFi->mvInvLevelSigma2Lines[kLUn.octave];
	    eFir->setInformation(Eigen::Matrix3d::Identity()*invSigma2);
	    
	    if(bRobust)
	    {
		g2o::RobustKernelCauchy* rkFir = new g2o::RobustKernelCauchy;
		eFir->setRobustKernel(rkFir);
		rkFir->setDelta(thHuber2DLines);
	    }
	    
	    eFir->fx = pKFi->fx;
	    eFir->fy = pKFi->fy;
	    eFir->cx = pKFi->cx;
	    eFir->cy = pKFi->cy;

	    optimizer.addEdge(eFir);
	    vpEdgesMonoLines.push_back(eFir);
	    
	    
	    EdgeSE3ProjectXYZLines* eEnd = new EdgeSE3ProjectXYZLines();

	    eEnd->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(idEnd)));
	    eEnd->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKFi->mnId)));
	    eEnd->setMeasurement(obs);
	    //const float invSigma2 = pKFi->mvInvLevelSigma2Lines[kLUn.octave];
	    eEnd->setInformation(Eigen::Matrix3d::Identity()*invSigma2);

	    if(bRobust)
	    {
		g2o::RobustKernelCauchy* rkEnd = new g2o::RobustKernelCauchy;
		eEnd->setRobustKernel(rkEnd);
		rkEnd->setDelta(thHuber2DLines);
	    }

	    eEnd->fx = pKFi->fx;
	    eEnd->fy = pKFi->fy;
	    eEnd->cx = pKFi->cx;
	    eEnd->cy = pKFi->cy;

	    optimizer.addEdge(eEnd);
	    vpEdgesMonoLines.push_back(eEnd);
		    
	    vpEdgeKFMonoLines.push_back(pKFi);		    		 	    		    
	    vpMapLineEdgeMono.push_back(pML);
            
	}
    }
    
    optimizer.initializeOptimization();
    optimizer.optimize(nIterations);
    
    vector<pair<KeyFrame*,MapPoint*> > vToErasePoints;
    vToErasePoints.reserve(vpEdgesMonoPoints.size());
    for(size_t i=0, iend=vpEdgesMonoPoints.size(); i<iend;i++)
    {
        g2o::EdgeSE3ProjectXYZ* e = vpEdgesMonoPoints[i];
        MapPoint* pMP = vpMapPointEdgeMono[i];

        if(pMP->isBad())
            continue;

        if(e->chi2()>5.991 || !e->isDepthPositive())
        {
            KeyFrame* pKFi = vpEdgeKFMonoPoints[i];
            vToErasePoints.push_back(make_pair(pKFi,pMP));
        }
    }
    
    vector<pair<KeyFrame*,MapLine*> > vToEraseLines;
    vToEraseLines.reserve(vpEdgesMonoLines.size());
    
     for(size_t i=0, iend=vpEdgesMonoLines.size(); i<iend;i=i+2)
    {
        EdgeSE3ProjectXYZLines* eFir = vpEdgesMonoLines[i];
	EdgeSE3ProjectXYZLines* eEnd = vpEdgesMonoLines[i+1];
	
        MapLine* pML = vpMapLineEdgeMono[i / 2];
	
	if(pML->isBad())
            continue;
	 if(eFir->chi2()+eEnd->chi2()>5.991 || !eFir->isDepthPositive() || !eEnd->isDepthPositive())
        {
	    KeyFrame* pKFi = vpEdgeKFMonoLines[i / 2];
	    vToEraseLines.push_back(make_pair(pKFi,pML));
	}
    }

    if(!vToErasePoints.empty())
    {
        for(size_t i=0;i<vToErasePoints.size();i++)
        {
            KeyFrame* pKFi = vToErasePoints[i].first;
            MapPoint* pMPi = vToErasePoints[i].second;
            pKFi->EraseMapPointMatch(pMPi);
            pMPi->EraseObservation(pKFi);
        }
    }

    if(!vToEraseLines.empty())
    {
        for(size_t i=0;i<vToEraseLines.size();i++)
        {
            KeyFrame* pKFi = vToEraseLines[i].first;
            MapLine* pMLi = vToEraseLines[i].second;
            pKFi->EraseMapLineMatch(pMLi);
            pMLi->EraseObservation(pKFi);
        }
    }
    
    // Recover optimized data

    //Keyframes
    for(size_t i=0; i<vpKFs.size(); i++)
    {
        KeyFrame* pKF = vpKFs[i];
	
	g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
        g2o::SE3Quat SE3quat = vSE3->estimate();
	
	pKF->SetPose(Converter::toCvMat(SE3quat));	
    }
    
    //Points
    for(size_t i=0; i<vpMP.size(); i++)
    {
        MapPoint* pMP = vpMP[i];

         if(pMP->isBad())
             continue;
        g2o::VertexSBAPointXYZ* vPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(pMP->mnId+maxKFid+1));
       
        pMP->SetWorldPos(Converter::toCvMat(vPoint->estimate()));
        pMP->UpdateNormalAndDepth();
        
    }
    
    //Lines
    for(size_t i=0; i<vpML.size(); i++)
    {

        MapLine* pML = vpML[i];

         if(pML->isBad())
             continue;

        g2o::VertexSBAPointXYZ* vFirPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(2*pML->mnId+maxMapPointAndKFid+1));
	cv::Mat FirPWorldPos = Converter::toCvMat(vFirPoint->estimate());
        pML->SetFirPWorldPos(FirPWorldPos);

	g2o::VertexSBAPointXYZ* vEndPoint = static_cast<g2o::VertexSBAPointXYZ*>(optimizer.vertex(2*pML->mnId+1+maxMapPointAndKFid+1));
	cv::Mat EndPWorldPos = Converter::toCvMat(vEndPoint->estimate());
        pML->SetEndPWorldPos(EndPWorldPos);

	cv::Mat MidPWorldPos(3,1,CV_32F);
	for(int i=0;i<3;i++)
            MidPWorldPos.at<float>(i)=
            0.5*( FirPWorldPos.at<float>(i)+EndPWorldPos.at<float>(i) );	
	pML->SetMidPWorldPos(MidPWorldPos);

        pML->UpdateNormalAndDepth();
        
    }
}

} //namespace ORB_SLAM
