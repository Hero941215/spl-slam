# spl-slam
**A visual SLAM using both point and line features. The system is built on the basis of [orb-slam2](https://github.com/raulmur/ORB_SLAM2)**

** **

## 0. Features
The system uses line midpoint information to improve several modules, including front-end line feature matching, system initialization, and so on.

## 1. Prerequisites

Prequisites are identical to [orb-slam2](https://github.com/raulmur/ORB_SLAM2). The library related to line features has been integrated into Thirdparty and can be used directly.

## 3. Run
### 3.1. **run with hilti dataset**

roslaunch robotrun tight_slam_ouster_indoor.launch (hilti-2021)

roslaunch robotrun tight_slam_pandar_indoor.launch (hilti-2022)

### 3.2. **run with UrbanNav dataset**

roslaunch robotrun tight_slam_velodyne_outdoor.launch

### 3.3. **run with livox mid-360**

roslaunch robotrun tight_slam_mid360_indoor.launch

## 4. Acknowledgments

Thanks for LOAM(J. Zhang and S. Singh. LOAM: Lidar Odometry and Mapping in Real-time), [FAST-LIO2](https://github.com/hku-mars/FAST_LIO), [BALM2](https://github.com/hku-mars/BALM).
