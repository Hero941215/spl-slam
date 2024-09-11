echo "Building ROS nodes"

cd Examples/ROS/PL-SLAM
mkdir build
cd build
cmake .. 
make -lpthread
