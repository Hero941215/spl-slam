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

#ifndef FRAMEDRAWER_H
#define FRAMEDRAWER_H

#include "Tracking.h"
#include "MapPoint.h"
#include "MapLine.h"
#include "Map.h"

#include<opencv2/core/core.hpp>
#include<opencv2/features2d/features2d.hpp>

#include<mutex>


namespace PL_SLAM
{

class Tracking;
class Viewer;

class FrameDrawer
{
public:
    FrameDrawer(Map* pMap);

    // Update info from the last processed frame.
    void Update(Tracking *pTracker);
    
    void UpdateBoth(Tracking *pTracker);

    // Draw last processed frame.
    cv::Mat DrawFrame();
    
    cv::Mat DrawFrameBoth();

protected:

    void DrawTextInfo(cv::Mat &im, int nState, cv::Mat &imText);
    
    void DrawTextInfoBoth(cv::Mat &im, int nState, cv::Mat &imText);

    // Info of the frame to be drawn
    cv::Mat mIm;
    int N;
    int NL;
    vector<cv::KeyPoint> mvCurrentKeys;
    vector<KeyLine> mvCurrentKeyLines;
    vector<bool> mvbMap, mvbVO;
    vector<bool> mvbMapLines, mvbVOLines;
    
    bool mbOnlyTracking;
    
    int mnTracked, mnTrackedVO;
    int mnTrackedLines, mnTrackedVOLines;
    
    vector<cv::KeyPoint> mvIniKeys;
    vector<KeyLine> mvIniKeyLines;
    
    vector<int> mvIniMatches;
    vector<int> mvIniMatchesLines;
    
    int mState;

    Map* mpMap;

    std::mutex mMutex;
};

} //namespace ORB_SLAM

#endif // FRAMEDRAWER_H
