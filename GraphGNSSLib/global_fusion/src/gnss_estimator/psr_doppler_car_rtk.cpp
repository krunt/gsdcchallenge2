/*******************************************************
 * Copyright (C) 2019, Intelligent Positioning and Navigation Lab, Hong Kong Polytechnic University
 * 
 * This file is part of GraphGNSSLib.
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Weisong Wen (weisong.wen@connect.polyu.hk)
 * 
 * Main fucntions: pseudorange/Doppler/Carrier-phase fusion using factor graph optimization 
 * input: pseudorange (psr_), Doppler velocity (doppler_) and carrier-phase (car_) from GPS/BeiDou.
 * output: position of the GNSS receiver, ambiguity of carrier-phase 
 * Date: 2020/11/28
 *******************************************************/

// std inputs and outputs, fstream
#include <iostream>
#include <string>  
#include <fstream>
#include<sstream>
#include <stdlib.h>
#include <iomanip>

// math
#include <math.h>
//time 
#include <time.h>
//algorithm 
#include <algorithm>

// google eigen
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include<Eigen/Core>

// google implements commandline flags processing.
#include <gflags/gflags.h>
// google loging tools
#include <glog/logging.h>
// ros
#include <ros/ros.h>
#include <rosbag/bag.h>
#include <rosbag/view.h>
/* Reference from NovAtel GNSS/INS */
#include <novatel_msgs/INSPVAX.h> // novatel_msgs/INSPVAX
#include "gnss_tools.h"
#include <nlosExclusion/GNSS_Raw_Array.h>

#include <geometry_msgs/Point32.h>
#include <stdio.h>
#include <queue>
#include <map>
#include <queue>
#include <mutex>
#include <thread>
#include <nav_msgs/Odometry.h>

#include <nav_msgs/Path.h>

#include <ceres/ceres.h>
#include <ceres/rotation.h>

#include "../tic_toc.h"

// allign 
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>

#include <sensor_msgs/NavSatFix.h>

// rtklib
#include <stdarg.h>
#include "../../RTKLIB/src/rtklib.h"

 /* factor graph related head file */
#include "../../include/gnss_estimator/psr_doppler_car_rtk.h"

// static double rand01()
// {
//     return (double)rand() / RAND_MAX;
// }

// static double randm11()
// {
//     auto ret = rand01();
//     return ret * 2 - 1;
// }

std::vector<std::string> strSplit(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

class psr_doppler_car_rtk
{
    ros::NodeHandle nh;

    /* ros subscriber */
    ros::Publisher pub_WLSENU, pub_FGOENU, pub_global_path, pub_fgo_llh;
    std::map<double, nlosExclusion::GNSS_Raw_Array> gnss_raw_map;
    // std::map<double, nav_msgs::Odometry> doppler_map;
    std::map<double, nlosExclusion::GNSS_Raw_Array> station_gnss_raw_map;

    GNSS_Tools m_GNSS_Tools; // utilities  

    /* subscriber */
    std::unique_ptr<message_filters::Subscriber<nlosExclusion::GNSS_Raw_Array>> gnss_raw_array_sub;
    std::unique_ptr<message_filters::Subscriber<nav_msgs::Odometry>> doppler_sub;
    std::unique_ptr<message_filters::Subscriber<nlosExclusion::GNSS_Raw_Array>> station_gnss_raw_array_sub;
    std::unique_ptr<message_filters::TimeSynchronizer<nlosExclusion::GNSS_Raw_Array, nlosExclusion::GNSS_Raw_Array, nav_msgs::Odometry>> syncdoppler2GNSSRaw;

    /* thread lock for data safe */
    std::mutex m_gnss_raw_mux;

    /* thread for data processing */
    // std::thread optimizationThread;

    int gnss_frame = 0;
    Eigen::Matrix<double, 3,1> ENU_ref;
    int slidingWindowSize = 10000000; // ? epoch measurements 150 100000
    bool hasNewData = false;

    /* latest state in ENU */
    Eigen::Matrix<double ,3,1> FGOENULatest;

    /* path in ENU */
    nav_msgs::Path fgo_path;

    std::string inputSolPath, outputSolPath;

private:
    // std::unique_ptr<factor_graph> factor_graph_ptr_; // factor graph ptr
    FactorGraph factor_graph;
    
    bool initializeFactorGraph(ceres::Solver::Options& options) 
    {
        /* initialize the factor graph size */
        factor_graph.setWindowSize(slidingWindowSize);

        /* set up ceres-solver options */
        factor_graph.setupSolverOptions(options);

        /* set up loss functions (Huber, Cauchy)*/
        factor_graph.setupLossFunction("Huber");
        // factor_graph.setupLossFunction("Cauchy");
    }

    void parseRtkSolution(const std::string& rtkSolPath)
    {
        std::fstream fs(rtkSolPath.c_str(), std::ios_base::in);
        std::string line;
        for (int i = 0; i < 25; ++i)
            std::getline(fs, line);

        std::vector<std::pair<double, Eigen::MatrixXd>> rtkSol;

        // srand(32);
        while (std::getline(fs, line))
        {
            std::istringstream iss(line);
            int idx;
            double gpst, lat, lon, heigh;
            if (!(iss >> idx >> gpst >> lat >> lon >> heigh)) { break; } // error

            // Eigen::MatrixXd llh;
            // llh.resize(3, 1);
            // Eigen::MatrixXd ecef;
            // ecef.resize(3, 1);
            // llh(0) = lon;
            // llh(1) = lat;
            // llh(2) = heigh;
            // ecef = m_GNSS_Tools.llh2ecef(llh);

            Eigen::MatrixXd ecef;
            ecef.resize(3, 1);

            double pos[3] = { lat * D2R, lon * D2R, heigh };
            double ecefa[3];

            pos2ecef(pos, ecefa);

            ecef(0) = ecefa[0];
            ecef(1) = ecefa[1];
            ecef(2) = ecefa[2];

            // double scale = 0.5;
            // ecef(0) += scale * randm11();
            // ecef(1) += scale * randm11();
            // ecef(2) += scale * randm11();

            rtkSol.push_back({gpst, ecef});
        }

        factor_graph.setInitialSolution(rtkSol);

        LOG(INFO) << "parse initial solution with size: " << rtkSol.size() << std::endl;
    }


public:
    psr_doppler_car_rtk(const std::string& rtkInputSolPath, const std::string& fgoOutputSolPath)
    {      
        parseRtkSolution(rtkInputSolPath);

        inputSolPath = rtkInputSolPath;
        outputSolPath = fgoOutputSolPath;

        factor_graph.setInputSolPath(inputSolPath);
        factor_graph.setOutputSolPath(outputSolPath);
        
        /* reference point for ENU calculation */
        ENU_ref<< ref_lon, ref_lat, ref_alt;

    }


    /**
	 * @brief perform factor graph optimization
	 * @param none
	 * @return none
	 */
    void solveOptimization()
    {
        if(factor_graph.getDataStreamSize()<=3 )
        {
            return;
        }

        // int ntrials = 1;
        // int seed_initial = 32;
        // double fcost_best = 1e20;
        // for (int seed = seed_initial; seed < seed_initial+ntrials; ++seed)
        {
            /* define the problem */
            ceres::Problem problem;
            ceres::Solver::Options options;
            ceres::Solver::Summary summary;

            /* start clock for factor graph optimization */
            TicToc OptTime;

            /* initialize factor graph */
            initializeFactorGraph(options);

            /* get data stream size */
            factor_graph.getDataStreamSize();

            /* setup memory for position state */
            factor_graph.setupPositionStateMemory();

            /* setup memory for ambiguity state */
            factor_graph.setConstARStateMemoryAndInitialize();

            /* initialize factor graph parameters */
            factor_graph.initializeFactorGraphParas();

            /* initialize the newly added states */
            factor_graph.initializeNewlyAddedGraph(0);

            /* add parameter blocks */
            factor_graph.addParameterBlocksToGraph(problem);

            /* add Doppler FACTORS to factor graph */
            factor_graph.addMotionFactors(problem);

            // factor_graph.addDopplerFactors(problem);

            factor_graph.addAdrFactors(problem);

            // factor_graph.addPseudorangeFactors(problem);

            // /* add double-differenced pseudorange/Carrier-phase FACTORS */
            // factor_graph.addDDPsrCarFactors(problem);

            /* solve the factor graph (float solution) */
            factor_graph.solveFactorGraphFloatSolution(problem, options, summary);

            /* save graph state to variables */
            factor_graph.saveGraphStateToVector();

            // /* solve the ambiguity resolution of current epoch */
            // factor_graph.solveAmbiguityResolutionFixedSolution();

            // /* set reference point for ENU calculation */
            // factor_graph.setupReferencePoint();

            // /* get the path in factor graph */
            // FGOENULatest = factor_graph.getLatestStateENU();

            // /* print the lastest float state in factor graph */
            // factor_graph.printLatestFloatStateENU();

            // /* print the lastest fixed state in factor graph */
            // factor_graph.printLatestFixedStateENU();

            /* publish the path from FGO */ 
            // fgo_path = factor_graph.getPathENU(fgo_path);
            // pub_global_path.publish(fgo_path);

            /* remove the data outside sliding window */
            // factor_graph.removeStatesOutsideSlidingWindow();

            /** */
            // std::cout << "OptTime-> "<< OptTime.toc()<< std::endl;
        }
    }




    /**
   * @brief gnss raw msg and doppler msg callback
   * @param gnss raw msg and doppler msg
   * @return void
   @ 
   */
   void gnssraw_doppler_msg_callback(const nlosExclusion::GNSS_Raw_ArrayConstPtr& gnss_msg, const nlosExclusion::GNSS_Raw_ArrayConstPtr& station_gnss_msg, const nav_msgs::OdometryConstPtr& doppler_msg)
    {
        // m_gnss_raw_mux.lock();
        // hasNewData = true;
        // gnss_frame++;
        // double time0 = gnss_msg->GNSS_Raws[0].GNSS_time;
        // double time1 = station_gnss_msg->GNSS_Raws[0].GNSS_time;
        // double time_frame = doppler_msg->pose.pose.position.x;

        // // std::cout<<"gnss time0 " <<time0 <<std::endl; 
        // // std::cout<<"doppler time_frame " <<time_frame <<std::endl;
        // // std::cout<<"station time1 " <<time1 <<std::endl;

        // /* save the  */
        // if(checkValidEpoch(time_frame) && m_GNSS_Tools.checkRepeating(*gnss_msg))
        // {
        //     if(gnss_msg->GNSS_Raws.size())
        //     {
        //         doppler_map[time_frame] = *doppler_msg;
        //         gnss_raw_map[time_frame] = *gnss_msg;
        //         station_gnss_raw_map[time_frame] = *station_gnss_msg;

        //         factor_graph.input_gnss_raw_data(*gnss_msg, time_frame);
        //         factor_graph.input_station_data(*station_gnss_msg, time_frame);
        //         factor_graph.input_doppler_data(*doppler_msg, time_frame);

        //         Eigen::MatrixXd eWLSSolutionECEF = m_GNSS_Tools.WeightedLeastSquare(
        //                                 m_GNSS_Tools.getAllPositions(*gnss_msg),
        //                                 m_GNSS_Tools.getAllMeasurements(*gnss_msg),
        //                                 *gnss_msg, "WLS");
        //         Eigen::Matrix<double ,3,1> WLSENU;
        //         WLSENU = m_GNSS_Tools.ecef2enu(ENU_ref, eWLSSolutionECEF);
        //         LOG(INFO) << "WLSENU -> "<< std::endl << WLSENU;

        //         nav_msgs::Odometry odometry;
        //         odometry.header.frame_id = "map";
        //         odometry.child_frame_id = "map";
        //         odometry.pose.pose.position.x = WLSENU(0);
        //         odometry.pose.pose.position.y = WLSENU(1);
        //         odometry.pose.pose.position.z = WLSENU(2);
        //         pub_WLSENU.publish(odometry);
        //     }
        // }
        
        // /* release the lock */
        // m_gnss_raw_mux.unlock();
    }

    void addEpochMessages(
        const nlosExclusion::GNSS_Raw_ArrayConstPtr& gnss_msg, 
        const nlosExclusion::GNSS_Raw_ArrayConstPtr& rover_gnss_msg, 
        const nlosExclusion::GNSS_Raw_ArrayConstPtr& station_gnss_msg, 
        const nav_msgs::OdometryConstPtr& doppler_msg,
        const nav_msgs::OdometryConstPtr& adr_msg)
    {
        // double time_frame = doppler_msg->pose.pose.position.x;
        double time_frame = gnss_msg->GNSS_Raws[0].GNSS_time;
        factor_graph.input_gnss_raw_data(*gnss_msg, time_frame);
        factor_graph.input_rover_data(*rover_gnss_msg, time_frame);
        factor_graph.input_station_data(*station_gnss_msg, time_frame);
        factor_graph.input_doppler_data(*doppler_msg, time_frame);
        factor_graph.input_adrvel_data(*adr_msg, time_frame);
    }


    std::vector<nlosExclusion::GNSS_Raw> firstRecArray;

    void addEpochMessages(
        const nlosExclusion::GNSS_Raw_ArrayConstPtr& rover_gnss_msg, 
        const nlosExclusion::GNSS_Raw_ArrayConstPtr& station_gnss_msg)
    {
        double time_frame = rover_gnss_msg->GNSS_Raws[0].GNSS_time;
        firstRecArray.push_back(rover_gnss_msg->GNSS_Raws[0]);
        factor_graph.input_rover_data(*rover_gnss_msg, time_frame);
        factor_graph.input_station_data(*station_gnss_msg, time_frame);
    }

    ~psr_doppler_car_rtk()
    {
    }

    void dumpVelocityToFile(const std::string& fname)
    {
        // char path[PATH_MAX];
        // getcwd(path,sizeof(path));

        auto items = strSplit(fname, "/");
        assert(items.size() >= 4);

        auto outName = fname.substr(0, fname.size()-items.back().size()) + "doppler_velocity.csv";

        const int GPS_TO_UTC = 315964782;

        auto tripId=items[items.size()-4]+"/"+items[items.size()-3];

        FILE* fd = fopen(outName.c_str(), "w");
        assert(fd);
        fprintf(fd, "tripId,UnixTimeMillis,vx,vy,vz\n");
        for (const auto& item: firstRecArray) {
            long unixTime = (item.GNSS_week * 7 * 24 * 3600 + item.GNSS_time + GPS_TO_UTC) * 1000;
            fprintf(fd, "%s,%ld,%f,%f,%f\n", tripId.c_str(), unixTime, item.rover_vx, item.rover_vy, item.rover_vz);
        }
        
        fclose(fd);
    }
};

int main(int argc, char **argv)
{
    FLAGS_logtostderr = 1;  // output to console
    google::InitGoogleLogging(argv[0]); // init the google logging
    google::ParseCommandLineFlags(&argc, &argv, true); // parseCommandLineFlags 
    ros::init(argc, argv, "psr_doppler_car_rtk_node"); 
    ROS_INFO("\033[1;32m----> psr_doppler_car_rtk_node Started.\033[0m"); 
    // ...
    // ros::spin();

    // const char* solPath = "/home/alexey/smartphone-decimeter-2022/data/train/2020-08-06-US-MTV-1/GooglePixel4/supplemental/gnss_log_rtklib.pos";
    // const char* solPath = "/home/alexey/smartphone-decimeter-2022/data/train/2021-03-16-US-MTV-1/GooglePixel4XL/supplemental/gnss_log_rtklib.pos";
    // const char* solPath = "/home/alexey/smartphone-decimeter-2022/data/train/2021-03-16-US-MTV-1/GooglePixel4XL/supplemental/gnss_log_rtklib_gt.pos";
    // const char* solPath = "/home/alexey/smartphone-decimeter-2022/data/train/2020-07-24-US-MTV-2/GooglePixel4XL/supplemental/gnss_log_rtklib.pos";
    
    std::string inputFile, outputFile;
    ros::param::get("inputFile", inputFile);
	ros::param::get("outputFile", outputFile);

    double stationParam[3];
	ros::param::get("station_x", stationParam[0]);
	ros::param::get("station_y", stationParam[1]);
	ros::param::get("station_z", stationParam[2]);

    bool dumpSpeed = false;
    ros::param::get("dumpSpeed", dumpSpeed);

	station_x = stationParam[0];
	station_y = stationParam[1];
	station_z = stationParam[2];

    psr_doppler_car_rtk psr_doppler_car_rtk(inputFile, outputFile);

    rosbag::Bag bag("/home/alexey/rtklib.msg");
    // rosbag::View dopplerView(bag, rosbag::TopicQuery("doppler"));
    // rosbag::View adrvelView(bag, rosbag::TopicQuery("adrvel"));
    // rosbag::View rawView(bag, rosbag::TopicQuery("raw"));
    rosbag::View roverView(bag, rosbag::TopicQuery("rover"));
    rosbag::View stationView(bag, rosbag::TopicQuery("station"));

    // auto adrIt = adrvelView.begin();
    // auto dopIt = dopplerView.begin();
    // auto rawIt = rawView.begin();
    auto roverIt = roverView.begin();
    auto stationIt = stationView.begin();

    int msgCount = 0;
    for (;  stationIt != stationView.end() && roverIt != roverView.end(); 
        ++stationIt, ++roverIt)
    {
        // const nlosExclusion::GNSS_Raw_ArrayConstPtr gnssMsg = (*rawIt).instantiate<nlosExclusion::GNSS_Raw_Array>();
        const nlosExclusion::GNSS_Raw_ArrayConstPtr stationMsg = (*stationIt).instantiate<nlosExclusion::GNSS_Raw_Array>();
        const nlosExclusion::GNSS_Raw_ArrayConstPtr roverMsg = (*roverIt).instantiate<nlosExclusion::GNSS_Raw_Array>();
        // const nav_msgs::OdometryConstPtr dopplerMsg = (*dopIt).instantiate<nav_msgs::Odometry>();
        // const nav_msgs::OdometryConstPtr adrvelMsg = (*adrIt).instantiate<nav_msgs::Odometry>();

        // psr_doppler_car_rtk.addEpochMessages(gnssMsg, roverMsg, stationMsg, dopplerMsg, adrvelMsg);
        // psr_doppler_car_rtk.addEpochMessages(gnssMsg, stationMsg, adrvelMsg);
         psr_doppler_car_rtk.addEpochMessages(roverMsg, stationMsg);

        // std::cout << (double)gnssMsg->GNSS_Raws[0].GNSS_time << std::endl;
        // printf("GNSS_time=%f\n", gnssMsg->GNSS_Raws[0].GNSS_time);

        // for(int i =0; i < gnssMsg->GNSS_Raws.size(); i++) // for weighted least square
        // {
        //     std::cout << gnssMsg->GNSS_Raws[i].prn_satellites_index
        //         << " " << gnssMsg->GNSS_Raws[i].sat_pos_x
        //         << " " << gnssMsg->GNSS_Raws[i].sat_pos_y
        //         << " " << gnssMsg->GNSS_Raws[i].sat_pos_z
        //         << std::endl;
        // }

        // break;
        // if (msgCount > 300)
        //     break;
            // psr_doppler_car_rtk.solveOptimization();

        ++msgCount;
    }

    std::cout << "message count read: " << msgCount << std::endl;

    if (dumpSpeed) {
        psr_doppler_car_rtk.dumpVelocityToFile(outputFile);

        bag.close();
        return 0;
    }

    psr_doppler_car_rtk.solveOptimization();

    bag.close();

    return 0;
}
