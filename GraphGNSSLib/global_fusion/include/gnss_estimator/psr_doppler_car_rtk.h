/*******************************************************
 * Copyright (C) 2019, Intelligent Positioning and Navigation Lab, Hong Kong Polytechnic University
 * 
 * This file is part of GraphGNSSLib.
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Weisong Wen (weisong.wen@connect.polyu.hk)
 *******************************************************/

#define D2R 3.1415926/180.0
#include <nlosExclusion/GNSS_Raw_Array.h>
// google implements commandline flags processing.
#include <gflags/gflags.h>
// google loging tools
#include <glog/logging.h>
#include "../gnss_tools.h"

#include "pseudorange_factor.h"
#include "doppler_factor.hpp"
#include "carrier_phase_factor.h"

#define enanle_Gaussian 1
#define state_size 9 // x,y,z, (clk_gps,clk_BeiDou are removed by DD)
#define velocity_size 4
#define ar_size 15  // size of variable of ambiguity 15
#define com_size state_size+ar_size // combine the state and the ambiguity in a array

#define enable_doppler_factor 1
#define doppler_vel_var 0.001 //0.01
#define enable_pseudorange_factor 0
#define enable_DD_factor 1 // enable DD measurement factor
#define enable_static_MM 0 // enable zero velocity motion model?

#define enable_ambiguity_resolution 1 // enable Ambiguity Resolution
#define use_fixed_cov 0 // use fixed covariance for DD measurement?

GNSS_Tools m_GNSS_Tools; // utilities

static double rand01()
{
    return (double)rand() / RAND_MAX;
}

static double randm11()
{
    auto ret = rand01();
    return ret * 2 - 1;
}


static std::vector<std::string> strSplit(const std::string &s, char sep) 
{
std::vector<std::string> ret; 
int from = 0; 
while (1) {
    int nxt = s.find(sep, from);
    if (nxt == std::string::npos) {
        auto to_add = s.substr(from);
        if (!to_add.empty())
            ret.push_back(to_add);
        break;
    }    
    auto to_add = s.substr(from, nxt - from);
    if (!to_add.empty())
        ret.push_back(to_add);
    from = nxt + 1; 
}    
return ret; 
}

static std::string strJoin(const std::vector<std::string> &vec, char sep) 
{
if (vec.empty()) return "";
std::string ret = vec[0];
for (int i = 1; i < vec.size(); ++i) {
    ret += sep; 
    ret += vec[i];
}    
return ret; 
}

  template <typename T>
  inline static std::string num2string(T num) 
  {
    std::stringstream ss;
    ss << std::setprecision(10) << num; 
    return ss.str();
  }

  template <typename T>
  inline static T string2num(const std::string &s)
  {
    T ret; 
    std::stringstream ss;
    ss << s;
    ss >> ret; 
    return ret; 
  }

  template <>
  inline double string2num<double>(const std::string &s)
  {
    return strtod(s.c_str(), NULL);
  }

  static std::string bloatFile(const std::string &path)
  {
    std::ifstream s(path);
    if (!s.is_open())
      return "";
    s.seekg(0, std::ios_base::end);
    int fileSize = s.tellg();
    s.seekg(0, std::ios_base::beg);
    std::string res; 
    res.reserve(fileSize);
    res.assign((std::istreambuf_iterator<char>(s)),
      std::istreambuf_iterator<char>());
    return res; 
  }


class FactorGraph{
public:
    /* continuous data stream */
    std::map<double, nlosExclusion::GNSS_Raw_Array> rover_gnss_raw_map, station_gnss_raw_map;

    int d_swpos = 0;

    /* Ceres solver object */
    // ceres::Problem problem;
    // ceres::Solver::Options options;
    // ceres::Solver::Summary summary;
    ceres::LossFunction *loss_function = nullptr;

    /* size of factor graph */
    int sizeOfFactorGraph = 10000000;

    /* position state array of factor graph */
    std::vector<double*> state_array;

    std::vector<double*> velocity_array;

    /* ambiguity state array of factor graph */
    std::vector<double*> ar_state_array;

    std::vector<double> switch_array;

    /* array save the num of ambiguity unknowns for each epoch */
    std::vector<int*> ar_state_num;

    /* gps second array */
    std::vector<int> gps_sec_array;
    
    /* position state saved in vector pair */
    std::vector<std::pair<double, std::pair<Eigen::Vector3d,Eigen::Vector3d>>> Ps;

    /* ambiguity state saved in vector pair <time, <PRN, ambiguity>> */
    std::vector<std::pair<double, std::vector<std::pair<int, double>>>> AR;

    /* size of the last factor graph optimization */
    int lastFactorGraphSize = 0;

    /* fixed variance of doppler measurements */
    double var = 0.6;

    /* reference point for ENU calculation */
    Eigen::MatrixXd ENULlhRef;

    /* measurements size */
    int measSize = 0; 

    std::vector<std::vector<Eigen::Vector3d>> adrDirectionsPerSatellite;
    std::vector<std::vector<Eigen::Vector3d>> losVectorPerSatellite;

    /* parameters */
    int numOfPsrFactors = 0;
    int numOfDopplerFactors = 0;
    int numOfStates =0;

    /* latest GNSS-RTK solution with LAMBDA */
    Eigen::Matrix<double, 3,1> fixedStateGNSSRTK;
    int fixed_cnt = 0;

    std::vector<std::pair<double, Eigen::MatrixXd>> mPosSolution;

    std::string inputSolPath, outputSolPath;
    

public:

    void setInputSolPath(const std::string& i_inputSolPath)
    {
        inputSolPath = i_inputSolPath;
    }
    void setOutputSolPath(const std::string& i_outputSolPath)
    {
        outputSolPath = i_outputSolPath;
    }

    /* input gnss raw (pseudorange/carrier-phase) data  */
    bool input_gnss_raw_data(nlosExclusion::GNSS_Raw_Array GNSS_data, double timestamp)
    {
    }

    /* input Doppler data  */
    bool input_doppler_data(nav_msgs::Odometry dopplerData, double timestamp)
    {
    }

    bool input_adrvel_data(nav_msgs::Odometry adrData, double timestamp)
    {
    }

    bool input_rover_data(nlosExclusion::GNSS_Raw_Array GNSS_data, double timestamp)
    {
        if(timestamp<0) return false;
        else 
        {
            rover_gnss_raw_map[timestamp] = GNSS_data;
            return true;
        }
    }


    /* input GNSS data from station  */
    bool input_station_data(nlosExclusion::GNSS_Raw_Array GNSS_data, double timestamp)
    {
        if(timestamp<0) return false;
        else 
        {
            station_gnss_raw_map[timestamp] = GNSS_data;
            return true;
        }
    }

    void setInitialSolution(const std::vector<std::pair<double, Eigen::MatrixXd>>& initSolution)
    {
        mPosSolution = initSolution;
    }

    /* input gnss doppler data */
    bool setWindowSize(int windowSize)
    {
        sizeOfFactorGraph = windowSize;
        return true;
    }

    /* clear data stream */
    bool clearDataStream()
    {
        station_gnss_raw_map.clear();
        return true;
    }

    /* set up ceres-solver options */
    bool setupSolverOptions(ceres::Solver::Options& options)
    {
        options.use_nonmonotonic_steps = true;
        options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
        options.trust_region_strategy_type = ceres::TrustRegionStrategyType::DOGLEG;
        options.dogleg_type = ceres::DoglegType::SUBSPACE_DOGLEG;
        options.minimizer_progress_to_stdout = true;
        options.num_threads = 8;
        options.max_num_iterations = 500; //258;
        // options.parameter_tolerance = 1e-7;
        return true;
    }

    /* set up Loss functions options */
    bool setupLossFunction(std::string loss)
    {
        if(loss=="Huber")
            loss_function = new ceres::HuberLoss(1.0);
        else 
        {
            loss_function = new ceres::CauchyLoss(1.0);
        }
        return true;
    }

    /* get data stream size */  
    int getDataStreamSize()
    {
        measSize = rover_gnss_raw_map.size();
        return measSize;
    }

    bool initializeFactorGraphParas()
    {
        numOfPsrFactors = 0;
        numOfDopplerFactors = 0;
        numOfStates =0;
    }

    /* setup state size */
    bool setupPositionStateMemory()
    {
        state_array.reserve(measSize);
        ar_state_num.reserve(measSize);
        velocity_array.reserve(measSize);
        int length = measSize;
        LOG(INFO) << "length" << length << std::endl;

        for(int i = 0; i < length;i++)
        {
            /* ECEF_x, ECEF_y, ECEF_z */
            state_array[i]  = new double[state_size]; //
            velocity_array[i] = new double[velocity_size];
            ar_state_num[i] = new int[1];
        }

        return true;
    }

    /* setup const ambiguity state size 
    */
    bool setConstARStateMemoryAndInitialize()
    {
        /* allocate memory to the ambiguity state */
        int length = measSize;
        gps_sec_array.reserve(length);
        ar_state_array.reserve(length);

        for(int i = 0; i < length;i++)
        {
            /* ECEF_x, ECEF_y, ECEF_z */
            ar_state_array[i] = new double[ar_size]; //
            for(int j = 0; j < ar_size; j++)
            {
                ar_state_array[i][j] = 0;
            }
        }

        {
            int total_eqs = 0;
            int length = measSize;
            auto iter_sm = station_gnss_raw_map.begin();
            for(int k = 0;  k < length; k++,iter_sm++) // tranverse the whole station gnss measurements map
            {
                nlosExclusion::GNSS_Raw_Array st_gnss_data = (iter_sm->second);
                int sv_cnt = st_gnss_data.GNSS_Raws.size();
                double t = iter_sm->first;
                
                /* find user end gnss data with closest time */
                nlosExclusion::GNSS_Raw_Array closest_gnss_data;
                findClosestEpoch(t, rover_gnss_raw_map, closest_gnss_data);

                if (closest_gnss_data.GNSS_Raws.empty())
                    continue;

                total_eqs += sv_cnt;
            }

            total_eqs += rover_gnss_raw_map.size();

            switch_array.resize(3*total_eqs);
        }

        return true;
    }

    Eigen::MatrixXd getInterpPosition(double gpst)
    {
        assert(!mPosSolution.empty());
        if (gpst < mPosSolution.front().first)
            return mPosSolution.front().second;
        if (gpst > mPosSolution.back().first)
            return mPosSolution.back().second;
        for (int i = 1; i < mPosSolution.size(); ++i)
        {
            auto gpst0 = mPosSolution[i-1].first;
            auto gpst1 = mPosSolution[i].first;
            if (gpst >= gpst0 && gpst <= gpst1)
            {
                double dt = gpst1 - gpst0;
                double dt0 = gpst - gpst0;
                
                Eigen::MatrixXd ecef;
                ecef.resize(3, 1);

                ecef = ((dt - dt0) / dt)  * mPosSolution[i-1].second + (dt0 / dt) * mPosSolution[i].second;
                return ecef;
            }
        }
        assert(0);
    }

    /* initialize the newly added state using WLS*/
    bool initializeNewlyAddedGraph(int seed)
    {
        

        int length = measSize;

        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter;
        iter = rover_gnss_raw_map.begin();

        srand(seed);
        for(int i = 0; i <length; i++,iter++)
        {
            // if(i>=(lastFactorGraphSize-1))
            {
                nlosExclusion::GNSS_Raw_Array gnss_data = (iter->second);
                // Eigen::MatrixXd eWLSSolutionECEF = m_GNSS_Tools.WeightedLeastSquare(
                //                     m_GNSS_Tools.getAllPositions(gnss_data),
                //                     m_GNSS_Tools.getAllMeasurements(gnss_data),
                //                     gnss_data, "WLS");
                state_array[i][0] = 0; //eWLSSolutionECEF(0);
                state_array[i][1] = 0; //eWLSSolutionECEF(1); 
                state_array[i][2] = 0; //eWLSSolutionECEF(2);

                state_array[i][3] = 0;
                state_array[i][4] = 0;
                state_array[i][5] = 0;
                state_array[i][6] = 0;
                state_array[i][7] = 0;
                state_array[i][8] = 0;

                // state_array[i][0] = eWLSSolutionECEF(0);
                // state_array[i][1] = eWLSSolutionECEF(1); 
                // state_array[i][2] = eWLSSolutionECEF(2);

                double scale = 0.25;
                auto posEcef = getInterpPosition(gnss_data.GNSS_Raws[0].GNSS_time);

                state_array[i][0] = posEcef(0) + scale * randm11();
                state_array[i][1] = posEcef(1) + scale * randm11();
                state_array[i][2] = posEcef(2) + scale * randm11();

                for (int j = 0; j < velocity_size; ++j)
                    velocity_array[i][j] = 0.0;

                velocity_array[i][0] = gnss_data.GNSS_Raws[0].rover_vx + scale * randm11();
                velocity_array[i][1] = gnss_data.GNSS_Raws[0].rover_vy + scale * randm11();
                velocity_array[i][2] = gnss_data.GNSS_Raws[0].rover_vz + scale * randm11();

                // std::cout << gnss_data.GNSS_Raws[0].GNSS_time << " " << posEcef(0) << " " << posEcef(1) << " " << posEcef(2) << std::endl;
            }

            std::fill(switch_array.begin(), switch_array.end(), 30.0); //10.0);
        }

        losVectorPerSatellite.resize(length);

        iter = rover_gnss_raw_map.begin();
        for(int i = 0; i <length; i++, iter++)
        {
            losVectorPerSatellite[i].resize(MAXSAT);

            double pos[3],E[9],a[3],e[3];
            ecef2pos(state_array[i],pos); xyz2enu(pos,E);

            nlosExclusion::GNSS_Raw_Array gnss_data = (iter->second);
            for (const auto& item: gnss_data.GNSS_Raws)
            {
                double cosel=cos(item.elevation * D2R);
                a[0]=sin(item.azimuth * D2R)*cosel;
                a[1]=cos(item.azimuth * D2R)*cosel;
                a[2]=sin(item.elevation * D2R);
                matmul("TN",3,1,3,1.0,E,a,0.0,e);

                losVectorPerSatellite[i][item.prn_satellites_index] = Eigen::Vector3d(e[0], e[1], e[2]).normalized();
            }
        }
        

        // std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator prev_iter;
        // prev_iter = rover_gnss_raw_map.begin();
        // iter = prev_iter;
        // iter++;

        // for(int i = 1; i <length; i++,prev_iter=iter,iter++) {
        //     double dt = iter->first - prev_iter->first;
        //     for (auto& item: prev_iter->second.GNSS_Raws) {
        //         item.rover_vx = (state_array[i][0] - state_array[i-1][0]) / dt;
        //         item.rover_vy = (state_array[i][1] - state_array[i-1][1]) / dt;
        //         item.rover_vz = (state_array[i][2] - state_array[i-1][2]) / dt;
        //     }
        // }

        return true;
    }

    /* add parameter blocks */
    bool addParameterBlocksToGraph(ceres::Problem& problem)
    {
        int length = measSize;
        for(int i = 0; i <length; i++)
        {
            /* add parameter block for position state (ECEF_x, ECEF_y, ECEF_z) */
            problem.AddParameterBlock(state_array[i], state_size);

            /* add parameter block for ambiguity state */
            problem.AddParameterBlock(ar_state_array[i], ar_size);

            problem.AddParameterBlock(velocity_array[i], velocity_size);
        }

        for (int i = 0; i < switch_array.size(); ++i)
            problem.AddParameterBlock(&switch_array[i], 1);

        return true;
    }

    // /* fix the first parameter block */
    // bool fixFirstState(bool flag, ceres::Problem& problem)
    // {
    //     /* fixed the first state only after first optimization */
    //     if(flag && Ps.size()>10)
    //     {
    //         problem.SetParameterBlockConstant(state_array[0]);
    //         problem.SetParameterBlockConstant(ar_state_array[0]);
    //     }    
    //     return true;
    // }

/*
    bool addDopplerFactors(ceres::Problem& problem)
    {
        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter_pr;
        iter_pr = rover_gnss_raw_map.begin();
        int length = measSize;
        int fidx = 0;
        for(int m = 0;  m < length; m++,iter_pr++) // 
        {
            nlosExclusion::GNSS_Raw_Array gnss_data = (iter_pr->second);
            int sv_cnt = gnss_data.GNSS_Raws.size();
            double t = gnss_data.GNSS_Raws[0].GNSS_time;

            for(int i =0; i < sv_cnt; i++)
            {
                double s_g_x = 0, s_g_y = 0,s_g_z = 0, var = 1;
                double s_g_vx = 0, s_g_vy = 0,s_g_vz = 0;

                if (gnss_data.GNSS_Raws[i].freqidx != fidx)
                    continue;

                auto sat_id = gnss_data.GNSS_Raws[i].prn_satellites_index;
                int sat_sys = satsys(sat_id, NULL);
                if (sat_sys == SYS_NONE)
                    continue;

                s_g_x = gnss_data.GNSS_Raws[i].sat_pos_x;
                s_g_y = gnss_data.GNSS_Raws[i].sat_pos_y;
                s_g_z = gnss_data.GNSS_Raws[i].sat_pos_z;

                s_g_vx = gnss_data.GNSS_Raws[i].sat_pos_vx;
                s_g_vy = gnss_data.GNSS_Raws[i].sat_pos_vy;
                s_g_vz = gnss_data.GNSS_Raws[i].sat_pos_vz;

                auto losVector = losVectorPerSatellite[m][sat_id];

                if (gnss_data.GNSS_Raws[i].doppler == 0)
                    continue;
                if (losVector == Eigen::Vector3d(0,0,0))
                    continue;
                // if (gnss_data.GNSS_Raws[i].visable & 1)
                //     continue;

                auto* velFactor = new velocityFactor(gnss_data.GNSS_Raws[i].sat_clk_drift_err,
                                                losVector,
                                                gnss_data.GNSS_Raws[i].doppler,
                                                s_g_x, s_g_y, s_g_z, 
                                                s_g_vx, s_g_vy, s_g_vz);
                
                double residuals[1];
                velFactor->operator()(state_array[m], velocity_array[m], residuals);

                if (abs(residuals[0]) > 250)
                {
                    continue;
                }

                // std::cout << "vel-res," << std::fixed << std::setprecision(10) << residuals[0] << ",0,0" << std::endl;


                ceres::CostFunction* ps_function = new ceres::AutoDiffCostFunction<velocityFactor, 1 
                                                                , state_size, velocity_size>(velFactor);

                auto ID = problem.AddResidualBlock(ps_function, loss_function, state_array[m], velocity_array[m]);
            }            
        }

        return true;
    }
    */

    bool addMotionFactors(ceres::Problem& problem)
    {
        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter_pr, iter_pr_next;
        iter_pr = rover_gnss_raw_map.begin();
        int length = measSize;
        for(int m = 0;  m < length; m++,iter_pr++)
        {
            nlosExclusion::GNSS_Raw_Array gnss_prev;
            findClosestEpoch(iter_pr->first, rover_gnss_raw_map, gnss_prev);

            iter_pr_next = iter_pr;
            iter_pr_next++;

            if(iter_pr_next == rover_gnss_raw_map.end())
                continue;

            double delta_t = iter_pr_next->first - iter_pr->first;
            double v_x_i = iter_pr->second.GNSS_Raws[0].rover_vx;
            double v_y_i = iter_pr->second.GNSS_Raws[0].rover_vy;
            double v_z_i = iter_pr->second.GNSS_Raws[0].rover_vz;

            // double v_x_i = gnss_prev.GNSS_Raws[0].rover_vx;
            // double v_y_i = gnss_prev.GNSS_Raws[0].rover_vy;
            // double v_z_i = gnss_prev.GNSS_Raws[0].rover_vz;

            // double v_norm = sqrt(v_x_i*v_x_i + v_y_i*v_y_i + v_z_i*v_z_i);
            // if (v_norm < 10) {
            //     // v_x_i = v_y_i = v_z_i = 0;
            //     continue;
            // }


            double var = doppler_vel_var; // 
            // double var = 0.6;
            Eigen::Vector3d var_vec(var,var,var);

            auto* dopFactor = new dopplerFactor(v_x_i, v_y_i, v_z_i, delta_t,  var_vec);

            double residuals[3];
            dopFactor->operator()(state_array[m], state_array[m+1], velocity_array[m], &switch_array[d_swpos], residuals);

            // std::cout << "dop-res," << std::fixed << std::setprecision(10) << residuals[0] << "," << residuals[1] << "," << residuals[2] << std::endl;

            // if (abs(residuals[0]) > 200 || abs(residuals[1]) > 300 || abs(residuals[2]) > 300)
            // {
            //     continue;
            // }

            ceres::CostFunction* doppler_function = new ceres::AutoDiffCostFunction<dopplerFactor, 3 
                                                    , state_size,state_size,velocity_size,1>(dopFactor);
            problem.AddResidualBlock(doppler_function, loss_function, state_array[m],
                state_array[m+1], velocity_array[m], &switch_array[d_swpos++]);

            // auto* accFactor = new accelFactor(delta_t);
            // ceres::CostFunction* accel_function = new ceres::AutoDiffCostFunction<accelFactor, 2
            //                                         , velocity_size,velocity_size>(accFactor);
            // problem.AddResidualBlock(accel_function, loss_function, velocity_array[m], velocity_array[m+1]);
        }

        return true;
    }

    bool addAdrFactors(ceres::Problem& problem)
    {
        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter_pr, iter_pr_next;
        iter_pr = rover_gnss_raw_map.begin();
        int length = measSize;
        int skipped = 0, added = 0;
        for(int m = 0;  m < length; m++,iter_pr++)
        {
            nlosExclusion::GNSS_Raw_Array gnss_prev = iter_pr->second;

            iter_pr_next = iter_pr;
            iter_pr_next++;

            if(iter_pr_next == rover_gnss_raw_map.end())
                continue;

            nlosExclusion::GNSS_Raw_Array gnss_curr = iter_pr_next->second;

            int fidx = 0;
            for (auto& prev_item: gnss_prev.GNSS_Raws) {
                if (prev_item.freqidx != fidx) continue;

                double delta_t = iter_pr_next->first - iter_pr->first;

                double sat_id = prev_item.prn_satellites_index;

                nlosExclusion::GNSS_Raw curr_item;

                int nfound_mask = 0;
                for (auto it = gnss_curr.GNSS_Raws.begin(); it != gnss_curr.GNSS_Raws.end(); ++it) {
                    if (it->freqidx != fidx) continue;

                    if (it->prn_satellites_index == prev_item.prn_satellites_index) {
                        curr_item = *it;
                        nfound_mask |= 1;
                        break;
                    }
                }

                if (nfound_mask != 1)
                    continue;

                auto losVectorPrev = losVectorPerSatellite[m][sat_id];
                auto losVectorCurr = losVectorPerSatellite[m+1][sat_id];

                if ((prev_item.visable & 1) || (curr_item.visable & 1) || 
                    losVectorPrev == Eigen::Vector3d(0,0,0) || losVectorCurr == Eigen::Vector3d(0,0,0))
                {
                    double s_g_x = 0, s_g_y = 0,s_g_z = 0, var = 1;
                    double s_g_vx = 0, s_g_vy = 0,s_g_vz = 0;

                    int sat_sys = satsys(sat_id, NULL);
                    if (sat_sys == SYS_NONE) {
                        skipped++;
                        continue;
                    }

                    s_g_x = prev_item.sat_pos_x;
                    s_g_y = prev_item.sat_pos_y;
                    s_g_z = prev_item.sat_pos_z;

                    s_g_vx = prev_item.sat_pos_vx;
                    s_g_vy = prev_item.sat_pos_vy;
                    s_g_vz = prev_item.sat_pos_vz;

                    if (prev_item.doppler == 0 || losVectorPrev == Eigen::Vector3d(0,0,0)) {
                        skipped++;
                        continue;
                    }

                    auto* velFactor = new velocityFactor(prev_item.sat_clk_drift_err,
                                                    losVectorPrev,
                                                    prev_item.doppler,
                                                    s_g_x, s_g_y, s_g_z, 
                                                    s_g_vx, s_g_vy, s_g_vz);

                    double residuals[1];
                    velFactor->operator()(state_array[m], velocity_array[m], &switch_array[d_swpos], residuals);

                    if (abs(residuals[0]) > 250)
                    {
                        skipped++;
                        continue;
                    }

                    ceres::CostFunction* ps_function = new ceres::AutoDiffCostFunction<velocityFactor, 1 
                                                , state_size, velocity_size, 1>(velFactor);

                    problem.AddResidualBlock(ps_function, loss_function, state_array[m], velocity_array[m], &switch_array[d_swpos++]);

                    continue;
                }

                auto* func = new CarrierPhaseDeltaAdrFactor(delta_t, prev_item, curr_item, losVectorPrev, losVectorCurr,
                    curr_item.lamda, curr_item.phase_var);

                double residuals[1];
                func->operator()(state_array[m], state_array[m+1], velocity_array[m], &switch_array[d_swpos], residuals);

                // std::cout << "adr-res," << std::fixed << std::setprecision(10) << residuals[0] << ",0,0" << std::endl;
                // std::cout << "adr-state-res," << (state_array[m+1][0] - state_array[m][0]) << "," << (state_array[m+1][1] - state_array[m][1])
                //         << "," << (state_array[m+1][2] - state_array[m][2]) << std::endl;

                if (abs(residuals[0]) > 750) // 1e6 * 5)
                {
                    skipped++;
                    continue;
                }

                ceres::CostFunction* phase_function = new ceres::AutoDiffCostFunction<CarrierPhaseDeltaAdrFactor, 1, 
                    state_size, state_size, velocity_size, 1>(func);
                problem.AddResidualBlock(phase_function, loss_function, 
                    state_array[m], state_array[m+1], velocity_array[m], &switch_array[d_swpos++]);
                added++;
            }
        }

        LOG(INFO) << "adr skipped: " << skipped << " total: " << (skipped + added);

        return true;
    }


    /* add pseudorange FACTORS */
    bool addPseudorangeFactors(ceres::Problem& problem)
    {
        /* add pseudorange factor */
        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter_pr;
        iter_pr = rover_gnss_raw_map.begin();
        int length = measSize;
        for(int m = 0;  m < length; m++,iter_pr++) // 
        {
            nlosExclusion::GNSS_Raw_Array gnss_data = (iter_pr->second);
            // gps_sec_array[m] = int(gnss_data.GNSS_Raws[0].GNSS_time);
            // state_gps_sec_vec.push_back(int(gnss_data.GNSS_Raws[0].GNSS_time));
            // MatrixXd weight_matrix; //goGPS weighting
            // weight_matrix = m_GNSS_Tools.cofactorMatrixCal_WLS(gnss_data, "WLS"); //goGPS
            // std::cout << "weight_matrix-> "<<weight_matrix<<std::endl;
            int sv_cnt = gnss_data.GNSS_Raws.size();
            double t = gnss_data.GNSS_Raws[0].GNSS_time;

            int factor_index = -1;
            for(int i =0; i < sv_cnt; i++)
            {
                factor_index++;

                double s_g_x = 0, s_g_y = 0,s_g_z = 0, var = 1;
                double pseudorange = 0;
                // if(m_GNSS_Tools.PRNisGPS(gnss_data.GNSS_Raws[i].prn_satellites_index)) sat_sys = "GPS";
                // else sat_sys = "BeiDou";

                int sat_sys = satsys(gnss_data.GNSS_Raws[i].prn_satellites_index, NULL);
                if (sat_sys == SYS_NONE)
                    continue;
                int sat_sys_idx = -1;
                // sat_sys_idx = ffs(sat_sys >> 1) - 1;
                for (int i = 0; i < 8; ++i)
                    if (sat_sys & (1 << i)) {
                        sat_sys_idx = i;
                        break;
                    }
                
                if (sat_sys_idx == -1)
                    continue;

                // printf("%d,%d\n", sat_sys, sat_sys_idx);

                s_g_x = gnss_data.GNSS_Raws[i].sat_pos_x;
                s_g_y = gnss_data.GNSS_Raws[i].sat_pos_y;
                s_g_z = gnss_data.GNSS_Raws[i].sat_pos_z;

                // double pos_var = 2;
                double pos_var = 2*1e-1; // gives 9.552
                // auto pos_var = sqrt(gnss_data.GNSS_Raws[i].pos_var);

                pseudorange = gnss_data.GNSS_Raws[i].pseudorange;

                pseudorange -= gnss_data.GNSS_Raws[i].err_tropo;
                pseudorange -= gnss_data.GNSS_Raws[i].err_iono;

                double ele = gnss_data.GNSS_Raws[i].elevation;
                double snr = gnss_data.GNSS_Raws[i].snr;

                pos_var = pos_var / sin(ele * D2R);

                // pos_var += 0.01 * pow(10,0.1*MAX(52-snr*4*SNR_UNIT,0));

                // printf("pos_var=%f ele=%f\n", pos_var, ele);

                ceres::CostFunction* ps_function = new ceres::AutoDiffCostFunction<pseudorangeFactor, 1 
                                                                , state_size>(new 
                                                                pseudorangeFactor(sat_sys_idx, s_g_x, s_g_y, s_g_z, 
                                                                    pseudorange, pos_var));

                // assert(sw_pos < switch_array.size());

                auto ID = problem.AddResidualBlock(ps_function, loss_function, state_array[m]); //, &switch_array[sw_pos]);
                numOfPsrFactors++;
            }            
        }
       
        return true;
    }


    /* add double-differenced pseudorange/Carrier-phase FACTORS */
    bool addDDPsrCarFactors(ceres::Problem& problem)
    {
        // for (int f = 0; f < NFREQ; ++f)
        for (int f = 0; f < 1; ++f)
        {
            std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter_sm; // station gnss measurements map iterator
            int length = measSize;
            iter_sm = station_gnss_raw_map.begin();
            for(int k = 0;  k < length; k++,iter_sm++) // tranverse the whole station gnss measurements map
            {
                nlosExclusion::GNSS_Raw_Array st_gnss_data = (iter_sm->second);
                int sv_cnt = st_gnss_data.GNSS_Raws.size();
                double t = iter_sm->first;
                
                /* find user end gnss data with closest time */
                nlosExclusion::GNSS_Raw_Array closest_gnss_data;
                findClosestEpoch(t, rover_gnss_raw_map, closest_gnss_data);

                if (closest_gnss_data.GNSS_Raws.empty())
                    continue;

                gps_sec_array[k] = int(closest_gnss_data.GNSS_Raws[0].GNSS_time);

                /* get the dd measurements between
                * 1. st_gnss_data from station
                * 2. closest_gnss_data from user end
                */

                /* tranverse station gnss data at a epoch */
                int carrier_phase_index = 0;
                for(int q = 0; q < sv_cnt; q++)
                {
                    if (st_gnss_data.GNSS_Raws[q].freqidx != f) continue;
                    
                    double sat_id = st_gnss_data.GNSS_Raws[q].prn_satellites_index;
                    /*u_master_sv: user to master 
                    *u_iSV: user to ith satellite
                    *r_master_sv: reference to master satellite
                    */
                    nlosExclusion::GNSS_Raw u_master_sv, u_iSV, r_master_sv;
                    
                    /* find the master satellite from the user end */
                    if(findMasterSatellite(sat_id, closest_gnss_data, u_master_sv, u_iSV, f))
                    {
                        if(u_iSV.elevation==0)
                        {
                            LOG(INFO) << "satellite with zero elevation angle---";
                            LOG(INFO) << "satellite with pseudorange---"<<u_iSV.pseudorange;
                        }
                        /* find the satellite from station gnss iwth the same id with master satellite */
                        findSatellitewithSameId(u_master_sv.prn_satellites_index, st_gnss_data, r_master_sv, f);

                        DDMeasurement DD_measurement;
                        DD_measurement.u_master_SV = u_master_sv;
                        DD_measurement.u_iSV = u_iSV;

                        DD_measurement.r_master_SV = r_master_sv;
                        DD_measurement.r_iSV = st_gnss_data.GNSS_Raws[q];

                        if (u_master_sv.pos_var == 0 || u_master_sv.phase_var == 0)
                            continue;

                        double dpos_var = std::max(1e-4, std::min(u_master_sv.pos_var, 1e4));
                        double dphase_var = std::max(1e-4, std::min(u_master_sv.phase_var, 1e4));

                        Eigen::Vector3d base_pose(station_x, station_y, station_z);
                        if (1)
                        {
                            auto* ddFactor = new DDpseudorangeVSVFactor(DD_measurement, base_pose, dpos_var);
                            
                            double residuals[1];
                            ddFactor->operator()(state_array[k], &switch_array[d_swpos], residuals);

                            if (abs(residuals[0]) > 500)
                            {
                                continue;
                            }

                            // std::cout << "dd-res," << std::fixed << std::setprecision(10) << residuals[0] << ",0,0" << std::endl;

                            ceres::CostFunction* dd_pr_function = new ceres::AutoDiffCostFunction<DDpseudorangeVSVFactor, 1 
                                                                , state_size, 1>(ddFactor);

                            assert(d_swpos < switch_array.size());

                            auto ID = problem.AddResidualBlock(dd_pr_function, loss_function, state_array[k], &switch_array[d_swpos]);

                            d_swpos++;
                        }

                        // if(checkCarrierPhaseConsistency(DD_measurement)) 
                        // if (checkCarrierPhaseConsistencyCycleSlip(DD_measurement))
                        // (carrier_phase_index<10))
                        // if (0)
                        if (0)
                        { 
                            ceres::CostFunction* dd_cp_function = new ceres::AutoDiffCostFunction<DDCarrierPhaseFactor_DDBias, 1 
                                                            , state_size,ar_size>(new 
                                                            DDCarrierPhaseFactor_DDBias(DD_measurement, base_pose,carrier_phase_index, 
                                                            dphase_var));
                            auto IDs = problem.AddResidualBlock(dd_cp_function, loss_function, state_array[k],ar_state_array[k]); 

                            carrier_phase_index++;
                            // LOG(INFO)<<"add double-difference carrier phase factor->"<<carrier_phase_index;
                        }
                        else
                        {
                            // LOG(INFO)<<"no carrier-phase measurement";
                        }
                    }
                }

                // LOG(INFO) << "sv_cnt/carrier_phase_index-> "<<sv_cnt<<"  "<<carrier_phase_index;
                ar_state_num[k][0] = carrier_phase_index;
                
            }
        }
    }


    /* solve the factor graph */
    double solveFactorGraphFloatSolution(ceres::Problem& problem,
        ceres::Solver::Options& options, ceres::Solver::Summary& summary)
    {
        /* solve the problem*/
        ceres::Solve(options, &problem, &summary);

        std::cout << summary.FullReport() << std::endl;
        // std::cout << summary.BriefReport() << std::endl;

        return summary.final_cost;
    }

    /* save graph state to vector for next solving */
    bool saveGraphStateToVector()
    {
        /* save the size of current factor graph */
        // lastFactorGraphSize = measSize;

        /* get time from data stream */
        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter_pr;
        iter_pr = rover_gnss_raw_map.begin();
        int length = measSize;

        Ps.clear();

        /* tranverse the stateArray */
        for(int m = 0;  m < length; m++,iter_pr++) // 
        {
            nlosExclusion::GNSS_Raw_Array gnss_data = (iter_pr->second);
            double time = gnss_data.GNSS_Raws[0].GNSS_time;
            double prn = gnss_data.GNSS_Raws[0].prn_satellites_index;

            auto val = std::make_pair(time, 
                    std::make_pair(Eigen::Vector3d(state_array[m][0],state_array[m][1],state_array[m][2]),
                    Eigen::Vector3d(velocity_array[m][0],velocity_array[m][1],velocity_array[m][2]))
                    );

            /* if the state vector is empty, override */
            if(Ps.size()==0)
            {
                Ps.push_back(val);
                // Clocks.push_back(std::make_pair(time, Eigen::Vector2d(state_array[m][3],state_array[m][4])));
            }
            /* if the state vector is NOT empty, update */
            else
            {
                bool findTimeKey = false;
                for(int i = 0; i < Ps.size(); i++)
                {
                    if(time == Ps[i].first)
                    {
                        Ps[i] = val;
                        // Clocks[i] = std::make_pair(time, Eigen::Vector2d(state_array[m][3],state_array[m][4]));
                        findTimeKey = true;
                    }
                }
                /* new time frame, add to state vector*/
                if(findTimeKey==false)
                {
                    Ps.push_back(val);
                    // Clocks.push_back(std::make_pair(time, Eigen::Vector2d(state_array[m][3],state_array[m][4])));
                } 
            }
        }

        {
            // std::ifstream infile(inputSolPath.c_str());
            auto contents = bloatFile(inputSolPath);
            auto lines = strSplit(contents, '\n');
            FILE* wfd = fopen(outputSolPath.c_str(), "wt");

            int s_ix = 0;

            // std::string line;
            for (int i = 0; i < 25; ++i)
            {
                auto line = lines[s_ix++];
                if (line[line.size() - 1] == '\r')
                    line.resize(line.size() - 1);
                line += '\n';
                fwrite(line.c_str(), 1, line.size(), wfd);
            }

            //2115 512144.435   37.629524771 -122.422781146     9.7881   2  15   0.5565   0.6210   0.9971   0.2210  -0.5925   0.4939 -15.57    0.0


            int ix = 0;
            while (s_ix < lines.size()) //std::getline(infile, line))
            {
                auto line = lines[s_ix++];

                if (line[line.size() - 1] == '\r')
                    line.resize(line.size() - 1);

                if (ix >= Ps.size())
                {
                    line += '\n';
                    fwrite(line.c_str(), 1, line.size(), wfd);
                    continue;
                }

                auto fields = strSplit(line, ' ');

                double curr = string2num<double>(fields[1]);

                double tol = 0.03;
                while (ix < Ps.size() && Ps[ix].first < curr - tol) { ix++; } 

                if (ix >= Ps.size())
                {
                    line += '\n';
                    fwrite(line.c_str(), 1, line.size(), wfd);
                    continue;
                }

                if (fabs(Ps[ix].first - curr) < tol)
                {
                    // Eigen::MatrixXd ecef;
                    // ecef.resize(3, 1);
                    // Eigen::MatrixXd llh;
                    // llh.resize(3, 1);
                    // ecef(0) = Ps[ix].second(0);
                    // ecef(1) = Ps[ix].second(1);
                    // ecef(2) = Ps[ix].second(2);
                    // llh = m_GNSS_Tools.ecef2llh(ecef);

                    Eigen::MatrixXd llh;
                    llh.resize(3, 1);

                    double pos[3];
                    double ecefa[3] = { Ps[ix].second.first(0), Ps[ix].second.first(1), Ps[ix].second.first(2) };
                    double vel_ecef[3] = { Ps[ix].second.second(0), Ps[ix].second.second(1), Ps[ix].second.second(2) };

                    ecef2pos(ecefa, pos);

            //     fprintf(wfd, "%f,%f,%f,%f\n", item.first, llh(1), llh(0), llh(2));

                    fields[3] = num2string(pos[1] * R2D); //llh(0));
                    fields[2] = num2string(pos[0] * R2D); //llh(1));
                    fields[4] = num2string(pos[2]); //llh(2));

                    fields[10] = num2string(vel_ecef[0]);
                    fields[11] = num2string(vel_ecef[1]);
                    fields[12] = num2string(vel_ecef[2]);

                    auto to_write = strJoin(fields, ' ') + '\n';
                    fwrite(to_write.c_str(), 1, to_write.size(), wfd);

                    ix++;
                    continue;
                }

                if (Ps[ix].first > curr)
                {
                    line += '\n';
                    fwrite(line.c_str(), 1, line.size(), wfd);
                }
            }

            fclose(wfd);

            // fprintf(wfd, "gpst,lat,lon,heigh\n");
            // for (const auto& item: Ps)
            // {
            //     Eigen::MatrixXd ecef;
            //     ecef.resize(3, 1);
            //     Eigen::MatrixXd llh;
            //     llh.resize(3, 1);
            //     ecef(0) = item.second(0);
            //     ecef(1) = item.second(1);
            //     ecef(2) = item.second(2);
            //     llh = m_GNSS_Tools.ecef2llh(ecef);

            //     fprintf(wfd, "%f,%f,%f,%f\n", item.first, llh(1), llh(0), llh(2));
            // }

        }

        return true;
    }
    
    /* set up the reference point for ENU calculation */
    bool setupReferencePoint()
    {
        /* reference point for ENU calculation */
        ENULlhRef.resize(3,1);
        ENULlhRef<< ref_lon, ref_lat, ref_alt;
        return true;
    }

    /* get the latest state in ENU */
    Eigen::Matrix<double ,3,1> getLatestStateENU()
    {
        int length = measSize;
        Eigen::Matrix<double ,3,1> FGOENU;
        Eigen::Matrix<double, 3,1> state;
        state<< state_array[length-1][0], 
                state_array[length-1][1], 
                state_array[length-1][2];
        FGOENU = m_GNSS_Tools.ecef2enu(ENULlhRef, state);    
        return FGOENU;
    }

     /* print the latest state in ENU */
    bool printLatestFloatStateENU()
    {
        int length = measSize;
        Eigen::Matrix<double ,3,1> FGOENU;
        Eigen::Matrix<double, 3,1> state;
        state<< state_array[length-1][0], 
                state_array[length-1][1], 
                state_array[length-1][2];
        FGOENU = m_GNSS_Tools.ecef2enu(ENULlhRef, state);    
        std::cout << "FGOENU-> "<< FGOENU<< std::endl;  
        return true;
    }

    /* print the latest state in ENU */
    bool printLatestFixedStateENU()
    {
        int length = measSize;
        Eigen::Matrix<double ,3,1> fixedFGOENU;
        fixedFGOENU = m_GNSS_Tools.ecef2enu(ENULlhRef, fixedStateGNSSRTK);    
        std::cout << "fixedFGOENU-> "<< fixedFGOENU<< std::endl;  
        return true;
    }

    /* get the path of FGO in ENU */
    nav_msgs::Path getPathENU(nav_msgs::Path& fgo_path)
    {
        int length = measSize;
        Eigen::Matrix<double ,3,1> FGOENU;
        Eigen::Matrix<double, 3,1> state;
        fgo_path.poses.clear();
        fgo_path.header.frame_id = "map";
        for(int i = 0; i < length;i++)
        {
            state<< state_array[i][0], 
                    state_array[i][1], 
                    state_array[i][2];
            FGOENU = m_GNSS_Tools.ecef2enu(ENULlhRef, state);  
            geometry_msgs::PoseStamped pose_stamped;
            pose_stamped.header.stamp = ros::Time::now();
            pose_stamped.header.frame_id = "map";
            pose_stamped.pose.position.x = FGOENU(0);
            pose_stamped.pose.position.y = FGOENU(1);
            pose_stamped.pose.position.z = 10;
            fgo_path.poses.push_back(pose_stamped);
            // std::cout << "pose_stamped- FGO-> "<< std::endl<< pose_stamped;
        }
              
        return fgo_path;
    }

   /**
   * @brief maintain sliding window slidingWindowSize
   * @param gnss raw msg and doppler msg
   * @return void
   @ 
   */

    // void removeStatesOutsideSlidingWindow()
    // {
    //     int numElementsToRemove = 0;
        
    //     /* sliding window gnss raw pseudorange*/
    //     numElementsToRemove = gnss_raw_map.size() - sizeOfFactorGraph;
    //     if(numElementsToRemove<=0) return;
    //     auto i = gnss_raw_map.begin();
    //     while (i != gnss_raw_map.end() && numElementsToRemove > 0)
    //     {
    //         i = gnss_raw_map.erase(i);
    //         --numElementsToRemove;
    //     }

    //     /* sliding window station gnss raw pseudorange*/
    //     numElementsToRemove = station_gnss_raw_map.size() - sizeOfFactorGraph;
    //     if(numElementsToRemove<=0) return;
    //     auto s = station_gnss_raw_map.begin();
    //     while (s != station_gnss_raw_map.end() && numElementsToRemove > 0)
    //     {
    //         s = station_gnss_raw_map.erase(s);
    //         --numElementsToRemove;
    //     }
    // }


    /* free memory */
    bool freeStateMemory()
    {
        int length = measSize;
        for(int i = 0; i < length;i++)
        {
            free(state_array[i]);
        }     
        return true;
    }

};


/* check the valid epoch based on gps time span*/
bool checkValidEpoch(double gps_sec)
{
    return true;
    // if((gps_sec >= start_gps_sec) && (gps_sec <=end_gps_sec))
    // {
    //     return true;
    // }
    // else return false;
}