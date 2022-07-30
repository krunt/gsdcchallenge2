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

/* doppler factor*/
struct dopplerFactor
{
    dopplerFactor(
        double v_x, double v_y, double v_z,
        double delta_t, Eigen::Vector3d var_vector)
                : v_x(v_x), v_y(v_y), v_z(v_z),
                delta_t(delta_t), var_vector(var_vector){}

    template <typename T>
    bool operator()(const T* state_i, const T* state_j, const T* vel_state, const T* switchvar, T* residuals) const
    {
        T est_v_x = (state_j[0] - state_i[0])/ delta_t;
        T est_v_y = (state_j[1] - state_i[1])/ delta_t;
        T est_v_z = (state_j[2] - state_i[2])/ delta_t;

        // est_pseudorange = sqrt(delta_x+ delta_y + delta_z);

        auto switch_act_func = [](T val) {
            return T(1.0);
            // if (val < T(0.0)) return T(0.0);
            // if (val > T(1.0)) return T(1.0);
            // return val;
        };

        // residuals[0] = switch_act_func(switchvar[0]) * (est_v_x - T(v_x)) / T(var_vector(0));
        // residuals[1] = switch_act_func(switchvar[0]) * (est_v_y - T(v_y)) / T(var_vector(1));
        // residuals[2] = switch_act_func(switchvar[0]) * (est_v_z - T(v_z)) / T(var_vector(2));

        residuals[0] = switch_act_func(switchvar[0]) * (est_v_x - T(vel_state[0])) / T(var_vector(0));
        residuals[1] = switch_act_func(switchvar[0]) * (est_v_y - T(vel_state[1])) / T(var_vector(1));
        residuals[2] = switch_act_func(switchvar[0]) * (est_v_z - T(vel_state[2])) / T(var_vector(2));

        return true;
    }

    double v_x, v_y, v_z, var;
    Eigen::Vector3d var_vector;
    double delta_t;
    std::string sat_sys; // satellite system
};

struct accelFactor
{
    accelFactor(double delta_t): delta_t(delta_t) {}

    template <typename T>
    bool operator()(const T* vel_state_i, const T* vel_state_j, T* residuals) const
    {
        T est_a_x = (vel_state_j[0] - vel_state_i[0])/ delta_t;
        T est_a_y = (vel_state_j[1] - vel_state_i[1])/ delta_t;
        T est_a_z = (vel_state_j[2] - vel_state_i[2])/ delta_t;

        T est_a_sqr = est_a_x*est_a_x + est_a_y*est_a_y + est_a_z*est_a_z;

        // T thresh = T(5)*T(5);
        T thresh = T(10)*T(10);
        residuals[0] = est_a_sqr > thresh ? (est_a_sqr-thresh)*1e2 : T(0); //est_a_sqr);
        // residuals[0] = est_a_sqr;

        T vel_sqr = vel_state_i[0]*vel_state_i[0]+vel_state_i[1]*vel_state_i[1]+vel_state_i[2]*vel_state_i[2];
        // residuals[1] = vel_sqr < T(1*1) ? (vel_sqr*1e2) : T(0);

        residuals[1] = T(0);

        return true;
    }

    double delta_t;
};

/* doppler factor*/
struct dopplerConstraint
{
    typedef ceres::DynamicAutoDiffCostFunction<dopplerConstraint, 10>
      dopplerDynaCostFunction;
    dopplerConstraint(double v_x, double v_y, double v_z, double delta_t, Eigen::Vector3d var_vector, int keyIndex)
                :v_x(v_x),v_y(v_y), v_z(v_z), delta_t(delta_t), var_vector(var_vector),keyIndex(keyIndex){}

    template <typename T>
    bool operator()(T const* const* state, T* residuals) const
    {
        /* state[keyIndex][0]   -> state_i
           state[keyIndex+1][0] -> state_j 
        */
        T est_v_x = (state[keyIndex+1][0] - state[keyIndex][0])/ delta_t;
        T est_v_y = (state[keyIndex+1][1] - state[keyIndex][1])/ delta_t;
        T est_v_z = (state[keyIndex+1][2] - state[keyIndex][2])/ delta_t;

        // est_pseudorange = sqrt(delta_x+ delta_y + delta_z);

        residuals[0] = (est_v_x - T(v_x)) / T(var_vector(0));
        residuals[1] = (est_v_y - T(v_y)) / T(var_vector(1));
        residuals[2] = (est_v_z - T(v_z)) / T(var_vector(2));

        return true;
    }

    // Factory method to create a CostFunction from a DDCarrierPhaseConstraint to
    // conveniently add to a ceres problem.
    static dopplerDynaCostFunction* Create(double v_x, double v_y, double v_z, double delta_t, Eigen::Vector3d var_vector, int keyIndex, std::vector<double*>* state_array, std::vector<double*>* pose_parameter_blocks, std::vector<int*>* ar_state_num) {
        
        dopplerConstraint* constraint = new dopplerConstraint(
            v_x, v_y, v_z, delta_t, var_vector, keyIndex);
        
        dopplerDynaCostFunction* cost_function = new dopplerDynaCostFunction(constraint);
        
        pose_parameter_blocks->clear();
        // double a[5] = {1,2,3,4,5};
        // parameter_blocks->push_back(a);
        // parameter_blocks->push_back(&((*state_array)[keyIndex]));
        // parameter_blocks->push_back(state_array[keyIndex]);
        
        /* push back one more state ("+1") corresponding to state_j */
        for(int i = 0; i <(keyIndex+1+1); i++)
        {
            pose_parameter_blocks->push_back((*state_array)[i]);
            cost_function->AddParameterBlock(3 + (*ar_state_num)[i][0]);
        }
        // std::cout << "parameter_blocks.size()-> " << parameter_blocks->size() << std::endl;
        // cost_function->AddParameterBlock(1);
        // cost_function->AddParameterBlock(5);
        
        cost_function->SetNumResiduals(3);
        return (cost_function);
    }

    double v_x, v_y, v_z, var;
    Eigen::Vector3d var_vector;
    double delta_t;
    std::string sat_sys; // satellite system
    int keyIndex;
};