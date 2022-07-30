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


struct CarrierPhaseDeltaAdrFactor
{
    CarrierPhaseDeltaAdrFactor(
        double delta_t,
        const nlosExclusion::GNSS_Raw& i_prev,  
        const nlosExclusion::GNSS_Raw& i_curr,
        const Eigen::Vector3d& i_prev_los, 
        const Eigen::Vector3d& i_curr_los, 
        double lambda,
        double var):
            delta_t(delta_t),
            d_prev(i_prev),d_curr(i_curr),
            d_prev_los(i_prev_los), 
            d_curr_los(i_curr_los),
            d_lambda(lambda),
            dd_var(var)
    {}

    template <typename T>
    bool operator()(const T* prev_state, const T* state, const T* prev_velocity_state, const T* switchvar, T* residuals) const
    {
        double lambda = d_lambda;

        T delta_adr = T(d_curr.carrier_phase - d_prev.carrier_phase) * lambda;

        T delta_d = T(d_curr_los.dot(Eigen::Vector3d(d_curr.sat_pos_x, d_curr.sat_pos_y, d_curr.sat_pos_z))
         - d_prev_los.dot(Eigen::Vector3d(d_prev.sat_pos_x, d_prev.sat_pos_y, d_prev.sat_pos_z)));

        T delta_g = (d_curr_los.x() * prev_state[0] + d_curr_los.y() * prev_state[1] + d_curr_los.z() * prev_state[2])
            - (d_prev_los.x() * prev_state[0] + d_prev_los.y() * prev_state[1] + d_prev_los.z() * prev_state[2]);

        T delta_eru = 
                    (d_curr_los.x() * prev_velocity_state[0]
                     + d_curr_los.y() * prev_velocity_state[1]
                     + d_curr_los.z() * prev_velocity_state[2]) * T(delta_t);

        T delta_rbias = (state[3] - prev_state[3]);

        auto switch_act_func = [](T val) {
            return T(1.0);
            // if (val < T(0.0)) return T(0.0);
            // if (val > T(1.0)) return T(1.0);
            // return T(val);
            // return val;
            if (val>T(30))
                return T(1.0);
            if (val<T(-30))
                return T(0.0);  
            return T(1.0) / (T(1.0) + exp(-T(val)));
            // return T(1);
            // return val;
        };


        residuals[0] = switch_act_func(switchvar[0]) * 
            T(delta_adr - delta_d + delta_g + delta_eru - delta_rbias) / T(1e-2*dd_var + 1e-7);
        residuals[0] += 1e2 * abs(T(switch_act_func(switchvar[0]) - T(1.0)));

        // residuals[0] = T(delta_move - delta_adr) / T(1e-1 + 1e-7);

        // printf("p,%f,%f,%f,%f,%f,%f\n", delta_move, rover_delta, sat_delta, delta_adr, dd_var, residuals[0]);

        return true;
    }

    double delta_t;
    nlosExclusion::GNSS_Raw d_prev, d_curr;
    Eigen::Vector3d d_prev_los, d_curr_los;
    double dd_var = 1.0;
    double d_lambda = 1.0;
};

#if 0
struct CarrierPhaseDeltaAdrFactor
{
    CarrierPhaseDeltaAdrFactor(const nlosExclusion::GNSS_Raw& i_prev,  
        const nlosExclusion::GNSS_Raw& i_curr,
        const Eigen::Vector3d& i_dir, double lambda,
        double var):d_prev(i_prev),d_curr(i_curr), d_lambda(lambda),
            d_dir(i_dir), dd_var(var)
    {}

    template <typename T>
    bool operator()(const T* prev_state, const T* state, T* residuals) const
    {
        double OMGE_ = 7.2921151467E-5;
        double CLIGHT_ = 299792458.0;

        T prev_dist;
        T delta_x1 = pow((prev_state[0] - d_prev.sat_pos_x),2);
        T delta_y1 = pow((prev_state[1] - d_prev.sat_pos_y),2);
        T delta_z1 = pow((prev_state[2] - d_prev.sat_pos_z),2);
        prev_dist = sqrt(delta_x1+ delta_y1 + delta_z1);
        prev_dist = prev_dist + OMGE_ * (d_prev.sat_pos_x*prev_state[1]-d_prev.sat_pos_y*prev_state[0])/CLIGHT_;
        
        T curr_dist;
        delta_x1 = pow((state[0] - d_curr.sat_pos_x),2);
        delta_y1 = pow((state[1] - d_curr.sat_pos_y),2);
        delta_z1 = pow((state[2] - d_curr.sat_pos_z),2);
        curr_dist = sqrt(delta_x1+ delta_y1 + delta_z1);
        curr_dist = curr_dist + OMGE_ * (d_curr.sat_pos_x*state[1]-d_curr.sat_pos_y*state[0])/CLIGHT_;

        T delta_move = curr_dist - prev_dist;
        
        T delta_adr = T(d_curr.carrier_phase - d_prev.carrier_phase) * d_lambda;

        residuals[0] = T(delta_move - delta_adr) / T(1e-3*dd_var + 1e-7);

        // printf("p,%f,%f,%f,%f,%f,%f\n", delta_move, rover_delta, sat_delta, delta_adr, dd_var, residuals[0]);

        return true;
    }

    nlosExclusion::GNSS_Raw d_prev, d_curr;
    Eigen::Vector3d d_dir;
    double dd_var = 1.0;
    double d_lambda = 1.0;
};
#endif


#if 0
struct CarrierPhaseDeltaAdrFactor
{
    CarrierPhaseDeltaAdrFactor(const nlosExclusion::GNSS_Raw& i_prev, 
        const nlosExclusion::GNSS_Raw& i_curr,
        const nlosExclusion::GNSS_Raw& i_ref_prev,
        const nlosExclusion::GNSS_Raw& i_ref_curr)
            :d_prev(i_prev),d_curr(i_curr),
            d_ref_prev(i_ref_prev), d_ref_curr(i_ref_curr),
            d_var(i_curr.phase_var+i_ref_curr.phase_var),
            d_lambda(i_ref_prev.lamda)
    {}

    template <typename T>
    bool operator()(const T* prev_state, const T* state, T* residuals) const
    {
        double OMGE_ = 7.2921151467E-5;
        double CLIGHT_ = 299792458.0;

        T prev_dist0;
        T delta_x1 = pow((prev_state[0] - d_prev.sat_pos_x),2);
        T delta_y1 = pow((prev_state[1] - d_prev.sat_pos_y),2);
        T delta_z1 = pow((prev_state[2] - d_prev.sat_pos_z),2);
        prev_dist0 = sqrt(delta_x1+ delta_y1 + delta_z1);
        prev_dist0 += OMGE_ * (d_prev.sat_pos_x*prev_state[1]-d_prev.sat_pos_y*prev_state[0])/CLIGHT_;
        
        T prev_dist1;
        delta_x1 = pow((prev_state[0] - d_ref_prev.sat_pos_x),2);
        delta_y1 = pow((prev_state[1] - d_ref_prev.sat_pos_y),2);
        delta_z1 = pow((prev_state[2] - d_ref_prev.sat_pos_z),2);
        prev_dist1 = sqrt(delta_x1+ delta_y1 + delta_z1);
        prev_dist1 += OMGE_ * (d_ref_prev.sat_pos_x*prev_state[1]-d_ref_prev.sat_pos_y*prev_state[0])/CLIGHT_;

        T curr_dist0;
        delta_x1 = pow((state[0] - d_curr.sat_pos_x),2);
        delta_y1 = pow((state[1] - d_curr.sat_pos_y),2);
        delta_z1 = pow((state[2] - d_curr.sat_pos_z),2);
        curr_dist0 = sqrt(delta_x1+ delta_y1 + delta_z1);
        curr_dist0 += OMGE_ * (d_curr.sat_pos_x*state[1]-d_curr.sat_pos_y*state[0])/CLIGHT_;

        T curr_dist1;
        delta_x1 = pow((state[0] - d_ref_curr.sat_pos_x),2);
        delta_y1 = pow((state[1] - d_ref_curr.sat_pos_y),2);
        delta_z1 = pow((state[2] - d_ref_curr.sat_pos_z),2);
        curr_dist1 = sqrt(delta_x1+ delta_y1 + delta_z1);
        curr_dist1 += OMGE_ * (d_ref_curr.sat_pos_x*state[1]-d_ref_curr.sat_pos_y*state[0])/CLIGHT_;

        T delta_move = (curr_dist1 - curr_dist0) - (prev_dist1 - prev_dist0);
        
        T delta_adr = T(d_ref_curr.carrier_phase - d_curr.carrier_phase) * d_lambda
            - T(d_ref_prev.carrier_phase - d_prev.carrier_phase) * d_lambda;

        residuals[0] = T(delta_move - delta_adr) / T(1e-6); // T(1e-8*d_var + 1e-7);
        // printf("p,%f,%f,%f,%f,%f,%f\n", delta_move, rover_delta, sat_delta, delta_adr, dd_var, residuals[0]);

        return true;
    }

    nlosExclusion::GNSS_Raw d_prev, d_curr;
    nlosExclusion::GNSS_Raw d_ref_prev, d_ref_curr;
    double d_var = 1.0;
    double d_lambda = 1.0;
};
#endif

/* DDCarrierPhaseFactor factor (fix the bug: repeating estimation of the ambiguity between master satellite and user, station )
* fix the bug: only estimate the double-differenced carrier-phase bias
* M: master satellite
* i: satellite to be double-difference based on master satellite
* r: reference satellite
* u: use end satellite (GNSS receiver)
* state: 
*   0~4: x,y,z, clc_gps, clc_beidou
*   5~8: gps master satellite to reference, user. beidou master satellite to reference, user. 
*/
struct DDCarrierPhaseFactor_DDBias
{
    DDCarrierPhaseFactor_DDBias(DDMeasurement dd_measurement, Eigen::Vector3d base_pos, int index, double var)
                :dd_measurement(dd_measurement),base_pos(base_pos),index(index),dd_var(var) {}

    template <typename T>
    bool operator()(const T* state, const T* ar_state, T* residuals) const
    {
        /* get the position of station*/
        Eigen::Vector3d pose_r = base_pos;
        GNSS_Tools m_GNSS_Tools1; // utilities
        double OMGE_ = 7.2921151467E-5;
        double CLIGHT_ = 299792458.0;

        double lambda = dd_measurement.r_master_SV.lamda;
        // lambda = 0.192039;
        // LOG(INFO) << "CP lambda-> " << lambda;

        /* satellite position*/
        Eigen::Vector3d u_pose_m(dd_measurement.u_master_SV.sat_pos_x, dd_measurement.u_master_SV.sat_pos_y,dd_measurement.u_master_SV.sat_pos_z);
        double var_u2m = m_GNSS_Tools1.getVarofCp_ele_SNR(dd_measurement.u_master_SV);

        Eigen::Vector3d u_pose_i(dd_measurement.u_iSV.sat_pos_x, dd_measurement.u_iSV.sat_pos_y,dd_measurement.u_iSV.sat_pos_z); 
        double var_u2i = m_GNSS_Tools1.getVarofCp_ele_SNR(dd_measurement.u_iSV);

        Eigen::Vector3d r_pose_m(dd_measurement.r_master_SV.sat_pos_x, dd_measurement.r_master_SV.sat_pos_y,dd_measurement.r_master_SV.sat_pos_z);
        double var_r2m = m_GNSS_Tools1.getVarofCp_ele_SNR(dd_measurement.r_master_SV);
        
        Eigen::Vector3d r_pose_i(dd_measurement.r_iSV.sat_pos_x, dd_measurement.r_iSV.sat_pos_y,dd_measurement.r_iSV.sat_pos_z); 
        double var_r2i = m_GNSS_Tools1.getVarofCp_ele_SNR(dd_measurement.r_iSV);

        T est_p_r2m;
        est_p_r2m = T(m_GNSS_Tools1.getDistanceFrom2Points(pose_r, r_pose_m));
        est_p_r2m = est_p_r2m + OMGE_ * (r_pose_m(0)*pose_r(1)-r_pose_m(1)*pose_r(0))/CLIGHT_;

        T est_p_r2i = T(m_GNSS_Tools1.getDistanceFrom2Points(pose_r, r_pose_i));
        est_p_r2i = est_p_r2i + OMGE_ * (r_pose_i(0)*pose_r(1)-r_pose_i(1)*pose_r(0))/CLIGHT_;

        T est_p_u2m;
        T delta_x1 = pow((state[0] - u_pose_m(0)),2);
        T delta_y1 = pow((state[1] - u_pose_m(1)),2);
        T delta_z1 = pow((state[2] - u_pose_m(2)),2);
        est_p_u2m = sqrt(delta_x1+ delta_y1 + delta_z1);
        est_p_u2m = est_p_u2m + OMGE_ * (u_pose_m(0)*state[1]-u_pose_m(1)*state[0])/CLIGHT_;
        
        T est_p_u2i;
        T delta_x2 = pow((state[0] - u_pose_i(0)),2);
        T delta_y2 = pow((state[1] - u_pose_i(1)),2);
        T delta_z2 = pow((state[2] - u_pose_i(2)),2);
        est_p_u2i = sqrt(delta_x2+ delta_y2 + delta_z2);
        est_p_u2i = est_p_u2i + OMGE_ * (u_pose_i(0)*state[1]-u_pose_i(1)*state[0])/CLIGHT_;

        T est_DD_cp = (est_p_u2i - est_p_r2i) - (est_p_u2m - est_p_r2m);
        est_DD_cp = est_DD_cp + ar_state[index] * lambda;   

        T p_r2m = T(dd_measurement.r_master_SV.carrier_phase * lambda);
        T p_r2i = T(dd_measurement.r_iSV.carrier_phase * lambda);
        T p_u2m = T(dd_measurement.u_master_SV.carrier_phase * lambda);
        T p_u2i = T(dd_measurement.u_iSV.carrier_phase * lambda);
        T DD_cp = (p_u2i - p_r2i) - (p_u2m - p_r2m);

        double var = var_u2m + var_u2i + var_r2m + var_r2i;
        var = var /4.0;

        var = dd_var;
        // if ((dd_measurement.u_master_SV.visable & 1) || (dd_measurement.u_iSV.visable & 1))
        //     var = 1e20;

        #define use_fixed_cov 1
        #if use_fixed_cov
        residuals[0] = (est_DD_cp - DD_cp) / T(1e-4*var + 1e-7); //T(pow(0.004, 2)); 
        // residuals[0] = (est_DD_cp - DD_cp) / T(var + 1e-7); //T(pow(0.004, 2)); 
        // printf("cp: %f,%f,%f,%f,%f\n", est_DD_cp, DD_cp, (est_DD_cp - DD_cp), var, residuals[0]);
        // printf("cp: %f,%f,%f\n", est_DD_cp, DD_cp, (est_DD_cp - DD_cp));
        // residuals[0] = (est_DD_cp - DD_cp) / T(pow(0.004, 2)); 
        #else
        // residuals[0] = (est_DD_cp - DD_cp) / T(var_u2i);
        residuals[0] = (est_DD_cp - DD_cp) / T(var);
        #endif

        return true;
    }

    DDMeasurement dd_measurement;
    Eigen::Vector3d base_pos;
    int index; // index of the ambiguity in one epoch
    double dd_var = 1;
};


/* DDCarrierPhaseConstraint with DynamicAutoDiffCostFunction from Ceres-solver (fix the bug: repeating estimation of the ambiguity between master satellite and user, station )
* fix the bug: only estimate the double-differenced carrier-phase bias
* M: master satellite
* i: satellite to be double-difference based on master satellite
* r: reference satellite
* u: use end satellite (GNSS receiver)
* state: 
*   0~4: x,y,z, clc_gps, clc_beidou
*   5~8: gps master satellite to reference, user. beidou master satellite to reference, user. 
*/
struct DDCarrierPhaseConstraint
{
    typedef ceres::DynamicAutoDiffCostFunction<DDCarrierPhaseConstraint, 10>
      DDCarrierPhaseDynaCostFunction;
    DDCarrierPhaseConstraint(DDMeasurement dd_measurement, Eigen::Vector3d base_pos, int keyIndex, int index)
                :dd_measurement(dd_measurement),base_pos(base_pos),keyIndex(keyIndex),index(index){}

    template <typename T>
    bool operator()(T const* const* state,  T* residuals) const
    {
        /* get the position of station*/
        Eigen::Vector3d pose_r = base_pos;
        GNSS_Tools m_GNSS_Tools1; // utilities
        double OMGE_ = 7.2921151467E-5;
        double CLIGHT_ = 299792458.0;

        double lambda = dd_measurement.r_master_SV.lamda;
        // lambda = 0.192039;
        // LOG(INFO) << "CP lambda-> " << lambda;

        /* satellite position*/
        Eigen::Vector3d u_pose_m(dd_measurement.u_master_SV.sat_pos_x, dd_measurement.u_master_SV.sat_pos_y,dd_measurement.u_master_SV.sat_pos_z);
        double var_u2m = m_GNSS_Tools1.getVarofCp_ele_SNR(dd_measurement.u_master_SV);

        Eigen::Vector3d u_pose_i(dd_measurement.u_iSV.sat_pos_x, dd_measurement.u_iSV.sat_pos_y,dd_measurement.u_iSV.sat_pos_z); 
        double var_u2i = m_GNSS_Tools1.getVarofCp_ele_SNR(dd_measurement.u_iSV);

        Eigen::Vector3d r_pose_m(dd_measurement.r_master_SV.sat_pos_x, dd_measurement.r_master_SV.sat_pos_y,dd_measurement.r_master_SV.sat_pos_z);
        double var_r2m = m_GNSS_Tools1.getVarofCp_ele_SNR(dd_measurement.r_master_SV);
        
        Eigen::Vector3d r_pose_i(dd_measurement.r_iSV.sat_pos_x, dd_measurement.r_iSV.sat_pos_y,dd_measurement.r_iSV.sat_pos_z); 
        double var_r2i = m_GNSS_Tools1.getVarofCp_ele_SNR(dd_measurement.r_iSV);

        T est_p_r2m;
        est_p_r2m = T(m_GNSS_Tools1.getDistanceFrom2Points(pose_r, r_pose_m));
        est_p_r2m = est_p_r2m + OMGE_ * (r_pose_m(0)*pose_r(1)-r_pose_m(1)*pose_r(0))/CLIGHT_;

        T est_p_r2i = T(m_GNSS_Tools1.getDistanceFrom2Points(pose_r, r_pose_i));
        est_p_r2i = est_p_r2i + OMGE_ * (r_pose_i(0)*pose_r(1)-r_pose_i(1)*pose_r(0))/CLIGHT_;

        T est_p_u2m;
        T delta_x1 = pow((state[keyIndex][0] - u_pose_m(0)),2);
        T delta_y1 = pow((state[keyIndex][1] - u_pose_m(1)),2);
        T delta_z1 = pow((state[keyIndex][2] - u_pose_m(2)),2);
        est_p_u2m = sqrt(delta_x1+ delta_y1 + delta_z1);
        est_p_u2m = est_p_u2m + OMGE_ * (u_pose_m(0)*state[keyIndex][1]-u_pose_m(1)*state[keyIndex][0])/CLIGHT_;
        
        T est_p_u2i;
        T delta_x2 = pow((state[keyIndex][0] - u_pose_i(0)),2);
        T delta_y2 = pow((state[keyIndex][1] - u_pose_i(1)),2);
        T delta_z2 = pow((state[keyIndex][2] - u_pose_i(2)),2);
        est_p_u2i = sqrt(delta_x2+ delta_y2 + delta_z2);
        est_p_u2i = est_p_u2i + OMGE_ * (u_pose_i(0)*state[keyIndex][1]-u_pose_i(1)*state[keyIndex][0])/CLIGHT_;
        
        /* keyIndex: epoch index, index: ambiguity index in one epoch */
        T est_DD_cp = (est_p_u2i - est_p_r2i) - (est_p_u2m - est_p_r2m);
        est_DD_cp = est_DD_cp + state[keyIndex][3 + index] * lambda;   

        T p_r2m = T(dd_measurement.r_master_SV.carrier_phase * lambda);
        T p_r2i = T(dd_measurement.r_iSV.carrier_phase * lambda);
        T p_u2m = T(dd_measurement.u_master_SV.carrier_phase * lambda);
        T p_u2i = T(dd_measurement.u_iSV.carrier_phase * lambda);
        T DD_cp = (p_u2i - p_r2i) - (p_u2m - p_r2m);

        double var = var_u2m + var_u2i + var_r2m + var_r2i;
        var = var /4.0;

        #if use_fixed_cov
        residuals[0] = (est_DD_cp - DD_cp) / T(pow(0.004, 2)); 
        #else
        // residuals[0] = (est_DD_cp - DD_cp) / T(var_u2i);
        residuals[0] = (est_DD_cp - DD_cp) / T(var);
        #endif

        return true;
    }

    // Factory method to create a CostFunction from a DDCarrierPhaseConstraint to
    // conveniently add to a ceres problem.
    static DDCarrierPhaseDynaCostFunction* Create(DDMeasurement dd_measurement, Eigen::Vector3d base_pos, int keyIndex, int index, std::vector<double*>* state_array, std::vector<double*>* pose_parameter_blocks, std::vector<int*>* ar_state_num) {
        
        DDCarrierPhaseConstraint* constraint = new DDCarrierPhaseConstraint(
            dd_measurement, base_pos, keyIndex, index);
        
        DDCarrierPhaseDynaCostFunction* cost_function = new DDCarrierPhaseDynaCostFunction(constraint);
        
        pose_parameter_blocks->clear();
        // double a[5] = {1,2,3,4,5};
        // parameter_blocks->push_back(a);
        // parameter_blocks->push_back(&((*state_array)[keyIndex]));
        // parameter_blocks->push_back(state_array[keyIndex]);
        
        for(int i = 0; i <(keyIndex+1); i++)
        {
            pose_parameter_blocks->push_back((*state_array)[i]);
            cost_function->AddParameterBlock(3 + (*ar_state_num)[i][0]);
        }
        // std::cout << "parameter_blocks.size()-> " << parameter_blocks->size() << std::endl;
        // cost_function->AddParameterBlock(1);
        // cost_function->AddParameterBlock(5);
        
        cost_function->SetNumResiduals(1);
        return (cost_function);
    }

    DDMeasurement dd_measurement;
    Eigen::Vector3d base_pos;
    int index; // index of the ambiguity in one epoch
    int keyIndex; // index of epoch

};

/**
	 * @brief find master satellites from user end gnss data,
	 * @param satellite id, gnss data array (one epoch), gnss data (one satellite)
	 * @return none
	 */
    bool findMasterSatellite(int sat_id, nlosExclusion::GNSS_Raw_Array user_gnss_data, 
        nlosExclusion::GNSS_Raw& master_sv, nlosExclusion::GNSS_Raw& i_sv,
        int freqidx_filter=-1)
    {   
        /* get navigation satellite system */
        int sys=satsys(sat_id,NULL);   
        if(sys==SYS_GPS)
        {
            // LOG(INFO) << "GPS Satellite  ";
        }
        else if(sys==SYS_CMP)
        {
            // LOG(INFO) << "BeiDou Satellite   ";
        }
        else
        {
            // LOG(INFO) << "Unknow!!!!! Satellite   ";
        }
        
        int sv_cnt = user_gnss_data.GNSS_Raws.size(); 
        int same_satsystem_cnt = 0;
        double master_sv_id = -1;
        double max_elevation = -1;
        for(int i = 0; i < sv_cnt; i++)
        {
            if (freqidx_filter != -1 && user_gnss_data.GNSS_Raws[i].freqidx != freqidx_filter) continue;

            if(sys == (satsys(user_gnss_data.GNSS_Raws[i].prn_satellites_index,NULL))) // same satellite system
            {
                same_satsystem_cnt++;
                if(user_gnss_data.GNSS_Raws[i].elevation >= max_elevation)
                {
                    max_elevation = user_gnss_data.GNSS_Raws[i].elevation;
                    master_sv_id = i;
                }

                // same satellite system and same satellite id
                if(sat_id ==int(user_gnss_data.GNSS_Raws[i].prn_satellites_index))
                {
                    i_sv = user_gnss_data.GNSS_Raws[i];
                    i_sv.elevation = user_gnss_data.GNSS_Raws[i].elevation;
                    // if(i_sv.elevation==0)
                    // {
                    //     LOG(INFO) << "satellite with zero elevation angle---";
                    // }
                }
            }
        }

        if(same_satsystem_cnt>=2) // at least 2 satellites with same sat system
        {
            master_sv = user_gnss_data.GNSS_Raws[master_sv_id];
            // LOG(INFO) << "elevation of master satellite " << max_elevation;
            /*check if you have find the satellite with same ID in user end*/
            if(sat_id==master_sv.prn_satellites_index || (i_sv.pseudorange<10)) // the sat_id is same as master satellite (double-difference should not be done between master and master)
            {
                // LOG(INFO) << "Warning!!! master satellite is itself!!!!!";
                return false;
            }
            return true;
        }
        else 
            return false;

    }

    /**
	 * @brief find closest gnss measurement from user end gnss data
	 * @param time, gnss map from user end, gnss data (one epoch)
	 * @return none
	*/
    void findClosestEpoch(double t, std::map<double, nlosExclusion::GNSS_Raw_Array> gnss_raw_map, 
        nlosExclusion::GNSS_Raw_Array& cloest_epoch_gnss)
    {
        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator gnss_iter;
        gnss_iter = gnss_raw_map.begin();
        int length = gnss_raw_map.size();
        double time_diff = 100000;
        for(int i = 0;  i < length; i++,gnss_iter++) // initialize
        {
            if((fabs(t - gnss_iter->first))<time_diff)
            {
                // with smallest time difference
                time_diff = fabs(t - gnss_iter->first);
                cloest_epoch_gnss = gnss_iter->second;
            }
            

        }
        // LOG(INFO) << "time_diff " << time_diff;
    }

    /**
	 * @brief find satellite with same id
	 * @param id, gnss data in one epoch, one satellite
	 * @return none
	*/
    void findSatellitewithSameId(double id, nlosExclusion::GNSS_Raw_Array gnss_data, 
        nlosExclusion::GNSS_Raw& same_id_sv, int freqidx_filter=-1)
    {
        int length = gnss_data.GNSS_Raws.size();
        for(int i = 0; i < length; i++)
        {
            if (freqidx_filter != -1 && gnss_data.GNSS_Raws[i].freqidx != freqidx_filter) continue;
            if(gnss_data.GNSS_Raws[i].prn_satellites_index == id)
            {
                same_id_sv = gnss_data.GNSS_Raws[i];
            }
        }
    }

    /**
	 * @brief check the wheather all the consided satellite have carrier-phase
	 * @param gnss double-difference measurements
	 * @return true or false
	*/
    bool checkCarrierPhaseConsistency(DDMeasurement DD_measurement)
    {
        /* the master satellite is found from user end*/
        nlosExclusion::GNSS_Raw u_master_SV;
        nlosExclusion::GNSS_Raw u_iSV;

        nlosExclusion::GNSS_Raw r_master_SV;
        nlosExclusion::GNSS_Raw r_iSV;
        if(DD_measurement.u_master_SV.carrier_phase<10) return false;
        else if(DD_measurement.u_iSV.carrier_phase<10) return false;
        else if(DD_measurement.r_master_SV.carrier_phase<10) return false;
        else if(DD_measurement.r_iSV.carrier_phase<10) return false;
        else return true;
        
    }

    /**
	 * @brief check the wheather all the consided satellite have carrier-phase 
     * or have cycle slip?
	 * @param gnss double-difference measurements
	 * @return true or false
	*/
    bool checkCarrierPhaseConsistencyCycleSlip(DDMeasurement DD_measurement)
    {
        /* the master satellite is found from user end*/
        nlosExclusion::GNSS_Raw u_master_SV;
        nlosExclusion::GNSS_Raw u_iSV;

        nlosExclusion::GNSS_Raw r_master_SV;
        nlosExclusion::GNSS_Raw r_iSV;
        if((DD_measurement.u_master_SV.carrier_phase<10) || (DD_measurement.u_master_SV.visable>0)) return false;
        else if((DD_measurement.u_iSV.carrier_phase<10) || (DD_measurement.u_iSV.visable>0)) return false;

        else if((DD_measurement.r_master_SV.carrier_phase<10)||(DD_measurement.r_master_SV.visable>0)) return false;
        else if((DD_measurement.r_iSV.carrier_phase<10)||(DD_measurement.r_iSV.visable>0)) return false;
        else return true;
        
    }