/*******************************************************
 * Copyright (C) 2019, Intelligent Positioning and Navigation Lab, Hong Kong Polytechnic University
 * 
 * This file is part of GraphGNSSLib.
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Weisong Wen (weisong.wen@connect.polyu.hk)
 * Function: decode the RINEX file, output GNSS raw measurements via ros topics
 * Date: 2020/11/27
 *******************************************************/

#include <stdarg.h>
#include "../../RTKLIB/src/rtklib.h"
#include <ros/ros.h>
#include <stdio.h>
#include <assert.h>
#define ENACMP 1 // enable BeiDou

extern void postposRegisterPub(ros::NodeHandle &n);
extern void rtkposRegisterPub(ros::NodeHandle &n);
extern void pntposRegisterPub(ros::NodeHandle &n);

int main(int argc, char **argv)
{
	
    ros::init(argc, argv, "gnss_preprocessor_node");
	ros::NodeHandle nh("~");
	ROS_INFO("\033[1;32m----> gnss_preprocessor Started.\033[0m");

	/* get setup parameters */
	int mode, nf, soltype;
	std::string roverMeasureFile, baseMeasureFile, navFile;
	std::string out_folder;
	nh.param("mode",   mode, 2);
	nh.param("nf",     nf, 2);
	nh.param("soltype",soltype, 0);
	ros::param::get("roverMeasureFile", roverMeasureFile);
	ros::param::get("baseMeasureFile", baseMeasureFile);
	ros::param::get("navFile", navFile);
	ros::param::get("out_folder", out_folder);

	double stationParam[3];
	ros::param::get("station_x", stationParam[0]);
	ros::param::get("station_y", stationParam[1]);
	ros::param::get("station_z", stationParam[2]);

	station_x = stationParam[0];
	station_y = stationParam[1];
	station_z = stationParam[2];

	// std::cout << roverMeasureFile << ", " << baseMeasureFile << ", " << navFile << std::endl;
	// std::cout << station_x << ", " << station_y << ", " << station_z << std::endl;

	/* flag for state */
    int n=0,i,stat;

	/* input node handle */
	postposRegisterPub(nh);
	rtkposRegisterPub(nh);
	pntposRegisterPub(nh);

	/* processing time setting */
	double ti=0.0;						// processing interval  (s) (0:all)
	double tu=0.0;						// unit time (s) (0:all)
	gtime_t ts={0},te={0};
	ts.time=0;							// start time (0:all)
	te.time=0;							// end time (0:all)

	/* options */
	prcopt_t prcopt = prcopt_default;	// processing option
	solopt_t solopt = solopt_default;	// output solution option
	filopt_t filopt = {""};	            // file option
	prcopt.mode = mode;			// Kinematic RTK
	// prcopt.mode = PMODE_SINGLE;			// SPP
	prcopt.navsys = SYS_GPS | SYS_GAL | SYS_GLO;              // use all satellites system
	// prcopt.navsys = SYS_ALL;
	// prcopt.navsys = SYS_GPS | SYS_CMP; 
	// prcopt.navsys = SYS_GPS;
	prcopt.nf = nf;						// frequency (1:L1,2:L1+L2,3:L1+L2+L5) 
	prcopt.soltype = soltype;					// 0:forward,1:backward,2:combined
	prcopt.elmin = 15.0*D2R;				// elevation mask (rad)	

	prcopt.snrmask.ena[0] = 1;
	prcopt.snrmask.ena[1] = 1;
	for (int j = 0; j < 9; ++j) prcopt.snrmask.mask[0][j] = 24;
	for (int j = 0; j < 9; ++j) prcopt.snrmask.mask[1][j] = 34;
	for (int j = 0; j < 9; ++j) prcopt.snrmask.mask[2][j] = 24;

	prcopt.tidecorr = 0;					// earth tide correction (0:off,1-:on) 
	prcopt.posopt[4] = 0;               // use RAIM FDE (qmo)  1
	prcopt.tropopt = TROPOPT_SAAS;        // troposphere option: Saastamoinen model
	prcopt.ionoopt = IONOOPT_BRDC;		// ionosphere option: Broad cast
	prcopt.sateph = EPHOPT_BRDC;			// ephemeris option: broadcast ephemeris

	prcopt.posopt[0] = 0;
	prcopt.posopt[1] = 0;
	prcopt.posopt[2] = 0;
	prcopt.posopt[3] = 0;
	prcopt.posopt[4] = 0;
	prcopt.posopt[5] = 0;

	prcopt.dynamics = 1;

	prcopt.modear = 0;					// AR mode (0:off,1:continuous,2:instantaneous,3:fix and hold)
	prcopt.glomodear = 0;
	prcopt.bdsmodear = 0;
	prcopt.arfilter = 1;

	prcopt.thresar[0] = 3;
	prcopt.thresar[5] = 3;
	prcopt.thresar[6] = 3;
	prcopt.thresar[1] = 0.05;
	prcopt.thresar[2] = 0;
	prcopt.thresar[3] = 1e-09;
	prcopt.thresar[4] = 1e-05;

	prcopt.varholdamb = 0.1;
	prcopt.gainholdamb = 0.01;

	prcopt.minlock = 5;
	prcopt.minfixsats = 4;
	prcopt.minholdsats = 5;
	prcopt.mindropsats = 10;
	prcopt.minfix = 10;
	prcopt.armaxiter = 1;
	prcopt.elmaskar = 15;
	prcopt.elmaskhold = 15;
	prcopt.maxout = 4;
	prcopt.maxtdiff = 30;
	prcopt.syncsol = 0;
	prcopt.thresslip = 0.1;
	prcopt.thresdop = 5;
	prcopt.maxinno = 1;
	prcopt.maxgdop = 30;
	prcopt.niter = 1;

	prcopt.eratio[0] = 300;
	prcopt.eratio[1] = 300;
	prcopt.eratio[2] = 300;
	prcopt.err[1] = 0; //0.003;
	prcopt.err[2] = 0; //0.003;
	prcopt.err[3] = 0;
	prcopt.err[4] = 1;
	prcopt.err[5] = 52;
	prcopt.err[6] = 0.003;
	prcopt.err[7] = 0.1;
	
	prcopt.std[0] = 30;
	prcopt.std[1] = 0.03;
	prcopt.std[2] = 0.3;

	prcopt.prn[3] = 3;
	prcopt.prn[4] = 1;
	prcopt.prn[0] = 0.01;
	prcopt.prn[1] = 0.001;
	prcopt.prn[2] = 0.0001;
	prcopt.prn[5] = 0;
	prcopt.sclkstab = 5e-12;

	solopt.outopt = 1;					// output processing options (0:no,1:yes)
	solopt.timef = 0;						// time format (0:sssss.s,1:yyyy/mm/dd hh:mm:ss.s)
	solopt.timeu = 3;						// time digits under decimal point
	solopt.sep[0] = ',';					// field separator
	solopt.sstat= 0;						// solution statistics level (0:off,1:states,2:residuals)
	solopt.trace = 0;						// debug trace level (0:off,1-5:debug)
	solopt.sstat = 0;						// get the solution file
	solopt.posf = SOLF_LLH;
	solopt.height = 0;

	char *rov="",*base="";
	char infile_[10][1024]={""}, *infile[10];
	char outfile[1024];

	/* set input files */
	for (i=0;i<10;i++) infile[i]=infile_[i];

	strcpy(infile[n++],strdup(roverMeasureFile.c_str()));
	strcpy(infile[n++],strdup(baseMeasureFile.c_str()));
	strcpy(infile[n++],strdup(navFile.c_str()));
	// strcpy(infile[n++],strdup(GPSEmpFile.c_str()));

	/* if you use the RTK mode, specify the position of the station (only used by RTKLIB)
	 * following is an example position of the base HKSC in Hong Kong */
	prcopt.rb[0] = station_x;			// base position for relative mode {x,y,z} (ecef) (m)
	prcopt.rb[1] = station_y;			// base position for relative mode {x,y,z} (ecef) (m)
	prcopt.rb[2] = station_z;			// base position for relative mode {x,y,z} (ecef) (m)

	/* set output files */
	strcpy(outfile, strdup(out_folder.c_str()));

	/* decode the RINEX files*/
	// while(ros::ok())
	{
		/* decode the RINEX file positioning */
		stat=postpos(ts,te,ti,tu,&prcopt,&solopt,&filopt,infile,n,outfile,rov,base);

		printf("\n");
		if(stat==0){
			ROS_INFO("\033[1;32m----> gnss_preprocessor Finished.\033[0m");
		}
		else if(stat>0){
			ROS_INFO("\033[1;32m----> gnss_preprocessor Error!!!.\033[0m");
		}
		else if(stat==-1){
			ROS_INFO("\033[1;32m----> gnss_preprocessor Aborted!!!.\033[0m");
		}
		// system("pause");
	}
	// while(ros::ok())
	// {}
    // ros::spin();
	ros::spinOnce();
    return 0;
}

/*****dummy application functions for shared library*****/
extern int showmsg(char *format,...) {
	va_list arg;
	char buff[1024];
	if (*format) {
		va_start(arg,format);
		vsprintf(buff,format,arg);
		va_end(arg);
		printf("%s\n",buff);
	}
	return 0;	
}
extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}
