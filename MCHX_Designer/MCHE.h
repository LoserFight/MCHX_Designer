#pragma once
#include"rei.h"
#include "segment.h"
#include "louverFin.h"
#include "microTube.h"
#include<cmath>
#include<vector>
#include<iostream>
#include<fstream>
#include<valarray>

//#include "matrix.h"
//#include "mclcppclass.h"
//#include <mclmcrrt.h>
//#include"NewtonToSolveGdp.h" // 修改为你自己的工程

struct air_con {
	double v_in;
	double T;
	double P;
	double w;
};

class MCHE
{
	int m_segment;
	Segment* m_HE;
	double m_q_air;
	double m_q_ref;
	double P_air_in;
	double P_ref_in;


	//double m_T_hot_in;
	//double m_T_cold_in;
	//double m_T_hot_out;
	//double m_T_cold_out;



	const double lambda_316L = 15.1;
	const double Sigma_316L = 117E6;
	const double Pa_design = 10E6;
	
	double Delta_p_air;
	double Delta_p_ref;

	vector<UINT> passes;

	int sum_tubes;
	microTube* m_tube;//microTube para
	louverFin* m_fin;



	refri refrigerant_in;
	refri refrigerant_out;

	refri refrigerant_tmp_in;
	refri refrigerant_tmp_out;


	//homogenous
	double a_B = 1.0;
	double a_n1=1.0;
	double a_n2 = 0.0;
	double a_n3 = 1.0;

	////lockhart-Martinelli
	//double a_B = 0.28;
	//double a_n1 = 0.64;
	//double a_n2 = 0.007;
	//double a_n3 = 0.36;

	double total_L;
	double total_W;
	double total_H;
	double total_volume;

	double Afr_air;
	double Ac_air;
	double Af;//area of fin;
	double At;//area of tube;
	double Dh_air;

	double Ac_tube;

public:
	double q_re;//mass flow rate kg/s
	air_con air_in;
	double Q_total;
	double DP_re;
	double DelSat;
	int blockeds = 1;
	int bbk = 1;
	double reChargeL;
	double reChargeV;
	double reChargeVL;
	double HE_volume_tube;
	double HE_volume_fin;



	MCHE(vector<UINT>& p,int segNum,const refri & refiIn,const air_con &airin);
	MCHE( vector<UINT>& p, int segNum, const refri& refiIn, const air_con& airin, const microTube& tdata, const louverFin& fdata);

	double linearInd(const double& a, const double& b, const double& ai, const double& bi, const double& i) {
		return a + (b - a) * (i - ai) / (bi - ai);
	}

	~MCHE();
	//void heatbalance(Segment* a);
	void geometryDesign();
	//bool check_P_design();
	//void set_volume();
	//void LMTD_cal_KA(bool (MCHE::* rec)(Segment&));

	bool solve(int i);
	double get_dp_air() { return Delta_p_air; }
	double get_dp_ref() { return Delta_p_ref; }
	double get_HE_L() { return total_L; }
	double get_HE_W() { return total_W; }
	double get_HE_H() { return total_H; }
	double get_HE_volume() { return total_volume; }
	//Correlations
	bool Gnielinski(Segment& a, const double G, double& h, double& f);

	bool Friedel(Segment& a, double G_r, double& dp, double& Re_L);
	bool Blasius(double Re, double& f);
	bool BlasiusOnly(double Re, double& f);

	bool condShah2016(Segment& a, double G_t, double& hi);

	//air
	bool WangChang_louver(double& f, double& ho);
	double WangChang_fin_eta(double h0);
	double classic(double h0);
	bool jacobi_louver(double& f, double& ho, const refri& a_in, const refri& a_out);
	double KeCal(double sig);
	double KcCal(double sig);

	double alpha(const refri & r);
	//eps-NTU
	double AirLessNTU(double Cx, double NTU);

	double NTsolveG(double R_T, double R_L,double R_a, double n1, double n2,double n3, double Dp,double G0);

	double TubeLessNTU(double Cx, double NTU);

	void savePressure();
	void saveResult();
};

