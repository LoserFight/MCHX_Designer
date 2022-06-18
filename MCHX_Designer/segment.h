#pragma once
#include "rei.h"
struct Segment {

	double Q;
	double L;

	//refrigerant
	//double h_rei_in;//j/kgk
	double h_rei_out;//j/kgk
	//double P_rei_in;//pa
	double P_rei_out;//pa
	double m_s; //kg/s
	double quality;
	double Rs_La;
	double Rs_Tu;
	double Rs_a;
	//refri out;

	double Re;
	double G_refri;
	double x;
	double alpha;


	//air
	//double P_air_in;
	double P_air_out;
	double h_air_out;
	double T_air_out;
	double ht_air;
	double v_air;//m/s fr
	//double W;

	//Segment(const char* a):out(a){}

};