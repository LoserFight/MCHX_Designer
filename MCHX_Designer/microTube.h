#pragma once
class microTube
{
public:
	double tube_length;
	double tube_thickness;
	double tube_width;
	double vertical_space;
	double covered_lengh;
	
	int tube_ports;
	double port_width;
	double port_height;
	double port_Ac;
	double Dh;
	double P;
	double k= 237.000;


	microTube(double length,double t,double w,double p,double covered,double Pw,double Ph,int ports);
	void geoCalculate();
	void geoCalculateCircle();
	void specific(double cDh, double Ac);
	~microTube() {};

};

