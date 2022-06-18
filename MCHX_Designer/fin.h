#pragma once
class fin
{
public:
	double FPI;
	double finTickness;
	double finSpacing;
	double k = 180.0;
	fin(double xFPI,double xfinTickness);
	~fin() {};
	void geoCalculate();

};

