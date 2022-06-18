#include "microTube.h"

microTube::microTube(double length, double t, double w, double p, double covered, double Pw, double Ph, int ports)
{
	tube_length = length;
	tube_thickness = t;
	tube_width=w;
	vertical_space=p;
	covered_lengh = covered;

	tube_ports=ports;
	port_width=Pw;
	port_height=Ph;


	geoCalculate();
}

void microTube::geoCalculate()
{
	port_Ac = tube_ports * port_height * port_width;
	Dh = 2.0 * port_height * port_width/(port_height+ port_width);
	P = (port_height + port_width) * 2.0;
}

void microTube::geoCalculateCircle()
{
	port_Ac = tube_ports * 3.1415 * port_width * port_width/4.0;
	Dh = port_width;
	P = 3.1415 * Dh;
}

void microTube::specific(double cDh, double Ac)
{
	port_Ac = tube_ports * Ac;
	Dh = cDh;
	P = 4.0 * Ac / Dh;
}
