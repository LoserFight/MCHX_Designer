#include "louverFin.h"



louverFin::louverFin(double xFPI, double xfinTickness, double xtheta, double xLh, double xFd, double xLl, double xLp, double xFl):fin(xFPI, xfinTickness)
{
	Fp = finSpacing;
    theta=xtheta;//degree
    Lh=xLh;
    Fd=xFd;
    Ll=xLl;
    Lp=xLp;
    Fl= xFl;

}
