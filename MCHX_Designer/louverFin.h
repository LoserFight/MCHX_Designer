#pragma once
#include "fin.h"
class louverFin :
    public fin
{
public:
    double theta;//degree
    double Lh;
    double Fd;
    double Ll;
    double Lp;
    double Fl;
    //double sigf;//finTickness
    double Fp;//same with finspacing
    louverFin(double xFPI, double xfinTickness,double xtheta,double xLh,double xFd,double xLl,double xLp,double xFl);
};

