#ifndef INPUTPARA_H
#define INPUTPARA_H
#include"MCHE.h"
#include"rei.h"

struct InPut_data{
    air_con air_in;
    int segms;
    vector<UINT> pas;
    refri re;
    microTube t;
    louverFin f;
    int tube_s;
    InPut_data():re("r134a.fld"),t(2.323,
        0.0013,
        0.016,
        0.0096,
        0.0,
        0.0009, 0.001, 15),
        f(18.0, 0.00005, 28.99, 0.001 * std::tan(28.99 / 180.0 * 3.14)/2.0,t.tube_width, 0.008, 0.001,t.vertical_space-t.tube_thickness){
    }
};
extern InPut_data sed;


#endif // INPUTPARA_H
