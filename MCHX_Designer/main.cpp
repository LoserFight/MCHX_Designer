#include"rei.h"
#include <stdlib.h>
#include <stdio.h>
#include<vector>
#include"MCHE.h"
#include"InputPara.h"
#include<regex>
#include "matrix.h"
#include "mclcppclass.h"
#include <mclmcrrt.h>
EXTERN_C bool mclInitializeApplication(const char** options, size_t count);
InPut_data sed;
int main(int argc,char* argv[])
{
    if (!mclInitializeApplication(NULL, 0))
    {
        cout << "Could not initialize the application.\n";
        exit(0);
    }

    if (!NewtonToSolveGdpInitialize())
    {

        cout << "Could not initialize NewtonToSolveGdp!" << endl;

        exit(0);

    }
    //refri test("air.ppf");
    //test.tem = 300.0;
    //test.pressure = 101325;
    ////test.quality = 0.1;
    //test.tpProperty();
    //cout << test.quality << endl;

    //vector<UINT> pas = { 69,26 };
    //refri re("r134a.fld");
    //re.tem = 353.15;
    //re.pressure = 1317905.491;
    //re.tpProperty();
    //re.qm = 0.25;
    //air_con a_in;
    //a_in.P = 101325;
    //a_in.T = 303.15;
    //a_in.v_in = 2.0;
    //a_in.w = 0.4;
    //for (int i = 0; i < argc; i++) {
    //    std::cout << argv[i] << endl;
    //}

    sed.air_in.P = atof(argv[1]);
    sed.air_in.T = atof(argv[2]);
    sed.air_in.v_in = atof(argv[3]);

    sed.t.tube_length = atof(argv[4]);
    sed.t.tube_width = atof(argv[5]);
    sed.t.tube_ports = atoi(argv[6]);
    sed.t.tube_thickness = atof(argv[7]);
    sed.t.vertical_space = atof(argv[8]);
    sed.t.port_height = atof(argv[9]);
    sed.t.port_width = atof(argv[10]);
    sed.t.geoCalculate();

    sed.f.FPI = atof(argv[11]);

    sed.f.Fd = atof(argv[12]);
    sed.f.finTickness = atof(argv[13]);
    sed.f.Ll = atof(argv[14]);
    sed.f.Lp = atof(argv[15]);
    sed.f.theta = atof(argv[16]);
    sed.f.Lh = sed.f.Lp * std::tan(sed.f.theta / 180.0 * 3.14) / 2.0;
    sed.f.Fl = sed.t.vertical_space - sed.t.tube_thickness;
    sed.f.geoCalculate();
    string reName(argv[17]);
    reName += ".fld";
    sed.re = refri(reName.c_str());

    sed.re.tem = atof(argv[18]);
    sed.re.pressure = atof(argv[19]);
    sed.re.qm = atof(argv[20]);
    sed.re.tpProperty();
    sed.segms = atoi(argv[21]);

    string passes(argv[22]);

    std::regex patten(",");
    std::sregex_token_iterator first{ passes.begin(),passes.end(), patten, -1 }, last;
    vector<string> sPa = { first,last };
    for (auto strr : sPa) {
        sed.pas.push_back(stoi(strr));
        sed.tube_s += stoi(strr);
    }




    MCHE* HE;
    HE=new MCHE(sed.pas, sed.segms, sed.re, sed.air_in, sed.t, sed.f);

    HE->geometryDesign();
    HE->bbk = 1;
    if (!HE->solve(2)) {
        cout << "The input mass flow value is too large." << endl;
        return -1;
    }
    HE->savePressure();
    HE->saveResult();
    cout<<"Heat load:"<<HE->Q_total<<endl;
    cout << "Refri pressure drop:" << HE->DP_re << endl;
    cout << "Refri Charge:  L:" << HE->reChargeL << endl;

    cout << "Refri Charge:  VL:" << HE->reChargeVL << endl;

    cout << "Refri Charge:  V:" << HE->reChargeV << endl;
    cout << "outlet refri deltaSatT: " << HE->DelSat << endl;
    mclTerminateApplication();
    NewtonToSolveGdpTerminate();
    return 0;
}