#pragma once

#include <iostream>
#include <iomanip>
#include <math.h>
#include <windows.h>
#include <string>

using namespace std;

//const char refrigerant[] = "r22.fld";     //制冷剂

/*----------------------------------------------------结构体：制冷剂、空气、水----------------------------------------------------*/
/*制冷剂性质求解结构体*/
typedef struct refri
{
	
	/* REFPROP函数单位 t-K, p-kPa, d-mol/L, h-J/mol, s-J/(mol.K) vis- uPa.s cp-J/(mol.K) th_con-W/(m.K) */
	double tem, pressure, p_cri, molmass;                  //温度 压力 临界压力 摩尔质量
	double quality, density, densityLiq, densityGas;       //干度 密度
	double viscosity, viscosityLiq, viscosityGas;		   //动力粘度
	double enthalpy, entropy, k_com;					   //焓 熵 绝热指数
	double HeatCapacity_p, HeatCapacity_v, thermalConduct;  //定压/定容比热容 热导率
	double Tension;											//表面张力
	double qm, heatFlux;                                   //流量 热流量

	/*临时参数（REFPROP函数使用）*/
	typedef void(__stdcall* fp_SETUPdllTYPE)(long&, char*, char*, char*, long&, char*, long, long, long, long); fp_SETUPdllTYPE SETUPdll;
	typedef void(__stdcall* fp_SETMIXdllTYPE)(char*, char*, char*, long&, char*, double*, long&, char*, long, long, long, long, long); fp_SETMIXdllTYPE SETMIXdll;
	typedef void(__stdcall* fp_INFOdllTYPE)(long&, double&, double&, double&, double&, double&, double&, double&, double&, double&, double&); fp_INFOdllTYPE INFOdll;
	typedef void(__stdcall* fp_PHFLSHdllTYPE)(double&, double&, double*, double&, double&, double&, double&, double*, double*, double&, double&, double&, double&, double&, double&, long&, char*, long); fp_PHFLSHdllTYPE PHFLSHdll;
	typedef void(__stdcall* fp_PSFLSHdllTYPE)(double&, double&, double*, double&, double&, double&, double&, double*, double*, double&, double&, double&, double&, double&, double&, long&, char*, long); fp_PSFLSHdllTYPE PSFLSHdll;
	typedef void(__stdcall* fp_TPFLSHdllTYPE)(double&, double&, double*, double&, double&, double&, double*, double*, double&, double&, double&, double&, double&, double&, double&, long&, char*, long); fp_TPFLSHdllTYPE TPFLSHdll;
	
	typedef void(__stdcall* fp_PQFLSHdllTYPE)(double&, double&, double*, double&, double&, double&, double&, double*, double*, double&, double&, double&, double&, double&, double&, long&, char*, long); fp_PQFLSHdllTYPE PQFLSHdll;
	//typedef void(__stdcall* fp_TQFLSHdllTYPE)(double&, double&, double*, double&, double&, double&, double*, double*, double&, double&, double&, double&, double&, double&, double&, long&, char*, long); fp_TQFLSHdllTYPE TQFLSHdll;

	
	typedef void(__stdcall* fp_SATPdllTYPE)(double&, double*, long&, double&, double&, double&, double*, double*, long&, char*, long); fp_SATPdllTYPE SATPdll;
	typedef void(__stdcall* fp_SATTdllTYPE)(double&, double*, long&, double&, double&, double&, double*, double*, long&, char*, long); fp_SATTdllTYPE SATTdll;
	typedef void(__stdcall* fp_TRNPRPdllTYPE)(double&, double&, double*, double&, double&, long&, char*, long); fp_TRNPRPdllTYPE TRNPRPdll;
	typedef void(__stdcall* fp_SURFTdllTYPE)(double&, double&, double*, double&, long&, char*, long); fp_SURFTdllTYPE SURFTdll;

	double x[20], xliq[20], xvap[20], f[20]; long ierr;
	char hf[255 * 20], hrf[3 + 1], herr[255 + 1], hfmix[255 + 1];
	char refrigerant[255 * 20];
	double t, p, dl, dv, d, q, e, h, s, cv, cp, w, b, c, eta, tcx;
	
	/*构造函数（初始化REFPROP模块）*/
	refri(const char * a)
	{
		double wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas;
		long info_index = 1; long ii = 1;
		HINSTANCE RefpropdllInstance = LoadLibrary("REFPRP64.DLL"); //加载DLL

		/*可能用到的REFPROP函数*/
		SETUPdll = (fp_SETUPdllTYPE)GetProcAddress(RefpropdllInstance, "SETUPdll");
		SETMIXdll = (fp_SETMIXdllTYPE)GetProcAddress(RefpropdllInstance, "SETMIXdll");
		INFOdll = (fp_INFOdllTYPE)GetProcAddress(RefpropdllInstance, "INFOdll");
		PHFLSHdll = (fp_PHFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "PHFLSHdll");
		
		PSFLSHdll = (fp_PSFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "PSFLSHdll");
		TPFLSHdll = (fp_TPFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "TPFLSHdll");

		//TQFLSHdll = (fp_TQFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "TQFLSHdll");
		//PQFLSHdll = (fp_PQFLSHdllTYPE)GetProcAddress(RefpropdllInstance, "PQFLSHdll");
		SATPdll = (fp_SATPdllTYPE)GetProcAddress(RefpropdllInstance, "SATPdll");
		SATTdll = (fp_SATTdllTYPE)GetProcAddress(RefpropdllInstance, "SATTdll");
		TRNPRPdll = (fp_TRNPRPdllTYPE)GetProcAddress(RefpropdllInstance, "TRNPRPdll");

		SURFTdll = (fp_SURFTdllTYPE)GetProcAddress(RefpropdllInstance, "SURFTdll");

		strcpy_s(refrigerant, a);
		strcpy_s(hf, refrigerant);
		strcpy_s(hfmix, "hmx.bnc");
		strcpy_s(hrf, "DEF");
		strcpy_s(herr, "Ok");
		x[0] = 1.0;
		SETUPdll(ii, hf, hfmix, hrf, ierr, herr, 255 * 20, 255, 3, 255);
		INFOdll(info_index, wm, ttp, tnbp, tc, pc, dc, zc, acf, dip, rgas);
		if (ierr != 0) printf("%s\n", herr);
		molmass = wm; p_cri = 1000 * pc;
	}

	/*初始化setup制冷剂求解REFPROP */
	void setuprefri()
	{
		long ii = 1;
		strcpy_s(hf, refrigerant);
		strcpy_s(hfmix, "hmx.bnc");
		strcpy_s(hrf, "DEF");
		strcpy_s(herr, "Ok");
		SETUPdll(ii, hf, hfmix, hrf, ierr, herr, 255 * 20, 255, 3, 255);
	}

	/*根据 p 饱和态求解（i=1 液 i=2 气）*/
	void pSatProperty(long i)
	{
		setuprefri();
		p = 0.001 * pressure;
		SATPdll(p, x, i, t, dl, dv, xliq, xvap, ierr, herr, 255);
		if (i == 1) t -= 0.001; else if (i == 2) t += 0.001;
		TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		if (i == 1) density = dl * molmass; else if (i == 2) density = dv * molmass;
		tem = t; pressure = 1000 * p;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta;
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass;
		k_com = cp / cv; thermalConduct = tcx;
	}

	/*根据 T 饱和态求解（i=1 液 i=2 气）*/
	void tSatProperty(long i)
	{
		setuprefri();
		t = tem;
		SATTdll(t, x, i, p, dl, dv, xliq, xvap, ierr, herr, 255);
		if (i == 1) t -= 0.001; else if (i == 2) t += 0.001;
		TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		if (i == 1) density = dl * molmass; else if (i == 2) density = dv * molmass;
		tem = t; pressure = 1000 * p;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta;
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass;
		k_com = cp / cv; thermalConduct = tcx;
	}

	/*根据 T p 求解状态参数*/
	void tpProperty()
	{
		setuprefri();
		t = tem; p = 0.001 * pressure;
		TPFLSHdll(t, p, x, d, dl, dv, xliq, xvap, q, e, h, s, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		tem = t; pressure = 1000 * p; density = d * molmass; quality = q;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta;
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass;
		k_com = cp / cv; thermalConduct = tcx;
	}

	/*根据 p h 求解状态参数*/
	void phProperty()
	{
		setuprefri();
		p = 0.001 * pressure; h = molmass * enthalpy;
		PHFLSHdll(p, h, x, t, d, dl, dv, xliq, xvap, q, e, s, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		tem = t; pressure = 1000 * p; quality = q;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta;
		density = d * molmass; densityLiq = dl * molmass; densityGas = dv * molmass;
		//if (quality > 0) density = densityGas * quality + densityLiq * (1 - quality);
		if (quality > 0) density = 1.0 / (quality / densityGas + (1.0 - quality) / densityLiq);
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass; k_com = cp / cv;
	}

	/*根据 p s 求解状态参数*/
	void psProperty()
	{
		setuprefri();
		p = 0.001 * pressure; s = molmass * entropy;
		PSFLSHdll(p, s, x, t, d, dl, dv, xliq, xvap, q, e, h, cv, cp, w, ierr, herr, 255);
		TRNPRPdll(t, d, x, eta, tcx, ierr, herr, 255);
		tem = t; pressure = 1000 * p; quality = q;
		enthalpy = h / molmass; entropy = s / molmass; viscosity = 1e-6 * eta;
		density = d * molmass; densityLiq = dl * molmass; densityGas = dv * molmass;
		//if (quality > 0) density = densityGas * quality + densityLiq * (1 - quality);
		if (quality > 0) density = 1.0 / (quality / densityGas + (1.0 - quality) / densityLiq);
		HeatCapacity_p = cp / molmass; HeatCapacity_v = cv / molmass; k_com = cp / cv;
	}

	/*根据 p q 求解状态参数*/
	void pqProperty()
	{

		double tmp_q= quality;
		pSatProperty(1L);
		double h_l = enthalpy;
		pSatProperty(2L);
		double h_v = enthalpy;
		enthalpy = tmp_q * h_v + (1.0 - tmp_q) * h_l;
		phProperty();
	}

	//求液体表面张力，在获得其他参数后再使用
	void surTension() {
		SURFTdll(t, dl, x, Tension, ierr, herr, 255);

	}


}refri;