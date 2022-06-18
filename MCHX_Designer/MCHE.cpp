#include "MCHE.h"

//EXTERN_C bool mclInitializeApplication(const char** options, size_t count);
////#pragma comment(lib,"NewtonToSolveGdp.lib")

bool MCHE::Gnielinski(Segment& a,const double G,double &h,double &f)
{
	double Nu_hot;
	double f_hot;
	double Dh = m_tube->Dh;

	//refrigerant_tmp_out.pressure = a.P_rei_out;
	//refrigerant_tmp_out.enthalpy = a.h_rei_out;
	//refrigerant_tmp_out.phProperty();


	//double t_hot_m = (a.T_hot_in + a.T_hot_out) / 2.0;
	double rou_hot = (refrigerant_tmp_in.density+refrigerant_tmp_out.density)/2.0;
	double u_hot = (refrigerant_tmp_in.viscosity + refrigerant_tmp_out.viscosity) / 2.0;
	double lambda_hot = (refrigerant_tmp_in.thermalConduct + refrigerant_tmp_out.thermalConduct) / 2.0;
	double cp_hot = (refrigerant_tmp_in.HeatCapacity_p ) *1000;
	//double Pr_hot = u_hot * cp_hot * 1000 / lambda_hot;
	double Pr_hot = u_hot * cp_hot  / lambda_hot;
	double Re_t;
	Re_t = G * Dh / u_hot;
	f_hot = Re_t <= 2300.0 ? 15.767 / Re_t : 0.25 * pow((1.0 / (1.8 * log10(Re_t) - 1.5)), 2.0);
	if (Re_t <= 2300.0) {
		Nu_hot = 4.089;

	}
	else if (Re_t < 5000.0) {

		Nu_hot = f_hot / 2.0 * (5000.0 - 1000.0) * Pr_hot / (1.0 + 12.7 * (pow(Pr_hot, 2.0 / 3.0) - 1) * sqrt(f_hot / 2.0));
		Nu_hot = linearInd(4.089, Nu_hot, 2300, 5000, Re_t);
	}
	else if (Re_t < 5e6) {

		Nu_hot = f_hot / 2.0 * (Re_t - 1000) * Pr_hot / (1.0 + 12.7 * (pow(Pr_hot, 2.0 / 3.0) - 1) * sqrt(f_hot / 2.0));
	}
	else {
		std::cout << "Re int hot channel may be too large:\n Re_hot:" << Re_t << std::endl;
		return false;
	}
	bool isMicro = true;
	if (isMicro) {
		Nu_hot = Nu_hot * (1 + 7.6E-5 * (1.0 - pow(m_tube->Dh / 0.001164, 2.0)));
	}


	h = Nu_hot * lambda_hot / Dh;
	f = f_hot;
	//a.K_i = 1.0 / (1.0 / a.h_hot + 1.0 / a.h_cold + t / lambda_316L);

	return true;
}

bool MCHE::Friedel(Segment& a,double G_r,double &dp,double &Re_LL)
{
	double G = G_r;//
	refri vapor(refrigerant_tmp_in.hf);
	vapor.tem = (refrigerant_tmp_in.tem + refrigerant_tmp_in.tem) / 2.0;//
	vapor.tSatProperty(2L);
	refri Liq(refrigerant_tmp_in.hf);
	Liq.tem = (refrigerant_tmp_in.tem + refrigerant_tmp_in.tem) / 2.0;
	Liq.tSatProperty(1L);

	double u_g = (vapor.viscosity + vapor.viscosity) / 2.0;
	double u_l= (Liq.viscosity+ Liq.viscosity) / 2.0;
	double rou_g = (vapor.density + vapor.density) / 2.0;
	double rou_l = (Liq.density + Liq.density) / 2.0;

	refrigerant_tmp_in.surTension();
	double sig = refrigerant_tmp_in.Tension;
	double Re_G = G * m_tube->Dh / u_g;
	double Re_L = G * m_tube->Dh / u_l;
	double f_GO, f_LO;
	Re_LL = Re_L;
	//BlasiusOnly(Re_G, f_GO);
	//BlasiusOnly(Re_L, f_LO);
	Blasius(Re_G, f_GO);
	Blasius(Re_L, f_LO);
	double x = (refrigerant_tmp_in.quality + refrigerant_tmp_in.quality) / 2;
	double rou_ave = 1.0 / (x / rou_g + (1.0 - x) / rou_l);
	//double alpha = pow(1.0 + a_B * pow((1 - x) / x, a_n1) * pow(u_l / u_g, a_n2) * pow(rou_g / rou_l, a_n3), -1.0);
	double A1 = pow(1.0 - x, 2.0) + x * x * (rou_l / rou_g * f_GO / f_LO);
	double A2 = pow(x, 0.78) * pow(1.0 - x, 0.224);
	double A3 = pow(rou_l / rou_g, 0.91) * pow(u_g / u_l, 0.19) * pow(1.0 - u_g / u_l, 0.7);
	double Fr = G / 9.8 * G / rou_ave / rou_ave / m_tube->Dh;
	double We=G/rou_ave*G* m_tube->Dh/sig;


	double faiLO = A1+3.24*A2*A3/pow(Fr,0.045)/pow(We,0.035);
	double deltaP_LO=4.0*a.L/ m_tube->Dh* f_LO*G*G/2.0/rou_l;
	double deltaP_f = faiLO * deltaP_LO;
	dp = deltaP_f;
	return false;
}

bool MCHE::Blasius(double Re, double& f)
{
	if(Re>2300.0)
		f = 0.0791 * pow(Re, -0.25);//fanning
	else {
		f = 16 / Re;
	}
	return false;
}


bool MCHE::BlasiusOnly(double Re, double& f)
{

	f = 0.0791 * pow(Re, -0.25);//fanning

	return false;
}

bool MCHE::condShah2016(Segment& a,double G_t,double & hi)
{
	refri vapor(refrigerant_tmp_in.hf);
	vapor.tem =( refrigerant_tmp_in.tem+ refrigerant_tmp_out.tem)/2.0;
	vapor.tSatProperty(2L);
	refri Liq(refrigerant_tmp_in.hf);
	Liq.tem = (refrigerant_tmp_in.tem + refrigerant_tmp_out.tem) / 2.0;
	Liq.tSatProperty(1L);

	double G = G_t;//
	double u_g = (vapor.viscosity + vapor.viscosity) / 2;
	double u_l = (Liq.viscosity + Liq.viscosity) / 2;
	double rou_g = vapor.density;
	double rou_l = Liq.density;

	//double x = (refrigerant_tmp_in.quality + refrigerant_tmp_out.quality) / 2;
	double x = (refrigerant_tmp_in.quality + refrigerant_tmp_in.quality) / 2;
	refrigerant_tmp_in.surTension();

	double u = (refrigerant_tmp_in.viscosity + refrigerant_tmp_out.viscosity) / 2.0;
	//double lambda = (refrigerant_tmp_in.thermalConduct + refrigerant_tmp_out.thermalConduct) / 2.0;
	//double cp = (refrigerant_tmp_in.HeatCapacity_p + refrigerant_tmp_out.HeatCapacity_p) / 2.0;
	double sig = refrigerant_tmp_in.Tension;
	double Z = pow(1.0 / x - 1.0, 0.8) * (refrigerant_tmp_in.pressure- refrigerant_tmp_out.pressure);
	double Jg = x * G / (9.8 * m_tube->Dh * rou_g * (rou_l - rou_g));
	double We_GT = G * G * m_tube->Dh / rou_g / sig;


	double Pr_L = u_l * Liq.HeatCapacity_p / Liq.thermalConduct*1000.0;

	double Re_LT = G * m_tube->Dh / u_l;
	double Re_GT = G * m_tube->Dh / u_g;
	double Re_LO = Re_LT * (1 - x);

	double h_tp;// two - phase heat transfer coefficient[W m−2 C−1]
	if (We_GT > 100 && Jg > 0.98 * pow(Z + 0.263, -0.62)) {
		double h_LT = 0.023 * pow(Re_LT, 0.8) * pow(Pr_L, 0.4) * Liq.thermalConduct / m_tube->Dh;
		h_tp = h_LT * (1.0 + 1.128 * pow(x, 0.817) * pow(Liq.density / vapor.density, 0.3685)
			* pow(Liq.viscosity / vapor.viscosity, 0.2363)
			* pow(1.0 - vapor.viscosity / Liq.viscosity, 2.144) * pow(Pr_L, -0.1));

	}
	else if (Jg<0.95 * pow(1.254 + 2.27 * pow(Z, 1.249), -1.0)) {
		h_tp = 1.32 * pow(Re_LO, -0.3333) * 
			pow(Liq.density / Liq.viscosity * (Liq.density - vapor.density)
			/ Liq.viscosity * 9.8
			* pow(Liq.thermalConduct, 3), -0.3333);
	}
	else {
		double h_LT = 0.023 * pow(Re_LT, 0.8) * pow(Pr_L, 0.4) * Liq.thermalConduct / m_tube->Dh;
	double h1= h_LT * (1.0 + 1.128 * pow(x, 0.817) * pow(Liq.density / vapor.density, 0.3685)
		* pow(Liq.viscosity / vapor.viscosity, 0.2363)
		* pow(1.0 - vapor.viscosity / Liq.viscosity, 2.144) * pow(Pr_L, -0.1));
	double h2= 1.32 * pow(Re_LO, -0.3333) *
		pow(Liq.density / Liq.viscosity * (Liq.density - vapor.density)
			/ Liq.viscosity * 9.8
			* pow(Liq.thermalConduct, 3), -0.3333);
	h_tp = h1 + h2;


	
	}
	hi = h_tp;

	return false;
}

bool MCHE::WangChang_louver(double &f,double &ho)
{
	double Vc = air_in.v_in * Afr_air / Ac_air;

	//double Vc =1.985;

	double Re_Lp=1.165* Vc *m_fin->Lp/18.6*1E6;

	double Dh = Dh_air;
	double f1, f2, f3;
	if(Re_Lp>5000.0) return false;
	double j1 = pow(Re_Lp, -0.49) * pow(m_fin->theta / 90.0, 0.27) * pow(m_fin->Fp / m_fin->Lp, -0.14);
	double j2 = pow(m_fin->Fl / m_fin->Lp, -0.29) * pow(m_tube->tube_width / m_fin->Lp, -0.23);
	double j3= pow(m_fin->Ll / m_fin->Lp, 0.68) * pow(m_tube->vertical_space / m_fin->Lp, -0.28) * pow(m_fin->finTickness / m_fin->Lp, -0.05);
	double j = j1 * j2 * j3;
	ho = j * 1.165 * Vc * 1005 / pow(0.70, 2.0 / 3.0);

	if (Re_Lp < 150.0) {
		 f1 = 14.39 * pow(Re_Lp, -0.805 * m_fin->Fp / m_fin->Fl) * pow(log(1.0 + (m_fin->Fp / m_fin->Lp)), 3.04);
		 f2 = pow(log(0.9 + pow(m_fin->finTickness / m_fin->Fp, 0.48)), -1.435) * pow(Dh / m_fin->Lp, -3.01) * pow(log(0.5 * Re_Lp), -3.01);
		 f3 = pow(m_fin->Fp / m_fin->Ll, -0.308) * pow(m_fin->Fd / m_fin->Ll, -0.308) * pow(m_fin->theta, 0.35) * exp(-0.1167 * m_tube->vertical_space / m_tube->tube_thickness);
	}
	else {
		f1 = 4.97 * pow(Re_Lp, 0.6049 - 1.064 / pow(m_fin->theta, 0.2));
		double f12=	 pow(log(0.9 + pow(m_fin->finTickness / m_fin->Fp,0.5)), -0.527);
		f1 = f1 * f12;
		f2 = pow(log(0.3 * Re_Lp) * (Dh / m_fin->Lp), -2.966) * pow(m_fin->Fp / m_fin->Ll, -0.7931*(m_tube->vertical_space/(m_tube->vertical_space-m_tube->tube_thickness)));//
		f3 = pow(m_tube->vertical_space / m_tube->tube_thickness, -0.0446) * pow(log(1.2 + pow(m_fin->Lp / m_fin->Fp, 1.4)  ), -3.553) * pow(m_fin->theta, -0.477);
		//f3 = pow(m_tube->vertical_space / m_tube->tube_thickness, -0.0446) * (-3.553)*log(1.2 + pow(m_fin->Lp / m_fin->Fp, 1.4)) * pow(m_fin->theta, -0.477);

	}
	f = f1 * f2 * f3;
	return true;
}

double MCHE::WangChang_fin_eta(double h0)
{
	double m = sqrt(2.0 * h0 / m_fin->k / m_fin->finTickness);
	double XL = m_tube->vertical_space/2.0;
	double XM = m_tube->tube_width;
	double r = m_tube->tube_thickness;
	double rr = 1.28 * XM / r * pow(XL / XM - 0.2, 0.5);
	double fai = (rr - 1) * (1 + 0.35 * log(rr));
	double mrfai = m * r * fai;
	double nf = tanh(mrfai) / mrfai;
	return 1 - Af / (Af + At) * (1 - nf);
}

double MCHE::classic(double h0)
{
	double l1 = m_fin->Fl / 2.0 - m_fin->finTickness;
	double m_i = sqrt(2.0*h0/m_fin->k/ m_fin->finTickness*(1.0+ m_fin->finTickness/m_fin->Fl));
	double ml = m_i * l1;
	double nf = tanh(ml) / ml;


	return 1 - Af / (Af + At) * (1 - nf);
}

bool MCHE::jacobi_louver(double& f, double& ho,const refri & a_in, const refri& a_out)
{
	double c1 = 0.872;
	double c2 = 0.219;
	double c3 = -0.0881;
	double c4 = 0.149;
	double c5 = -0.259;
	double c6 = 0.54;
	double c7 = -0.902;
	double c8 = 2.62;
	double c9 = 0.301;
	double c10 = -0.458;
	double c11 = -0.00874;
	double c12 = 0.049;
	double c13 = 0.142;
	double c14 = -0.0065;

	double d1 = 3.69;
	double d2 = -0.256;
	double d3 = 0.904;
	double d4 = 0.2;
	double d5 = 0.733;
	double d6 = 0.648;
	double d7 =- 0.647;
	double d8 = 0.799;
	double d9 = -0.845;
	double d10 = 0.0013;
	double d11 = 0.126;

	double N_LB = 2;

	double Vc = air_in.v_in * Afr_air / Ac_air;

	//double Vc =1.985;

	double Re_Lp = 1.165 * Vc * m_fin->Lp / 18.6 * 1E6;
	double Fp = m_fin->Fp;
	double Fl = m_fin->Fl;
	double Ll = m_fin->Ll;
	double Lp = m_fin->Lp;
	double Fd = m_fin->Fd;
	double f_re = pow(Re_Lp * Fp / Lp, d9) + d10 * pow(Re_Lp, d11 * (m_fin->finTickness / Fp));
	double a = m_fin->theta / 180.0 * 3.1416;
	double f_cor = d1 * f_re * pow(N_LB, d2) * pow(Fp / Lp, d3) * sin(a + d4) * pow(1.0 - Fl / m_tube->vertical_space,d5)*pow(Ll/Fl,d6)*pow(m_fin->finTickness/Lp,d7)
					*pow(Fl/Fp,d8);

	f = f_cor;

	double j_re = pow(Re_Lp, c10 + c11 * cosh(Fp / Lp - 1.0));
	double j_low = 1.0 - sin(Lp / Fp * a) / (cosh(c12 * Re_Lp - c13 * Fd / N_LB / Fp));
	double j_louver = 1 - c14 * tan(a) * (Fd / N_LB / Fp) * cos(2 * 3.1416 * (Fp / Lp / tan(a) - 1.8));
	double j_cor = c1 * j_re * j_low * j_louver * pow(a, c2) * pow(N_LB, c3) * pow(Fl / Lp, c4) *
		pow(Fd / Fp, c5) * pow(Ll / Fl, c6) * pow(Fl / m_tube->vertical_space, c7) * pow(1.0 - m_fin->finTickness / Lp, c8)
		* pow(Lp / Fp, c9);

	double G = Vc * (a_in.density + a_out.density) / 2.0;
	double Cp = (a_in.HeatCapacity_p + a_out.HeatCapacity_p)*1000 / 2.0;
	double Pr_in = a_in.viscosity * a_in.HeatCapacity_p * 1000 / a_in.thermalConduct;
	double Pr_out= a_out.viscosity * a_out.HeatCapacity_p * 1000 / a_out.thermalConduct;
	double Pr = (Pr_in + Pr_out) / 2.0;

	ho = j_cor * G * Cp * pow(Pr, -2.0 / 3.0);

	return false;
}

double MCHE::KeCal(double sig)
{

	return -0.03354*pow(sig,3)+1.013*pow(sig,2)-1.949*sig+0.9702;
}

double MCHE::KcCal(double sig)
{
	return -0.14*pow(sig,3)-0.24*pow(sig,2)-0.047*sig+0.42;
}

double MCHE::alpha(const refri& r)
{
	double x = r.quality;
	if (x <= 0.0) return 0.001;
	refri vapor(r.hf);
	vapor.tem = r.tem ;
	vapor.tSatProperty(2L);
	refri Liq(r.hf);
	Liq.tem = r.tem;
	Liq.tSatProperty(1L);
	


	return pow(1.0+ a_B*pow((1.0-x)/x,a_n1)*pow(Liq.viscosity/vapor.viscosity,a_n2)*pow(vapor.density/Liq.density,a_n3),-1.0);
}

double MCHE::AirLessNTU(double Cx, double NTU)
{
	if (Cx < 0.001) {
		return exp(Cx * (1.0 - exp(-1.0 * NTU))) * (1.0 - exp(-1.0 * NTU));
	}
	else if(Cx>0.0) {
		return 1.0 / Cx * (1.0 - exp(-Cx * (1.0 - exp(-NTU))));

	}
	return 0.0;
}

double MCHE::NTsolveG(double R_T, double R_L,double R_a, double n1, double n2,double n3, double Dp,double G0)
{
	const double eps = 1E-10;
	double x = G0;
	while (true) {
		double f = R_T * pow(x, n1) + R_L * pow(x, n2)+ R_a*pow(x,n3) - Dp;
		double df = n1 * R_T * pow(x, n1 - 1.0) + n2*R_L * pow(x, n2 - 1.0)+ n3 * R_a * pow(x, n3 - 1.0);
		double nx = x - f / df;
		if (abs(nx - x) < eps) break;
		x = nx;
	}
	return x;
}

double MCHE::TubeLessNTU(double Cx, double NTU)
{
	return 1.0-exp((exp(-NTU*Cx)-1.0)/Cx);
}

void MCHE::savePressure()
{
	ofstream out("pressure.xls");
	ofstream out2("quality.xls");
	ofstream out3("Re.xls");
	ofstream out4("airDrop.xls");
	ofstream out5("Rs.xls");
	for (int i = 0; i <= m_segment; i++) {
		for (int j = 0; j < sum_tubes*blockeds; j++) {
			out << m_HE[j * (m_segment + 1) + i].P_rei_out << (j < sum_tubes* blockeds - 1 ? '\t' : '\n');
			double tmp = ((m_HE[j * (m_segment + 1) + i].quality > -0.00) && (m_HE[j * (m_segment + 1) + i].quality <= 1.00)) ? m_HE[j * (m_segment + 1) + i].quality
				: ((m_HE[j * (m_segment + 1) + i].quality < -0.00) ? -1.1 : 1.1);
			out2<< tmp<< (j < sum_tubes* blockeds - 1 ? '\t' : '\n');
			out3<< m_HE[j * (m_segment + 1) + i].Re<< (j < sum_tubes* blockeds - 1 ? '\t' : '\n');
			if((j+1)%blockeds==0) out4 << air_in.P-m_HE[j * (m_segment + 1) + i].P_air_out << (j < sum_tubes* blockeds - 1 ? '\t' : '\n');
			out5<< m_HE[j * (m_segment + 1) + i].Rs_La << (j < sum_tubes* blockeds - 1 ? '\t' : '\n');
		}
	}
	out.close();
	out2.close();
	out3.close();
	out4.close();
	out5.close();
	

}

void MCHE::saveResult()
{
	ofstream out("Result.txt");
	out << "HeatLoad/W\t" << Q_total<<std::endl;
	out << "RefrigerantPressureDrop/Pa\t" << DP_re << std::endl;
	out << "RefrigerantCharge/kg\t" << reChargeL+ reChargeVL + reChargeV<< std::endl;
	out << "DeltaSatT/K\t" << DelSat << std::endl;
	double Dp_air=0.0;
	double ht_air_ave=0.0;

		for (int jj = 0; jj < sum_tubes * blockeds; jj++) {
			for (int i = 0; i <= m_segment; i++) {
				if ((jj + 1) % blockeds == 0&&i!=0) Dp_air += air_in.P - m_HE[jj * (m_segment + 1) + i].P_air_out;
				ht_air_ave += m_HE[jj * (m_segment + 1) + i].ht_air;
			}
		}
	Dp_air /= double(sum_tubes) * m_segment;
	ht_air_ave /= double(sum_tubes) * blockeds * m_segment;

	out << "AverageDeltaP_air/Pa\t" << Dp_air << std::endl;
	out << "AverageHTC_air/(W/m2s)\t" << ht_air_ave << std::endl;
	out << "Volume\t" << HE_volume_tube + HE_volume_fin << std::endl;

	out.close();
}

MCHE::MCHE(vector<UINT>& p, int segNum, const refri& refiIn,const air_con& airin):
	refrigerant_in(refiIn), refrigerant_out(refiIn), refrigerant_tmp_in(refiIn), refrigerant_tmp_out(refiIn),
	air_in(airin)
{
	m_segment = segNum; 
	passes = std::move(p);
	sum_tubes = 0;
	for (size_t i = 0; i < passes.size(); i++) {
		sum_tubes += passes[i];
	}

	m_HE = new Segment[long(sum_tubes) * (m_segment + 1L)]; //(refrigerant_in.refrigerant);//microTube(double length,double t,double w,double p,double covered,double Pw,double Ph,int ports);	
	m_tube = new microTube(2.323,
		0.0013,
		0.016,
		0.0096,
		0.0,
		0.0009, 0.001, 15);
	
	double Lvtheta = 28.9974;
	//louverFin(double xFPI, double xfinTickness,double xtheta,double xLh,double xFd,double xLl,double xLp)
	m_fin = new louverFin(18.0, 0.00005, Lvtheta, 0.001 * std::tan(Lvtheta / 180.0 * 3.14)/2.0,
		m_tube->tube_width, 0.008, 0.001,m_tube->vertical_space-m_tube->tube_thickness);


	//geometryDesign();
	q_re = refrigerant_in.qm;
	reChargeVL=reChargeL = reChargeV = 0.0;

}

MCHE::MCHE( vector<UINT>& p, int segNum, const refri& refiIn, const air_con& airin, const microTube& tdata, const louverFin& fdata) :
	refrigerant_in(refiIn), refrigerant_out(refiIn), refrigerant_tmp_in(refiIn), refrigerant_tmp_out(refiIn),
	air_in(airin)
{
	m_segment = segNum;
	passes = std::move(p);
	sum_tubes = 0;
	for (size_t i = 0; i < passes.size(); i++) {
		sum_tubes += passes[i];
	}
	m_HE = new Segment[long(sum_tubes) * (m_segment + 1L)]; //(refrigerant_in.refrigerant);//microTube(double length,double t,double w,double p,double covered,double Pw,double Ph,int ports);	
	m_tube = new microTube(tdata);
	m_fin = new louverFin(fdata);
	q_re = refrigerant_in.qm;
}

MCHE::~MCHE()
{
	delete m_fin;
	delete m_tube;
	delete[] m_HE ;

}

void MCHE::geometryDesign()
{
	const double inchToM = 0.0254;
	double HE_L = m_tube->tube_length;
	double HE_H = m_tube->vertical_space * (sum_tubes- 1);
	int Nf = int(m_tube->tube_length / m_fin->finSpacing)+1;
	double thta = m_fin->theta / 180.0 * 3.1415;
	Afr_air = HE_L * HE_H;
	//double P = refrigerant_in.pressure;
	//double TensileStrength = 92E6;
	//double w = m_tube->port_width;
	//double H = m_tube->port_height;
	//double Np = m_tube->tube_ports;
	//double t2 = (m_tube->tube_thickness - H) / 2.0;
	//Ac_air = Afr_air - sum_tubes * (m_tube->tube_thickness * m_tube->tube_length +
	//	Nf *
	//	(
	//		m_fin->finTickness * (
	//			m_tube->vertical_space - m_tube->tube_thickness - m_fin->Ll
	//			)
	//		+ m_fin->Ll * m_fin->Lp * tan(thta)
	//		)
	//	);
	Ac_air = Afr_air - sum_tubes * (m_tube->tube_thickness * m_tube->tube_length + Nf * m_fin->finTickness * m_fin->Fl);

	double LdFactor = 0.8;
	Af = 2.0 * Nf * (m_fin->Fd * HE_H - m_tube->tube_thickness * m_tube->tube_width * sum_tubes);
	Af = 2.0 * Nf * (m_fin->Fd * m_fin->Fl * (sum_tubes-1.0));
	At = 2.0 * ((m_tube->tube_width + m_tube->tube_thickness)* m_tube->tube_length -Nf*m_fin->finTickness* m_tube->tube_width) * sum_tubes;
	Dh_air = 4.0 * Ac_air * m_fin->Fd / (Af + At);

	HE_volume_tube = sum_tubes*(m_tube->tube_thickness*m_tube->tube_width-m_tube->port_Ac)*m_tube->tube_length;
	HE_volume_fin = Af / 2.0 * m_fin->finTickness;

}

bool MCHE::solve(int iii)
{
	Q_total = 0.0;
	int dir[] = { 1,1 };
	//int dir[] = { 1,-1 };
	int locStart[] = { 1,1 };
	//int locStart[] = { 1,m_segment-1 };
	int dirInd = 0;
	refri tmp(refrigerant_in);


	refri air_tmp("air.ppf");
	refri air_tmp_in("air.ppf");

	double epsMax = 100;
	double deltaX = 0.1;
	air_tmp.pressure = air_in.P;
	air_tmp.tem = air_in.T;
	air_tmp.tpProperty();
	air_tmp_in.pressure = air_in.P;
	air_tmp_in.tem = air_in.T;
	air_tmp_in.tpProperty();
	double air_in_enthalpy=air_tmp.enthalpy;
	double air_density_i = air_tmp.density;
	double ho, f_air,fin_eta;
	double L_seg= m_tube->tube_length / m_segment;
	int ntube = 0;
	if (iii == 1) {
		int i = 0, j = 0;
		int last_j = j;
		for (size_t p = 0; p < passes.size(); p++) {
			double G = q_re / passes[p] / m_tube->port_Ac;
			for (int v_t = i; v_t < i+passes[p]; v_t++) {
				m_HE[v_t * (m_segment+1) + locStart[dirInd] - dir[dirInd]].P_rei_out = refrigerant_tmp_out.pressure;
				m_HE[v_t * (m_segment+1) + locStart[dirInd] - dir[dirInd]].h_rei_out = refrigerant_tmp_out.enthalpy;
				m_HE[v_t * (m_segment + 1) + locStart[dirInd] - dir[dirInd]].quality = refrigerant_tmp_out.quality;
				m_HE[v_t * (m_segment+1) + locStart[dirInd] - dir[dirInd]].G_refri = G;
			}

			ntube+= passes[p];
			for (; i < ntube; i++) {

				for (j=locStart[dirInd]; j <= m_segment && j >= 0; j += dir[dirInd]) {
					refrigerant_tmp_in.pressure = m_HE[i * (m_segment + 1) + j- dir[dirInd]].P_rei_out;
					refrigerant_tmp_in.enthalpy = m_HE[i * (m_segment + 1) + j - dir[dirInd]].h_rei_out;
					refrigerant_tmp_in.phProperty();
					
					refrigerant_tmp_out = refrigerant_tmp_in;
					air_tmp.pressure = air_tmp_in.pressure;
					air_tmp.tem = air_tmp_in.tem;
					air_tmp.tpProperty();

					if (refrigerant_tmp_in.quality >= 1.0) {
					//gas
						double eps = 2000;
						while (eps> epsMax) {
							jacobi_louver(f_air, ho, air_tmp_in, air_tmp);
							fin_eta= classic(ho);
							double hi;
							double fi;
							double Re_t = G * m_tube->Dh*2.0 / (refrigerant_tmp_in.viscosity+ refrigerant_tmp_out.viscosity);
							Gnielinski(m_HE[i * (m_segment + 1) + j], G, hi, fi);

							Blasius(Re_t, fi);

							double m_r = G * m_tube->port_Ac;
							double Cr = m_r * refrigerant_tmp_in.HeatCapacity_p * 1000;
							double m_a = air_in.v_in * air_tmp_in.density * L_seg * m_tube->vertical_space;
							double Ca = m_a * (air_tmp_in.HeatCapacity_p + air_tmp.HeatCapacity_p) / 2.0 * 1000;
							double Cmin = min(Cr, Ca);

							double rou_m_t = (refrigerant_tmp_in.densityGas + refrigerant_tmp_out.densityGas) / 2.0;
							//4.17
							double delta_re_P = 4.0 * L_seg / m_tube->Dh * fi * 0.5 * G * G / rou_m_t +
								0.5 * (G * G) * (1.0 / refrigerant_tmp_out.densityGas - 1.0 / refrigerant_tmp_in.densityGas);

							double Gc_air = air_in.v_in * air_density_i * Afr_air / Ac_air;
							double sigA = Ac_air / Afr_air;
							double Ke = KeCal(sigA);
							double Kc = KcCal(sigA);
							double air_density_o = air_tmp.density;
							//double delta_air_P = 4.0 * f_air * m_fin->Fd / Dh_air*pow(air_in.v_in * Afr_air / Ac_air,2.0);
							double delta_air_P;
							delta_air_P = (At + Af) / Ac_air * pow(Gc_air, 2.0) / air_density_i * f_air / 2.0;
							delta_air_P += pow(Gc_air, 2.0) / 2.0 / air_density_i * (1.0 - sigA * sigA + Kc);
							delta_air_P -= pow(Gc_air, 2.0) / 2.0 / air_density_o * (1.0 - sigA * sigA -Ke);
							delta_air_P += pow(Gc_air, 2.0) * (1.0/ air_density_o -1.0/ air_density_i);
							double KA = 1.0 / (1.0 / (fin_eta * ho * (Af + At) / (double(m_segment * sum_tubes))) + 1.0 / (hi * L_seg * m_tube->P * m_tube->tube_ports)+(m_tube->tube_thickness-m_tube->port_height)/(2.0*m_tube->k * L_seg * m_tube->P * m_tube->tube_ports));
							double NTU = KA / Cmin;
							double epsilon;
							if (Ca < Cr) {
								epsilon = AirLessNTU(Ca / Cr, NTU);
							}
							else {
								epsilon = TubeLessNTU(Cr / Ca, NTU);
							}
							double Q = epsilon * Cmin * (refrigerant_tmp_in.tem - air_in.T);
							refrigerant_tmp_out.enthalpy = refrigerant_tmp_in.enthalpy - Q / (q_re / passes[p]) / 1000.0;
							refrigerant_tmp_out.pressure = refrigerant_tmp_in.pressure - delta_re_P;
							refrigerant_tmp_out.phProperty();
							air_tmp.enthalpy= air_in_enthalpy + Q / (air_in.v_in * air_density_i * Afr_air/ (double(m_segment * sum_tubes))) / 1000.0;
							air_tmp.pressure = air_in.P - delta_air_P;
							air_tmp.phProperty();
							eps = abs(m_HE[i * (m_segment + 1) + j].P_rei_out - refrigerant_tmp_out.pressure);
							m_HE[i * (m_segment + 1) + j].P_rei_out = refrigerant_tmp_out.pressure;
							m_HE[i * (m_segment + 1) + j].h_rei_out = refrigerant_tmp_out.enthalpy;
							m_HE[i * (m_segment + 1) + j].quality = refrigerant_tmp_out.quality;
							m_HE[i * (m_segment + 1) + j].Re = Re_t;
							m_HE[i * (m_segment + 1) + j].P_air_out = air_tmp.pressure;
							m_HE[i * (m_segment + 1) + j].h_air_out = air_tmp.enthalpy;
							m_HE[i * (m_segment + 1) + j].T_air_out = air_tmp.tem;
						}
						last_j = j;
						reChargeV += L_seg * m_tube->port_Ac* (refrigerant_tmp_out.density+refrigerant_tmp_in.density)/2.0;

					}
					else if(refrigerant_tmp_in.quality<=0.0) {
						//liq
						double eps = 2000;
						while (eps  > epsMax) {
							jacobi_louver(f_air, ho, air_tmp_in, air_tmp);
							fin_eta = classic(ho);
							double hi;
							double fi;
							double Re_t = G * m_tube->Dh / (refrigerant_tmp_in.viscosity);
							Gnielinski(m_HE[i * (m_segment + 1) + j], G, hi, fi);

							Blasius(Re_t, fi);
							double rou_m_t = (refrigerant_tmp_in.densityLiq + refrigerant_tmp_out.densityLiq) / 2.0;
							double delta_re_P = 4.0 * L_seg / m_tube->Dh * fi * 0.5 * G * G / rou_m_t +
								0.5 * (G * G) * (1 / refrigerant_tmp_out.densityLiq - 1 / refrigerant_tmp_in.densityLiq);

							//air pres drop
							double Gc_air = air_in.v_in * air_density_i * Afr_air / Ac_air;
							double sigA = Ac_air / Afr_air;
							double Ke = KeCal(sigA);
							double Kc = KcCal(sigA);
							double air_density_o = air_tmp.density;
							double delta_air_P = (At + Af) / Ac_air * pow(Gc_air, 2.0) / air_density_i * f_air / 2.0;
							delta_air_P += pow(Gc_air, 2.0) / 2.0 / air_density_i * (1.0 - sigA * sigA + Kc);
							delta_air_P -= pow(Gc_air, 2.0) / 2.0 / air_density_o * (1.0 - sigA * sigA - Ke);
							delta_air_P += pow(Gc_air, 2.0) * (1.0 / air_density_o - 1.0 / air_density_i);

							//double delta_air_P=

							double m_r = G * m_tube->port_Ac;
							double Cr = m_r * refrigerant_tmp_in.HeatCapacity_p * 1000;
							double m_a = air_in.v_in * air_tmp_in.density * L_seg * m_tube->vertical_space;
							double Ca = m_a * (air_tmp_in.HeatCapacity_p+air_tmp.HeatCapacity_p)/2.0*1000;
							double Cmin = min(Cr, Ca);
							double KA = 1.0 / (1.0 / (fin_eta * ho * (Af + At) / (double(m_segment * sum_tubes))) + 1.0 / (hi * L_seg * m_tube->P * m_tube->tube_ports)+(m_tube->tube_thickness - m_tube->port_height) / (2.0 * m_tube->k * L_seg * m_tube->P * m_tube->tube_ports));
							double NTU = KA / Cmin;
							double epsilon;
							if (Ca < Cr) {
								epsilon = AirLessNTU(Ca / Cr, NTU);
							}
							else {
								epsilon = TubeLessNTU(Cr / Ca, NTU);
							}
							double Q = epsilon * Cmin * (refrigerant_tmp_in.tem - air_in.T);

							refrigerant_tmp_out.enthalpy = refrigerant_tmp_in.enthalpy - Q / (q_re / passes[p]) / 1000.0;
							refrigerant_tmp_out.pressure = refrigerant_tmp_in.pressure - delta_re_P;
							refrigerant_tmp_out.phProperty();
							air_tmp.enthalpy = air_in_enthalpy + Q / (air_in.v_in * air_density_i * Afr_air / (double(m_segment * sum_tubes))) / 1000.0;
							air_tmp.pressure = air_in.P - delta_air_P;
							air_tmp.phProperty();

							eps = abs(m_HE[i * (m_segment + 1) + j].P_rei_out - refrigerant_tmp_out.pressure);
							m_HE[i * (m_segment + 1) + j].P_rei_out = refrigerant_tmp_out.pressure;
							m_HE[i * (m_segment + 1) + j].h_rei_out = refrigerant_tmp_out.enthalpy;
							m_HE[i * (m_segment + 1) + j].quality = refrigerant_tmp_out.quality;
							m_HE[i * (m_segment + 1) + j].Re = Re_t;

							m_HE[i * (m_segment + 1) + j].P_air_out = air_tmp.pressure;
							m_HE[i * (m_segment + 1) + j].h_air_out = air_tmp.enthalpy;
							m_HE[i * (m_segment + 1) + j].T_air_out = air_tmp.tem;
						}
						last_j = j;

						reChargeL += L_seg * m_tube->port_Ac * (refrigerant_tmp_out.density + refrigerant_tmp_in.density)/2.0;

					
					}
					else {
					
					double eps =2000;
					while (eps > epsMax) {
						jacobi_louver(f_air, ho, air_tmp_in, air_tmp);
						fin_eta = classic(ho);
						double hi;
						double delta_re_P;
						m_HE[i * (m_segment + 1) + j].L = L_seg;
						double Rel;
						Friedel(m_HE[i * (m_segment + 1) + j], G, delta_re_P, Rel);
						//dp
						double ai= alpha(refrigerant_tmp_in);
						double ao= alpha(refrigerant_tmp_out);
						delta_re_P = delta_re_P +
							G * G*(pow(refrigerant_tmp_out.quality, 2.0) / (ao * refrigerant_tmp_out.densityGas) +
								pow(1.0 - refrigerant_tmp_out.quality, 2.0) / (1.0 - ao) / refrigerant_tmp_out.densityLiq
								- (
									pow(refrigerant_tmp_in.quality, 2.0) / (ai * refrigerant_tmp_in.densityGas) +
									pow(1.0 - refrigerant_tmp_in.quality, 2.0) / (1.0 - ai) / refrigerant_tmp_in.densityLiq

									));

						//air pres drop
						double Gc_air = air_in.v_in * air_density_i * Afr_air / Ac_air;
						double sigA = Ac_air / Afr_air;
						double Ke = KeCal(sigA);
						double Kc = KcCal(sigA);
						double air_density_o = air_tmp.density;
						double delta_air_P = (At + Af) / Ac_air * pow(Gc_air, 2.0) / air_density_i * f_air / 2.0;
						delta_air_P += pow(Gc_air, 2.0) / 2.0 / air_density_i * (1.0 - sigA * sigA + Kc);
						delta_air_P -= pow(Gc_air, 2.0) / 2.0 / air_density_o * (1.0 - sigA * sigA - Ke);
						delta_air_P += pow(Gc_air, 2.0) * (1.0 / air_density_o - 1.0 / air_density_i);


						refrigerant_tmp_out.pressure = refrigerant_tmp_in.pressure - delta_re_P;
						//refrigerant_tmp_out.tpProperty();

						condShah2016(m_HE[i * (m_segment + 1) + j], G, hi);

						double m_r = G * m_tube->port_Ac;
						//double Cr = m_r * refrigerant_tmp_in.HeatCapacity_p;
						double m_a = air_in.v_in * air_tmp_in.density * L_seg * m_tube->vertical_space;
						double Ca = m_a * (air_tmp_in.HeatCapacity_p + air_tmp.HeatCapacity_p) / 2.0 * 1000;
						double Cmin = Ca;
						double KA = 1.0 / (1.0 / (fin_eta * ho * (Af + At) / (double(m_segment * sum_tubes))) + 1.0 / (hi * L_seg * m_tube->P * m_tube->tube_ports)+ 
							(m_tube->tube_thickness - m_tube->port_height) / (2.0 * m_tube->k * L_seg * m_tube->P * m_tube->tube_ports));
						double NTU = KA / Cmin;
						double epsilon;

						epsilon = AirLessNTU(0.0, NTU);

						double Q = epsilon * Cmin * (refrigerant_tmp_in.tem - air_in.T);
						refrigerant_tmp_out.enthalpy = refrigerant_tmp_in.enthalpy - Q / (q_re / passes[p]) / 1000.0;
						refrigerant_tmp_out.pressure = refrigerant_tmp_in.pressure - delta_re_P;
						refrigerant_tmp_out.phProperty();

						eps = abs(m_HE[i * (m_segment + 1) + j].P_rei_out - refrigerant_tmp_out.pressure);

						m_HE[i * (m_segment + 1) + j].P_rei_out = refrigerant_tmp_out.pressure;
						m_HE[i * (m_segment + 1) + j].h_rei_out = refrigerant_tmp_out.enthalpy;
						m_HE[i * (m_segment + 1) + j].P_air_out = air_in.P - delta_air_P;
						m_HE[i * (m_segment + 1) + j].h_air_out = air_in_enthalpy + Q / (air_density_i * air_in.v_in * (Afr_air / (sum_tubes * m_segment)));
						m_HE[i * (m_segment + 1) + j].quality = refrigerant_tmp_out.quality;
						m_HE[i * (m_segment + 1) + j].Re = 0.0;
						air_tmp.pressure = air_in.P - delta_air_P;
						air_tmp.enthalpy = air_in_enthalpy + Q / (air_density_i * air_in.v_in * (Afr_air / (sum_tubes * m_segment))) / 1000.0;
						air_tmp.phProperty();

						m_HE[i * (m_segment + 1) + j].T_air_out = air_tmp.tem;
					}
					double ai = alpha(refrigerant_tmp_in);
					double ao = alpha(refrigerant_tmp_out);
					if (refrigerant_tmp_in.quality > 0.05) {
						reChargeVL += L_seg * m_tube->port_Ac * (refrigerant_tmp_out.densityLiq + refrigerant_tmp_in.densityLiq) * (1.0 - (ao + ai) / 2.0) / 2.0;
						reChargeVL += L_seg * m_tube->port_Ac * (refrigerant_tmp_out.densityGas + refrigerant_tmp_in.densityGas) * (ao + ai) / 4.0;
					}
					else { 
						reChargeL += L_seg * m_tube->port_Ac * (refrigerant_tmp_out.density + refrigerant_tmp_in.density) / 2.0; 
					}
						last_j = j;
					
					}

				}
				
				
			}
			dirInd = 1 - dirInd;
		
		}
	
	}
	else if (iii == 2) {

		int i = 0, j = 0;
		int last_j = j;
		delete[] m_HE;
		blockeds =bbk;
		double maxEps=300;
		//std::valarray<double> G_dis(blockeds);
		m_HE = new Segment[long(sum_tubes*blockeds) * (m_segment  + 1L)];

		for (size_t p = 0; p < passes.size(); p++) {
			double G = q_re / passes[p] / m_tube->port_Ac;
			std::valarray<double> G_dis(G, blockeds);
			std::valarray<double> Dp_refri(1.0, blockeds);
			for (int v_t = i; v_t < i + passes[p] * blockeds; v_t++) {
				m_HE[v_t * (m_segment + 1) + locStart[dirInd] - dir[dirInd]].P_rei_out = refrigerant_tmp_out.pressure;
				m_HE[v_t * (m_segment + 1) + locStart[dirInd] - dir[dirInd]].h_rei_out = refrigerant_tmp_out.enthalpy;
				m_HE[v_t * (m_segment + 1) + locStart[dirInd] - dir[dirInd]].quality = refrigerant_tmp_out.quality;
				m_HE[v_t * (m_segment + 1) + locStart[dirInd] - dir[dirInd]].G_refri = G;
			}

			ntube += passes[p];

			air_tmp.pressure = air_in.P;
			air_tmp.tem = air_in.T;
			air_tmp.tpProperty();

			double PressureIn = refrigerant_tmp_out.pressure;
			double meanP = 4000;
			int count = 0;
			double reChargeV_t, reChargeVL_t, reChargeL_t;
			do {
				
				reChargeV_t = reChargeVL_t = reChargeL_t = 0.0;
				for (int b = 0; b < blockeds; b++) {
					
					for (j = locStart[dirInd]; j <= m_segment && j >= 0; j += dir[dirInd]) {
						refrigerant_tmp_in.pressure = m_HE[(i+b) * (m_segment + 1) + j - dir[dirInd]].P_rei_out;
						refrigerant_tmp_in.enthalpy = m_HE[(i+b) * (m_segment + 1) + j - dir[dirInd]].h_rei_out;
						refrigerant_tmp_in.phProperty();

						refrigerant_tmp_out = refrigerant_tmp_in;
						if (b == 0) {
							air_tmp_in.pressure = air_in.P;
							air_tmp_in.tem = air_in.T;
						}
						else {
							air_tmp_in.pressure = m_HE[(i + b - 1) * (m_segment + 1) + j].P_air_out;
							air_tmp_in.tem = m_HE[(i + b - 1) * (m_segment + 1) + j].T_air_out;
						}

						air_tmp_in.tpProperty();
						m_HE[(i + b) * (m_segment + 1) + j].P_rei_out = 0.0;
						if (refrigerant_tmp_in.quality >= 1.0) {
							//gas
							double eps = 2000;
							
							while (eps > epsMax) {
								jacobi_louver(f_air, ho, air_tmp_in, air_tmp);
								fin_eta = classic(ho);
								double hi;
								double fi;
								double Re_t = G_dis[b] * m_tube->Dh * 2.0 / (refrigerant_tmp_in.viscosity + refrigerant_tmp_out.viscosity);
								Gnielinski(m_HE[i * (m_segment + 1) + j], G, hi, fi);

								Blasius(Re_t, fi);

								double m_r = G_dis[b] * m_tube->port_Ac / blockeds;
								double Cr = m_r * refrigerant_tmp_in.HeatCapacity_p * 1000;
								double m_a = air_in.v_in * air_tmp_in.density * L_seg * m_tube->vertical_space;
								double Ca = m_a * (air_tmp_in.HeatCapacity_p + air_tmp.HeatCapacity_p) / 2.0 * 1000;
								double Cmin = min(Cr, Ca);

								double rou_m_t = (refrigerant_tmp_in.densityGas + refrigerant_tmp_out.densityGas) / 2.0;
								//4.17
								double delta_re_P = 4.0 * L_seg / m_tube->Dh * fi * 0.5 * G_dis[b] * G_dis[b] / rou_m_t;
									

								//if (Re_t > 2300) {
								//	double c = 1.75;
								//	m_HE[(i + b) * (m_segment + 1) + j].Rs_Tu = delta_re_P / pow(G_dis[b], c);
								//	m_HE[(i + b) * (m_segment + 1) + j].Rs_La = 0.0;
								//}
								//else {
								//	double c = 1.0;
								//	m_HE[(i + b) * (m_segment + 1) + j].Rs_La = delta_re_P / pow(G_dis[b], c);
								//	m_HE[(i + b) * (m_segment + 1) + j].Rs_Tu = 0.0;
								//}
								//m_HE[(i + b) * (m_segment + 1) + j].Rs_a = 0.5  * (1.0 / refrigerant_tmp_out.densityGas - 1.0 / refrigerant_tmp_in.densityGas);
								delta_re_P+= 0.5 * (G_dis[b] * G_dis[b]) * (1.0 / refrigerant_tmp_out.densityGas - 1.0 / refrigerant_tmp_in.densityGas);

								double Gc_air = air_in.v_in * air_density_i * Afr_air / Ac_air;
								double sigA = Ac_air / Afr_air;
								double Ke = KeCal(sigA);
								double Kc = KcCal(sigA);
								double air_density_o = air_tmp.density;
								//double delta_air_P = 4.0 * f_air * m_fin->Fd / Dh_air*pow(air_in.v_in * Afr_air / Ac_air,2.0);
								double delta_air_P = (At + Af) / blockeds / Ac_air * pow(Gc_air, 2.0) / air_density_i * f_air / 2.0;
								//delta_air_P = (At + Af) / Ac_air * pow(air_in.v_in * air_density_i, 2.0) / air_density_i * f_air / 2.0;
								delta_air_P += pow(Gc_air, 2.0) / 2.0 / air_density_i * (1.0 - sigA * sigA + Kc) / blockeds;
								delta_air_P -= pow(Gc_air, 2.0) / 2.0 / air_density_o * (1.0 - sigA * sigA - Ke) / blockeds;
								delta_air_P += pow(Gc_air, 2.0) * (1.0 / air_density_o - 1.0 / air_density_i) / blockeds;





								double KA = 1.0 / (1.0 / (fin_eta * ho * (Af + At) / blockeds / (double(m_segment * sum_tubes))) + 1.0 / (hi * L_seg / blockeds * m_tube->P * m_tube->tube_ports) + (m_tube->tube_thickness - m_tube->port_height) / (2.0 * m_tube->k * L_seg * m_tube->P * m_tube->tube_ports / blockeds));
								double NTU = KA / Cmin;
								double epsilon;
								if (Ca < Cr) {
									epsilon = AirLessNTU(Ca / Cr, NTU);
								}
								else {
									epsilon = TubeLessNTU(Cr / Ca, NTU);
								}
								double Q = epsilon * Cmin * (refrigerant_tmp_in.tem - air_tmp_in.tem);
								if (Q <= 0.0) return false;
								refrigerant_tmp_out.enthalpy = refrigerant_tmp_in.enthalpy - Q / (m_r) / 1000.0;
								refrigerant_tmp_out.pressure = refrigerant_tmp_in.pressure - delta_re_P;
								refrigerant_tmp_out.phProperty();
								air_tmp.pressure = air_tmp_in.pressure - delta_air_P;
								air_tmp.enthalpy = air_tmp_in.enthalpy + Q / (air_density_i * air_in.v_in * (Afr_air / (sum_tubes * m_segment))) / 1000.0;
								air_tmp.phProperty();

								eps = abs(m_HE[(i + b) * (m_segment + 1) + j].P_rei_out - refrigerant_tmp_out.pressure);
								m_HE[(i + b) * (m_segment + 1) + j].P_rei_out = refrigerant_tmp_out.pressure;
								m_HE[(i + b) * (m_segment + 1) + j].h_rei_out = refrigerant_tmp_out.enthalpy;
								m_HE[(i + b) * (m_segment + 1) + j].quality = refrigerant_tmp_out.quality;
								m_HE[(i + b) * (m_segment + 1) + j].Re = Re_t;
								m_HE[(i + b) * (m_segment + 1) + j].P_air_out = air_tmp.pressure;
								m_HE[(i + b) * (m_segment + 1) + j].h_air_out = air_tmp.enthalpy;
								m_HE[(i + b) * (m_segment + 1) + j].T_air_out = air_tmp.tem;
								m_HE[(i + b) * (m_segment + 1) + j].ht_air = ho;
								if (Re_t > 2300) {
									double c = 1.75;
									m_HE[(i + b) * (m_segment + 1) + j].Rs_Tu = delta_re_P / pow(G_dis[b], c);
									m_HE[(i + b) * (m_segment + 1) + j].Rs_La = 0.0;
								}
								else {
									double c = 1.0;
									m_HE[(i + b) * (m_segment + 1) + j].Rs_La = delta_re_P / pow(G_dis[b], c);
									m_HE[(i + b) * (m_segment + 1) + j].Rs_Tu = 0.0;
								}
							}
							last_j = j;
							reChargeV_t += passes[p] * L_seg * m_tube->port_Ac / blockeds * (refrigerant_tmp_out.density + refrigerant_tmp_in.density) / 2.0;
							//refrigerant_tmp_in = refrigerant_tmp_out;

						}
						else if (refrigerant_tmp_in.quality <= 0.0) {
							//liq
							double eps = 2000;
							while (eps > epsMax) {
								jacobi_louver(f_air, ho, air_tmp_in, air_tmp);
								fin_eta = classic(ho);
								double hi;
								double fi;
								double Re_t = G_dis[b] * m_tube->Dh / (refrigerant_tmp_in.viscosity);
								Gnielinski(m_HE[i * (m_segment + 1) + j], G, hi, fi);

								Blasius(Re_t, fi);
								double rou_m_t = (refrigerant_tmp_in.densityLiq + refrigerant_tmp_out.densityLiq) / 2.0;
								double delta_re_P = 4.0 * L_seg / m_tube->Dh * fi * 0.5 * G_dis[b] * G_dis[b] / rou_m_t;
								//if (Re_t > 2300) {
								//	double c = 1.75;
								//	m_HE[(i + b) * (m_segment + 1) + j].Rs_Tu = delta_re_P / pow(G_dis[b], c);
								//	m_HE[(i + b) * (m_segment + 1) + j].Rs_La = 0.0;
								//}
								//else {
								//	double c = 1.0;
								//	m_HE[(i + b) * (m_segment + 1) + j].Rs_La = delta_re_P / pow(G_dis[b], c);
								//	m_HE[(i + b) * (m_segment + 1) + j].Rs_Tu = 0.0;
								//}
								//m_HE[(i + b) * (m_segment + 1) + j].Rs_a = 0.5 * (1.0 / refrigerant_tmp_out.densityLiq - 1.0 / refrigerant_tmp_in.densityLiq);
								delta_re_P+=0.5 * (G_dis[b] * G_dis[b]) * (1.0/ refrigerant_tmp_out.densityLiq - 1.0 / refrigerant_tmp_in.densityLiq);

								//air pres drop
								double Gc_air = air_in.v_in * air_density_i * Afr_air / Ac_air;
								double sigA = Ac_air / Afr_air;
								double Ke = KeCal(sigA);
								double Kc = KcCal(sigA);
								double air_density_o = air_tmp.density;
								//double delta_air_P = 4.0 * f_air * m_fin->Fd / Dh_air*pow(air_in.v_in * Afr_air / Ac_air,2.0);
								double delta_air_P = (At + Af) / blockeds / Ac_air * pow(Gc_air, 2.0) / air_density_i * f_air / 2.0;
								//delta_air_P = (At + Af) / Ac_air * pow(air_in.v_in * air_density_i, 2.0) / air_density_i * f_air / 2.0;
								delta_air_P += pow(Gc_air, 2.0) / 2.0 / air_density_i * (1.0 - sigA * sigA + Kc) / blockeds;
								delta_air_P -= pow(Gc_air, 2.0) / 2.0 / air_density_o * (1.0 - sigA * sigA - Ke) / blockeds;
								delta_air_P += pow(Gc_air, 2.0) * (1.0 / air_density_o - 1.0 / air_density_i) / blockeds;

								//double delta_air_P=

								double m_r = G_dis[b] * m_tube->port_Ac / blockeds;
								double Cr = m_r * refrigerant_tmp_in.HeatCapacity_p * 1000;
								double m_a = air_in.v_in * air_tmp_in.density * L_seg * m_tube->vertical_space;
								double Ca = m_a * (air_tmp_in.HeatCapacity_p + air_tmp.HeatCapacity_p) / 2.0 * 1000;
								double Cmin = min(Cr, Ca);
								double KA = 1.0 / (1.0 / (fin_eta * ho * (Af + At) / blockeds / (double(m_segment * sum_tubes))) + 1.0 / (hi * L_seg / blockeds * m_tube->P * m_tube->tube_ports) + (m_tube->tube_thickness - m_tube->port_height) / (2.0 * m_tube->k * L_seg / blockeds * m_tube->P * m_tube->tube_ports));
								double NTU = KA / Cmin;
								double epsilon;
								if (Ca < Cr) {
									epsilon = AirLessNTU(Ca / Cr, NTU);
								}
								else {
									epsilon = TubeLessNTU(Cr / Ca, NTU);
								}
								double Q = epsilon * Cmin * (refrigerant_tmp_in.tem - air_tmp_in.tem);
								if (Q <= 0.0) return false;
								refrigerant_tmp_out.enthalpy = refrigerant_tmp_in.enthalpy - Q / (m_r) / 1000.0;
								refrigerant_tmp_out.pressure = refrigerant_tmp_in.pressure - delta_re_P;
								refrigerant_tmp_out.phProperty();
								air_tmp.pressure = air_tmp_in.pressure - delta_air_P;
								air_tmp.enthalpy = air_tmp_in.enthalpy + Q / (air_density_i * air_in.v_in * (Afr_air / (sum_tubes * m_segment))) / 1000.0;
								air_tmp.phProperty();

								eps = abs(m_HE[(i + b) * (m_segment + 1) + j].P_rei_out - refrigerant_tmp_out.pressure);
								m_HE[(i + b) * (m_segment + 1) + j].P_rei_out = refrigerant_tmp_out.pressure;
								m_HE[(i + b) * (m_segment + 1) + j].h_rei_out = refrigerant_tmp_out.enthalpy;
								m_HE[(i + b) * (m_segment + 1) + j].quality = refrigerant_tmp_out.quality;
								m_HE[(i + b) * (m_segment + 1) + j].Re = Re_t;
								m_HE[(i + b) * (m_segment + 1) + j].P_air_out = air_tmp.pressure;
								m_HE[(i + b) * (m_segment + 1) + j].h_air_out = air_tmp.enthalpy;
								m_HE[(i + b) * (m_segment + 1) + j].T_air_out = air_tmp.tem;
								m_HE[(i + b) * (m_segment + 1) + j].ht_air = ho;
								if (Re_t > 2300) {
									double c = 1.75;
									m_HE[(i + b) * (m_segment + 1) + j].Rs_Tu = delta_re_P / pow(G_dis[b], c);
									m_HE[(i + b) * (m_segment + 1) + j].Rs_La = 0.0;
								}
								else {
									double c = 1.0;
									m_HE[(i + b) * (m_segment + 1) + j].Rs_La = delta_re_P / pow(G_dis[b], c);
									m_HE[(i + b) * (m_segment + 1) + j].Rs_Tu = 0.0;
								}
							}
							last_j = j;

							reChargeL_t += passes[p] * L_seg * m_tube->port_Ac / blockeds * (refrigerant_tmp_out.density + refrigerant_tmp_in.density) / 2.0;
							//reChargeL += L_seg * m_tube->port_Ac * (refrigerant_tmp_out.density) ;
							//refrigerant_tmp_in = refrigerant_tmp_out;


						}
						else {

							double eps = 2000;
							while (eps > epsMax) {
								jacobi_louver(f_air, ho, air_tmp_in, air_tmp);
								fin_eta = classic(ho);
								double hi;
								double delta_re_P;
								double Re_t=G_dis[b] * m_tube->Dh / (refrigerant_tmp_in.viscosity);;
								m_HE[i * (m_segment + 1) + j].L = L_seg;
								double Rel;
								Friedel(m_HE[i * (m_segment + 1) + j], G_dis[b], delta_re_P, Rel);
								//dp
								double ai = alpha(refrigerant_tmp_in);
								double ao = alpha(refrigerant_tmp_out);

								delta_re_P = delta_re_P +
									G_dis[b] * G_dis[b] * (pow(refrigerant_tmp_out.quality, 2.0) / (ao * refrigerant_tmp_out.densityGas) +
										pow(1.0 - refrigerant_tmp_out.quality, 2.0) / (1.0 - ao) / refrigerant_tmp_out.densityLiq
										- (
											pow(refrigerant_tmp_in.quality, 2.0) / (ai * refrigerant_tmp_in.densityGas) +
											pow(1.0 - refrigerant_tmp_in.quality, 2.0) / (1.0 - ai) / refrigerant_tmp_in.densityLiq

											));


									//air pres drop
								double Gc_air = air_in.v_in * air_density_i * Afr_air / Ac_air;
								double sigA = Ac_air / Afr_air;
								double Ke = KeCal(sigA);
								double Kc = KcCal(sigA);
								double air_density_o = air_tmp.density;
								//double delta_air_P = 4.0 * f_air * m_fin->Fd / Dh_air*pow(air_in.v_in * Afr_air / Ac_air,2.0);
								double delta_air_P = (At + Af) / blockeds / Ac_air * pow(Gc_air, 2.0) / air_density_i * f_air / 2.0;
								//delta_air_P = (At + Af) / Ac_air * pow(air_in.v_in * air_density_i, 2.0) / air_density_i * f_air / 2.0;
								delta_air_P += pow(Gc_air, 2.0) / 2.0 / air_density_i * (1.0 - sigA * sigA + Kc) / blockeds;
								delta_air_P -= pow(Gc_air, 2.0) / 2.0 / air_density_o * (1.0 - sigA * sigA - Ke) / blockeds;
								delta_air_P += pow(Gc_air, 2.0) * (1.0 / air_density_o - 1.0 / air_density_i) / blockeds;


								refrigerant_tmp_out.pressure = refrigerant_tmp_in.pressure - delta_re_P;
								//refrigerant_tmp_out.tpProperty();

								condShah2016(m_HE[i * (m_segment + 1) + j], G_dis[b], hi);

								double m_r = G_dis[b] * m_tube->port_Ac / blockeds;
								double m_a = air_in.v_in * air_tmp_in.density * L_seg * m_tube->vertical_space;
								double Ca = m_a * (air_tmp_in.HeatCapacity_p + air_tmp.HeatCapacity_p) / 2.0 * 1000;
								double Cmin = Ca;
								double KA = 1.0 / (1.0 / (fin_eta * ho * (Af + At)/blockeds / (double(m_segment * sum_tubes))) + 1.0 / (hi * L_seg / blockeds * m_tube->P * m_tube->tube_ports) +
									(m_tube->tube_thickness - m_tube->port_height) / (2.0 * m_tube->k * L_seg * m_tube->P / blockeds * m_tube->tube_ports));
								double NTU = KA/ Cmin;
								double epsilon;

								epsilon = AirLessNTU(0.0, NTU);

								double Q = epsilon * Cmin * (refrigerant_tmp_in.tem - air_tmp_in.tem);
								if (Q <= 0.0) return false;
								refrigerant_tmp_out.enthalpy = refrigerant_tmp_in.enthalpy - Q / (m_r) / 1000.0;
								refrigerant_tmp_out.pressure = refrigerant_tmp_in.pressure - delta_re_P;
								refrigerant_tmp_out.phProperty();

								eps = abs(m_HE[(i + b) * (m_segment + 1) + j].P_rei_out - refrigerant_tmp_out.pressure);
								air_tmp.pressure = air_tmp_in.pressure - delta_air_P;
								air_tmp.enthalpy = air_tmp_in.enthalpy + Q / (air_density_i * air_in.v_in * (Afr_air / (sum_tubes * m_segment))) / 1000.0;
								air_tmp.phProperty();

								m_HE[(i + b) * (m_segment + 1) + j].P_rei_out = refrigerant_tmp_out.pressure;
								m_HE[(i + b) * (m_segment + 1) + j].h_rei_out = refrigerant_tmp_out.enthalpy;
								m_HE[(i + b) * (m_segment + 1) + j].quality = refrigerant_tmp_out.quality;
								m_HE[(i + b) * (m_segment + 1) + j].Re = Re_t;
								m_HE[(i + b) * (m_segment + 1) + j].P_air_out = air_tmp.pressure;
								m_HE[(i + b) * (m_segment + 1) + j].h_air_out = air_tmp.enthalpy;
								m_HE[(i + b) * (m_segment + 1) + j].T_air_out = air_tmp.tem;
								m_HE[(i + b) * (m_segment + 1) + j].ht_air = ho;
								if (Rel > 2300) {
									double c = 1.75;
									m_HE[(i + b) * (m_segment + 1) + j].Rs_Tu = delta_re_P / pow(G_dis[b], c);
									m_HE[(i + b) * (m_segment + 1) + j].Rs_La = 0.0;
								}
								else {
									double c = 1.0;
									m_HE[(i + b) * (m_segment + 1) + j].Rs_La = delta_re_P / pow(G_dis[b], c);
									m_HE[(i + b) * (m_segment + 1) + j].Rs_Tu = 0.0;
								}
							}
							double ai = alpha(refrigerant_tmp_in);
							double ao = alpha(refrigerant_tmp_out);
							//if (refrigerant_tmp_in.quality > 0.05) {
							//	reChargeVL_t += passes[p] * L_seg * m_tube->port_Ac / blockeds * (refrigerant_tmp_out.densityLiq + refrigerant_tmp_in.densityLiq) * (1.0 - (ao + ai) / 2.0) / 2.0;
							//	reChargeVL_t += passes[p] * L_seg * m_tube->port_Ac / blockeds * (refrigerant_tmp_out.densityGas + refrigerant_tmp_in.densityGas) * (ao + ai) / 4.0;
							//}
							//else {
							//	reChargeL_t += passes[p] * L_seg * m_tube->port_Ac * (refrigerant_tmp_out.density + refrigerant_tmp_in.density) / 2.0;
							//}
							if (refrigerant_tmp_in.quality > deltaX*0.5) {
								reChargeVL_t += passes[p] * L_seg * m_tube->port_Ac / blockeds * (refrigerant_tmp_out.densityLiq + refrigerant_tmp_in.densityLiq) * (1.0 - (ao + ai) / 2.0) / 2.0;
								reChargeVL_t += passes[p] * L_seg * m_tube->port_Ac / blockeds * (refrigerant_tmp_out.densityGas + refrigerant_tmp_in.densityGas) * (ao + ai) / 4.0;
							}
							else {
								reChargeL_t += passes[p] * L_seg * m_tube->port_Ac * (refrigerant_tmp_out.density + refrigerant_tmp_in.density) / 2.0;
							}

							deltaX = refrigerant_tmp_in.quality - refrigerant_tmp_out.quality;
							last_j = j;

						}

					}

					Dp_refri[b] = PressureIn- m_HE[(i + b) * (m_segment + 1) + last_j].P_rei_out;
				}
				//double meanDp = Dp_refri.sum() / blockeds;
				std::valarray<double> meanDp(Dp_refri.sum() / blockeds, blockeds);
				std::valarray<double> absDis=abs(Dp_refri - meanDp);
				maxEps =absDis[0];
				std::vector<double> R_stream(blockeds, 0.0);
				double sumR_stream=0.0;

				//2)Nd
				std::vector<double> Rs_L(blockeds, 0.0);
				std::vector<double> Rs_T(blockeds, 0.0);
				//std::vector<double> Rs_a(blockeds, 0.0);
				std::vector<double> G_calcu(blockeds, 0.0);
				double sumG_stream = 0.0;
				for (int s = locStart[dirInd]; s <= m_segment && s >= 0; s += dir[dirInd]) {
					for (int id = 0; id < blockeds; id++) {
						Rs_L[id] += m_HE[(i + id) * (m_segment + 1) + s].Rs_La;
						Rs_T[id] += m_HE[(i + id) * (m_segment + 1) + s].Rs_Tu;
					}
				
				}
				for (int id = 0; id < blockeds; id++) {
					maxEps = max(absDis[id], maxEps);
					G_calcu[id] = NTsolveG(Rs_T[id], Rs_L[id], 0.0, 1.75, 1.0, 2.0,meanDp[0], G_dis[id]);
					sumG_stream += G_calcu[id];
				}
				double releaseCo = 0.1;
				maxEps /= meanDp[0];
				for (int id = 0; id < blockeds; id++) {
					G_dis[id] = (1 - releaseCo) * G_dis[id] + releaseCo * G_calcu[id] / sumG_stream * q_re * blockeds / passes[p] / m_tube->port_Ac;
				}
				////////////////////////////////////////////////////////////////
				//3)Nd+1
				//double* Rs_L = new double[blockeds] {0.0};
				//double* Rs_T = new double[blockeds] {0.0};
				//std::vector<double> G_calcu(blockeds, 0.0);
				//double sumG_stream = 0.0;
				//for (int s = locStart[dirInd]; s <= m_segment && s >= 0; s += dir[dirInd]) {
				//	for (int id = 0; id < blockeds; id++) {
				//		Rs_L[id] += m_HE[(i + id) * (m_segment + 1) + s].Rs_La;
				//		Rs_T[id] += m_HE[(i + id) * (m_segment + 1) + s].Rs_Tu;
				//	}

				//}
				//mwArray x0(1, blockeds+1, mxDOUBLE_CLASS);
				//mwArray R_L(1, blockeds, mxDOUBLE_CLASS);
				//mwArray R_T(1, blockeds, mxDOUBLE_CLASS);
				//mwArray sumG(1, 1, mxDOUBLE_CLASS);
				//mwArray epsd(1, 1, mxDOUBLE_CLASS);
				//sumG(1, 1) = q_re * blockeds / passes[p] / m_tube->port_Ac;
				//epsd(1, 1) = 1E-6;
				//x0(1, blockeds + 1)= meanDp[0];

				//for (int id = 0; id < blockeds; id++) {
				//	maxEps = max(absDis[id], maxEps);
				//}
				//R_L.SetData(Rs_L, blockeds);
				//R_T.SetData(Rs_T, blockeds);
				//mwArray r(1, blockeds+1, mxDOUBLE_CLASS);
				//mwArray Fx(1, blockeds+1, mxDOUBLE_CLASS);
				//mwArray R(1,2*(blockeds+1), mxDOUBLE_CLASS);
				//NewtonToSolveGdp(3, r, Fx, R, x0, R_L, R_T, sumG, epsd);

				//double releaseCo = 0.5;
				//maxEps /= meanDp[0];
				////std::cout << r << endl;
				//double* Gn = new double[blockeds+1] {0.0};
				//r.GetData(Gn, r.NumberOfElements());
				//for (int id = 0; id < blockeds; id++) {
				//	G_dis[id] = (1 - releaseCo) * G_dis[id] + releaseCo * Gn[id];
				//}
				//delete[] Rs_L;
				//delete[] Rs_T;
				//delete[] Gn;
				////////////////////////////////////////////////////////////////
				// 
				 ////1)I=U/R
				/*double cstC = 1.7;
				for (int id = 0; id < blockeds; id++) {
					maxEps = max(absDis[id], maxEps);
					R_stream[id] =Dp_refri[id] / pow(G_dis[id], cstC);
					sumR_stream += pow(1.0/R_stream[id],1.0/cstC);
				}
				double releaseCo = 1;
				maxEps /= meanDp[0];
				
					for (int id = 0; id < blockeds; id++) {

						G_dis[id] = (1 - releaseCo) * G_dis[id] + releaseCo * pow(1.0 / R_stream[id], 1.0 / cstC) / sumR_stream * q_re * blockeds / passes[p] / m_tube->port_Ac;
					}*/
				
				std::cout << maxEps << std::endl;
				if (++count >= 10) break;
				
			} while (maxEps>0.01);
			reChargeV += reChargeV_t;
			reChargeVL += reChargeVL_t;
			reChargeL+=reChargeL_t;
			double h_out = 0.0;
			double P_out = 0.0;
			for (int db = 0; db < blockeds; db++) {
				h_out+=m_HE[(i + db) * (m_segment + 1) + last_j].h_rei_out*G_dis[db]*m_tube->port_Ac/blockeds;
				P_out+= m_HE[(i + db) * (m_segment + 1) + last_j].P_rei_out;
			}
			h_out =h_out/ (q_re/ passes[p]);
			P_out /= blockeds;
				dirInd = 1 - dirInd;

				for (int start = i + blockeds; start <i+ passes[p] * blockeds; start += blockeds) {
					memcpy(&m_HE[(start + 0) * (m_segment + 1) + 0], &m_HE[(i + 0) * (m_segment + 1) + 0], sizeof(Segment)* (m_segment+1)* blockeds);
				}
				i=i+ blockeds*passes[p];


		}

	}

		refrigerant_out = refrigerant_tmp_out;
		Q_total = abs(refrigerant_out.enthalpy - refrigerant_in.enthalpy) * 1000 * q_re;
		DP_re = abs(refrigerant_out.pressure - refrigerant_in.pressure);
		double t = refrigerant_out.tem;
		refrigerant_out.pSatProperty(1);
		DelSat = t - refrigerant_out.tem;
		return true;
	}	


