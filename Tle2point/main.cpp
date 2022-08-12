//
// main.cpp 
// This sample code demonstrates how to use the C++ classes in order
// to determine satellite position and look angles.
//
// Copyright (c) 2003-2013 Michael F. Henry
//
// 01/2013
//

#include "stdafx.h"
#include "tle2point.h"
#include "cNoradSDP4.h"

#include <stdio.h>

#include "coreLib.h"
#include "cOrbit.h"



// This namespace contains all the OrbitTools classes; be sure
// to use this namespace.
using namespace std;
using namespace Zeptomoby::OrbitTools;
using namespace TLE2Point;

/////////////////////////////////////////////////////////////////////////////
// Test routine to output position and velocity information
void PrintPosVel(const cTle &tle)
{
   cOrbit       orbit(tle);
   vector<cEci> Pos;

   // Calculate position, velocity
   // mpe = "minutes past epoch"
   for (int mpe = 0; mpe <= (360 * 4); mpe += 360)
   {
      // Get the position of the satellite at time "mpe"
      cEciTime eci = orbit.GetPosition(mpe);
    
      // Push the coordinates object onto the end of the vector.
      Pos.push_back(eci);
   }

}


/* 数组求逆：仅限3 * 3的数组
输入参数：
A：求逆前的二维数组
输出参数：
D：求逆后的二维数组
*/
bool inv(double A[3][3], double D[3][3])
{
	int n = 3;
	double C[3][6] = { 0 };
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			C[i][j] = A[i][j];
	for (int i = 0; i < n; i++)
		C[i][i + n] = 1;
	for (int k = 0; k < n; k++)
	{
		double max = abs(C[k][k]);
		int ii = k;
		for (int m = k + 1; m < n; m++)
			if (max < abs(C[m][k]))
			{
				max = abs(C[m][k]);
				ii = m;
			}
		for (int m = k; m < 2 * n; m++)
		{
			if (ii == k) break;
			double c;
			c = C[k][m];
			C[k][m] = C[ii][m];
			C[ii][m] = c;
		}
		if (C[k][k] != 1)
		{
			double bs = C[k][k];
			if (bs == 0)
			{
				return false;
			}
			C[k][k] = 1;
			for (int p = k + 1; p < n * 2; p++)
			{
				C[k][p] /= bs;
			}
		}
		for (int q = k + 1; q < n; q++)
		{
			double bs = C[q][k];
			for (int p = k; p < n * 2; p++)
			{
				C[q][p] -= bs * C[k][p];
			}
		}
	}
	for (int q = n - 1; q > 0; q--)
	{
		for (int k = q - 1; k > -1; k--)
		{
			double bs = C[k][q];
			for (int m = k + 1; m < 2 * n; m++)
			{
				C[k][m] -= bs * C[q][m];
			}
		}
	}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			D[i][j] = C[i][j + n];
	return true;
}
/*两个向量叉积：仅限大小为3的数组
输入参数：
A：一维数组，大小为3
B：一维数组，大小为3
输出参数：
C：一维数组，大小为3
*/
void cross(double A[], double B[], double C[3])
{
	double x = A[0] * B[2] - A[2] * B[1];//计算三阶行列式
	double y = A[2] * B[0] - A[0] * B[2];
	double z = A[0] * B[1] - A[1] * B[0];
	C[0] = x;
	C[1] = y;
	C[2] = z;
}
//两个向量叉积：仅限大小为3 * 1的数组
/*
输入参数：
A：二维数组，大小为3 * 1
B：二维数组，大小为3 * 1
输出参数：
C：二维数组，大小为3 * 1
*/
void cross(double A[3][1], double B[3][1], double C[3][1])
{
	double x = A[1][0] * B[2][0] - A[2][0] * B[1][0];//计算三阶行列式
	double y = A[2][0] * B[0][0] - A[0][0] * B[2][0];
	double z = A[0][0] * B[1][0] - A[1][0] * B[0][0];
	C[0][0] = x;
	C[1][0] = y;
	C[2][0] = z;
}

//两个3 * 3的数组乘积，输出3 * 3的数组
void mult33(double a[3][3], double b[3][3], double c[3][3])
{
	int i, j, k;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			c[i][j] = 0;
			for (k = 0; k < 3; k++)
			{
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}

//两个二维数组乘积：输入分别为3 * 3和3 * 1，输出3 * 1
void mult31(double a[3][3], double b[3][1], double c[3][1])
{
	int i, j, k;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 1; j++)
		{
			c[i][j] = 0;
			for (k = 0; k < 3; k++)
			{
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}

//两个数组乘积：输入分别为1 * 3和3 * 1，输出数值
double mult11(double a[3], double b[3][1])
{
	int i, j, k;
	double c = 0;

	for (i = 0; i < 3; i++)
	{
		c += a[i] * b[i][0];
	}
	return c;
}

//矩阵转置: 1 * 3 转 3 * 1
template<typename T >
void transpose(T a[3], T at[3][1])
{
	for (int i = 0; i < 3; i++)
	{
		at[i][0] = a[i];
	}
}

///矩阵转置: 3 * 1 转 1 * 3
template<typename T >
void transpose(T a[3][1], T at[3])
{
	int m = 3;
	int n = 1;
	int i, j, t;

	at[0] = a[0][0];
	at[1] = a[1][0];
	at[2] = a[2][0];
}

double computeHg(TIMEDATE m_time)  
{

	int year, month, day, hour, minute;
	double second;
	year = m_time.year;
	month = m_time.month;
	day = m_time.day;
	hour = m_time.hour;
	minute = m_time.minute;
	second = m_time.second;

	double Hg;
	double m_jdc, m_Tu, GMST, m_dT, Om, L, eps, delta_psi;

	//Calculate the Greenwich Mean Sidereal Time Angle
	m_jdc = 367 * year
		- floor(7 * (year + floor((month + 9) / 12)) / 4)
		+ floor(275 * month / 9) + day +
		(hour + minute / 60.0 + second / 3600.0) / 24.0 - 730531.5;
	m_Tu = m_jdc / 36525.0;
	GMST = 280.46061837 + 360.98564736629 * m_jdc + 0.000388 * m_Tu * m_Tu;
	GMST = fmod(GMST, 360.0);
	GMST = (GMST / 360.0) * 2 * Pai;

	//Considering precession, Nutation and 月亮轨道对黄道平均升交点的黄经
	m_dT = m_jdc - 2451545.0;
	Om = 125.04452 - 1934.136261 * m_dT;
	L = 280.4665 + 36000.7698 * m_dT;
	eps = 23.439 - 0.0000004 * m_dT;
	delta_psi = -0.000319 * sin(Om) - 0.000024 * sin(2 * L);

	//The correction to GMST in seconds of time is given by    
	GMST = GMST + delta_psi * cos(eps) * (WE * 3600);
	Hg = GMST;

	return Hg;
}



// 矩阵A与矩阵B相乘
void MatrixAB(double A[3][3], double B[3], double C[3])
{
	int i, j;
	for (i = 0; i < 3; i++)
	{
		C[i] = 0;
		for (j = 0; j < 3; j++)
			C[i] += A[i][j] * B[j];
	}

}



void getXYZFromLOSandEarth(double xo, double yo, double zo, double LOS[3][1], double H, double XYZ[3])
{
	double a = Re + H;
	double b = Rp + H;
	double c = a / b;
	double AA = pow(LOS[0][0], 2.0) + pow(LOS[1][0], 2) + c * c * pow(LOS[2][0], 2);
	double BB = xo * LOS[0][0] + yo * LOS[1][0] + c * c * zo * LOS[2][0];
	double CC = xo * xo + yo * yo + c * c * zo * zo - a * a;
	double tt = (-BB - sqrt(BB * BB - AA * CC)) / AA;
	XYZ[0] = xo + tt * (LOS[0][0]);
	XYZ[1] = yo + tt * (LOS[1][0]);
	XYZ[2] = zo + tt * (LOS[2][0]);
}

/*******************************
函数说明
根据地固坐标系下的位置矢量计算地理经纬度
输入量
XYZG：目标点在地固坐标系下的坐标
flag: 标志（0：假设高程为零；1：利用直接法计算地理经纬度及高度）
输出量
lat, long : 地理经纬度longitude and latitude
********************************/
void getBL(double XYZG[3], double flag, double& lat, double& lon, double h)
{
	double e = sqrt((Re * Re - Rp * Rp) / (Re * Re));
	if (flag == -1)
	{
		double X = XYZG[0];
		double Y = XYZG[1];
		double Z = XYZG[2];
		double L = atan2(Y, X);
		double p = -2 * Z / sqrt(X * X + Y * Y);
		double q = 1 / (1 - e * e) + Z * Z / (X * X + Y * Y)
			- Re * Re / (X * X + Y * Y) * pow(e, 4) / (1 - pow(e, 4));
		double r = -Z / sqrt(X * X + Y * Y) * 2 / (1 - e * e);
		double s = Z * Z / (X * X + Y * Y) / (1 - e * e);
		double PP = -q / 2;
		double QQ = (p * r - 4 * s) / 4;
		double RR = (-p) * (-p) * s / 8 + q * s / 2 - r * r / 8;
		double pp = QQ - PP * PP / 3;
		double qq = pow(PP, 3) * 2 / 27 - PP * QQ / 3 + RR;
		double kk = pow((-qq / 2 + sqrt(qq * qq / 4 + pow(pp, 3) / 27))
			, (1 / 3)) + pow((-qq / 2 - sqrt(qq * qq / 4 + pow(pp, 3) / 27))
				, (1 / 3)) - PP / 3;
		double bb = -sqrt(kk * kk - s);
		double aa = (p * kk - r) / 2 / bb;
		double TT1 = (-(p / 2 + aa) + sqrt(((p / 2 + aa), 2) - 4 * (kk + bb))) / 2;
		double TT2 = (-(p / 2 + aa) - sqrt(((p / 2 + aa), 2) - 4 * (kk + bb))) / 2;
		double TT3 = (-(p / 2 - aa) + sqrt(((p / 2 - aa), 2) - 4 * (kk - bb))) / 2;
		double TT4 = (-(p / 2 - aa) - sqrt(((p / 2 - aa), 2) - 4 * (kk - bb))) / 2;

		double B = atan(TT4);
		double N = Re / sqrt(1 - e * e * pow((sin(B)), 2));
		h = X / cos(B) / cos(L) - N;
		lat = B;
		lon = L;
	}
	else {
		h = flag;
		double a = Re + h;
		double b = Rp + h;
		double c = sqrt((a * a - b * b) / (a * a));
		double B = atan2(XYZG[2], (1 - c * c) * sqrt(pow(XYZG[0], 2)
			+ pow(XYZG[1], 2)));
		double L = atan2(XYZG[1], XYZG[0]);
		lat = B;
		lon = L;
	}
}

void get2BL(double* polar, const POSITION& Pos)
{
	double tmp1, tmp2;
	tmp1 = (Pos.x * Pos.x + Pos.y * Pos.y) * (1 - ECCENT) * (1 - ECCENT);
	tmp2 = Pos.z / sqrt(tmp1);
	polar[0] = sqrt(Pos.x * Pos.x + Pos.y * Pos.y + Pos.z * Pos.z);
	polar[1] = atan2(Pos.y, Pos.x) * 180.0 / Pai;
	polar[2] = atan(tmp2) * 180.0 / Pai;
}

//求解视点矢量和地球的交点
void getXYZO(double satellite[3], double Los[3], double H, double xyzo[3])
{
	double a, b, c;
	a = Re + H;
	b = Rp + H;
	c = a / b;
	double AA, BB, CC, tt;
	AA = Los[0] * Los[0] + Los[1] * Los[1] + c * c * Los[2] * Los[2];
	BB = satellite[0] * Los[0] + satellite[1] * Los[1] + c * c * satellite[2] * Los[2];
	CC = satellite[0] * satellite[0] + satellite[1] * satellite[1]
		+ c * c * satellite[2] * satellite[2] - Re * Re;

	tt = (-BB - sqrt(BB * BB - AA * CC)) / AA;

	xyzo[0] = satellite[0] + tt * Los[0];
	xyzo[1] = satellite[1] + tt * Los[1];
	xyzo[2] = satellite[2] + tt * Los[2];


}

double GetJDE(int year, double day)
{
	// 1582 A.D.: 10 days removed from calendar
	// 3000 A.D.: Arbitrary error checking limit
	//assert((year > 1582) && (year < 3000));
	//assert((day >= 1.0) && (day < 367.0));

	// Now calculate Julian date

	year--;

	// Centuries are not leap years unless they divide by 400
	int A = (year / 100);
	int B = 2 - A + (A / 4);

	double NewYears = (int)(365.25 * year) +
		(int)(30.6001 * 14) +
		1720994.5 + B;  // 1720994.5 = Oct 30, year -1

	double m_Date = NewYears + day;
	return m_Date;
}


//位置由J2000到WGS84
POSITION TLE2Point::RJ2000toECF(AREAPARA areapara)
{
	POSITION SatPosECF, Sartemp;
	double theta = areapara.GMST;
	//TIMEDATE mytime = { 2021, 8 ,16, 0,0 , 0 };
	//TIMEDATE mytime = { 2019, 3 ,17, 3,21 , 52 };
	/*theta = computeHg(mytime);*/
	//double poseMat[9];// = { 0,0,0,0,0,0,0,0,0 }
	//QMGen(poseMat, areapara);
	//double postempmat[3][3] = { { poseMat[0] ,poseMat[1] ,poseMat[2] } ,{ poseMat[3] ,poseMat[4] ,poseMat[5] } ,{ poseMat[6] ,poseMat[7] ,poseMat[8] } };
	double postempmat[3][3] = { {1,0,0 } ,{ 0,1,0, } ,{ 0,0,1 } };
	//double Sarmat[3] = { SARpoint.x , SARpoint.y , SARpoint.z };
	//double Sarmattemp[3];
	//MatrixA0B(postempmat, Sarmat, Sarmattemp);
	Sartemp.x = areapara.X;
	Sartemp.y = areapara.Y;
	Sartemp.z = areapara.Z;
	SatPosECF.x = cos(theta) * Sartemp.x + sin(theta) * Sartemp.y;
	SatPosECF.y = -sin(theta) * Sartemp.x + cos(theta) * Sartemp.y;
	SatPosECF.z = Sartemp.z;
	return SatPosECF;
}

//速度J2000转换为地固坐标系
POSITION TLE2Point::VJ2000toECF(AREAPARA areapara)
{
	POSITION SatVelECF;
	double Hg = areapara.GMST;
	/*Hg = computeHg(mytime);*/


	//SatVelECF.x = areapara.Vx - WE * areapara.Vy;
	//SatVelECF.y = areapara.Vy - WE * areapara.Vx;
	//SatVelECF.z = areapara.Vz;

	double B_ma[3][3]; double invB_ma[3][3];


	B_ma[0][0] = cos(Hg);
	B_ma[0][1] = sin(Hg);
	B_ma[0][2] = 0;
	B_ma[1][0] = -sin(Hg);
	B_ma[1][1] = cos(Hg);
	B_ma[1][2] = 0;
	B_ma[2][0] = 0;
	B_ma[2][1] = 0;
	B_ma[2][2] = 1;

	double VB_ma[3][3];
	//vB_ma[0][0] = -sin(Hg)*WE;
	//vB_ma[0][1] = cos(Hg)*WE;
	//vB_ma[0][2] = 0;
	//vB_ma[1][0] = -cos(Hg)*WE;
	//vB_ma[1][1] = -sin(Hg)*WE;
	//vB_ma[1][2] = 0;
	//vB_ma[2][0] = 0;
	//vB_ma[2][1] = 0;
	//vB_ma[2][2] = 0;

	VB_ma[0][0] = -sin(Hg) * WE;
	VB_ma[0][1] = cos(Hg) * WE;
	VB_ma[0][2] = 0;
	VB_ma[1][0] = -cos(Hg) * WE;
	VB_ma[1][1] = -sin(Hg) * WE;
	VB_ma[1][2] = 0;
	VB_ma[2][0] = 0;
	VB_ma[2][1] = 0;
	VB_ma[2][2] = 0;

	double vj2000[3] = { areapara.Vx, areapara.Vy, areapara.Vz };
	double j2000[3] = { areapara.X, areapara.Y, areapara.Z };
	double vwgs841[3], vwgs842[3];
	MatrixAB(VB_ma, j2000, vwgs841);




	MatrixAB(B_ma, vj2000, vwgs842);


	double ECF[3];
	ECF[0] = vwgs841[0] + vwgs842[0];
	ECF[1] = vwgs841[1] + vwgs842[1];
	ECF[2] = vwgs841[2] + vwgs842[2];

	SatVelECF.x = ECF[0];
	SatVelECF.y = ECF[1];
	SatVelECF.z = ECF[2];

	return SatVelECF;
}

void getXYZGfroXYZO(double xyzo[3], int epochYear, double epochDay, double xyzg[3],double GMST)
{

	double JED, T;
	JED = GetJDE(epochYear, epochDay);
	T = (JED - 2451545) / 36525;

	//1.0岁差矩阵
	double fi, Z, cita;
	double D_ma[3][3];
	fi = (2306.2181 * T + 0.30188 * T * T + 0.017998 * T * T * T) * Pai / (180 * 60 * 60);
	Z = (2306.2181 * T + 1.09468 * T * T + 0.018203 * T * T * T) * Pai / (180 * 60 * 60);
	cita = (2004.3109 * T - 0.42665 * T * T - 0.041833 * T * T * T) * Pai / (180 * 60 * 60);
	D_ma[0][0] = cos(Z) * cos(cita) * cos(fi) - sin(Z) * sin(fi);
	D_ma[0][1] = -cos(Z) * cos(cita) * sin(fi) - sin(Z) * cos(fi);
	D_ma[0][2] = -cos(Z) * sin(cita);
	D_ma[1][0] = sin(Z) * cos(cita) * cos(fi) + cos(Z) * sin(fi);
	D_ma[1][1] = -sin(Z) * cos(cita) * sin(fi) + cos(Z) * cos(fi);
	D_ma[1][2] = -sin(Z) * sin(cita);
	D_ma[2][0] = sin(cita) * cos(fi);
	D_ma[2][1] = -sin(cita) * sin(fi);
	D_ma[2][2] = cos(cita);


	//2.0 章动矩阵
	double womiga_m, Ms, Mss, F, D, fi_delta
		, yipuxilong_delta, yipuxilong_average, yipuxilong;
	double C_ma[3][3];
	womiga_m = (450160.280 - (5 * 1296000 + 482890.539) * T
		+ 7.455 * T * T + 0.008 * T * T * T) * Pai / (180 * 60 * 60);
	Ms = (1287099.804 + (99 * 1296000 + 1292581.244) * T
		- 0.577 * T * T - 0.012 * T * T * T) * Pai / (180 * 60 * 60);
	Mss = (485866.733 + (1325 * 1296000 + 715922.633) * T
		+ 31.310 * T * T + 0.064 * T * T * T) * Pai / (180 * 60 * 60);
	F = (335778.877 + (1342 * 1296000 + 295263.137) * T
		- 13.257 * T * T + 0.011 * T * T * T) * Pai / (180 * 60 * 60);
	D = (1072261.307 + (1236 * 1296000 + 1105601.328) * T
		- 6.891 * T * T + 0.019 * T * T * T) * Pai / (180 * 60 * 60);
	fi_delta = ((-17.1996 - 0.1742 * T) * sin(womiga_m)
		+ (-1.3187 - 0.00016 * T) * sin(2 * F - 2 * D + 2 * womiga_m)
		- 0.2274 * sin(2 * F + 2 * womiga_m) + 0.2062 * sin(2 * womiga_m) + 0.1426 * sin(Ms)) * Pai / (180 * 60 * 60);
	yipuxilong_delta = (9.2025 * cos(womiga_m)
		+ 0.5736 * cos(2 * F - 2 * D + 2 * womiga_m) + 0.0977 * cos(2 * F + 2 * womiga_m)
		- 0.0895 * cos(2 * womiga_m)) * Pai / (180 * 60 * 60);
	yipuxilong_average = (84381.448 - 46.8150 * T - 0.00059 * T * T + 0.001813 * T * T * T) * Pai / (180 * 60 * 60);
	yipuxilong = yipuxilong_average + yipuxilong_delta;
	C_ma[0][0] = cos(fi_delta);
	C_ma[0][1] = -sin(fi_delta) * cos(yipuxilong_average);
	C_ma[0][2] = -sin(fi_delta) * sin(yipuxilong_average);
	C_ma[1][0] = cos(yipuxilong) * sin(fi_delta);
	C_ma[1][1] = cos(yipuxilong) * cos(fi_delta) * cos(yipuxilong_average)
		+ sin(yipuxilong) * sin(yipuxilong_average);
	C_ma[1][2] = cos(yipuxilong) * cos(fi_delta) * sin(yipuxilong_average)
		- sin(yipuxilong) * cos(yipuxilong_average);
	C_ma[2][0] = sin(yipuxilong) * sin(fi_delta);
	C_ma[2][1] = sin(yipuxilong) * cos(fi_delta) * cos(yipuxilong_average)
		- cos(yipuxilong) * sin(yipuxilong_average);
	C_ma[2][2] = sin(yipuxilong) * cos(fi_delta) * sin(yipuxilong_average)
		+ cos(yipuxilong) * cos(yipuxilong_average);


	//3.0 格林威治角旋转矩
	double Hg;
	double B_ma[3][3];
	Hg = GMST;
	B_ma[0][0] = cos(Hg);
	B_ma[0][1] = sin(Hg);
	B_ma[0][2] = 0;
	B_ma[1][0] = -sin(Hg);
	B_ma[1][1] = cos(Hg);
	B_ma[1][2] = 0;
	B_ma[2][0] = 0;
	B_ma[2][1] = 0;
	B_ma[2][2] = 1;

	//4.0 极移忽略
	double j2000[3];
	double tmp[3], tmp1[3];

	j2000[0] = xyzo[0];
	j2000[1] = xyzo[1];
	j2000[2] = xyzo[2];

	MatrixAB(D_ma, j2000, tmp);
	MatrixAB(C_ma, tmp, tmp1);
	MatrixAB(B_ma, tmp1, xyzg);


}

void mygetpoint(int epochYear, double epochDay,  AREAPARA areapara
	, POSITION TSatPosECF, POSITION TSatVelECF, EARTHPOINT& earthpoint)

{
	double sitaX = 0;
	double sitaY = 0;
	double Losc[3]; double AA[3][3]; double Losd[3];
	sitaX = 0;
	Losc[0] = 0;
	Losc[1] = 0;
	Losc[2] = 1;
	//		printf("Losc:%f,%f,%f\n",Losc[0],Losc[1],Losc[2]);

	sitaY = 0;
	AA[0][0] = cos(sitaY); AA[0][1] = 0; AA[0][2] = sin(sitaY);
	AA[1][0] = 0; AA[1][1] = 1; AA[1][2] = 0;
	AA[2][0] = -sin(sitaY); AA[2][1] = 0; AA[2][2] = cos(sitaY);
	MatrixAB(AA, Losc, Losd);
	//		printf("Losd:%f,%f,%f\n",Losd[0],Losd[1],Losd[2]);


	//1.四元数转欧拉角
	double yaw;
	double roll;
	double pitch;
	double x = areapara.Q1;
	double y = areapara.Q2;
	double z = areapara.Q3;
	double w = areapara.Q0;

	// roll (x-axis rotation)
	double sinr_cosp = +2.0 * (w * x + y * z);
	double cosr_cosp = +1.0 - 2.0 * (x * x + y * y);
	roll = atan2(sinr_cosp, cosr_cosp);

	// pitch (y-axis rotation)
	double sinp = +2.0 * (w * y - z * x);
	if (fabs(sinp) >= 1)
		pitch = copysign(M_PI / 2, sinp); // use 90 degrees if out of range
	else
		pitch = asin(sinp);


	// yaw (z-axis rotation)
	double siny_cosp = +2.0 * (w * z + x * y);
	double cosy_cosp = +1.0 - 2.0 * (y * y + z * z);
	yaw = atan2(siny_cosp, cosy_cosp);
	double a = 2;
	yaw = yaw + a;
	//    return yaw;

	//2.0 根据姿态角建立旋转矩阵



	double Mp[3][3], Mr[3][3], My[3][3];
	My[0][0] = cos(yaw); My[0][1] = -sin(yaw); My[0][2] = 0;
	My[1][0] = sin(yaw); My[1][1] = cos(yaw); My[1][2] = 0;
	My[2][0] = 0; My[2][1] = 0; My[2][2] = 1;

	Mr[0][0] = 1; Mr[0][1] = 0; Mr[0][2] = 0;
	Mr[1][0] = 0; Mr[1][1] = cos(roll); Mr[1][2] = -sin(roll);
	Mr[2][0] = 0; Mr[2][1] = sin(roll); Mr[2][2] = cos(roll);

	Mp[0][0] = cos(pitch); Mp[0][1] = 0; Mp[0][2] = sin(pitch);
	Mp[1][0] = 0; Mp[1][1] = 1; Mp[1][2] = 0;
	Mp[2][0] = -sin(pitch); Mp[2][1] = 0; Mp[2][2] = cos(pitch);

	//double Lose[3], Lostmp1[3][3], Lostmp2[3][3];
	//mult33(Mp, My, Lostmp1);
	//mult33(Lostmp1, Mr, Lostmp2);
	double Ae_roll[3][3] = { { 1, 0, 0 },{ 0, cos(roll), sin(roll) },{ 0, -sin(roll), cos(roll) } };
	double Ae_pitch[3][3] = { { cos(pitch), 0, -sin(pitch) },{ 0, 1, 0 },{ sin(pitch), 0, cos(pitch) } };
	double Ae_yaw[3][3] = { { cos(yaw), sin(yaw), 0 },{ -sin(yaw), cos(yaw), 0 },{ 0, 0, 1 } };
	double Aer[3][3] = { 0 };
	double AerTemp[3][3] = { 0 };
	mult33(Ae_yaw, Ae_pitch, AerTemp);
	mult33(AerTemp, Ae_roll, Aer);
	double AerInv[3][3] = { 0 };
	inv(Aer, AerInv);
	//MatrixAB(Lostmp2, Losd, Lose);
	double Lose[3], Lostmp1[3], Lostmp2[3];

	//MatrixAB(Mp, Losd, Lostmp1);
	//MatrixAB(My, Lostmp1, Lostmp2);
	//MatrixAB(Mr, Lostmp2, Lose);
	MatrixAB(AerInv, Losd, Lose);
	//MatrixAB(My, Lostmp1, Lostmp2);
	//MatrixAB(Mr, Lostmp2, Lose);
	//		printf("Lose:%f,%f,%f\n",Lose[0],Lose[1],Lose[2]);

	//
	double xyzo[6], xyz0temp1[3], xyzotemp2[3];
	xyz0temp1[0] = areapara.X; xyz0temp1[1] = areapara.Y; xyz0temp1[2] = areapara.Z;
	//xyzo[3] = areapara.Vx; xyzo[4] = areapara.Vy; xyzo[5] = areapara.Vz;
	//xyzo[0] = areapara.X; xyzo[1] = areapara.Y; xyzo[2] = areapara.Z;
	//xyzo[3] = areapara.Vx; xyzo[4] = areapara.Vy; xyzo[5] = areapara.Vz;
	/*xyzo[0] = TSatPosECF.x; xyzo[1] = TSatPosECF.y; xyzo[2] = TSatPosECF.z;*/
	xyzo[3] = TSatVelECF.x; xyzo[4] = TSatVelECF.y; xyzo[5] = TSatVelECF.z;
	//xyzg[0] = areapara.X; xyzg[1] = areapara.Y; xyzg[2] = areapara.Z;
	//xyzg[3] = areapara.Vx; xyzg[4] = areapara.Vy; xyzg[5] = areapara.Vz;

	getXYZGfroXYZO(xyz0temp1, epochYear,  epochDay, xyzotemp2, areapara.GMST);

	
	//getXYZOfroXYZG(xyzg, m_line_timedate, xyzo);
	//getXYZOfroXYZG(xyzg, m_line_timedate, xyzo);
	xyzo[0] = xyzotemp2[0];
	xyzo[1] = xyzotemp2[1];
	xyzo[2] = xyzotemp2[2];


	//		printf("xyzo:%f,%f,%f,%f,%f,%f\n",xyzo[0],xyzo[1],xyzo[2],xyzo[3],xyzo[4],xyzo[5]);


	//3.0 从姿态坐标系到惯性坐标系  
	double XX[3], YY[3], ZZ[3];
	double MXYZ[3][3], MXYZT[3][3];
	double rr; rr = sqrt(xyzo[0] * xyzo[0] + xyzo[1] * xyzo[1] + xyzo[2] * xyzo[2]);
	ZZ[0] = -xyzo[0] / rr;
	ZZ[1] = -xyzo[1] / rr;
	ZZ[2] = -xyzo[2] / rr;
	double tmp1, tmp2, tmp3, tmp;
	tmp1 = xyzo[4] * ZZ[2] - xyzo[5] * ZZ[1];
	tmp2 = xyzo[5] * ZZ[0] - xyzo[3] * ZZ[2];
	tmp3 = xyzo[3] * ZZ[1] - xyzo[4] * ZZ[0];
	tmp = sqrt(tmp1 * tmp1 + tmp2 * tmp2 + tmp3 * tmp3);
	YY[0] = -tmp1 / tmp;
	YY[1] = -tmp2 / tmp;
	YY[2] = -tmp3 / tmp;
	XX[0] = ZZ[1] * YY[2] - ZZ[2] * YY[1];
	XX[1] = ZZ[2] * YY[0] - ZZ[0] * YY[2];
	XX[2] = ZZ[0] * YY[1] - ZZ[1] * YY[0];

	MXYZ[0][0] = -XX[0]; MXYZ[0][1] = YY[0]; MXYZ[0][2] = ZZ[0];
	MXYZ[1][0] = -XX[1]; MXYZ[1][1] = YY[1]; MXYZ[1][2] = ZZ[1];
	MXYZ[2][0] = -XX[2]; MXYZ[2][1] = YY[2]; MXYZ[2][2] = ZZ[2];


	double Loso[3];
	MatrixAB(MXYZ, Lose, Loso);
	//		printf("Loso:%f,%f,%f\n",Loso[0],Loso[1],Loso[2]);

	//4.0 视点矢量和地球模型的交点
	double sate[3];
	sate[0] = xyzo[0]; sate[1] = xyzo[1]; sate[2] = xyzo[2];
	double P_xyzo[3];
	double HS;
	HS = 0;
	getXYZO(sate, Loso, HS, P_xyzo);
	//getXYZO(sate, Loso, m_earth, HS, P_xyzo);
	//		printf("P_xyzo:%f,%f,%f\n",P_xyzo[0],P_xyzo[1],P_xyzo[2]);

	/*double P_xyzg[3];
	POSITION j2000point;
	j2000point.x = P_xyzo[0];
	j2000point.y = P_xyzo[1];
	j2000point.z = P_xyzo[2];
	RJ2000toECF(areapara, j2000point);
	P_xyzg[0] =  j2000point.x ;
	P_xyzg[1] =  j2000point.y;
	P_xyzg[2] =  j2000point.z;
*/
//getXYZGfroXYZO(P_xyzo, m_line_timedate, P_xyzg);
//		printf("P_xyzg:%f,%f,%f\n",P_xyzg[0],P_xyzg[1],P_xyzg[2]);


//5.0 获得经纬度


	double hh = 0;

	getBL(P_xyzo, 0, earthpoint.lat, earthpoint.lon, hh);
	earthpoint.lat = earthpoint.lat * 180.0 / Pai;
	earthpoint.lon = earthpoint.lon * 180.0 / Pai;
	//cout << "oulalon：" << earthpoint.lon << '\n' << "oulalat：" << earthpoint.lat << endl;
	//printf("BL:%f,%f\n", LatLonH[0] * 180 / PI, LatLonH[1] * 180 / PI);


}




//////////////////////////////////////////////////////////////////////////////
void TLE2Point::PointFromTle(string str1, string str2, string str3, EARTHPOINT &earthpoint, double minutes)
{
   // Test SGP4
  

   cTle tleSGP4(str1, str2, str3);

   PrintPosVel(tleSGP4);

   //printf("\n");



   cTle tleSDP4(str1, str2, str3);

   PrintPosVel(tleSDP4);

   //printf("\nExample output:\n");

   // Example: Define a location on the earth, then determine the look-angle
   // to the SDP4 satellite defined above.

   // Create an Orbit object using the SDP4 TLE object.
   

   // Get the location of the satellite from the Orbit object. The 
   // earth-centered inertial information is placed into eciSDP4.
   // Here we ask for the location of the satellite 90 minutes after
   // the TLE epoch.
   cOrbit orbitSDP4(tleSDP4);
   cEciTime eciSDP4 = orbitSDP4.GetPosition(minutes);
   AREAPARA areapara;
   areapara.Vx = eciSDP4.Velocity().m_x;
   areapara.Vy = eciSDP4.Velocity().m_y;
   areapara.Vz = eciSDP4.Velocity().m_z;

   areapara.X = eciSDP4.Position().m_x;
   areapara.Y = eciSDP4.Position().m_y;
   areapara.Z = eciSDP4.Position().m_z;
   areapara.Q0 = 1;
   areapara.Q1 = 0;
   areapara.Q2 = 0;
   areapara.Q3 = 0;

   
   //Zeptomoby::OrbitTools::cOrbit m_Orbit;

   int    epochYear = (int)tleSGP4.GetField(cTle::FLD_EPOCHYEAR);
   double epochDay = tleSGP4.GetField(cTle::FLD_EPOCHDAY);
   cJulian jd = orbitSDP4.Epoch();

   areapara.GMST = jd.ToGmst();
   POSITION TSatPosECF, TSatVelECF;
   TSatPosECF = RJ2000toECF(areapara);
   TSatVelECF = VJ2000toECF(areapara);

   mygetpoint(epochYear,  epochDay, areapara
	   ,  TSatPosECF,  TSatVelECF,  earthpoint);


   earthpoint.hight = (sqrt(SQUARE(areapara.X) + SQUARE(areapara.Y) + SQUARE(areapara.Z)) - 6377.830)*1000;



   // Now create a site object. Site objects represent a location on the 
   // surface of the earth. Here we arbitrarily select a point on the
   // equator.
   
}

int main()
{

	string str1 = "SGP4 Test";
	string str2 = "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0     8";
	string str3 = "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518   105";
	EARTHPOINT  earthpoint;
	PointFromTle(str1, str2, str3, earthpoint,60);
	std::cout << "lat:" << earthpoint.lat << "*****lon:" << earthpoint.lon << "********hight:" << earthpoint.hight << endl;

}