//
// stdafx.h 
//

#include <iostream>
#define WIN32_LEAN_AND_MEAN   // Exclude rarely-used stuff from Windows headers
#include <stdio.h>
#include <tchar.h>
#include <time.h>
#include <Math.h>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>

using namespace std;
namespace TLE2Point {
#define  WE         (7.2921151467E-5)
#define SQUARE(x) ((x)*(x))
#define M_PI       3.14159265358979323846   // pi
	const double Pai = 3.1415926535897932384626433832795;
	const double Re = 6378137.0;
	const double Rp = 6356752.314245;//6356800.0;
	const int gps_sample = 7;
	const int quar_sample = 5;
	const double DEG = 0.01745329251994;	//1�ȵĻ���
	const double R = 6378137;				//����볤��
	const double RS = 6356752.314245215709;	//�����̾�
	const double ECCENT = 0.00669437999013;	//������ƽ��e^2
	const double e1 = sqrt((6378137.0 * 6378137.0 - 6356752.31 * 6356752.31) / (6378137.0 * 6378137.0));
	const double e2 = sqrt((6378137.0 * 6378137.0 - 6356752.31 * 6356752.31) / (6356752.31 * 6356752.31));
	const double VOL = 300000000;			//����
	typedef struct _areapara {
		double GMST;			// GreenwichSiderealTime
		double X;				// λ��,��λm
		double Y;
		double Z;
		double Vx;				// �ٶ�,��λm/sec
		double Vy;
		double Vz;
		double Q0;				// ��̬��Ԫ��
		double Q1;
		double Q2;
		double Q3;
		int Num;				// ���Ƿ�Χ��������
		int Resolution;			// ��Ʒֱ��ʣ���λm
		double theta;			// ��
	}AREAPARA;

	//��ά����ṹ�壨������λ�á��ٶȣ�
	typedef struct Position
	{
		double x;
		double y;
		double z;
		bool flag;//λ����Ч��־
	}POSITION;
	typedef struct _earthpoint {
		double lat;
		double lon;
		double hight;
	}EARTHPOINT;
	typedef struct _TIMEDATE
	{
		int year;
		int month;
		int day;
		int hour;
		int minute;
		double second;
	}TIMEDATE;

#ifdef __cplusplus
	extern "C"
#endif
	{

		__declspec(dllexport) POSITION RJ2000toECF(AREAPARA areapara);
		__declspec(dllexport) POSITION VJ2000toECF(AREAPARA areapara);
		__declspec(dllexport) void PointFromTle(string str1, string str2, string str3, EARTHPOINT& earthpoint,  double minutes);
#ifdef __cplusplus
	}
#endif



}

