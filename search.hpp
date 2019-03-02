#ifndef PLUME_SEARCH_HPP
#define PLUME_SEARCH_HPP

#include<cmath>
#include<vector>
#include<stdlib.h>
#define DEG2RAD 0.0174532925

struct Data
{
	int no;
	double x;
	double y;
	double direction;
	double facing;
	double pm25;
	double difference;
};

// 建立一筆飛行量測資料
Data newData(int no,double x,double y,double direction,double facing,double pm25,double difference);

// 排放源觀測誤差
double calcErrDis(std::vector<Data> &log,const int xs,const int ys,const int interval);

// 最近觀測的最高濃度
double windowHighest(std::vector<Data> &log,const int windowSize);

#endif
