#include "search.hpp"
Data newData(int no,double x, double y,double direction,double facing,double pm25,double difference){
	Data d;
	d.x = x;
	d.y = y;
	d.no = no;
	d.direction = direction;
	d.facing = facing;
	d.pm25 =pm25;
	d.difference =difference;
	return d;
}

// 藉由飛行紀錄計算最大濃度與排放源距離差
double calcErrDis(std::vector<Data> &log,const int xs,const int ys,const int interval){
	double max = 0;
	int max_ind = 0;
	for(int i=0;i<log.size();i++){
		if(log[i].pm25>max){
			max = log[i].pm25;
			max_ind = i;		
		}
	}
	double err_dis = sqrt((log[max_ind].x-xs)*(log[max_ind].x-xs)+(log[max_ind].y-ys)*(log[max_ind].y-ys));
	return err_dis;
	
}

// 取得最近windowSize次觀測中的最高濃度
double windowHighest(std::vector<Data> &log,const int windowSize){
	if(log.size()==0){return -1;}
	double highest = log[log.size()-1].pm25;
	for(int i=0;i<windowSize;i++){
		if(i>log.size()){
			break;
		}else if(log[log.size()-1-i].pm25>highest){
			highest = log[log.size()-1-i].pm25;
		}
	}
	return highest;
}
