#ifndef GAUSSIAN_DISPERSION_MODEL
#define GAUSSIAN_DISPERSION_MODEL
#include <cstdlib>
#include <stdio.h>
#include <sstream>
#include <random>
#include <chrono>
#include <iostream>
#include <math.h>
#define PI 3.1415926 
#include "statistics.hpp"

// 高斯煙流模型
class GaussianDispersionModel{
private:
	double _x,_y,_z;
	double _Q,_u,_dir;
	double _xs,_ys,_H;
	int _st;
	double _bg;
	bool _noise;
	double _stddev;
	unsigned _seed;
	std::default_random_engine _gen;

public:
	GaussianDispersionModel();
    
    // 設定煙流源頭座標
    void setSource(double xs,double ys,double H);
	
    // 設定無人機座標
    void setCurLoc(double x,double y,double z);
	
    // 設定風向,風速,穩定係數
    void setWind(double u,double dir,int st);
	
    // 設定排放率
    void setRate(double Q);
	
    // 設定背景濃度值
    void setBackGround(double n);
	
    // 將背景加入雜訊
    void setNoise(bool bg);
	
    // 設定雜訊標準差
    void setStddev(double std);
	
    // 取得無人機座標
    double getX();
	double getY();
	double getZ();

    // (棄用)使用python產生煙流模型濃度
	double evaluatePy();
	
    // 使用C產生煙流模型濃度
    double evaluate();
	double evaluate(int times);
	
    // 計算風力擴散係數
    void calc_sigmas(int cate,double x1,double &sig_y,double &sig_z);

};

#endif

