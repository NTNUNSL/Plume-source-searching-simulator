#include "gdm.hpp"
GaussianDispersionModel::GaussianDispersionModel(){
	_x=100;
	_y=0;
	_z=70;
	_Q=40;
	_u=5;
	_dir=270;
	_xs=0;
	_ys=0;
	_H=50;
	_st=1;
	_bg=0;
	_noise=false;
	_stddev=0.1;
	_seed = std::chrono::system_clock::now().time_since_epoch().count();
	_gen = std::default_random_engine(_seed);
}

// 利用風場穩定度與下風處位置計算擴散係數,使用查表法
void GaussianDispersionModel::calc_sigmas(int cate,double x1,double &sig_y,double &sig_z){
    double x = abs(x1);
    double a,b,c,d;
    if(cate==1)
    {
        if(x<100 && x>0){a=122.8;b=0.9447;}
        else if(x>=100 && x<150){a=158.08;b=1.0542;}
        else if(x>=150 && x<200){a=170.22;b=1.0932;}
        else if(x>=200 && x<250){a=179.52;b=1.1262;}
        else if(x>=250 && x<300){a=217.410;b=1.2644;}
        else if(x>=300 && x<400){a=258.89;b=1.4094;}
        else if(x>=400 && x<500){a=346.75;b=1.7283;}
        else if(x>=500 && x<3110){a=453.85;b=2.1166;}
        else if(x>=3110){a=453.85;b=2.1166;}
        c=24.167;
        d=2.5334;
    }
    else if(cate==2)
    {
        if(x<200 && x>0){a=90.673;b=0.93198;}
        else if(x>=200 && x<400){a=98.483;b=0.98332;}
        else if(x>=400){a=109.3;b=1.0971;}
        c=18.333;
        d=1.8096;
    }
    else if(cate==3)
    {
        a=61.141;
        b=0.91465;
        c=12.5;
        d=1.0857;
    }
    else if(cate==4)
    {
        if(x<300 && x>0){a=34.459;b=0.86974;}
        else if(x>=300 && x<1000){a=32.093;b=0.81066;}
        else if(x>=1000 && x<3000){a=32.093;b=0.64403;}
        else if(x>=3000 && x<10000){a=33.504;b=0.60486;}
        else if(x>=10000 && x<30000){a=36.650;b=0.56589;}
        else if(x>=30000){a=44.053;b=0.51179;}
        c=8.333;
        d=0.72382;
    }
    else if(cate==5)
    {
        if(x<100 && x>0){a=24.26;b=0.8366;}
        else if(x>=100 && x<300){a=23.331;b=0.81956;}
        else if(x>=300 && x<1000){a=21.628;b=0.7566;}
        else if(x>=1000 && x<2000){a=21.628;b=0.63077;}
        else if(x>=2000 && x<4000){a=22.534;b=0.57154;}
        else if(x>=4000 && x<10000){a=24.703;b=0.50527;}
        else if(x>=10000 && x<20000){a=26.970;b=0.46713;}
        else if(x>=20000 && x<40000){a=35.420;b=0.37615;}
        else if(x>=40000){a=47.618;b=0.29592;}
        c=6.25;
        d=0.54287;
    }
    else if(cate==6)
    {
        if(x<200 && x>0){a=15.209;b=0.81158;}
        else if(x>=200 && x<700){a=14.457;b=0.78407;}
        else if(x>=700 && x<1000){a=13.953;b=0.68465;}
        else if(x>=1000 && x<2000){a=13.953;b=0.63227;}
        else if(x>=2000 && x<3000){a=14.823;b=0.54503;}
        else if(x>=3000 && x<7000){a=16.187;b=0.46490;}
        else if(x>=7000 && x<15000){a=17.836;b=0.41507;}
        else if(x>=15000 && x<30000){a=22.651;b=0.32681;}
        else if(x>=30000 && x<60000){a=27.074;b=0.27436;}
        else if(x>=60000){a=34.219;b=0.21716;}
        c=4.1667;
        d=0.36191;
    }
    sig_z = a*pow(x/1000,b);
    if(sig_z>5000){sig_z=5000;}

    double theta=0.017453293*(c-d*log(abs(x+1e-15)/1000));
    sig_y=465.11628*x/1000*tan(theta);
}

// 多次量測平均值
double GaussianDispersionModel::evaluate(int times){
    double s = 0;
    for(int i =0;i<times;i++){
        s = s + evaluate();
    }
    s = s/times;
    return s;
}

// 實作高斯煙流模型並計算無人機座標的濃度
double GaussianDispersionModel::evaluate(){
	double Q=_Q;
	double u=_u;
	double dir1=_dir;
	double x=_x;
	double y=_y;
	double z=_z;
	double xs=_xs;
	double ys=_ys;
	double H=_H;
	int st=_st;
    
	double u1 = u;
    double x1 = x-xs;
    double y1 = y-ys;

    double wx = u1*sin((dir1-180)*PI/180);
    double wy = u1*cos((dir1-180)*PI/180);
    double dot_product=wx*x1+wy*y1;
    double magnitudes=u1*sqrt(x1*x1+y1*y1); 
    double subtended=acos(dot_product/(magnitudes+1e-15));
    double hypotenuse=sqrt(x1*x1+y1*y1);
    double downwind=cos(subtended)*hypotenuse;
    double crosswind=sin(subtended)*hypotenuse;

    double sig_y,sig_z;
    calc_sigmas(st,downwind,sig_y,sig_z);

    double C=Q/(2*PI*u1*sig_y*sig_z)
        * exp(-(crosswind*crosswind)/(2*sig_y*sig_y))
        *(exp(-(z-H)*(z-H)/(2*sig_z*sig_z)) +
        exp(-(z+H)*(z+H)/(2*sig_z*sig_z)) )*1e6;

	if(_noise)
	{
        double stddevLow = 0.05;
        double stddevHigh = 0.1;
        double stddev = 0;
        double maxDist = 100;
        double distance = sqrt((x-xs)*(x-xs)+(y-ys)*(y-ys));

        if(distance <= maxDist){
            stddev = stddevHigh - (distance/maxDist)*(stddevHigh-stddevLow);
        }else{
            stddev = stddevLow;
        }
		std::normal_distribution<double> dist(0, C*stddev);
		C+=dist(_gen);
	}

    if(_bg!=0){
		C+=_bg;
	}

    if(_noise)
    {
		std::normal_distribution<double> dist(0, C*_stddev);
		C+=dist(_gen);
    }

	if(C<0)
	{
		C==0;
	}

    return C;
}

// (棄用)使用現成的python程式產生高斯煙流模型並計算無人機座標濃度
double GaussianDispersionModel::evaluatePy(){
	double PM25_concentration;
	std::ostringstream cmd;
	cmd<<"python3 gdm.py "<<_Q<<' '<<_u<<' '<<_dir<<' '<<_x<<' '<<_y<<' '<<_z<<' '<<_xs<<' '<<_ys<<' '<<_H<<' '<<_st;
	std::string tmp = cmd.str();
	//std::cout<<tmp.c_str()<<std::endl;

	FILE *fp = popen(tmp.c_str(),"r");
	fscanf(fp,"%lf",&PM25_concentration);
	pclose(fp);

	if(_bg!=0){
		PM25_concentration+=_bg;
	}

	if(_noise)
	{
		std::normal_distribution<double> dist(0, PM25_concentration*_stddev);
		PM25_concentration+=dist(_gen);
	}

	if(PM25_concentration<0)
	{
		PM25_concentration==0;
	}

	return PM25_concentration;
}
void GaussianDispersionModel::setSource(double xs,double ys,double H){
	_xs = xs;
	_ys = ys;
	_H = H;
}
void GaussianDispersionModel::setCurLoc(double x,double y,double z){
	_x = x;
	_y = y;
	_z = z;
}
void GaussianDispersionModel::setWind(double u,double dir,int st){
	_u = u;
	_dir = dir;
	_st = st;
}
void GaussianDispersionModel::setRate(double Q){
	_Q = Q;
}
void GaussianDispersionModel::setBackGround(double bg)
{
	_bg = bg;
}
void GaussianDispersionModel::setNoise(bool n)
{
	_noise = n;
}
void GaussianDispersionModel::setStddev(double std)
{
	_stddev = std;
}
double GaussianDispersionModel::getX(){return _x;}
double GaussianDispersionModel::getY(){return _y;}
double GaussianDispersionModel::getZ(){return _z;}
