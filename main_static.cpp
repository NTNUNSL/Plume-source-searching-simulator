#include "gdm.hpp"
#include "zigzag.hpp"
#include "search.hpp"
#include <vector>
#include <unistd.h>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#define DEG2RAD 0.0174532925

using namespace std;

// 無動態調整的靜態飛行法
int main(int argc,char *argv[])
{
	srand(time(0));
	//cout<<"algorithm: static"<<endl;

	// Program parameters
	int loop = 1000;
	bool saveLog = true;
	int counter = 0;
	int success = 0;
	vector<Data> log;
	double prePM25 = 0;
	double curPM25 = 0;
	int err_dis_interval = 2;
	int lostCount = 0;
	double maxPM25 = 0;	
	double flight_speed = 10;
	double sample_duration = 10;
	double flight_height = 0;
	statis stat;

	// GDM parameters
	int scale = 100;
	int xs = 0;
	int ys = 0;
	int H = 0;
	double Q = 0.1;
	int u = 5;
	double windDir = 270;
	int stableLevel = 3;
	double stddev = 0.01;
	bool noise = true;
	double background = 20;

	// Algorithm parameters
	double direction = rand() % 360;
	double windDis = 5;
	int far_c = argc>1?atoi(argv[1]):20;
	int close_c = argc>2?atoi(argv[2]):200;
	int area = argc>3?atoi(argv[3]):3600;
	double r_m = argc>4?atoi(argv[4]):10;
	double r_s = argc>5?atoi(argv[5]):3;
	double time_limit = argc>6?atoi(argv[6]):1200;

	double far_distance = sqrt(area)/r_m;
	double close_distance = r_s;
	double far_step = far_distance/r_s;
	double close_step = 1;
	//double beta = argc>7?atof(argv[7]):0.48;
	int windowSize = 5;


	for (int t = 0; t < loop; t++)
	{
		// Initialize
		xs = -30 + (rand()%10);
		ys = 10 + (-1)*(rand()%20);
		Offset offset;
		int sample = 0;
		double length = 0;
		log.clear();
		prePM25 = 0;
		curPM25 = 0;
		counter = 0;
		lostCount = 0;
		maxPM25 = 0;
		double flight_time = 0;
		double temp_max_pm25 = 0;
		// GDM
		GaussianDispersionModel gdm;
		gdm.setSource(xs, ys, H);
		gdm.setCurLoc(0, 0, flight_height);
		gdm.setWind(u, windDir, stableLevel);
		gdm.setRate(Q);
		gdm.setNoise(noise);
		gdm.setStddev(stddev);
		gdm.setBackGround(background);
		// Dynamic algorithm
		direction = rand() % 360;

		// Start algorithm
		// Evaluate wind
		direction = windDir;
		// Choose better start position
		double dx1 = windDis * sin((direction+90) * DEG2RAD);
		double dy1 = windDis * cos((direction+90) * DEG2RAD);
		double dx2 = windDis * sin((direction-90) * DEG2RAD);
		double dy2 = windDis * cos((direction-90) * DEG2RAD);
		double pmRight,pmLeft;
		curPM25 = gdm.evaluate(10);
		log.push_back(newData(counter,0,0,direction,direction,curPM25,0));
		gdm.setCurLoc(dx1, dy1, H);
        length+=windDis;
		pmRight = gdm.evaluate(10);
		log.push_back(newData(++counter,dx1,dy1,direction,direction,pmRight,0));
		gdm.setCurLoc(dy2, dy2, H);
        length+=windDis*2;
		pmLeft = gdm.evaluate(10);
		log.push_back(newData(++counter,dx2,dy2,direction,direction,pmLeft,0));
		if(curPM25>pmRight && curPM25>pmLeft){
			gdm.setCurLoc(0, 0, H);
            length+=windDis;
			log.push_back(newData(++counter,0,0,direction,direction,curPM25,0));
		}else if(pmRight>pmLeft){
			gdm.setCurLoc(dx1, dy1, H);
            length+=windDis*2;
			log.push_back(newData(++counter,dx1,dy1,direction,direction,pmRight,0));
		}else{
			gdm.setCurLoc(dx2, dy2, H);
			log.push_back(newData(++counter,dx2,dy2,direction,direction,pmLeft,0));
		}
		sample+=4;
		// Start Looping
		while (true)
		{
			if(lostCount>=10 || flight_time>time_limit )
			{
				double err_dis = calcErrDis(log,xs,ys,err_dis_interval);
				stat.sample.push_back(sample);
				stat.length.push_back(length);
				stat.err.push_back(err_dis);
				stat.time.push_back(flight_time);
				break;
			}

			double originX = gdm.getX();
			double originY = gdm.getY();
			prePM25 = curPM25;
			curPM25 = gdm.evaluate(10);
			sample ++;

			// Log
			log.push_back(newData(counter,originX,originY,direction,direction,curPM25,curPM25-prePM25));
			double highest = windowHighest(log,windowSize);

			// Adjust distances
			//double ratio = abs((highest-far_c)/(close_c-far_c));
			//if(ratio>1){ratio=1;}
			//ratio = (std::exp(ratio*beta)-1)/(std::exp(1*beta)-1);  // Do exponential scale
			double distance = far_step;// - ratio*(far_step-close_step);
			double distanceLarge = far_distance;// - ratio*(far_distance-close_distance);


			// Sample and determin next direction
			double x1, x2, x3, y1, y2, y3, pm1, pm2, pm3;
			x1 = distance * sin(direction * DEG2RAD);
			y1 = distance * cos(direction * DEG2RAD);
			x2 = distance * sin((direction + 90) * DEG2RAD);
			y2 = distance * cos((direction + 90) * DEG2RAD);
			x3 = distance * sin((direction + 270) * DEG2RAD);
			y3 = distance * cos((direction + 270) * DEG2RAD);
			gdm.setCurLoc(originX + x1, originY + y1, gdm.getZ());
			length += sqrt(x1*x1+y1*y1);
			pm1 = gdm.evaluate(10) - curPM25;
			sample++;
			log.push_back(newData(counter,originX+x1,originY+y1,direction,direction,curPM25+pm1,pm1));
			gdm.setCurLoc(originX + x2, originY + y2, gdm.getZ());
			length += sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
			pm2 = gdm.evaluate(10) - curPM25;
			sample++;
			log.push_back(newData(counter,originX+x2,originY+y2,direction,direction,curPM25+pm2,pm2));
			gdm.setCurLoc(originX + x3, originY + y3, gdm.getZ());
			length += sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2));
			pm3 = gdm.evaluate(10) - curPM25;
			sample++;
			log.push_back(newData(counter,originX+x3,originY+y3,direction,direction,curPM25+pm3,pm3));
			if (pm1 > pm2 && pm2 > pm3)
			{
				direction += 90 * pm2 / (pm1 + pm2);
			}
			else if (pm1 > pm3 && pm3 > pm2)
			{
				direction -= 90 * pm3 / (pm1 + pm3);
			}
			else if (pm2 > pm1 && pm1 > pm3)
			{
				direction += 90;
			}
			else if (pm2 > pm3 && pm3 > pm1)
			{
				direction += 120;
			}
			else if (pm3 > pm1 && pm1 > pm2)
			{
				direction -= 90;
			}
			else if (pm3 > pm2 && pm2 > pm1)
			{
				direction -= 120;
			}
			direction = fmod(direction, 360);

			// Large move
			double goX = distanceLarge * sin(direction * DEG2RAD);
			double goY = distanceLarge * cos(direction * DEG2RAD);
			gdm.setCurLoc(originX + goX, originY + goY, gdm.getZ());
			length += sqrt(goX*goX+goY*goY);
			counter++;
			flight_time = sample*sample_duration+length/flight_speed;
			temp_max_pm25 = curPM25;
			temp_max_pm25 = pm1>temp_max_pm25?pm1:temp_max_pm25;
			temp_max_pm25 = pm2>temp_max_pm25?pm2:temp_max_pm25;
			temp_max_pm25 = pm3>temp_max_pm25?pm3:temp_max_pm25;
			if(temp_max_pm25 > maxPM25 ){
				maxPM25 = temp_max_pm25;
				lostCount = 0;
			}else if(temp_max_pm25 > 0.1*maxPM25){
				lostCount = 0;
			}else{
				lostCount++;
			}
		}
	}

	// Show result
	/*
	cout<<"loop counts: "<<stat.sample.size()<<"/"<<loop<<endl;
	cout<<"avg sample: "<<stat.mean(stat.sample)<<" std_dev: "<<stat.std_dev(stat.sample)<<endl;
	cout<<"avg length: "<<stat.mean(stat.length)<<" std_dev: "<<stat.std_dev(stat.length)<<endl;
	cout<<"avg time: "<<stat.mean(stat.time)<<" std_dev: "<<stat.std_dev(stat.time)<<endl;
	cout<<"avg error distance: "<<stat.mean(stat.err)<<" std_dev: "<<stat.std_dev(stat.err)<<endl;
	*/
	
    for(int i=0;i<stat.sample.size();i++){    
        cout<<argv[1]<<","<<argv[2]<<","<<"static"<<","<<argv[4]<<","<<argv[5]<<","<<argv[3]<<","<<argv[6];
        cout<<","<<stat.sample[i]<<","<<stat.length[i]<<","<<stat.time[i]<<","<<stat.err[i]<<endl;
    }
	//cout<<","<<stat.mean(stat.sample)<<","<<stat.std_dev(stat.sample)<<","<<stat.mean(stat.length)<<","<<stat.std_dev(stat.length);
	//cout<<","<<stat.mean(stat.time)<<","<<stat.std_dev(stat.time)<<","<<stat.mean(stat.err)<<","<<stat.std_dev(stat.err)<<endl;

	if(saveLog)
	{
		ofstream f;
		f.open("log.csv");
		f<<Q<<','<<u<<','<<windDir<<','<<xs<<','<<ys<<','<<H<<','<<stableLevel<<"\n";
		for(int i=0;i<log.size();i++){
			f<<log[i].x<<','<<log[i].y<<','<<log[i].direction<<','<<log[i].facing<<',';
			f<<log[i].pm25<<','<<log[i].difference<<"\n";
		}
		f.close();
	}

	return 0;
}
