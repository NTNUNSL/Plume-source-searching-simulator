#include "gdm.hpp"
#include "zigzag.hpp"
#include "search.hpp"
#include <vector>
#include <unistd.h>
#include <fstream>
#include <stdlib.h>
#include <time.h>

using namespace std;

int main(int argc,char *argv[])
{
	srand(time(0));
	//cout<<"algorithm: space-filling"<<endl;

	// Program parameters
	int loop = 1000;
	bool saveLog = true;
	int counter = 0;
	int success = 0;
	vector<Data> log;
	double prePM25 = 0;
	double curPM25 = 0;
	double err_dis_interval = 2;
	int lostCount = 0;
	double maxPM25 = 0;
	double flight_speed = 10;
	double sample_duration = 10;
	double flight_height = 0;
	int max_flight_time = 30;
	statis stat;

	// Algorithm parameters
	double spacing = 1; // 行進間格距離
	double area = 3600;
	double side = sqrt(area); // 邊長

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

	for(int t=0;t<loop;t++)
	{
		// Initialize
		xs = -30 + (rand()%10);
		ys = 10 + (-1)*(rand()%20);
		Offset offset;
		int sample = 0;
		double length = 0;
		double flight_time = 0;
		log.clear();
		prePM25 = 0;
		curPM25 = 0;
		counter = 0;
		lostCount = 0;
		maxPM25 = 0;
		// GDM
		GaussianDispersionModel gdm;
		gdm.setSource(xs,ys,H);
		gdm.setCurLoc(0,0,flight_height);
		gdm.setWind(u,windDir,stableLevel);
		gdm.setRate(Q);
		gdm.setNoise(noise);
		gdm.setStddev(stddev);
		gdm.setBackGround(background);

		while(true){
			if(flight_time>=1200 || lostCount >= 10 ){
				double err_dis = calcErrDis(log,xs,ys,err_dis_interval);
				stat.length.push_back(length);
				stat.sample.push_back(sample);
				stat.err.push_back(err_dis);
				stat.time.push_back(flight_time);
				break;
			}

			int y = 0;
			// move
			for(int x=0;x>-side;){
				for(;y>(-1)*side/2;y-=spacing){
					gdm.setCurLoc(x,y,gdm.getZ());
					prePM25 = curPM25;
					curPM25 = gdm.evaluate(10);
					log.push_back(newData(counter,x,y,0,0,curPM25,curPM25-prePM25));
					length += spacing;
					sample++;
				}
                y+=spacing;
				x-=spacing;
				length+=spacing;
				gdm.setCurLoc(x,y,gdm.getZ());
				prePM25 = curPM25;
				curPM25 = gdm.evaluate(10);
				log.push_back(newData(counter,x,y,0,0,curPM25,curPM25-prePM25));
				for(;y<side/2;y+=spacing){
					gdm.setCurLoc(x,y,gdm.getZ());
					prePM25 = curPM25;
					curPM25 = gdm.evaluate(10);
					log.push_back(newData(counter,x,y,0,0,curPM25,curPM25-prePM25));
					length += spacing;
					sample++;
				}
                y-=spacing;
				x-=spacing;
				length+=spacing;
				gdm.setCurLoc(x,y,gdm.getZ());
				prePM25 = curPM25;
				curPM25 = gdm.evaluate(10);
				log.push_back(newData(counter,x,y,0,0,curPM25,curPM25-prePM25));
				for(;y>0;y-=spacing){
					gdm.setCurLoc(x,y,gdm.getZ());
					prePM25 = curPM25;
					curPM25 = gdm.evaluate(10);
					log.push_back(newData(counter,x,y,0,0,curPM25,curPM25-prePM25));
					length += spacing;
					sample++;
				}
			}
			

			flight_time = sample*sample_duration+length/flight_speed;

			if(curPM25>maxPM25){
				maxPM25 = curPM25;
				lostCount = 0;
			}else if(curPM25>0.1*maxPM25){
				lostCount = 0;
			}else{
				lostCount++;
			}
		}
	}


	// Show result
    for(int i=0;i<stat.sample.size();i++){
        cout<<stat.sample[i]<<","<<stat.length[i]<<","<<stat.time[i]<<","<<stat.err[i]<<endl;
    }

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

