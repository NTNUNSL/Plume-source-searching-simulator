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
	//cout<<"algorithm: greedy"<<endl;

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
	statis stat;

	// Algorithm parameters
	double z_direction = rand() % 360;
	double z_drift = argc>1?atoi(argv[1]):45; // degree
	double z_distance = argc>2?atof(argv[2]):1; // meters
	double z_trim = argc>3?atoi(argv[3]):15; // trim degree
	int max_flight_time = argc>4?atoi(argv[4]):1200; // flight time limit

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

	//cout<<"source: "<<xs<<","<<ys<<" wind_direction: "<<windDir<<endl;

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
		// Zigzag
		double trimAngle = 0;
		z_direction = windDir;
		Zigzag zigzag(z_direction,z_drift,z_distance);

		while(true){
			if(flight_time>=max_flight_time || lostCount >= 10 ){
				double err_dis = calcErrDis(log,xs,ys,err_dis_interval);
				success++;
				stat.length.push_back(length);
				stat.sample.push_back(sample);
				stat.err.push_back(err_dis);
				stat.time.push_back(flight_time);
				break;
			}

			// move
			offset = zigzag.genNextOffset();
			gdm.setCurLoc(gdm.getX()+offset.x,gdm.getY()+offset.y,gdm.getZ());
			length += sqrt(offset.x*offset.x+offset.y*offset.y);

			// measure PM2.5 and trim angle
			prePM25 = curPM25;
			curPM25 = gdm.evaluate(10);
			sample ++;

			// Log
			Data d;
			d.no = counter;
			d.x = gdm.getX();
			d.y = gdm.getY();
			d.direction = zigzag.getDirection();
			d.facing = zigzag.getDirectionDrifted();
			d.pm25 = curPM25;
			d.difference = curPM25 - prePM25;
			log.push_back(d);

			// Trim angle for next step
			if(counter>=1){
				if(zigzag.curStat()=='r')
				{
					if(d.difference>log[counter-1].difference){
						trimAngle = z_trim;
					}else{
						trimAngle = -z_trim;
					}
				}
				else if (zigzag.curStat() == 'l')
				{
					if(d.difference>log[counter-1].difference){
						trimAngle = -z_trim;
					}else{
						trimAngle = z_trim;
					}
				}
				zigzag.trim(trimAngle);
			}
			flight_time = sample*sample_duration+length/flight_speed;
			counter++;

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
        cout<<argv[1]<<","<<argv[2]<<","<<argv[3];
        cout<<","<<stat.sample[i]<<","<<stat.length[i]<<","<<stat.time[i]<<","<<stat.err[i]<<endl;
    }
	//cout<<"loop counts: "<<success<<"/"<<loop<<endl;
	//cout<<"avg sample: "<<stat.mean(stat.sample)<<" std_dev: "<<stat.std_dev(stat.sample)<<endl;
	//cout<<"avg length: "<<stat.mean(stat.length)<<" std_dev: "<<stat.std_dev(stat.length)<<endl;
	//cout<<"avg time: "<<stat.mean(stat.time)<<" std_dev: "<<stat.std_dev(stat.time)<<endl;
	//cout<<"avg error distance: "<<stat.mean(stat.err)<<" std_dev: "<<stat.std_dev(stat.err)<<endl;

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

