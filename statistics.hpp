#ifndef SEARCH_STATISTICS_HPP
#define SEARCH_STATISTICS_HPP
#include<cmath>
#include<vector>
using namespace std;

// 統計模型, 方便計算飛行結果
class statis{
public:
	vector<double> length;
	vector<double> sample;
	vector<double> time;
	vector<double> err;

	double sum(vector<double> data);
	double mean(vector<double> data);
	double variance(vector<double> data);
	double std_dev(vector<double> data);	
	double get95confidence(vector<double> data);
	    
};


#endif
