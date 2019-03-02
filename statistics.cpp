#include "statistics.hpp"
double statis::sum(vector<double> data){
	double sum=0;
	for(int i=0;i<data.size();i++){
		sum += data[i];
	}
	return sum;
}
double statis::mean(vector<double> data){
	return statis::sum(data)/data.size();
}
double statis::variance(vector<double> data){
	double mean = statis::mean(data);
	double vari_sum = 0;
	for(int i=0;i<data.size();i++){
		vari_sum += (data[i]-mean)*(data[i]-mean);
	}
	return vari_sum/data.size();

}
double statis::std_dev(vector<double> data){
	return sqrt(variance(data));
}

