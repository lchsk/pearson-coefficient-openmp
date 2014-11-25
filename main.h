#ifndef MAIN_H
#define MAIN_H

#include <ctime>

class Config
{
public:
    Config();

    int input_length;


};

class Stats
{
    public:
        double mean_time;
        double stddev_time;
        double pearson_time;
        double overall_time;

        double compute_diff(double time1, double time2);
        double compute_diff(clock_t time1, clock_t time2);
        void compute_overall_time();
};

class Data
{
public:
    double* a;
    double* b;
    Config& conf;

    Data(Config& p_c);
    ~Data();

private:
    
};

#endif