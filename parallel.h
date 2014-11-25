#ifndef PARALLEL_H
#define PARALLEL_H

#include "main.h"

enum class ParallelType{
    PARALLEL_FOR,
    REDUCTION,
    INDICES,
    PTHREADS
};

// class PData
// {
// public:
//     PData(double* p_array, int p_index, int p_length)
//     : array(p_array), index(p_index), length(p_length)
//     {
//         global = 0;
//     };

//     PData(){};

//     double* array;
//     int index;
//     int length;
//     double global;
// };

typedef struct t_pdata
{
    double* array;
    int length;
    int index;
    double global;
} t_pdata;


class Parallel
{
private:
    Data& data;
    Stats stats;
public:
    Parallel(Data& p_data);

// Parallel for
double mean_for(double* p_array);
double stddev_for(double* p_array, double p_mean);
double pearson_for(
    double* p_array_a, 
    double* p_array_b, 
    double p_mean_a, 
    double p_mean_b, 
    double p_std_dev_a, 
    double p_std_dev_b
);

// Reduction
double mean_red(double* p_array);
double stddev_red(double* p_array, double p_mean);
double pearson_red(
    double* p_array_a, 
    double* p_array_b, 
    double p_mean_a, 
    double p_mean_b, 
    double p_std_dev_a, 
    double p_std_dev_b
);

// Indices
double mean_ind(double* p_array);
double stddev_ind(double* p_array, double p_mean);
double pearson_ind(
    double* p_array_a, 
    double* p_array_b, 
    double p_mean_a, 
    double p_mean_b, 
    double p_std_dev_a, 
    double p_std_dev_b
);

// Indices with pthreads
double mean_ph(double* p_array);
double stddev_ph(double* p_array, double p_mean);
double pearson_ph(
    double* p_array_a, 
    double* p_array_b, 
    double p_mean_a, 
    double p_mean_b, 
    double p_std_dev_a, 
    double p_std_dev_b
);

void set_threads(int p_threads);

// to start calculations
double run_parallel_pearson(ParallelType p_type, int p_threads);

int threads;

};


void* mean_ph_worker(void* d);




#endif