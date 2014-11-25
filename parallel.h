#ifndef PARALLEL_H
#define PARALLEL_H

#include "main.h"

enum class ParallelType{
    // OpenMP
    PARALLEL_FOR,   // parallel for with critical section
    REDUCTION,      // parallel for with reduction
    INDICES,        // parallel for with calculating indices of the array
    
    // pthreads
    PTHREADS,

    // C++11 multithreading
    CPP11
};

// structure for exchanging data in pthreads version
typedef struct t_pdata
{
    double* array;
    double* array_b;
    int length;
    int index;
    double global;
    double mean_a;
    double mean_b;
} t_pdata;

// pthreads worker functions
void* mean_ph_worker(void* d);
void* stddev_ph_worker(void* d);
void* pearson_ph_worker(void* d);

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
double dispatcher_ph(
    double* p_array, 
    int what, 
    double p_mean_a, 
    double p_mean_b,
    double p_stddev_a,
    double p_stddev_b,
    double* p_array_b
);

// C++11 Multithreading
double mean_cpp11(double* p_array, int p_index, int p_length, double* global);
double stddev_cpp11(double* p_array, int p_index, int p_length, double p_mean, double* global);
double pearson_cpp11(double* p_array_a, double* p_array_b, int p_index, int p_length, double p_mean_a, double p_mean_b, double* global);
double dispatcher_cpp11(
    double* p_array, 
    int what, 
    double p_mean_a, 
    double p_mean_b,
    double p_stddev_a,
    double p_stddev_b,
    double* p_array_b
);


void set_threads(int p_threads);

// to start calculations
double run_parallel_pearson(ParallelType p_type, int p_threads);

int threads;

};

#endif