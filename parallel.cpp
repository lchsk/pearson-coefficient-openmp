#include "definitions.h"
#include "parallel.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "main.h"

Parallel::Parallel(Data& p_data) : data(p_data)
{
    threads = 2;
}

void Parallel::set_threads(int p_threads)
{
    threads = p_threads;
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
}

double 
Parallel::mean_red(double* p_array)
{
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < data.conf.input_length; ++i)
    {
        sum += p_array[i];
    }

    return sum / data.conf.input_length;
}

double 
Parallel::stddev_red(double* p_array, double p_mean)
{
    double nominator = 0;

    #pragma omp parallel for reduction(+:nominator)
        for (int i = 0; i < data.conf.input_length; ++i)
        {
            nominator += (p_array[i] - p_mean) * (p_array[i] - p_mean);
        }

    return sqrt(nominator / data.conf.input_length);
}

double 
Parallel::pearson_red(
    double* p_array_a, 
    double* p_array_b, 
    double p_mean_a, 
    double p_mean_b, 
    double p_std_dev_a, 
    double p_std_dev_b
)
{
    double sum = 0;

    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < data.conf.input_length; ++i)
    {
        sum += ((p_array_a[i] - p_mean_a) * (p_array_b[i] - p_mean_b));
    }

    return (sum / data.conf.input_length) / (p_std_dev_a * p_std_dev_b);
}


double 
Parallel::mean_ind(double* p_array)
{
    double global = 0;

    #pragma omp parallel shared(global)
    {
        int i = omp_get_thread_num() * (data.conf.input_length / omp_get_num_threads());
        int per = data.conf.input_length / omp_get_num_threads();
        
        //printf("Thread %d here, start=%d, per=%d\n", omp_get_thread_num(), i, per);

        int c;
        double local = 0;

        for (i, c = 0; c < per; c++, i++)
        {
            local += p_array[i];
            //if (c > 1 && c % 250000 == 0)
            //    printf("T %d in the middle, c=%d, value: %.10f\n", omp_get_thread_num(), c, local);
        }
        global += local;
    }

    //printf("reduced: %.16f, \n", global / data.conf.input_length);
    //double ret = global / data.conf.input_length;

    return global / data.conf.input_length;
}

double 
Parallel::stddev_ind(double* p_array, double p_mean)
{
    double nominator = 0;

    #pragma omp parallel shared(nominator)
    {
        int i = omp_get_thread_num() * (data.conf.input_length / omp_get_num_threads());
        int per = data.conf.input_length / omp_get_num_threads();

        int c;
        double local = 0;

        for (i, c = 0; c < per; c++, i++)
        {
            local += (p_array[i] - p_mean) * (p_array[i] - p_mean);
        }
        nominator += local;
    }

    return sqrt(nominator / data.conf.input_length);
}

double 
Parallel::pearson_ind(
    double* p_array_a, 
    double* p_array_b, 
    double p_mean_a, 
    double p_mean_b, 
    double p_std_dev_a, 
    double p_std_dev_b
)
{
    double sum = 0;

    #pragma omp parallel shared(sum)
    {
        int i = omp_get_thread_num() * (data.conf.input_length / omp_get_num_threads());
        int per = data.conf.input_length / omp_get_num_threads();

        int c;
        double local = 0;

        for (i, c = 0; c < per; c++, i++)
            local += ((p_array_a[i] - p_mean_a) * (p_array_b[i] - p_mean_b));
        
        sum += local;
    }

    return (sum / data.conf.input_length) / (p_std_dev_a * p_std_dev_b);
}



double 
Parallel::mean_for(double* p_array)
{
    double sum = 0;

    #pragma omp parallel for shared(sum)
        for (int i = 0; i < data.conf.input_length; ++i)
        {
            #pragma omp critical
            sum += p_array[i];
        }

    return sum / data.conf.input_length;
}

double 
Parallel::stddev_for(double* p_array, double p_mean)
{
    double nominator = 0;

    #pragma omp parallel for shared(nominator)
        for (int i = 0; i < data.conf.input_length; ++i)
        {
            #pragma omp critical
            //nominator += pow((p_array[i] - p_mean), 2.0);
            nominator += (p_array[i] - p_mean) * (p_array[i] - p_mean);
        }

    return sqrt(nominator / data.conf.input_length);  
}

double 
Parallel::pearson_for(
    double* p_array_a, 
    double* p_array_b, 
    double p_mean_a, 
    double p_mean_b, 
    double p_std_dev_a, 
    double p_std_dev_b
)
{
    double sum = 0;

    #pragma omp parallel for shared(sum)
    for (int i = 0; i < data.conf.input_length; ++i)
    {
        #pragma omp critical
        sum += ((p_array_a[i] - p_mean_a) * (p_array_b[i] - p_mean_b));
    }

    return (sum / data.conf.input_length) / (p_std_dev_a * p_std_dev_b);
}

double 
Parallel::run_parallel_pearson(ParallelType p_type, int p_threads)
{
    printf("\n\tParallel algorithm now running...\n");
    printf("\tNumber of values: %d\n", data.conf.input_length);
    set_threads(p_threads);

    if (p_type == ParallelType::PARALLEL_FOR)
    {
        printf("\tUsing parallel for and critical section...\n");
    }
    else if (p_type == ParallelType::REDUCTION)
    {
        printf("\tUsing parallel for and reduction...\n");
    }
    else if (p_type == ParallelType::INDICES)
    {
        printf("\tUsing pragma parallel with calculated array indices...\n");
    }

    clock_t start, end;
    double diff;
    double mean_a, mean_b;
    double stddev_a, stddev_b;
    double pearson;

    start = clock();

    if (p_type == ParallelType::PARALLEL_FOR)
    {
        mean_a = mean_for(data.a);
        mean_b = mean_for(data.b);
    }
    else if (p_type == ParallelType::REDUCTION)
    {
        mean_a = mean_red(data.a);
        mean_b = mean_red(data.b);
    }
    else if (p_type == ParallelType::INDICES)
    {
        mean_a = mean_ind(data.a);
        mean_b = mean_ind(data.b);
    }
    
    end = clock();
    stats.mean_time = stats.compute_diff(start, end);
    printf("\tParallel mean computed in %*.*f ms\n", RESULT_LENGTH, FLOAT_PRECISION, stats.mean_time);

    #pragma omp barrier
    start = clock();
    if (p_type == ParallelType::PARALLEL_FOR)
    {
        stddev_a = stddev_for(data.a, mean_a);
        stddev_b = stddev_for(data.b, mean_b);
    }
    else if (p_type == ParallelType::REDUCTION)
    {
        stddev_a = stddev_red(data.a, mean_a);
        stddev_b = stddev_red(data.b, mean_b);
    }
    else if (p_type == ParallelType::INDICES)
    {
        stddev_a = stddev_ind(data.a, mean_a);
        stddev_b = stddev_ind(data.b, mean_b);
    }
    end = clock();

    stats.stddev_time = stats.compute_diff(start, end);
    printf("\tParallel stddev computed in %*.*f ms\n", RESULT_LENGTH, FLOAT_PRECISION, stats.stddev_time);

    #pragma omp barrier
    start = clock();
    if (p_type == ParallelType::PARALLEL_FOR)
    {
        pearson = pearson_for(data.a, data.b, mean_a, mean_b, stddev_a, stddev_b);
    }
    else if (p_type == ParallelType::REDUCTION)
    {
        pearson = pearson_red(data.a, data.b, mean_a, mean_b, stddev_a, stddev_b);
    }
    else if (p_type == ParallelType::INDICES)
    {
        pearson = pearson_ind(data.a, data.b, mean_a, mean_b, stddev_a, stddev_b);
    }

    end = clock();
    stats.pearson_time = stats.compute_diff(start, end);
    printf("\tParallel Pearson coefficient computed in %*.*f ms\n", RESULT_LENGTH, FLOAT_PRECISION, stats.pearson_time);

    stats.compute_overall_time();
    printf("\tOverall time: %*.*f ms\n", RESULT_LENGTH, FLOAT_PRECISION, stats.overall_time);

    printf("\n\tResults:\n");

    printf("\tMean X: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, mean_a);
    printf("\tMean Y: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, mean_b);

    printf("\tStd dev X: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, stddev_a);
    printf("\tStd dev Y: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, stddev_b);

    printf("\tPearson coefficient: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, pearson);
}