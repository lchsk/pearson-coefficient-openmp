#include "definitions.h"
#include "parallel.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "main.h"
#include <pthread.h>
#include <thread>

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

/* -------------------
REDUCTION 
--------------------*/
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
// ----------------------------

/* -------------------
INDICES
--------------------*/
double 
Parallel::mean_ind(double* p_array)
{
    double global = 0;

    #pragma omp parallel shared(global)
    {
        int i = omp_get_thread_num() * (data.conf.input_length / omp_get_num_threads());
        int per = data.conf.input_length / omp_get_num_threads();
        
        int c;
        double local = 0;

        for (i, c = 0; c < per; c++, i++)
        {
            local += p_array[i];
        }
        global += local;
    }

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
// ------------------------------

/* -------------------
PHTREADS
--------------------*/
void* 
mean_ph_worker(void* d)
{
    t_pdata* p = (t_pdata*) d;
    double local_sum = 0;

    for (int i = 0, index = p->index; i < p->length; i++, index++)
    {
        local_sum += p->array[index];
    }

    p->global = local_sum;
}

void* 
stddev_ph_worker(void* d)
{
    t_pdata* p = (t_pdata*) d;
    double local_sum = 0;

    for (int i = 0, index = p->index; i < p->length; i++, index++)
    {
        local_sum += ((p->array[index] - p->mean_a) * (p->array[index] - p->mean_a));
    }

    p->global = local_sum;
}

void* 
pearson_ph_worker(void* d)
{
    t_pdata* p = (t_pdata*) d;
    double local_sum = 0;

    for (int i = 0, index = p->index; i < p->length; i++, index++)
    {
        local_sum += ((p->array[index] - p->mean_a) * (p->array_b[index] - p->mean_b));
    }

    p->global = local_sum;
}

double 
Parallel::dispatcher_ph(
    double* p_array, 
    int p_what, 
    double p_mean_a,
    double p_mean_b,
    double p_stddev_a,
    double p_stddev_b,
    double* p_array_b
)
{
    int ret;
    pthread_t* thr = new pthread_t[threads - 1];
    t_pdata* pdat = new t_pdata[threads];
    int per_thread = data.conf.input_length / threads;

    for (int i = 0; i < threads - 1; ++i)
    {
        pdat[i].array = p_array;
        pdat[i].array_b = p_array_b;
        pdat[i].index = i * per_thread;
        pdat[i].length = per_thread;
        pdat[i].mean_a = p_mean_a;
        pdat[i].mean_b = p_mean_b;

        switch(p_what)
        {
            case 0:
                ret = pthread_create(&thr[i], NULL, mean_ph_worker, (void*) &pdat[i]);
                break;
            case 1: 
                ret = pthread_create(&thr[i], NULL, stddev_ph_worker, (void*) &pdat[i]);
                break;
            case 2: 
                ret = pthread_create(&thr[i], NULL, pearson_ph_worker, (void*) &pdat[i]);
                break;
        }
        
        if (ret)
        {
            fprintf(stderr,"Error - pthread_create() return code: %d\n", ret);
            exit(EXIT_FAILURE);
        }
    }

    // make master thread busy as well
    pdat[threads-1].index = pdat[threads-2].index + per_thread;
    pdat[threads-1].length = per_thread;
    pdat[threads-1].array = p_array;
    pdat[threads-1].array_b = p_array_b;
    pdat[threads-1].mean_a = p_mean_a;
    pdat[threads-1].mean_b = p_mean_b;

    switch(p_what)
    {
        case 0:
            mean_ph_worker((void*) &pdat[threads-1]);
            break;
        case 1: 
            stddev_ph_worker((void*) &pdat[threads-1]);
            break;
        case 2: 
            pearson_ph_worker((void*) &pdat[threads-1]);
            break;
    }

    double global = 0;

    for (int i = 0; i < threads - 1; ++i)
    {
        pthread_join(thr[i], NULL);
        global += pdat[i].global;
    }

    global += pdat[threads-1].global;

    delete[] pdat;
    delete[] thr;

    switch(p_what)
    {
        case 0:
            return global / data.conf.input_length;
        case 1: 
            return sqrt(global / data.conf.input_length);
        case 2:
            return (global / data.conf.input_length) / (p_stddev_a * p_stddev_b);
    }
}
// ------------------------------


/* -------------------
C++11 Multithreading
--------------------*/
double 
Parallel::mean_cpp11(double* p_array, int p_index, int p_length, double* global)
{
    double local_sum = 0;

    for (int i = 0, index = p_index; i < p_length; i++, index++)
    {
        local_sum += p_array[index];
    }

    *global = *global + local_sum;
}

double 
Parallel::stddev_cpp11(double* p_array, int p_index, int p_length, double p_mean, double* global)
{
    double local_sum = 0;

    for (int i = 0, index = p_index; i < p_length; i++, index++)
    {
        local_sum += ((p_array[index] - p_mean) * (p_array[index] - p_mean));
    }

    *global = *global + local_sum;
}

double 
Parallel::pearson_cpp11(double* p_array_a, double* p_array_b, int p_index, int p_length, double p_mean_a, double p_mean_b, double* global)
{
    double local_sum = 0;

    for (int i = 0, index = p_index; i < p_length; i++, index++)
    {
        local_sum += ((p_array_a[index] - p_mean_a) * (p_array_b[index] - p_mean_b));
    }

    *global = *global + local_sum;
}


double 
Parallel::dispatcher_cpp11(
    double* p_array, 
    int p_what, 
    double p_mean_a,
    double p_mean_b,
    double p_stddev_a,
    double p_stddev_b,
    double* p_array_b
)
{
    std::thread* thr = new std::thread[threads - 1];
    int per_thread = data.conf.input_length / threads;
    double global_mean = 0;
    double global_stddev = 0;
    double global_pearson = 0;

    for (int i = 0; i < threads - 1; ++i)
    {
        switch(p_what)
        {
            case 0:
                thr[i] = std::thread(&Parallel::mean_cpp11, this, p_array, i * per_thread, per_thread, &global_mean);
                break;
            case 1: 
                thr[i] = std::thread(&Parallel::stddev_cpp11, this, p_array, i * per_thread, per_thread, p_mean_a, &global_stddev);
                break;
            case 2: 
                thr[i] = std::thread(&Parallel::pearson_cpp11, this, p_array, p_array_b, i * per_thread, per_thread, p_mean_a, p_mean_b, &global_pearson);
                break;
        }
    }

    // make master thread busy as well
    switch(p_what)
    {
        case 0:
            mean_cpp11(p_array, data.conf.input_length - per_thread, per_thread, &global_mean);
            break;
        case 1:
            stddev_cpp11(p_array, data.conf.input_length - per_thread, per_thread, p_mean_a, &global_stddev); 
            break;
        case 2: 
            pearson_cpp11(p_array, p_array_b, data.conf.input_length - per_thread, per_thread, p_mean_a, p_mean_b, &global_pearson);
            break;
    }

    for (int i = 0; i < threads - 1; ++i)
    {
        thr[i].join();
    }

    delete[] thr;

    switch(p_what)
    {
        case 0:
            return global_mean / data.conf.input_length;
        case 1: 
            return sqrt(global_stddev / data.conf.input_length);
        case 2:
            return (global_pearson / data.conf.input_length) / (p_stddev_a * p_stddev_b);
    }
}
//----------------------------------

/* -------------------
PARALLEL FOR/CRITICAL SECTION
--------------------*/
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
    if ( ! PRINT_REPORT_RESULTS)
    {
        printf("\n\tParallel algorithm now running...\n");
        printf("\tNumber of values: %d\n", data.conf.input_length);
    }
    set_threads(p_threads);

    if ( ! PRINT_REPORT_RESULTS)
    {
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
        else if (p_type == ParallelType::PTHREADS)
        {
            printf("\tUsing p_threads...\n");
        }
        else if (p_type == ParallelType::CPP11)
        {
            printf("\tUsing C++11 multithreading...\n");
        }
    }

    clock_t start, end;
    double s, e;
    double diff;
    double mean_a, mean_b;
    double stddev_a, stddev_b;
    double pearson;

    s = omp_get_wtime();
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
    else if (p_type == ParallelType::PTHREADS)
    {
        mean_a = dispatcher_ph(data.a, 0, 0, 0, 0, 0, NULL);
        mean_b = dispatcher_ph(data.b, 0, 0, 0, 0, 0, NULL);
    }
    else if (p_type == ParallelType::CPP11)
    {
        mean_a = dispatcher_cpp11(data.a, 0, 0, 0, 0, 0, NULL);
        mean_b = dispatcher_cpp11(data.b, 0, 0, 0, 0, 0, NULL);
    }
    
    e = omp_get_wtime();
    stats.mean_time = stats.compute_diff(s, e);
    if ( ! PRINT_REPORT_RESULTS)
        printf("\tParallel mean computed in %*.*f ms\n", RESULT_LENGTH, FLOAT_PRECISION, stats.mean_time);

    #pragma omp barrier
    s = omp_get_wtime();
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
    else if (p_type == ParallelType::PTHREADS)
    {
        stddev_a = dispatcher_ph(data.a, 1, mean_a, 0, 0, 0, NULL);
        stddev_b = dispatcher_ph(data.b, 1, mean_b, 0, 0, 0, NULL);
    }
    else if (p_type == ParallelType::CPP11)
    {
        stddev_a = dispatcher_cpp11(data.a, 1, mean_a, 0, 0, 0, NULL);
        stddev_b = dispatcher_cpp11(data.b, 1, mean_b, 0, 0, 0, NULL);
    }
    e = omp_get_wtime();

    stats.stddev_time = stats.compute_diff(s, e);
    if ( ! PRINT_REPORT_RESULTS)
        printf("\tParallel stddev computed in %*.*f ms\n", RESULT_LENGTH, FLOAT_PRECISION, stats.stddev_time);

    #pragma omp barrier
    s = omp_get_wtime();
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
    else if (p_type == ParallelType::PTHREADS)
    {
        pearson = dispatcher_ph(data.a, 2, mean_a, mean_b, stddev_a, stddev_b, data.b);
    }
    else if (p_type == ParallelType::CPP11)
    {
        pearson = dispatcher_ph(data.a, 2, mean_a, mean_b, stddev_a, stddev_b, data.b);
    }
    e = omp_get_wtime();
    stats.pearson_time = stats.compute_diff(s, e);
    if ( ! PRINT_REPORT_RESULTS)
        printf("\tParallel Pearson coefficient computed in %*.*f ms\n", RESULT_LENGTH, FLOAT_PRECISION, stats.pearson_time);

    stats.compute_overall_time();
    if ( ! PRINT_REPORT_RESULTS)
        printf("\tOverall time: %*.*f ms\n", RESULT_LENGTH, FLOAT_PRECISION, stats.overall_time);

    if (PRINT_REPORT_RESULTS)
    {
        printf("%*.*f\t%*.*f\t%*.*f\t%*.*f\n", 
            RESULT_LENGTH, FLOAT_PRECISION, stats.mean_time,
            RESULT_LENGTH, FLOAT_PRECISION, stats.stddev_time,
            RESULT_LENGTH, FLOAT_PRECISION, stats.pearson_time,
            RESULT_LENGTH, FLOAT_PRECISION, stats.overall_time
        );
    }

    if ( ! PRINT_REPORT_RESULTS)
    {
        printf("\n\tResults:\n");
        printf("\tMean X: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, mean_a);
        printf("\tMean Y: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, mean_b);

        printf("\tStd dev X: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, stddev_a);
        printf("\tStd dev Y: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, stddev_b);

        printf("\tPearson coefficient: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, pearson);        
    }
}