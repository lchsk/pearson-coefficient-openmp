#include <iostream>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <string>
#include <string.h>
#include <omp.h>
#include "main.h"
#include "definitions.h"
#include "serial.h"
#include "parallel.h"

using namespace std;

int 
main (int argc, char* argv[])
{
    
    Config conf;
    conf.input_length = 1000000;
    
    Data d(conf);

    Serial s(d);
    Parallel p(d);

    ParallelType type;
    int threads = 0;

    if (argc > 2)
    {
        if (strcmp(argv[1], "for") == 0)
        {
            type = ParallelType::PARALLEL_FOR;
        }
        else if (strcmp(argv[1], "red") == 0)
        {
            type = ParallelType::REDUCTION;
        }
        else if (strcmp(argv[1], "ind") == 0)
        {
            type = ParallelType::INDICES;
        }
        else if (strcmp(argv[1], "pth") == 0)
        {
            type = ParallelType::PTHREADS;
        }
        else if (strcmp(argv[1], "c11") == 0)
        {
            type = ParallelType::CPP11;
        }

        threads = atoi(argv[2]);
    }
    else
    {
        printf("\n\tRUNNING WITH DEFAULT SETTINGS:\n\t2 threads (parallel for/reduction)\n\n");
        threads = 2;
        type = ParallelType::REDUCTION;
    }

    s.run_serial_pearson();

    if (PRINT_REPORT_RESULTS)
    {
        for (int i = 0; i < REPORT_ITERATIONS; ++i)
            p.run_parallel_pearson(type, threads);
    }
    else
        p.run_parallel_pearson(type, threads);

    return 0;
}

Config::Config()
{
    input_length = 100;
}

double
Stats::compute_diff(clock_t time1, clock_t time2)
{
    return ((double)time2 - (double)time1) * 1000.0 / CLOCKS_PER_SEC;
}

double
Stats::compute_diff(double time1, double time2)
{
    return (time2 - time1) * 1000;
}

void 
Stats::compute_overall_time()
{
    overall_time = mean_time + stddev_time + pearson_time;
}

Data::Data(Config& p_c) : conf(p_c)
{
    a = new double[conf.input_length];
    b = new double[conf.input_length];

    for (int i = 0; i < conf.input_length; ++i)
    {
        a[i] = sin(i);
        b[i] = sin(i + 1);
    }
}

Data::~Data()
{
    delete[] a;
    delete[] b;
}