#include <iostream>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <string>
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

    s.run_serial_pearson();
    p.run_parallel_pearson(ParallelType::INDICES, 2);

    // s_time = omp_get_wtime();
    // double meanX = p.parallel_mean(input_data.a);
    // double meanY = p.parallel_mean(input_data.b);
    // e_time = omp_get_wtime();

    // printf("Mean X: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, meanX);
    // printf("Mean Y: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, meanY);

    // printf("Duration = %lf ms\n", (e_time - s_time) * 1000);

    // double stddevX = p.parallel_stddev(input_data.a, meanX);
    // double stddevY = p.parallel_stddev(input_data.b, meanY);

    // printf("Std dev X: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, stddevX);
    // printf("Std dev Y: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, stddevY);

    // double pearson = p.parallel_pearson(input_data.a, input_data.b, meanX, meanY, stddevX, stddevY);
    
    // printf("Pearson coefficient: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, pearson);

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