#include "serial.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "definitions.h"
#include "parallel.h"
#include "main.h"

Serial::Serial(Data& p_data) : data(p_data)
{
}

double 
Serial::compute_mean(double* p_array)
{
    double sum = 0;

    int i = 0;
    for (i; i < data.conf.input_length; ++i)
    {
        sum = sum + p_array[i];
    }

    return sum / data.conf.input_length;
}

double 
Serial::compute_std_dev(double* p_array, double p_mean)
{
    double nominator = 0;

    int i = 0;
    for (i; i < data.conf.input_length; ++i)
    {
        nominator += (p_array[i] - p_mean) * (p_array[i] - p_mean);
    }

    return sqrt(nominator / data.conf.input_length);
}

double 
Serial::compute_pearson(
    double* p_array_a, 
    double* p_array_b, 
    double p_mean_a, 
    double p_mean_b, 
    double p_std_dev_a, 
    double p_std_dev_b
    )
{
    int i = 0;
    double sum = 0;
    for (i; i < data.conf.input_length; ++i)
    {
        sum += ((p_array_a[i] - p_mean_a) * (p_array_b[i] - p_mean_b));
    }

    return (sum / data.conf.input_length) / (p_std_dev_a * p_std_dev_b);
}

double 
Serial::run_serial_pearson()
{
    clock_t start, end;
    double diff;

    printf("\tSerial algorithm now running...\n");
    printf("\tNumber of values: %d\n", data.conf.input_length);

    start = clock();

    double mean_a = compute_mean(data.a);
    double mean_b = compute_mean(data.b);

    end = clock();
    stats.mean_time = stats.compute_diff(start, end);
    printf("\tSerial mean computed in %*.*f ms\n", RESULT_LENGTH, FLOAT_PRECISION, stats.mean_time);

    start = clock();
    double stddev_a = compute_std_dev(data.a, mean_a);
    double stddev_b = compute_std_dev(data.b, mean_b);
    end = clock();

    stats.stddev_time = stats.compute_diff(start, end);
    printf("\tSerial stddev computed in %*.*f ms\n", RESULT_LENGTH, FLOAT_PRECISION, stats.stddev_time);

    start = clock();
    double pearson = compute_pearson(data.a, data.b, mean_a, mean_b, stddev_a, stddev_b);
    end = clock();
    stats.pearson_time = stats.compute_diff(start, end);
    printf("\tSerial Pearson coefficient computed in %*.*f ms\n", RESULT_LENGTH, FLOAT_PRECISION, stats.pearson_time);

    stats.compute_overall_time();
    printf("\tOverall time: %*.*f ms\n", RESULT_LENGTH, FLOAT_PRECISION, stats.overall_time);

    printf("\n\tResults:\n");

    printf("\tMean X: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, mean_a);
    printf("\tMean Y: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, mean_b);

    printf("\tStd dev X: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, stddev_a);
    printf("\tStd dev Y: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, stddev_b);

    printf("\tPearson coefficient: %*.*f\n", RESULT_LENGTH, FLOAT_PRECISION, pearson);

    return 0;
}