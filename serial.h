#ifndef SERIAL_H
#define SERIAL_H

#include "main.h"
#include "definitions.h"

class Serial
{
public:
    Serial(Data& p_data);

    double compute_mean(double* p_array);
    double compute_std_dev(double* p_array, double p_mean);
    double compute_pearson(
        double* p_array_a, 
        double* p_array_b, 
        double p_mean_a, 
        double p_mean_b, 
        double p_std_dev_a, 
        double p_std_dev_b
    );
    double run_serial_pearson();

private:
    Data& data;
    Stats stats;
};

#endif