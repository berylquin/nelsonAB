#ifndef NELSONC_H
#define NELSONC_H

#define _USE_MATH_DEFINES
#include "math.h"
#include "stdio.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "random"

class nelsonc
{
public:
    nelsonc(int runs);

    void rstart();
    double gauss(double mean, double stddev);
    void schritt();

    int n;
     double start[1];
     double x0;
     double pos[1];
    double r=0.0;
    double a=2.0;
    double s=20.0;
    double dt=0.05;
    double ti=0.0;
    double B=0.01;

};

#endif // NELSONC_H
