//
//  gz.h
//  kurs
//
//  Created by BlackFox on 28.05.15.
//  Copyright (c) 2015 BlackFox. All rights reserved.
//

#ifndef kurs_gz_h
#define kurs_gz_h
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "finite-element.h"

bool converge(double *xk, double *xkp)
{
    double norm = 0;
    for (int i = 0; i < m; i++)
    {
        norm += (xk[i] - xkp[i])*(xk[i] - xkp[i]);
    }
    if(sqrt(norm) >= eps)
        return false;
    return true;
}

void GaussZ(double** a, double* b, double* x, double e)
{
    double*p = (double*)malloc(sizeof(double)*m);
    do
    {
        for (int i = 0; i < m; i++)
            p[i] = x[i];
        
        for (int i = 0; i < m; i++)
        {
            double var = 0;
            for (int j = 0; j < i; j++)
                var += (a[j][i] * x[j]);
            for (int j = i + 1; j < m; j++)
                var += (a[j][i] * p[j]);
            x[i] = (b[i] - var) / a[i][i];
        }
    }
    while (!converge(x, p));
}
#endif
