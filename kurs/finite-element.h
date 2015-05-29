//
//  mf.h
//  kurs
//
//  Created by BlackFox on 25.05.15.
//  Copyright (c) 2015 BlackFox. All rights reserved.
//

#ifndef kurs_finiteElement_h
#define kurs_finiteElement_h

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
using namespace::std;

struct coord
{
    double x = 0.0;
    double y = 0.0;
    int num = 0;
};
double** Kend, *Fend, *Xend;
double L = 2, M = 4, eps = 0.001;//L - длина прямоугольника, M - ширина
int n1 = 6, n2 = 6, m = n1*n2, n = 2*(n1-1)*(n2-1);
double h = L/(n1-1);

double lmb(double x, double y){
    return 1;//x*x+4*x*y-3;
}

double q(coord x){
    return 1;//2*x.x*x.y+2;
}

void MultiplyMrtxOmK(double ** a, double** b, double** c) { //умножение матрицы омега на Kij
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < 3; col++) {
            for (int inner = 0; inner < 3; inner++) {
                c[row][col] += a[row][inner] * b[inner][col];
            }
        }
    }
}
void MultiplyMrtx(double ** a, double** b, double** c) { //произведение Omega*Kij*Omega^t
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < m; col++) {
            for (int inner = 0; inner < m; inner++) {
                c[row][col] += a[row][inner] * b[col][inner];
            }
        }
    }
}
void MultiplyMrtxVector(double ** a, double* b, double* c) { //произведение матрицы на вектор
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < 3; i++) {
            c[j] += a[j][i] * b[i];
        }
    }
}

void element(coord x1, coord x2, coord x3){ //построение матрицы Kij и сливание с матрицей теплопроводности
    coord triangle[3] = {x1, x2, x3};
    double l1 = x3.x-x2.x;
    double l2 = x1.x-x3.x;
    double l3 = x2.x-x1.x;
    double h1 = x2.y-x3.y;
    double h2 = x3.y-x1.y;
    double h3 = x1.y-x2.y;
    
    double** k = (double**)malloc(sizeof(double*)*3);
    double** tmp1 = (double**)malloc(sizeof(double*)*m);
    double** tmp2 = (double**)malloc(sizeof(double*)*m);
    double** omeg = (double**)malloc(sizeof(double*)*m);
    double* tmpf = (double*)malloc(sizeof(double)*m);
    for (int j = 0; j < m; j++) {
        omeg[j] = (double*)malloc(sizeof(double)*3);
        tmp1[j] = (double*)malloc(sizeof(double)*m);
        tmp2[j] = (double*)malloc(sizeof(double)*m);
    };
    for (int j = 0; j < 3; j++) {
        k[j] =(double*)malloc(sizeof(double)*3);
    }
    k[0][0] = h1*h1+l1*l1*lmb(x1.x, x1.y)/(2*h*h);//матрица кэ на треугольнике
    k[1][0] = h1*h2+l1*l2*lmb(x1.x, x1.y)/(2*h*h);
    k[2][0] = h1*h3+l1*l3*lmb(x1.x, x1.y)/(2*h*h);
    k[0][1] = h1*h2+l1*l2*lmb(x1.x, x1.y)/(2*h*h);
    k[1][1] = h2*h2+l2*l2*lmb(x1.x, x1.y)/(2*h*h);
    k[2][1] = h2*h3+l2*l3*lmb(x1.x, x1.y)/(2*h*h);
    k[0][2] = h1*h3+l1*l3*lmb(x1.x, x1.y)/(2*h*h);
    k[1][2] = h2*h3+l2*l3*lmb(x1.x, x1.y)/(2*h*h);
    k[2][2] = h3*h3+l3*l3*lmb(x1.x, x1.y)/(2*h*h);
    
    double* f = (double*)malloc(sizeof(double)*3);
    f[0] = (2*q(x1)+q(x2)+q(x3))*h*h/24;
    f[1] = (q(x1)+2*q(x2)+q(x3))*h*h/24;
    f[2] = (q(x1)+q(x2)+2*q(x3))*h*h/24;
    
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            tmp1[i][j] = 0.0;
            tmp2[i][j] = 0.0;
            tmpf[j] = 0.0;
        }
    }
    for (int j = 0; j < 3; j++) {
        for (int i = 0;  i < m; i++) {
            if (i+1 == triangle[j].num) omeg[i][j] = 1.0;
            else omeg[i][j] = 0.0;
        }
    }
    MultiplyMrtxOmK(omeg, k, tmp1);
    MultiplyMrtx(tmp1, omeg, tmp2);
    MultiplyMrtxVector(omeg, f, tmpf);
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++){
            cout << k[i][j] <<'\t';
            if(i == 2) cout <<endl;
        }
    }
    cout <<endl;
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            Kend[i][j]+=tmp2[i][j];
        }
        Fend[j]+=tmpf[j];
    }
}

void squareTwoEl(coord x1, coord x2, coord x3, coord x4){ //разбиваем прямоугольную область на два треугольника и вып функцию выше
    element(x1, x2, x3);
    element(x4, x2, x3);
}

void Allelements(coord ** x){ //р
    Kend = (double**)malloc(sizeof(double*)*m);
    Fend = (double*)malloc(sizeof(double)*m);
    Xend = (double*)malloc(sizeof(double)*m);
    for (int j = 0; j < m; j++) {
        Kend[j] =(double*)malloc(sizeof(double)*m);
    }
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            Kend[i][j] = 0.0;
        }
        Fend[j] = 0.0;
        Xend[j] = 0.0;
    }
    for (int j = 0; j < n2-1; j++) {
        for (int i = 0; i < n1-1; i++) {
            squareTwoEl(x[i][j], x[i+1][j], x[i][j+1], x[i+1][j+1]);
        }
    }
}

int CheckSymmetric(double** fun){
    for (int j = 0; j < m; j++) {
        for (int i = 0; i< m; i++) {
            if (fun[i][j] != fun[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

#endif
