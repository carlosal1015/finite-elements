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
#include <math.h>
using namespace::std;

struct coord
{
    double x;
    double y;
    int num;
    coord()
    {
        x = 0.0;
        y = 0.0;
        num = 0;
    }
    coord(double xx, double yy)
    {
        x = xx;
        y = yy;
        num = 0;
    }
};
const double PI = acos(-1.0);
const double L = 50, M = 50, eps = 6.5;//L - длина прямоугольника, M - ширина
const int n1 = 250, n2 = 250, m = n1*n2, n = 2*(n1-1)*(n2-1);
const double hx = L/(n1-1), hy = M/(n2-1);
double *Kend, *Fend, *Xend;

double lmb(double x, double y){
    return 1;
    return x*x+4*x*y-3;
}

double q(coord x){
    // return -2.0*cos(2.0*x.x);
    return 50.0;
    return 0;
    return 2.0*cos(2.0*x.x);
}
double g(coord x){ // краевые условия на верхней и нижней границе области
    // return cos(x.x)*cos(x.x);
    // return sin(x.x)*sin(x.x) + 2.0*x.x*x.y;
    return (x.y == 0 ? 30 : 220);
}

void element(coord x1, coord x2, coord x3){ //построение матрицы Kij и сливание с матрицей теплопроводности
    double k[3][3];
    coord triangle[3] = {x1, x2, x3};
    double l1 = x3.x-x2.x;
    double l2 = x1.x-x3.x;
    double l3 = x2.x-x1.x;
    double h1 = x2.y-x3.y;
    double h2 = x3.y-x1.y;
    double h3 = x1.y-x2.y;
    
    k[0][0] = (h1*h1+l1*l1)*lmb(x1.x, x1.y)/(2*hx*hy);//матрица кэ на треугольнике
    k[1][0] = (h1*h2+l1*l2)*lmb(x1.x, x1.y)/(2*hx*hy);
    k[2][0] = (h1*h3+l1*l3)*lmb(x1.x, x1.y)/(2*hx*hy);
    k[0][1] = (h1*h2+l1*l2)*lmb(x1.x, x1.y)/(2*hx*hy);
    k[1][1] = (h2*h2+l2*l2)*lmb(x1.x, x1.y)/(2*hx*hy);
    k[2][1] = (h2*h3+l2*l3)*lmb(x1.x, x1.y)/(2*hx*hy);
    k[0][2] = (h1*h3+l1*l3)*lmb(x1.x, x1.y)/(2*hx*hy);
    k[1][2] = (h2*h3+l2*l3)*lmb(x1.x, x1.y)/(2*hx*hy);
    k[2][2] = (h3*h3+l3*l3)*lmb(x1.x, x1.y)/(2*hx*hy);
    
    double f[3];
    f[0] = (2*q(x1)+q(x2)+q(x3))*hx*hy/24;
    f[1] = (q(x1)+2*q(x2)+q(x3))*hx*hy/24;
    f[2] = (q(x1)+q(x2)+2*q(x3))*hx*hy/24;
    
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
            int t1 = triangle[i].num-1, t2 = triangle[j].num-1;
            if (t2 < n1 || t2 >= m-n1 || k[j][i] == 0.0) continue;
            if (t1 == t2-n1 || t1 == t2+n1) Kend[(t2-n1)*4] += k[j][i]*0.5;
            else if (t1 == t2-1) Kend[(t2-n1)*4 + 1] += k[j][i];
            else if (t1 == t2)   Kend[(t2-n1)*4 + 2] += k[j][i];
            else if (t1 == t2+1) Kend[(t2-n1)*4 + 3] += k[j][i];
            else { // если мы всё-таки пропустили какой-нибудь случай
                cerr << "ERROR!" << endl;
                exit(0);
            }
        }
        Fend[triangle[j].num-1]+=f[j];
    }
}

void squareTwoEl(coord x1, coord x2, coord x3, coord x4){ //разбиваем прямоугольную область на два треугольника и вып функцию выше
    element(x1, x2, x3);
    element(x4, x2, x3);
}

void Allelements(coord ** x){ //р
    for (int j = 0; j < m; j++) {
        Fend[j] = 0.0;
        Xend[j] = 0.0;
    }
    for (int j = 0; j < n2-1; j++) {
        for (int i = 0; i < n1-1; i++) {
            squareTwoEl(x[i][j], x[i+1][j], x[i][j+1], x[i+1][j+1]);
        }
    }
    for (int i = 0; i < n1; i++) {
        Fend[i] = g(coord(i*hx, 0));
    }
    for (int i = 0; i < n1; i++) {
        Fend[m-n1+i] = g(coord(i*hx, M));
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
