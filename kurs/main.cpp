//
//  main.cpp
//  kurs
//
//  Created by BlackFox on 25.05.15.
//  Copyright (c) 2015 BlackFox. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "finite-element.h"
#include "gz.h"
#include "fom.h"

using namespace::std;
int main(void) {
    clock_t t1 = clock();
    coord** x= (coord**)malloc(sizeof(coord*)*n1);
    for (int j = 0; j < n1; j++) {
        x[j] =(coord*)malloc(sizeof(coord)*n2);
    }
    Kend = (double*)malloc(sizeof(double*)*(m-n1*2)*4);
    Fend = (double*)malloc(sizeof(double)*m);
    Xend = (double*)malloc(sizeof(double)*m);
    for (int i = 0; i < 4*(m-n1*2); i++) {
        Kend[i] = 0.0;
    }
    int count = 1;
    for (int j = 0; j < n2; j++) {
        for (int i = 0; i < n1; i++){
            x[i][j].x = i*hx;
            x[i][j].y = j*hy;
            x[i][j].num = count;
            count++;
        }
    }
    Allelements(x);
    fom(Kend, Fend, Xend, m);
    for (int i = 0; i < m; i++){
        cout << fixed << setprecision(5) << (i%n1)*hx << " " << (i/n1)*hy << " " << Xend[i] << endl;
    }
    cerr << endl<< "Time: " << (double) (clock()-t1) / (double)CLOCKS_PER_SEC << endl;
    for (int j = 0; j < n1; j++) {
        free(x[j]);
    }
    free(x);
    free(Kend);
    free(Fend);
    free(Xend);
    return 0;
}
