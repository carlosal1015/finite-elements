//
//  main.cpp
//  kurs
//
//  Created by BlackFox on 25.05.15.
//  Copyright (c) 2015 BlackFox. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "finite-element.h"
#include "gz.h"
using namespace::std;
int main(void) {
    clock_t t1 = clock();
    coord** x= (coord**)malloc(sizeof(coord*)*n2);
    for (int j = 0; j < n2; j++) {
        x[j] =(coord*)malloc(sizeof(coord)*n2);
    }
    int count = 1;
    for (int j = 0; j < n2; j++) {
        for (int i = 0; i < n1; i++){
            x[i][j].x = i*h;
            x[i][j].y = j*h;
            x[i][j].num = count;
            count++;
        }
    }
    Allelements(x);
    cout << CheckSymmetric(Kend);
    cout <<'\t'<< m<<endl<< endl;
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++){
            cout << Kend[i][j] <<'\t';
            if(i == m-1) cout <<"\t\t"<<Fend[j]<<endl;
        }
    }
    cout <<endl;
    GaussZ(Kend, Fend, Xend, eps);
    for (int i = 0; i < m; i++){
        cout << Xend[i] <<endl;
    }
    cout << endl<< "Time: " << (double) (clock()-t1) / (double)CLOCKS_PER_SEC << endl;
    return 0;
}
