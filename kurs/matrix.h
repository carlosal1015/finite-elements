#include <cmath>
#include <cstdio>
#include <cstdlib>

inline double mult_vector_to_vector(double* v1, double* v2, int n) {
    double res = 0.0;
    for (int i = 0; i < n; i++) {
        res += v1[i]*v2[i];
    }
    return res;
}

inline double vector_norm_2(double* v, int n) {
    return sqrt(mult_vector_to_vector(v, v, n));
}

inline double* mult_matrix_to_vector(double** A, double* v, int n) {
    double* w = new double[n];
    for (int i = 0; i < n; i++) {
        w[i] = mult_vector_to_vector(A[i], v, n);
    }
    return w;
}

inline double* mult_matrix_to_vector_inplace(double** A, double* v, int n, double* w) {
    for (int i = 0; i < n; i++) {
        w[i] = mult_vector_to_vector(A[i], v, n);
    }
    return w;
}

inline double* mult_ndiag_matrix_to_vector(double* A, double* b, int n)  {
    double* v = new double[n];
    for (int i = 0; i < n1; i++) {
        v[i] = b[i];  // потому что первые n1 строк -- единичная диагональ
    }
    for (int i = n1; i < n-n1; i++) {
        v[i] = A[(i-n1)*4] * (b[i-n1] + b[i+n1]) + A[(i-n1)*4 + 1] * b[i-1]
             + A[(i-n1)*4 + 2] * b[i] + A[(i-n1)*4 + 3] * b[i+1];
    }
    for (int i = n-n1; i < n; i++) {
        v[i] = b[i];  // потому что последние n1 строк -- единичная диагональ
    }
    return v;
}

inline double* mult_ndiag_matrix_to_vector_inplace(double* A, double* b, int n, double* v) {
    for (int i = 0; i < n1; i++) {
        v[i] = b[i];  // потому что первые т1 строк -- единичная диагональ
    }
    for (int i = n1; i < n-n1; i++) {
        v[i] = A[(i-n1)*4] * (b[i-n1] + b[i+n1]) + A[(i-n1)*4 + 1] * b[i-1]
             + A[(i-n1)*4 + 2] * b[i] + A[(i-n1)*4 + 3] * b[i+1];
    }
    for (int i = n-n1; i < n; i++) {
        v[i] = b[i];  // потому что последние т1 строк -- единичная диагональ
    }
    return v;
}

inline double* mult_vector_to_scalar(double* v, double c, int n) {
    double* w = new double[n];
    for (int i = 0; i < n; i++) {
        w[i] = c * v[i];
    }
    return w;
}

inline double* mult_vector_to_scalar_inplace(double* v, double c, int n) {
    for (int i = 0; i < n; i++) {
        v[i] = c * v[i];
    }
    return v;
}

inline double* add_vector_to_vector(double* v1, double* v2, int n) {
    double* w = new double[n];
    for (int i = 0; i < n; i++) {
        w[i] = v1[i] + v2[i];
    }
    return w;
}

inline double* sub_vector_from_vector(double* v1, double* v2, int n) {
    double* w = new double[n];
    for (int i = 0; i < n; i++) {
        w[i] = v1[i] - v2[i];
    }
    return w;
}

inline double* sub_vector_from_vector_inplace(double* v1, double* v2, int n) {
    for (int i = 0; i < n; i++) {
        v1[i] = v1[i] - v2[i];
    }
    return v1;
}

inline double* add_vector_to_vector_inplace(double* v1, double* v2, int n) {
    for (int i = 0; i < n; i++) {
        v1[i] = v1[i] + v2[i];
    }
    return v1;
}

double** read_MatrixMarket(const char* filename, int & n) {
    double** A = NULL;
    FILE* f = fopen(filename, "rt");
    ssize_t read;
    size_t len;
    char* line = NULL;
    n = -1;
    if (f == NULL) {
        fprintf(stderr, "File '%s' not found.", filename);
        exit(1);
    }
    while ((read = getline(&line, &len, f)) != -1) {
        if (line[0] != '%') {
            if (n != -1) {
                int i,j;
                double x;
                sscanf(line, "%d %d %lf", &i, &j, &x);
                A[i-1][j-1] = x;
            } else {
                int nn, mm, ll;
                sscanf(line, "%d %d %d", &nn, &mm, &ll);
                n = (nn >= mm ? nn : mm);
                A = new double*[n];
                for (int i = 0; i < n; i++) {
                    A[i] = new double[n];
                }
            }
        }
    }
    return A;
}

double* LUSolve(double** A, double* b, int n, bool transpose=false) {
    double* tmp = new double[n];
    double* x = new double[n];
    for (int i = 0; i < n; i++) {
        double div = (!transpose ? 1 : A[i][i]);
        tmp[i] = b[i]/div;
        for (int j = 0; j < i; j++) {
            tmp[i] -= A[i][j]*tmp[j]/div;
        }
    }
    for (int i = n-1; i >= 0; i--) {
        double div = (!transpose ? A[i][i] : 1);
        x[i] = tmp[i]/div;
        for (int j = i+1; j < n; j++) {
            x[i] -= A[i][j]*x[j]/div;
        }
    }
    delete[] tmp;
    return x;
}
