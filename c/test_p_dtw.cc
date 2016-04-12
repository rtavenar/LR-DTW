#include <iostream>

#include "p_dtw.h"

#define N 50

int main(int argc, char **argv)
{
    mat_t s_x, s_y, mat_probas;
    mat3d_t probas;

    s_x.resize(N);
    for(size_t i=0; i<N; ++i)
        s_x[i].resize(1);

    s_y.resize(N);
    for(size_t i=0; i<N; ++i)
        s_y[i].resize(1);

    p_dtw(s_x, s_y, 10., &probas);
    p_dtw_backtrace(probas, &mat_probas);

    for(size_t i=0; i<N; ++i) {
        for(size_t j=0; j<N; ++j) {
            std::cout << mat_probas[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}