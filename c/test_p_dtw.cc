#include <iostream>

#include "p_dtw.h"

#define N 50

int main(int argc, char **argv)
{
    mat_t s_x, s_y;
    mat3d_t probas;

    s_x.resize(N);
    for(size_t i=0; i<N; ++i)
        s_x[i].resize(1);

    s_y.resize(N);
    for(size_t i=0; i<N; ++i)
        s_y[i].resize(1);

    std::cout << p_dtw(s_x, s_y, 10., &probas) << std::endl;
    return 0;
}