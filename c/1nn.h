#include <vector>
#include <algorithm>
#include<iostream>

#include "types.h"

vec_int_t dtw_reg_1nn_clf(const mat3d_t train_data, const mat3d_t test_data, const vec_int_t train_labels,
                          const float gamma, const bool entropy_regularized);

