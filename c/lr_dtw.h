#include <vector>
#include <cmath>

#include "types.h"

#define UP 0
#define RIGHT 1
#define DIAGONAL 2

float d(const vec_t x, const vec_t y);
void cdist(const mat_t s_x, const mat_t s_y, mat_t *distances);
vec_t get_probas_formula(float cost_up, float cost_right, float cost_diagonal, const float gamma,
                         const bool entropy_regularized, const bool normalized_costs);
float lr_dtw(const mat_t s_x, const mat_t s_y, const float gamma, mat3d_t *probas, const bool entropy_regularized,
             const bool normalized_costs);
void lr_dtw_backtrace(const mat3d_t probas, mat_t *mat_probas);
