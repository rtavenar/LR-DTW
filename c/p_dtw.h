#include <vector>
#include <cmath>

#define UP 0
#define RIGHT 1
#define DIAGONAL 2

typedef std::vector< std::vector<float> > mat_t;
typedef std::vector< std::vector< std::vector<float> > > mat3d_t;
typedef std::vector<float> vec_t;

float d(const vec_t x, const vec_t y);
void cdist(const mat_t s_x, const mat_t s_y, mat_t *distances);
vec_t get_probas_formula(const float cost_up, const float cost_right, const float cost_diagonal, const float gamma,
                         const bool entropy_regularized);
float p_dtw(const mat_t s_x, const mat_t s_y, const float gamma, mat3d_t *probas, const bool entropy_regularized);
void p_dtw_backtrace(const mat3d_t probas, mat_t *mat_probas);
