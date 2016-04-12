#include <vector>
#include <cmath>
#include <iostream>
#include "p_dtw.h"

float d(const vec_t x, const vec_t y) {
    float dst = 0.;
    for(size_t i=0; i<x.size(); ++i) {
        float diff = (x[i] - y[i]);
        dst += diff * diff;
    }
    return sqrt(dst);
}

void cdist(const mat_t s_x, const mat_t s_y, mat_t *distances) {
    (*distances).resize(s_x.size());
    for(size_t i=0; i<s_x.size(); ++i) {
        (*distances)[i].resize(s_y.size());
        for(size_t j=0; j<s_y.size(); ++j) {
            (*distances)[i][j] = d(s_x[i], s_y[j]);
        }
    }
}

vec_t get_probas_formula(const float cost_up, const float cost_right, const float cost_diagonal, const float gamma) {
    vec_t probas;
    probas.resize(3);
    std::fill(probas.begin(), probas.end(), 0.);

    if(gamma < 1e-12) {
        float min_val = std::min(std::min(cost_up, cost_right), cost_diagonal);
        if(cost_up == min_val)
            probas[UP] = 1.;
        if(cost_right == min_val)
            probas[RIGHT] = 1.;
        if(cost_diagonal == min_val)
            probas[DIAGONAL] = 1.;
        float s = 0;
        for(size_t i=0; i<3; ++i)
            s += probas[i];
        for(size_t i=0; i<3; ++i)
            probas[i] /= s;
    } else {
        float p_up = 1. / 3. * (1. + (cost_right + cost_diagonal - 2. * cost_up) / (2. * gamma));
        float p_right = 1. / 3. * (1. + (cost_up + cost_diagonal - 2. * cost_right) / (2. * gamma));
        if(p_up >= 0. && p_right >= 0. && p_up + p_right <= 1.) {
            probas[UP] = p_up;
            probas[RIGHT] = p_right;
            probas[DIAGONAL] = 1. - p_up - p_right;
        } else if(p_up < 0.) {
            p_right = .5 * (1. + (cost_diagonal - cost_right) / (2 * gamma));
            if(p_right <= 1. && p_right >= 0.) {
                probas[UP] = 0.;
                probas[RIGHT] = p_right;
                probas[DIAGONAL] = 1. - p_right;
            } else if(p_right < 0.) {
                probas[UP] = 0.;
                probas[RIGHT] = 0.;
                probas[DIAGONAL] = 1.;
            } else {
                probas[UP] = 0.;
                probas[RIGHT] = 1.;
                probas[DIAGONAL] = 0.;
            }
        } else if(p_right < 0.) {
            p_up = .5 * (1. + (cost_diagonal - cost_up) / (2 * gamma));
            if(p_up <= 1. && p_up >= 0.) {
                probas[UP] = p_up;
                probas[RIGHT] = 0.;
                probas[DIAGONAL] = 1. - p_up;
            } else if(p_up < 0.) {
                probas[UP] = 0.;
                probas[RIGHT] = 0.;
                probas[DIAGONAL] = 1.;
            } else {
                probas[UP] = 1.;
                probas[RIGHT] = 0.;
                probas[DIAGONAL] = 0.;
            }
        } else {
            p_up = .5 * (1. + (cost_right - cost_up) / (2 * gamma));
            if(p_up <= 1. && p_up >= 0.) {
                probas[UP] = p_up;
                probas[RIGHT] = 1. - p_up;
                probas[DIAGONAL] = 0.;
            } else if(p_up < 0.) {
                probas[UP] = 0.;
                probas[RIGHT] = 1.;
                probas[DIAGONAL] = 0.;
            } else {
                probas[UP] = 1.;
                probas[RIGHT] = 0.;
                probas[DIAGONAL] = 0.;
            }
        }
    }
    return probas;
}

float p_dtw(const mat_t s_x, const mat_t s_y, const float gamma, mat3d_t *probas) {
    mat_t distances, mat_cost;

    cdist(s_x, s_y, &distances);
    (*probas).resize(s_x.size());
    mat_cost.resize(s_x.size());
    for(size_t i=0; i<s_x.size(); ++i) {
        (*probas)[i].resize(s_y.size());
        mat_cost[i].resize(s_y.size());
        std::fill(mat_cost[i].begin(), mat_cost[i].end(), 0.);
    }

    for(size_t i=1; i<s_x.size(); ++i) {
        (*probas)[i][0].resize(3);
        std::fill((*probas)[i][0].begin(), (*probas)[i][0].end(), 0.);
        (*probas)[i][0][UP] = 1.;
        mat_cost[i][0] = mat_cost[i - 1][0] + distances[i][0];
    }
    for(size_t j=1; j<s_y.size(); ++j) {
        (*probas)[0][j].resize(3);
        std::fill((*probas)[0][j].begin(), (*probas)[0][j].end(), 0.);
        (*probas)[0][j][RIGHT] = 1.;
        mat_cost[0][j] = mat_cost[0][j - 1] + distances[0][j];
    }

    for(size_t i=1; i<s_x.size(); ++i) {
        for(size_t j=1; j<s_y.size(); ++j) {
            (*probas)[i][j] = get_probas_formula(mat_cost[i - 1][j], mat_cost[i][j - 1], mat_cost[i - 1][j - 1], gamma);
            mat_cost[i][j] = (*probas)[i][j][UP] * mat_cost[i - 1][j] + (*probas)[i][j][RIGHT] * mat_cost[i][j - 1] +
                (*probas)[i][j][DIAGONAL] * mat_cost[i - 1][j - 1] + distances[i][j];
        }
    }

    return mat_cost[s_x.size() - 1][s_y.size() - 1];
}

void p_dtw_backtrace(const mat3d_t probas, mat_t *mat_probas) {
    (*mat_probas).resize(probas.size());
    for(size_t i=0; i<(*mat_probas).size(); ++i)
        (*mat_probas)[i].resize(probas[i].size());
    (*mat_probas)[probas.size() - 1][probas[0].size() - 1] = 1.;
    for(int i=probas.size() - 2; i>=0; --i) {
        size_t m = (*mat_probas)[i].size();
        (*mat_probas)[i][m - 1] = (*mat_probas)[i + 1][m - 1] * probas[i + 1][m - 1][UP];
    }
    for(int j=probas[0].size() - 2; j>=0; --j) {
        size_t n = (*mat_probas).size();
        (*mat_probas)[n - 1][j] = (*mat_probas)[n - 1][j + 1] * probas[n - 1][j + 1][RIGHT];
    }
    for(int i=probas.size() - 2; i>=0; --i) {
        for(int j=probas[i].size() - 2; j>=0; --j) {
            (*mat_probas)[i][j] = (*mat_probas)[i + 1][j] * probas[i + 1][j][UP] +
                                  (*mat_probas)[i][j + 1] * probas[i][j + 1][RIGHT] +
                                  (*mat_probas)[i + 1][j + 1] * probas[i + 1][j + 1][DIAGONAL];
        }
    }
}
