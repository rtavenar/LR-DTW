#include "1nn.h"
#include "p_dtw.h"

vec_int_t dtw_reg_1nn_clf(const mat3d_t train_data, const mat3d_t test_data, const vec_int_t train_labels,
                          const float gamma, const bool entropy_regularized) {
    vec_int_t predicted_labels;
    for(size_t i=0; i<test_data.size(); ++i) {
        //std::cout << i << std::endl;
        vec_t distances;
        int idx;
        for(size_t j=0; j<train_data.size(); ++j) {
            mat3d_t useless_probas;
            float cur_dist = p_dtw(test_data[i], train_data[j], gamma, &useless_probas, entropy_regularized);
            distances.push_back(cur_dist);
            useless_probas.clear();
        }
        idx = std::min_element(distances.begin(), distances.end()) - distances.begin();
        predicted_labels.push_back(train_labels[idx]);
        distances.clear();
    }

    return predicted_labels;
}