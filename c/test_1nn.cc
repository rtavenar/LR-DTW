#include <iostream>

#include "ucr.h"
#include "1nn.h"
#include "types.h"

int main(int argc, char **argv) {
    vec_int_t labels_train, labels_test, predicted_labels;
    mat3d_t data_train, data_test;
    float gamma;
    bool entropy_regularized;

    if(argc != 5) {
        std::cout << "usage: " << argv[0] << " fname_train fname_test gamma entropy_regularization[0/1]" << std::endl;
        return -1;
    }
    read_ucr_file(argv[1], &labels_train, &data_train);
    read_ucr_file(argv[2], &labels_test, &data_test);

    gamma = std::stof(argv[3]);
    entropy_regularized = std::stoi(argv[4]);

    predicted_labels = dtw_reg_1nn_clf(data_train, data_test, labels_train, gamma, entropy_regularized);
    float n_cc = 0.;
    for(size_t i=0; i<predicted_labels.size(); ++i)
        if(predicted_labels[i] == labels_test[i])
            n_cc += 1.;
    std::cout << argv[2] << ";" << gamma << ";" << entropy_regularized << ";" << n_cc / predicted_labels.size()
              << std::endl;

    return 0;
}