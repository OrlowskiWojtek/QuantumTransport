#include "utils/include/basis_funcs.hpp"
#include "utils/include/params.hpp"

#include <eigen3/Eigen/Eigen>
#include <cassert>
#include <cmath>
#include <fstream>

QDBasis::QDBasis() {
    p = QDParams::getInstance();
    calculate_indexes();
}

void QDBasis::calculate_indexes() {
    nodes_x.resize(p->N);
    nodes_y.resize(p->N);

    for (int k = 0; k < p->N; k++) {
        int i = static_cast<int>(k) / static_cast<int>(p->n);
        int j = static_cast<int>(k) % static_cast<int>(p->n);

        nodes_x[k] = p->xmin + p->dx * i;
        nodes_y[k] = p->xmin + p->dx * j;
    }
}

double QDBasis::base_gaussian(double x, double y, int k) {
    assert(k < nodes_x.size() && k < nodes_y.size());

    return 1. / std::pow((p->alpha_x * M_PI), 0.25) * std::exp(-std::pow(x - nodes_x[k], 2) / (2 * p->alpha_x)) 
    * 1. / std::pow((p->alpha_y * M_PI), 0.25) *
           std::exp(-std::pow(y - nodes_y[k], 2) / (2. * p->alpha_y));
}

void QDBasis::base_gaussian_to_file(int k) {
    double step = (p->xmax - p->xmin) / 100;

    std::ofstream file("../data/basis_k_" + std::to_string(k) + ".dat");

    file << "# xmax:" << p->xmax << "\n";
    for (double x = p->xmin; x < p->xmax; x += step) {
        for (double y = p->xmin; y < p->xmax; y += step) {
            file << base_gaussian(x, y, k) << "\t";
        }
        file << "\n";
    }
    file.close();
}

double QDBasis::S(int k, int l){
    return std::exp(- std::pow(nodes_x[k] - nodes_x[l], 2) / (4*p->alpha_x) 
                    - std::pow(nodes_y[k] - nodes_y[l], 2) / (4*p->alpha_y));
}

double QDBasis::K(int k, int l){
    return -1. / (2. * p->reducted_mass) * 
            (((std::pow(nodes_x[k] - nodes_x[l], 2) - 2 * p->alpha_x)/ (4 * std::pow(p->alpha_x, 2)))
            + ((std::pow(nodes_y[k] - nodes_y[l], 2) - 2 * p->alpha_y) / (4 * std::pow(p->alpha_y, 2)))) * S(k,l);
}

double QDBasis::V(int k, int l){
    return p->reducted_mass / (2.) * 
    ((std::pow(p->omegax, 2) * (std::pow(nodes_x[k]+nodes_x[l],2) + 2 * p->alpha_x ) / 4)
    +
    (std::pow(p->omegay, 2) * (std::pow(nodes_y[k]+nodes_y[l],2) + 2 * p->alpha_y ) / 4)) * S(k,l);
}

double QDBasis::H(int k, int l){
    return K(k,l) + V(k,l);
}

Eigen::MatrixXd QDBasis::get_H_matrix(){
    Eigen::MatrixXd h_matrix(p->N, p->N);

    for(int k = 0; k < p->N; k++){
        for(int l = 0; l < p->N; l++){
            h_matrix(k,l) = H(k,l);
        }
    }

    return h_matrix;
}


Eigen::MatrixXd QDBasis::get_S_matrix(){
    Eigen::MatrixXd s_matrix(p->N, p->N);

    for(int k = 0; k < p->N; k++){
        for(int l = 0; l < p->N; l++){
            s_matrix(k,l) = S(k,l);
        }
    }

    return s_matrix;
}
