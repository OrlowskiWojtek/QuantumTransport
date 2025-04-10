// basis but also noding system
// indexing from 

#ifndef QD_BASIS_HPP
#define QD_BASIS_HPP

#include "params.hpp"
#include <vector>
#include <eigen3/Eigen/Eigen>

class QDBasis{
private:
    // nodes;
    std::vector<double> nodes_x;
    std::vector<double> nodes_y;

    QDParams* p;

    double H(int k, int l);
    double S(int k, int l);
    double K(int k, int l);
    double V(int k, int l);
public:
    QDBasis();

    void base_gaussian_to_file(int k); // in some easy to read in julia format     
    double base_gaussian(double x, double y, int k); // k -> node 

    Eigen::MatrixXd get_H_matrix();
    Eigen::MatrixXd get_S_matrix();

    void calculate_indexes();
    
};

#endif
