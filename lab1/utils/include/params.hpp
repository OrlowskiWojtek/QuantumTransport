#ifndef QD_PARAMS_HPP
#define QD_PARAMS_HPP

#include <cmath>
#include <functional>

class QDParams {
private:
    static QDParams *instance;
    QDParams() {}

public:
    QDParams(const QDParams &) = delete;
    QDParams &operator=(const QDParams &) = delete;
    static QDParams *getInstance() {
        if (!instance) {
            instance = new QDParams();
        }

        return instance;
    }

    void set_omega_x(double omx_meV);
    void set_omega_y(double omy_meV);

    double hartree_to_meV = 27211.6;
    double bohr_radius = 0.0529;

    double omegax = 80 / hartree_to_meV;
    double omegay = 200 / hartree_to_meV;
    double dx = 1. / bohr_radius;
    double reducted_mass = 0.24;
    double alpha_x = 1. / (reducted_mass * omegax);
    double alpha_y = 1. / (reducted_mass * omegay);

    std::function<double(double, double)> potential = [&](double x, double y) {
        return 0.5 * reducted_mass * (std::pow(omegax * x, 2) + std::pow(omegay * y, 2));
    };

    int n = 9;
    int N = n * n;        // n x n -> to be set
    double xmin = - (n-1) * dx / 2.;
    double xmax = (n-1) * dx / 2.;
};

#endif
