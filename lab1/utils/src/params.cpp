#include "utils/include/params.hpp"

QDParams *QDParams::instance = nullptr;

void QDParams::set_omega_x(double omx) {
    omx = omx / hartree_to_meV;
    this->omegax = omx;
    this->alpha_x = 1. / (omx * reducted_mass);
}

void QDParams::set_omega_y(double omy) {
    omy = omy / hartree_to_meV;
    this->omegay = omy;
    this->alpha_y = 1. / (omy * reducted_mass);
}
