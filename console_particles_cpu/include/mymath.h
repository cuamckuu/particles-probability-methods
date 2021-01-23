#ifndef MYMATH_H
#define MYMATH_H

#include <cmath>

using ldouble = double;

const ldouble PI = std::acos(-1.0);
const ldouble E  = std::exp(1.0);

const ldouble light_speed = 2.99792458e8;
const ldouble c = light_speed;

const ldouble m_e = 9.10938356e-31;
const ldouble q_e = 1.6021766209e-19;

const ldouble eps0 = 8.854188e-12;
const ldouble K = 1 / (4 * PI * eps0);

#endif // MYMATH_H
