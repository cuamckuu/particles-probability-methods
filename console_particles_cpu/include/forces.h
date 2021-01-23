#ifndef FORCES_H
#define FORCES_H

#include <vector>

#include "vec3.h"
#include "particle.h"

using ldoubleArr2D = std::vector<std::vector<ldouble>>;

Vec3 ConvertToXYZ(const Vec3 &v, const ldouble teta);

Vec3 F_Coulomb(const Particle &p1, const Particle &p2);

Vec3 F_Cherenkov_2(const Particle & p1, const Particle & p2, const int nu_, const int n_, const ldoubleArr2D &k, const ldoubleArr2D &Amp);

#endif // FORCES_H
