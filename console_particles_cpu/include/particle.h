#ifndef PARTICLE_H
#define PARTICLE_H

#include "vec3.h"
#include "mymath.h"

class Particle {
public:
    ldouble q = 0.0; // Заряд частицы
    ldouble W = 0.0; // Энергия частицы

    Vec3 pos;
    Vec3 pulse;

    void recountEnergy();
    Particle() = default;
    Particle(const ldouble q, const Vec3 &r, const Vec3 &p);
    Particle(const ldouble q, const ldouble W, const Vec3 &r, const Vec3 &p);
};

#endif // PARTICLE_H
