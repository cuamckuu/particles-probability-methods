#ifndef VEC3_H
#define VEC3_H

#include <iostream>
#include <cmath>
#include "mymath.h"

class Vec3 {
public:
    ldouble x, y, z;

    Vec3();
    Vec3(const Vec3 &v); // Copy constructor
    Vec3(ldouble x, ldouble y, ldouble z);
    ldouble magnitudeSqr() const;
    ldouble magnitude() const;

    ldouble dot(Vec3 &v) const;
    Vec3 cross(Vec3 &v) const;
    Vec3 norm() const;
    Vec3 proj(Vec3 &v) const;

    ldouble comp(Vec3 &v) const;
    ldouble diff_angle(Vec3 &v) const;

    Vec3 rotate(ldouble theta, Vec3 normal) const;

    Vec3& operator=(const Vec3 &v); // Copy assigment
    Vec3 operator+(const Vec3 &v);
    Vec3 operator-(const Vec3 &v);

    Vec3 operator*(const ldouble a);
    Vec3 operator/(const ldouble a);

    Vec3& operator+=(const Vec3 &v);
    Vec3& operator-=(const Vec3 &v);
    Vec3& operator*=(const ldouble num);

    friend std::ostream& operator<< (std::ostream& stream, const Vec3& v);
};

#endif // VEC3_H
