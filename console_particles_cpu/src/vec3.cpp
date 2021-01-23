#include "vec3.h"

Vec3::Vec3()
    :x(0.0)
    ,y(0.0)
    ,z(0.0)
{}

Vec3::Vec3(const Vec3 &v)
    :x(v.x)
    ,y(v.y)
    ,z(v.z)
{}

Vec3::Vec3(ldouble x, ldouble y, ldouble z)
    :x(x)
    ,y(y)
    ,z(z)
{}

ldouble Vec3::magnitudeSqr() const {
    return std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2);
}

ldouble Vec3::magnitude() const {
    return std::sqrt(magnitudeSqr());
}

ldouble Vec3::dot(Vec3 &v) const {
    return (x * v.x) + (y * v.y) + (z * v.z);
}

Vec3 Vec3::cross(Vec3 &v) const {
    Vec3 cross;
    cross.x = (y * v.z) - (z * v.y);
    cross.y = (z * v.x) - (x * v.z);
    cross.z = (x * v.y) - (y * v.x);

    return cross;
}

Vec3 Vec3::norm() const {
    ldouble mag = magnitude();

    Vec3 norm(x/mag, y/mag, z/mag);

    return norm;
}

Vec3 Vec3::proj(Vec3 &v) const {
    Vec3 v_n = v.norm();
    ldouble dotPr = dot(v_n);

    Vec3 proj(dotPr*v_n.x, dotPr*v_n.y, dotPr*v_n.z);

    return proj;
}

ldouble Vec3::comp(Vec3 &v) const {
    Vec3 v_n = v.norm();
    return dot(v_n);
}

ldouble Vec3::diff_angle(Vec3 &v) const {
    return acos(dot(v) / (magnitude() * v.magnitude()));
}

Vec3 Vec3::rotate(ldouble theta, Vec3 normal) const {
    ldouble mod = normal.magnitude();
    ldouble modSqr = normal.magnitudeSqr();

    ldouble cos_theta = cos(theta);
    ldouble sin_theta = sin(theta);

    Vec3 rotate;
    rotate.x = (cos_theta + (1 - cos_theta)*pow(normal.x / mod, 2))*x
        + ((1 - cos_theta)*normal.x*normal.y / modSqr - sin_theta * normal.z / mod)*y
        + ((1 - cos_theta)*normal.x*normal.z + sin_theta * normal.y / mod)*z;

    rotate.y = ((1 - cos_theta)*normal.x*normal.y / modSqr + sin_theta * normal.z / mod)*x
        + (cos_theta + (1 - cos_theta)*pow(normal.y / mod, 2))*y
        + ((1 - cos_theta)*normal.y*normal.z / modSqr - sin_theta * normal.x / mod)*z;

    rotate.z = ((1 - cos_theta)*normal.x*normal.z / modSqr - sin_theta * normal.y / mod)*x
        + ((1 - cos_theta)*normal.y*normal.z / modSqr + sin_theta * normal.x / mod)*y
        + (cos_theta + (1 - cos_theta)*pow(normal.z / mod, 2))*z;

    return rotate;
}

Vec3& Vec3::operator=(const Vec3 &v) {
    x = v.x;
    y = v.y;
    z = v.z;

    return *this;
} // Copy assigment

Vec3 Vec3::operator+(const Vec3 &v) {
    Vec3 sum(*this);

    sum.x += v.x;
    sum.y += v.y;
    sum.z += v.z;

    return sum;
}

Vec3 Vec3::operator-(const Vec3 &v) {
    Vec3 diff(*this);

    diff.x -= v.x;
    diff.y -= v.y;
    diff.z -= v.z;

    return diff;
}

Vec3 Vec3::operator*(const ldouble a) {
    Vec3 mul(*this);

    mul.x *= a;
    mul.y *= a;
    mul.z *= a;

    return mul;
}

Vec3 Vec3::operator/(const ldouble a) {
    Vec3 div(*this);

    div.x /= a;
    div.y /= a;
    div.z /= a;

    return div;
}

Vec3& Vec3::operator+=(const Vec3 &v) {
    x += v.x;
    y += v.y;
    z += v.z;

    return *this;
}

Vec3& Vec3::operator-=(const Vec3 &v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;

    return *this;
}

Vec3& Vec3::operator*=(const ldouble num) {
    x *= num;
    y *= num;
    z *= num;

    return *this;
}

std::ostream& operator<< (std::ostream& stream, const Vec3& v) {
    stream << "Vec3([" << v.x << ", " << v.y << ", " << v.z << "])";

    return stream;
}
