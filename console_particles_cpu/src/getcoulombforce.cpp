#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <cmath>
#include <vector>
#include <array>
#include <cstdlib>
#include <iomanip>
#include <limits>

#include "forces.h"
#include "vec3.h"
#include "particle.h"

std::vector<Particle> read_particles_from_file(std::string filename) {
    std::vector<Particle> particles;

    std::ifstream fin(filename);
    fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip headers

    while (fin) {
        ldouble q;
        ldouble W;
        ldouble x, y, z;
        ldouble px, py, pz;

        fin >> q >> W >> x >> y >> z >> px >> py >> pz;

        if(!fin)
            break;

        particles.emplace_back(
            q,
            W,
            Vec3(x, y, z),
            Vec3(px, py, pz)
        );
    }

    return particles;
}

Vec3 get_coulomb_force(Particle p1, Particle p2) {
    ldouble gamma = (p1.W / (m_e * std::pow(c, 2)));

    ldouble r0 = std::sqrt(std::pow(p1.pos.x, 2) + std::pow(p1.pos.y, 2));
    ldouble teta0 = 0;
    if(r0 > 1e-20) {
        if(p1.pos.y >= 0) {
            teta0 = std::acos(p1.pos.x / r0); // Angle from [0; pi]
        } else {
            teta0 = 2*PI - std::acos(p1.pos.x / r0);
        }
    }

    ldouble r = std::sqrt(std::pow(p2.pos.x, 2) + std::pow(p2.pos.y, 2));
    ldouble teta = 0;
    if(r > 1e-20) {
        if(p2.pos.y >= 0) {
            teta = std::acos(p2.pos.x / r); // Angle from [0; pi]
        } else {
            teta = 2*PI - std::acos(p2.pos.x / r);
        }
    }

    ldouble dteta_sin = std::sin(teta - teta0);
    ldouble dteta_cos = std::cos(teta - teta0);

    ldouble z_hat = ((p1.pos.z - p2.pos.z) * gamma);
    ldouble denominator = std::pow(
        r*r + r0*r0 - 2*r*r0*dteta_cos + z_hat*z_hat,
        3.0/2.0
    );

    if(denominator < 1e-20)
        std::cout << "p1 " << p1.pos << " p2 " << p2.pos << denominator << std::endl;

    ldouble Q = (p1.q * p2.q);
    ldouble F_r    = (Q * (r - r0 * dteta_cos)) / (gamma * denominator);
    ldouble F_teta = (Q * (r0 * dteta_sin)) / (gamma * denominator);

    ldouble teta_sin = std::sin(teta);
    ldouble teta_cos = std::cos(teta);

    ldouble Fx = K * (F_r * teta_cos - F_teta * teta_sin);
    ldouble Fy = K * (F_r * teta_sin + F_teta * teta_cos);
    ldouble Fz = K * (Q * z_hat / denominator);

    return Vec3(Fx, Fy, Fz);
}

void save_vec_to_file(std::string filename, Vec3& vec) {
    std::ofstream fout(filename);

    fout << std::scientific;
    fout.precision(std::numeric_limits<double>::digits10+2);

    fout << "x y z" << std::endl;
    fout << vec.x << " " << vec.y << " " << vec.z << std::endl;
}

int main() {
    std::vector<Particle> particles = read_particles_from_file("data.csv");

    Vec3 force = get_coulomb_force(particles[0], particles[1]);

    save_vec_to_file("res.csv", force);

    return 0;
}
