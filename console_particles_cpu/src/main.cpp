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

void save_particles_to_file(std::string filename, const std::vector<Particle> &particles) {
    std::ofstream fout(filename);

    fout << "q w x y z px py pz" << std::endl;

    for (size_t i = 0; i < particles.size(); i++) {
        const Particle &p = particles[i];

        fout << std::scientific;
        fout.precision(std::numeric_limits<double>::digits10+2);

        fout << p.q << " " << p.W \
             << " " << p.pos.x << " " << p.pos.y << " " << p.pos.z \
             << " " << p.pulse.x << " " << p.pulse.y << " " << p.pulse.z \
             << std::endl;
    }
}

int main(int argc, char *argv[]) {
    if(argc != 2) {
        std::cout << "Usage: " << argv[0] << " Rc dt" << std::endl;
    }
    ldouble Rc = std::atof(argv[1]);
    ldouble dt = std::atof(argv[2]);

    std::vector<Particle> particles = read_particles_from_file("data.csv");

    std::vector<Vec3> forces(particles.size());

    for (size_t i = 0; i < particles.size(); i++) {
        const Particle &p1 = particles[i];

        bool p1_out_of_bound = (
            std::pow(p1.pos.x, 2) + std::pow(p1.pos.y, 2) > Rc*Rc
        );
        if (p1_out_of_bound) {
            continue;
        }

        for (size_t j = i+1; j < particles.size(); j++) {
            const Particle &p2 = particles[j];

            bool p2_out_of_bound = (
                std::pow(p2.pos.x, 2) + std::pow(p2.pos.y, 2) > Rc*Rc
            );
            if (p2_out_of_bound) {
                continue;
            }

            Vec3 force = get_coulomb_force(p1, p2);
            forces[i] += force;
            forces[j] -= force;
        }
    }

    for (size_t i = 0; i < particles.size(); i++) {
        particles[i].pulse += forces[i] * dt;
        particles[i].pos += particles[i].pulse * (c*c / particles[i].W) * dt;
    }

    save_particles_to_file("res.csv", particles);

    return 0;
}
