#include "particle.h"

Particle::Particle(const ldouble q, const Vec3 &pos, const Vec3 &pulse)
    :q(q)
    ,pos(pos)
    ,pulse(pulse)
{
    recountEnergy();
}

Particle::Particle(const ldouble q, const ldouble W, const Vec3 &pos, const Vec3 &pulse)
    :q(q)
    ,W(W)
    ,pos(pos)
    ,pulse(pulse)
{}

void Particle::recountEnergy() {
    W = light_speed * std::sqrt(pulse.magnitudeSqr() + m_e * m_e * light_speed * light_speed);
}
