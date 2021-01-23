#include "mymath.h"
#include <cmath>

#ifndef BESSEL_H
#define BESSEL_H

inline ldouble I0(ldouble x) {
    return ldouble(std::cyl_bessel_i(0, x));
}

inline ldouble I1(ldouble x) {
    return ldouble(std::cyl_bessel_i(1, x));
}

inline ldouble In(ldouble x, int n) {
    return ldouble(std::cyl_bessel_i(n, x));
}

/*
inline ldouble I0(ldouble x) {
    ldouble t_2, t_1, Bess, Bess1;

    if (abs(x) <= 3.75) {
        t_2 = pow(x / 3.75, 2.0);
        Bess = 0.0360768 + t_2 * 0.0045813;
        Bess = 0.2659732 + t_2 * Bess;
        Bess = 1.2067492 + t_2 * Bess;
        Bess = 3.0899424 + t_2 * Bess;
        Bess = 3.5156229 + t_2 * Bess;

        return 1.0 + t_2 * Bess;
    } else {
        t_1 = 3.75 / x;
        Bess1 = 1.0 / sqrtl(x);
        Bess1 = Bess1 * expl(x);
        Bess = -0.01647633 + t_1 * 0.00392377;
        Bess = 0.02635537 + t_1 * Bess;
        Bess = -0.02057706 + t_1 * Bess;
        Bess = 0.00916281 + t_1 * Bess;
        Bess = -0.00157565 + t_1 * Bess;
        Bess = 0.00225319 + t_1 * Bess;
        Bess = 0.01328592 + t_1 * Bess;
        Bess = 0.39894228 + t_1 * Bess;

        return Bess1 * Bess;
    }
}

inline ldouble I1(ldouble x) {
    ldouble t2, t_1, Bess, Bess1;

    if (abs(x) <= 3.75) {
        t2 = pow(x / 3.75, 2.0);
        Bess = 0.00301532 + t2 * 0.00032411;
        Bess = 0.02658733 + t2 * Bess;
        Bess = 0.15084934 + t2 * Bess;
        Bess = 0.51498869 + t2 * Bess;
        Bess = 0.87890594 + t2 * Bess;
        Bess = 1.0 / 2.0 + t2 * Bess;

        return x * Bess;
    } else {
        t_1 = 3.75 / x;
        Bess1 = 1.0 / sqrtl(x);
        Bess1 = Bess1 * expl(x);
        Bess = 0.01787654 - t_1 * 0.00420059;
        Bess = -0.02895312 + t_1 * Bess;
        Bess = 0.02282967 + t_1 * Bess;
        Bess = -0.01031555 + t_1 * Bess;
        Bess = 0.00163801 + t_1 * Bess;
        Bess = -0.00362018 + t_1 * Bess;
        Bess = -0.03988024 + t_1 * Bess;
        Bess = 0.39894228 + t_1 * Bess;

        return Bess1 * Bess;
    }
}

inline ldouble In(ldouble x, int m) {
    if (m == 0) return I0(x);
    if (m == 1) return I1(x);
    if (x == 0) return 0;

    int iacc = 40; // Increase to enhance accuracy

    const ldouble kBigPositive = 1.e10;
    const ldouble kBigNegative = 1.e-10;

    ldouble tox = 2.0 / abs(x);
    ldouble bip = 0;
    ldouble bim;
    ldouble bi = 1;
    ldouble result = 0;

    int mm = 2 * (m + (int)round(sqrt(iacc * m))); // May be round should be Integer
    for (int j = mm; j >= 1; j--) {
        bim = bip + j * tox * bi;
        bip = bi;
        bi = bim;

        // Renormalise to prevent overflows
        if (abs(bi) > kBigPositive) {
            result *= kBigNegative;
            bip    *= kBigNegative;
            bi     *= kBigNegative;
        }

        if (j == m) {
            result = bip;
        }
    }

    result *= I0(x) / bi; // Normalise with BesselI0(x)

    if ((x < 0) && (m % 2 == 1)) {
        result = -result;
    }

    return result;
}
*/

inline ldouble InD(ldouble x, int n) {
    if (n == 0) {
        return I1(x);
    }

    return (In(x, n - 1) + In(x, n + 1)) / ldouble(2.0);
}

#endif // BESSEL_H
