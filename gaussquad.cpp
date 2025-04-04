#include "gaussquad.hpp"
#include <iostream>
#include <algorithm>
//#include <omp.h>

double solve(int n, double shift, const std::vector<double>& a, const std::vector<double>& b) {
    double t = a[0] - shift;
    for (int i = 1; i < n - 1; ++i) {
        t = a[i] - shift - (b[i - 1] * b[i - 1]) / t;
    }
    return ONE / t;
}

void custom_gauss_rule(int n, std::vector<double>& a, std::vector<double>& b, std::vector<double>& w,
                       double muzero, char endpt, double lo, double hi, int& info) {
    double g, t1;
    if (endpt == 'L') {
        a[n - 1] = (n == 1) ? lo : solve(n, lo, a, b) * b[n - 2] * b[n - 2] + lo;
    } else if (endpt == 'R') {
        a[n - 1] = (n == 1) ? hi : solve(n, hi, a, b) * b[n - 2] * b[n - 2] + hi;
    } else if (endpt == 'B') {
        if (n == 1) {
            info = -6;
            return;
        }
        g = solve(n, lo, a, b);
        t1 = (lo - hi) / (solve(n, hi, a, b) - g);
        b[n - 2] = std::sqrt(t1);
        a[n - 1] = lo + g * t1;
    }

    st_eigenproblem(n, a, b, w, info);
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        w[i] = muzero * w[i] * w[i];
    }
}

void recur_coeffs(int icode, int n, double alpha, double beta,
                  std::vector<double>& a, std::vector<double>& b,
                  double& muzero, int& info) {
    info = 0;
    double pi = FOUR * std::atan(ONE);

    switch (icode) {
        case LEGENDRE:
            muzero = TWO;
            for (int i = 0; i < n; ++i) {
                a[i] = ZERO;
                double ri = i + 1;
                b[i] = ri / std::sqrt(FOUR * ri * ri - ONE);
            }
            break;
        case CHEBYSHEV_FIRST:
            muzero = pi;
            for (int i = 0; i < n; ++i) {
                a[i] = ZERO;
                b[i] = HALF;
            }
            b[0] = std::sqrt(HALF);
            break;
        case CHEBYSHEV_SECOND:
            muzero = HALF * pi;
            for (int i = 0; i < n; ++i) {
                a[i] = ZERO;
                b[i] = HALF;
            }
            break;
        case HERMITE:
            muzero = std::sqrt(pi);
            for (int i = 0; i < n; ++i) {
                a[i] = ZERO;
                b[i] = std::sqrt(HALF * (i + 1));
            }
            break;
        case JACOBI: {
            if (alpha <= -ONE || beta <= -ONE) {
                info = -1;
                return;
            }
            double ab = alpha + beta;
            double abi = TWO + ab;
            muzero = std::pow(TWO, ab + ONE) *
                     std::tgamma(alpha + ONE) * std::tgamma(beta + ONE) / std::tgamma(abi);
            a[0] = (beta - alpha) / abi;
            b[0] = std::sqrt(FOUR * (ONE + alpha) * (ONE + beta) /
                            ((abi + ONE) * abi * abi));
            double a2b2 = beta * beta - alpha * alpha;
            for (int i = 1; i < n; ++i) {
                abi = TWO * (i + 1) + ab;
                a[i] = a2b2 / ((abi - TWO) * abi);
                b[i] = std::sqrt(FOUR * (i + 1) * (i + 1 + alpha) *
                                 (i + 1 + beta) * (i + 1 + ab) /
                                 ((abi * abi - ONE) * abi * abi));
            }
            break;
        }
        case LAGUERRE:
            if (alpha <= -ONE) {
                info = -1;
                return;
            }
            muzero = std::tgamma(alpha + ONE);
            for (int i = 0; i < n; ++i) {
                a[i] = TWO * (i + 1) - ONE + alpha;
                b[i] = std::sqrt((i + 1) * (i + 1 + alpha));
            }
            break;
        default:
            info = -1;
    }
}

void st_eigenproblem(int n, std::vector<double>& d, std::vector<double>& e,
                     std::vector<double>& z, int& info) {
    info = 0;
    z[0] = ONE;
    std::fill(z.begin() + 1, z.end(), ZERO);
    if (n == 1) return;

    e[n - 1] = ZERO;
    for (int l = 0; l < n; ++l) {
        for (int j = 0; j < MAXITS; ++j) {
            int m = l;
            while (m < n - 1 && std::abs(e[m]) > EPSILON * (std::abs(d[m]) + std::abs(d[m + 1])))
                ++m;
            if (m == l) break;
            if (j == MAXITS - 1) {
                info = l + 1;
                return;
            }

            double g = (d[l + 1] - d[l]) / (TWO * e[l]);
            double r = std::sqrt(g * g + ONE);
            g = d[m] - d[l] + e[l] / (g + std::copysign(r, g));
            double s = ONE, c = ONE, p = ZERO;

            for (int i = m - 1; i >= l; --i) {
                double f = s * e[i];
                double b = c * e[i];
                if (std::abs(f) < std::abs(g)) {
                    s = f / g;
                    r = std::sqrt(s * s + ONE);
                    e[i + 1] = g * r;
                    c = ONE / r;
                    s *= c;
                } else {
                    c = g / f;
                    r = std::sqrt(c * c + ONE);
                    e[i + 1] = f * r;
                    s = ONE / r;
                    c *= s;
                }
                g = d[i + 1] - p;
                r = (d[i] - g) * s + TWO * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;

                double temp = z[i + 1];
                z[i + 1] = s * z[i] + c * temp;
                z[i] = c * z[i] - s * temp;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = ZERO;
        }
    }

    for (int i = 0; i < n - 1; ++i) {
        int k = i;
        double p = d[i];
        for (int j = i + 1; j < n; ++j) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            std::swap(d[i], d[k]);
            std::swap(z[i], z[k]);
        }
    }
}

void orthonormal_polynomials(int n, int m, const std::vector<double>& x,
                             const std::vector<double>& a, const std::vector<double>& b,
                             double muzero, std::vector<std::vector<double>>& p) {
    double c = ONE / std::sqrt(muzero);
    double rb = ONE / b[0];

    #pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        p[i][0] = c;
        p[i][1] = rb * (x[i] - a[0]) * c;
    }

    for (int j = 2; j <= n; ++j) {
        double rbj = ONE / b[j - 1];
        #pragma omp parallel for
        for (int i = 0; i < m; ++i) {
            p[i][j] = rbj * ((x[i] - a[j - 1]) * p[i][j - 1] - b[j - 2] * p[i][j - 2]);
        }
    }
}

void gauss_rule(int icode, int n, std::vector<double>& t, std::vector<double>& w,
                std::vector<double>& work, double alpha, double beta,
                char endpt, int& info) {
    info = 0;
    double muzero, lo = -ONE, hi = ONE;

    if (n < 1) { info = -2; return; }
    if ((icode == JACOBI || icode == LAGUERRE) && alpha <= -ONE) { info = -6; return; }
    if (icode == JACOBI && beta <= -ONE) { info = -7; return; }
    if (!(endpt == 'N' || endpt == 'B' || endpt == 'L' || endpt == 'R')) { info = -8; return; }

    recur_coeffs(icode, n, alpha, beta, t, work, muzero, info);
    if (info != 0) return;

    if (n == 1 && endpt == 'B') { info = -8; return; }

    if (icode == LAGUERRE && endpt == 'L') lo = ZERO;
    if ((icode == LAGUERRE && (endpt == 'R' || endpt == 'B')) ||
        (icode == HERMITE && endpt != 'N')) {
        info = -8;
        return;
    }

    custom_gauss_rule(n, t, work, w, muzero, endpt, lo, hi, info);
}

