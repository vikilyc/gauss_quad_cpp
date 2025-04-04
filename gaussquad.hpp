#ifndef GAUSSQUAD_HPP
#define GAUSSQUAD_HPP

#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <stdexcept>

enum RuleCode {
    LEGENDRE = 0,
    CHEBYSHEV_FIRST,
    CHEBYSHEV_SECOND,
    JACOBI,
    LAGUERRE,
    HERMITE
};

const int MAXITS = 30;
const double ZERO = 0.0;
const double HALF = 0.5;
const double ONE  = 1.0;
const double TWO  = 2.0;
const double FOUR = 4.0;
const double EPSILON = std::numeric_limits<double>::epsilon();

void gauss_rule(int icode, int n, std::vector<double>& t, std::vector<double>& w,
                std::vector<double>& work, double alpha, double beta,
                char endpt, int& info);

void recur_coeffs(int icode, int n, double alpha, double beta,
                  std::vector<double>& a, std::vector<double>& b,
                  double& muzero, int& info);

void custom_gauss_rule(int n, std::vector<double>& a, std::vector<double>& b,
                       std::vector<double>& w, double muzero, char endpt,
                       double lo, double hi, int& info);

void st_eigenproblem(int n, std::vector<double>& d, std::vector<double>& e,
                     std::vector<double>& z, int& info);

void orthonormal_polynomials(int n, int m, const std::vector<double>& x,
                             const std::vector<double>& a, const std::vector<double>& b,
                             double muzero, std::vector<std::vector<double>>& p);

double solve(int n, double shift, const std::vector<double>& a, const std::vector<double>& b);

#endif // GAUSSQUAD_HPP

