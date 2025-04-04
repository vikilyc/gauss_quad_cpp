#include "gaussquad.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

std::string rule_name(int code) {
    switch (code) {
        case LEGENDRE: return "Legendre";
        case CHEBYSHEV_FIRST: return "Chebyshev I";
        case CHEBYSHEV_SECOND: return "Chebyshev II";
        case JACOBI: return "Jacobi(alpha=0.5, beta=0.5)";
        case LAGUERRE: return "Laguerre(alpha=0.5)";
        case HERMITE: return "Hermite";
        default: return "Unknown";
    }
}

// Reference integrals of sin(x) * w(x)
double true_integral(int rule) {
    switch (rule) {
        case LEGENDRE: return 0.0;                         // ∫_{-1}^{1} sin(x) dx = 0
        case CHEBYSHEV_FIRST: return 0.0;                  // Symmetric odd function
        case CHEBYSHEV_SECOND: return 0.0;
        case JACOBI: return 0.0;
        case LAGUERRE: return 1.0 / (1.0 * 1.0 + 1.0);      // ∫₀^∞ sin(x) * e^(-x) dx ≈ 0.5
        case HERMITE: return 0.0;                          // Odd function with symmetric weight
        default: return 0.0;
    }
}

int main() {
    std::vector<int> ns = {5, 10, 20};
    std::vector<int> rules = {LEGENDRE, CHEBYSHEV_FIRST, CHEBYSHEV_SECOND, JACOBI, LAGUERRE, HERMITE};

    for (int rule : rules) {
        std::cout << "Rule: " << rule_name(rule) << std::endl;
        double true_val = true_integral(rule);
        for (int n : ns) {
            int info = 0;
            std::vector<double> t(n), w(n), work(n);
            double alpha = 0.5, beta = 0.5;
            char endpt = 'N';

            gauss_rule(rule, n, t, w, work, alpha, beta, endpt, info);

            if (info != 0) {
                std::cerr << "  [n=" << n << "] Error in gauss_rule: info = " << info << std::endl;
                continue;
            }

            double integral = 0.0;
            for (int i = 0; i < n; ++i) {
                integral += w[i] * std::sin(t[i]);
            }

            double err = std::abs(integral - true_val);
            std::cout << "  n = " << std::setw(2) << n
                      << ", approx = " << std::setw(16) << std::setprecision(10) << integral
                      << ", error = "  << std::setw(12) << err << std::endl;
        }
        std::cout << std::endl;
    }

    return 0;
}
