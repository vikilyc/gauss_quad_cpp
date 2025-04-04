# Gaussian Quadrature Module

This project provides a modern and parallelized implementation of Gaussian-type quadrature rule generation, adapted from historic Fortran 77 code and now available in both modern Fortran and C++ with OpenMP.

---

## 📜 History

- **Original Fortran 77**: January 20, 1975, Stanford University  
- **Modified**: December 21, 1983, by Eric Grosse  
- **Fortran 95 Version**: November 2005, by Bill McLean  
- **C++ Version**: April 2025, by ChatGPT and Yuchen Liu

---

## 🎯 Purpose

This module computes the **nodes** `t(j)` and **weights** `w(j)` for Gaussian quadrature rules with pre-assigned weight functions, used to approximate integrals of the form:

\\[
\\int_a^b f(x) w(x) dx \\approx \\sum_{j=1}^{n} f(t_j) w_j
\\]

> Note: `w(x)` (the continuous weight function) and `w_j` (discrete quadrature weights) are unrelated.

Gaussian quadrature is especially effective on infinite intervals with appropriate weight functions, where standard numerical methods may fail.

---

## 🔢 Supported Weight Functions

| Code | Weight Function \( w(x) \)                          | Interval         | Orthogonal Polynomials |
|------|------------------------------------------------------|------------------|------------------------|
| 0    | \( 1 \)                                              | \((-1, 1)\)      | Legendre              |
| 1    | \( \\frac{1}{\\sqrt{1 - x^2}} \)                     | \((-1, 1)\)      | Chebyshev (First Kind)|
| 2    | \( \\sqrt{1 - x^2} \)                                | \((-1, 1)\)      | Chebyshev (Second Kind)|
| 3    | \( (1 - x)^\\alpha (1 + x)^\\beta \)                | \((-1, 1)\)      | Jacobi                |
| 4    | \( x^\\alpha e^{-x} \)                               | \((0, \\infty)\) | Laguerre              |
| 5    | \( e^{-x^2} \)                                       | \((-\\infty, \\infty)\) | Hermite       |

---

## 📚 References

1. **Golub, G. H., and Welsch, J. H.**  
   _Calculation of Gaussian Quadrature Rules_,  
   Mathematics of Computation, Vol. 23, April 1969, pp. 221–230.

2. **Golub, G. H.**  
   _Some Modified Matrix Eigenvalue Problems_,  
   SIAM Review, Vol. 15, April 1973, pp. 318–334.

3. **Stroud and Secrest**  
   _Gaussian Quadrature Formulas_,  
   Prentice-Hall, Englewood Cliffs, N.J., 1966.
