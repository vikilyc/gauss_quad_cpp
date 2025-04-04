  ! This module is a revision of the Fortran 77 code gaussq.f
  !
  !   Original version 20 Jan 1975 from Stanford
  !   Modified 21 Dec 1983 by Eric Grosse
  !   f95 version Nov 2005 by Bill McLean
  !   c++ version Apr 2025 by ChatGPT and Yuchen Liu
  !
  ! The routines are used to compute the nodes t(j) and weights
  ! w(j) for Gaussian-type quadrature rules with pre-assigned
  ! nodes.  These quantities are used when one wishes to approximate
  !
  !                 / b
  !                 |
  !                 | f(x) w(x) dx
  !                 |
  !                 / a
  !
  ! by
  !
  !                  n
  !                 ---
  !                 \
  !                  |  f(t(j)) * w(j).
  !                 /
  !                 ---
  !                 j=1
  !
  ! (Note w(x) and w(j) have no connection with each other.)
  ! Here w(x) is one of six possible non-negative weight
  ! functions (listed below), and f(x) is the
  ! function to be integrated.  Gaussian quadrature is particularly
  ! useful on infinite intervals (with appropriate weight
  ! functions), since then other techniques often fail.
  !
  ! Associated with each weight function w(x) is a set of
  ! orthogonal polynomials.  The nodes t(j) are just the zeroes
  ! of the proper n-th degree polynomial.
  !
  ! References:
  !
  !      1.  Golub, G. H., and Welsch, J. H., Calculation of gaussian
  !          quadrature rules, Mathematics of Computation 23 (april,
  !          1969), pp. 221-230.
  !      2.  Golub, G. H., Some modified matrix eigenvalue problems,
  !          Siam Review 15 (april, 1973), pp. 318-334 (section 7).
  !      3.  Stroud and Secrest, Gaussian Quadrature Formulas, Prentice-
  !          Hall, Englewood Cliffs, N.J., 1966.
