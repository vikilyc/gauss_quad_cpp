program test_sin
  use gaussquad_module
  implicit none

  integer :: i, n, icode, info
  real(kind=WP), allocatable :: t(:), w(:), work(:)
  real(kind=WP) :: alpha, beta, result, exact
  character(len=1) :: endpt
  character(len=20) :: rulename

  do icode = 0, 5
    select case(icode)
    case(0)
      rulename = "Legendre"
      exact = 0.0_WP
    case(1)
      rulename = "Chebyshev I"
      exact = 0.0_WP
    case(2)
      rulename = "Chebyshev II"
      exact = 0.0_WP
    case(3)
      rulename = "Jacobi"
      exact = 0.0_WP
    case(4)
      rulename = "Laguerre"
      exact = 0.5_WP
    case(5)
      rulename = "Hermite"
      exact = 0.0_WP
    end select

    print *, "Rule: ", trim(rulename)
    do n = 5, 20, 5
      allocate(t(n), w(n), work(n))
      alpha = 0.5_WP
      beta  = 0.5_WP
      endpt = 'N'

      call gauss_rule(icode, n, t, w, work, alpha, beta, endpt, info)
      if (info /= 0) then
        print *, "  Error in gauss_rule, info = ", info
        cycle
      end if

      result = 0.0_WP
      do i = 1, n
        result = result + w(i) * sin(t(i))
      end do

      print '(a, i2, a, f18.12, a, es12.4)', "  n = ", n, ": approx = ", result, ", error = ", abs(result - exact)

      deallocate(t, w, work)
    end do
    print *
  end do

end program test_sin
