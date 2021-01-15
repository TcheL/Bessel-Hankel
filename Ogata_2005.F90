! References:
!     Ogata, 2005. A numerical integration formula based on the Bessel
!   functions. Publ. RIMS, 41, 949-970.
!     shixun22, 2020. hankel (Python code for fast computation of Hankel
!   transforms using Ogata 2005 method). https://github.com/shixun22/hankel

program main
  use example_Functions
  use Bessel_Function
  implicit none

  ! Hankel transform, https://en.wikipedia.org/wiki/Hankel_transform#Definition
  !   $ F(r) = \int_0^\infty f(\lambda) J_\nu(r \lambda) \lambda d\lambda $

  integer, parameter :: nRoots = 500, nSample = 400
  real(kind = myKind), parameter :: pi = 3.141592653589793_myKind, &
    & steph = 1.0d-5, delta = 1.0d-3
  real(kind = myKind) :: t(nRoots), x(nRoots), xi(nRoots), &
    & psi(nRoots), psi_d(nRoots), omega(nRoots)

  integer :: iExample, i, j, nu = 0
  real(kind = myKind) :: r, numRes, anaRes, relatError
  character(len = :), allocatable :: fmtStr

  procedure(real(kind = myKind)), pointer :: pFunc2exampleFunction => null()

  ! ===== to prepare =====
  write(*, "(A)", advance = "NO") "> which example you want to check (4-10)? "
  read(*, *) iExample

  select case(iExample)
    case(4)
      pFunc2exampleFunction => exampleFunc1
      nu = 0
    case(5)
      pFunc2exampleFunction => exampleFunc2
      nu = 0
    case(6)
      pFunc2exampleFunction => exampleFunc3
      nu = 0
    case(7)
      pFunc2exampleFunction => exampleFunc4
      nu = 1
    case(8)
      pFunc2exampleFunction => exampleFunc5
      nu = 1
    case(9)
      pFunc2exampleFunction => exampleFunc6
      nu = 1
    case(10)
      pFunc2exampleFunction => exampleFunc7
      nu = 1
    case default
      stop "<<< Illegal input option, please have a check."
  end select

  if(myKind == kind(1.0e0)) then
    fmtStr = "(I5, 1P, 1X, E14.7, 1X, E14.7, 1X, 0P, F9.2, A)"
  else
    fmtStr = "(I5, 1P, 1X, E23.16, 1X, E23.16, 1X, 0P, F9.2, A)"
  end if

  ! ===== to calculate =====
  call jvn(nu, nRoots, xi)
  xi = xi/pi
  t = steph*xi

  ! equation (5.1) in (Ogata, 2005)
  psi = t*tanh(pi/2.0_myKind*sinh(t))
  x = pi/steph*psi

  psi_d = cosh(pi/2.0_myKind*sinh(t))
  psi_d = tanh(pi/2.0_myKind*sinh(t)) + t*(pi/2.0_myKind*cosh(t)/(psi_d*psi_d))

  omega = bessel_yn(nu, pi*xi)/bessel_jn(nu + 1, pi*xi)

  do i = 1, nSample, 1
    r = i*delta
    numRes = 0.0_myKind
    do j = 1, nRoots, 1
      ! modified equation (5.2) in (Ogata, 2005)
      numRes = numRes + 1.0_myKind/(r*r)*( pi*omega(j) &
        & *funcKernel(pFunc2exampleFunction, "lambda", x(j)/r) &
        & *bessel_jn(nu, x(j))*psi_d(j) )*x(j)
    end do
    anaRes = funcKernel(pFunc2exampleFunction, 'r', r)
    relatError = (numRes - anaRes)/anaRes*100.0_myKind
    write(*, fmtStr) i, numRes, anaRes, relatError, '%'
  end do

  contains

    real(kind = myKind) function funcKernel(pFunc, varName, var) result(r)
      real(kind = myKind) :: pFunc
      character(len = *), intent(in) :: varName
      real(kind = myKind), intent(in) :: var
      if(varName == "lambda") then
        r = pFunc(varName, var)/var
      else if(varName == 'r') then
        r = pFunc(varName, var)
      end if
    end function funcKernel

end program
