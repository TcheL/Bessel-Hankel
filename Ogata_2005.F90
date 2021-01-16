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

  integer, parameter :: nRoots = 1000
  real(kind = myKind), parameter :: pi = 3.141592653589793d0
  real(kind = myKind), parameter :: steph = 1.0d-7

  integer :: iExample, nSample, nu
  real(kind = myKind) :: rlogMin, rlogMax, rlogStep
  real(kind = myKind) :: t(nRoots), x(nRoots), xi(nRoots), &
    & psi(nRoots), psi_d(nRoots), omega(nRoots)

  integer :: fileID, i, j
  real(kind = myKind) :: r, numRes, anaRes, relErr
  character(len = :), allocatable :: fmtStr

  procedure(real(kind = myKind)), pointer :: pFunc2exampleFunction => null()
  namelist /input/ iExample, nSample, rlogMin, rlogMax

  ! ===== to prepare =====
  open(newunit = fileID, file = 'input.nml', status = 'old')
    read(fileID, nml = input)
  close(fileID)
  rlogStep = (rlogMax - rlogMin)/(nSample - 1)

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
    fmtStr = "(1P, 3(E16.7E3, 1X), 0P, F9.2, A)"
  else
    fmtStr = "(1P, 3(E25.16E3, 1X), 0P, F9.2, A)"
  end if

  ! ===== to calculate =====
  call jvn(nu, nRoots, xi)
  xi = xi/pi
  t = steph*xi

  ! equation (5.1) in (Ogata, 2005)
  psi = t*tanh(pi/2.0d0*sinh(t))
  x = pi/steph*psi

  psi_d = cosh(pi/2.0d0*sinh(t))
  psi_d = tanh(pi/2.0d0*sinh(t)) + t*(pi/2.0d0*cosh(t)/(psi_d*psi_d))

  omega = bessel_yn(nu, pi*xi)/bessel_jn(nu + 1, pi*xi)

  do i = 1, nSample, 1
    r = 10**(rlogMin + rlogStep*(i - 1))
    numRes = 0.0d0
    do j = 1, nRoots, 1
      ! modified equation (5.2) in (Ogata, 2005)
      numRes = numRes + 1.0d0/(r*r)*( pi*omega(j) &
        & *funcKernel(pFunc2exampleFunction, "lambda", x(j)/r) &
        & *bessel_jn(nu, x(j))*psi_d(j) )*x(j)
    end do
    anaRes = funcKernel(pFunc2exampleFunction, 'r', r)
    relErr = (numRes - anaRes)/anaRes*100.0d0
    write(*, fmtStr) r, numRes, anaRes, relErr, '%'
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
