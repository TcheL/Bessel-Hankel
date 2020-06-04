#define KV 8

! References:
!     Ogata, 2005. A numerical integration formula based on the Bessel
!   functions. Publ. RIMS, 41, 949-970.
!     Guptasarma D. and Singh, 1997. New digital linear filters for Hankel J_0
!   and J_1 tranforms. Geophysical Prospecting, 45, 745-762.
!     shixun22, 2020. hankel (Python code for fast computation of Hankel
!   transforms using Ogata 2005 method). https://github.com/shixun22/hankel

program main
  use Bessel_Function
  implicit none

  ! Hankel transform, https://en.wikipedia.org/wiki/Hankel_transform#Definition
  !   $ F(r) = \int_0^\infty f(\lambda) J_nu(r \lambda) \lambda d\lambda $

  integer, parameter :: nRoots = 500, nSample = 400
  real(kind = KV), parameter :: pi = 3.141592653589793d0, steph = 1.0d-5, &
    & delta = 1.0d-3
  real(kind = KV) :: t(nRoots), x(nRoots), xi(nRoots), &
    & psi(nRoots), psi_d(nRoots), omega(nRoots)

  integer :: i, j, nu = 0
  real(kind = KV) :: r, res, funcKernel

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
    r = i*delta
    res = 0.0d0
    do j = 1, nRoots, 1
      ! modified equation (5.2) in (Ogata, 2005)
      res = res + 1.0d0/(r*r)*(pi*omega(j)*funcKernel("lambda", x(j)/r) &
        & *bessel_jn(nu, x(j))*psi_d(j))*x(j)
    end do
    print*, i, res, funcKernel('r', r)
  end do
end program

real(kind = KV) function funcKernel(varName, var) result(f)
  ! equation (6) or (8) in (Guptasarma and Singh, 1997)
  character(len = *), intent(in) :: varName
  real(kind = KV), intent(in) :: var
  if(varName == "lambda") then
    f = exp( - var )
  else if(varName == 'r') then
    f = 1.0d0/sqrt((1.0d0 + var*var)**3) ! for nu = 0
    ! f = var/sqrt((1.0d0 + var*var)**3) ! for nu = 1
  end if
end function funcKernel
