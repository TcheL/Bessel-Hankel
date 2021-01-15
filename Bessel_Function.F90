module Bessel_Function

  ! Reference:
  !   mjyzo.f90, J-P Moreau, Paris, (http://www.jpmoreau.fr):
  !     http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/mjyzo_f90.txt
  !   Properties of Bessel function:
  !     $ Z_v'(x) = \frac{v}{x} Z_v(x) - Z_{v+1}(x) $
  !     $ Z_v''(x) = (\frac{v^2}{x^2} - 1) Z_v(x) - \frac{1}{x} Z_v'(x) $
  !   First zero of Bessel function for large orders:
  !     equations (9.5.14-17) on page 371 of "Handbook of Mathematical
  !     Functions with Formulas, Graphs, and Mathematical Tables" (Milton
  !     Abramowitz and Irene A. Stegun, 1970)
  !   Fixed point iteration to solve zeros:
  !     https://en.wikipedia.org/wiki/Fixed-point_iteration#Applications
  !   Halley's method to solve zeros:
  !     https://en.wikipedia.org/wiki/Halley's_method#Method

  implicit none
  public

  integer, parameter, private :: KV = kind(1.0d0)
  real(kind = KV), parameter :: eps = 1.0d-9

  contains

#ifdef MOREAU

    real(kind = KV) function Jv(nu, x, of, obs, osu, obj) result(J)
      ! Bessel function of the 1st kind
      integer, intent(in) :: nu
      real(kind = KV), intent(in) :: x
      real(kind = KV), intent(out), optional :: of, obs, osu, obj(nu + 2)
      integer :: nt, mt, k
      real(kind = KV) :: f0, f1, f, bs, su, bj(nu + 2)
      do nt = 1, 900, 1
        mt = int(0.5d0*log10(6.28d0*nt) - nt*log10(1.36d0*abs(x)/nt))
        if(mt > 20) exit
      end do
      bs = 0.0d0
      f0 = 0.0d0
      f1 = tiny(1.0d0)
      ! Original: f1 = 1.0d-35
      su = 0.0d0
      bj = 0.0d0
      do k = nt, 0, - 1
        f = 2.0d0*(k + 1.0d0)*f1/x - f0
        if(k <= nu + 1) bj(k + 1) = f
        if(mod(k, 2) == 0) then
          bs = bs + 2.0d0*f
          if(k /= 0) su = su + ( - 1)**(k/2)*f/k
        end if
        f0 = f1
        f1 = f
      end do
      do k = 0, nu + 1, 1
        bj(k + 1) = bj(k + 1)/(bs - f)
      end do
      J = bj(nu + 1)
      if(present(of )) of  = f
      if(present(obs)) obs = bs
      if(present(osu)) osu = su
      if(present(obj)) obj = bj
    end function Jv

    real(kind = KV) function Yv(nu, x, oby) result(Y)
      ! Bessel function of the 2nd kind
      integer, intent(in) :: nu
      real(kind = KV), intent(in) :: x
      real(kind = KV), intent(out), optional :: oby(nu + 2)
      integer :: k
      real(kind = KV) :: ec, e0, s1, f0, f1, by(nu + 2)
      real(kind = KV) :: J, f, su, bs, bj(nu + 2)
      J = Jv(nu, x, of = f, obs = bs, osu = su, obj = bj)
      ec = 0.5772156649015329d0
      e0 = 0.3183098861837907d0
      s1 = 2.0d0*e0*(log(x/2.0d0)+ ec)*bj(1)
      f0 = s1 - 8.0d0*e0*su/(bs - f)
      f1 = (bj(2)*f0 - 2.0d0*e0/x)/bj(1)
      by(1) = f0
      by(2) = f1
      do k = 2, nu + 1, 1
        f = 2.0d0*(k - 1.0d0)*f1/x - f0
        by(k + 1) = f
        f0 = f1
        f1 = f
      end do
      Y = by(nu + 1)
      if(present(oby)) oby = by
    end function Yv

    real(kind = KV) function Jv1(nu, x, oJ) result(J1)
      ! The derivative of the 1st order of Bessel function of the 1st kind
      integer, intent(in) :: nu
      real(kind = KV), intent(in) :: x
      real(kind = KV), intent(out), optional :: oJ
      real(kind = KV) :: J, bj(nu + 2)
      J = Jv(nu, x, obj = bj)
      J1 = - bj(nu + 2) + nu*bj(nu + 1)/x
      if(present(oJ)) oJ = J
    end function Jv1

    real(kind = KV) function Yv1(nu, x, oY) result(Y1)
      ! The derivative of the 1st order of Bessel function of the 2nd kind
      integer, intent(in) :: nu
      real(kind = KV), intent(in) :: x
      real(kind = KV), intent(out), optional :: oY
      real(kind = KV) :: Y, by(nu + 2)
      Y = Yv(nu, x, oby = by)
      Y1 = - by(nu + 2) + nu*by(nu + 1)/x
      if(present(oY)) oY = Y
    end function Yv1

#else

    real(kind = KV) function Jv(nu, x) result(J)
      ! Bessel function of the 1st kind
      integer, intent(in) :: nu
      real(kind = KV), intent(in) :: x
      J = bessel_jn(nu, x)
    end function Jv

    real(kind = KV) function Yv(nu, x) result(Y)
      ! Bessel function of the 2nd kind
      integer, intent(in) :: nu
      real(kind = KV), intent(in) :: x
      Y = bessel_yn(nu, x)
    end function Yv

    real(kind = KV) function Jv1(nu, x, oJ) result(J1)
      ! The derivative of the 1st order of Bessel function of the 1st kind
      integer, intent(in) :: nu
      real(kind = KV), intent(in) :: x
      real(kind = KV), intent(out), optional :: oJ
      real(kind = KV) :: J
      J = bessel_jn(nu, x)
      J1 = nu/x*J - bessel_jn(nu + 1, x)
      if(present(oJ)) oJ = J
    end function Jv1

    real(kind = KV) function Yv1(nu, x, oY) result(Y1)
      ! The derivative of the 1st order of Bessel function of the 2nd kind
      integer, intent(in) :: nu
      real(kind = KV), intent(in) :: x
      real(kind = KV), intent(out), optional :: oY
      real(kind = KV) :: Y
      Y = bessel_yn(nu, x)
      Y1 = nu/x*Y - bessel_yn(nu + 1, x)
      if(present(oY)) oY = Y
    end function Yv1

#endif

    real(kind = KV) function Jv2(nu, x, oJ, oJ1) result(J2)
      ! The derivative of the 2nd order of Bessel function of the 1st kind
      integer, intent(in) :: nu
      real(kind = KV), intent(in) :: x
      real(kind = KV), intent(out), optional :: oJ, oJ1
      real(kind = KV) :: J, J1
      J1 = Jv1(nu, x, oJ = J)
      J2 = (nu*nu/(x*x) - 1.0d0)*J - J1/x
      if(present(oJ )) oJ  = J
      if(present(OJ1)) oJ1 = J1
    end function Jv2

    real(kind = KV) function Yv2(nu, x, oY, oY1) result(Y2)
      ! The derivative of the 2nd order of Bessel function of the 2nd kind
      integer, intent(in) :: nu
      real(kind = KV), intent(in) :: x
      real(kind = KV), intent(out), optional :: oY, oY1
      real(kind = KV) :: Y, Y1
      Y1 = Yv1(nu, x, oY = Y)
      Y2 = (nu*nu/(x*x) - 1.0d0)*Y - Y1/x
      if(present(oY )) oY  = Y
      if(present(oY1)) oY1 = Y1
    end function Yv2

#ifndef HALLEY

    subroutine jvn(nu, n, zeros)
      ! The zeros of Bessel function of the 1st kind
      integer, intent(in) :: nu, n
      real(kind = KV), intent(out) :: zeros(n)
      integer :: k
      real(kind = KV) :: x, x0, J, J1
      zeros = 0.0d0
      if(nu <= 20) then
        x = 2.82141d0 + 1.15859d0*nu
      else
        x = nu**(1.0d0/3)
        x = nu + 1.85576d0*x + 1.03315d0/x
      end if
      k = 0
      do while(.true.)
        x0 = x
        J1 = Jv1(nu, x, oJ = J)
        ! Fixed-point iteration
        x = x - J/J1
        if(abs(x - x0) > eps) cycle
        k = k + 1
        zeros(k) = x
        if(k >= n) exit
        x = x + 3.1416d0 + (0.0972d0 + 0.0679d0*nu - 0.000354d0*nu*nu)/k
      end do
    end subroutine jvn

    subroutine yvn(nu, n, zeros)
      ! The zeros of Bessel function of the 2nd kind
      integer, intent(in) :: nu, n
      real(kind = KV), intent(out) :: zeros(n)
      integer :: k
      real(kind = KV) :: x, x0, Y, Y1
      zeros = 0.0d0
      if(nu <= 20) then
        x = 1.19477d0 + 1.08933d0*nu
      else
        x = nu**(1.0d0/3)
        x = nu + 0.93158d0*x + 0.26035d0/x
      end if
      k = 0
      do while(.true.)
        x0 = x
        Y1 = Yv1(nu, x, oY = Y)
        ! Fixed-point iteration
        x = x - Y/Y1
        if(abs(x - x0) > eps) cycle
        k = k + 1
        zeros(k) = x
        if(k >= n) exit
        x = x + 3.1416d0 + (0.3120d0 + 0.0852d0*nu - 0.000403d0*nu*nu)/k
      end do
    end subroutine yvn

#else

    subroutine jvn(nu, n, zeros)
      ! The zeros of Bessel function of the 1st kind
      integer, intent(in) :: nu, n
      real(kind = KV), intent(out) :: zeros(n)
      integer :: k
      real(kind = KV) :: x, x0, J, J1, J2
      zeros = 0.0d0
      if(nu <= 20) then
        x = 2.82141d0 + 1.15859d0*nu
      else
        x = nu**(1.0d0/3)
        x = nu + 1.85576d0*x + 1.03315d0/x
      end if
      k = 0
      do while(.true.)
        x0 = x
        J2 = Jv2(nu, x, oJ = J, oJ1 = J1)
        ! Halley's method
        x = x - 2.0d0*J*J1/(2.0d0*J1*J1 - J*J2)
        if(abs(x - x0) > eps) cycle
        k = k + 1
        zeros(k) = x
        if(k >= n) exit
        x = x + 3.1416d0 + (0.0972d0 + 0.0679d0*nu - 0.000354d0*nu*nu)/k
      end do
    end subroutine jvn

    subroutine yvn(nu, n, zeros)
      ! The zeros of Bessel function of the 2nd kind
      integer, intent(in) :: nu, n
      real(kind = KV), intent(out) :: zeros(n)
      integer :: k
      real(kind = KV) :: x, x0, Y, Y1, Y2
      zeros = 0.0d0
      if(nu <= 20) then
        x = 1.19477d0 + 1.08933d0*nu
      else
        x = nu**(1.0d0/3)
        x = nu + 0.93158d0*x + 0.26035d0/x
      end if
      k = 0
      do while(.true.)
        x0 = x
        Y2 = Yv2(nu, x, oY = Y, oY1 = Y1)
        ! Halley's method
        x = x - 2.0d0*Y*Y1/(2.0d0*Y1*Y1 - Y*Y2)
        if(abs(x - x0) > eps) cycle
        k = k + 1
        zeros(k) = x
        if(k >= n) exit
        x = x + 3.1416d0 + (0.3120d0 + 0.0852d0*nu - 0.000403d0*nu*nu)/k
      end do
    end subroutine yvn

#endif

    subroutine jvn1(nu, n, zeros)
      ! The zeros of the derivative of the 1st order of Bessel function
      !   of the 1st kind
      integer, intent(in) :: nu, n
      real(kind = KV), intent(out) :: zeros(n)
      integer :: k
      real(kind = KV) :: x, x0, J1, J2
      zeros = 0.0d0
      if(nu <= 20) then
        x = 0.961587d0 + 1.07703d0*nu
      else
        x = nu**(1.0d0/3)
        x = nu + 0.80861d0*x + 0.07249d0/x
      end if
      if(nu == 0) x = 3.8317d0
      k = 0
      do while(.true.)
        x0 = x
        J2 = Jv2(nu, x, oJ1 = J1)
        ! Fixed-point iteration
        x = x - J1/J2
        if(abs(x - x0) > eps) cycle
        k = k + 1
        zeros(k) = x
        if(k >= n) exit
        x = x + 3.1416d0 + (0.4955d0 + 0.0915d0*nu - 0.000435d0*nu*nu)/k
      end do
    end subroutine jvn1

    subroutine yvn1(nu, n, zeros)
      ! The zeros of the derivative of the 1st order of Bessel function
      !   of the 2nd kind
      integer, intent(in) :: nu, n
      real(kind = KV), intent(out) :: zeros(n)
      integer :: k
      real(kind = KV) :: x, x0, Y1, Y2
      zeros = 0.0d0
      if(nu <= 20) then
        x = 2.67257d0 + 1.16099d0*nu
      else
        x = nu**(1.0d0/3)
        x = nu + 1.8211d0*x + 0.94001d0/x
      end if
      k = 0
      do while(.true.)
        x0 = x
        Y2 = Yv2(nu, x, oY1 = Y1)
        ! Fixed-point iteration
        x = x - Y1/Y2
        if(abs(x - x0) > eps) cycle
        k = k + 1
        zeros(k) = x
        if(k >= n) exit
        x = x + 3.1416d0 + (0.1970d0 + 0.0643d0*nu - 0.000286d0*nu*nu)/k
      end do
    end subroutine yvn1

end module Bessel_Function

!===============================================================================
!*                      An example for using this module                       *
!===============================================================================
!|  program main
!|    use Bessel_Function
!|    implicit none
!|  
!|    integer, parameter :: KV = kind(1.0d0)
!|    integer, parameter :: nu = 1, n = 10
!|    real(kind = KV) :: J0(n), J1(n), Y0(n), Y1(n)
!|    integer :: i
!|    call jvn (nu, n, J0)
!|    call jvn1(nu, n, J1)
!|    call yvn (nu, n, Y0)
!|    call yvn1(nu, n, Y1)
!|    do i = 1, n, 1
!|      write(*, *) i, J0(i), J1(i), Y0(i), Y1(i)
!|    end do
!|  
!|  end program
!===============================================================================
