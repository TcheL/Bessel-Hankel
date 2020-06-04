#ifdef FAST
#define _HANKJ0TRAN(func, arg) HankelJ0TransformFast(func, arg)
#define _HANKJ1TRAN(func, arg) HankelJ1TransformFast(func, arg)
#else
#define _HANKJ0TRAN(func, arg) HankelJ0Transform(func, arg)
#define _HANKJ1TRAN(func, arg) HankelJ1Transform(func, arg)
#endif

! Reference:
!   Guptasarma D. and Singh, 1997. New digital linear filters for Hankel J_0
! and J_1 tranforms. Geophysical Prospecting, 45, 745-762.

module exampleFunctions

  ! examples from (Guptasarma and Singh, 1997)

  implicit none
  public

  integer, parameter :: myKind = kind(1.0d0)
  real(kind = myKind), private :: c = 1.0_myKind, alpha = 1.0_myKind

  contains

    real(kind = myKind) function exampleFunc1(varName, var) result(func)
      ! equation (4) in (Guptasarma and Singh, 1997)
      character(len = *), intent(in) :: varName
      real(kind = myKind), intent(in) :: var
      if(varName == "lambda") then
        func = exp( - c*var)
      else if(varName == 'r') then
        func = 1.0_myKind/sqrt(c*c + var*var)
      else
        func = 0.0_myKind
      end if
    end function exampleFunc1

    real(kind = myKind) function exampleFunc2(varName, var) result(func)
      ! equation (5) in (Guptasarma and Singh, 1997)
      character(len = *), intent(in) :: varName
      real(kind = myKind), intent(in) :: var
      if(varName == "lambda") then
        func = var*exp( - c*var*var)
      else if(varName == 'r') then
        func = 1.0_myKind/(2*c)*exp( - var*var/(4*c))
      else
        func = 0.0_myKind
      end if
    end function exampleFunc2

    real(kind = myKind) function exampleFunc3(varName, var) result(func)
      ! equation (6) in (Guptasarma and Singh, 1997)
      character(len = *), intent(in) :: varName
      real(kind = myKind), intent(in) :: var
      if(varName == "lambda") then
        func = var*exp( - c*var)
      else if(varName == 'r') then
        func = c/sqrt((c*c + var*var)**3)
      else
        func = 0.0_myKind
      end if
    end function exampleFunc3

    real(kind = myKind) function exampleFunc4(varName, var) result(func)
      ! equation (7) in (Guptasarma and Singh, 1997)
      character(len = *), intent(in) :: varName
      real(kind = myKind), intent(in) :: var
      if(varName == "lambda") then
        func = var*exp( - c*var) &
           & + alpha*var*var*exp( - c*var*var)
      else if(varName == 'r') then
        func = var/sqrt((c*c + var*var)**3) + alpha*var*exp( - var*var/(4*c))/(4*c*c)
      else
        func = 0.0_myKind
      end if
    end function exampleFunc4

    real(kind = myKind) function exampleFunc5(varName, var) result(func)
      ! equation (8) in (Guptasarma and Singh, 1997)
      character(len = *), intent(in) :: varName
      real(kind = myKind), intent(in) :: var
      if(varName == "lambda") then
        func = var*exp( - c*var)
      else if(varName == 'r') then
        func = var/sqrt((c*c + var*var)**3)
      else
        func = 0.0_myKind
      end if
    end function exampleFunc5

    real(kind = myKind) function exampleFunc6(varName, var) result(func)
      ! equation (9) in (Guptasarma and Singh, 1997)
      character(len = *), intent(in) :: varName
      real(kind = myKind), intent(in) :: var
      if(varName == "lambda") then
        func = var*var*exp( - c*var*var)
      else if(varName == 'r') then
        func = var*exp( - var*var/(4*c))/(4*c*c)
      else
        func = 0.0_myKind
      end if
    end function exampleFunc6

    real(kind = myKind) function exampleFunc7(varName, var) result(func)
      ! equation (10) in (Guptasarma and Singh, 1997)
      character(len = *), intent(in) :: varName
      real(kind = myKind), intent(in) :: var
      if(varName == "lambda") then
        func = exp( - c*var)
      else if(varName == 'r') then
        func = (sqrt(var*var + c*c) - c)/(var*sqrt(var*var + c*c))
      else
        func = 0.0_myKind
      end if
    end function exampleFunc7

end module exampleFunctions

module HankelTransform

  ! Hankel transform, i.e. equation (1) in (Guptasarma and Singh, 1997):
  !   $ f(r) = \int_0^\infty K(\lambda) J_i(r \lambda) d\lambda $

  use exampleFunctions, only: myKind
  implicit none
  public

  contains

    real(kind = myKind) function HankelJ0TransformFast(extFunc, r) result(f)
      ! eqations (2) and (3), and table 1 in (Guptasarma and Singh, 1997)
      real(kind = myKind) :: extFunc
      real(kind = myKind), intent(in) :: r
      integer :: i
      integer, parameter :: n = 61
      real(kind = myKind) :: lambda, K(n)
      real(kind = myKind) :: a = -5.08250000000d+00, s = 1.16638303862d-01
      real(kind = myKind) :: W(n) = [ &
&  3.30220475766d-04, -1.18223623458d-03,  2.01879495264d-03, -2.13218719891d-03, &
&  1.60839063172d-03, -9.09156346708d-04,  4.37889252738d-04, -1.55298878782d-04, &
&  7.98411962729d-05,  4.37268394072d-06,  3.94253441247d-05,  4.02675924344d-05, &
&  5.66053344653d-05,  7.25774926389d-05,  9.55412535465d-05,  1.24699163157d-04, &
&  1.63262166579d-04,  2.13477133718d-04,  2.79304232173d-04,  3.65312787897d-04, &
&  4.77899413107d-04,  6.25100170825d-04,  8.17726956451d-04,  1.06961339341d-03, &
&  1.39920928148d-03,  1.83020380399d-03,  2.39417015791d-03,  3.13158560774d-03, &
&  4.09654426763d-03,  5.35807925630d-03,  7.00889482693d-03,  9.16637526490d-03, &
&  1.19891721272d-02,  1.56755740646d-02,  2.04953856060d-02,  2.67778388247d-02, &
&  3.49719672729d-02,  4.55975312615d-02,  5.93498881451d-02,  7.69179091244d-02, &
&  9.91094769804d-02,  1.26166963993d-01,  1.57616825575d-01,  1.89707800260d-01, &
&  2.13804195282d-01,  2.08669340316d-01,  1.40250562745d-01, -3.65385242807d-02, &
& -2.98004010732d-01, -4.21898149249d-01,  5.94373771266d-02,  5.29621428353d-01, &
& -4.41362405166d-01,  1.90355040550d-01, -6.19966386785d-02,  1.87255115744d-02, &
& -5.68736766738d-03,  1.68263510609d-03, -4.38587145792d-04,  8.59117336292d-05, &
& -9.15853765160d-06                                                              &
        & ]
      do i = 1, n, 1
        lambda = 10.0_myKind**(a + (i - 1)*s)/r
        K(i) = extFunc("lambda", lambda)
      end do
      f = dot_product(K, W)/r
    end function HankelJ0TransformFast

    real(kind = myKind) function HankelJ0Transform(extFunc, r) result(f)
      ! eqations (2) and (3), and table 2 in (Guptasarma and Singh, 1997)
      real(kind = myKind) :: extFunc
      real(kind = myKind), intent(in) :: r
      integer :: i
      integer, parameter :: n = 120
      real(kind = myKind) :: lambda, K(n)
      real(kind = myKind) :: a = -8.38850000000d+00, s = 9.0422646867d-02
      real(kind = myKind) :: W(n) = [ &
&  9.62801364263d-07, -5.02069203805d-06,  1.25268783953d-05, -1.99324417376d-05, &
&  2.29149033546d-05, -2.04737583809d-05,  1.49952002937d-05, -9.37502840980d-06, &
&  5.20156955323d-06, -2.62939890538d-06,  1.26550848081d-06, -5.73156151923d-07, &
&  2.76281274155d-07, -1.09963734387d-07,  7.38038330280d-08, -9.31614600001d-09, &
&  3.87247135578d-08,  2.10303178461d-08,  4.10556513877d-08,  4.13077946246d-08, &
&  5.68828741789d-08,  6.59543638130d-08,  8.40811858728d-08,  1.01532550003d-07, &
&  1.26437360082d-07,  1.54733678097d-07,  1.91218582499d-07,  2.35008851918d-07, &
&  2.89750329490d-07,  3.56550504341d-07,  4.39299297826d-07,  5.40794544880d-07, &
&  6.66136379541d-07,  8.20175040653d-07,  1.01015545059d-06,  1.24384500153d-06, &
&  1.53187399787d-06,  1.88633707689d-06,  2.32307100992d-06,  2.86067883258d-06, &
&  3.52293208580d-06,  4.33827546442d-06,  5.34253613351d-06,  6.57906223200d-06, &
&  8.10198829111d-06,  9.97723263578d-06,  1.22867312381d-05,  1.51305855976d-05, &
&  1.86329431672d-05,  2.29456891669d-05,  2.82570465155d-05,  3.47973610445d-05, &
&  4.28521099371d-05,  5.27705217882d-05,  6.49856943660d-05,  8.00269662180d-05, &
&  9.85515408752d-05,  1.21361571831d-04,  1.49454562334d-04,  1.84045784500d-04, &
&  2.26649641428d-04,  2.79106748890d-04,  3.43716968725d-04,  4.23267056591d-04, &
&  5.21251001943d-04,  6.41886194381d-04,  7.90483105615d-04,  9.73420647376d-04, &
&  1.19877439042d-03,  1.47618560844d-03,  1.81794224454d-03,  2.23860214971d-03, &
&  2.75687537633d-03,  3.39471308297d-03,  4.18062141752d-03,  5.14762977308d-03, &
&  6.33918155348d-03,  7.80480111772d-03,  9.61064602702d-03,  1.18304971234d-02, &
&  1.45647517743d-02,  1.79219149417d-02,  2.20527911163d-02,  2.71124775541d-02, &
&  3.33214363101d-02,  4.08864842127d-02,  5.01074356716d-02,  6.12084049407d-02, &
&  7.45146949048d-02,  9.00780900611d-02,  1.07940155413d-01,  1.27267746478d-01, &
&  1.46676027814d-01,  1.62254276550d-01,  1.68045766353d-01,  1.52383204788d-01, &
&  1.01214136498d-01, -2.44389126667d-03, -1.54078468398d-01, -3.03214415655d-01, &
& -2.97674373379d-01,  7.93541259524d-03,  4.26273267393d-01,  1.00032384844d-01, &
& -4.94117404043d-01,  3.92604878741d-01, -1.90111691178d-01,  7.43654896362d-02, &
& -2.78508428343d-02,  1.09992061155d-02, -4.69798719697d-03,  2.12587632706d-03, &
& -9.81986734159d-04,  4.44992546836d-04, -1.89983519162d-04,  7.31024164292d-05, &
& -2.40057837293d-05,  6.23096824846d-06, -1.12363896552d-06,  1.04470606055d-07  &
        & ]
      do i = 1, n, 1
        lambda = 10.0_myKind**(a + (i - 1)*s)/r
        K(i) = extFunc("lambda", lambda)
      end do
      f = dot_product(K, W)/r
    end function HankelJ0Transform

    real(kind = myKind) function HankelJ1TransformFast(extFunc, r) result(f)
      ! eqations (2) and (3), and table 3 in (Guptasarma and Singh, 1997)
      real(kind = myKind) :: extFunc
      real(kind = myKind), intent(in) :: r
      integer :: i
      integer, parameter :: n = 47
      real(kind = myKind) :: lambda, K(n)
      real(kind = myKind) :: a = -3.05078187595d+00, s = 1.10599010095d-01
      real(kind = myKind) :: W(n) = [ &
&  3.17926147465d-06, -9.73811660718d-06,  1.64866227408d-05, -1.81501261160d-05, &
&  1.87556556369d-05, -1.46550406038d-05,  1.53799733803d-05, -6.95628273934d-06, &
&  1.41881555665d-05,  3.41445665537d-06,  2.13941715512d-05,  2.34962369042d-05, &
&  4.84340283290d-05,  7.33732978590d-05,  1.27703784430d-04,  2.08120025730d-04, &
&  3.49803898913d-04,  5.79107814687d-04,  9.65887918451d-04,  1.60401273703d-03, &
&  2.66903777685d-03,  4.43111590040d-03,  7.35631696247d-03,  1.21782796293d-02, &
&  2.01097829218d-02,  3.30096953061d-02,  5.37143591532d-02,  8.60516613299d-02, &
&  1.34267607144d-01,  2.00125033067d-01,  2.74027505792d-01,  3.18168749246d-01, &
&  2.41655667461d-01, -5.40549161658d-02, -4.46912952135d-01, -1.92231885629d-01, &
&  5.52376753950d-01, -3.57429049025d-01,  1.41510519002d-01, -4.61421935309d-02, &
&  1.48273761923d-02, -5.07479209193d-03,  1.83829713749d-03, -6.67742804324d-04, &
&  2.21277518118d-04, -5.66248732755d-05,  7.88229202853d-06                      &
        & ]
      do i = 1, n, 1
        lambda = 10.0_myKind**(a + (i - 1)*s)/r
        K(i) = extFunc("lambda", lambda)
      end do
      f = dot_product(K, W)/r
    end function HankelJ1TransformFast

    real(kind = myKind) function HankelJ1Transform(extFunc, r) result(f)
      ! eqations (2) and (3), and table 4 in (Guptasarma and Singh, 1997)
      real(kind = myKind) :: extFunc
      real(kind = myKind), intent(in) :: r
      integer :: i
      integer, parameter :: n = 140
      real(kind = myKind) :: lambda, K(n)
      real(kind = myKind) :: a = -7.91001919000d+00, s = 8.7967143957d-02
      real(kind = myKind) :: W(n) = [ &
& -6.76671159511d-14,  3.39808396836d-13, -7.43411889153d-13,  8.93613024469d-13, &
& -5.47341591896d-13, -5.84920181906d-14,  5.20780672883d-13, -6.92656254606d-13, &
&  6.88908045074d-13, -6.39910528298d-13,  5.82098912530d-13, -4.84912700478d-13, &
&  3.54684337858d-13, -2.10855291368d-13,  1.00452749275d-13,  5.58449957721d-15, &
& -5.67206735175d-14,  1.09107856853d-13, -6.04067500756d-14,  8.84512134731d-14, &
&  2.22321981827d-14,  8.38072239207d-14,  1.23647835900d-13,  1.44351787234d-13, &
&  2.94276480713d-13,  3.39965995918d-13,  6.17024672340d-13,  8.25310217692d-13, &
&  1.32560792613d-12,  1.90949961267d-12,  2.93458179767d-12,  4.33454210095d-12, &
&  6.55863288798d-12,  9.78324910827d-12,  1.47126365223d-11,  2.20240108708d-11, &
&  3.30577485691d-11,  4.95377381480d-11,  7.43047574433d-11,  1.11400535181d-10, &
&  1.67052734516d-10,  2.50470107577d-10,  3.75597211630d-10,  5.63165204681d-10, &
&  8.44458166896d-10,  1.26621795331d-09,  1.89866561359d-09,  2.84693620927d-09, &
&  4.26886170263d-09,  6.40104325574d-09,  9.59798498616d-09,  1.43918931885d-08, &
&  2.15798696769d-08,  3.23584600810d-08,  4.85195105813d-08,  7.27538583183d-08, &
&  1.09090191748d-07,  1.63577866557d-07,  2.45275193920d-07,  3.67784458730d-07, &
&  5.51470341585d-07,  8.26916206192d-07,  1.23991037294d-06,  1.85921554669d-06, &
&  2.78777669034d-06,  4.18019870272d-06,  6.26794044911d-06,  9.39858833064d-06, &
&  1.40925408889d-05,  2.11312291505d-05,  3.16846342900d-05,  4.75093313246d-05, &
&  7.12354794719d-05,  1.06810848460d-04,  1.60146590551d-04,  2.40110903628d-04, &
&  3.59981158972d-04,  5.39658308918d-04,  8.08925141201d-04,  1.21234066243d-03, &
&  1.81650387595d-03,  2.72068483151d-03,  4.07274689463d-03,  6.09135552241d-03, &
&  9.09940027636d-03,  1.35660714813d-02,  2.01692550906d-02,  2.98534800308d-02, &
&  4.39060697220d-02,  6.39211368217d-02,  9.16763946228d-02,  1.28368795114d-01, &
&  1.73241920046d-01,  2.19830379079d-01,  2.51193131178d-01,  2.32380049895d-01, &
&  1.17121080205d-01, -1.17252913088d-01, -3.52148528535d-01, -2.71162871370d-01, &
&  2.91134747110d-01,  3.17192840623d-01, -4.93075681595d-01,  3.11223091821d-01, &
& -1.36044122543d-01,  5.12141261934d-02, -1.90806300761d-02,  7.57044398633d-03, &
& -3.25432753751d-03,  1.49774676371d-03, -7.24569558272d-04,  3.62792644965d-04, &
& -1.85907973641d-04,  9.67201396593d-05, -5.07744171678d-05,  2.67510121456d-05, &
& -1.40667136728d-05,  7.33363699547d-06, -3.75638767050d-06,  1.86344211280d-06, &
& -8.71623576811d-07,  3.61028200288d-07, -1.05847108097d-07, -1.51569361490d-08, &
&  6.67633241420d-08, -8.33741579804d-08,  8.31065906136d-08, -7.53457009758d-08, &
&  6.48057680299d-08, -5.37558016587d-08,  4.32436265303d-08, -3.37262648712d-08, &
&  2.53558687098d-08, -1.81287021528d-08,  1.20228328586d-08, -7.10898040664d-09, &
&  3.53667004588d-09, -1.36030600198d-09,  3.52544249042d-10, -4.53719284366d-11  &
        & ]
      do i = 1, n, 1
        lambda = 10.0_myKind**(a + (i - 1)*s)/r
        K(i) = extFunc("lambda", lambda)
      end do
      f = dot_product(K, W)/r
    end function HankelJ1Transform

end module HankelTransform

program main

  use exampleFunctions
  use HankelTransform
  implicit none

  integer, parameter :: nSample = 400
  integer :: iExample, i
  real(kind = myKind) :: delta, r
  real(kind = myKind) :: numerSolution(nSample), analySolution(nSample), &
    & relatError(nSample)

  delta = 1.0d-03 ! 1.0d+00
  write(*, "(A)", advance = "NO") "> which example you want to check (4-10)? "
  read(*, *) iExample

  select case(iExample)
    case(4)
      do i = 1, nSample, 1
        r = i*delta
        numerSolution(i) = _HANKJ0TRAN(exampleFunc1, r)
        analySolution(i) = exampleFunc1('r', r)
      end do
    case(5)
      do i = 1, nSample, 1
        r = i*delta
        numerSolution(i) = _HANKJ0TRAN(exampleFunc2, r)
        analySolution(i) = exampleFunc2('r', r)
      end do
    case(6)
      do i = 1, nSample, 1
        r = i*delta
        numerSolution(i) = _HANKJ0TRAN(exampleFunc3, r)
        analySolution(i) = exampleFunc3('r', r)
      end do
    case(7)
      do i = 1, nSample, 1
        r = i*delta
        numerSolution(i) = _HANKJ1TRAN(exampleFunc4, r)
        analySolution(i) = exampleFunc4('r', r)
      end do
    case(8)
      do i = 1, nSample, 1
        r = i*delta
        numerSolution(i) = _HANKJ1TRAN(exampleFunc5, r)
        analySolution(i) = exampleFunc5('r', r)
      end do
    case(9)
      do i = 1, nSample, 1
        r = i*delta
        numerSolution(i) = _HANKJ1TRAN(exampleFunc6, r)
        analySolution(i) = exampleFunc6('r', r)
      end do
    case(10)
      do i = 1, nSample, 1
        r = i*delta
        numerSolution(i) = _HANKJ1TRAN(exampleFunc7, r)
        analySolution(i) = exampleFunc7('r', r)
      end do
    case default
      stop "<<< Illegal input option, please have a check."
  end select
        
  relatError = (numerSolution - analySolution)/analySolution*100_myKind
! relatError = (analySolution - numerSolution)/numerSolution*100_myKind

  do i = 1, 400, 1
    if(myKind == kind(1.0e0)) then
      write(*, "(I5, 1P, 1X, E14.7, 1X, E14.7, 1X, 0P, F9.2, A)") i, &
        & numerSolution(i), analySolution(i), relatError(i), '%'
    else if(myKind == kind(1.0d0)) then
      write(*, "(I5, 1P, 1X, E23.16, 1X, E23.16, 1X, 0P, F9.2, A)") i, &
        & numerSolution(i), analySolution(i), relatError(i), '%'
    end if
  end do

end program
