! Reference:
!   Guptasarma D. and Singh, 1997. New digital linear filters for Hankel J_0
! and J_1 tranforms. Geophysical Prospecting, 45, 745-762.

program main

  use example_Functions
  use Hankel_Transform
  implicit none

  integer, parameter :: nSample = 400
  integer :: iExample, i
  real(kind = myKind) :: delta, r
  real(kind = myKind) :: numerSolution(nSample), analySolution(nSample), &
    & relatError(nSample)
  character(len = :), allocatable :: fmtStr

  procedure(real(kind = myKind)), pointer :: pFunc2HankelTransform => null()
  procedure(real(kind = myKind)), pointer :: pFunc2exampleFunction => null()

  delta = 1.0d-03 ! 1.0d+00
  write(*, "(A)", advance = "NO") "> which example you want to check (4-10)? "
  read(*, *) iExample

  select case(iExample)
    case(4)
      pFunc2exampleFunction => exampleFunc1
    case(5)
      pFunc2exampleFunction => exampleFunc2
    case(6)
      pFunc2exampleFunction => exampleFunc3
    case(7)
      pFunc2exampleFunction => exampleFunc4
    case(8)
      pFunc2exampleFunction => exampleFunc5
    case(9)
      pFunc2exampleFunction => exampleFunc6
    case(10)
      pFunc2exampleFunction => exampleFunc7
    case default
      stop "<<< Illegal input option, please have a check."
  end select

  select case(iExample)
#ifndef FAST
    case(4:6)
      pFunc2HankelTransform => HankelJ0Transform
    case(7:10)
      pFunc2HankelTransform => HankelJ1Transform
#else
    case(4:6)
      pFunc2HankelTransform => HankelJ0TransformFast
    case(7:10)
      pFunc2HankelTransform => HankelJ1TransformFast
#endif
  end select

  do i = 1, nSample, 1
    r = i*delta
    numerSolution(i) = pFunc2HankelTransform(pFunc2exampleFunction, r)
    analySolution(i) = pFunc2exampleFunction('r', r)
  end do
        
  relatError = (numerSolution - analySolution)/analySolution*100.0_myKind
! relatError = (analySolution - numerSolution)/numerSolution*100.0_myKind


  if(myKind == kind(1.0e0)) then
    fmtStr = "(I5, 1P, 1X, E14.7, 1X, E14.7, 1X, 0P, F9.2, A)"
  else
    fmtStr = "(I5, 1P, 1X, E23.16, 1X, E23.16, 1X, 0P, F9.2, A)"
  end if

  do i = 1, 400, 1
    write(*, fmtStr) i, numerSolution(i), analySolution(i), relatError(i), '%'
  end do

end program
