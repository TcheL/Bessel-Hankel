! Reference:
!   Guptasarma D. and Singh, 1997. New digital linear filters for Hankel J_0
! and J_1 tranforms. Geophysical Prospecting, 45, 745-762.

program main

  use example_Functions
  use Hankel_Transform
  implicit none

  integer :: iExample, nSample
  real(kind = myKind) :: rlogMin, rlogMax, rlogStep

  integer :: fileID, i
  real(kind = myKind) :: r, numRes, anaRes, relErr
  character(len = :), allocatable :: fmtStr

  procedure(real(kind = myKind)), pointer :: pFunc2HankelTransform => null()
  procedure(real(kind = myKind)), pointer :: pFunc2exampleFunction => null()
  namelist /input/ iExample, nSample, rlogMin, rlogMax

  open(newunit = fileID, file = 'input.nml', status = 'old')
    read(fileID, nml = input)
  close(fileID)
  rlogStep = (rlogMax - rlogMin)/(nSample - 1)

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

  if(myKind == kind(1.0e0)) then
    fmtStr = "(1P, E14.7, 1X, E14.7, 1X, E14.7, 1X, 0P, F9.2, A)"
  else
    fmtStr = "(1P, E23.16, 1X, E23.16, 1X, E23.16, 1X, 0P, F9.2, A)"
  end if

  do i = 1, nSample, 1
    r = 10**(rlogMin + rlogStep*(i - 1))
    numRes = pFunc2HankelTransform(pFunc2exampleFunction, r)
    anaRes = pFunc2exampleFunction('r', r)
    relErr = (numRes - anaRes)/anaRes*100.0d0
    write(*, fmtStr) r, numRes, anaRes, relErr, '%'
  end do

end program
