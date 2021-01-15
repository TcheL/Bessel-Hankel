! Reference:
!   Guptasarma D. and Singh, 1997. New digital linear filters for Hankel J_0
! and J_1 tranforms. Geophysical Prospecting, 45, 745-762.

module example_Functions

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

end module example_Functions

