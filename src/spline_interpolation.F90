module spline_interpolation
	use H5_rdwt
	implicit none
	!This module calculates U kernel using 1) spline interpolation 2) gauss chebyshev quadrature	
contains
	
	! utable contains information of velocity
	! gtable contains information of k1, k2, lambda
	! Note original spline interpolation file is at the 'Spline_Final' Folder
	! File that includes all the process that Pyeon did (including trapezoidal method, adaptive method.. etc) is at spline_interpolation_07_26.F90 file 
	subroutine vcalc(utable, gtable, gamma, gammap)
	! This subroutine gets velocity data from col_f_module and calculateds k1,k2, lambda
	IMPLICIT NONE
	Double precision, parameter		:: cvel = 299792458
	Double precision	:: gamma		!Lorentz factor
	Double precision	:: gammap
	Double precision	:: ur, urp, uz, uzp	!momentum per unit mass in r,z direction
	Double precision	:: lambdacalc		
	Double precision	:: k1calc, k2calc, tcalc
	Double precision, Dimension(4)	:: utable
	Double precision, Dimension(3)	:: gtable

	! From col_f_module.F90 get the velocity information
	
	ur=utable(1)
	urp=utable(2)
	uz=utable(3)
	uzp=utable(4)

	gamma=sqrt(1+((ur**2+uz**2) / cvel**2))
	gammap=sqrt(1+((urp**2+uzp**2) / cvel**2))
	
	! lambda, k1, k2, t calculation
	lambdacalc = (cvel**2) * ( (gamma**2 * gammap**2) + (uz**2 * uzp**2 / cvel**4) - (2*gamma*gammap*uz*uzp / cvel**2) -1 )
	k1calc = (cvel**2/lambdacalc) * ( (2*ur*urp*uz*uzp / cvel**4) - (2*gamma*gammap*ur*urp / cvel**2) )
	k2calc = (cvel**2/lambdacalc) * ( ur**2 * urp**2 / cvel**4 )
	tcalc = abs( 1 + k1calc + k2calc ) ! the case T~~0 is very rare, and for further calculation, T should be larger than 0 so, use abs to make T always > 0  - Jinkyu

	gtable(1) = tcalc
	gtable(2) = k1calc	
	gtable(3) = lambdacalc

	end subroutine vcalc




	! Gauss-Chebyshev quadrature subroutine
	subroutine gauss_chebyshev_quad(ur,urp,uz,uzp,n, NUperpperp,NUperpparl,NUparlparl,NUperpperpprime,NUparlperpprime,T)
	! This subroutine calculates U kernel using gauss_chebyshev quadrature
	  implicit none
	  real (kind=8) :: NUperpperp,NUperpparl,NUparlparl,NUperpperpprime,NUparlperpprime
	  real(kind=8) :: ur, urp, uz, uzp
	  real(kind=8) :: lambda, T, k1, k2
	  real(kind=8) :: gamma, gammap
	  integer, intent(in) :: n ! Number of quadrature points
	  real(kind=8) :: result ! Result of the integration
	  real(kind=8) :: a, b, c, d, x, f, f2
	  integer :: i
	  integer :: option = 0	 	! option1 Default , option2 Singularity Control
	  real(kind=8) :: cvel=299792458D0
	  Double precision, Dimension(3)	:: gtable			! gtable = T=1, k1=2, lambda=3	  
	  Double precision, Dimension(4)	:: utable			! utable = ur=1, urp=2, uz=3, uzp=4
	  Double precision, Dimension(5)	:: stable			! stable = Sx1=1 Sx2=2 Sx3=3 Sx4=4 Sx5=5
	  Real					:: tinit, tfin, tspend		! calculation time 

	!%%%%%%%%%%%%

	  utable(1) = ur
	  utable(2) = urp
	  utable(3) = uz
	  utable(4) = uzp

	  call vcalc(utable, gtable, gamma, gammap)

	  T= gtable(1)
	  k1= gtable(2)
	  k2= T-k1-1
	  lambda= gtable(3)

	  a = -1.0d0 ! Integration interval start
	  b = 1.0d0 ! Integration interval end
	  c = (a + b) / 2.0d0 ! Midpoint of the interval
	  d = (b - a) / 2.0d0 ! Width of the interval

	  if (1+k1+k2<1.5D-2) then  ! The value 1.5D-2 -> not clear, need further consideration
	      option = 2
	  else
	      option = 1		  
	  endif 

	
	  if (option.eq.1) then
	  ! option1: No Singularity Control (Default)
		  result = 0.0d0 ! Initialize the result	  	 
		  do i = 1, n
		    x = c + d * cos((2*i - 1) * acos(-1.0d0) / (2*n)) ! Chebyshev points
		    f = funUperpperpG(x, ur, urp, lambda, k1, k2, cvel) ! Evaluate the integrand at the quadrature point
		    f2= f*acos(-1.0d0)/n
		    result = result + f2 ! Accumulate the result
		  end do

	 	  
		  !write(*,*) 'Uperpperp calculation time, No S C', tfin-tinit

		  NUperpperp = (4.0 * acos(-1.0) * (1.0 / (gamma * gammap)))*(d * result)

		  result= 0.0d0  ! Initialize the result

		  do i = 1, n
		    x = c + d * cos((2*i - 1) * acos(-1.0d0) / (2*n)) ! Chebyshev points
		    f = funUperpparlG(x, ur, urp, uz, uzp, lambda, k1, k2, cvel) ! Evaluate the integrand at the quadrature point
		    f2= f*acos(-1.0d0)/n
		    result = result + f2 ! Accumulate the result
		  end do

		  NUperpparl = (4.0 * acos(-1.0) * (1.0 / (gamma * gammap)))*d * result ! Scale the result by the interval width

		  result= 0.0d0  ! Initialize the result

		  do i = 1, n
		    x = c + d * cos((2*i - 1) * acos(-1.0d0) / (2*n)) ! Chebyshev points
		    f = funUparlparlG(x, uz, uzp, lambda, k1, k2, cvel) ! Evaluate the integrand at the quadrature point
		    f2= f*acos(-1.0d0)/n
		    result = result + f2 ! Accumulate the result
		  end do
		 
		  NUparlparl = (4.0 * acos(-1.0) * (1.0 / (gamma * gammap)))*d * result ! Scale the result by the interval width

		  result= 0.0d0  ! Initialize the result

		  do i = 1, n
		    x = c + d * cos((2*i - 1) * acos(-1.0d0) / (2*n)) ! Chebyshev points
		    f = funUperpperpprimeG(x, ur, urp, lambda, k1, k2, cvel) ! Evaluate the integrand at the quadrature point
		    f2= f*acos(-1.0d0)/n
		    result = result + f2 ! Accumulate the result
		  end do
		  
		  NUperpperpprime = (4.0 * acos(-1.0) * (1.0 / (gamma * gammap)))*d * result ! Scale the result by the interval width

		  result= 0.0d0  ! Initialize the result

		  do i = 1, n
		    x = c + d * cos((2*i - 1) * acos(-1.0d0) / (2*n)) ! Chebyshev points
		    f = funUparlperpprimeG(x, ur, urp, uz, uzp, lambda, k1, k2, cvel) ! Evaluate the integrand at the quadrature point
		    f2= f*acos(-1.0d0)/n
		    result = result + f2 ! Accumulate the result
		  end do

		  NUparlperpprime = (4.0 * acos(-1.0) * (1.0 / (gamma * gammap)))*d * result ! Scale the result by the interval width

	  elseif (option.eq.2) then
	  ! option2 : Singularity Control
		  
		  result = 0.0d0 ! Initialize the result	  	 
		  do i = 1, n
		    x = c + d * cos((2*i - 1) * acos(-1.0d0) / (2*n)) ! Chebyshev points
		    f = funUperpperpS(x, ur, urp, lambda, k1, k2, cvel) ! Evaluate the integrand at the quadrature point
		    f2= f*acos(-1.0d0)/n
		    result = result + f2 ! Accumulate the result
		  end do


		  ! Scale the result by the interval width
		  NUperpperp = (4.0 * acos(-1.0) * (1.0 / (gamma * gammap)))* &
		  (d * result + & ! adding anayltic integration value of singularity 	 
		  (1-(ur**2.0d0)/(cvel**2.0d0)-(urp**2.0d0)/(cvel**2.0d0))* &
		  (lambda**(-1.0d0/2.0d0))*((1.0d0) * sqrt(2.0d0) * asinh(sqrt((2-2*k2)/t))/ sqrt(1-k2)) + &
		  (-ur**2.0d0-urp**2.0d0)* &
		  (lambda**(-3.0d0/2.0d0))*((1.0d0) * sqrt(2.0d0) / ((t**(3.0d0/2.0d0)) * sqrt(1.0d0 + t / (2.0d0 - k2)) * sqrt((1.0d0 - k2) / t))) + &
		  (2.0d0*ur*urp/(cvel**2))* &
		  (lambda**(-1.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)) * sqrt(2.0d0) * asinh(sqrt((2-2*k2)/t))/ sqrt(1-k2))+  &
		  (2.0d0*ur*urp)* &
		  (lambda**(-3.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)) * sqrt(2.0d0) / ((t**(3.0d0/2.0d0)) * sqrt(1.0d0 + t / (2.0d0 - k2)) * sqrt((1.0d0 - k2) / t))) &
		  )

		  result= 0.0d0  ! Initialize the result

		  do i = 1, n
		    x = c + d * cos((2*i - 1) * acos(-1.0d0) / (2*n)) ! Chebyshev points
		    f = funUperpparlS(x, ur, urp, uz, uzp, lambda, k1, k2, cvel) ! Evaluate the integrand at the quadrature point
		    f2= f*acos(-1.0d0)/n
		    result = result + f2 ! Accumulate the result
		  end do

		  NUperpparl = (4.0 * acos(-1.0) * (1.0 / (gamma * gammap)))* &
		  (d * result + & ! adding anayltic integration value of singularity	 
		  (-(ur*uz)/(cvel**2.0d0)-(urp*uzp)/(cvel**2.0d0))* &
		  (lambda**(-1.0d0/2.0d0))*((1.0d0) * sqrt(2.0d0) * asinh(sqrt((2-2*k2)/t))/ sqrt(1-k2)) + &
		  (-ur*uz-urp*uzp)* &
		  (lambda**(-3.0d0/2.0d0))*((1.0d0) * sqrt(2.0d0) / ((t**(3.0d0/2.0d0)) * sqrt(1.0d0 + t / (2.0d0 - k2)) * sqrt((1.0d0 - k2) / t))) + &
		  ((uz*urp)/(cvel**2.0d0)+(ur*uzp)/(cvel**2.0d0))* &
		  (lambda**(-1.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)) * sqrt(2.0d0) * asinh(sqrt((2-2*k2)/t))/ sqrt(1-k2))+  &
		  (uz*urp+ur*uzp)* &
		  (lambda**(-3.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)) * sqrt(2.0d0) / ((t**(3.0d0/2.0d0)) * sqrt(1.0d0 + t / (2.0d0 - k2)) * sqrt((1.0d0 - k2) / t))) &
		  )

		  result= 0.0d0  ! Initialize the result

		  do i = 1, n
		    x = c + d * cos((2*i - 1) * acos(-1.0d0) / (2*n)) ! Chebyshev points
		    f = funUparlparlS(x, uz, uzp, lambda, k1, k2, cvel) ! Evaluate the integrand at the quadrature point
		    f2= f*acos(-1.0d0)/n
		    result = result + f2 ! Accumulate the result
		  end do
		  NUparlparl = (4.0 * acos(-1.0) * (1.0 / (gamma * gammap)))* &
		  (d * result + & ! adding anayltic integration value of singularity	 
		  (1.0d0-(uz**2.0d0)/(cvel**2.0d0)-(uzp**2.0d0)/(cvel**2.0d0))* &
		  (lambda**(-1.0d0/2.0d0))*((1.0d0) * sqrt(2.0d0) * asinh(sqrt((2-2*k2)/t))/ sqrt(1-k2)) + &
		  (-uz**2.0d0-uzp**2.0d0)* &
		  (lambda**(-3.0d0/2.0d0))*((1.0d0) * sqrt(2.0d0) / ((t**(3.0d0/2.0d0)) * sqrt(1.0d0 + t / (2.0d0 - k2)) * sqrt((1.0d0 - k2) / t))) + &
		  (2.0d0*uz*uzp/(cvel**2))* &
		  (lambda**(-1.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)) * sqrt(2.0d0) * asinh(sqrt((2-2*k2)/t))/ sqrt(1-k2))+  &
		  (2.0d0*uz*uzp)* &
		  (lambda**(-3.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)) * sqrt(2.0d0) / ((t**(3.0d0/2.0d0)) * sqrt(1.0d0 + t / (2.0d0 - k2)) * sqrt((1.0d0 - k2) / t))) &
		  )

		  result= 0.0d0  ! Initialize the result

		  do i = 1, n
		    x = c + d * cos((2*i - 1) * acos(-1.0d0) / (2*n)) ! Chebyshev points
		    f = funUperpperpprimeS(x, ur, urp, lambda, k1, k2, cvel) ! Evaluate the integrand at the quadrature point
		    f2= f*acos(-1.0d0)/n
		    result = result + f2 ! Accumulate the result
		  end do

		  NUperpperpprime = (4.0 * acos(-1.0) * (1.0 / (gamma * gammap)))* &
		  (d * result + & ! adding anayltic integration value of singularity	 
		  (1.0d0-(ur**2.0d0)/(cvel**2.0d0)-(urp**2.0d0)/(cvel**2.0d0))* &
		  (lambda**(-1.0d0/2.0d0))*((1.0d0) * sqrt(2.0d0) * asinh(sqrt((2-2*k2)/t))/ sqrt(1-k2)) + &
		  (-ur**2.0d0-urp**2.0d0)* &
		  (lambda**(-3.0d0/2.0d0))*((1.0d0) * sqrt(2.0d0) / ((t**(3.0d0/2.0d0)) * sqrt(1.0d0 + t / (2.0d0 - k2)) * sqrt((1.0d0 - k2) / t))) + &
		  ((2.0d0*ur*urp)/(cvel**2.0d0))* &
		  (lambda**(-1.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)) * sqrt(2.0d0) * asinh(sqrt((2-2*k2)/t))/ sqrt(1-k2))+  &
		  (2.0d0*ur*urp)* &
		  (lambda**(-3.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)) * sqrt(2.0d0) / ((t**(3.0d0/2.0d0)) * sqrt(1.0d0 + t / (2.0d0 - k2)) * sqrt((1.0d0 - k2) / t))) &
		  )
		  
		  result= 0.0d0  ! Initialize the result

		  do i = 1, n
		    x = c + d * cos((2*i - 1) * acos(-1.0d0) / (2*n)) ! Chebyshev points
		    f = funUparlperpprimeS(x, ur, urp, uz, uzp, lambda, k1, k2, cvel) ! Evaluate the integrand at the quadrature point
		    f2= f*acos(-1.0d0)/n
		    result = result + f2 ! Accumulate the result
		  end do

		  NUparlperpprime = (4.0 * acos(-1.0) * (1.0 / (gamma * gammap)))* &
		  (d * result + & ! adding anayltic integration value of singularity	 
		  ((-urp*uzp)/(cvel**2.0d0)+(-ur*uz)/(cvel**2.0d0))* &
		  (lambda**(-1.0d0/2.0d0))*((1.0d0) * sqrt(2.0d0) * asinh(sqrt((2-2*k2)/t))/ sqrt(1-k2)) + &
		  (-urp*uzp-ur*uz)* &
		  (lambda**(-3.0d0/2.0d0))*((1.0d0) * sqrt(2.0d0) / ((t**(3.0d0/2.0d0)) * sqrt(1.0d0 + t / (2.0d0 - k2)) * sqrt((1.0d0 - k2) / t))) + &
		  ((uz*urp)/(cvel**2.0d0)+(ur*uzp)/(cvel**2.0d0))* &
		  (lambda**(-1.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)) * sqrt(2.0d0) * asinh(sqrt((2-2*k2)/t))/ sqrt(1-k2))+  &
		  (uz*urp+ur*uzp)* &
		  (lambda**(-3.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)) * sqrt(2.0d0) / ((t**(3.0d0/2.0d0)) * sqrt(1.0d0 + t / (2.0d0 - k2)) * sqrt((1.0d0 - k2) / t))) &
		  )

	  endif 



	end subroutine gauss_chebyshev_quad



!! 					 Function Define					!!
!================================================================================================!
!!												!!
!!												!!
!!					       Uperpperp					!!
!!												!!
!!												!!
!================================================================================================!


!==========================Equation for Singularity Cancellation Gauss Chebyshev=================!
!Below function is written as following order
!(coefficients)*(lambda value)*(r)*(singularity term)   => Here, r^2= 1+ w^2/c^2  (Definition)
!funUperpperpS -> S denotes for Singularity control
	function funUperpperpS(x, ur, urp, lambda, k1, k2, cvel)
	  implicit none
	  real(kind=8), intent(in) :: x, ur, urp, lambda, k1, k2, cvel
	  real(kind=8) :: funUperpperpS

	  funUperpperpS = ((lambda**(-3.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-3.0d0/2.0d0)) + &
	    ((lambda**(-1.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-1.0d0/2.0d0)) / cvel**2)) * &
	      (lambda * (1.0d0+k1*x+k2*x**2) - ur**2 - urp**2*x**2 + &
	      2.0d0*ur*urp*x*(sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))) - &  !! Singularity Cancellation
	      (1.0d0-(ur**2.0d0)/(cvel**2.0d0)-(urp**2.0d0)/(cvel**2.0d0))* &
	      (lambda**(-1.0d0/2.0d0))*(1)*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**1)))*(sqrt(1-x**2))- &
	      (-ur**2.0d0-urp**2.0d0)* &
	      (lambda**(-3.0d0/2.0d0))*(1)*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**3)))*(sqrt(1-x**2))- &
	      (2.0d0*ur*urp/(cvel**2))* & 
    	      (lambda**(-1.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)))*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**1)))*(sqrt(1-x**2))- &
	      (2.0d0*ur*urp)* &
	      (lambda**(-3.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)))*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**3)))*(sqrt(1-x**2))      
	end function funUperpperpS

!==================================Equation for Gauss Chebyshev==================================!
! G denotes for Gauss Chebyshev	
	function funUperpperpG(x, ur, urp, lambda, k1, k2, cvel)
	  implicit none
	  real(kind=8), intent(in) :: x, ur, urp, lambda, k1, k2, cvel
	  real(kind=8) :: funUperpperpG

	  funUperpperpG = ((lambda**(-3.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-3.0d0/2.0d0)) + &
	    ((lambda**(-1.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-1.0d0/2.0d0)) / cvel**2)) * &
	      (lambda * (1.0d0+k1*x+k2*x**2) - ur**2 - urp**2*x**2 + &
	      2.0d0*ur*urp*x*(sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))) 
	end function funUperpperpG

!==================================Equation for other integation=================================!
! If you need function for other integration method trapezoidal, simpson ... etc, use this function. (Note that sqrt(1-x^2) is multiplied)
	function funUperpperp(x, ur, urp, lambda, k1, k2, cvel)
	  implicit none
	  real(kind=8), intent(in) :: x, ur, urp, lambda, k1, k2, cvel
	  real(kind=8) :: funUperpperp
	 
 	funUperpperp = (1.0d0/sqrt(1.0d0-x**2)) * &
	     ((lambda**(-3.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-3.0d0/2.0d0)) + &
	      ((lambda**(-1.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-1.0d0/2.0d0)) / cvel**2)) * &
	      (lambda * (1.0d0+k1*x+k2*x**2) - ur**2 - urp**2*x**2 + &
	      2.0d0*ur*urp*x*(sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2)))
	end function funUperpperp
!================================================================================================!
!!												!!
!!												!!
!!					       Uperpparl					!!
!!												!!
!!												!!
!================================================================================================!

!==========================Equation for Singularity Cancellation Gauss Chebyshev=================!
	function funUperpparlS(x, ur, urp, uz, uzp, lambda, k1, k2, cvel)
	  implicit none

	  real(kind=8), intent(in) :: x, ur, urp, uz, uzp, lambda, k1, k2, cvel
	  real(kind=8) :: funUperpparlS

	funUperpparlS = ((lambda**(-1.5d0) * (1.0d0+k1*x+k2*x**2)**(-1.5d0) + &
               lambda**(-0.5d0) * (1.0d0+k1*x+k2*x**2)**(-0.5d0) / cvel**2) * &
              (-ur*uz - x*urp*uzp + &
              (sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2)) * &
              (x*uz*urp + ur*uzp)))- &  !! Singularity Cancellation
	      (-(ur*uz)/(cvel**2.0d0)-(urp*uzp)/(cvel**2.0d0))* &
	      (lambda**(-1.0d0/2.0d0))*(1)*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**1)))*(sqrt(1-x**2)) - & 
	      (-ur*uz-urp*uzp)* &
	      (lambda**(-3.0d0/2.0d0))*(1)*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**3)))*(sqrt(1-x**2)) - &
	      ((uz*urp)/(cvel**2.0d0)+(ur*uzp)/(cvel**2.0d0))* & 
   	      (lambda**(-1.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)))*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**1)))*(sqrt(1-x**2))- &
	      (uz*urp+ur*uzp)* &
	      (lambda**(-3.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)))*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**3)))*(sqrt(1-x**2)) 
	end function funUperpparlS
!==================================Equation for Gauss Chebyshev==================================!	
	function funUperpparlG(x, ur, urp, uz, uzp, lambda, k1, k2, cvel)
	  implicit none

	  real(kind=8), intent(in) :: x, ur, urp, uz, uzp, lambda, k1, k2, cvel
	  real(kind=8) :: funUperpparlG
	funUperpparlG = ((lambda**(-1.5d0) * (1.0d0+k1*x+k2*x**2)**(-1.5d0) + &
              lambda**(-0.5d0) * (1.0d0+k1*x+k2*x**2)**(-0.5d0) / cvel**2) * &
              (-ur*uz - x*urp*uzp + &
              (sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2)) * &
              (x*uz*urp + ur*uzp)))
	end function funUperpparlG
!==================================Equation for other integation=================================!
	function funUperpparl(x, ur, urp, uz, uzp, lambda, k1, k2, cvel)
	  implicit none

	  real(kind=8), intent(in) :: x, ur, urp, uz, uzp, lambda, k1, k2, cvel
	  real(kind=8) :: funUperpparl
	funUperpparl = (1.0d0/sqrt(1.0d0-x**2)) * &
              ((lambda**(-1.5d0) * (1.0d0+k1*x+k2*x**2)**(-1.5d0) + &
               lambda**(-0.5d0) * (1.0d0+k1*x+k2*x**2)**(-0.5d0) / cvel**2) * &
              (-ur*uz - x*urp*uzp + &
              (sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2)) * &
              (x*uz*urp + ur*uzp)))

	end function funUperpparl

!================================================================================================!
!!												!!
!!												!!
!!					       Uparlparl					!!
!!												!!
!!												!!
!================================================================================================!


!==========================Equation for Singularity Cancellation Gauss Chebyshev=================!
	function funUparlparlS(x, uz, uzp, lambda, k1, k2, cvel)
	  implicit none
	  real(kind=8), intent(in) :: x, uz, uzp, lambda, k1, k2, cvel
	  real(kind=8) :: funUparlparlS
	  funUparlparlS = (((lambda**(-3.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2.0d0)**(-3.0d0/2.0d0))) + &
	    (((lambda**(-1.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2.0d0)**(-1.0d0/2.0d0)))/cvel**2)) * &
	    (lambda*(1.0d0+k1*x+k2*x**2) - uz**2 - uzp**2 + &
	    2.0d0*(sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*(uz*uzp)) - &  !! Singularity Cancellation
	      (1.0d0-(uz**2.0d0)/(cvel**2.0d0)-(uzp**2.0d0)/(cvel**2.0d0))* &
	      (lambda**(-1.0d0/2.0d0))*(1)*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**1)))*(sqrt(1-x**2))- &
	      (-uz**2.0d0-uzp**2.0d0)* &
	      (lambda**(-3.0d0/2.0d0))*(1)*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**3)))*(sqrt(1-x**2))- &
	      (2.0d0*uz*uzp/(cvel**2))* & 
    	      (lambda**(-1.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)))*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**1)))*(sqrt(1-x**2))- &
	      (2.0d0*uz*uzp)* &
	      (lambda**(-3.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)))*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**3)))*(sqrt(1-x**2)) 
	end function funUparlparlS     
!==================================Equation for Gauss Chebyshev==================================!
	function funUparlparlG(x, uz, uzp, lambda, k1, k2, cvel)
	  implicit none
	  real(kind=8), intent(in) :: x, uz, uzp, lambda, k1, k2, cvel
	  real(kind=8) :: funUparlparlG
	 funUparlparlG = (((lambda**(-3.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2.0d0)**(-3.0d0/2.0d0))) + &
	    (((lambda**(-1.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2.0d0)**(-1.0d0/2.0d0)))/cvel**2)) * &
	    (lambda*(1.0d0+k1*x+k2*x**2) - uz**2 - uzp**2 + &
	    2.0d0*(sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*(uz*uzp))
	end function funUparlparlG
!==================================Equation for other integation=================================!
	function funUparlparl(x, uz, uzp, lambda, k1, k2, cvel)
	  implicit none
	  real(kind=8), intent(in) :: x, uz, uzp, lambda, k1, k2, cvel
	  real(kind=8) :: funUparlparl
	  funUparlparl = (1.0d0/sqrt(1d0-x**2)) * &
	    (((lambda**(-3.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2.0d0)**(-3.0d0/2.0d0))) + &
	    (((lambda**(-1.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2.0d0)**(-1.0d0/2.0d0)))/cvel**2)) * &
	    (lambda*(1.0d0+k1*x+k2*x**2) - uz**2 - uzp**2 + &
	    2.0d0*(sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*(uz*uzp))
	end function funUparlparl


!================================================================================================!
!!												!!
!!												!!
!!					  Uperpperpprime					!!
!!												!!
!!												!!
!================================================================================================!

	
!==========================Equation for Singularity Cancellation Gauss Chebyshev=================!	
	function funUperpperpprimeS(x, ur, urp, lambda, k1, k2, cvel)
	  implicit none
	  real(kind=8), intent(in) :: x, ur, urp, lambda, k1, k2, cvel
	  real(kind=8) :: funUperpperpprimeS	
	funUperpperpprimeS =(((lambda**(-3.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-3.0d0/2.0d0))) + &
	    (((lambda**(-1.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-1.0d0/2.0d0)))/cvel**2)) * &
	    ((x)*((lambda)*(1.0d0+k1*x+k2*x**2)- &
	    ur**2+2.0d0*(sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*(x)*(ur*urp)-(x**2)*(urp**2)) + &
	    (1.0d0-x**2)*((sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*(ur*urp)-(x)*(urp**2)))- &  !! Singularity Cancellation
	      (1.0d0-(ur**2.0d0)/(cvel**2.0d0)-(urp**2.0d0)/(cvel**2.0d0))* &
	      (lambda**(-1.0d0/2.0d0))*(1)*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**1)))*(sqrt(1-x**2))- &
	      (-ur**2.0d0-urp**2.0d0)* &
	      (lambda**(-3.0d0/2.0d0))*(1)*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**3)))*(sqrt(1-x**2))- &
	      ((2.0d0*ur*urp)/(cvel**2.0d0))* & 
    	      (lambda**(-1.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)))*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**1)))*(sqrt(1-x**2))- &
	      (2.0d0*ur*urp)* &
	      (lambda**(-3.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)))*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**3)))*(sqrt(1-x**2))           
	end function funUperpperpprimeS

!==================================Equation for Gauss Chebyshev==================================!
	function funUperpperpprimeG(x, ur, urp, lambda, k1, k2, cvel)
	  implicit none
	  real(kind=8), intent(in) :: x, ur, urp, lambda, k1, k2, cvel
	  real(kind=8) :: funUperpperpprimeG	
	funUperpperpprimeG =(((lambda**(-3.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-3.0d0/2.0d0))) + &
	    (((lambda**(-1.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-1.0d0/2.0d0)))/cvel**2)) * &
	    ((x)*((lambda)*(1.0d0+k1*x+k2*x**2)- &
	    ur**2+2.0d0*(sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*(x)*(ur*urp)-(x**2)*(urp**2)) + &
	    (1.0d0-x**2)*((sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*(ur*urp)-(x)*(urp**2)))
	end function funUperpperpprimeG
!==================================Equation for other integation=================================!
	function funUperpperpprime(x, ur, urp, lambda, k1, k2, cvel)
	  implicit none
	  real(kind=8), intent(in) :: x, ur, urp, lambda, k1, k2, cvel
	  real(kind=8) :: funUperpperpprime
	funUperpperpprime = (1.0d0/sqrt(1.0d0-x**2)) * &
	    (((lambda**(-3.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-3.0d0/2.0d0))) + &
	    (((lambda**(-1.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-1.0d0/2.0d0)))/cvel**2)) * &
	    ((x)*((lambda)*(1.0d0+k1*x+k2*x**2)- &
	    ur**2+2.0d0*(sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*(x)*(ur*urp)-(x**2)*(urp**2)) + &
	    (1.0d0-x**2)*((sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*(ur*urp)-(x)*(urp**2)))

	end function funUperpperpprime

!================================================================================================!
!!												!!
!!												!!
!!					  Uparlperpprime					!!
!!												!!
!!												!!
!================================================================================================!

!==========================Equation for Singularity Cancellation Gauss Chebyshev=================!
	function funUparlperpprimeS(x, ur, urp, uz, uzp, lambda, k1, k2, cvel)
	  implicit none
	  real(kind=8), intent(in) :: x, ur, urp, uz, uzp, lambda, k1, k2, cvel
	  real(kind=8) :: funUparlperpprimeS

	  funUparlperpprimeS = (((lambda**(-3.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-3.0d0/2.0d0))) + &
	    (((lambda**(-1.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-1.0d0/2.0d0)))/cvel**2)) * &
	    ((sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*(1.0d0-x**2)*(uz*urp) - &
	    (1.0d0-x**2)*(urp*uzp) + x*(-ur*uz-x*urp*uzp + &
	    (sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*((x)*(uz*urp) + (ur*uzp))))- &  !! Singularity Cancellation
	      ((-urp*uzp)/(cvel**2.0d0)+(-ur*uz)/(cvel**2.0d0))* &
	      (lambda**(-1.0d0/2.0d0))*(1)*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**1)))*(sqrt(1-x**2))- &
	      (-urp*uzp-ur*uz)* &
	      (lambda**(-3.0d0/2.0d0))*(1)*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**3)))*(sqrt(1-x**2))- &
	      ((uz*urp)/(cvel**2.0d0)+(ur*uzp)/(cvel**2.0d0))* & 
    	      (lambda**(-1.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)))*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**1)))*(sqrt(1-x**2))- &
	      (uz*urp+ur*uzp)* &
	      (lambda**(-3.0d0/2.0d0))*(sqrt(1+lambda*(1+k1+k2)/(cvel**2)))*(1.0d0 / sqrt((2.d0) * (-(x - 1.d0) * ((k2 - 1.d0) * (x - 1.d0) + (1+k1+k2))**3)))*(sqrt(1-x**2))    
	end function funUparlperpprimeS
!==================================Equation for Gauss Chebyshev==================================!
	function funUparlperpprimeG(x, ur, urp, uz, uzp, lambda, k1, k2, cvel)
	  implicit none
	  real(kind=8), intent(in) :: x, ur, urp, uz, uzp, lambda, k1, k2, cvel
	  real(kind=8) :: funUparlperpprimeG

	  funUparlperpprimeG = (((lambda**(-3.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-3.0d0/2.0d0))) + &
	    (((lambda**(-1.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-1.0d0/2.0d0)))/cvel**2)) * &
	    ((sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*(1.0d0-x**2)*(uz*urp) - &
	    (1.0d0-x**2)*(urp*uzp) + x*(-ur*uz-x*urp*uzp + &
	    (sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*((x)*(uz*urp) + (ur*uzp))))
	end function funUparlperpprimeG
!==================================Equation for other integation=================================!
	function funUparlperpprime(x, ur, urp, uz, uzp, lambda, k1, k2, cvel)
	  implicit none
	  real(kind=8), intent(in) :: x, ur, urp, uz, uzp, lambda, k1, k2, cvel
	  real(kind=8) :: funUparlperpprime

	  funUparlperpprime = (1.0d0/sqrt(1.0d0-x**2)) * &
	    (((lambda**(-3.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-3.0d0/2.0d0))) + &
	    (((lambda**(-1.0d0/2.0d0)) * ((1.0d0+k1*x+k2*x**2)**(-1.0d0/2.0d0)))/cvel**2)) * &
	    ((sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*(1.0d0-x**2)*(uz*urp) - &
	    (1.0d0-x**2)*(urp*uzp) + x*(-ur*uz-x*urp*uzp + &
	    (sqrt(1.0d0+lambda*(1.0d0+k1*x+k2*x**2)/cvel**2))*((x)*(uz*urp) + (ur*uzp))))

	end function funUparlperpprime



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPLINE BEGIN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPLINE BEGIN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------------------------------------------------------------------------------------------------------------------------------
	subroutine spline2D(Tarray, k1array , rdata, adata, bdata, cdata, ddata, target, Sinterp)
	implicit none
!	Character(len=*)	:: filenamerd						! filename to read data
!	Character(len=*)	:: dsetrdname						! dataset name to read at file
!	Character(len=*)	:: dsetwrname						! dataset name to write at file	
!	Character(len=*)	:: cname						! name of  C coefficient for T direction
!	Character(len=*)	:: dname						! name of  D coefficient for T direction
!	Character(len=*)	:: bname						! name of  B coefficient for T direction
	integer		,  parameter 	:: k1 = 101, T = 501			  	! number of points
	Double precision,  parameter	:: hk = 0.02D0, hT = 0.002D0			! grid interval
	Double precision 		:: Tarray(1,T), k1array(1,k1)			! grid array for interpolarion
	integer 			:: l						! iteration number
	Double precision		:: rdata(T,k1), cdata(T,k1)			! dataset and basic spline coefficent buffer
	Double precision		:: ddata(T-1,k1), bdata(T-1,k1), adata(T-1,k1)	! spline coefficient buffer
	integer				:: idx1, idx2					! grid index
	Double precision		:: xp1, xp2					! x-x(index) = displacement from on-grid point	
	Double precision		:: yt1(4)					! interpolated value at T dir
	Double precision		:: yt2						! interpolated value at T, k1 dir
	Double precision		:: target(3)					! interpolation target grid for verifying interpolation 
	Double precision, intent(out)	:: Sinterp					! interpolation target result for verifying interpolation
	Real				:: tinit, tfin, tspend				! calculation time 

	!%%%%%%%%%%%%
		call cpu_time(tinit)
	!%%%%%%%%%%%%

 	! READ DATA
!	call h5_read2(filenamerd, dsetrdname, rdata)      			 	! 
!	call h5_cread2(cdata, cname, bdata, bname, ddata, dname)			! spline coefficients in T direction
!	adata= rdata(1:T-1,:)								! calculate spline coefficient a in T direction
!	call h5_gread2(filenamerd, "Tarray", "k1array", Tarray, k1array)		! grid data 

	!%%%%%%%%%%%%
!		call cpu_time(tfin)
!		tspend=tfin-tinit
		!write(*,*) "2D readdata time="
		!write(*,*) tspend
	!%%%%%%%%%%%%


	!%%%%%%%%%%%%
		call cpu_time(tinit)
	!%%%%%%%%%%%%


	!do l = 1, 140
	! find 2D index
	call find_idx(hT, T, target(1), Tarray, idx1, xp1)
	call find_idx(hk, k1, target(2), k1array, idx2, xp2)
	!write(*,*) 'Tindex',idx1
	!write(*,*) 'k1index',idx2

	! interpolate in T dim using table
	
		! assing interpolation target matrix (4*4) and interpolate in 1D : result = yt1[4]
	if (idx1 .NE. T) then
		if (idx2 <= k1-3) then
			yt1 = ddata(idx1,idx2:idx2+3)*xp1**3 &
			    + cdata(idx1,idx2:idx2+3)*xp1**2 &
			    + bdata(idx1,idx2:idx2+3)*xp1 &
			    + adata(idx1,idx2:idx2+3)
		else
			yt1 = ddata(idx1,k1-3:k1)*xp1**3 &
			    + cdata(idx1,k1-3:k1)*xp1**2 &
			    + bdata(idx1,k1-3:k1)*xp1 &
			    + adata(idx1,k1-3:k1)
		end if
	else
		if (idx2 <= k1-3) then
			yt1 = rdata(T,idx2:idx2+3)
		else
			yt1 = rdata(T,k1-3:k1)
		end if
	end if


	! interpolate in k1 dim, result : yt2
		call interpolate3(k1, idx2, xp2, yt1, hk, yt2)	

	! assign interpolation result to result array 
		Sinterp = yt2
	!end do	


	!%%%%%%%%%%%%
		call cpu_time(tfin)
		tspend=tfin-tinit
		!write(*,*) "spline2D time="
		!write(*,*) tspend
	!%%%%%%%%%%%%

	end subroutine spline2D

!-----------------------------------------------------------------------------------------------------------------------------------------------------
	!subroutine spline3D(filenamerd, dsetrdname, dsetwrname, cname, dname, bname, target, Sinterp)
	subroutine spline3D(Tarray, k1array, loglambdaarray, rdata, adata, bdata, cdata, ddata, target, Sinterp)
	implicit none
!	Character(len=*)		:: filenamerd									! filename to read data
!	Character(len=*)		:: dsetrdname									! dataset name to read at file
!	Character(len=*)		:: dsetwrname									! dataset name to write at file	
!	Character(len=*)		:: cname									! name of  C coefficient for T direction
!	Character(len=*)		:: dname									! name of  D coefficient for T direction
!	Character(len=*)		:: bname									! name of  B coefficient for T direction
	integer		,  parameter 	:: k1 = 101, T = 501, loglambda = 136   				! number of point
	Double precision,  parameter	:: hk = 0.02D0, hT = 0.002D0, hl = 0.2D0				! grid interval
	Double precision 		:: Tarray(1,T), k1array(1,k1), loglambdaarray(1,loglambda)		! point array for interpolarion
	integer 			:: l									! iteration number
	Double precision		:: rdata(T,k1,loglambda), cdata(T,k1,loglambda)				! dataset buffer
	Double precision		:: ddata(T-1,k1,loglambda), bdata(T-1,k1,loglambda), adata(T-1,k1,loglambda)	! dataset buffer
	integer				:: idx1, idx2, idx3							! grid index
	Double precision		:: xp1, xp2, xp3							! x-x(index)	
	Double precision		:: yt1(4,4)								! interpolated value at T dir
	Double precision		:: yt2(4)								! interpolated value at T, k1 dir
	Double precision		:: yt3									! final interpolated value
	Double precision		:: target(3)								! interpolation target grid
	Double precision, intent(out)	:: Sinterp								! interpolation target result
	Real				:: tinit, tfin, tspend							! calc cpu time 



	!%%%%%%%%%%%%
		call cpu_time(tinit)
	!%%%%%%%%%%%%
 	! READ DATA
!	call h5_read3(filenamerd, dsetrdname, rdata)									! read dataset
!	call h5_cread3(cdata, cname, bdata, bname, ddata, dname)							! read spline coefficients in T direction
!	adata(:,:,:)= rdata(1:T-1,:,:)											! calculate spline coefficient a in T direction
!	call h5_gread3(filenamerd, "Tarray", "k1array", "lambdaarraylog", Tarray, k1array, loglambdaarray)		! read grid data 
	!%%%%%%%%%%%%
!		call cpu_time(tfin)
!		tspend=tfin-tinit
		!write(*,*) "3D reading data time="
		!write(*,*) tspend
	!%%%%%%%%%%%%



	!%%%%%%%%%%%%
		call cpu_time(tinit)
	!%%%%%%%%%%%%
	! find 3D index
!do l = 1, 140
	call find_idx(hT, T, target(1), Tarray, idx1, xp1)
	call find_idx(hk, k1, target(2), k1array, idx2, xp2)
	call find_idx(hl, loglambda, log(target(3)), loglambdaarray, idx3, xp3)


! interpolate in T dim

	! assing interpolation target matrix (4*4*4) and interpolate in 1D : result = yt1[4*4]
	if (idx1 /= T) then 
		if (idx2 <= k1-3 .AND. idx3 <= loglambda-3) then
			yt1 = ddata(idx1,idx2:idx2+3,idx3:idx3+3)*xp1**3 &
			    + cdata(idx1,idx2:idx2+3,idx3:idx3+3)*xp1**2 &
			    + bdata(idx1,idx2:idx2+3,idx3:idx3+3)*xp1 &
			    + adata(idx1,idx2:idx2+3,idx3:idx3+3)

		else if (idx2 <= k1-3 .AND. idx3 > loglambda-3) then
			yt1 = ddata(idx1,idx2:idx2+3,loglambda-3:loglambda)*xp1**3 &
			    + cdata(idx1,idx2:idx2+3,loglambda-3:loglambda)*xp1**2 &
			    + bdata(idx1,idx2:idx2+3,loglambda-3:loglambda)*xp1 &
			    + adata(idx1,idx2:idx2+3,loglambda-3:loglambda)

		else if (idx2 > k1-3 .AND. idx3 <= loglambda-3) then
			yt1 = ddata(idx1,k1-3:k1,idx3:idx3+3)*xp1**3 &
			    + cdata(idx1,k1-3:k1,idx3:idx3+3)*xp1**2 &
			    + bdata(idx1,k1-3:k1,idx3:idx3+3)*xp1 &
			    + adata(idx1,k1-3:k1,idx3:idx3+3)

		else
			yt1 = ddata(idx1,k1-3:k1,loglambda-3:loglambda)*xp1**3 &
			    + cdata(idx1,k1-3:k1,loglambda-3:loglambda)*xp1**2 &
			    + bdata(idx1,k1-3:k1,loglambda-3:loglambda)*xp1 &
			    + adata(idx1,k1-3:k1,loglambda-3:loglambda)
		end if
	else
		if (idx2 <= k1-3 .AND. idx3 <= loglambda-3) then
			yt1 = rdata(T,idx2:idx2+3,idx3:idx3+3)

		else if (idx2 <= k1-3 .AND. idx3 > loglambda-3) then
			yt1 = rdata(T,idx2:idx2+3,loglambda-3:loglambda)

		else if (idx2 > k1-3 .AND. idx3 <= loglambda-3) then
			yt1 = rdata(T,k1-3:k1,idx3:idx3+3)

		else
			yt1 = rdata(T,k1-3:k1,loglambda-3:loglambda)
		end if

	end if


! interpolate in k1 dim, result : yt2[4]
	call interpolate2(k1, idx2, xp2, yt1, hk, yt2)	

! interpolate in loglambda dim, result : yt3
	call interpolate3(loglambda, idx3, xp3, yt2, hl, yt3)

	Sinterp = yt3
!end do	


!%%%%%%%%%%%%
	call cpu_time(tfin)
	tspend=tfin-tinit
	!write(*,*) "spline3D time="
	!write(*,*) tspend
!%%%%%%%%%%%%
	end subroutine spline3D



	! cubic spline interpolation in second direction
	subroutine interpolate2(n2, idx, xp, yt1, hk, yt2)
		implicit none
		integer		,    intent(in) :: n2, idx					! number of point (n2 = 2d n3 = 3d)
		Double precision,    intent(in) :: xp, hk					
		Double precision,    intent(in) :: yt1(4,4)					! 1st interpolated data array
		Double precision,    intent(out):: yt2(4)					! 2nd interpolated data array
		Double precision 		:: a(3,4), b(3,4), c(4,4), d(3,4)		! spline coefficients

		c(1,:)=0
		c(4,:)=0
		c(2,:)=(3/(4*hk**2))*(yt1(1,:)-2.25*yt1(2,:)+1.5*yt1(3,:)-0.25*yt1(4,:))
		c(3,:)=(3/(4*hk**2))*(yt1(2,:)-2*yt1(3,:)+yt1(4,:))

		a = yt1(1:3,:)
		b = (yt1(2:4,:)-yt1(1:3,:))/hk-hk*(2*c(1:3,:)+c(2:4,:))/3
		d = (c(2:4,:)-c(1:3,:))/(3*hk)

		if (idx .LE. n2-3) then
			yt2(:) = d(1,:)*xp**3 + c(1,:)*xp**2 + b(1,:)*xp + a(1,:)
		else
			if (idx .EQ. n2-2) then
				yt2(:) = d(2,:)*xp**3 + c(2,:)*xp**2 + b(2,:)*xp + a(2,:)
			else if (idx .EQ. n2-1) then
				yt2(:) = d(3,:)*xp**3 + c(3,:)*xp**2 + b(3,:)*xp + a(3,:)
			else
				yt2(:) = yt1(4,:)
			end if
		end if
	end subroutine interpolate2



	! cubic spline interpolation in 3rd dim	(at fixed T, k1 , interpolation in loglambda grid)
	subroutine interpolate3(n3, idx, xp, yt2, hl, yt3)
		implicit none
		integer		,    intent(in) :: n3, idx			! number of point (n2 = 2d n3 = 3d)
		Double precision,    intent(in) :: xp, hl				
		Double precision,    intent(in) :: yt2(4)			! 2nd interpolated data array
		Double precision,    intent(out):: yt3				! final interpolated data array
		Double precision 		:: a(3), b(3), c(4), d(3)	! spline coefficients
		integer				:: k					

		! basic calculation
		c = calc_c(yt2, hl)
		a = calc_a(yt2)
		b = calc_b(yt2, c, hl)
		d = calc_d(c, hl)

		! interpolation
		if (idx .GE. 1 .AND. idx .LE. n3-3) then
			yt3 = d(1)*xp**3 + c(1)*xp**2 + b(1)*xp + a(1)
		else if (idx .EQ. n3-2) then
			yt3 = d(2)*xp**3 + c(2)*xp**2 + b(2)*xp + a(2)
		else if (idx .EQ. n3-1) then
			yt3 = d(3)*xp**3 + c(3)*xp**2 + b(3)*xp + a(3)
		else
			yt3 = yt2(4)
		end if

	end subroutine interpolate3


	! :return    real(8)     c(4): c coefficient
	! Use Tomas algorithm at uniform grid in [4^n] Matrix
	function calc_c(rp, h) result(c)
		implicit none
		Double precision,    intent(in) :: h
		Double precision,    intent(in) :: rp(4)
		Double precision 		:: c(4)

			c(1)=0
			c(2)=(3/(4*h**2))*(rp(1)-2.25*rp(2)+1.5*rp(3)-0.25*rp(4))
			c(3)=(3/(4*h**2))*(rp(2)-2*rp(3)+rp(4))
			c(4)=0
	end function calc_c


	! :return    real(8) a(3): a coefficient
	function calc_a(rp) result(a)
		implicit none
		Double precision,    intent(in) :: rp(4)
		Double precision 		:: a(3)

		a = rp(1:3)

	end function calc_a


	! :return    real(8) b(3): b coefficient
	function calc_b(rp, c, h) result(b)
		implicit none
		Double precision,    intent(in) :: h
		Double precision,    intent(in) :: c(4)
		Double precision,    intent(in) :: rp(4)
		Double precision 		:: b(3)

		b = (rp(2:4)-rp(1:3))/h-h*(2*c(1:3)+c(2:4))/3

	end function calc_b


	! :return    real(8)     d(3): d coefficient
	function calc_d(c, h) result(d)
		implicit none
		Double precision,    intent(in) :: h
		Double precision,    intent(in) :: c(4)
		Double precision 		:: d(3)

		d = (c(2:4)-c(1:3))/(3*h)

	end function calc_d
	
	
	
	! :return   integer idx: index
	! :return   real(8) xp: x-x(index)
	subroutine find_idx(h, n, x, xarr, idx, xp)
		integer		,  intent(out)	:: idx
		Double precision,  intent(out)	:: xp
		integer 	,  intent(in)	:: n
		Double precision,  intent(in)   :: h
		Double precision,  intent(in)   :: x
		Double precision,  intent(in)   :: xarr(1,n)

		idx = int((x-xarr(1,1))/h)+1
		xp  = x - xarr(1,idx)
	end subroutine find_idx
	


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! gyro averaged U calculation
	subroutine Ucalc(ur, urp, uz, uzp, Uperpperp, Uperpparl, Uparlparl, Uperpperpprime, Uparlperpprime, &
			 Tarray, k1array, loglambdaarray, rdata2D, adata2D, bdata2D, cdata2D, ddata2D, &
			 rdata3D, adata3D, bdata3D, cdata3D, ddata3D, k1value, Tvalue)
		IMPLICIT NONE
	Double precision	:: Ninterp
	Double precision	:: N11
	Double precision	:: N12
	Double precision	:: N13
	Double precision	:: N14
	Double precision	:: N15
	Double precision	:: N21
	Double precision	:: N22
	Double precision	:: N23
	Double precision	:: N24
	Double precision	:: N25
!	Double precision	:: N31
	Double precision	:: N32
	Double precision	:: N33
	Double precision	:: N34
	Double precision	:: N35
	Double precision	:: N42
	Double precision	:: N43
	Double precision	:: ur, urp, uz, uzp, k1value, Tvalue
	Double precision, parameter		:: cvel = 299792458		! light velocity
	Double precision, parameter		:: pi = 3.1415926535
	Double precision	:: gamma			
	Double precision	:: gammap
	Double precision, Dimension(4)	:: utable			! utable = ur=1, urp=2, uz=3, uzp=4
	Double precision, Dimension(3)	:: gtable			! gtable = T=1, k1=2, lambda=3
	Double precision, Dimension(5)	:: stable			! stable = Sx1=1 Sx2=2 Sx3=3 Sx4=4 Sx5=5
	Double precision	:: Uperpperp  , Uperpparl, Uparlparl, Uperpperpprime, Uparlperpprime
	Real				:: tinit, tfin, tspend				! calculation time 
	!!=====================================================================================
	integer		,  parameter 	:: k1 = 101, T = 501, loglambda= 136	 				! number of points
	Double precision 		:: Tarray(1,T), k1array(1,k1),  loglambdaarray(1,loglambda)		! grid array for interpolarion
	!Double precision		:: rdata(T,k1), cdata(T,k1)						! dataset and basic spline coefficent buffer
	!Double precision		:: ddata(T-1,k1), bdata(T-1,k1), adata(T-1,k1)				! spline coefficient buffer
	!!=====================================================================================
	Double precision		:: rdata2D(T,k1,10), cdata2D(T,k1,10)		! order of 8 -> 11,12,13,21,22,23,32,33,42,43
	Double precision		:: ddata2D(T-1,k1,10), bdata2D(T-1,k1,10), adata2D(T-1,k1,10)	! spline coefficient buffer
	!!=====================================================================================
	Double precision		:: rdata3D(T,k1,loglambda,6), cdata3D(T,k1,loglambda,6)	! order of 6 -> 14,15,24,25,34,35
	Double precision		:: ddata3D(T-1,k1,loglambda,6), bdata3D(T-1,k1,loglambda,6), adata3D(T-1,k1,loglambda,6)	
 
	utable(1)= ur
	utable(2)= urp
	utable(3)= uz
	utable(4)= uzp
	! if you want to get U kernel with spline remove the below comment
	!call vcalc(utable, gtable, gamma, gammap, stable)


	!write(*,*) "Check1"
	!write(*,*) utable(3)
	!write(*,*) "Check1"
	!%%%%%%%%%%%%
		call cpu_time(tinit)
	!%%%%%%%%%%%%

	!2D interpolation
	k1value = gtable(2)
	Tvalue = gtable(1)

	call spline2D(Tarray, k1array, rdata2D(:,:,1), adata2D(:,:,1), bdata2D(:,:,1), cdata2D(:,:,1), ddata2D(:,:,1), gtable, Ninterp)
	N11 = (Ninterp+stable(1))*gtable(3)**(0.5)
	!write(*,*) 'N11',N11
	call spline2D(Tarray, k1array, rdata2D(:,:,2), adata2D(:,:,2), bdata2D(:,:,2), cdata2D(:,:,2), ddata2D(:,:,2), gtable, Ninterp)
	N12 = (Ninterp+stable(2))*gtable(3)**(-0.5)
	!write(*,*) 'N12',N12
	call spline2D(Tarray, k1array, rdata2D(:,:,3), adata2D(:,:,3), bdata2D(:,:,3), cdata2D(:,:,3), ddata2D(:,:,3), gtable, Ninterp)
	N13 = (Ninterp+stable(3))*gtable(3)**(-1.5)
	!write(*,*) 'N13',N13
	call spline2D(Tarray, k1array, rdata2D(:,:,4), adata2D(:,:,4), bdata2D(:,:,4), cdata2D(:,:,4), ddata2D(:,:,4), gtable, Ninterp)
	N21 = (Ninterp+stable(1))*gtable(3)**(0.5)
	!write(*,*) 'N21',N21
	call spline2D(Tarray, k1array, rdata2D(:,:,5), adata2D(:,:,5), bdata2D(:,:,5), cdata2D(:,:,5), ddata2D(:,:,5), gtable, Ninterp)
	N22 = (Ninterp+stable(2))*gtable(3)**(-0.5)
	!write(*,*) 'N22',N22
	call spline2D(Tarray, k1array, rdata2D(:,:,6), adata2D(:,:,6), bdata2D(:,:,6), cdata2D(:,:,6), ddata2D(:,:,6), gtable, Ninterp)
	N23 = (Ninterp+stable(3))*gtable(3)**(-1.5)
	!write(*,*) 'N23',N23
!!	call spline2D("UperpperpCheck.h5", "Nint31", "Ninterp31", "C31array", "D31array", "B31array", gtable, Ninterp)
!!	N31 = (Ninterp+stable(1))*gtable(1)**(-0.5)

	call spline2D(Tarray, k1array, rdata2D(:,:,7), adata2D(:,:,7), bdata2D(:,:,7), cdata2D(:,:,7), ddata2D(:,:,7), gtable, Ninterp)
	N32 = (Ninterp+stable(2))*gtable(3)**(-0.5)
	!write(*,*) 'N32',N32
	call spline2D(Tarray, k1array, rdata2D(:,:,8), adata2D(:,:,8), bdata2D(:,:,8), cdata2D(:,:,8), ddata2D(:,:,8), gtable, Ninterp)
	N33 = (Ninterp+stable(3))*gtable(3)**(-1.5)
	!write(*,*) 'N33',N33
	call spline2D(Tarray, k1array, rdata2D(:,:,9), adata2D(:,:,9), bdata2D(:,:,9), cdata2D(:,:,9), ddata2D(:,:,9), gtable, Ninterp)
	N42= (Ninterp+stable(2))*gtable(3)**(-0.5)
	!write(*,*) 'N42',N42	
	call spline2D(Tarray, k1array, rdata2D(:,:,10), adata2D(:,:,10), bdata2D(:,:,10), cdata2D(:,:,10), ddata2D(:,:,10), gtable, Ninterp)
	N43 = (Ninterp+stable(3))*gtable(3)**(-1.5)
	!write(*,*) 'N43',N43
	!3D interpolation

	call spline3D(Tarray, k1array, loglambdaarray, rdata3D(:,:,:,1), adata3D(:,:,:,1), bdata3D(:,:,:,1), cdata3D(:,:,:,1), ddata3D(:,:,:,1), gtable, Ninterp)
	N14 = Ninterp+stable(4)
	!write(*,*) 'N14',N14
	call spline3D(Tarray, k1array, loglambdaarray, rdata3D(:,:,:,2), adata3D(:,:,:,2), bdata3D(:,:,:,2), cdata3D(:,:,:,2), ddata3D(:,:,:,2), gtable, Ninterp)
	N15 = Ninterp+stable(5)
	!write(*,*) 'N15',N15
	call spline3D(Tarray, k1array, loglambdaarray, rdata3D(:,:,:,3), adata3D(:,:,:,3), bdata3D(:,:,:,3), cdata3D(:,:,:,3), ddata3D(:,:,:,3), gtable, Ninterp)
	N24 = Ninterp+stable(4)
	!write(*,*) 'N24',N24
	call spline3D(Tarray, k1array, loglambdaarray, rdata3D(:,:,:,4), adata3D(:,:,:,4), bdata3D(:,:,:,4), cdata3D(:,:,:,4), ddata3D(:,:,:,4), gtable, Ninterp)
	N25 = Ninterp+stable(5)
	!write(*,*) 'N25',N25
	call spline3D(Tarray, k1array, loglambdaarray, rdata3D(:,:,:,5), adata3D(:,:,:,5), bdata3D(:,:,:,5), cdata3D(:,:,:,5), ddata3D(:,:,:,5), gtable, Ninterp)
	N34 = Ninterp+stable(4)
	!write(*,*) 'N34',N34
	call spline3D(Tarray, k1array, loglambdaarray, rdata3D(:,:,:,6), adata3D(:,:,:,6), bdata3D(:,:,:,6), cdata3D(:,:,:,6), ddata3D(:,:,:,6), gtable, Ninterp)
	N35 = Ninterp+stable(5)
	!write(*,*) 'N35',N35
	Uperpperp = ( (4*pi) / (gamma*gammap) ) * ( N11 / cvel**2 &
						  + N12 * ( 1 - ( utable(1) / cvel )**2 ) &
						  - N13 * utable(1)**2 &
						  + N24 * 2 * utable(1) * utable(2) / cvel**2 &
						  + N25 * 2 * utable(1) * utable(2) &
						  - N32 * ( utable(2) / cvel )**2 &
						  - N33 * utable(2)**2 )
	
	!write(*,*) "Check1"
	!write(*,*) Uperpperp
	!write(*,*) "Check1"

	Uperpparl = ( (4*pi) / (gamma*gammap) ) * ( - N13 * utable(1) * utable(3) &
						    - N23 * utable(2) * utable(4) &
						    + N25 * utable(3) * utable(2) &
						    + N15 * utable(1) * utable(4) &
						    - N12 * utable(1) * utable(3) / cvel**2 &
						    - N22 * utable(2) * utable(4) / cvel**2 &
						    + N24 * utable(3) * utable(2) / cvel**2 &
						    + N14 * utable(1) * utable(4) / cvel**2 )
	
	
	Uparlparl = ( (4*pi) / (gamma*gammap) ) * ( N11 / cvel**2 &
						  + N12 * ( 1 - utable(3)**2 / cvel**2 - utable(4)**2 / cvel**2 ) &
						  - N13 * ( utable(3)**2 + utable(4)**2 ) &
						  + N14 * 2 * utable(3) * utable(4) / cvel**2 &
						  + N15 * 2 * utable(3) * utable(4) )
	
	


	Uperpperpprime = ( (4*pi) / (gamma*gammap) ) * ( N22 * 1 &
						       - N23 * utable(1)**2 &
						       + N35 * 2 * utable(1) * utable(2) &
						       - N43 * utable(2)**2 &
						       + N15 * utable(1) * utable(2) &
						       - N35 * utable(1) * utable(2) &
						       - N23 * utable(2)**2 &
						       + N43 * utable(2)**2 &
						       + N21 * 1 / cvel**2 &
						       - N22 * utable(1)**2  / cvel**2 &
						       + N34 * 2 * utable(1) * utable(2) / cvel**2 &
						       - N42 * utable(2)**2  / cvel**2 &
						       + N14 * utable(1) * utable(2) / cvel**2 &
						       - N34 * utable(1) * utable(2) / cvel**2 &
						       - N22 * utable(2)**2  / cvel**2 &
						       + N42 * utable(2)**2  / cvel**2 )
	

	Uparlperpprime = ( (4*pi) / (gamma*gammap) ) * ( N15 * utable(3) * utable(2) &
						       - N35 * utable(3) * utable(2) &
						       - N13 * utable(2) * utable(4) &
						       + N33 * utable(2) * utable(4) &
						       - N23 * utable(1) * utable(3) &
						       - N33 * utable(2) * utable(4) &
						       + N35 * utable(2) * utable(3) &
						       + N25 * utable(1) * utable(4) &
						       + N14 * utable(3) * utable(2) / cvel**2 &
						       - N34 * utable(3) * utable(2) / cvel**2 &
						       - N12 * utable(2) * utable(4) / cvel**2 &
						       + N32 * utable(2) * utable(4) / cvel**2 &
						       - N22 * utable(1) * utable(3) / cvel**2 &
						       - N32 * utable(2) * utable(4) / cvel**2 &
						       + N32 * utable(2) * utable(3) / cvel**2 &
						       + N24 * utable(1) * utable(4) / cvel**2 )
	
	end subroutine Ucalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPLINE FINISH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SPLINE FINISH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module spline_interpolation
