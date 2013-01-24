! Fourth-Order Runge-Kutta Algorithm to solve ODE system: X'(t) = F(X,t)
! Last modified by GuuD on Jan. 9, 2013

! Module:
! M_qsort: Quicksort
! M_rk4_parameter
! M_derivative_linear_ODE
! M_rk4: Main computing.

module M_qsort!{{{
    implicit none
    private
    public :: S_qsort_real
    contains
    recursive subroutine S_qsort_real ( A )!{{{
        implicit none
        ! Argument
        real(kind=8) , intent(inout) , dimension(:) :: A
        ! Local
        integer(kind=4) :: marker
        ! Core
        if ( size(A) > 1) then
            call Partition ( A, marker )
            call S_qsort_real ( A ( : marker-1 ) )
            call S_qsort_real ( A ( marker : ) )
        end if
    end subroutine S_qsort_real!}}}
    subroutine Partition ( A , marker )!{{{
        implicit none
        ! Argument
        real(kind=8) , intent(inout) , dimension(:) :: A
        integer(kind=4) , intent(out) :: marker
        ! Local
        integer(kind=4) :: i
        integer(kind=4) :: j
        real(kind=8) :: temp
        real(kind=8) :: x
        ! Initialization
        x = A(1)
        i = 0
        j = size(A) + 1
        ! Loop
        Find_Marker: do
            j = j-1
            do
                if ( A(j) .le. x ) exit
                j = j-1
            end do
            i = i+1
            do
                if ( A(i) >= x ) exit
                i = i+1
            end do
            if ( i < j ) then
                temp = A(i)
                A(i) = A(j)
                A(j) = temp
            elseif (i == j) then
                marker = i+1
                return
            else
                marker = i
                return
            endif
        end do Find_Marker
    end subroutine Partition!}}}
end module M_qsort!}}}

module M_rk4_parameter!{{{
    implicit none
    private
    public :: Nequation
    protected :: Nequation
    public :: Ntime
    protected :: Ntime
    public :: timestep
    protected :: timestep
    public :: TIMEPOINT
    protected :: TIMEPOINT
    public :: S_set_Nequation
    public :: S_set_Ntime
    public :: S_set_timestep
    public :: S_set_TIMEPOINT
    ! Variable
    integer(kind=4) :: Nequation
    integer(kind=4) :: Ntime
    real(kind=8) :: timestep
    real(kind=8) , dimension(:) , pointer :: TIMEPOINT => null()
    contains
    subroutine S_set_Nequation ( v )!{{{
        implicit none
        integer(kind=4) , intent(in) :: v
        if ( v .le. 0 ) then
            write(*,*) 'ERROR: Number of equation must be positive.'
            stop
        end if
        Nequation = v
        return
    end subroutine S_set_Nequation!}}}
    subroutine S_set_Ntime ( v )!{{{
        implicit none
        integer(kind=4) , intent(in) :: v
        if ( v .le. 0 ) then
            write(*,*) 'ERROR: Number of time-point must be positive.'
            stop
        end if
        Ntime = v
        return
    end subroutine S_set_Ntime!}}}
    subroutine S_set_timestep ( v )!{{{
        implicit none
        real(kind=8) , intent(in) :: v
        integer(kind=4) :: i
        if ( v .le. 0 ) then
            write(*,*) 'ERROR: Time-step for RK4 must be positive.'
            stop
        end if
        timestep = v
        return
    end subroutine S_set_timestep!}}}
    subroutine S_set_TIMEPOINT ( V )!{{{
        use M_qsort ,&
            only : S_qsort_real
        implicit none
        real(kind=8) , dimension(:) , intent(inout) , target :: V
        integer(kind=4) :: i
        if ( size(V) .ne. Ntime ) then
            write(*,*) 'ERROR: Size of Time-point array not equal to Number &
            of time-point.'
            stop
        end if
        Ordering_Check: do i=1,size(V)-1
            if ( V(i) .eq. V(i+1) ) then
                write(*,*) 'ERROR: Time-point array contains identical items.'
                stop
            end if
            if ( V(i) .gt. V(i+1) ) then
                write(*,*) 'WARNING: Time-point array not in ascending order. &
                    Sorting performed.'
                call S_qsort_real ( V )
                exit
            end if
        end do Ordering_Check
        Identity_Check: do i=1,size(V)-1
            if ( V(i) .eq. V(i+1) ) then
                write(*,*) 'ERROR: Time-point array contains identical items.'
                stop
            end if
        end do Identity_Check
        TIMEPOINT => V
        return
    end subroutine S_set_TIMEPOINT!}}}
end module M_rk4_parameter!}}}

module M_derivative_linear_ODE!{{{
    implicit none
    private
    public :: Nnonzero
    protected :: Nnonzero
!    public :: NONZERO_ROW
!    protected :: NONZERO_ROW
!    public :: NONZERO_COLUMN
!    protected :: NONZERO_COLUMN
!    public :: NONZERO_VALUE
!    protected :: NONZERO_VALUE
    public :: S_set_Nnonzero
    public :: S_set_NONZERO_ROW
    public :: S_set_NONZERO_COLUMN
    public :: S_set_NONZERO_VALUE
    public :: S_set_INTERCEPT
    public :: S_derivative_linear_ODE
    ! Variable
    integer(kind=4) :: Nnonzero
    integer(kind=4) , dimension(:) , pointer :: NONZERO_ROW => null()
    integer(kind=4) , dimension(:) , pointer :: NONZERO_COLUMN => null()
    real(kind=8) , dimension(:) , pointer :: NONZERO_VALUE => null()
    real(kind=8) , dimension(:) , pointer :: INTERCEPT => null()
    contains
    subroutine S_set_Nnonzero ( v )!{{{
        use M_rk4_parameter ,&
            only : Nequation
        implicit none
        integer(kind=4) , intent(in) :: v
        Nnonzero = v
        if ( v .le. 0 .or. v .ge. Nequation**2 ) then
            write(*,*) 'WARNING: Parameter matrix is full.'
            Nnonzero = Nequation**2
        end if
        return
    end subroutine S_set_Nnonzero!}}}
    subroutine S_set_NONZERO_ROW ( V )!{{{
        use M_rk4_parameter ,&
            only : Nequation
        implicit none
        integer(kind=4) , dimension(:) , intent(in) , target :: V
        if ( size(V) .ne. Nnonzero ) then
            write(*,*) 'ERROR: Size of Nonzero-row-position array not &
                equal to Number of nonzero entries.'
            stop
        end if
        if ( minval(V) .le. 0 .or. maxval(V) .gt. Nequation ) then
            write(*,*) 'Inappropriate row number in &
                Nonzero-row-position array.'
            write(*,*) 'Minimal value:' , minval(V)
            write(*,*) 'Maximal value:' , maxval(V)
            stop
        end if
        NONZERO_ROW => V
!        write(*,*) size(NONZERO_ROW) , NONZERO_ROW
        return
    end subroutine S_set_NONZERO_ROW!}}}
    subroutine S_set_NONZERO_COLUMN ( V )!{{{
        use M_rk4_parameter ,&
            only : Nequation
        implicit none
        integer(kind=4) , dimension(:) , intent(in) , target :: V
        if ( size(V) .ne. Nnonzero ) then
            write(*,*) 'ERROR: Size of Nonzero-column-position array not &
                equal to Number of nonzero entries.'
            stop
        end if
        if ( minval(V) .le. 0 .or. maxval(V) .gt. Nequation ) then
            write(*,*) 'Inappropriate column number in &
                Nonzero-column-position array.'
            write(*,*) 'Minimal value:' , minval(V)
            write(*,*) 'Maximal value:' , maxval(V)
            stop
        end if
        NONZERO_COLUMN => V
        return
    end subroutine S_set_NONZERO_COLUMN!}}}
    subroutine S_set_NONZERO_VALUE ( V )!{{{
        implicit none
        real(kind=8) , dimension(:) , intent(in) , target :: V
        if ( size(V) .ne. Nnonzero ) then
            write(*,*) 'ERROR: Size of Nonzero-entry-value array not &
                equal to Number of nonzero entries.'
            write(*,*) 'Array size:' , size(V)
            stop
        end if
        NONZERO_VALUE => V
        return
    end subroutine S_set_NONZERO_VALUE!}}}
    subroutine S_set_INTERCEPT ( V )!{{{
        use M_rk4_parameter ,&
            only : Nequation
        implicit none
        real(kind=8) , dimension(:) , intent(in) , target :: V
        if ( size(V) .ne. Nequation ) then
            write(*,*) 'ERROR: Size of intercept array not &
                equal to Number of equation.'
            stop
        end if
        INTERCEPT => V
        return
    end subroutine S_set_INTERCEPT!}}}
    subroutine S_derivative_linear_ODE ( D , X , t )!{{{
        use M_rk4_parameter ,&
            only : Nequation
        implicit none
        real(kind=8) , intent(out) , dimension(:) :: D
        real(kind=8) , intent(in) , dimension(:) :: X
        real(kind=8) , intent(in) :: t
        integer(kind=4) :: i
        if ( size(D) .ne. Nequation ) then
            write(*,*) 'ERROR: Size of derivative output not &
                equal to Number of equation.'
            stop
        end if
        if ( size(X) .ne. Nequation ) then
            write(*,*) 'ERROR: Size of derivative input not &
                equal to Number of equation.'
            stop
        end if
        D = 0
        do i=1,Nnonzero
            D (NONZERO_ROW(i)) = D (NONZERO_ROW(i)) +&
                NONZERO_VALUE(i) * X (NONZERO_COLUMN(i))
        end do
        D = D + INTERCEPT
        return
    end subroutine S_derivative_linear_ODE!}}}
end module M_derivative_linear_ODE!}}}

module M_rk4!{{{
    implicit none
    private
    public :: S_set_INITIAL
    public :: S_rk4
    ! Variable
    real(kind=8) , dimension(:) , pointer :: INITIAL => null()
    contains
    subroutine S_set_INITIAL ( V )!{{{
        use M_rk4_parameter ,&
            only : Nequation
        implicit none
        real(kind=8) , dimension(:) , intent(in) , target :: V
        if ( size(V) .ne. Nequation ) then
            write(*,*) 'ERROR: Size of initial-value array not &
                equal to Number of equation.'
            stop
        end if
        INITIAL => V
        return
    end subroutine S_set_INITIAL!}}}
    subroutine S_rk4 ( ODERESULT , f )!{{{
        use M_rk4_parameter ,&
            only : Nequation , Ntime , timestep , TIMEPOINT
        ! Argument
        implicit none
        real(kind=8) , intent(out) , dimension(:,:) :: ODERESULT
        external :: f
        interface
            subroutine f ( D , X , t )
                implicit none
                real(kind=8) , intent(in) :: t
                real(kind=8) , intent(in) , dimension(:) :: X
                real(kind=8) , intent(out) , dimension(:) :: D
            end subroutine f
        end interface
        ! Local Variable
        integer(kind=4) :: i
        integer(kind=4) :: j
        integer(kind=4) :: Nstep
        real(kind=8) :: step
        real(kind=8) :: currenttime
        real(kind=8) , allocatable , dimension(:,:) , save :: TEMP
            ! Column 1: x(t_{i+1}) intermediate value.
            ! Column 2: x(t_i).
            ! Column 3: derivative evaluating argument.
            ! Column 4: k_1,k_2,k_3,k_4 in RK4
        if ( size(ODERESULT,1) .ne. Nequation &
            .or. size(ODERESULT,2) .ne. Ntime ) then
            write(*,*) 'ERROR: Size of ODE-result array not equal &
                to Number of (Equation,Time-point).'
            stop
        end if
        if ( .not. ALLOCATED ( TEMP ) ) then
            ALLOCATE ( TEMP ( Nequation , 4 ) )
        end if
        ! Iteration
        ODERESULT(:,1) = INITIAL
        TEMP(:,1) = INITIAL
        TEMP(:,2) = INITIAL
        currenttime = TIMEPOINT(1)
        do i = 2 , Ntime
            step = TIMEPOINT(i) - TIMEPOINT(i-1)
            Nstep = CEILING ( step / timestep )
            step = step / Nstep
            ! Inner Loop
            do j = 1 , Nstep
                if ( j == Nstep ) then
                    step = TIMEPOINT(i) - currenttime
                end if
                CALL f ( TEMP(:,4) , TEMP(:,2) , currenttime )
                TEMP(:,4) = TEMP(:,4) * step
                TEMP(:,1) = TEMP(:,1) + TEMP(:,4) / 6
                TEMP(:,3) = TEMP(:,2) + TEMP(:,4) / 2
                currenttime = currenttime + step / 2
                CALL f ( TEMP(:,4) , TEMP(:,3) , currenttime )
                TEMP(:,4) = TEMP(:,4) * step
                TEMP(:,1) = TEMP(:,1) + TEMP(:,4) / 3
                TEMP(:,3) = TEMP(:,2) + TEMP(:,4) / 2
                CALL f ( TEMP(:,4) , TEMP(:,3) , currenttime )
                TEMP(:,4) = TEMP(:,4) * step
                TEMP(:,1) = TEMP(:,1) + TEMP(:,4) / 3
                TEMP(:,3) = TEMP(:,2) + TEMP(:,4)
                currenttime = currenttime + step / 2
                CALL f ( TEMP(:,4) , TEMP(:,3) , currenttime )
                TEMP(:,4) = TEMP(:,4) * step
                TEMP(:,1) = TEMP(:,1) + TEMP(:,4) / 6
                TEMP(:,2) = TEMP(:,1)
            end do
            ODERESULT(:,i) = TEMP(:,1)
        end do
        return
    end subroutine S_rk4!}}}
end module M_rk4!}}}
module M_least_squares!{{{
    implicit none
    private
    public :: S_set_Nobservation
    public :: S_set_OBSERVATION_EQUATION
    public :: S_set_OBSERVATION_TIME
    public :: S_set_OBSERVATION_VALUE
    public :: S_set_WEIGHT
    public :: S_ODE_linear_ls
    public :: S_ODE_linear_ls_lmdif
    ! Variable
    integer(kind=4) :: Nobservation
    integer(kind=4) , dimension(:) , pointer :: OBSERVATION_EQUATION
    integer(kind=4) , dimension(:) , allocatable :: OBSERVATION_TIME
    real(kind=8) , dimension(:) , pointer :: OBSERVATION_VALUE
    real(kind=8) , dimension(:) , pointer :: WEIGHT
    contains
    subroutine S_set_Nobservation ( v )!{{{
        implicit none
        integer(kind=4) , intent(in) :: v
        if ( v .le. 0 ) then
            write(*,*) 'ERROR: Number of observations must be positive.'
            stop
        end if
        Nobservation = v
        return
    end subroutine S_set_Nobservation!}}}
    subroutine S_set_OBSERVATION_EQUATION ( V )!{{{
        use M_rk4_parameter ,&
            only : Nequation
        implicit none
        integer(kind=4) , dimension(:) , intent(in) , target :: V
        if ( size(V) .ne. Nobservation ) then
            write(*,*) 'ERROR: Size of observation-equation array not equal &
                to number of observations.'
            stop
        end if
        if ( minval(V) .lt. 1 .or. maxval(V) .gt. Nequation ) then
            write(*,*) 'ERROR: Inappropriate value in &
                observations-equation array.'
            stop
        end if
        OBSERVATION_EQUATION => V
        return
    end subroutine S_set_OBSERVATION_EQUATION!}}}
    subroutine S_set_OBSERVATION_TIME ( V )!{{{
        use M_rk4_parameter ,&
            only : Ntime , TIMEPOINT
        implicit none
        real(kind=8) , dimension(:) , intent(in) , target :: V
        integer(kind=4) :: i
        integer(kind=4) :: j
        if ( size(V) .ne. Nobservation ) then
            write(*,*) 'ERROR: Size of observation-equation array not equal &
                to number of observations.'
            stop
        end if
        if ( allocated(OBSERVATION_TIME) .and. &
            size(OBSERVATION_TIME) .ne.  Nobservation ) then
            deallocate ( OBSERVATION_TIME )
        end if
        if ( .not. allocated(OBSERVATION_TIME) ) then
            allocate ( OBSERVATION_TIME (Nobservation) , stat = i )
            if ( i .ne. 0 ) then
                write(*,*) 'ERROR: Unable to allocate memory for &
                    observation-time array.'
                stop
            end if
        end if
        do i = 1 , Nobservation
            Inner: do j = 1 , Ntime
                if ( abs(V(i)-TIMEPOINT(j)) .lt. 1e-3 ) then
                    OBSERVATION_TIME(i) = j
                    exit Inner
                end if
                if ( V(i) .lt. TIMEPOINT(j) - 1e-3 ) then
                    write(*,*) 'ERROR: Item in observation-time array not &
                        contained in time-point array.'
                    stop
                end if
            end do Inner
        end do
        return
    end subroutine S_set_OBSERVATION_TIME!}}}
    subroutine S_set_OBSERVATION_VALUE ( V )!{{{
        implicit none
        real(kind=8) , dimension(:) , intent(in) , target :: V
        if ( size(V) .ne. Nobservation ) then
            write(*,*) 'ERROR: Size of observation-equation array not equal &
                to number of observations.'
            stop
        end if
        OBSERVATION_VALUE => V
        return
    end subroutine S_set_OBSERVATION_VALUE!}}}
    subroutine S_set_WEIGHT ( V )!{{{
        implicit none
        real(kind=8) , dimension(:) , intent(in) , target :: V
        if ( size(V) .ne. Nobservation ) then
            write(*,*) 'ERROR: Size of weight array not equal &
                to number of observations.'
            stop
        end if
        WEIGHT => V
        return
    end subroutine S_set_WEIGHT!}}}
    subroutine S_ODE_linear_ls ( Y , X )!{{{
        use M_rk4_parameter ,&
            only : Nequation , Ntime! , timestep , TIMEPOINT
        use M_derivative_linear_ODE ,&
            only : Nnonzero , S_set_NONZERO_VALUE , S_set_INTERCEPT ,&
                S_derivative_linear_ODE
        use M_rk4 ,&
            only : S_set_INITIAL , S_rk4
        implicit none
        ! Argument
        real(kind=8) , dimension(:) , intent(out) :: Y
        real(kind=8) , dimension(:) , intent(in) :: X
        ! Local
        real(kind=8) , dimension(:,:) , allocatable :: ODERESULT
        integer(kind=4) :: i
!        real(kind=8) , dimension(:) , allocatable , save :: T_WEIGHT
!        if ( .not. allocated(T_WEIGHT) ) allocate ( T_WEIGHT(Nequation) )
!        T_WEIGHT = (/8,16,250,120,250000,500000,120000,35000,120000,30000/)
        ! Size check
        if ( size(Y) .ne. Nobservation ) then
            write(*,*) 'ERROR: Size of least-squares function output &
                variable not equal to number of observations.'
            write(*,*) 'Output variable size:' , size(Y)
            write(*,*) 'Number of observation:' , Nobservation
            stop
        end if
        if ( size(X) .ne. (Nnonzero+2*Nequation) ) then
            write(*,*) 'ERROR: Incorrect size of &
                least-squares input variable size.'
            stop
        end if
        ! Initialization
        allocate ( ODERESULT ( Nequation , Ntime ) , stat = i )
        if ( i .ne. 0 ) then
            write(*,*) 'ERROR: Unable to allocate memory for &
                ODE-result array.'
            stop
        end if
        call S_set_INITIAL ( X ( 1:Nequation ) )
        call S_set_INTERCEPT ( X ( (Nequation+1):(2*Nequation) ) )
        call S_set_NONZERO_VALUE ( X ( (2*Nequation+1): ) )
        ! RK4 to solve ODE
        call S_rk4 ( ODERESULT , S_derivative_linear_ODE )
        ! Output assignment
        do i = 1 , Nobservation
            Y(i) = ( OBSERVATION_VALUE(i) - &
                ODERESULT(OBSERVATION_EQUATION(i),OBSERVATION_TIME(i)) ) &
                * WEIGHT(i)
!                / T_WEIGHT(OBSERVATION_EQUATION(i))
        end do
        return
    end subroutine S_ODE_linear_ls!}}}
    subroutine S_ODE_linear_ls_lmdif ( m , n , X , Y , i )!{{{
        implicit none
        ! Argument
        INTEGER, INTENT(IN)      :: m, n
        REAL (8), INTENT(IN)    :: x(:)
        REAL (8), INTENT(OUT)   :: Y(:)
        INTEGER, INTENT(IN OUT)  :: i
!        integer(kind=4) , intent(in) :: m
!        integer(kind=4) , intent(in) :: n
!        real(kind=8) , dimension(:) , intent(in) :: X
!        real(kind=8) , dimension(:) , intent(out) :: Y
!        integer(kind=4) , intent(inout) :: i
        if ( i .lt. 0 ) then
            write(*,*) 'ERROR: User abortion.'
            stop
        end if
        call S_ODE_linear_ls ( Y , X )
        if ( i .eq. 0 ) then
            write(*,*) 'Current RSS: ' , norm2(Y)**2
        end if
        return
    end subroutine S_ODE_linear_ls_lmdif!}}}
end module M_least_squares!}}}
MODULE M_Levenberg_Marquardt
! MINPACK routines which are used by both LMDIF & LMDER

IMPLICIT NONE
!INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER, PARAMETER :: dp = 8

PRIVATE
!PUBLIC :: dp, lmdif1, lmdif, lmder1, lmder, enorm
public :: lmdif

CONTAINS


SUBROUTINE lmdif1(fcn, m, n, x, fvec, tol, info, iwa)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-11  Time: 00:51:44

! N.B. Arguments WA & LWA have been removed.

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(OUT)     :: fvec(:)
REAL (dp), INTENT(IN)      :: tol
INTEGER, INTENT(OUT)       :: info
INTEGER, INTENT(OUT)       :: iwa(:)

! EXTERNAL fcn

INTERFACE
  SUBROUTINE fcn(m, n, x, fvec, iflag)
    IMPLICIT NONE
!    INTEGER, PARAMETER       :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, PARAMETER       :: dp = 8
    INTEGER, INTENT(IN)      :: m, n
    REAL (dp), INTENT(IN)    :: x(:)
    REAL (dp), INTENT(OUT)   :: fvec(:)
    INTEGER, INTENT(IN OUT)  :: iflag
  END SUBROUTINE fcn
END INTERFACE

!     **********

!     subroutine lmdif1

!     the purpose of lmdif1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of the
!     levenberg-marquardt algorithm. this is done by using the more
!     general least-squares solver lmdif. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.

!     the subroutine statement is

!       subroutine lmdif1(fcn, m, n, x, fvec, tol, info, iwa)

!     where

!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.

!         subroutine fcn(m, n, x, fvec, iflag)
!         integer m, n, iflag
!         REAL (dp) x(n), fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end

!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmdif1.
!         in this case set iflag to a negative integer.

!       m is a positive integer input variable set to the number of functions.

!       n is a positive integer input variable set to the number
!         of variables.  n must not exceed m.

!       x is an array of length n.  on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.

!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.

!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at most tol.

!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.

!         info = 0  improper input parameters.

!         info = 1  algorithm estimates that the relative error
!                   in the sum of squares is at most tol.

!         info = 2  algorithm estimates that the relative error
!                   between x and the solution is at most tol.

!         info = 3  conditions for info = 1 and info = 2 both hold.

!         info = 4  fvec is orthogonal to the columns of the
!                   jacobian to machine precision.

!         info = 5  number of calls to fcn has reached or
!                   exceeded 200*(n+1).

!         info = 6  tol is too small. no further reduction in
!                   the sum of squares is possible.

!         info = 7  tol is too small. no further improvement in
!                   the approximate solution x is possible.

!       iwa is an integer work array of length n.

!       wa is a work array of length lwa.

!       lwa is a positive integer input variable not less than
!         m*n+5*n+m.

!     subprograms called

!       user-supplied ...... fcn

!       minpack-supplied ... lmdif

!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more

!     **********
INTEGER   :: maxfev, mode, nfev, nprint
REAL (dp) :: epsfcn, ftol, gtol, xtol, wa(2*n), fjac(m,n)
REAL (dp), PARAMETER :: factor = 100._dp, zero = 0.0_dp

info = 0

!     check the input parameters for errors.

IF (n <= 0 .OR. m < n .OR. tol < zero) GO TO 10

!     call lmdif.

maxfev = 200*(n + 1)
ftol = tol
xtol = tol
gtol = zero
epsfcn = zero
mode = 1
nprint = 0
CALL lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, wa,   &
           mode, factor, nprint, info, nfev, fjac, iwa, wa(n+1:))
IF (info == 8) info = 4

10 RETURN

!     last card of subroutine lmdif1.

END SUBROUTINE lmdif1



SUBROUTINE lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,  &
                 diag, mode, factor, nprint, info, nfev, fjac, ipvt, qtf)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 12:45:59

! N.B. Arguments LDFJAC, WA1, WA2, WA3 & WA4 have been removed.

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(OUT)     :: fvec(:)
REAL (dp), INTENT(IN)      :: ftol
REAL (dp), INTENT(IN)      :: xtol
REAL (dp), INTENT(IN OUT)  :: gtol
INTEGER, INTENT(IN OUT)    :: maxfev
REAL (dp), INTENT(IN OUT)  :: epsfcn
REAL (dp), INTENT(OUT)     :: diag(:)
INTEGER, INTENT(IN)        :: mode
REAL (dp), INTENT(IN)      :: factor
INTEGER, INTENT(IN)        :: nprint
INTEGER, INTENT(OUT)       :: info
INTEGER, INTENT(OUT)       :: nfev
REAL (dp), INTENT(OUT)     :: fjac(:,:)    ! fjac(ldfjac,n)
INTEGER, INTENT(OUT)       :: ipvt(:)
REAL (dp), INTENT(OUT)     :: qtf(:)

! EXTERNAL fcn

INTERFACE
  SUBROUTINE fcn(m, n, x, fvec, iflag)
    IMPLICIT NONE
!    INTEGER, PARAMETER       :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, PARAMETER       :: dp = 8
    INTEGER, INTENT(IN)      :: m, n
    REAL (dp), INTENT(IN)    :: x(:)
    REAL (dp), INTENT(OUT)   :: fvec(:)
    INTEGER, INTENT(IN OUT)  :: iflag
  END SUBROUTINE fcn
END INTERFACE

!     **********

!     subroutine lmdif

!     the purpose of lmdif is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.

!     the subroutine statement is

!       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
!                        diag,mode,factor,nprint,info,nfev,fjac,
!                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)

!     where

!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.

!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         REAL (dp) x(:),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end

!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmdif.
!         in this case set iflag to a negative integer.

!       m is a positive integer input variable set to the number
!         of functions.

!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.

!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.

!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.

!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.

!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.

!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.

!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least
!         maxfev by the end of an iteration.

!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.

!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.

!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.

!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.

!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.

!       info is an integer output variable.  If the user has terminated
!         execution, info is set to the (negative) value of iflag.
!         See description of fcn.  Otherwise, info is set as follows.

!         info = 0  improper input parameters.

!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.

!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.

!         info = 3  conditions for info = 1 and info = 2 both hold.

!         info = 4  the cosine of the angle between fvec and any column of
!                   the Jacobian is at most gtol in absolute value.

!         info = 5  number of calls to fcn has reached or exceeded maxfev.

!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.

!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.

!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.

!       nfev is an integer output variable set to the number of calls to fcn.

!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that

!                t     t           t
!               p *(jac *jac)*p = r *r,

!         where p is a permutation matrix and jac is the final calculated
!         Jacobian.  Column j of p is column ipvt(j) (see below) of the
!         identity matrix. the lower trapezoidal part of fjac contains
!         information generated during the computation of r.

!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.

!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.

!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.

!       wa1, wa2, and wa3 are work arrays of length n.

!       wa4 is a work array of length m.

!     subprograms called

!       user-supplied ...... fcn

!       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac

!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod

!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more

!     **********
INTEGER   :: i, iflag, iter, j, l
REAL (dp) :: actred, delta, dirder, epsmch, fnorm, fnorm1, gnorm,  &
             par, pnorm, prered, ratio, sum, temp, temp1, temp2, xnorm
REAL (dp) :: wa1(n), wa2(n), wa3(n), wa4(m)
REAL (dp), PARAMETER :: one = 1.0_dp, p1 = 0.1_dp, p5 = 0.5_dp,  &
                        p25 = 0.25_dp, p75 = 0.75_dp, p0001 = 0.0001_dp, &
                        zero = 0.0_dp

!     epsmch is the machine precision.

epsmch = EPSILON(zero)

info = 0
iflag = 0
nfev = 0

!     check the input parameters for errors.

IF (n <= 0 .OR. m < n .OR. ftol < zero .OR. xtol < zero .OR. gtol < zero  &
    .OR. maxfev <= 0 .OR. factor <= zero) GO TO 300
IF (mode /= 2) GO TO 20
DO  j = 1, n
  IF (diag(j) <= zero) GO TO 300
END DO

!     evaluate the function at the starting point and calculate its norm.

20 iflag = 1
CALL fcn(m, n, x, fvec, iflag)
nfev = 1
IF (iflag < 0) GO TO 300
fnorm = enorm(m, fvec)

!     initialize levenberg-marquardt parameter and iteration counter.

par = zero
iter = 1

!     beginning of the outer loop.


!        calculate the jacobian matrix.

30 iflag = 2
CALL fdjac2(fcn, m, n, x, fvec, fjac, iflag, epsfcn)
nfev = nfev + n
IF (iflag < 0) GO TO 300

!        if requested, call fcn to enable printing of iterates.

IF (nprint <= 0) GO TO 40
iflag = 0
IF (MOD(iter-1,nprint) == 0) CALL fcn(m, n, x, fvec, iflag)
IF (iflag < 0) GO TO 300

!        compute the qr factorization of the jacobian.

40 CALL qrfac(m, n, fjac, .true., ipvt, wa1, wa2)

!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.

IF (iter /= 1) GO TO 80
IF (mode == 2) GO TO 60
DO  j = 1, n
  diag(j) = wa2(j)
  IF (wa2(j) == zero) diag(j) = one
END DO

!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.

60 wa3(1:n) = diag(1:n)*x(1:n)
xnorm = enorm(n, wa3)
delta = factor*xnorm
IF (delta == zero) delta = factor

!        form (q transpose)*fvec and store the first n components in qtf.

80 wa4(1:m) = fvec(1:m)
DO  j = 1, n
  IF (fjac(j,j) == zero) GO TO 120
  sum = DOT_PRODUCT( fjac(j:m,j), wa4(j:m) )
  temp = -sum/fjac(j,j)
  DO  i = j, m
    wa4(i) = wa4(i) + fjac(i,j)*temp
  END DO
  120 fjac(j,j) = wa1(j)
  qtf(j) = wa4(j)
END DO

!        compute the norm of the scaled gradient.

gnorm = zero
IF (fnorm == zero) GO TO 170
DO  j = 1, n
  l = ipvt(j)
  IF (wa2(l) == zero) CYCLE
  sum = zero
  DO  i = 1, j
    sum = sum + fjac(i,j)*(qtf(i)/fnorm)
  END DO
  gnorm = MAX(gnorm, ABS(sum/wa2(l)))
END DO

!        test for convergence of the gradient norm.

170 IF (gnorm <= gtol) info = 4
IF (info /= 0) GO TO 300

!        rescale if necessary.

IF (mode == 2) GO TO 200
DO  j = 1, n
  diag(j) = MAX(diag(j), wa2(j))
END DO

!        beginning of the inner loop.

!           determine the levenberg-marquardt parameter.

200 CALL lmpar(n, fjac, ipvt, diag, qtf, delta, par, wa1, wa2)

!           store the direction p and x + p. calculate the norm of p.

DO  j = 1, n
  wa1(j) = -wa1(j)
  wa2(j) = x(j) + wa1(j)
  wa3(j) = diag(j)*wa1(j)
END DO
pnorm = enorm(n, wa3)

!           on the first iteration, adjust the initial step bound.

IF (iter == 1) delta = MIN(delta, pnorm)

!           evaluate the function at x + p and calculate its norm.

iflag = 1
CALL fcn(m, n, wa2, wa4, iflag)
nfev = nfev + 1
IF (iflag < 0) GO TO 300
fnorm1 = enorm(m, wa4)

!           compute the scaled actual reduction.

actred = -one
IF (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

!           compute the scaled predicted reduction and
!           the scaled directional derivative.

DO  j = 1, n
  wa3(j) = zero
  l = ipvt(j)
  temp = wa1(l)
  DO  i = 1, j
    wa3(i) = wa3(i) + fjac(i,j)*temp
  END DO
END DO
temp1 = enorm(n,wa3)/fnorm
temp2 = (SQRT(par)*pnorm)/fnorm
prered = temp1**2 + temp2**2/p5
dirder = -(temp1**2 + temp2**2)

!           compute the ratio of the actual to the predicted reduction.

ratio = zero
IF (prered /= zero) ratio = actred/prered

!           update the step bound.

IF (ratio <= p25) THEN
  IF (actred >= zero) temp = p5
  IF (actred < zero) temp = p5*dirder/(dirder + p5*actred)
  IF (p1*fnorm1 >= fnorm .OR. temp < p1) temp = p1
  delta = temp*MIN(delta,pnorm/p1)
  par = par/temp
ELSE
  IF (par /= zero .AND. ratio < p75) GO TO 260
  delta = pnorm/p5
  par = p5*par
END IF

!           test for successful iteration.

260 IF (ratio < p0001) GO TO 290

!           successful iteration. update x, fvec, and their norms.

DO  j = 1, n
  x(j) = wa2(j)
  wa2(j) = diag(j)*x(j)
END DO
fvec(1:m) = wa4(1:m)
xnorm = enorm(n, wa2)
fnorm = fnorm1
iter = iter + 1

!           tests for convergence.

290 IF (ABS(actred) <= ftol .AND. prered <= ftol .AND. p5*ratio <= one) info = 1
IF (delta <= xtol*xnorm) info = 2
IF (ABS(actred) <= ftol .AND. prered <= ftol  &
    .AND. p5*ratio <= one .AND. info == 2) info = 3
IF (info /= 0) GO TO 300

!           tests for termination and stringent tolerances.

IF (nfev >= maxfev) info = 5
IF (ABS(actred) <= epsmch .AND. prered <= epsmch  &
    .AND. p5*ratio <= one) info = 6
IF (delta <= epsmch*xnorm) info = 7
IF (gnorm <= epsmch) info = 8
IF (info /= 0) GO TO 300

!           end of the inner loop. repeat if iteration unsuccessful.

IF (ratio < p0001) GO TO 200

!        end of the outer loop.

GO TO 30

!     termination, either normal or user imposed.

300 IF (iflag < 0) info = iflag
iflag = 0
IF (nprint > 0) CALL fcn(m, n, x, fvec, iflag)
RETURN

!     last card of subroutine lmdif.

END SUBROUTINE lmdif



!SUBROUTINE lmder1(fcn, m, n, x, fvec, fjac, tol, info, ipvt)
! 
!! Code converted using TO_F90 by Alan Miller
!! Date: 1999-12-09  Time: 12:45:54
!
!! N.B. Arguments LDFJAC, WA & LWA have been removed.
!
!INTEGER, INTENT(IN)        :: m
!INTEGER, INTENT(IN)        :: n
!REAL (dp), INTENT(IN OUT)  :: x(:)
!REAL (dp), INTENT(IN OUT)  :: fvec(:)
!REAL (dp), INTENT(IN OUT)  :: fjac(:,:)    ! fjac(ldfjac,n)
!REAL (dp), INTENT(IN)      :: tol
!INTEGER, INTENT(OUT)       :: info
!INTEGER, INTENT(IN OUT)    :: ipvt(:)
!
!
!! EXTERNAL fcn
!
!INTERFACE
!  SUBROUTINE fcn(m, n, x, fvec, fjac, iflag)
!    IMPLICIT NONE
!!    INTEGER, PARAMETER       :: dp = SELECTED_REAL_KIND(12, 60)
!    INTEGER, PARAMETER       :: dp = 8
!    INTEGER, INTENT(IN)      :: m, n
!    REAL (dp), INTENT(IN)    :: x(:)
!    REAL (dp), INTENT(OUT)   :: fvec(:)
!    REAL (dp), INTENT(OUT)   :: fjac(:,:)
!    INTEGER, INTENT(IN OUT)  :: iflag
!  END SUBROUTINE fcn
!END INTERFACE
!
!!     **********
!
!!     subroutine lmder1
!
!!     The purpose of lmder1 is to minimize the sum of the squares of
!!     m nonlinear functions in n variables by a modification of the
!!     levenberg-marquardt algorithm.  This is done by using the more
!!     general least-squares solver lmder.  The user must provide a
!!     subroutine which calculates the functions and the jacobian.
!
!!     the subroutine statement is
!
!!       subroutine lmder1(fcn, m, n, x, fvec, fjac, tol, info, ipvt)
!
!!     where
!
!!       fcn is the name of the user-supplied subroutine which
!!         calculates the functions and the jacobian.  fcn must
!!         be declared in an interface statement in the user
!!         calling program, and should be written as follows.
!
!!         subroutine fcn(m, n, x, fvec, fjac, iflag)
!!         integer   :: m, n, ldfjac, iflag
!!         REAL (dp) :: x(:), fvec(:), fjac(:,:)
!!         ----------
!!         if iflag = 1 calculate the functions at x and
!!         return this vector in fvec. do not alter fjac.
!!         if iflag = 2 calculate the jacobian at x and
!!         return this matrix in fjac. do not alter fvec.
!!         ----------
!!         return
!!         end
!
!!         the value of iflag should not be changed by fcn unless
!!         the user wants to terminate execution of lmder1.
!!         in this case set iflag to a negative integer.
!
!!       m is a positive integer input variable set to the number of functions.
!
!!       n is a positive integer input variable set to the number
!!         of variables.  n must not exceed m.
!
!!       x is an array of length n. on input x must contain
!!         an initial estimate of the solution vector. on output x
!!         contains the final estimate of the solution vector.
!
!!       fvec is an output array of length m which contains
!!         the functions evaluated at the output x.
!
!!       fjac is an output m by n array. the upper n by n submatrix
!!         of fjac contains an upper triangular matrix r with
!!         diagonal elements of nonincreasing magnitude such that
!
!!                t     t           t
!!               p *(jac *jac)*p = r *r,
!
!!         where p is a permutation matrix and jac is the final calculated
!!         Jacobian.  Column j of p is column ipvt(j) (see below) of the
!!         identity matrix.  The lower trapezoidal part of fjac contains
!!         information generated during the computation of r.
!
!!       ldfjac is a positive integer input variable not less than m
!!         which specifies the leading dimension of the array fjac.
!
!!       tol is a nonnegative input variable. termination occurs
!!         when the algorithm estimates either that the relative
!!         error in the sum of squares is at most tol or that
!!         the relative error between x and the solution is at most tol.
!
!!       info is an integer output variable.  If the user has terminated
!!         execution, info is set to the (negative) value of iflag.
!!         See description of fcn.  Otherwise, info is set as follows.
!
!!         info = 0  improper input parameters.
!
!!         info = 1  algorithm estimates that the relative error
!!                   in the sum of squares is at most tol.
!
!!         info = 2  algorithm estimates that the relative error
!!                   between x and the solution is at most tol.
!
!!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!!         info = 4  fvec is orthogonal to the columns of the
!!                   jacobian to machine precision.
!
!!         info = 5  number of calls to fcn with iflag = 1 has reached 100*(n+1).
!
!!         info = 6  tol is too small.  No further reduction in
!!                   the sum of squares is possible.
!
!!         info = 7  tol is too small.  No further improvement in
!!                   the approximate solution x is possible.
!
!!       ipvt is an integer output array of length n. ipvt
!!         defines a permutation matrix p such that jac*p = q*r,
!!         where jac is the final calculated jacobian, q is
!!         orthogonal (not stored), and r is upper triangular
!!         with diagonal elements of nonincreasing magnitude.
!!         column j of p is column ipvt(j) of the identity matrix.
!
!!       wa is a work array of length lwa.
!
!!       lwa is a positive integer input variable not less than 5*n+m.
!
!!     subprograms called
!
!!       user-supplied ...... fcn
!
!!       minpack-supplied ... lmder
!
!!     argonne national laboratory. minpack project. march 1980.
!!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!!     **********
!INTEGER   :: maxfev, mode, nfev, njev, nprint
!REAL (dp) :: ftol, gtol, xtol, wa(2*n)
!REAL (dp), PARAMETER :: factor = 100._dp, zero = 0.0_dp
!
!info = 0
!
!!     check the input parameters for errors.
!
!IF ( n <= 0 .OR. m < n .OR. tol < zero ) GO TO 10
!
!!     call lmder.
!
!maxfev = 100*(n + 1)
!ftol = tol
!xtol = tol
!gtol = zero
!mode = 1
!nprint = 0
!CALL lmder(fcn, m, n, x, fvec, fjac, ftol, xtol, gtol, maxfev,  &
!           wa, mode, factor, nprint, info, nfev, njev, ipvt, wa(n+1:) )
!IF (info == 8) info = 4
!
!10 RETURN
!
!!     last card of subroutine lmder1.
!
!END SUBROUTINE lmder1
!
!
!
!SUBROUTINE lmder(fcn, m, n, x, fvec, fjac, ftol, xtol, gtol, maxfev, &
!                 diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf)
! 
!! Code converted using TO_F90 by Alan Miller
!! Date: 1999-12-09  Time: 12:45:50
!
!! N.B. Arguments LDFJAC, WA1, WA2, WA3 & WA4 have been removed.
!
!INTEGER, INTENT(IN)        :: m
!INTEGER, INTENT(IN)        :: n
!REAL (dp), INTENT(IN OUT)  :: x(:)
!REAL (dp), INTENT(OUT)     :: fvec(m)
!REAL (dp), INTENT(OUT)     :: fjac(:,:)    ! fjac(ldfjac,n)
!REAL (dp), INTENT(IN)      :: ftol
!REAL (dp), INTENT(IN)      :: xtol
!REAL (dp), INTENT(IN OUT)  :: gtol
!INTEGER, INTENT(IN OUT)    :: maxfev
!REAL (dp), INTENT(OUT)     :: diag(:)
!INTEGER, INTENT(IN)        :: mode
!REAL (dp), INTENT(IN)      :: factor
!INTEGER, INTENT(IN)        :: nprint
!INTEGER, INTENT(OUT)       :: info
!INTEGER, INTENT(OUT)       :: nfev
!INTEGER, INTENT(OUT)       :: njev
!INTEGER, INTENT(OUT)       :: ipvt(:)
!REAL (dp), INTENT(OUT)     :: qtf(:)
!
!INTERFACE
!  SUBROUTINE fcn(m, n, x, fvec, fjac, iflag)
!    IMPLICIT NONE
!!    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
!    INTEGER, PARAMETER :: dp = 8
!    INTEGER, INTENT(IN)     :: m, n
!    REAL (dp), INTENT(IN)   :: x(:)
!    REAL (dp), INTENT(OUT)  :: fvec(:)
!    REAL (dp), INTENT(OUT)  :: fjac(:,:)
!    INTEGER, INTENT(IN)     :: iflag
!  END SUBROUTINE fcn
!END INTERFACE
!
!
!!     **********
!
!!     subroutine lmder
!
!!     the purpose of lmder is to minimize the sum of the squares of
!!     m nonlinear functions in n variables by a modification of
!!     the levenberg-marquardt algorithm. the user must provide a
!!     subroutine which calculates the functions and the jacobian.
!
!!     the subroutine statement is
!
!!       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
!!                        maxfev,diag,mode,factor,nprint,info,nfev,
!!                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
!
!!     where
!
!!       fcn is the name of the user-supplied subroutine which
!!         calculates the functions and the jacobian. fcn must
!!         be declared in an external statement in the user
!!         calling program, and should be written as follows.
!
!!         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
!!         integer m,n,ldfjac,iflag
!!         REAL (dp) x(:),fvec(m),fjac(ldfjac,n)
!!         ----------
!!         if iflag = 1 calculate the functions at x and
!!         return this vector in fvec. do not alter fjac.
!!         if iflag = 2 calculate the jacobian at x and
!!         return this matrix in fjac. do not alter fvec.
!!         ----------
!!         return
!!         end
!
!!         the value of iflag should not be changed by fcn unless
!!         the user wants to terminate execution of lmder.
!!         in this case set iflag to a negative integer.
!
!!       m is a positive integer input variable set to the number
!!         of functions.
!
!!       n is a positive integer input variable set to the number
!!         of variables. n must not exceed m.
!
!!       x is an array of length n. on input x must contain
!!         an initial estimate of the solution vector. on output x
!!         contains the final estimate of the solution vector.
!
!!       fvec is an output array of length m which contains
!!         the functions evaluated at the output x.
!
!!       fjac is an output m by n array. the upper n by n submatrix
!!         of fjac contains an upper triangular matrix r with
!!         diagonal elements of nonincreasing magnitude such that
!
!!                t     t           t
!!               p *(jac *jac)*p = r *r
!
!!         where p is a permutation matrix and jac is the final calculated
!!         jacobian.  Column j of p is column ipvt(j) (see below) of the
!!         identity matrix.  The lower trapezoidal part of fjac contains
!!         information generated during the computation of r.
!
!!       ldfjac is a positive integer input variable not less than m
!!         which specifies the leading dimension of the array fjac.
!
!!       ftol is a nonnegative input variable.  Termination occurs when both
!!         the actual and predicted relative reductions in the sum of squares
!!         are at most ftol.   Therefore, ftol measures the relative error
!!         desired in the sum of squares.
!
!!       xtol is a nonnegative input variable. termination
!!         occurs when the relative error between two consecutive
!!         iterates is at most xtol. therefore, xtol measures the
!!         relative error desired in the approximate solution.
!
!!       gtol is a nonnegative input variable.  Termination occurs when the
!!         cosine of the angle between fvec and any column of the jacobian is
!!         at most gtol in absolute value.  Therefore, gtol measures the
!!         orthogonality desired between the function vector and the columns
!!         of the jacobian.
!
!!       maxfev is a positive integer input variable.  Termination occurs when
!!         the number of calls to fcn with iflag = 1 has reached maxfev.
!
!!       diag is an array of length n.  If mode = 1 (see below), diag is
!!         internally set.  If mode = 2, diag must contain positive entries
!!         that serve as multiplicative scale factors for the variables.
!
!!       mode is an integer input variable.  if mode = 1, the
!!         variables will be scaled internally.  if mode = 2,
!!         the scaling is specified by the input diag.  other
!!         values of mode are equivalent to mode = 1.
!
!!       factor is a positive input variable used in determining the
!!         initial step bound. this bound is set to the product of
!!         factor and the euclidean norm of diag*x if nonzero, or else
!!         to factor itself. in most cases factor should lie in the
!!         interval (.1,100.).100. is a generally recommended value.
!
!!       nprint is an integer input variable that enables controlled printing
!!         of iterates if it is positive.  In this case, fcn is called with
!!         iflag = 0 at the beginning of the first iteration and every nprint
!!         iterations thereafter and immediately prior to return, with x, fvec,
!!         and fjac available for printing.  fvec and fjac should not be
!!         altered.  If nprint is not positive, no special calls of fcn with
!!         iflag = 0 are made.
!
!!       info is an integer output variable.  If the user has terminated
!!         execution, info is set to the (negative) value of iflag.
!!         See description of fcn.  Otherwise, info is set as follows.
!
!!         info = 0  improper input parameters.
!
!!         info = 1  both actual and predicted relative reductions
!!                   in the sum of squares are at most ftol.
!
!!         info = 2  relative error between two consecutive iterates
!!                   is at most xtol.
!
!!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!!         info = 4  the cosine of the angle between fvec and any column of
!!                   the jacobian is at most gtol in absolute value.
!
!!         info = 5  number of calls to fcn with iflag = 1 has reached maxfev.
!
!!         info = 6  ftol is too small.  No further reduction in
!!                   the sum of squares is possible.
!
!!         info = 7  xtol is too small.  No further improvement in
!!                   the approximate solution x is possible.
!
!!         info = 8  gtol is too small.  fvec is orthogonal to the
!!                   columns of the jacobian to machine precision.
!
!!       nfev is an integer output variable set to the number of
!!         calls to fcn with iflag = 1.
!
!!       njev is an integer output variable set to the number of
!!         calls to fcn with iflag = 2.
!
!!       ipvt is an integer output array of length n.  ipvt
!!         defines a permutation matrix p such that jac*p = q*r,
!!         where jac is the final calculated jacobian, q is
!!         orthogonal (not stored), and r is upper triangular
!!         with diagonal elements of nonincreasing magnitude.
!!         column j of p is column ipvt(j) of the identity matrix.
!
!!       qtf is an output array of length n which contains
!!         the first n elements of the vector (q transpose)*fvec.
!
!!       wa1, wa2, and wa3 are work arrays of length n.
!
!!       wa4 is a work array of length m.
!
!!     subprograms called
!
!!       user-supplied ...... fcn
!
!!       minpack-supplied ... dpmpar,enorm,lmpar,qrfac
!
!!       fortran-supplied ... ABS,MAX,MIN,SQRT,mod
!
!!     argonne national laboratory. minpack project. march 1980.
!!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!!     **********
!INTEGER   :: i, iflag, iter, j, l
!REAL (dp) :: actred, delta, dirder, epsmch, fnorm, fnorm1, gnorm,  &
!             par, pnorm, prered, ratio, sum, temp, temp1, temp2, xnorm
!REAL (dp) :: wa1(n), wa2(n), wa3(n), wa4(m)
!REAL (dp), PARAMETER :: one = 1.0_dp, p1 = 0.1_dp, p5 = 0.5_dp,  &
!                        p25 = 0.25_dp, p75 = 0.75_dp, p0001 = 0.0001_dp, &
!                        zero = 0.0_dp
!
!!     epsmch is the machine precision.
!
!epsmch = EPSILON(zero)
!
!info = 0
!iflag = 0
!nfev = 0
!njev = 0
!
!!     check the input parameters for errors.
!
!IF (n <= 0 .OR. m < n .OR. ftol < zero .OR. xtol < zero .OR. gtol < zero  &
!    .OR. maxfev <= 0 .OR. factor <= zero) GO TO 300
!IF (mode /= 2) GO TO 20
!DO  j = 1, n
!  IF (diag(j) <= zero) GO TO 300
!END DO
!
!!     evaluate the function at the starting point and calculate its norm.
!
!20 iflag = 1
!CALL fcn(m, n, x, fvec, fjac, iflag)
!nfev = 1
!IF (iflag < 0) GO TO 300
!fnorm = enorm(m, fvec)
!
!!     initialize levenberg-marquardt parameter and iteration counter.
!
!par = zero
!iter = 1
!
!!     beginning of the outer loop.
!
!!        calculate the jacobian matrix.
!
!30 iflag = 2
!CALL fcn(m, n, x, fvec, fjac, iflag)
!njev = njev + 1
!IF (iflag < 0) GO TO 300
!
!!        if requested, call fcn to enable printing of iterates.
!
!IF (nprint <= 0) GO TO 40
!iflag = 0
!IF (MOD(iter-1,nprint) == 0) CALL fcn(m, n, x, fvec, fjac, iflag)
!IF (iflag < 0) GO TO 300
!
!!        compute the qr factorization of the jacobian.
!
!40 CALL qrfac(m, n, fjac, .true., ipvt, wa1, wa2)
!
!!        on the first iteration and if mode is 1, scale according
!!        to the norms of the columns of the initial jacobian.
!
!IF (iter /= 1) GO TO 80
!IF (mode == 2) GO TO 60
!DO  j = 1, n
!  diag(j) = wa2(j)
!  IF (wa2(j) == zero) diag(j) = one
!END DO
!
!!        on the first iteration, calculate the norm of the scaled x
!!        and initialize the step bound delta.
!
!60 wa3(1:n) = diag(1:n)*x(1:n)
!xnorm = enorm(n,wa3)
!delta = factor*xnorm
!IF (delta == zero) delta = factor
!
!!        form (q transpose)*fvec and store the first n components in qtf.
!
!80 wa4(1:m) = fvec(1:m)
!DO  j = 1, n
!  IF (fjac(j,j) == zero) GO TO 120
!  sum = DOT_PRODUCT( fjac(j:m,j), wa4(j:m) )
!  temp = -sum/fjac(j,j)
!  DO  i = j, m
!    wa4(i) = wa4(i) + fjac(i,j)*temp
!  END DO
!  120 fjac(j,j) = wa1(j)
!  qtf(j) = wa4(j)
!END DO
!
!!        compute the norm of the scaled gradient.
!
!gnorm = zero
!IF (fnorm == zero) GO TO 170
!DO  j = 1, n
!  l = ipvt(j)
!  IF (wa2(l) == zero) CYCLE
!  sum = zero
!  DO  i = 1, j
!    sum = sum + fjac(i,j)*(qtf(i)/fnorm)
!  END DO
!  gnorm = MAX(gnorm,ABS(sum/wa2(l)))
!END DO
!
!!        test for convergence of the gradient norm.
!
!170 IF (gnorm <= gtol) info = 4
!IF (info /= 0) GO TO 300
!
!!        rescale if necessary.
!
!IF (mode == 2) GO TO 200
!DO  j = 1, n
!  diag(j) = MAX(diag(j), wa2(j))
!END DO
!
!!        beginning of the inner loop.
!
!!           determine the levenberg-marquardt parameter.
!
!200 CALL lmpar(n, fjac, ipvt, diag, qtf, delta, par, wa1, wa2)
!
!!           store the direction p and x + p. calculate the norm of p.
!
!DO  j = 1, n
!  wa1(j) = -wa1(j)
!  wa2(j) = x(j) + wa1(j)
!  wa3(j) = diag(j)*wa1(j)
!END DO
!pnorm = enorm(n, wa3)
!
!!           on the first iteration, adjust the initial step bound.
!
!IF (iter == 1) delta = MIN(delta,pnorm)
!
!!           evaluate the function at x + p and calculate its norm.
!
!iflag = 1
!CALL fcn(m, n, wa2, wa4, fjac, iflag)
!nfev = nfev + 1
!IF (iflag < 0) GO TO 300
!fnorm1 = enorm(m, wa4)
!
!!           compute the scaled actual reduction.
!
!actred = -one
!IF (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2
!
!!           compute the scaled predicted reduction and
!!           the scaled directional derivative.
!
!DO  j = 1, n
!  wa3(j) = zero
!  l = ipvt(j)
!  temp = wa1(l)
!  DO  i = 1, j
!    wa3(i) = wa3(i) + fjac(i,j)*temp
!  END DO
!END DO
!temp1 = enorm(n,wa3)/fnorm
!temp2 = (SQRT(par)*pnorm)/fnorm
!prered = temp1**2 + temp2**2/p5
!dirder = -(temp1**2 + temp2**2)
!
!!           compute the ratio of the actual to the predicted reduction.
!
!ratio = zero
!IF (prered /= zero) ratio = actred/prered
!
!!           update the step bound.
!
!IF (ratio <= p25) THEN
!  IF (actred >= zero) temp = p5
!  IF (actred < zero) temp = p5*dirder/(dirder + p5*actred)
!  IF (p1*fnorm1 >= fnorm .OR. temp < p1) temp = p1
!  delta = temp*MIN(delta, pnorm/p1)
!  par = par/temp
!ELSE
!  IF (par /= zero .AND. ratio < p75) GO TO 260
!  delta = pnorm/p5
!  par = p5*par
!END IF
!
!!           test for successful iteration.
!
!260 IF (ratio < p0001) GO TO 290
!
!!           successful iteration. update x, fvec, and their norms.
!
!DO  j = 1, n
!  x(j) = wa2(j)
!  wa2(j) = diag(j)*x(j)
!END DO
!fvec(1:m) = wa4(1:m)
!xnorm = enorm(n,wa2)
!fnorm = fnorm1
!iter = iter + 1
!
!!           tests for convergence.
!
!290 IF (ABS(actred) <= ftol .AND. prered <= ftol .AND. p5*ratio <= one) info = 1
!IF (delta <= xtol*xnorm) info = 2
!IF (ABS(actred) <= ftol .AND. prered <= ftol  &
!    .AND. p5*ratio <= one .AND. info == 2) info = 3
!IF (info /= 0) GO TO 300
!
!!           tests for termination and stringent tolerances.
!
!IF (nfev >= maxfev) info = 5
!IF (ABS(actred) <= epsmch .AND. prered <= epsmch  &
!    .AND. p5*ratio <= one) info = 6
!IF (delta <= epsmch*xnorm) info = 7
!IF (gnorm <= epsmch) info = 8
!IF (info /= 0) GO TO 300
!
!!           end of the inner loop. repeat if iteration unsuccessful.
!
!IF (ratio < p0001) GO TO 200
!
!!        end of the outer loop.
!
!GO TO 30
!
!!     termination, either normal or user imposed.
!
!300 IF (iflag < 0) info = iflag
!iflag = 0
!IF (nprint > 0) CALL fcn(m, n, x, fvec, fjac, iflag)
!RETURN
!
!!     last card of subroutine lmder.
!
!END SUBROUTINE lmder


SUBROUTINE lmpar(n, r, ipvt, diag, qtb, delta, par, x, sdiag)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 12:46:12

! N.B. Arguments LDR, WA1 & WA2 have been removed.

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: r(:,:)
INTEGER, INTENT(IN)        :: ipvt(:)
REAL (dp), INTENT(IN)      :: diag(:)
REAL (dp), INTENT(IN)      :: qtb(:)
REAL (dp), INTENT(IN)      :: delta
REAL (dp), INTENT(OUT)     :: par
REAL (dp), INTENT(OUT)     :: x(:)
REAL (dp), INTENT(OUT)     :: sdiag(:)


!     **********

!     subroutine lmpar

!     given an m by n matrix a, an n by n nonsingular diagonal
!     matrix d, an m-vector b, and a positive number delta,
!     the problem is to determine a value for the parameter
!     par such that if x solves the system

!           a*x = b ,     sqrt(par)*d*x = 0 ,

!     in the least squares sense, and dxnorm is the euclidean
!     norm of d*x, then either par is zero and

!           (dxnorm-delta) <= 0.1*delta ,

!     or par is positive and

!           abs(dxnorm-delta) <= 0.1*delta .

!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization, with column pivoting, of a. that is, if
!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!     columns, and r is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then lmpar expects
!     the full upper triangle of r, the permutation matrix p,
!     and the first n components of (q transpose)*b. on output
!     lmpar also provides an upper triangular matrix s such that

!            t   t                   t
!           p *(a *a + par*d*d)*p = s *s .

!     s is employed within lmpar and may be of separate interest.

!     only a few iterations are generally needed for convergence
!     of the algorithm. if, however, the limit of 10 iterations
!     is reached, then the output par will contain the best
!     value obtained so far.

!     the subroutine statement is

!       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag, wa1,wa2)

!     where

!       n is a positive integer input variable set to the order of r.

!       r is an n by n array. on input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         on output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix s.

!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.

!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that a*p = q*r. column j of p
!         is column ipvt(j) of the identity matrix.

!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.

!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.

!       delta is a positive input variable which specifies an upper
!         bound on the euclidean norm of d*x.

!       par is a nonnegative variable. on input par contains an
!         initial estimate of the levenberg-marquardt parameter.
!         on output par contains the final estimate.

!       x is an output array of length n which contains the least
!         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
!         for the output par.

!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix s.

!       wa1 and wa2 are work arrays of length n.

!     subprograms called

!       minpack-supplied ... dpmpar,enorm,qrsolv

!       fortran-supplied ... ABS,MAX,MIN,SQRT

!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more

!     **********
INTEGER   :: i, iter, j, jm1, jp1, k, l, nsing
REAL (dp) :: dxnorm, dwarf, fp, gnorm, parc, parl, paru, sum, temp
REAL (dp) :: wa1(n), wa2(n)
REAL (dp), PARAMETER :: p1 = 0.1_dp, p001 = 0.001_dp, zero = 0.0_dp

!     dwarf is the smallest positive magnitude.

dwarf = TINY(zero)

!     compute and store in x the gauss-newton direction. if the
!     jacobian is rank-deficient, obtain a least squares solution.

nsing = n
DO  j = 1, n
  wa1(j) = qtb(j)
  IF (r(j,j) == zero .AND. nsing == n) nsing = j - 1
  IF (nsing < n) wa1(j) = zero
END DO

DO  k = 1, nsing
  j = nsing - k + 1
  wa1(j) = wa1(j)/r(j,j)
  temp = wa1(j)
  jm1 = j - 1
  DO  i = 1, jm1
    wa1(i) = wa1(i) - r(i,j)*temp
  END DO
END DO

DO  j = 1, n
  l = ipvt(j)
  x(l) = wa1(j)
END DO

!     initialize the iteration counter.
!     evaluate the function at the origin, and test
!     for acceptance of the gauss-newton direction.

iter = 0
DO  j = 1, n
  wa2(j) = diag(j)*x(j)
END DO
dxnorm = enorm(n, wa2)
fp = dxnorm - delta
IF (fp <= p1*delta) GO TO 220

!     if the jacobian is not rank deficient, the newton
!     step provides a lower bound, parl, for the zero of
!     the function.  Otherwise set this bound to zero.

parl = zero
IF (nsing < n) GO TO 120
DO  j = 1, n
  l = ipvt(j)
  wa1(j) = diag(l)*(wa2(l)/dxnorm)
END DO
DO  j = 1, n
  sum = zero
  jm1 = j - 1
  DO  i = 1, jm1
    sum = sum + r(i,j)*wa1(i)
  END DO
  wa1(j) = (wa1(j) - sum)/r(j,j)
END DO
temp = enorm(n,wa1)
parl = ((fp/delta)/temp)/temp

!     calculate an upper bound, paru, for the zero of the function.

120 DO  j = 1, n
  sum = zero
  DO  i = 1, j
    sum = sum + r(i,j)*qtb(i)
  END DO
  l = ipvt(j)
  wa1(j) = sum/diag(l)
END DO
gnorm = enorm(n,wa1)
paru = gnorm/delta
IF (paru == zero) paru = dwarf/MIN(delta,p1)

!     if the input par lies outside of the interval (parl,paru),
!     set par to the closer endpoint.

par = MAX(par,parl)
par = MIN(par,paru)
IF (par == zero) par = gnorm/dxnorm

!     beginning of an iteration.

150 iter = iter + 1

!        evaluate the function at the current value of par.

IF (par == zero) par = MAX(dwarf, p001*paru)
temp = SQRT(par)
wa1(1:n) = temp*diag(1:n)
CALL qrsolv(n, r, ipvt, wa1, qtb, x, sdiag)
wa2(1:n) = diag(1:n)*x(1:n)
dxnorm = enorm(n, wa2)
temp = fp
fp = dxnorm - delta

!        if the function is small enough, accept the current value
!        of par. also test for the exceptional cases where parl
!        is zero or the number of iterations has reached 10.

IF (ABS(fp) <= p1*delta .OR. parl == zero .AND. fp <= temp  &
    .AND. temp < zero .OR. iter == 10) GO TO 220

!        compute the newton correction.

DO  j = 1, n
  l = ipvt(j)
  wa1(j) = diag(l)*(wa2(l)/dxnorm)
END DO
DO  j = 1, n
  wa1(j) = wa1(j)/sdiag(j)
  temp = wa1(j)
  jp1 = j + 1
  DO  i = jp1, n
    wa1(i) = wa1(i) - r(i,j)*temp
  END DO
END DO
temp = enorm(n,wa1)
parc = ((fp/delta)/temp)/temp

!        depending on the sign of the function, update parl or paru.

IF (fp > zero) parl = MAX(parl,par)
IF (fp < zero) paru = MIN(paru,par)

!        compute an improved estimate for par.

par = MAX(parl, par+parc)

!        end of an iteration.

GO TO 150

!     termination.

220 IF (iter == 0) par = zero
RETURN

!     last card of subroutine lmpar.

END SUBROUTINE lmpar



SUBROUTINE qrfac(m, n, a, pivot, ipvt, rdiag, acnorm)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 12:46:17

! N.B. Arguments LDA, LIPVT & WA has been removed.

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: a(:,:)
LOGICAL, INTENT(IN)        :: pivot
INTEGER, INTENT(OUT)       :: ipvt(:)
REAL (dp), INTENT(OUT)     :: rdiag(:)
REAL (dp), INTENT(OUT)     :: acnorm(:)


!     **********

!     subroutine qrfac

!     this subroutine uses householder transformations with column
!     pivoting (optional) to compute a qr factorization of the
!     m by n matrix a. that is, qrfac determines an orthogonal
!     matrix q, a permutation matrix p, and an upper trapezoidal
!     matrix r with diagonal elements of nonincreasing magnitude,
!     such that a*p = q*r. the householder transformation for
!     column k, k = 1,2,...,min(m,n), is of the form

!                           t
!           i - (1/u(k))*u*u

!     where u has zeros in the first k-1 positions. the form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding linpack subroutine.

!     the subroutine statement is

!       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)

!     where

!       m is a positive integer input variable set to the number of rows of a.

!       n is a positive integer input variable set to the number
!         of columns of a.

!       a is an m by n array. on input a contains the matrix for
!         which the qr factorization is to be computed.  on output
!         the strict upper trapezoidal part of a contains the strict
!         upper trapezoidal part of r, and the lower trapezoidal
!         part of a contains a factored form of q (the non-trivial
!         elements of the u vectors described above).

!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.

!       pivot is a logical input variable.  if pivot is set true,
!         then column pivoting is enforced.  if pivot is set false,
!         then no column pivoting is done.

!       ipvt is an integer output array of length lipvt.  ipvt
!         defines the permutation matrix p such that a*p = q*r.
!         column j of p is column ipvt(j) of the identity matrix.
!         if pivot is false, ipvt is not referenced.

!       lipvt is a positive integer input variable. if pivot is false,
!         then lipvt may be as small as 1. if pivot is true, then
!         lipvt must be at least n.

!       rdiag is an output array of length n which contains the
!         diagonal elements of r.

!       acnorm is an output array of length n which contains the
!         norms of the corresponding columns of the input matrix a.
!         If this information is not needed, then acnorm can coincide
!         with rdiag.

!       wa is a work array of length n. if pivot is false, then wa
!         can coincide with rdiag.

!     subprograms called

!       minpack-supplied ... dpmpar,enorm

!       fortran-supplied ... MAX,SQRT,MIN

!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more

!     **********
INTEGER   :: i, j, jp1, k, kmax, minmn
REAL (dp) :: ajnorm, epsmch, sum, temp, wa(n)
REAL (dp), PARAMETER :: one = 1.0_dp, p05 = 0.05_dp, zero = 0.0_dp

!     epsmch is the machine precision.

epsmch = EPSILON(zero)

!     compute the initial column norms and initialize several arrays.

DO  j = 1, n
  acnorm(j) = enorm(m,a(1:,j))
  rdiag(j) = acnorm(j)
  wa(j) = rdiag(j)
  IF (pivot) ipvt(j) = j
END DO

!     reduce a to r with householder transformations.

minmn = MIN(m,n)
DO  j = 1, minmn
  IF (.NOT.pivot) GO TO 40
  
!        bring the column of largest norm into the pivot position.
  
  kmax = j
  DO  k = j, n
    IF (rdiag(k) > rdiag(kmax)) kmax = k
  END DO
  IF (kmax == j) GO TO 40
  DO  i = 1, m
    temp = a(i,j)
    a(i,j) = a(i,kmax)
    a(i,kmax) = temp
  END DO
  rdiag(kmax) = rdiag(j)
  wa(kmax) = wa(j)
  k = ipvt(j)
  ipvt(j) = ipvt(kmax)
  ipvt(kmax) = k
  
!        compute the householder transformation to reduce the
!        j-th column of a to a multiple of the j-th unit vector.
  
  40 ajnorm = enorm(m-j+1, a(j:,j))
  IF (ajnorm == zero) CYCLE
  IF (a(j,j) < zero) ajnorm = -ajnorm
  DO  i = j, m
    a(i,j) = a(i,j)/ajnorm
  END DO
  a(j,j) = a(j,j) + one
  
!        apply the transformation to the remaining columns
!        and update the norms.
  
  jp1 = j + 1
  DO  k = jp1, n
    sum = zero
    DO  i = j, m
      sum = sum + a(i,j)*a(i,k)
    END DO
    temp = sum/a(j,j)
    DO  i = j, m
      a(i,k) = a(i,k) - temp*a(i,j)
    END DO
    IF (.NOT.pivot .OR. rdiag(k) == zero) CYCLE
    temp = a(j,k)/rdiag(k)
    rdiag(k) = rdiag(k)*SQRT(MAX(zero, one-temp**2))
    IF (p05*(rdiag(k)/wa(k))**2 > epsmch) CYCLE
    rdiag(k) = enorm(m-j, a(jp1:,k))
    wa(k) = rdiag(k)
  END DO
  rdiag(j) = -ajnorm
END DO
RETURN

!     last card of subroutine qrfac.

END SUBROUTINE qrfac



SUBROUTINE qrsolv(n, r, ipvt, diag, qtb, x, sdiag)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 12:46:21

! N.B. Arguments LDR & WA has been removed.

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: r(:,:)
INTEGER, INTENT(IN)        :: ipvt(:)
REAL (dp), INTENT(IN)      :: diag(:)
REAL (dp), INTENT(IN)      :: qtb(:)
REAL (dp), INTENT(OUT)     :: x(:)
REAL (dp), INTENT(OUT)     :: sdiag(:)


!     **********

!     subroutine qrsolv

!     given an m by n matrix a, an n by n diagonal matrix d,
!     and an m-vector b, the problem is to determine an x which
!     solves the system

!           a*x = b ,     d*x = 0 ,

!     in the least squares sense.

!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization, with column pivoting, of a. that is, if
!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!     columns, and r is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then qrsolv expects
!     the full upper triangle of r, the permutation matrix p,
!     and the first n components of (q transpose)*b. the system
!     a*x = b, d*x = 0, is then equivalent to

!                  t       t
!           r*z = q *b ,  p *d*p*z = 0 ,

!     where x = p*z. if this system does not have full rank,
!     then a least squares solution is obtained. on output qrsolv
!     also provides an upper triangular matrix s such that

!            t   t               t
!           p *(a *a + d*d)*p = s *s .

!     s is computed within qrsolv and may be of separate interest.

!     the subroutine statement is

!       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)

!     where

!       n is a positive integer input variable set to the order of r.

!       r is an n by n array. on input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         on output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix s.

!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.

!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that a*p = q*r. column j of p
!         is column ipvt(j) of the identity matrix.

!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.

!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.

!       x is an output array of length n which contains the least
!         squares solution of the system a*x = b, d*x = 0.

!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix s.

!       wa is a work array of length n.

!     subprograms called

!       fortran-supplied ... ABS,SQRT

!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more

!     **********
INTEGER   :: i, j, jp1, k, kp1, l, nsing
REAL (dp) :: COS, cotan, qtbpj, SIN, sum, TAN, temp, wa(n)
REAL (dp), PARAMETER :: p5 = 0.5_dp, p25 = 0.25_dp, zero = 0.0_dp

!     copy r and (q transpose)*b to preserve input and initialize s.
!     in particular, save the diagonal elements of r in x.

DO  j = 1, n
  DO  i = j, n
    r(i,j) = r(j,i)
  END DO
  x(j) = r(j,j)
  wa(j) = qtb(j)
END DO

!     eliminate the diagonal matrix d using a givens rotation.

DO  j = 1, n
  
!        prepare the row of d to be eliminated, locating the
!        diagonal element using p from the qr factorization.
  
  l = ipvt(j)
  IF (diag(l) == zero) CYCLE
  sdiag(j:n) = zero
  sdiag(j) = diag(l)
  
!        the transformations to eliminate the row of d
!        modify only a single element of (q transpose)*b
!        beyond the first n, which is initially zero.
  
  qtbpj = zero
  DO  k = j, n
    
!           determine a givens rotation which eliminates the
!           appropriate element in the current row of d.
    
    IF (sdiag(k) == zero) CYCLE
    IF (ABS(r(k,k)) < ABS(sdiag(k))) THEN
      cotan = r(k,k)/sdiag(k)
      SIN = p5/SQRT(p25 + p25*cotan**2)
      COS = SIN*cotan
    ELSE
      TAN = sdiag(k)/r(k,k)
      COS = p5/SQRT(p25 + p25*TAN**2)
      SIN = COS*TAN
    END IF
    
!           compute the modified diagonal element of r and
!           the modified element of ((q transpose)*b,0).
    
    r(k,k) = COS*r(k,k) + SIN*sdiag(k)
    temp = COS*wa(k) + SIN*qtbpj
    qtbpj = -SIN*wa(k) + COS*qtbpj
    wa(k) = temp
    
!           accumulate the tranformation in the row of s.
    
    kp1 = k + 1
    DO  i = kp1, n
      temp = COS*r(i,k) + SIN*sdiag(i)
      sdiag(i) = -SIN*r(i,k) + COS*sdiag(i)
      r(i,k) = temp
    END DO
  END DO
  
!        store the diagonal element of s and restore
!        the corresponding diagonal element of r.
  
  sdiag(j) = r(j,j)
  r(j,j) = x(j)
END DO

!     solve the triangular system for z. if the system is
!     singular, then obtain a least squares solution.

nsing = n
DO  j = 1, n
  IF (sdiag(j) == zero .AND. nsing == n) nsing = j - 1
  IF (nsing < n) wa(j) = zero
END DO

DO  k = 1, nsing
  j = nsing - k + 1
  sum = zero
  jp1 = j + 1
  DO  i = jp1, nsing
    sum = sum + r(i,j)*wa(i)
  END DO
  wa(j) = (wa(j) - sum)/sdiag(j)
END DO

!     permute the components of z back to components of x.

DO  j = 1, n
  l = ipvt(j)
  x(l) = wa(j)
END DO
RETURN

!     last card of subroutine qrsolv.

END SUBROUTINE qrsolv



FUNCTION enorm(n,x) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 12:45:34

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: x(:)
REAL (dp)              :: fn_val


!     **********

!     function enorm

!     given an n-vector x, this function calculates the euclidean norm of x.

!     the euclidean norm is computed by accumulating the sum of squares in
!     three different sums.  The sums of squares for the small and large
!     components are scaled so that no overflows occur.  Non-destructive
!     underflows are permitted.  Underflows and overflows do not occur in the
!     computation of the unscaled sum of squares for the intermediate
!     components.  The definitions of small, intermediate and large components
!     depend on two constants, rdwarf and rgiant.  The main restrictions on
!     these constants are that rdwarf**2 not underflow and rgiant**2 not
!     overflow.  The constants given here are suitable for every known computer.

!     the function statement is

!       REAL (dp) function enorm(n,x)

!     where

!       n is a positive integer input variable.

!       x is an input array of length n.

!     subprograms called

!       fortran-supplied ... ABS,SQRT

!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more

!     **********
INTEGER   :: i
REAL (dp) :: agiant, floatn, s1, s2, s3, xabs, x1max, x3max
REAL (dp), PARAMETER :: one = 1.0_dp, zero = 0.0_dp, rdwarf = 3.834E-20_dp,  &
                        rgiant = 1.304E+19_dp

s1 = zero
s2 = zero
s3 = zero
x1max = zero
x3max = zero
floatn = n
agiant = rgiant/floatn
DO  i = 1, n
  xabs = ABS(x(i))
  IF (xabs > rdwarf .AND. xabs < agiant) GO TO 70
  IF (xabs <= rdwarf) GO TO 30
  
!              sum for large components.
  
  IF (xabs <= x1max) GO TO 10
  s1 = one + s1*(x1max/xabs)**2
  x1max = xabs
  GO TO 20

  10 s1 = s1 + (xabs/x1max)**2

  20 GO TO 60
  
!              sum for small components.
  
  30 IF (xabs <= x3max) GO TO 40
  s3 = one + s3*(x3max/xabs)**2
  x3max = xabs
  GO TO 60

  40 IF (xabs /= zero) s3 = s3 + (xabs/x3max)**2

  60 CYCLE
  
!           sum for intermediate components.
  
  70 s2 = s2 + xabs**2
END DO

!     calculation of norm.

IF (s1 == zero) GO TO 100
fn_val = x1max*SQRT(s1 + (s2/x1max)/x1max)
GO TO 120

100 IF (s2 == zero) GO TO 110
IF (s2 >= x3max) fn_val = SQRT(s2*(one + (x3max/s2)*(x3max*s3)))
IF (s2 < x3max) fn_val = SQRT(x3max*((s2/x3max) + (x3max*s3)))
GO TO 120

110 fn_val = x3max*SQRT(s3)

120 RETURN

!     last card of function enorm.

END FUNCTION enorm



SUBROUTINE fdjac2(fcn, m, n, x, fvec, fjac, iflag, epsfcn)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 12:45:44

! N.B. Arguments LDFJAC & WA have been removed.

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
REAL (dp), INTENT(IN)      :: fvec(m)
REAL (dp), INTENT(OUT)     :: fjac(:,:)    ! fjac(ldfjac,n)
INTEGER, INTENT(INOUT)        :: iflag
REAL (dp), INTENT(IN)      :: epsfcn

INTERFACE
  SUBROUTINE fcn(m, n, x, fvec, iflag)
    IMPLICIT NONE
!    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, PARAMETER :: dp = 8
    INTEGER, INTENT(IN)     :: m, n
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: fvec(:)
    INTEGER, INTENT(IN OUT)     :: iflag
  END SUBROUTINE fcn
END INTERFACE

!     **********

!     subroutine fdjac2

!     this subroutine computes a forward-difference approximation
!     to the m by n jacobian matrix associated with a specified
!     problem of m functions in n variables.

!     the subroutine statement is

!       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)

!     where

!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.

!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         REAL (dp) x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end

!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of fdjac2.
!         in this case set iflag to a negative integer.

!       m is a positive integer input variable set to the number of functions.

!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.

!       x is an input array of length n.

!       fvec is an input array of length m which must contain the
!         functions evaluated at x.

!       fjac is an output m by n array which contains the
!         approximation to the jacobian matrix evaluated at x.

!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.

!       iflag is an integer variable which can be used to terminate
!         the execution of fdjac2.  see description of fcn.

!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine precision.

!       wa is a work array of length m.

!     subprograms called

!       user-supplied ...... fcn

!       minpack-supplied ... dpmpar

!       fortran-supplied ... ABS,MAX,SQRT

!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more

!     **********
INTEGER   :: i, j
REAL (dp) :: eps, epsmch, h, temp, wa(m)
REAL (dp), PARAMETER :: zero = 0.0_dp

!     epsmch is the machine precision.

epsmch = EPSILON(zero)

eps = SQRT(MAX(epsfcn, epsmch))
DO  j = 1, n
  temp = x(j)
  h = eps*ABS(temp)
  IF (h == zero) h = eps
  x(j) = temp + h
  CALL fcn(m, n, x, wa, iflag)
  IF (iflag < 0) EXIT
  x(j) = temp
  DO  i = 1, m
    fjac(i,j) = (wa(i) - fvec(i))/h
  END DO
END DO

RETURN

!     last card of subroutine fdjac2.

END SUBROUTINE fdjac2


END MODULE M_Levenberg_Marquardt
subroutine S_lmdif_ODE_linear &
    ( Nequation , Ntime , timestep , TIMEPOINT , Nnonzero , NONZERO_ROW ,&
    NONZERO_COLUMN , Nobservation , OBSERVATION_EQUATION ,&
    OBSERVATION_TIME , OBSERVATION_VALUE , WEIGHT ,&
    tol , max_func , dif_step , Nprint , X , Y )
    use M_rk4_parameter ,&
        only : S_set_Nequation , S_set_Ntime , S_set_timestep ,&
            S_set_TIMEPOINT
    use M_derivative_linear_ODE ,&
        only : S_set_Nnonzero , S_set_NONZERO_ROW , S_set_NONZERO_COLUMN
    use M_least_squares ,&
        only : S_set_Nobservation , S_set_OBSERVATION_EQUATION ,&
            S_set_OBSERVATION_TIME , S_set_OBSERVATION_VALUE ,&
            S_set_WEIGHT , S_ODE_linear_ls_lmdif
    use M_Levenberg_Marquardt ,&
        only : lmdif
    implicit none
    ! Argument
    integer(kind=4) , intent(in) :: Nequation
    integer(kind=4) , intent(in) :: Ntime
    real(kind=8) , intent(in) :: timestep
    real(kind=8) , dimension(Ntime) , intent(inout) :: TIMEPOINT
    integer(kind=4) , intent(in) :: Nnonzero
    integer(kind=4) , dimension(Nnonzero) , intent(in) :: NONZERO_ROW
    integer(kind=4) , dimension(Nnonzero) , intent(in) :: NONZERO_COLUMN
    integer(kind=4) , intent(in) :: Nobservation
    integer(kind=4) , dimension(Nobservation) , intent(in) :: OBSERVATION_EQUATION
    real(kind=8) , dimension(Nobservation) , intent(in) :: OBSERVATION_TIME
    real(kind=8) , dimension(Nobservation) , intent(in) :: OBSERVATION_VALUE
    real(kind=8) , dimension(Nobservation) , intent(in) :: WEIGHT
    real(kind=8) , intent(inout) :: tol
    integer(kind=4) , intent(inout) :: max_func
    real(kind=8) , intent(inout) :: dif_step
    integer(kind=4) , intent(in) :: Nprint
    real(kind=8) , dimension(Nnonzero+2*Nequation) , intent(inout) :: X
    real(kind=8) , dimension(Nobservation) , intent(out) :: Y
    ! LMDIF argument
    integer(kind=4) :: m
    integer(kind=4) :: n
    real(kind=8) , dimension(:) , allocatable :: diag
    integer(kind=4) :: mode = 1
    real(kind=8) :: factor = 100.0
    integer(kind=4) :: info
    integer(kind=4) :: Nfunc
    real(kind=8) , dimension(:,:) , allocatable :: fjac
    integer(kind=4) , dimension(:) , allocatable :: ipvt
    real(kind=8) , dimension(:) , allocatable :: qtf
    ! F77 LMDIF argument
!    integer(kind=4) :: ldfjac
!    real(kind=8) , dimension(:) , allocatable :: w1
!    real(kind=8) , dimension(:) , allocatable :: w2
!    real(kind=8) , dimension(:) , allocatable :: w3
!    real(kind=8) , dimension(:) , allocatable :: w4
    ! Local
    integer(kind=4) :: i
    ! Initialization
    CALL S_set_Nequation ( Nequation )
    CALL S_set_Ntime ( Ntime )
    CALL S_set_timestep ( timestep )
    CALL S_set_TIMEPOINT ( TIMEPOINT )
    CALL S_set_Nnonzero ( Nnonzero )
    call S_set_NONZERO_ROW ( NONZERO_ROW )
    call S_set_NONZERO_COLUMN ( NONZERO_COLUMN )
    CALL S_set_Nobservation ( Nobservation )
    call S_set_OBSERVATION_EQUATION ( OBSERVATION_EQUATION )
    call S_set_OBSERVATION_TIME ( OBSERVATION_TIME )
    call S_set_OBSERVATION_VALUE ( OBSERVATION_VALUE )
    call S_set_WEIGHT ( WEIGHT )
    ! LMDIF array initialization
    m = size(Y)
!    write(*,*) 'm=',m
    n = size(X)
    allocate ( diag (n) , stat = i )
    if ( i .ne. 0 ) then
        write(*,*) 'Unable to allocate memory for LMDIF-initialization.'
        stop
    end if
    diag = 1.0
    allocate ( fjac (m,n) , stat = i )
    if ( i .ne. 0 ) then
        write(*,*) 'Unable to allocate memory for LMDIF-initialization.'
        stop
    end if
    allocate ( ipvt (n) , stat = i )
    if ( i .ne. 0 ) then
        write(*,*) 'Unable to allocate memory for LMDIF-initialization.'
        stop
    end if
    allocate ( qtf (n) , stat = i )
    if ( i .ne. 0 ) then
        write(*,*) 'Unable to allocate memory for LMDIF-initialization.'
        stop
    end if
    write(*,*) 'Initialization Completed.  Now iteration starts.'
    call lmdif ( S_ODE_linear_ls_lmdif , m , n , X , Y , tol , tol , tol, &
        max_func , dif_step , diag , mode , factor , Nprint , info , Nfunc ,&
        fjac , ipvt , qtf )
    ! F77 LMDIF
!    ldfjac = m
!    allocate ( w1 (n) , stat = i )
!    if ( i .ne. 0 ) then
!        write(*,*) 'Unable to allocate memory for LMDIF-initialization.'
!        stop
!    end if
!    allocate ( w2 (n) , stat = i )
!    if ( i .ne. 0 ) then
!        write(*,*) 'Unable to allocate memory for LMDIF-initialization.'
!        stop
!    end if
!    allocate ( w3 (n) , stat = i )
!    if ( i .ne. 0 ) then
!        write(*,*) 'Unable to allocate memory for LMDIF-initialization.'
!        stop
!    end if
!    allocate ( w4 (m) , stat = i )
!    if ( i .ne. 0 ) then
!        write(*,*) 'Unable to allocate memory for LMDIF-initialization.'
!        stop
!    end if
!    call lmdif ( S_ODE_linear_ls_lmdif , m , n , X , Y , tol , tol , tol, &
!        max_func , dif_step , diag , mode , factor , Nprint , info , Nfunc ,&
!        fjac , ldfjac , ipvt , qtf , w1 , w2 , w3 , w4 )
    ! Output
    write(*,*) 'Iteration terminates because:'
    select case ( info )
        case ( 0 )
            write(*,*) 'Improper input parameters.'
        case ( 1 )
            write(*,*) 'Both actual and predicted relative reductions &
                in the object function value are smaller than &
                given tolerance.'
        case ( 2 )
            write(*,*) 'Relative error between two consecutive iterates &
                is smaller than given tolerance.'
        case ( 3 )
            write(*,*) 'Actual and predicted relative reductions &
                in the object function value , &
                and relative error between two consecutive iterates &
                are all smaller than given tolerance.'
        case ( 4 )
            write(*,*) 'The cosine of the angle between object function &
                and any column of the Jacobian matrix is smaller than &
                given tolerance in absolute value.'
        case ( 5 )
            write(*,*) 'Maximal number of object function evaluation &
                is reached.'
        case ( 6:8 )
            write(*,*) 'Given tolerance may be too small.  &
                No further reduction is possible.'
        case default
            write(*,*) 'Unknown mistake.'
    end select
    write(*,*) 'Total number of object function evaluation:' , Nfunc
    write(*,*) 'Final residual sum of squares:' , norm2(Y)**2
    return
end subroutine S_lmdif_ODE_linear
