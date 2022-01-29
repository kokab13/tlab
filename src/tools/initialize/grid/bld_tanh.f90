#include "types.h"

!########################################################################
!# The grid spacing dx/ds is given by a hyperbolic tangent function
!#
!# dx/ds = f_0 + [(f_1-f_0)/2]*[ 1 + TANH[(s-s_1)/(2*delta_1)] ]
!#             + [(f_2-f_0)/2]*[ 1 + TANH[(s-s_2)/(2*delta_2)] ]
!#             + ...
!########################################################################
SUBROUTINE BLD_TANH(idir, x, imax, scalex, work)
  USE GRID_LOCAL
  IMPLICIT NONE

  TINTEGER, INTENT(IN   ) :: idir, imax
  TREAL,    INTENT(INOUT) :: x(imax), scalex, work(imax)

  ! -----------------------------------------------------------------------
  TREAL st(3), f(3), delta(3)    ! superposition of up to 3 modes, each with 3 parameters
  TINTEGER i, iloc, iseg, im
  TREAL ds

  ! #######################################################################
  iseg = 1                      ! Assumes just 1 segment

  iloc = 1
  IF ( g_build(idir)%mirrored ) iloc = imax/2 ! mirrored case; first point in array is imax/2

  ! create uniform reference grid s
  ds = scalex /M_REAL(imax-iloc)
  DO i=iloc,imax
    x(i) = M_REAL(i-iloc) *ds
  ENDDO

  ! -----------------------------------------------------------------------
  ! define local variables for readability below; strides of 3
  st    = g_build(idir)%vals(1::3,iseg)     ! transition point in uniform grid
  f     = g_build(idir)%vals(2::3,iseg)     ! ratio dx(i)/dx_0
  delta = g_build(idir)%vals(3::3,iseg)

  ! create mapping from grid s to grid x
  work = C_0_R
  DO im = 1,3                           ! 3 modes at most
    IF ( ABS(delta(im)) > C_0_R ) THEN
      work(iloc:imax) = work(iloc:imax) + (f(im)-C_1_R)*delta(im)*LOG(EXP((x(iloc:imax)-st(im))/delta(im))+C_1_R)
    END IF
  END DO
  x = x +work

  ! correct value of scale
  x(:) = x(:) - x(iloc)
  scalex = x(imax)

  RETURN
END SUBROUTINE BLD_TANH
