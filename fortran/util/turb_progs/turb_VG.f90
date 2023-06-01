    !	Copyright 2011 Johns Hopkins University
    !
    !  Licensed under the Apache License, Version 2.0 (the "License");
    !  you may not use this file except in compliance with the License.
    !  You may obtain a copy of the License at
    !
    !      http://www.apache.org/licenses/LICENSE-2.0
    !
    !  Unless required by applicable law or agreed to in writing, software
    !  distributed under the License is distributed on an "AS IS" BASIS,
    !  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    !  See the License for the specific language governing permissions and
    !  limitations under the License.


    program TurbTest

    use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_float, c_f_pointer
    implicit none

    type, bind(c) :: thresholdinfo
        INTEGER(c_int) :: x, y, z
        REAL(c_float) :: value
    end type thresholdinfo

    INTEGER, PARAMETER :: RP=4 ! Number of bytes for REALs (single precision)

    ! ---- Temporal Interpolation Options ----
    INTEGER, PARAMETER :: NoTInt = 0   ! No temporal interpolation
    INTEGER, PARAMETER :: PCHIPInt = 1 ! Piecewise cubic Hermit interpolation in time

    ! ---- Spatial Interpolation Flags for GetVelocity & GetVelocityAndPressure ----
    INTEGER, PARAMETER :: NoSInt = 0 ! No spatial interpolation
    INTEGER, PARAMETER :: Lag4 = 4   ! 4th order Lagrangian interpolation in space
    INTEGER, PARAMETER :: Lag6 = 6   ! 6th order Lagrangian interpolation in space
    INTEGER, PARAMETER :: Lag8 = 8   ! 8th order Lagrangian interpolation in space

    ! ---- Spatial Differentiation & Interpolation Flags for GetVelocityGradient & GetPressureGradient ----
    INTEGER, PARAMETER :: FD4NoInt = 40 ! 4th order finite differential scheme for grid values, no spatial interpolation
    INTEGER, PARAMETER :: FD6NoInt = 60 ! 6th order finite differential scheme for grid values, no spatial interpolation
    INTEGER, PARAMETER :: FD8NoInt = 80 ! 8th order finite differential scheme for grid values, no spatial interpolation
    INTEGER, PARAMETER :: FD4Lag4 = 44  ! 4th order finite differential scheme for grid values, 4th order Lagrangian interpolation in space

    ! ---- Spline interpolation and differentiation Flags for getVelocity,
    !      getPressure, getVelocityGradient, getPressureGradient,
    !      getVelocityHessian, getPressureHessian
    INTEGER, PARAMETER :: M1Q4 = 104   ! Splines with smoothness 1 (3rd order) over 4 data points. Not applicable for Hessian.
    INTEGER, PARAMETER :: M2Q8 = 208   ! Splines with smoothness 2 (5th order) over 8 data points.
    INTEGER, PARAMETER :: M2Q14 = 214  ! Splines with smoothness 2 (5th order) over 14 data points.

    !
    ! Choose which dataset to use in this query
    ! Currently, only valid datasets are:
    !   'isotropic1024coarse', 'isotropic1024fine', 'mhd1024', 'channel', 'mixing' and 'isotropic4096'
    !
    CHARACTER(*), PARAMETER :: dataset = 'channel' // CHAR(0)

    !
    ! Access key
    CHARACTER(*), PARAMETER :: authkey = 'edu.cornell.gjr68-f5b59202' // CHAR(0)
    INTEGER, PARAMETER :: nprt = 100, nts = 1000
    REAL(RP), PARAMETER :: pi = 3.1415926535897932384626433_RP

    REAL(RP) :: u(3, nprt), ut(3, nprt), dU(9, nprt), &
                x(3, nprt), xt(3, nprt), dt, time
    INTEGER :: cts, rc, i, n, s
    INTEGER, ALLOCATABLE :: seed(:)
    REAL, ALLOCATABLE :: seedR(:)
    REAL(KIND=8), ALLOCATABLE :: Gtmp(:,:,:,:)

!   Declare the return type of the turblib functions as INTEGER.
!   This is required for custom error handling (see the README).
    INTEGER :: getvelocity, getvelocitygradient

    ALLOCATE(Gtmp(nts, 3, 3, nprt))

!   Seed
    s = 5965965
    CALL RANDOM_SEED(SIZE=n)
    ALLOCATE(seed(n), seedR(n))
    CALL RANDOM_NUMBER(seedR)
    seed = s + INT(n*seedR)
    CALL RANDOM_SEED(PUT=seed)


    dt = 0.02
    time = 0.01_RP

!   Intialize the gSOAP runtime.
!   This is required before any WebService routines are called.
!   
    CALL soapinit()

    ! Enable exit on error.  See README for details.
    CALL turblibSetExitOnError(1)

    DO i = 1,nprt
        CALL RANDOM_NUMBER(x(1,i))
        CALL RANDOM_NUMBER(x(2,i))
        CALL RANDOM_NUMBER(x(3,i))
        x(1,i) = x(1,i)*8_RP*pi
        x(2,i) = x(2,i)*2_RP-1
        x(3,i) = x(3,i)*3_RP*pi
    ENDDO


!   Velocity
    DO cts = 1, nts
!       Initial velocity
        rc = getvelocity(authkey, dataset, time, Lag6, PCHIPInt, nprt, x, u)

!       Time stepped
        xt = x + dt*u

!       Manage periodicity (only in x and z)
        DO i = 1, nprt
            IF(xt(1,i) .gt. 8*pi) xt(1,i) = xt(1,i) - 8*pi
            IF(xt(1,i) .lt. 0) xt(1,i) = xt(1,i) + 8*pi
            IF(xt(3,i) .gt. 3*pi) xt(3,i) = xt(3,i) - 3*pi
            IF(xt(3,i) .lt. 0) xt(3,i) = xt(3,i) + 3*pi
        ENDDO

!       Step
        time = time + dt

!       End point vel
        rc = getvelocity(authkey, dataset, time, Lag6, PCHIPInt, nprt, xt, ut)

!       Heun's
        x = x + dt*(ut + u)*0.5

!       Manage periodicity (only in x and z)
        DO i = 1, nprt
            IF(x(1,i) .gt. 8*pi) x(1,i) = x(1,i) - 8*pi
            IF(x(1,i) .lt. 0) x(1,i) = x(1,i) + 8*pi
            IF(x(3,i) .gt. 3*pi) x(3,i) = x(3,i) - 3*pi
            IF(x(3,i) .lt. 0) x(3,i) = x(3,i) + 3*pi
        ENDDO

!       Newest velgrad
        rc = getvelocitygradient(authkey, dataset, time, FD4Lag4, PCHIPInt, nprt, x, dU)

        print *, 'Timestep: ', cts, ' Time: ', time
        Gtmp(cts, 1, 1:3, :) = dU(1:3, :)
        Gtmp(cts, 2, 1:3, :) = dU(4:6, :)
        Gtmp(cts, 3, 1:3, :) = dU(7:9, :)
    ENDDO

    print *, "Writing binary file"
    ! Binary write
    OPEN(1,FILE = 'VG_turb.bin', ACCESS = 'stream', ACTION = 'write', FORM="unformatted")
    WRITE(1) nts, nprt
    WRITE(1) Gtmp(1:nts,:,:,1:nprt)
    CLOSE(1)
    
!   Destroy the gSOAP runtime.
!   No more WebService routines may be called.
    CALL soapdestroy()

    end program TurbTest

