

!   Old/unused vaiable
    ! REAL(KIND = 8) :: cq, normev1, normev2, mab(2,2), mabt(2,2), &
    !                   mabp(2,2), dmab2(2,2),q1, q2, dq(2)
! !           Alternate way, staying in Cartesian


                ! cell%gnR(:,:,i,j)  = gn
                ! cell%dgtR(:,:,i,j) = dgt
                ! cell%dgpR(:,:,i,j) = dgp


! !           First get cartesian derivatives of tau via chain rule
!             tauders(1,:,:) = taut*c1(1) + taup*c2(1)
!             tauders(2,:,:) = taut*c1(2) + taup*c2(2)
!             tauders(3,:,:) = taut*c1(3) + taup*c2(3)

! !           Surface divergence
!             cell%ff(1,i,j) = -SUM(Prj*tauders(:,:,1))
!             cell%ff(2,i,j) = -SUM(Prj*tauders(:,:,2))
!             cell%ff(3,i,j) = -SUM(Prj*tauders(:,:,3))


! !           This is a way to just get the tensions directly in coontravariant form. The tensions
! !           work fine, but for some reason the derivatives are incorrect. I also didn't try that hard
! !           to get it to work. Might be nice to have, but not essential (Walter 2010 Eq. 13)

!             tauab = 0.5D0*(B/(es(1)*es(2))*(I1 + 1)*cell%gnR(:,:,i,j) &
!             + es(1)*es(2)*(C*I2 - B)*gn)

!          dtauabt = (-B/(2D0*es(1)*es(1)*es(2))*(I1 + 1D0)*cell%gnR(:,:,i,j) + 0.5D0*es(2)*(C*I2 - B)*gn)*es1t & 
!                  + (-B/(2D0*es(1)*es(2)*es(2))*(I1 + 1D0)*cell%gnR(:,:,i,j) + 0.5D0*es(1)*(C*I2 - B)*gn)*es2t &
!                  + 0.5D0*B/(es(1)*es(2))*cell%gnR(:,:,i,j)*I1t & ! Strain invariants
!                  + 0.5D0*es(1)*es(2)*C*gn*I2t &
!                  + 0.5D0*B/(es(1)*es(2))*(I1 + 1D0)*cell%dgtR(:,:,i,j) & ! Tensors
!                  + 0.5D0*es(1)*es(2)*(I2 - B)*dgt
!          dtauabp = (-B/(2D0*es(1)*es(1)*es(2))*(I1 + 1D0)*cell%gnR(:,:,i,j) + 0.5D0*es(2)*(C*I2 - B)*gn)*es1p & 
!                  + (-B/(2D0*es(1)*es(2)*es(2))*(I1 + 1D0)*cell%gnR(:,:,i,j) + 0.5D0*es(1)*(C*I2 - B)*gn)*es2p &
!                  + 0.5D0*B/(es(1)*es(2))*cell%gnR(:,:,i,j)*I1p & ! Strain invariants
!                  + 0.5D0*es(1)*es(2)*C*gn*I2p &
!                  + 0.5D0*B/(es(1)*es(2))*(I1 + 1D0)*cell%dgpR(:,:,i,j) & ! Tensors
!                  + 0.5D0*es(1)*es(2)*(I2 - B)*dgp



!           This is the old way I did bending. It was mostly just a qualitative measure
!           of bending, where the moment varied linearly with the difference in curvature
!           between the current and reference state. I don't use it, and the bending model
!           that I do use below is better I believe. The only catch is that the bending model
!           below requires the use of a "spontaneous curvature", a material property that says
!           how much a lipid bilayer bends without any forces. This is a material property,
!           but one person uses it as a varying property so there is no stress in the biconcave
!           resting shape.
!             mab = -(k - cell%kR(i,j))*gn
!             mabt = -(kt - cell%kdR(1,i,j))*gn - (k - cell%kR(i,j))*dgt
!             mabp = -(kp - cell%kdR(2,i,j))*gn - (k - cell%kR(i,j))*dgp
            
!             dmab2(1,1) =-(kt2 - cell%kd2R(1,i,j))*gn(1,1) - 2D0*(kt - cell%kdR(1,i,j))*dgt(1,1) &
!                        - (k - cell%kR(i,j))*dgt2(1,1)
!             dmab2(1,2) =-(ktp - cell%kd2R(3,i,j))*gn(1,2) - (kt - cell%kdR(1,i,j))*dgp(1,2) &
!                        - (kp  - cell%kdR(2,i,j))*dgt(1,2) - (k  - cell%kR(i,j))*dgtp(1,2)
!             dmab2(2,1) =-(ktp - cell%kd2R(3,i,j))*gn(2,1) - (kt - cell%kdR(1,i,j))*dgp(2,1) &
!                        - (kp  - cell%kdR(2,i,j))*dgt(2,1) - (k  - cell%kR(i,j))*dgtp(2,1)
!             dmab2(2,2) =-(kp2 - cell%kd2R(2,i,j))*gn(2,2) - 2D0*(kp - cell%kdR(2,i,j))*dgp(2,2) &
!                        - (k - cell%kR(i,j))*dgp2(2,2)

! !           Transverse shears and partials
!             q1 = Eb*(mabt(1,1) + mabp(2,1) + mab(1,1)*(2D0*c111 + c221) &
!                + mab(2,1)*(2D0*c112 + c222) + mab(1,2)*c112 + mab(2,2)*c122)
!             q2 = Eb*(mabt(1,2) + mabp(2,2) + mab(1,2)*(c111 + 2D0*c221) &
!                + mab(2,2)*(c112 + 2D0*c222) + mab(1,1)*c211 + mab(2,1)*c221)
            
! !           Partials of q
!             dq(1) = Eb*(dmab2(1,1) + dmab2(2,1) &
!                   + mabt(1,1)*(2D0*c111 + c221) + mab(1,1)*(2D0*c111t + c221t) &
!                   + mabt(2,1)*(2D0*c112 + c222) + mab(2,1)*(2D0*c112t + c222t) &
!                   + mabt(1,2)*c112 + mab(1,2)*c112t + mabt(2,2)*c122 + mab(2,2)*c122t)
            
!             dq(2) = Eb*(dmab2(1,2) + dmab2(2,2) &
!                   + mabp(1,2)*(c111 + 2D0*c221) + mab(1,2)*(c111p + 2D0*c221p) &
!                   + mabp(2,2)*(c112 + 2D0*c222) + mab(2,2)*(c112p + 2D0*c222p) &
!                   + mabp(1,1)*c211 + mab(1,1)*c211p + mabp(2,1)*c221 + mab(2,1)*c221p)
!
! !           Covariant divergence of q 
!             cq = dq(1) + dq(2) + c111*q1 + c112*q2 + c221*q1 + c222*q2
