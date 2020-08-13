SUBROUTINE conab(sgm,sge,sgx)

!-------------------------------------------------------------
!
!  Construct the A and B matrices that make the gamma matrix
!  A and B matrices specific for an ordinate (n) and cell (x,y,z)
! 
!-------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: sgm, sge, sgx
INTEGER :: t, u, v, ieq, tt, col, mltx, mlty, mltz, mltxx, mltyy, mltzz

! Initialize the amat and bmat values
amat = 0.0
bmat = 0.0

! Set ieq
ieq = 0

! Begin constructing the matrices by advancing through the equations
! From the AHOT-N equations
DO t = 0, lambda
   mltx = sgm**t
   mltxx = ((-1)**t)*mltx 
   DO u = 0, lambda
      mlty = sge**u
      mltyy = ((-1)**u)*mlty
      DO v = 0, lambda
         mltz = sgx**v
         mltzz = ((-1)**v)*mltz
    
         ! Start with the 000 equation, then 001, 002, ..., 010, ..., 100, 
         ieq = ieq + 1
      
         ! amat contribution from total interaction
         amat(ieq,ieq) = amat(ieq,ieq) + 1.0
      
         ! contribution from x-dir summation
         DO tt = MOD((t+1),2), (t-1), 2
            col = ordsq*tt + order*u + v + 1
            amat(ieq,col) = amat(ieq,col) - sgm*(2.0*tt+1.0)/ex
         END DO
      
         ! Contribution from y-dir summation
         DO tt = MOD((u+1),2), (u-1), 2
            col = ordsq*t + order*tt + v + 1
            amat(ieq,col) = amat(ieq,col) - sge*(2.0*tt+1.0)/ey
         END DO

         ! Contribution from z-dir summation
         DO tt = MOD((v+1),2), (v-1), 2
            col = ordsq*t + order*u + tt + 1
            amat(ieq,col) = amat(ieq,col) - sgx*(2.0*tt+1.0)/ez
         END DO
      
         ! Contribution from outgoing z flux (amat) and incoming z flux (bmat)
         col = ordcb + order*t + u + 1
         amat(ieq,col) = amat(ieq,col) + mltz/(2.0*ez)
         bmat(ieq,col) = bmat(ieq,col) + mltzz/(2.0*ez)

         ! Contribution from outgoing y flux (amat) and incoming y flux (bmat)
         col = ordcb + ordsq + order*t + v + 1
         amat(ieq,col) = amat(ieq,col) + mlty/(2.0*ey)
         bmat(ieq,col) = bmat(ieq,col) + mltyy/(2.0*ey)
      
         ! Contribution from outgoing x flux (amat) and incoming x flux (bmat)
         col = ordcb + 2*ordsq + order*u + v + 1
         amat(ieq,col) = amat(ieq,col) + mltx/(2.0*ex)
         bmat(ieq,col) = bmat(ieq,col) + mltxx/(2.0*ex)
      
         ! bmat contribution from scattering and fixed source
         bmat(ieq,ieq) = bmat(ieq,ieq) + 1.0
      
      ! Finished with the AHOT-N equations
      END DO
   END DO
END DO

! Contributions from the WDD equations
! z-direction
DO t = 0, lambda
   DO u = 0, lambda
      ieq = ieq + 1
      ! Contributions to amat from even summations
      DO tt = 0, lambda, 2
         col = ordsq*t + order*u + tt + 1
         amat(ieq,col) = amat(ieq,col) + (2.0*tt + 1.0)
      END DO
      ! Contributions to amat from odd summations
      DO tt = 1, lambda, 2
         col = ordsq*t + order*u + tt + 1
         amat(ieq,col) = amat(ieq,col) + sgx*gamma*(2.0*tt+1.0)
      END DO
      ! Contributions from outgoing flux
      col = ordcb + order*t + u + 1
      amat(ieq,col) = amat(ieq,col) - (1.0+gamma)/2.0
      ! Contribution to bmat from incoming flux
      bmat(ieq,col) = bmat(ieq,col) + (1.0-gamma)/2.0
   ! Done with z-direction
   END DO
END DO
! y-direction
DO t = 0, lambda
   DO v = 0, lambda
      ieq = ieq + 1
      ! Contributions to amat from even summations
      DO tt = 0, lambda, 2
         col = ordsq*t + order*tt + v + 1
         amat(ieq,col) = amat(ieq,col) + (2.0*tt + 1.0)
      END DO
      ! Contributions to amat from odd summations
      DO tt = 1, lambda, 2
         col = ordsq*t + order*tt + v + 1
         amat(ieq,col) = amat(ieq,col) + sge*beta*(2.0*tt+1.0)
      END DO
      ! Contribution from outgoing flux
      col = ordcb + ordsq + order*t + v + 1
      amat(ieq,col) = amat(ieq,col) - (1.0+beta)/2.0
      ! Contribution to bmat from incoming flux
      bmat(ieq,col) = bmat(ieq,col) + (1.0-beta)/2.0
   ! Done with y-direction
   END DO
END DO
! x-direction
DO u = 0, lambda
   DO v = 0, lambda
      ieq = ieq + 1
      ! Contributions to amat from even summations
      DO tt = 0, lambda, 2
         col = ordsq*tt + order*u + v + 1
         amat(ieq,col) = amat(ieq,col) + (2.0*tt + 1.0)
      END DO
      ! Contributions to amat from odd summations
      DO tt = 1, lambda, 2
         col = ordsq*tt + order*u + v + 1
         amat(ieq,col) = amat(ieq,col) + sgm*alpha*(2.0*tt+1.0)
      END DO
      ! Contribution from outgoing flux
      col = ordcb + 2*ordsq + order*u + v + 1
      amat(ieq,col) = amat(ieq,col) - (1.0+alpha)/2.0
      ! Contribution to bmat from incoming flux
      bmat(ieq,col) = bmat(ieq,col) + (1.0-alpha)/2.0
   ! Done with x-direction
   END DO
END DO

RETURN
END SUBROUTINE conab
