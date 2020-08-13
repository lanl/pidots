SUBROUTINE symtrz(g,symmvec)

!---------------------------------------------------------------
!
! Symmetrizes the (I-J_phi) matrix
!  Works for tpose = 1 or 0 (on or off)
!
!---------------------------------------------------------------

USE invar
USE solvar
IMPLICIT NONE
INTEGER, INTENT(IN) :: g
INTEGER :: i, j, k, m, t1, t2, ieq, ieqo, l, ll, tt
REAL*8 :: symm
REAL*8, DIMENSION(neq), INTENT(OUT) :: symmvec

IF (tpose == 1) THEN
   ! Symmetrize jmat and put q into sv
   DO k = 1, nz
      DO j = 1, ny
         DO i = 1, nx
            m = mat(i,j,k)
            t1 = ((i-1) + (j-1)*nx + (k-1)*nx*ny)*nmom
            DO l = 0, anord
               symm = sigs(l,m,g,g)*dx(i)*dy(j)*dz(k)*((-1.0)**l)
               t2 = l**2 + 1
               ieq = t1 + t2
               symmvec(ieq) = symm
               DO tt = 1, neq
                  jmat(tt,ieq) = symmvec(ieq)*jmat(tt,ieq)
               END DO
               DO ll = 1, l
                  t2 = l**2 + 2*ll
                  ieq = t1 + t2
                  ieqo = ieq + 1
                  symmvec(ieq) = 2.0*symm           ! even moment
                  DO tt = 1, neq
                     jmat(tt,ieq) = symmvec(ieq)*jmat(tt,ieq)
                  END DO
                  symmvec(ieqo) = 2.0*symm          ! odd moment
                  DO tt = 1, neq
                     jmat(tt,ieqo) = symmvec(ieqo)*jmat(tt,ieqo)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
ELSE
   DO k = 1, nz
      DO j = 1, ny
         DO i = 1, nx
            m = mat(i,j,k)
            t1 = ((i-1) + (j-1)*nx + (k-1)*nx*ny)*nmom
            DO l = 0, anord
               symm = sigs(l,m,g,g)*dx(i)*dy(j)*dz(k)*((-1.0)**l)
               t2 = l**2 + 1
               ieq = t1 + t2
               symmvec(ieq) = symm
               DO tt = 1, neq
                  jmat(ieq,tt) = symmvec(ieq)*jmat(ieq,tt)
               END DO
               DO ll = 1, l
                  t2 = l**2 + 2*ll
                  ieq = t1 + t2
                  ieqo = ieq + 1
                  symmvec(ieq) = 2.0*symm           ! even moment
                  DO tt = 1, neq
                     jmat(ieq,tt) = symmvec(ieq)*jmat(ieq,tt)
                  END DO
                  symmvec(ieqo) = 2.0*symm          ! odd moment
                  DO tt = 1, neq
                     jmat(ieqo,tt) = symmvec(ieqo)*jmat(ieqo,tt)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
END IF

RETURN
END SUBROUTINE symtrz
