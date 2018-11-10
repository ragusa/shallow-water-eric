MODULE smoothing
  USE matrix_type
  TYPE(matrice_bloc), PUBLIC                 :: mat
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: test_in,test_out

CONTAINS

  SUBROUTINE apply_smoothing(test_in,test_out)
    USE st_matrix
    USE mesh_handling
    USE fem_s_M
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:)  :: test_in
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: test_out
    REAL(KIND=8)                :: neighbors_sum
    INTEGER      :: i,j,p,n
    CALL st_csr(mesh%jj, mat%ia, mat%ja)
    ALLOCATE(mat%aa(SIZE(mat%ja)))
    mat%aa = 0.d0
    neighbors_sum = 0.d0
    DO i = 1, mesh%np
      ! WRITE(*,*) 'this is ith node ', i
      ! WRITE(*,*) 'value at original ith node', test_in(i)
      ! WRITE(*,*) 'this is how many neighbors', (mat%ia(i+1) - 1) - mat%ia(i) + 1
      DO p = mat%ia(i), mat%ia(i+1) - 1
         !WRITE(*,*) 'this is p ', p
         n = (mat%ia(i+1) - 1) - mat%ia(i) + 1
         j = mat%ja(p)
         !WRITE(*,*) 'this is a jth neighbor and its value', j, test_in(j)
         !test_out(i) = (n * test_in(i) + test_in(j))/(2.d0 * n)
         neighbors_sum = neighbors_sum + test_in(j)
      END DO
      ! WRITE(*,*) 'this is sum of neighbors  ', neighbors_sum
      ! WRITE(*,*) 'this is old way average           ', test_out(i)
      ! WRITE(*,*) 'this is new way average           ', (n*test_in(i) + neighbors_sum)/(2.d0 * n)
      test_out(i) = (n*test_in(i) + neighbors_sum)/(2.d0 * n)
      ! WRITE(*,*) '--------------------------  '
      !  STOP
      neighbors_sum = 0.d0
    END DO
  END SUBROUTINE apply_smoothing

END MODULE smoothing
