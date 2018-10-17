program extract
  IMPLICIT NONE
  REAL(KIND=8) :: err1, err10, err2, err20, erri, erri0
  CHARACTER(LEN=40) :: fich, stuff1, stuff2, stuff3, stuff4
  CHARACTER(LEN=1)  :: c
  CHARACTER(LEN=2)  :: cc
  INTEGER :: nb_fich, i, k, i0, i1
  open(unit=10,file='data_file',form='formatted',status='unknown')
  OPEN(unit=12,file='error',form='formatted',status='unknown')
  READ(10,*) nb_fich
  DO i = 1, nb_fich
     READ(10,*) fich
     OPEN(unit=11,file=fich,form='formatted',status='unknown')
     DO k= 1, 1
        READ(11,*) stuff1
     END DO
     READ(11,*) stuff1, stuff2, stuff3, stuff4, err1
     READ(11,*) stuff1, stuff2, stuff3, stuff4, err2
     READ(11,*) stuff1, stuff2, stuff3, stuff4, erri
     CLOSE(11)
     IF (i==1) THEN
        i0 = 100*2**(i-1)
        WRITE(12,'(I4,3(A,es8.2),A)') i0, ' & ', err1, ' & --  & ', err2, ' & --  & ', erri, ' & --  & \\ \hline'

     ELSE
        !backspace(12)
        !read(12,*) i0, c, err10, c, cc, c, err20, c, cc, c, erri0
        !write(*,*) i0, err10
        i0 = 100*2**(i-2)
        i1 = 2*i0
        WRITE(12,'(I4,3(A,es8.2,A,f4.2),A)') i1, ' & ', err1, ' & ', log(err10/err1)/log(i1/float(i0)), &
             ' & ', err2, ' & ', log(err20/err2)/log(i1/float(i0)), &
             ' & ', erri, ' & ', log(erri0/erri)/log(i1/float(i0)), ' \\ \hline'
     END IF
     i0 = i1
     err10 = err1
     err20 = err2
     erri0 = erri
  END DO
  
end program extract
