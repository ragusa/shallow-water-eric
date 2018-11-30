MODULE mesh_interpolation

  PUBLIC :: coupe, mesh_interp, inside, FE_interpolation
  PRIVATE

  !============================================================================!
  REAL(KIND=8) ::  eps_diam_ref = 1.d-12  ! 1.d-14 Controle d'erreur d'arrondi !
  !                pour determiner si r_old est un vertex du triangle m_old    !
  !============================================================================!

CONTAINS


  SUBROUTINE mesh_interp(mesh,new_mesh,field,new_field,new)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh, new_mesh
    REAL(KIND=8), DIMENSION(:),   INTENT(IN ):: field
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT):: new_field
    LOGICAL,                      INTENT(IN) :: new

    REAL(KIND=8), DIMENSION(2)               :: r_init, r_goal
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: m_goal
    REAL(KIND=8) ::  valeur
    INTEGER :: m, n, m_init
    LOGICAL :: ok
    LOGICAL, SAVE :: once = .false.

    IF (.NOT.new .AND. .NOT.once) THEN
       WRITE(*,*)' BUG mesh_interp: initialization not done'
       STOP
    END IF

    IF (new) THEN
       IF (ALLOCATED(m_goal)) DEALLOCATE(m_goal)
       ALLOCATE(m_goal(new_mesh%np))
       r_init = new_mesh%rr(:,1)
       DO m = 1, mesh%me ! search of element containing r
          IF (inside(mesh,m,r_init)) THEN
             m_init = m
             EXIT
          END IF
       END DO
       m_goal(1) = m_init
       r_init(1) = SUM(mesh%rr(1,mesh%jj(1:3,m_init)))/3
       r_init(2) = SUM(mesh%rr(2,mesh%jj(1:3,m_init)))/3

       DO n = 2, new_mesh%np
          r_goal = new_mesh%rr(:,n)
          CALL find_elem(mesh, r_init, m_init, r_goal, m_goal(n), ok)
          IF (.NOT.ok) THEN
             WRITE(*,*) ' On sort du domaine'
             RETURN
          END IF
          !r_init = r_goal
          r_init(1) = SUM(mesh%rr(1,mesh%jj(1:3,m_goal(n))))/3
          r_init(2) = SUM(mesh%rr(2,mesh%jj(1:3,m_goal(n))))/3
          m_init = m_goal(n)
          !write(*,*) ' m_goal', m_goal(n),  r_goal, inside(mesh,m_goal(n),r_goal)
       END DO
       once = .TRUE.
    END IF

    DO n = 1, new_mesh%np
       r_goal = new_mesh%rr(:,n)
       new_field(n) = SUM(field(mesh%jj(:,m_goal(n))) * ff(mesh,m_goal(n),r_goal))
       !write(*,*) '  new_field(n)',  field(mesh%jj(:,m_goal(n))), new_field(n)
    END DO

100 FORMAT(2(e11.5,2X))
    CLOSE(20)

  END SUBROUTINE mesh_interp

  SUBROUTINE coupe(mesh,r_0,r_1,n_point,field,file_name)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh
    REAL(KIND=8), DIMENSION(2),   INTENT(IN) :: r_0, r_1
    INTEGER,                      INTENT(IN) :: n_point
    REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: field
    CHARACTER(*),                 INTENT(IN) :: file_name

    REAL(KIND=8), DIMENSION(2)       :: r_init, r_goal
    REAL(KIND=8) :: s, d, valeur
    INTEGER :: m, n, m_init, m_goal
    LOGICAL :: ok

    OPEN (UNIT = 20, FILE = file_name, FORM = 'formatted', STATUS = 'unknown')
    WRITE (20,*) "$DATA=CURVE2D"
    WRITE(20,'(2(A,e11.5,2x,e11.5),A)') '%toplabel = " x0=',r_0,',  x1=',r_1,' " '
    WRITE(20,*) ' %xlabel = " " '
    WRITE(20,*) ' %ylabel = " " '
    WRITE(20,*) ' '


    d = SQRT(pd_scal(r_0-r_1,r_0-r_1))
    r_init = r_0
    DO m = 1, mesh%me ! search of element containing r

       IF (inside(mesh,m,r_init)) THEN
          m_init = m
          EXIT
       END IF

    END DO

    s = 0.d0
    valeur = SUM(field(mesh%jj(:,m_init)) * ff(mesh,m_init,r_init))
    WRITE(20,100) s, valeur

    DO n = 2, n_point
       s = d*(n-1.d0)/(n_point-1.d0)
       r_goal = r_0 + (n-1)*(r_1-r_0)/(n_point-1)
       CALL find_elem(mesh, r_init, m_init, r_goal, m_goal, ok)
       IF (.NOT.ok) THEN
          WRITE(*,*) ' La coupe sort du domaine'
          RETURN
       END IF
       r_init = r_goal
       m_init = m_goal
       write(*,*) ' m_goal',m_goal, ok,  r_goal
       valeur = SUM(field( mesh%jj(:,m_goal)) * ff(mesh,m_goal,r_goal))
       WRITE(20,100) s, valeur

    END DO

100 FORMAT(2(e11.5,2X))
    CLOSE(20)

  END SUBROUTINE coupe

  FUNCTION pd_scal(x,y)

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(2), INTENT(IN) :: x, y
    REAL(KIND=8)                           :: pd_scal

    pd_scal = x(1)*y(1) + x(2)*y(2)

  END FUNCTION pd_scal


  FUNCTION pd_vect(y)

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(2), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(2)             :: pd_vect

    pd_vect(1) = -y(2)
    pd_vect(2) = y(1)

  END FUNCTION pd_vect

  FUNCTION ff(mesh, m,r)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                        :: mesh
    INTEGER,        	        INTENT(IN) :: m
    REAL(KIND=8), DIMENSION(2), INTENT(IN) :: r

    REAL(KIND=8), DIMENSION(mesh%gauss%n_w):: ff

    REAL(KIND=8), DIMENSION(2) :: r1, r2, r3, &
         n1, n2, n3

    REAL(KIND=8) :: x, y, z, t, one=1, two=2, four=4, half  = 0.5d0

    r1 = mesh%rr(:,mesh%jj(1,m))
    r2 = mesh%rr(:,mesh%jj(2,m))
    r3 = mesh%rr(:,mesh%jj(3,m))

    n1 = pd_vect(r3-r2)   ! inward normal
    n2 = pd_vect(r1-r3)   ! inward normal
    n3 = pd_vect(r2-r1)   ! inward normal

    x = 1 - pd_scal(n2, r2 - r) / pd_scal(n2, r2 - r3)
    y = 1 - pd_scal(n3, r3 - r) / pd_scal(n3, r3 - r1)

    !IF (x<-1.d-15 .OR. y<-1.d-15 .OR. x-1.d0>1.d-15 .OR. y-1.d0>1.d-15 &
    !     .OR. 1-x-y<-1.d-15 .OR. 1-x-y-1.d0>1.d-15) THEN
    !   WRITE(*,*) ' BUG in ff'
    !   WRITE(*,*) x, y, 1-x-y,(x<-1.d-15),y<-1.d-15,x-1.d0>1.d-15,y> 1.d0,1-x-y<-1.d-15, 1-x-y>1.d0
    !   STOP
    !END IF

    IF (mesh%gauss%n_w==3) THEN
       ff(1) = 1-x-y
       ff(2) = x
       ff(3) = y
    ELSE IF (mesh%gauss%n_w==6) THEN
       ff(1) = (half - x - y) * (one - x - y) * two
       ff(2) = x * (x - half) * two
       ff(3) = y * (y - half) * two
       ff(4) = x * y * four
       ff(5) = y * (one - x - y) * four
       ff(6) = x * (one - x - y) * four
    ELSE
       WRITE(*,*) ' BUG in ff '
       STOP
    END IF


  END FUNCTION ff

    FUNCTION FE_interpolation(mesh, m,r) RESULT(ff)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                        :: mesh
    INTEGER,        	        INTENT(IN) :: m
    REAL(KIND=8), DIMENSION(2), INTENT(IN) :: r

    REAL(KIND=8), DIMENSION(mesh%gauss%n_w):: ff

    REAL(KIND=8), DIMENSION(2) :: r1, r2, r3, &
         n1, n2, n3

    REAL(KIND=8) :: x, y, z, t, one=1, two=2, four=4, half  = 0.5d0

    r1 = mesh%rr(:,mesh%jj(1,m))
    r2 = mesh%rr(:,mesh%jj(2,m))
    r3 = mesh%rr(:,mesh%jj(3,m))

    n1 = pd_vect(r3-r2)   ! inward normal
    n2 = pd_vect(r1-r3)   ! inward normal
    n3 = pd_vect(r2-r1)   ! inward normal

    x = 1 - pd_scal(n2, r2 - r) / pd_scal(n2, r2 - r3)
    y = 1 - pd_scal(n3, r3 - r) / pd_scal(n3, r3 - r1)

    !IF (x<-1.d-15 .OR. y<-1.d-15 .OR. x-1.d0>1.d-15 .OR. y-1.d0>1.d-15 &
    !     .OR. 1-x-y<-1.d-15 .OR. 1-x-y-1.d0>1.d-15) THEN
    !   WRITE(*,*) ' BUG in ff'
    !   WRITE(*,*) x, y, 1-x-y,(x<-1.d-15),y<-1.d-15,x-1.d0>1.d-15,y> 1.d0,1-x-y<-1.d-15, 1-x-y>1.d0
    !   STOP
    !END IF

    IF (mesh%gauss%n_w==3) THEN
       ff(1) = 1-x-y
       ff(2) = x
       ff(3) = y
    ELSE IF (mesh%gauss%n_w==6) THEN
       ff(1) = (half - x - y) * (one - x - y) * two
       ff(2) = x * (x - half) * two
       ff(3) = y * (y - half) * two
       ff(4) = x * y * four
       ff(5) = y * (one - x - y) * four
       ff(6) = x * (one - x - y) * four
    ELSE
       WRITE(*,*) ' BUG in ff '
       STOP
    END IF


  END FUNCTION FE_interpolation

  SUBROUTINE find_elem(mesh, r_init, m_init, r_goal, m_goal, ok, i_dom)
    !
    ! Trouve l'element auquel appartient le point r_goal.
    ! Retourne ok=true si le rayon reste dans le domaine,
    ! et retourne ok=false sinon. On va de r_init a r_goal en ligne droite.
    !
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                        :: mesh

    REAL(KIND=8), DIMENSION(2), INTENT(IN)  :: r_goal, r_init
    INTEGER,                    INTENT(IN)  :: m_init
    INTEGER,                    INTENT(OUT) :: m_goal
    LOGICAL,                    INTENT(OUT) :: ok
    INTEGER,  OPTIONAL,         INTENT(IN)  :: i_dom

    REAL(KIND=8), DIMENSION(2) :: r_old, r_next
    INTEGER :: m_old, m_next, face_old, face_next

    r_old = r_init
    m_old = m_init
    face_old = 0

    IF (inside(mesh,m_old,r_goal)) THEN
       ok = .true.
       m_goal = m_old
       RETURN
    END IF

    DO WHILE(.TRUE.)

       CALL find_next(mesh, r_goal, r_old, face_old, m_old, r_next, face_next, m_next)
       !TEST
       !write(*,*) ' m_next', m_next, ' m_old', m_old
       !TEST
       IF (PRESENT(i_dom)) THEN
          IF (m_next == 0  .OR.  mesh%i_d(m_next) /= i_dom) THEN
             ok = .FALSE.
             RETURN
          END IF
       ELSE
          IF (m_next == 0) THEN
             ok = .FALSE.
             RETURN
          END IF
       END IF
!June 11 2007
       IF ((m_next==m_old) .OR. inside(mesh,m_next,r_goal)) THEN ! To account for curved boundaries.
       !IF (inside(mesh,m_next,r_goal)) THEN
!June 11 2007
          m_goal = m_next
          ok = .TRUE.
          RETURN
       END IF

       m_old = m_next
       r_old = r_next
       face_old = face_next

    END DO

  END SUBROUTINE find_elem


  SUBROUTINE find_next(mesh, r_goal, r_old, face_old, m_old, r_next, face_next, m_next)

    ! Find the next element that is neighboring m_old
    ! and is in the direction of vector r - r_old.
    ! warning : r is asumed to be internal to the mesh.
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                        :: mesh

    REAL(KIND=8), DIMENSION(2), INTENT(IN)   :: r_goal
    REAL(KIND=8), DIMENSION(2), INTENT(INOUT):: r_old
    INTEGER,                    INTENT(IN)   :: face_old, m_old
    INTEGER,                    INTENT(OUT)  :: face_next, m_next
    REAL(KIND=8), DIMENSION(2), INTENT(OUT)  :: r_next

    REAL(KIND=8), DIMENSION(2,3) :: rv, n
    REAL(KIND=8), DIMENSION(2)   :: d
    REAL(KIND=8), DIMENSION(3)   :: s
    INTEGER,      DIMENSION(1)   :: dummy
    REAL(KIND=8)                 :: scal1, scal2, diam=0.d0, eps_diam
    INTEGER                      :: j, jp

    rv = mesh%rr(:,mesh%jj(1:3,m_old)) ! vertex de m_old
    DO j = 1, 3
       diam = MAX(diam,MAXVAL(rv(:,j)-r_old))
    END DO
    eps_diam = eps_diam_ref*diam
    !Test pour voir si r_old n'est pas un vertex du triangle
    DO j = 1, 3
       IF (MAXVAL(ABS(r_old-rv(:,j))) < eps_diam) THEN
          r_old = rv(:,j) + 1.d-10*((rv(:,1)+rv(:,2)+rv(:,3))/3 -r_old)
          EXIT
       END IF
    END DO


    n(:,1) = pd_vect(rv(:,3)-rv(:,2))   ! inward (non-normalized) normal
    n(:,2) = pd_vect(rv(:,1)-rv(:,3))
    n(:,3) = pd_vect(rv(:,2)-rv(:,1))

    d = r_goal - r_old

    DO j = 1, 3
       scal1 = pd_scal(n(:,j), d)
       jp = MODULO(j,3) + 1
       scal2 = pd_scal(n(:,j), rv(:,jp) - r_old)

       !TEST
       !write(*,*) ' AVANT: scal1 ',scal1, ' scal2', scal2
       !TEST

       IF (ABS(scal2) > ABS(scal1)) THEN
          scal1=1.d0
          scal2=1.d20
       END IF

       IF (scal2==0.d0 .and. scal1>0.d0) THEN
          scal1=1.d0
          scal2=1.d20
       END IF

       !TEST
       !write(*,*) ' APRES: scal1 ',scal1, ' scal2', scal2
       !TEST


       !IF (ABS(scal1) .LE. eps) THEN
       !   s(j) = 1.d20
       !ELSE
       s(j) = scal2/scal1
       !ENDIF

       IF (s(j) < 0.d0) s(j) = 1.d20

    END DO

    IF (face_old /= 0) s(face_old) = 1.d20

    dummy = MINLOC(s)
    !TEST
    !write(*,*) 'i On sort par la face ', dummy(1)
    !TEST

    IF (s(dummy(1)) == 1.d20) THEN
       WRITE(*,*) ' BUG : find_next, s(dummy(1)) == 1.d20'
       write(*,*) ' rv(1) ', rv(:,1)
       write(*,*) ' rv(2) ', rv(:,2)
       write(*,*) ' rv(3) ', rv(:,3)
       write(*,*) ' r_old', r_old
       write(*,*) ' r_goal', r_goal
       STOP
    END IF

    IF (dummy(1) == face_old) THEN
       WRITE(*,*) ' BUG : on tente de sortir par ou on est entre'
       WRITE(*,*) ' s ', s, 'face_old ', face_old, 'dummy(1) ', dummy(1)
       WRITE(*,*) ' r_goal ', r_goal, ' r_old ', r_old
       WRITE(*,*) ' m_old ', m_old
       STOP
    END IF

    r_next = r_old + s(dummy(1)) * d
    m_next = mesh%neigh(dummy(1),m_old)

    !TEST
    !write(*,*) ' s(dummy(1)) ', s(dummy(1))
    !TEST
    IF (m_next==0) THEN
       WRITE(*,*) ' BUG: Algorithm tries to get out of the domain'
       WRITE(*,*) ' I assume the boundary is curved and proceed, is that okay?'
       m_next = m_old  ! We stay here
       write(*,*) ' rv(1) ', rv(:,1)
       write(*,*) ' rv(2) ', rv(:,2)
       write(*,*) ' rv(3) ', rv(:,3)
       write(*,*) ' r_old', r_old
       write(*,*) ' r_goal', r_goal
       RETURN
       !STOP
    END IF

    DO j = 1, 3
       IF (mesh%neigh(j,m_next)==m_old) THEN
          face_next = j
          RETURN
       END IF
    END DO


    WRITE(*,*) ' BUG : face_next pas trouv''e', m_next, m_old
    WRITE(*,*) ' neigh(:,m_old) ',mesh%neigh(:,m_old)
    WRITE(*,*) ' neigh(:,m_next) ',mesh%neigh(:,m_next)
    WRITE(*,*) ' l''arete par laquelle on sort ', dummy(1)
    WRITE(*,*) r_goal
    WRITE(*,*) rv(:,1)
    WRITE(*,*) rv(:,2)
    WRITE(*,*) rv(:,3)
    WRITE(*,*) inside(mesh,m_old,r_goal)
    WRITE(*,*)pd_scal(pd_vect(rv(:,3)-rv(:,2)),r_goal-rv(:,2))
    WRITE(*,*)pd_scal(pd_vect(rv(:,1)-rv(:,3)),r_goal-rv(:,3))
    WRITE(*,*)pd_scal(pd_vect(rv(:,2)-rv(:,1)),r_goal-rv(:,1))

    WRITE(*,*) r_old
    WRITE(*,*) inside(mesh,m_old,r_old)


   OPEN (UNIT=20, FILE='test.plt', FORM='formatted', STATUS='unknown')

   WRITE (20, *) '$ DATA = CONTCURVE'
   WRITE (20, *) '% contstyle = 2'
   WRITE (20, *) '% nsteps = 50'
   WRITE (20, *) '% meshplot  = true'
   WRITE (20, *)


      WRITE (20, 100) rv(1,1), rv(2,1), 1., 1
      WRITE (20, 100) rv(1,2), rv(2,2), 1., 2
      WRITE (20, 100) rv(1,3), rv(2,3), 1., 3
      WRITE (20, *)

      WRITE (20, 100) rv(1,1), rv(2,1), 1., 1
      WRITE (20, 100) rv(1,2), rv(2,2), 1., 2
      WRITE (20, 100) r_goal(1), r_goal(2), 1., 4
      WRITE (20, *)


100 FORMAT(3(e11.5,3x),i5)
    STOP

  END SUBROUTINE find_next


  FUNCTION inside(mesh,m,r)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                        :: mesh
    INTEGER,                    INTENT(IN) :: m
    REAL(KIND=8), DIMENSION(2), INTENT(IN) :: r

    LOGICAL                    :: inside
    REAL(KIND=8), DIMENSION(2) :: r1, r2, r3, &
         n1, n2, n3
    REAL(KIND=8) :: eps_ref=0.d0, s1, s2, s3, epsilon=-1.d-12
    r1 = mesh%rr(:,mesh%jj(1,m))
    r2 = mesh%rr(:,mesh%jj(2,m))
    r3 = mesh%rr(:,mesh%jj(3,m))

    inside = .FALSE.

    n1 = pd_vect(r3-r2)
    s1 = pd_scal(n1, r - r2)
    !IF (s1 < eps_ref ) RETURN

    n2 = pd_vect(r1-r3)
    s2 = pd_scal(n2, r - r3)
    !IF (s2 < eps_ref ) RETURN

    n3 = pd_vect(r2-r1)
    s3 = pd_scal(n3, r - r1)
    !IF (s3 < eps_ref ) RETURN

    IF (MIN(s1,s2,s3)<epsilon) RETURN
    inside = .TRUE.

  END FUNCTION inside


END MODULE mesh_interpolation
