PROGRAM shallow_water
  USE space_dim
  USE input_data
  USE mesh_handling
  USE update
  USE boundary_conditions
  USE sub_plot
  USE fem_tn
  USE mesh_interpolation
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: rk, un, ui, uo, uinit
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: hmovie
  INTEGER                                   :: it, it_max, kit=0
  REAL(KIND=8)                              :: tps, to, q1=0.d0, q2, q3, dt_frame, t_frame=0.d0
  REAL(KIND=8)                              :: hmax0, norm, dt0
  INTEGER :: nb_frame=0, i, n
  CHARACTER(LEN=200) :: header
  CHARACTER(LEN=3)   :: frame
  LOGICAL :: once=.TRUE.
  inputs%syst_size=k_dim+1+2 !For shallow water+ fake GN
  CALL read_my_data('data')
  CALL construct_mesh
  CALL construct_matrices

  ALLOCATE(rk(inputs%syst_size,mesh%np),un(inputs%syst_size,mesh%np),&
       ui(inputs%syst_size,mesh%np),uo(inputs%syst_size,mesh%np),hmovie(mesh%np),&
       uinit(inputs%syst_size,mesh%np))
  inputs%time =0.d0
  CALL init(un)
  uinit = un

  hmax0 = MAXVAL(un(1,:))
  CALL plot_1d(mesh%rr(1,:), bath, 'bath.plt')
  CALL plot_1d(mesh%rr(1,:), uinit(1,:), 'hinit.plt')
  CALL plot_1d(mesh%rr(1,:), uinit(1,:)+bath, 'hpluszinit.plt')
  CALL plot_1d(mesh%rr(1,:), uinit(2,:), 'qxinit.plt')
  CALL plot_1d(mesh%rr(1,:), uinit(3,:), 'hetainit.plt')
  CALL plot_1d(mesh%rr(1,:), uinit(4,:), 'homegainit.plt')
  CALL COMPUTE_DT(un)

  nb_frame=0
  dt_frame = 1.d10
  t_frame = 1.d10

  WRITE(*,*) ' Mass1', SUM(un(1,:)*lumped)

  it_max = INT(inputs%Tfinal/inputs%dt)
  IF (it_max==0) THEN
     it_max = 1
     inputs%dt = inputs%Tfinal
  ELSE
     inputs%dt = inputs%Tfinal/it_max
  END IF
  dt0 = inputs%dt

  !DO it = 1, it_max
  DO WHILE(inputs%time<inputs%Tfinal)
     CALL COMPUTE_DT(un)
     write(*,*) 'time ', inputs%time, inputs%dt,' Mass1', SUM(un(1,:)*lumped)
     !IF (inputs%dt<0.01*dt0) THEN
     !   WRITE(*,*) ' Time step decresases too much'
     !   STOP
     !END IF

     to = inputs%time
     uo = un  !t
     !===Step 1
     CALL euler(uo,un) !t
     write(*,*) 'time ', inputs%time, inputs%dt, ' Mass2', SUM(uo(1,:)*lumped), SUM(un(1,:)*lumped)
     !if (ABS(SUM(uo(1,:)*lumped)-SUM(un(1,:)*lumped))/SUM(uo(1,:)*lumped).ge.1d-7) STOP
     CALL bdy(un,inputs%time+inputs%dt) !t+dt

     !===Step 2
     inputs%time=to+inputs%dt
     CALL euler(un,ui) !t+dt
     write(*,*) 'time ', inputs%time, inputs%dt, ' Mass3', SUM(un(1,:)*lumped), SUM(ui(1,:)*lumped)
     un = (3*uo+ ui)/4
     CALL bdy(un,inputs%time+inputs%dt/2) !t+dt/2
     write(*,*) 'time ', inputs%time, inputs%dt, ' Mass4', SUM(ui(1,:)*lumped), SUM(un(1,:)*lumped)

     !===Step 3
     inputs%time =  to + inputs%dt/2
     CALL euler(un,ui) !t+dt/2
     write(*,*) 'time ', inputs%time, inputs%dt, ' Mass5', SUM(un(1,:)*lumped), SUM(ui(1,:)*lumped)
     un = (uo+ 2*ui)/3
     CALL bdy(un,inputs%time+inputs%dt) !t+dt
     inputs%time = to + inputs%dt
     write(*,*) 'time ', inputs%time, inputs%dt, ' Mass6', SUM(ui(1,:)*lumped), SUM(un(1,:)*lumped)

     !===Monitor convergence
     SELECT CASE(inputs%type_test)
     CASE(11)
        hmovie = sol_anal(1,mesh%rr,inputs%time)
        CALL ns_l1 (mesh, hmovie-un(1,:), q3)
        CALL ns_l1 (mesh, hmovie, q2)
        WRITE(11,*) 1/inputs%time, q3/q2
     END SELECT

     IF (once) THEN
        tps = user_time()
        once=.FALSE.
     END IF

  END DO
  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step and dof', tps/(it_max*mesh%np), it_max

  CALL compute_errors

  CALL plot_1d(mesh%rr(1,:), un(1,:)+bath, 'HplusZ.plt')
  CALL plot_1d(mesh%rr(1,:), un(1,:)-sol_anal(1,mesh%rr,inputs%time), 'errh.plt')
  CALL plot_1d(mesh%rr(1,:), un(1,:), 'h.plt')
  CALL plot_1d(mesh%rr(1,:), un(1,:)**2, 'h_h.plt')
  CALL plot_1d(mesh%rr(1,:), un(2,:), 'qx.plt')
  CALL plot_1d(mesh%rr(1,:), un(3,:), 'h_eta.plt')
  CALL plot_1d(mesh%rr(1,:), un(4,:), 'h_w.plt')
  CALL plot_1d(mesh%rr(1,:), un(1,:)*velocity(1,:)**2/2, 'rho_e.plt')
  CALL plot_1d(mesh%rr(1,:), un(2,:)-sol_anal(2,mesh%rr,inputs%time), 'errqx.plt')

CONTAINS

  SUBROUTINE COMPUTE_DT(u0)
    USE input_data
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN) :: u0
    !CALL compute_dij(u0,.true.)
    CALL compute_dij(u0,.false.)
    inputs%dt = inputs%CFL*1/MAXVAL(ABS(dij%aa(diag))/lumped)
  END SUBROUTINE COMPUTE_DT

  SUBROUTINE bdy(uu,t)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: uu
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(2,mesh%np) :: vv
    INTEGER :: k
    IF (SIZE(h_js_D).NE.0)  uu(1,h_js_D)  = sol_anal(1,mesh%rr(:,h_js_D),t)
    IF (SIZE(ux_js_D).NE.0) uu(2,ux_js_D) = sol_anal(2,mesh%rr(:,ux_js_D),t)
    DO k = 3, inputs%syst_size
       IF (SIZE(h_js_D).NE.0)  uu(k,h_js_D)  = sol_anal(k,mesh%rr(:,h_js_D),t)
    END DO
  END SUBROUTINE bdy

  FUNCTION user_time() RESULT(time)
    IMPLICIT NONE
    REAL(KIND=8) :: time
    INTEGER :: count, count_rate, count_max
    CALL SYSTEM_CLOCK(COUNT, COUNT_RATE, COUNT_MAX)
    time = (1.d0*count)/count_rate
  END FUNCTION user_time

  SUBROUTINE compute_errors
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np) :: hh
    REAL(KIND=8), DIMENSION(k_dim,mesh%gauss%l_G*mesh%me):: rr_gauss
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me)  :: uexact
    REAL(KIND=8), DIMENSION(mesh%np)                 :: zero
    REAL(KIND=8) :: err, errb, norm, normb, waterh_ref = 0.d0
    SELECT CASE(inputs%type_test)
    CASE(3,4,5,6,7,11,12,13)
       IF (inputs%type_test==13) THEN
          waterh_ref = max_water_h
       END IF

       hh = un(1,:) - waterh_ref
       CALL r_gauss(rr_gauss)
       uexact = sol_anal(1,rr_gauss,inputs%time) - waterh_ref
       zero = 0.d0
       CALL ns_anal_l1(mesh, hh, uexact, err)
       CALL ns_anal_l1(mesh, zero, uexact, norm)
       WRITE(*,*) ' Relative L1 error on h (Gaussian)', err/norm!, err
       !CALL ns_l1 (mesh, hh-un(1,:), err)
       !CALL ns_l1 (mesh, hh, norm)
       !WRITE(*,*) ' Relative L1 error on h', err/norm

       CALL ns_anal_0(mesh, hh, uexact, err)
       CALL ns_anal_0(mesh, zero, uexact, norm)
       WRITE(*,*) ' Relative L2 error on h (Gaussian)', err/norm!, err

       hh = un(1,:)-sol_anal(1,mesh%rr,inputs%time)
       err = MAXVAL(ABS(hh))
       norm = MAXVAL(ABS(sol_anal(1,mesh%rr,inputs%time)-waterh_ref))
       WRITE(*,*) ' Relative Linfty error on h       ', err/norm, waterh_ref!, err

       uexact = sol_anal(2,rr_gauss,inputs%time)
       CALL ns_anal_l1(mesh, un(2,:), uexact, err)
       CALL ns_anal_l1(mesh, zero, uexact, norm)

       WRITE(*,*) ' Relative L1 error on Q (Gaussian)', (err)/(norm)!, (err)

       uexact = sol_anal(2,rr_gauss,inputs%time)
       CALL ns_anal_0(mesh, un(2,:), uexact, err)
       CALL ns_anal_0(mesh, zero, uexact, norm)
       WRITE(*,*) ' Relative L2 error on Q (Gaussian)', (err)/(norm)!, (err)

       hh = ABS(un(2,:)-sol_anal(2,mesh%rr,inputs%time))
       err = MAXVAL(hh)
       hh = ABS(sol_anal(2,mesh%rr,inputs%time))
       norm = MAXVAL(hh)
       WRITE(*,*) ' Relative Linfty error on Q       ', err/norm!, err

    CASE DEFAULT
       RETURN
    END SELECT

  END SUBROUTINE compute_errors

  SUBROUTINE ns_anal_mass_L1 (mesh, ff, f_anal,  t)
    USE Gauss_points
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: f_anal
    REAL(KIND=8),                 INTENT(OUT) :: t
    INTEGER ::  m, l, n, index
    REAL(KIND=8) :: fl, flindex
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    t = 0
    DO m = 1, me
       index = (m-1)*l_G
       DO l = 1, l_G; index = index + 1
          fl = 0
          DO n = 1, n_w
             fl = fl + ff(jj(n,m)) * ww(n,l)
          END DO
          IF (fl.GE. 0.5d0) THEN
             fl = 1.d0
          ELSE
             fl = 0.d0
          END IF
          IF (f_anal(index).GE. 0.5d0) THEN
             flindex = 1.d0
          ELSE
             flindex = 0.d0
          END IF

          t = t + ABS(fl-flindex) * rj(l,m)
       ENDDO
    ENDDO
  END SUBROUTINE ns_anal_mass_L1

  SUBROUTINE r_gauss(rr_gauss)
    USE Gauss_points
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT) :: rr_gauss
    INTEGER :: index, k, l, m, n
    REAL(KIND=8) :: s
    index = 0
    DO m = 1, mesh%me
       DO l = 1,  mesh%gauss%l_G; index = index + 1
          DO k = 1, mesh%gauss%k_d
             s = 0
             DO n=1, mesh%gauss%n_w
                s = s + mesh%gauss%ww(n,l)*mesh%rr(k,mesh%jj(n,m))
             END DO
             rr_gauss(k,index) = s
          END DO
       END DO
    END DO

  END SUBROUTINE r_gauss

  SUBROUTINE plot_1d(rr,un,file)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rr, un
    INTEGER :: n, unit=10
    CHARACTER(*) :: file
    OPEN(unit,FILE=TRIM(ADJUSTL(file)),FORM='formatted')
    WRITE(unit,*) '%toplabel='' '''
    WRITE(unit,*) '%xlabel='' '''
    WRITE(unit,*) '%ylabel='' '''
    WRITE(unit,*) '%ymax=', MAXVAL(un)
    WRITE(unit,*) '%ymin=', MINVAL(un)
    WRITE(unit,*) '%xyratio=1'
    WRITE(unit,*) '%mt=4'
    WRITE(unit,*) '%mc=2'
    WRITE(unit,*) '%lc=2'
    DO n = 1, SIZE(rr)
       WRITE(unit,*) rr(n), un(n)
    END DO
    CLOSE(unit)
  END SUBROUTINE plot_1d
END PROGRAM
