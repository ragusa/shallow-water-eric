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
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: rk, un, ui, uo
  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: hmovie, hetamovie, time_steps
  INTEGER                                   :: it, it_max, kit=0,counter
  REAL(KIND=8)                              :: tps, to, q1, q2, q3, dt_frame, t_frame=0.d0
  REAL(KIND=8)                              :: hmax0
  INTEGER :: nb_frame=200, i, n
  CHARACTER(LEN=200) :: header, etaHeader
  CHARACTER(LEN=3)   :: frame
  CALL read_my_data('data')
  CALL construct_mesh
  CALL construct_matrices
  inputs%syst_size=k_dim+1 !For shallow water

  IF (inputs%type_test==14) THEN
    nb_frame=14
  END IF

  IF (inputs%type_test==13 .OR. inputs%type_test==14 .OR. inputs%type_test==15 &
                  .OR. inputs%type_test==16 .OR. inputs%type_test==17 &
                    .OR. inputs%type_test==18 .OR. inputs%type_test==19) THEN
  inputs%syst_size=k_dim + 3 !For 2d SGN model: h, hu, hv, h*eta, hw. (Eric T.)
  END IF

  ALLOCATE(rk(inputs%syst_size,mesh%np),un(inputs%syst_size,mesh%np),&
       ui(inputs%syst_size,mesh%np),uo(inputs%syst_size,mesh%np),hmovie(mesh%np),hetamovie(mesh%np))
  inputs%time =0.d0
  CALL init(un)
  hmax0 = MAXVAL(un(1,:))
  IF (k_dim==1) THEN
    CALL plot_1d(mesh%rr(1,:), bath, 'bath.plt')
    CALL plot_1d(mesh%rr(1,:), un(1,:), 'hinit.plt')
    CALL plot_1d(mesh%rr(1,:), un(1,:)+bath, 'hpluszinit.plt')
    CALL plot_1d(mesh%rr(1,:), un(2,:), 'huinit.plt')
    CALL plot_1d(mesh%rr(1,:), un(3,:), 'hetainit.plt')
    CALL plot_1d(mesh%rr(1,:), un(4,:), 'hwinit.plt')
  ELSE
    CALL plot_scalar_field(mesh%jj, mesh%rr, bath, 'bath.plt')
    CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'hinit.plt')
  END IF
  ! IF (inputs%type_test==13 .OR. inputs%type_test==14) THEN
  ! CALL plot_scalar_field(mesh%jj, mesh%rr, un(4,:), 'h_eta_init.plt')
  ! CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:), 'h_u_init.plt')
  ! END IF

  CALL COMPUTE_DT(un)

  IF (inputs%type_test==8 .OR. inputs%type_test==5 .OR. inputs%type_test==9 &
  .OR. inputs%type_test==11 .OR. inputs%type_test==12 .OR. inputs%type_test==13 &
  .OR. inputs%type_test==14 .OR. inputs%type_test==15 &
        .OR. inputs%type_test==16 .OR. inputs%type_test==17 &
        .OR. inputs%type_test==18 .OR. inputs%type_test==19) THEN
     dt_frame = inputs%Tfinal/(nb_frame-1)
     CALL vtk_2d(mesh, bath, 10, 'bath.vtk')
  END IF
  ! output exact solution at Tfinal for soliton with flat bath
  IF (inputs%type_test==13) THEN
    CALL vtk_2d(mesh, sol_anal(1,mesh%rr,inputs%Tfinal),11,'hexact.vtk')
  END IF


  tps = user_time()
  it_max = INT(inputs%Tfinal/inputs%dt)
  IF (it_max==0) THEN
     it_max = 1
     inputs%dt = inputs%Tfinal
  ELSE
     inputs%dt = inputs%Tfinal/it_max
  END IF


  !DO it = 1, it_max
  counter = 0
  DO WHILE(inputs%time<inputs%Tfinal)
    counter = counter + 1
     !===Paraboloid seems to work with inputs%htiny=1d-4
     !inputs%htiny=aspect_ratio*inputs%gravity*inputs%dt**2/2
     !inputs%htiny = 10*aspect_ratio**2*MINVAL(lumped)/max_water_h
     !inputs%htiny = 80*aspect_ratio**2*MINVAL(lumped)/max_water_h
     !write(*,*) ' htiny', inputs%htiny
     IF (inputs%type_test==9) THEN
        inputs%htiny = 1.d-3
        inputs%dt=2.5d-2
     ELSE
        CALL COMPUTE_DT(un)
     END IF

     to = inputs%time
     uo = un  !t
     !===Step 1
     CALL euler(uo,un) !t
     CALL bdy(un,inputs%time+inputs%dt) !t+dt

     !===Step 2
     inputs%time=to+inputs%dt
     CALL euler(un,ui) !t+dt
     un = (3*uo+ ui)/4
     CALL bdy(un,inputs%time+inputs%dt/2) !t+dt/2

     !===Step 3
     inputs%time =  to + inputs%dt/2
     CALL euler(un,ui) !t+dt/2
     un = (uo+ 2*ui)/3
     CALL bdy(un,inputs%time+inputs%dt) !t+dt
     inputs%time = to + inputs%dt
     write(*,*) 'time ', inputs%time, inputs%dt
     !!!! outputing time steps to file!!!!
     IF (inputs%lambdaSGN > 0.d0) THEN
       open (unit = 7, file = "time_steps_GN.txt")
       write (7,*) inputs%dt
     ELSE
       open (unit = 7, file = "time_steps_SW.txt")
       write (7,*) inputs%dt
     END IF

     SELECT CASE(inputs%type_test)
     CASE(1,2,10)
        CALL ns_l1(mesh, un(2,:), q1)
        CALL ns_l1(mesh, un(3,:), q2)
        CALL ns_l1(mesh, un(1,:)*SQRT(un(1,:)), q3)
        WRITE(10,*) inputs%time, 0.5d0*(q1+q2)/(inputs%gravity*q3)
     END SELECT

     IF (inputs%type_test==9) THEN !Malpasset gauges
        IF (inputs%time.LE.inputs%dt) THEN
           DO n = 1, 12
              OPEN(10+n,FILE=TRIM(ADJUSTL(malpasset_file(n))), FORM='formatted')
              WRITE(*,*) TRIM(ADJUSTL(malpasset_title(n))), malpasset_rr(:,n)
              WRITE(10+n,*) TRIM(ADJUSTL(malpasset_title(n))), malpasset_rr(:,n)
           END DO
        END IF
        DO n = 1, 12
           WRITE(10+n,*) inputs%time, &
                SUM((un(1,mesh%jj(:,malpasset_m(n)))) &
                * FE_interpolation(mesh,malpasset_m(n),malpasset_rr(1:2,n))), &
                SUM((un(1,mesh%jj(:,malpasset_m(n)))+ bath(mesh%jj(:,malpasset_m(n)))) &
                * FE_interpolation(mesh,malpasset_m(n),malpasset_rr(1:2,n)))
        END DO
     END IF

     IF (inputs%type_test==16) THEN !seawall gauges
        IF (inputs%time.LE.inputs%dt) THEN
           DO n = 1, 7
              OPEN(20+n,FILE=TRIM(ADJUSTL(seawall_file(n))), FORM='formatted')
              ! WRITE(*,*) TRIM(ADJUSTL(seawall_title(n))), seawall_rr(:,n)
              ! WRITE(20+n,*) TRIM(ADJUSTL(seawall_title(n))), seawall_rr(:,n)
              !WRITE(*,*) TRIM(ADJUSTL(seawall_title(n))), '  t  ', 'h  ', 'hpz'
              !WRITE(20+n,*) TRIM(ADJUSTL(seawall_title(n))), '  t  ', 'h  ', 'hpz'
           END DO
        END IF
        DO n = 1, 7
           WRITE(20+n,*) inputs%time, &
                SUM(un(1,mesh%jj(:,seawall_m(n))) &
                * FE_interpolation(mesh,seawall_m(n),seawall_rr(1:2,n))), &
                SUM((un(1,mesh%jj(:,seawall_m(n)))+ bath(mesh%jj(:,seawall_m(n)))) &
                * FE_interpolation(mesh,seawall_m(n),seawall_rr(1:2,n)))
        END DO
     END IF

     IF (inputs%type_test==8 .OR. inputs%type_test==5 .OR. inputs%type_test==9 .OR. &
     inputs%type_test==11 .OR. inputs%type_test==12 .OR. inputs%type_test==13 &
     .OR. inputs%type_test==14 .OR. inputs%type_test==15 &
            .OR.  inputs%type_test==16 .OR. inputs%type_test==17 &
            .OR. inputs%type_test==18 .OR. inputs%type_test==19) THEN
        IF (0.d0 .LE. inputs%time) THEN
           IF (inputs%time.GE.t_frame-1.d-10) THEN
              kit=kit+1
              t_frame = t_frame+dt_frame
              DO i = 1, SIZE(bath)
                 IF (un(1,i).LE. 1.d-4*hmax0) THEN
                    hmovie(i) = -1.d-7*hmax0+bath(i) !0.32
                    hetamovie(i) = -1.d-7*hmax0+bath(i)
                 ELSE
                    hmovie(i) = un(1,i)+bath(i)
                    hetamovie(i) = un(4,i)+bath(i)
                 END IF
              END DO

              WRITE(frame,'(I3)') kit
              header = 'hpz_'//trim(adjustl(frame))//'.vtk'
              etaHeader = 'hetapz_'//trim(adjustl(frame))//'.vtk'
              CALL vtk_2d(mesh,hetamovie,13,etaHeader)
              CALL vtk_2d(mesh, hmovie, 10, header)
              header = 'h_'//trim(adjustl(frame))//'.vtk'
              etaHeader='heta_'//trim(adjustl(frame))//'.vtk'

              DO i = 1, SIZE(bath)
                 IF (un(1,i).LE. 1.d-4*hmax0) THEN
                    hmovie(i) = 0.d0
                    hetamovie(i) = 0.d0
                 ELSE
                    hmovie(i) = un(1,i)
                    hetamovie(i) = un(4,i)
                 END IF
              !hetamovie = un(4,:)
              END DO
              CALL vtk_2d(mesh, hmovie, 10, header)
              CALL vtk_2d(mesh,hetamovie,13,etaHeader)
              !CALL plot_scalar_field(mesh%jj, mesh%rr, hmovie, 'h_'//trim(adjustl(frame))//'.plt')
           END IF
        END IF
     END IF

  END DO

  tps = user_time() - tps
  WRITE(*,*) 'total time', tps, 'Time per time step and dof', tps/(it_max*mesh%np), it_max
  WRITE(*,*) 'total time in minutes', tps/60.d0
  WRITE(*,*) 'Total iterations', counter
  ALLOCATE(time_steps(counter))
  !OPEN(unit = 7, file = 'time_steps_SW.txt', status='old')
  !READ(7,*) time_steps
  !WRITE(*,*) time_steps(1)
  !WRITE(*,*) 'Average Time Step', SUM(time_steps(:))/SIZE(time_steps(:))
  WRITE(*,*) 'Number of elemets', mesh%me
  WRITE(*,*) 'Number of nodes  ', mesh%np
  !WRITE(*,*) 'Average Mesh Size', inputs%avgMeshSize

  CALL compute_errors
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:)+bath, 'HplusZ.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(1,:), 'h.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(2,:), 'qx.plt')
  CALL plot_scalar_field(mesh%jj, mesh%rr, un(3,:), 'qy.plt')
  ! IF (inputs%type_test==13 .OR. inputs%type_test==14) THEN
  !   CALL plot_scalar_field(mesh%jj,mesh%rr,un(4,:),'h_eta.plt')
  !   CALL plot_scalar_field(mesh%jj,mesh%rr,un(5,:),'h_w.plt')
  ! END IF

CONTAINS

  SUBROUTINE COMPUTE_DT(u0)
    USE input_data
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN) :: u0
    CALL compute_dij(u0)
    inputs%dt = inputs%CFL*1/MAXVAL(ABS(dij%aa(diag))/lumped)
  END SUBROUTINE COMPUTE_DT

  SUBROUTINE bdy(uu,t)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: uu
    REAL(KIND=8), INTENT(IN) :: t
    IF (SIZE(h_js_D).NE.0)  uu(1,h_js_D)  = sol_anal(1,mesh%rr(:,h_js_D),t)
    IF (SIZE(ux_js_D).NE.0) uu(2,ux_js_D) = sol_anal(2,mesh%rr(:,ux_js_D),t)
    IF (SIZE(uy_js_D).NE.0) uu(3,uy_js_D) = sol_anal(3,mesh%rr(:,uy_js_D),t)
    ! to add Boundary Conditions to eta and w ! actually need to fix this
    IF (inputs%type_test==13 .OR. inputs%type_test==14 .OR. inputs%type_test==15 &
                      .OR. inputs%type_test==16 .OR. inputs%type_test==17 &
        .OR. inputs%type_test==18  .OR. inputs%type_test==19) THEN
      IF (SIZE(ux_js_D).NE.0) uu(4,ux_js_D) = sol_anal(4,mesh%rr(:,ux_js_D),t)
      IF (SIZE(uy_js_D).NE.0) uu(5,uy_js_D) = sol_anal(5,mesh%rr(:,uy_js_D),t)
    END IF
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
    REAL(KIND=8) :: err, errb, norm, normb
    SELECT CASE(inputs%type_test)
    CASE(3,4,5,6,7,11,12,13)
       CALL r_gauss(rr_gauss)
       uexact = sol_anal(1,rr_gauss,inputs%time)
       zero = 0.d0
       CALL ns_anal_l1(mesh, un(1,:), uexact, err)
       CALL ns_anal_l1(mesh, zero, uexact, norm)
       WRITE(*,*) ' Total number of elements', mesh%me
       WRITE(*,*) ' Total number of points  ', mesh%np
       WRITE(*,*) ' Relative L1 error on h (Gaussian)', err/norm
       !hh = sol_anal(1,mesh%rr,inputs%time)
       !CALL ns_l1 (mesh, hh-un(1,:), err)
       !CALL ns_l1 (mesh, hh, norm)
       !WRITE(*,*) ' Relative L1 error on h', err/norm

       CALL ns_anal_0(mesh, un(1,:), uexact, err)
       CALL ns_anal_0(mesh, zero, uexact, norm)
       WRITE(*,*) ' Relative L2 error on h (Gaussian)', err/norm
       !hh = sol_anal(1,mesh%rr,inputs%time)
       !CALL ns_0 (mesh, hh-un(1,:), err)
       !CALL ns_0 (mesh, hh, norm)
       !WRITE(*,*) ' Relative L2 error on h', err/norm


       uexact = sol_anal(2,rr_gauss,inputs%time)
       CALL ns_anal_l1(mesh, un(2,:), uexact, err)
       CALL ns_anal_l1(mesh, zero, uexact, norm)
       uexact = sol_anal(3,rr_gauss,inputs%time)
       CALL ns_anal_l1(mesh, un(3,:), uexact, errb)
       CALL ns_anal_l1(mesh, zero, uexact, normb)
       !WRITE(*,*) ' Relative L1 error on q (Gaussian)', (err+errb)/(norm+normb)
       !hh = sol_anal(2,mesh%rr,inputs%time)
       !CALL ns_l1 (mesh, hh-un(2,:), err)
       !CALL ns_l1 (mesh, hh, norm)
       !hh = sol_anal(3,mesh%rr,inputs%time)
       !CALL ns_l1 (mesh, hh-un(3,:), errb)
       !CALL ns_l1 (mesh, hh, normb)
       !WRITE(*,*) ' Relative L1 error on q', (err+errb)/(norm+normb)

       uexact = sol_anal(2,rr_gauss,inputs%time)
       CALL ns_anal_0(mesh, un(2,:), uexact, err)
       CALL ns_anal_0(mesh, zero, uexact, norm)
       uexact = sol_anal(3,rr_gauss,inputs%time)
       CALL ns_anal_0(mesh, un(3,:), uexact, err)
       CALL ns_anal_0(mesh, zero, uexact, normb)
       !WRITE(*,*) ' Relative L2 error on q (Gaussian)', (err+errb)/(norm+normb)
       !hh = sol_anal(2,mesh%rr,inputs%time)
       !CALL ns_0 (mesh, hh-un(2,:), err)
       !CALL ns_0 (mesh, hh, norm)
       !hh = sol_anal(3,mesh%rr,inputs%time)
       !CALL ns_0 (mesh, hh-un(3,:), errb)
       !CALL ns_0 (mesh, hh, normb)
       !WRITE(*,*) ' Relative L2 error on q', (err+errb)/(norm+normb)

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
END PROGRAM shallow_water
