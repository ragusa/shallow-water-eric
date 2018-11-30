MODULE boundary_conditions
  !in 2D
  !5=Paraboloid
  !9=Malpaseset
  !10=Well balancing + friction
  !1D
  !11=1D accuracy test + Manning
  !12=1D paraboloid + linear friction
  !13=1D Solitary wave
  !14=1D Solitary wave Run UP
  !15=1D mGN steady state
  !16=1D Solitary wave Seawall
  !17=1D Periodic waves over submerged bar
  USE input_data
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE   :: bath
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: velocity
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: un_over_h
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE   :: one_over_h
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE   :: regul_h, limit_h
  REAL(KIND=8)                              :: max_water_h
  REAL(KIND=8), PARAMETER                   :: pi=3.141592653589793d0
  !===Malpasset
  INTEGER, PARAMETER              :: nb_of_gauges=29
  REAL(KIND=8), DIMENSION(2,nb_of_gauges)   :: malpasset_rr
  INTEGER,      DIMENSION(nb_of_gauges)     :: malpasset_m
  CHARACTER(LEN=20), DIMENSION(nb_of_gauges):: malpasset_file, malpasset_title
  !==Seawall Gauges
  REAL(KIND=8), DIMENSION(1,7)    :: seawall_rr
  INTEGER,      DIMENSION(7)      :: seawall_m
  CHARACTER(LEN=20), DIMENSION(7):: seawall_file, seawall_title
  !==bar Gauges
  REAL(KIND=8), DIMENSION(1,8)    :: bar_rr
  INTEGER,      DIMENSION(8)      :: bar_m
  CHARACTER(LEN=20), DIMENSION(8):: bar_file, bar_title
  INTEGER, DIMENSION(:), ALLOCATABLE :: new_to_old
  INTEGER, DIMENSION(1) :: vmax
  !===
CONTAINS

  FUNCTION compute_divide_by_h(un) RESULT(vv)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,SIZE(un,2)) :: vv
    INTEGER :: n
    one_over_h = compute_one_over_h(un(1,:))
    DO n = 1, mesh%np
       vv(1,n) = 1.d0
       vv(2:,n) = un(2:,n)*one_over_h(n)
    END DO
  END FUNCTION  compute_divide_by_h

  FUNCTION flux(un) RESULT(vv)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,k_dim,mesh%np):: vv
    REAL(KIND=8) :: hh, vel, reg
    INTEGER :: n
    SELECT CASE(inputs%type_test)
    CASE(13,14,15,16,17)
       DO n = 1, mesh%np
          vv(1:inputs%syst_size,1,n) = velocity(1,n)*un(1:inputs%syst_size,n)
       END DO
    CASE DEFAULT
       WRITE(*,*) ' BUG in flux'
       STOP
    END SELECT
  END FUNCTION flux

  FUNCTION compute_one_over_h(un) RESULT(vv)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(mesh%np) :: un
    REAL(KIND=8), DIMENSION(mesh%np) :: vv
    REAL(KIND=8) :: vel_square!
    INTEGER :: n
    LOGICAL, SAVE :: once=.TRUE.
    DO n = 1, mesh%np
       vv(n) = 2*un(n)/(un(n)**2+max(un(n),inputs%htiny)**2)
    END DO
  END FUNCTION compute_one_over_h

  SUBROUTINE init(un)
    USE mesh_handling
    USE fem_tn
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: un
    REAL(KIND=8), DIMENSION(mesh%np) :: aux, rr, x, gauss_height, sechSqd
    INTEGER :: i, k
    REAL(KIND=8) :: ray1, ray2, ray3, scale, x0, h0, a, hL, eta, omega, &
         hcone, rcone, scone, htop, radius, bx, q0, surface, h_l1, &
         h1, h2, slope, SS, z, q, L, cGN, cBer, k_wavenumber, term1, term2, term3
    !===For malpasset
    INTEGER, DIMENSION(3,mesh%me) :: jj_old
    REAL(KIND=8), DIMENSION(mesh%np) :: bath_old
    INTEGER :: m
    !===
    ALLOCATE(bath(mesh%np))
    ALLOCATE(one_over_h(mesh%np))
    ALLOCATE(regul_h(mesh%np), limit_h(mesh%np))
    ALLOCATE(velocity(k_dim,mesh%np))
    ALLOCATE(un_over_h(inputs%syst_size,mesh%np))
    SELECT CASE(inputs%type_test)

    CASE(13) !===Solitary wave
       inputs%gravity=10.d0
       bath = 0.d0

    CASE(14) ! modified hyperbolic SGN solitary wave run up (Eric T., 5/21/2018)
       inputs%gravity = 9.81d0
       h0 = 1.d0
       a = 0.28d0 ! amplitude
       slope = 1.d0 / 19.85d0
       inputs%max_water_h = a + h0
       inputs%Tfinal = inputs%Tfinal*SQRT(h0/inputs%gravity)
       DO i = 1, mesh%np
          bath(i) = MAX((mesh%rr(1,i) ) * slope, -h0)
       END DO

    CASE(15) !modified hyperbolic GN steady state
       ! initial constants go here
       inputs%gravity = 9.81d0
       a = 1.5d0
       L = 1.d0
       h0 = 2.d0
       q = 4.5d0
       inputs%max_water_h = h0 + a
       cBer = h0 + q**2 / (2.d0 * inputs%gravity * h0**2)

       ! gauss_height = h0 + a * EXP(-(x/L)**2)
       ! term1 = a * (L-x)*(L+x)
       ! term2 = EXP((x/L)**2) * h0 * (L**2 - 2.d0 * x**2)
       ! bath = cBer - gauss_height - q**2 / (2.d0 * inputs%gravity * gauss_height**2) &
       !       + 2.d0 * a * q**2 / (3.d0 * inputs%gravity * L**4) *  &
       !       1.d0 / (a + h0*EXP((x/L)**2))**2 * (term1 + term2)

       ! for bath
       x  = mesh%rr(1,:)
       k_wavenumber = SQRT(3.d0 * a/(4.d0 * h0**3)) ! wavenumber
       z = SQRT(3.d0 * a * h0) / (2.d0 * h0 * SQRT(h0 * (1.d0 + a)))
       L = 2.d0 / k_wavenumber * ACOSH(SQRT(1.d0 / 0.05d0)) ! wavelength of solitary wave
       sechSqd = (1.0d0/COSH( z*x ))**2.0d0

       term1 = (1.d0 + a)*(2.d0*inputs%gravity*(cber - h0)*h0**2 - q**2)
       term2 = 2*a*((1 + a)*inputs%gravity*h0**3 - q**2)
       term3 = 2*(1 + a)*inputs%gravity*h0**2

       bath = (term1 - term2 * sechSqd)/term3

    CASE(16) !===Seawall (Eric T., 11/16/2018)
      ! initial constants go here
      inputs%gravity = 9.81d0
      h0 = 0.2d0 ! reference height
      a = 0.35d0 ! amplitude
      max_water_h = a + h0
      inputs%max_water_h = max_water_h


      ! shift x 3 metres to match the shift in origin introduced by Hsiao and Lin
      x = mesh%rr(1,:) + 3.d0
      DO i = 1, mesh%np
        IF (x(i) .LT. 10.d0) THEN
          bath(i) = 0.d0
        ELSE IF (x(i) .GE. 13.6d0 .AND. x(i) .LE. 13.9d0) THEN
          bath(i) = 3.6d0/20.d0 + (x(i)-13.6d0)*0.076d0/(13.9d0 - 13.6d0)
        ELSE IF (x(i) .GE. 13.9d0 .AND. x(i) .LE. 13.948d0) THEN
          bath(i) = 3.6d0/20.d0 + 0.076d0
        ELSE IF (x(i) .GE. 13.948d0 .AND. x(i) .LE. 14.045d0) THEN
          bath(i) = 3.6d0/20.d0 + 0.076d0 - (x(i)-13.948d0)*(0.076d0-0.022d0)/(14.045d0 - 13.948d0)
        ELSE
          bath(i) = (x(i) - 10.d0)/20.d0
        END IF
      END DO
      bath = bath - h0 ! shift bath down by reference height
      !===Find Gauges
      CALL seawall_gauges

    CASE(17) !===Period Waves over bar
      inputs%gravity = 9.81d0
      h0 = 0.03d0
      inputs%max_water_h = h0

      ! define bathymetry here
      x = mesh%rr(1,:)
      DO i = 1, mesh%np
        IF (x(i) .LE. 6.d0) THEN
          bath(i) = -0.4d0
        ELSE IF (x(i) .GE. 6.d0 .AND. x(i) .LE. 12.d0) THEN
          bath(i) = -0.4d0 + (x(i)-6.d0)*1.d0/20.d0
        ELSE IF (x(i) .GE. 12.d0 .AND. x(i) .LE. 14.d0) THEN
          bath(i) = -0.1d0
        ELSE IF (x(i) .GE. 14.d0 .AND. x(i) .LE. 17.d0) THEN
          bath(i) = -0.1d0 - (x(i)-14.d0)*1.d0/10.d0
        ELSE
          bath(i) = -0.4d0
        END IF
      END DO
      !===Find Gauges
      CALL bar_gauges


    CASE DEFAULT
       WRITE(*,*) ' BUG in init'
       STOP
    END SELECT

    !===Initialization of un
    DO k = 1, inputs%syst_size
       un(k,:) = sol_anal(k,mesh%rr,inputs%time)
    END DO

    IF (inputs%type_test==13) THEN
       max_water_h = MINVAL(un(1,:))
    ELSE
       max_water_h = MAXVAL(un(1,:))
    END IF
    inputs%max_water_h = max_water_h

    !===Definition of htiny
    inputs%htiny=inputs%epsilon_htiny*max_water_h

    !===compute velocity
    un_over_h = compute_divide_by_h(un)
    velocity = un_over_h(2:2+k_dim-1,:)

    !===Define regul_h
    regul_h = 0.d0
    DO m = 1, mesh%me
       regul_h(mesh%jj(:,m)) = regul_h(mesh%jj(:,m)) + SUM(mesh%gauss%rj(:,m))
    END DO
    CALL ns_l1 (mesh, un(1,:), h_l1)
    surface = h_l1/max_water_h
    limit_h = inputs%epsilon_limit*max_water_h*(regul_h/surface)
    regul_h = inputs%epsilon_regul_h*max_water_h*(regul_h/surface)
    IF (max_water_h==0.d0) THEN
      WRITE(*,*) ' Forgot to initialize max_water_h', max_water_h
      STOP
    END IF
    write(*,*) ' limit_h/max_water_h' , MINVAL(limit_h)/max_water_h, MAXVAL(limit_h)/max_water_h, max_water_h
  END SUBROUTINE init

  SUBROUTINE seawall_gauges
    USE mesh_handling
    !USE mesh_interpolation
    IMPLICIT NONE
    INTEGER :: m, n
    ! these are gauges away from "dry land"
    ! these gauges measures change water_height + bath
    seawall_rr(1,1) = 5.9d0
    seawall_rr(1,2) = 7.6d0
    seawall_rr(1,3) = 9.64d0
    seawall_rr(1,4) = 10.462d0

    ! these are locations on "dry land"
    seawall_rr(1,5) = 10.732d0
    seawall_rr(1,6) = 11.005d0
    seawall_rr(1,7) = 11.12d0

    ! here we find the appropriate physical elements for gauge interpolation
    WRITE(*,*) '--------- seawall gauges stuff -------------'
    DO n = 1, 7
       write(*,*) 'n ', n
       DO m = 1, mesh%me
         IF (ABS(seawall_rr(1,n) - mesh%rr(1,m)) < 1.d-14 ) THEN
           !WRITE(*,*) 'im on a mesh node'
           !WRITE(*,*) m
           !seawall_m(n)=m
         ELSE  IF ((seawall_rr(1,n) .GE. mesh%rr(1,m)) .AND. (seawall_rr(1,n) .LE. mesh%rr(1,m+1))) THEN
          seawall_m(n)=m
          !WRITE(*,*) m
          !WRITE(*,*) seawall_rr(1,n)
          END IF
       END DO

       WRITE(*,*) mesh%rr(1,mesh%jj(:,seawall_m(n))), seawall_rr(1,n)
       WRITE(*,*) bath(mesh%jj(:,seawall_m(n))), SUM(bath(mesh%jj(:,seawall_m(n))) &
              * FE_interp_1d(mesh,seawall_m(n),seawall_rr(1,n)))
    END DO
    WRITE(*,*) seawall_m
    WRITE(*,*) '--------------------------------------------'
    !STOP
    ! create txt files here
    seawall_file(1)  ='WG1.txt'
    seawall_file(2)  ='WG3.txt'
    seawall_file(3)  ='WG10.txt'
    seawall_file(4)  ='WG22.txt'
    seawall_file(5)  ='WG28.txt'
    seawall_file(6)  ='WG37.txt'
    seawall_file(7)  ='WG40.txt'
  END SUBROUTINE seawall_gauges

  SUBROUTINE bar_gauges
    USE mesh_handling
    !USE mesh_interpolation
    IMPLICIT NONE
    INTEGER :: m, n
    ! these are gauges away from "dry land"
    ! these gauges measures change water_height + bath
    bar_rr(1,1) = 10.5d0
    bar_rr(1,2) = 12.5d0
    bar_rr(1,3) = 13.5d0
    bar_rr(1,4) = 14.5d0

    bar_rr(1,5) = 15.7d0
    bar_rr(1,6) = 17.3d0
    bar_rr(1,7) = 19.d0
    bar_rr(1,8) = 21.d0

    ! here we find the appropriate physical elements for gauge interpolation
    WRITE(*,*) '--------- bar gauges stuff -------------'
    DO n = 1, 8
       write(*,*) 'n ', n
       DO m = 1, mesh%me
         IF (ABS(bar_rr(1,n) - mesh%rr(1,m)) < 1.d-14 ) THEN
           !WRITE(*,*) 'im on a mesh node'
           !WRITE(*,*) m
           !seawall_m(n)=m
         ELSE  IF ((bar_rr(1,n) .GE. mesh%rr(1,m)) .AND. (bar_rr(1,n) .LE. mesh%rr(1,m+1))) THEN
          bar_m(n)=m
          !WRITE(*,*) m
          !WRITE(*,*) seawall_rr(1,n)
          END IF
       END DO

       WRITE(*,*) mesh%rr(1,mesh%jj(:,bar_m(n))), bar_rr(1,n)
       WRITE(*,*) bath(mesh%jj(:,bar_m(n))), SUM(bath(mesh%jj(:,bar_m(n))) &
              * FE_interp_1d(mesh,bar_m(n),bar_rr(1,n)))
    END DO
    WRITE(*,*) bar_m
    WRITE(*,*) '--------------------------------------------'
    !STOP
    ! create txt files here
    bar_file(1)  ='WG4.txt'
    bar_file(2)  ='WG5.txt'
    bar_file(3)  ='WG6.txt'
    bar_file(4)  ='WG7.txt'
    bar_file(5)  ='WG8.txt'
    bar_file(6)  ='WG9.txt'
    bar_file(7)  ='WG10.txt'
    bar_file(8)  ='WG11.txt'
  END SUBROUTINE bar_gauges

  FUNCTION sol_anal(k,rr,t) RESULT(vv) !CAREFUL HERE, rr comes from the calling subroutine
    USE mesh_handling
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k
    REAL(KIND=8), DIMENSION(:,:),  INTENT(IN) :: rr
    REAL(KIND=8),                  INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2)) :: vv, aux, x_coord, bath, gauss_height
    INTEGER :: i
    REAL(KIND=8) :: x, x0, speed, q0, hL, b, d, x1, x2, x3, &
         theta, Rcard, Scard, Tcard, Qcard, Dcard, tpio3, fpio3, a, omega, eta, h0, bernoulli, &
         xshock, h_pre_shock, h_post_shock, bath_shock, bathi, Ber_pre, Ber_post, &
         alpha, beta, chi, vel, xs, hcone, htop, radius, rcone, scone, bx, dd, h1, h2, &
         D_wave, htilde, slope, SS, z, c, sechSqd, hTildePrime, k_wavenumber, L, &
         cSW, cGN, q, cBer, term1,term2,term3, Tperiod
    !===Malpasset
    REAL(KIND=8), DIMENSION(SIZE(rr,2)) :: h_old
    !===
    SELECT CASE(inputs%type_test)

    CASE(13) !===Solitary wave
       h1 = 10.d0
       h2 = 11.d0
       x0 = 200
       dd = SQRT(inputs%gravity*h2)
       eta = SQRT(3.d0*(h2-h1)/(h2*h1**2))/2.d0
       aux = eta*(rr(1,:)-x0-dd*inputs%time)
       vv = h1 + (h2-h1)*(1.d0/COSH(aux)**2) !water height
       SELECT CASE(k)
       CASE(1)
          vv = h1 + (h2-h1)*(1.d0/COSH(aux)**2) !water height
       CASE(2)
          !IF (inputs%time<1.d-10) THEN
          vv = dd*(1-h1/vv)*vv !u_0*h_0
          !ELSE
          !   vv = 0.d0
          !END IF
       CASE(3)
          vv = (vv)**2  !h_0*eta_0
       CASE(4)
          aux = -eta*(h2-h1)*2*(1.d0/COSH(aux)**2)*tanh(aux) !dh/dx
          vv = (dd*(1-h1/vv)-dd)*aux*vv !(u0-d)dh0/dx * h0
       END SELECT

    CASE(14) !===Run up
       inputs%gravity = 9.81d0
       h0 = 1.d0
       a = 0.28d0 ! amplitude
       slope = 1.d0 / 19.85d0
       k_wavenumber = SQRT(3.d0 * a/(4.d0 * h0**3)) ! wavenumber
       z = SQRT(3.d0 * a * h0) / (2.d0 * h0 * SQRT(h0 * (1.d0 + a)))
       L = 2.d0 / k_wavenumber * ACOSH(SQRT(1.d0 / 0.05d0)) ! wavelength of solitary wave
       c = SQRT(inputs%gravity * (1.d0 + a) * h0)
       x0 = - h0/slope - L/2.d0  ! initial location of solitary wave

       SELECT CASE(k)
       CASE(1) ! h water height
          DO i = 1, SIZE(rr,2)
             bathi = MAX((rr(1,i)) * slope, -h0)
             sechSqd = (1.0d0/COSH( z*(rr(1,i)-x0-c*t)))**2.0d0
             htilde= a * h0 * sechSqd
             vv(i) = MAX(htilde - bathi,0.d0)
          END DO

       CASE(2) ! u*h component, u = c htilde/ (htilde + h0)
          DO i = 1, SIZE(rr,2)
             bathi = MAX((rr(1,i) ) * slope, -h0)
             sechSqd = (1.0d0/COSH( z*(rr(1,i)-x0-c*t)))**2.0d0
             htilde= a * h0 * sechSqd ! this is exact solitary wave
             vv(i) = MAX(htilde-bathi,0.d0)
             vv(i) =  vv(i) * c * htilde / (h0 + htilde)
          END DO

       CASE(3) ! eta*h component
          IF (t.LE.1.d-14) THEN
             DO i = 1, SIZE(rr,2)
                bathi = MAX((rr(1,i)) * slope, -h0)
                sechSqd = (1.0d0/COSH( z*(rr(1,i)-x0-c*t)))**2.0d0
                htilde= a * h0 * sechSqd
                vv(i) = MAX(htilde - bathi,0.d0)
                vv(i) = vv(i) * vv(i)
             END DO
          END IF
       CASE(4) ! w*h component of flow rate, which is -waterHeight^2 * div(velocity)
          IF (t.LE.1.d-14) THEN
             DO i = 1, SIZE(rr,2)
                bathi = MAX((rr(1,i) ) * slope, -h0)
                sechSqd = (1.0d0/COSH( z*(rr(1,i)-x0-c*t)))**2.0d0
                htilde= a * h0 * sechSqd ! this is exact solution
                hTildePrime = -2.d0 * z * htilde * TANH(z*(rr(1,i)-x0-c*t))
                vv(i) = MAX(htilde - bathi,0.d0)
                ! this is -waterHeight^2 * div(velocity)
                vv(i) = -vv(i)**2 * (c * h0 * hTildePrime /(h0 + htilde)**2)
             END DO
          END IF
       END SELECT

    CASE(15) ! mGN steady state
       ! initial constants go here
       a = 1.5d0
       L = 1.d0
       h0 = 2.d0
       q = 1.5d0
       cBer = h0 + q**2 / (2.d0 * inputs%gravity * h0**2)

       ! gauss_height = h0 + a * EXP(-(x_coord/L)**2)
       ! term1 = a * (L-x_coord)*(L+x_coord)
       ! term2 = EXP((x_coord/L)**2) * h0 * (L**2 - 2.d0 * x_coord**2)
       ! bath = cBer - gauss_height - q**2 / (2.d0 * inputs%gravity * gauss_height**2) &
       !       + 2.d0 * a * q**2 / (3.d0 * inputs%gravity * L**4) *  &
       !       1.d0 / (a + h0*EXP((x_coord/L)**2))**2 * (term1 + term2)

       ! for shallow water cardano's formula
       x0 = 0.d0
       tpio3=2*ACOS(-1.d0)/3
       fpio3=2*tpio3
       d = q**2 / (2.d0 * inputs%gravity)

       ! for bath
       x_coord = rr(1,:)

       k_wavenumber = SQRT(3.d0 * a/(4.d0 * h0**3)) ! wavenumber
       z = SQRT(3.d0 * a * h0) / (2.d0 * h0 * SQRT(h0 * (1.d0 + a)))
       L = 2.d0 / k_wavenumber * ACOSH(SQRT(1.d0 / 0.05d0)) ! wavelength of solitary wave
       sechSqd = (1.0d0/COSH( z*x ))**2.0d0

       term1 = (1.d0 + a)*(2.d0*inputs%gravity*(cber - h0)*h0**2 - q**2)
       term2 = 2*a*((1 + a)*inputs%gravity*h0**3 - q**2)
       term3 = 2*(1 + a)*inputs%gravity*h0**2


       SELECT CASE(k)
       CASE(1) ! h water height
         DO i = 1, SIZE(rr,2)
           sechSqd = (1.0d0/COSH( z*rr(1,i)))**2.0d0
           bathi = (term1 - term2 * sechSqd)/term3
           htilde= h0 + a * h0 * sechSqd
           vv(i) = htilde
        END DO
         ! ELSE
           ! DO i = 1, SIZE(rr,2)
           !   vv(i) = gauss_height(i)
           !   IF (inputs%lambda_bar < 1.d-10) THEN
           !     b = bath(i) - cBer
           !     d = q**2/(2*inputs%gravity)
           !     Rcard = (-27*d-2*b**3)/(54)
           !     Qcard = -b**2/9
           !     Dcard = Qcard**3+Rcard**2
           !     IF (Dcard.GE.0.d0) THEN
           !        !write(*,*) Rcard**2,Qcard**3,Qcard**3+Rcard**2
           !        !write(*,*) Rcard + sqrt(Dcard)
           !        Scard=-(-SQRT(Dcard)-Rcard)**(1.d0/3)
           !        Tcard=-(SQRT(Dcard)-Rcard)**(1.d0/3)
           !        x1 = Scard+Tcard-b/3.d0
           !        x2 = -(Scard+Tcard)/2.d0 - b/3.d0
           !        !write(*,*) x1**3+b*x1**2+d
           !        vv(i) = x2
           !        write(*,*) 'Dcard', Dcard, bath(i), x2
           !     ELSE
           !        !write(*,*) 'ratio', Rcard/SQRT(-Qcard**3), 'Dcard', Dcard, '-Qcard', -Qcard
           !        theta = ACOS(Rcard/SQRT(-Qcard**3))
           !        x1 = 2*SQRT(-Qcard)*cos(theta/3) - b/3
           !        !write(*,*) x1**3+b*x1**2+d
           !        x2 = 2*SQRT(-Qcard)*cos(theta/3+tpio3) - b/3
           !        x3 = 2*SQRT(-Qcard)*cos(theta/3+fpio3) - b/3
           !        !write (*,*) 'x1,x2,x3',x1, x2, x3
           !        !IF (rr(1,i)<x0) THEN
           !           vv(i) = MAX(x1,x2,x3)
           !        !ELSE
           !        !    vv(i) = MIN(x1,x3)
           !        ! END IF
           !        !write(*,*) 'vv(i)', vv(i)
           !     END IF
           !   END IF
           ! END DO
         ! END IF

       CASE(2) ! u*h component, hu = q = 6
         IF (t.LE.1.d-14) THEN
           DO i=1, SIZE(rr,2)
             vv(i) = q
           END DO
         ELSE
           DO i = 1, SIZE(rr,2)
             vv(i) = q
           END DO
         END IF

      CASE(3) ! eta*h component
        !DO i = 1, SIZE(rr,2)
        DO i = 1, SIZE(rr,2)
          sechSqd = (1.0d0/COSH( z*rr(1,i)))**2.0d0
          bathi  = (term1 - term2 * sechSqd)/term3
          htilde = h0 + a * h0 * sechSqd
          vv(i)  = htilde**2
        END DO
        !   vv(i) = gauss_height(i)**2
        !   IF (inputs%lambda_bar < 1.d-10) THEN
        !     b = bath(i) - cBer
        !     d = q**2/(2*inputs%gravity)
        !     Rcard = (-27*d-2*b**3)/(54)
        !     Qcard = -b**2/9
        !     Dcard = Qcard**3+Rcard**2
        !     IF (Dcard.GE.0.d0) THEN
        !        !write(*,*) Rcard**2,Qcard**3,Qcard**3+Rcard**2
        !        !write(*,*) Rcard + sqrt(Dcard)
        !        Scard=-(-SQRT(Dcard)-Rcard)**(1.d0/3)
        !        Tcard=-(SQRT(Dcard)-Rcard)**(1.d0/3)
        !        x1 = Scard+Tcard-b/3.d0
        !        x2 = -(Scard+Tcard)/2.d0 - b/3.d0
        !        !write(*,*) x1**3+b*x1**2+d
        !        vv(i) = x2
        !        write(*,*) 'Dcard', Dcard, bath(i), x2
        !     ELSE
        !        !write(*,*) 'ratio', Rcard/SQRT(-Qcard**3), 'Dcard', Dcard, '-Qcard', -Qcard
        !        theta = ACOS(Rcard/SQRT(-Qcard**3))
        !        x1 = 2*SQRT(-Qcard)*cos(theta/3) - b/3
        !        !write(*,*) x1**3+b*x1**2+d
        !        x2 = 2*SQRT(-Qcard)*cos(theta/3+tpio3) - b/3
        !        x3 = 2*SQRT(-Qcard)*cos(theta/3+fpio3) - b/3
        !        !write (*,*) 'x1,x2,x3',x1, x2, x3
        !        !IF (rr(1,i)<x0) THEN
        !           vv(i) = MAX(x1,x2,x3)**2
        !        !ELSE
        !        !    vv(i) = MIN(x1,x3)
        !        ! END IF
        !        !write(*,*) 'vv(i)', vv(i)
        !     END IF
        !   END IF
        ! END DO
       CASE(4) ! w*h component of flow rate, which is -waterHeight^2 * div(velocity)
         DO i = 1, SIZE(rr,2)
           vv(i) = 0.d0
         END DO
     END SELECT

    CASE(16) !===Seawall
      ! initial constants go here
      inputs%gravity = 9.81d0
      h0 = 0.2d0
      a = 0.35d0 ! amplitude
      k_wavenumber = SQRT(3.d0 * a/(4.d0 * h0**3)) ! wavenumber
      z = SQRT(3.d0 * a * h0) / (2.d0 * h0 * SQRT(h0 * (1.d0 + a)))
      L = 2.d0 / k_wavenumber * ACOSH(SQRT(1.d0 / 0.05d0)) ! wavelength of solitary wave
      c = SQRT(inputs%gravity * (1.d0 + a) * h0) ! wave speed
      x0 = 5.9d0 ! initial location of solitary wave

      ! define bathymetry here because it's easier
      ! need to shift the bathymetry 3 meters, see basilisk website
      x_coord = rr(1,:) + 3.d0
      DO i = 1, SIZE(rr,2)
        IF (x_coord(i) .LE. 10.d0) THEN
          bath(i) = 0.d0
        ELSE IF (x_coord(i) .GE. 13.6d0 .AND. x_coord(i) .LE. 13.9d0) THEN
          bath(i) = 3.6d0/20.d0 + (x_coord(i)-13.6d0)*0.076d0/(13.9d0 - 13.6d0)
        ELSE IF (x_coord(i) .GE. 13.9d0 .AND. x_coord(i) .LE. 13.948d0) THEN
          bath(i) = 3.6d0/20.d0 + 0.076d0
        ELSE IF (x_coord(i) .GE. 13.948d0 .AND. x_coord(i) .LE. 14.045d0) THEN
          bath(i) = 3.6d0/20.d0 + 0.076d0 - (x_coord(i)-13.948d0)*(0.076d0-0.022d0)/(14.045d0 - 13.948d0)
        ELSE
          bath(i) = (x_coord(i) - 10.d0)/20.d0
        END IF
      END DO
      bath = bath - h0 ! shift bath down by reference height

      SELECT CASE(k)
      CASE(1) ! h water height
          DO i = 1, SIZE(rr,2)
            sechSqd = (1.0d0/COSH( z*(rr(1,i)-x0-c*t)))**2.0d0
            htilde = a * h0 * sechSqd
            vv(i) = MAX(htilde - bath(i),0.d0)
          END DO
      CASE(2) ! h*u component, u = c htilde/ (htilde + h0)
        DO i = 1, SIZE(rr,2)
          sechSqd = (1.0d0/COSH( z*(rr(1,i)-x0-c*t)))**2.0d0
          htilde= a * h0 * sechSqd ! this is exact solitary wave
          vv(i) = MAX(htilde-bath(i),0.d0)
          vv(i) =  vv(i) * c * htilde / (h0 + htilde)
        END DO
      CASE(3) ! h*eta component, eta = h
           DO i = 1, SIZE(rr,2)
             sechSqd = (1.0d0/COSH( z*(rr(1,i)-x0-c*t)))**2.0d0
             htilde= a * h0 * sechSqd
             vv(i) = MAX(htilde - bath(i),0.d0)
             vv(i) = vv(i) * vv(i)
           END DO
      CASE(4) ! h*w component of flow rate, which is -waterHeight^2 * div(velocity)
           DO i = 1, SIZE(rr,2)
             sechSqd = (1.0d0/COSH( z*(rr(1,i)-x0-c*t)))**2.0d0
             htilde= a * h0 * sechSqd ! this is exact solution
             hTildePrime = -2.d0 * z * htilde * TANH(z*(rr(1,i)-x0-c*t))
             vv(i) = MAX(htilde - bath(i),0.d0)
             ! this is -waterHeight^2 * div(velocity)
             vv(i) = -vv(i)**2 * (c * h0 * hTildePrime /(h0 + htilde)**2)
           END DO
      END SELECT

    CASE(17) !===Period waves over bar
      ! initial constants go here
      inputs%gravity = 9.81d0
      a = 0.051d0
      Tperiod = 2.02d0

      ! define bathymetry here because it's easier
      x_coord = rr(1,:)
      DO i = 1, SIZE(rr,2)
        IF (x_coord(i) .LE. 6.d0) THEN
          bath(i) = -0.4d0
        ELSE IF (x_coord(i) .GE. 6.d0 .AND. x_coord(i) .LE. 12.d0) THEN
          bath(i) = -0.4d0 + (x_coord(i)-6.d0)*1.d0/20.d0
        ELSE IF (x_coord(i) .GE. 12.d0 .AND. x_coord(i) .LE. 14.d0) THEN
          bath(i) = -0.1d0
        ELSE IF (x_coord(i) .GE. 14.d0 .AND. x_coord(i) .LE. 17.d0) THEN
          bath(i) = -0.1d0 - (x_coord(i)-14.d0)*1.d0/10.d0
        ELSE
          bath(i) = -0.4d0
        END IF
      END DO

      SELECT CASE(k)
      CASE(1) ! h water height, defined to be 0 for BC at outlet
          DO i = 1, SIZE(rr,2)
            htilde = 0.d0
            vv(i) = MAX(htilde - bath(i),0.d0)
          END DO

      CASE(2) ! u*h component, here we have q (dot) n = a*sin()
          DO i = 1, SIZE(rr,2)
            IF (rr(1,i) .LE. 50.d0) THEN
              htilde = 0.4d0
              vv(i) = -htilde * a*SIN(2.d0 * pi * t/Tperiod)
            ELSE
              vv(i) = 0.d0
            END IF
          END DO

      CASE(3) ! eta*h component, initalized with heta = h^2
      !  IF (t .LE. 1d-14) THEN
          DO i = 1, SIZE(rr,2)
            htilde = 0.d0
            vv(i) = MAX(htilde - bath(i),0.d0)
            vv(i) = vv(i)*vv(i)
          END DO

      CASE(4) ! w*h component of flow rate, which is -waterHeight^2 * div(velocity)
           DO i = 1, SIZE(rr,2)
             vv(i) = 0.d0
           END DO
      END SELECT

    CASE DEFAULT
       WRITE(*,*) ' BUG in sol_anal'
       STOP
    END SELECT

  END FUNCTION sol_anal


  !SUBROUTINE compute_lambda_vacc(ul,ur,vell,velr,nij,lambda)
  SUBROUTINE compute_lambda_vacc(ul,ur,vell,velr,hhl,hhr,nij,lambda,if_dt)
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), INTENT(OUT) :: lambda
    REAL(KIND=8), DIMENSION(k_dim), INTENT(IN)  :: nij, vell, velr
    REAL(KIND=8), DIMENSION(inputs%syst_size), INTENT(IN)   :: ur, ul
    REAL(KIND=8), INTENT(IN)  :: hhl, hhr
    LOGICAL, INTENT(IN) :: if_dt
    REAL(KIND=8) :: fl, fr, ht, vl, vr, hl, hr, lbdl, lbdr, lbd_dry, sql, sqr, etal, etar
    REAL(KIND=8) :: fh, a, c, Delta, x0, hmin, hmax, vmin, vmax, sqrmin, sqrmax
    REAL(KIND=8) :: ovhl, ovhr
    hl = max(ul(1),1.d-16*max_water_h)
    hr = max(ur(1),1.d-16*max_water_h)
    ovhl=1.d0/hl
    ovhr=1.d0/hr
    vl =  SUM(vell*nij)
    vr =  SUM(velr*nij)
    sql = SQRT(inputs%gravity*ABS(hl))
    sqr = SQRT(inputs%gravity*ABS(hr))

    SELECT CASE(inputs%type_test)
    CASE(13,14,15,16,17)

       x0=(2.d0*SQRT(2.d0)-1.d0)**2
       IF(hl.LE.hr) THEN
          hmin = hl
          vmin = vl
          hmax = hr
          vmax = vr
       ELSE
          hmin = hr
          vmin = vr
          hmax = hl
          vmax = vl
       END IF

       !===Case 1
       fh = phi(x0*hmin,vl,hl,vr,hr)
       IF (0.d0.LE.fh) THEN
          ht = min(x0*hmin,(MAX(vl-vr+2*sql+2*sqr,0.d0))**2/(16*inputs%gravity))
          lbdl = vl - sql*SQRT((1+max((ht-hl)*ovhl/2,0.d0))*(1+max((ht-hl)*ovhl,0.d0)))
          lbdr = vr + sqr*SQRT((1+max((ht-hr)*ovhr/2,0.d0))*(1+max((ht-hr)*ovhr,0.d0)))
          !lambda = MAX(ABS(lbdl),ABS(lbdr))
          !RETURN
          GO TO 100
       END IF
       !===Case 2
       fh = phi(x0*hmax,vl,hl,vr,hr)
       IF (0.d0.LE.fh) THEN
          sqrmin = SQRT(hmin)
          sqrmax = SQRT(hmax)
          a = 1/(2.d0*SQRT(2.d0))
          c = -hmin*a -sqrmin*sqrmax + sqrmin*(vr-vl)/(2*sqrt(inputs%gravity))
          Delta = hmin-4*a*c
          IF (Delta<0.d0) THEN
             WRITE(*,*) ' BUG in compute_lambda_vacc'
             STOP
          END IF
          ht = min(x0*hmax,((-sqrmin+sqrt(Delta))/(2*a))**2)
          lbdl = vl - sql*SQRT((1+max((ht-hl)*ovhl/2,0.d0))*(1+max((ht-hl)*ovhl,0.d0)))
          lbdr = vr + sqr*SQRT((1+max((ht-hr)*ovhr/2,0.d0))*(1+max((ht-hr)*ovhr,0.d0)))
          !lambda = MAX(ABS(lbdl),ABS(lbdr))
          !RETURN
          GO TO 100
       END IF
       !===Case 3
       sqrmin = SQRT(hmin)
       sqrmax = SQRT(hmax)
       ht = sqrmin*sqrmax*(1+sqrt(2/inputs%gravity)*(vl-vr)/(sqrmin+sqrmax))
       !lbdl = vl - sql*SQRT((1+max((ht-hl)*ovhl/2,0.d0))*(1+max((ht-hl)*ovhl,0.d0)))
       !lbdr = vr + sqr*SQRT((1+max((ht-hr)*ovhr/2,0.d0))*(1+max((ht-hr)*ovhr,0.d0)))
       lbdl = vl - sql*SQRT((1+max((ht-hl)/(2*hl),0.d0))*(1+max((ht-hl)/hl,0.d0)))
       lbdr = vr + sqr*SQRT((1+max((ht-hr)/(2*hr),0.d0))*(1+max((ht-hr)/hr,0.d0)))

100    CONTINUE

       !===Correction for fake SGN, ERIC FIX
       etal = ul(inputs%syst_size-1)/hl
       etar = ur(inputs%syst_size-1)/hr
       IF (if_dt) THEN
         a = inputs%lambda_bar*inputs%gravity/hhl
         IF (etal < hl) THEN
           lbdl = min(lbdl,vl-sqrt(inputs%gravity*hl + (a/3)*(6.d0*hl+ 12.d0*(hl-etal))))
         ELSE
           lbdl = min(lbdl,vl-sqrt(inputs%gravity*hl + (a/3)*(6.d0*hl)))
         END IF
         IF (etar < hr) THEN
           lbdr = MAX(lbdr,vr-sqrt(inputs%gravity*hr + (a/3)*(6.d0*hr+ 12.d0*(hr-etar))))
         ELSE
           lbdr = MAX(lbdr,vr-sqrt(inputs%gravity*hr + (a/3)*(6.d0*hr)))
         END IF
       END IF
       ! IF (if_dt) THEN
       !   a = inputs%lambda_bar*inputs%gravity*(ul(inputs%syst_size-1)/hl)**2/hhl
       !   lbdl = min(lbdl,vl-sqrt(inputs%gravity*hl + (a/3)*ul(inputs%syst_size-1)**2/hl**4))
       !   a = inputs%lambda_bar*inputs%gravity*(ur(inputs%syst_size-1)/hr)**2/hhr
       !   lbdr = max(lbdl,vr+sqrt(inputs%gravity*hr + (a/3)*ur(inputs%syst_size-1)**2/hr**4))
       ! END IF
       !===end correction for fake SGN

       lambda = MAX(ABS(lbdl),ABS(lbdr))

    CASE DEFAULT
       WRITE(*,*) 'BUG in compute_lambda'
       STOP
    END SELECT
    RETURN
  END SUBROUTINE compute_lambda_vacc


  FUNCTION phi(h,ul,hl,ur,hr) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: h, ul, hl, ur, hr
    REAL(KIND=8)             :: vv, fl, fr, ovh, ovhl, ovhr
    IF (h>hl) THEN
       fl = (h-hl)*SQRT((inputs%gravity/2)*(h+hl)/(h*hl))
    ELSE
       fl = 2*(SQRT(inputs%gravity*h)-SQRT(inputs%gravity*hl))
    END IF
    IF (h>hr) THEN
       fr = (h-hr)*SQRT((inputs%gravity/2)*(h+hr)/(h*hr))
    ELSE
       fr = 2*(SQRT(inputs%gravity*h)-SQRT(inputs%gravity*hr))
    END IF
    vv = fl + fr + ur - ul
  END FUNCTION phi

  FUNCTION FE_interp_1d(mesh, m,r) RESULT(ff)
  USE def_type_mesh
  IMPLICIT NONE
  TYPE(mesh_type)                        :: mesh
  INTEGER,        	        INTENT(IN) :: m
  REAL(KIND=8), INTENT(IN) :: r

  REAL(KIND=8), DIMENSION(2):: ff

  REAL(KIND=8) :: r1, r2

  REAL(KIND=8) :: x

  r1 = mesh%rr(1,m)
  r2 = mesh%rr(1,m+1)

  x = (r - r1) / (r2 -r1)
  ff(1) = 1-x
  ff(2) = x

  END FUNCTION FE_interp_1d

END MODULE boundary_conditions
