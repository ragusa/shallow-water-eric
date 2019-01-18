MODULE input_data
  IMPLICIT NONE
  PUBLIC :: read_my_data
  TYPE my_data
     CHARACTER(len=200)             :: directory
     CHARACTER(len=200)             :: file_name
     CHARACTER(len=200)             :: viscous_type
     INTEGER                        :: nb_dom
     INTEGER, DIMENSION(:), POINTER :: list_dom
     INTEGER                        :: type_fe
     REAL(KIND=8)                   :: Tfinal
     REAL(KIND=8)                   :: CFL
     CHARACTER(LEN=15)              :: viscosity_type
     LOGICAL                        :: if_lumped
     LOGICAL                        :: if_alpha_limit
     LOGICAL                        :: if_FGN, if_FGN_update, want_movie
     CHARACTER(LEN=3)               :: limiter_type
     INTEGER                        :: type_test
     REAL(KIND=8)                   :: dt, time
     INTEGER                        :: h_nb_Dir_bdy, ux_nb_Dir_bdy, uy_nb_Dir_bdy
     INTEGER, DIMENSION(:), POINTER :: h_Dir_list, ux_Dir_list, uy_Dir_list
     REAL(KIND=8)                   :: gravity
     REAL(KIND=8)                   :: mannings, slope, discharge
     REAL(KIND=8)                   :: eta
     REAL(KIND=8)                   :: lambda_bar
     INTEGER                        :: nb_udotn_zero
     INTEGER                        :: nb_Dir_bdy
     INTEGER, DIMENSION(:), POINTER :: Dir_list
     INTEGER, DIMENSION(:), POINTER :: udotn_zero_list
     INTEGER                        :: syst_size
     INTEGER                        :: nb_frame
     REAL(KIND=8)                   :: htiny
     REAL(KIND=8)                   :: epsilon_regul_h, epsilon_limit
     REAL(KIND=8)                   :: epsilon_htiny, epsilon_pminus
     REAL(KIND=8)                   :: max_water_h
  END type my_data
  TYPE(my_data), PUBLIC  :: inputs
  PRIVATE
CONTAINS
  SUBROUTINE read_my_data(data_fichier)
    USE character_strings
    USE space_dim
    IMPLICIT NONE
    INTEGER, PARAMETER           :: in_unit=21
    CHARACTER(len=*), INTENT(IN) :: data_fichier
    LOGICAL :: okay
    inputs%epsilon_pminus   =-1.d-10 !Fct
    inputs%epsilon_htiny    = 1.d-10 !htiny
    inputs%epsilon_limit    = 1.d-10  !limiter + limit alpha
    inputs%epsilon_regul_h  = 1.d-10 !Velocity
    OPEN(UNIT = in_unit, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(in_unit, "===Name of directory for mesh file===")
    READ (in_unit,*) inputs%directory
    CALL read_until(in_unit, "===Name of mesh file===")
    READ (in_unit,*) inputs%file_name
    CALL read_until(in_unit, '===Number of subdomains in the mesh===')
    READ(in_unit,*) inputs%nb_dom
    ALLOCATE(inputs%list_dom(inputs%nb_dom))
    CALL read_until(21, '===List of subdomain in the mesh===')
    READ(in_unit,*) inputs%list_dom
    CALL read_until(21, '===Type of finite element===')
    READ(in_unit,*) inputs%type_fe
    CALL read_until(in_unit, "===Final time===")
    READ (in_unit,*) inputs%Tfinal
    CALL read_until(in_unit, "===CFL number===")
    READ (in_unit,*) inputs%CFL
    CALL read_until(in_unit, "===Viscosity type (galerkin, viscous, entropy_visc, fct)===")
    READ (in_unit,*) inputs%viscosity_type
    CALL read_until(in_unit, "===Type of first-order viscosity===")
    READ (in_unit,*) inputs%viscous_type
    CALL read_until(in_unit, "===Alpha limitation? (True/False)===")
    READ (in_unit,*) inputs%if_alpha_limit
    CALL read_until(in_unit, "===Lumping the mass matrix? (true/false)===")
    READ (in_unit,*) inputs%if_lumped
    CALL read_until(in_unit, "===Test case number===")
    READ (in_unit,*) inputs%type_test
    CALL find_string(in_unit, "===How many boundaries for u.n=0?===",okay)
    IF (okay) THEN
       READ (in_unit,*)  inputs%nb_udotn_zero
       IF (inputs%nb_udotn_zero>0) THEN
          CALL read_until(in_unit, "===List of boundarie for u.n=0?===")
          ALLOCATE(inputs%udotn_zero_list(inputs%nb_udotn_zero))
          READ (in_unit,*) inputs%udotn_zero_list
       ELSE
          inputs%nb_udotn_zero = 0
          ALLOCATE(inputs%udotn_zero_list(inputs%nb_udotn_zero))
       END IF
    ELSE
       inputs%nb_udotn_zero = 0
       ALLOCATE(inputs%udotn_zero_list(inputs%nb_udotn_zero))
    END IF
    CALL read_until(in_unit, "===How many Dirichlet boundaries for h?===")
    READ (in_unit,*)  inputs%h_nb_Dir_bdy
    CALL read_until(in_unit, "===List of Dirichlet boundaries for h?===")
    ALLOCATE(inputs%h_Dir_list(inputs%h_nb_Dir_bdy))
    READ (in_unit,*) inputs%h_Dir_list
    CALL read_until(in_unit, "===How many Dirichlet boundaries for ux?===")
    READ (in_unit,*)  inputs%ux_nb_Dir_bdy
    CALL read_until(in_unit, "===List of Dirichlet boundaries for ux?===")
    ALLOCATE(inputs%ux_Dir_list(inputs%ux_nb_Dir_bdy))
    READ (in_unit,*) inputs%ux_Dir_list
    IF (k_dim==2) THEN
       CALL read_until(in_unit, "===How many Dirichlet boundaries for uy?===")
       READ (in_unit,*)  inputs%uy_nb_Dir_bdy
       CALL read_until(in_unit, "===List of Dirichlet boundaries for uy?===")
       ALLOCATE(inputs%uy_Dir_list(inputs%uy_nb_Dir_bdy))
       READ (in_unit,*) inputs%uy_Dir_list
    END IF
    CALL find_string(in_unit, "===How many vtk frames?===",okay)
    IF (okay) THEN
       READ (in_unit,*)  inputs%nb_frame
    ELSE
       inputs%nb_frame = 0
    END IF

    CALL find_string(in_unit, "===FGN?===(.t.,.f.),okay===",okay)
    IF (okay) THEN
       READ (in_unit,*)  inputs%if_FGN
       CALL read_until(in_unit, "===Lambda_bar?===")
       READ (in_unit,*) inputs%lambda_bar
    ELSE
       inputs%if_FGN = .FALSE.
       inputs%lambda_bar = 0.d0
    END IF
    CALL find_string(in_unit, "===run old FGN update?===(.t.,.f.),okay===",okay)
    IF (okay) THEN
       READ (in_unit,*)  inputs%if_FGN_update
    ELSE
       inputs%if_FGN_update = .FALSE.
    END IF

    CALL find_string(in_unit, "===Want movie?===(.t.,.f.),okay===",okay)
    IF (okay) THEN
       READ (in_unit,*)  inputs%want_movie
    ELSE
       inputs%want_movie= .FALSE.
    END IF


    SELECT CASE(inputs%type_test)
    CASE(9,12,14,16)
       CALL read_until(in_unit, "===Mannings coefficient===")
       READ (in_unit,*) inputs%mannings
    CASE(10)
       CALL read_until(in_unit, "===Mannings coefficient===")
       READ (in_unit,*) inputs%mannings
       CALL read_until(in_unit, "===Slope===")
       READ (in_unit,*) inputs%slope
       CALL read_until(in_unit, "===Discharge===")
       READ (in_unit,*) inputs%discharge
    CASE(11)
       CALL read_until(in_unit, "===Mannings coefficient===")
       READ (in_unit,*) inputs%mannings
       CALL read_until(in_unit, "===Discharge===")
       READ (in_unit,*) inputs%discharge
    CASE DEFAULT
       inputs%mannings = 0.d0
    END SELECT

 CLOSE(in_unit)
END SUBROUTINE read_my_data
END MODULE input_data
