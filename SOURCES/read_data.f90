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
     LOGICAL                        :: if_fct_limit
     INTEGER                        :: type_test
     REAL(KIND=8)                   :: dt, time
     INTEGER                        :: h_nb_Dir_bdy, ux_nb_Dir_bdy, uy_nb_Dir_bdy
     INTEGER, DIMENSION(:), POINTER :: h_Dir_list, ux_Dir_list, uy_Dir_list
     REAL(KIND=8)                   :: gravity
     REAL(KIND=8)                   :: mannings
     REAL(KIND=8)                   :: eta
     INTEGER                        :: nb_Dir_bdy
     INTEGER, DIMENSION(:), POINTER :: Dir_list
     INTEGER                        :: syst_size
     REAL(KIND=8)                   :: htiny
     REAL(KIND=8)                   :: epsilon, max_water_h,epsilon_htiny
     REAL(KIND=8)                   :: lambdaSGN
     REAL(KIND=8)                   :: x1, x2, y1, y2
     REAL(KIND=8)                   :: refinement
  END type my_data
  TYPE(my_data), PUBLIC  :: inputs
  PRIVATE
CONTAINS
  SUBROUTINE read_my_data(data_fichier)
    USE character_strings
    IMPLICIT NONE
    INTEGER, PARAMETER           :: in_unit=21
    CHARACTER(len=*), INTENT(IN) :: data_fichier
    inputs%epsilon = 1.d-13
    inputs%epsilon_htiny = 1.d-10
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
    CALL read_until(in_unit, "===How many Dirichlet boundaries for uy?===")
    READ (in_unit,*)  inputs%uy_nb_Dir_bdy
    CALL read_until(in_unit, "===List of Dirichlet boundaries for uy?===")
    ALLOCATE(inputs%uy_Dir_list(inputs%uy_nb_Dir_bdy))
    READ (in_unit,*) inputs%uy_Dir_list

    SELECT CASE(inputs%type_test)
    CASE(9,10,14,15,16)
       CALL read_until(in_unit, "===Mannings coefficient===")
       READ (in_unit,*) inputs%mannings
    CASE DEFAULT
       inputs%mannings = 0.d0
    END SELECT

    SELECT CASE(inputs%type_test) ! for hyperbolic SGN model
    CASE(13,14,15,16,17,18,19)
       CALL read_until(in_unit, "===Lambda for SGN model===")
       READ (in_unit,*) inputs%lambdaSGN
    END SELECT

 CLOSE(in_unit)
END SUBROUTINE read_my_data
END MODULE input_data
