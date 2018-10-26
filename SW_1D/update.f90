MODULE update
  USE mesh_handling
  USE matrix_type
  USE space_dim
  USE input_data
  PUBLIC:: construct_matrices, euler, compute_dij
  TYPE(matrice_bloc), PUBLIC                  :: dij, betaij
  REAL(KIND=8), DIMENSION(:), POINTER, PUBLIC :: lumped
  REAL(KIND=8), DIMENSION(:), POINTER, PUBLIC :: mesh_size
  INTEGER, DIMENSION(:), POINTER, PUBLIC      :: diag
  REAL(KIND=8), DIMENSION(:), POINTER, PUBLIC :: fix_roundoff
  TYPE(matrice_bloc), PUBLIC                  :: mass, pre_mass, mc_minus_ml
  PRIVATE
  TYPE(matrice_bloc), DIMENSION(k_dim)  :: cij
  TYPE(matrice_bloc), DIMENSION(k_dim+3):: fctmat
  TYPE(matrice_bloc)                    :: muij, resij, dijL, muijL, lij
  INTEGER                               :: isolve, isolve_m
  LOGICAL                               :: if_roundoff_fix=.FALSE.
  INTEGER                               :: max_nb_pt_stencil
CONTAINS

  SUBROUTINE construct_matrices
    USE st_matrix
    USE mesh_handling
    USE fem_s_M
    IMPLICIT NONE
    INTEGER :: m, p, ni, nj, i, j, d
    REAL(KIND=8), DIMENSION(k_dim) :: x
    LOGICAL, DIMENSION(mesh%np) :: if_udotn_zero

    !===mass
    CALL st_csr(mesh%jj, mass%ia, mass%ja)
    ALLOCATE(mass%aa(SIZE(mass%ja)))
    mass%aa = 0.d0
    CALL qs_00_M (mesh, 1.d0, mass%ia, mass%ja, mass%aa)

    !===diag
    ALLOCATE(diag(mesh%np))
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          IF (i==mass%ja(p)) THEN
             diag(i) = p
             EXIT
          END IF
       END DO
    END DO

    !===fix_roundoff
    ALLOCATE(fix_roundoff(mesh%np))

    !===lumped
    ALLOCATE(lumped(mesh%np))
    DO i = 1, mesh%np
       lumped(i) = SUM(mass%aa(mass%ia(i):mass%ia(i+1)-1))
    END DO

    !===Lump mass matrix at boundary
    !DO ni = 1, SIZE(js_D)
    !   i = js_D(ni)
    !   mass%aa(mass%ia(i):mass%ia(i+1)-1)=0.d0
    !   mass%aa(diag(i)) = lumped(i)
    !END DO

    !===Mass - lumped
    CALL st_csr(mesh%jj, mc_minus_ml%ia, mc_minus_ml%ja)
    ALLOCATE(mc_minus_ml%aa(SIZE(mc_minus_ml%ja)))
    mc_minus_ml%aa = mass%aa
    mc_minus_ml%aa(diag) = mc_minus_ml%aa(diag) - lumped

    !===Mesh size (for fake SGN)
    ALLOCATE(mesh_size(mesh%np))
    IF (k_dim==1) THEN
       mesh_size = lumped
    ELSE
       mesh_size = SQRT(lumped)
    END IF

    !===Pre mass
    CALL st_csr(mesh%jj, pre_mass%ia, pre_mass%ja)
    ALLOCATE(pre_mass%aa(SIZE(pre_mass%ja)))
    DO i = 1, mesh%np
       pre_mass%aa(pre_mass%ia(i):pre_mass%ia(i+1)-1) = mass%aa(pre_mass%ia(i):pre_mass%ia(i+1)-1)/lumped(i)
    END DO

    !===dij
    CALL st_csr(mesh%jj, dij%ia, dij%ja)
    ALLOCATE(dij%aa(SIZE(dij%ja)))
    dij%aa = 0.d0

    !===betaij
    CALL st_csr(mesh%jj, betaij%ia, betaij%ja)
    ALLOCATE(betaij%aa(SIZE(betaij%ja)))
    IF (k_dim==2) THEN
       CALL compute_betaij
    ELSE
       betaij%aa=1.d0
    END IF

    !===muij
    CALL st_csr(mesh%jj, muij%ia, muij%ja)
    ALLOCATE(muij%aa(SIZE(muij%ja)))
    muij%aa = 0.d0

    !===dijL
    CALL st_csr(mesh%jj, dijL%ia, dijL%ja)
    ALLOCATE(dijL%aa(SIZE(dijL%ja)))
    dijL%aa = 0.d0

    !===muijL
    CALL st_csr(mesh%jj, muijL%ia, muijL%ja)
    ALLOCATE(muijL%aa(SIZE(muijL%ja)))
    muijL%aa = 0.d0

    !===fctmat
    DO d = 1, inputs%syst_size
       CALL st_csr(mesh%jj, fctmat(d)%ia, fctmat(d)%ja)
       ALLOCATE(fctmat(d)%aa(SIZE(fctmat(d)%ja)))
       fctmat(d)%aa = 0.d0
    END DO

    !===lij
    CALL st_csr(mesh%jj, lij%ia, lij%ja)
    ALLOCATE(lij%aa(SIZE(lij%ja)))
    lij%aa = 0.d0

    !===cij = \int_K \GRAD(\phi_j) \phi_i \dif x
    ! if_udotn_zero = .FALSE.
    ! if_udotn_zero(udotn_js_D) = .TRUE.   !ERIC
    ! if_udotn_zero(h_js_D) = .FALSE.
    DO d = 1, k_dim
       CALL st_csr(mesh%jj, cij(d)%ia, cij(d)%ja)
       ALLOCATE(cij(d)%aa(SIZE(cij(d)%ja)))
       cij(d)%aa = 0.d0
    END DO
    IF (inputs%type_test==9) THEN
       DO m = 1, mesh%me
          DO ni = 1, mesh%gauss%n_w
             i = mesh%jj(ni, m)
             DO nj = 1, mesh%gauss%n_w
                j = mesh%jj(nj, m)
                DO d = 1, k_dim
                   !TEST TEST TEST MALPASSET
                   x(d) = -SUM(mesh%gauss%dw(d,ni,:,m) * mesh%gauss%ww(nj,:)*mesh%gauss%rj(:,m))
                   !TEST TEST TEST MALPASSE
                END DO
                DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
                   IF (cij(1)%ja(p) == j) THEN
                      DO d = 1, k_dim
                         cij(d)%aa(p) = cij(d)%aa(p) + x(d)
                      END DO
                      EXIT
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ELSE
       DO m = 1, mesh%me
       DO ni = 1, mesh%gauss%n_w
          i = mesh%jj(ni, m)
          IF (if_udotn_zero(i)) THEN
             DO nj = 1, mesh%gauss%n_w
                j = mesh%jj(nj, m)
                DO d = 1, k_dim
                   x(d) = -SUM(mesh%gauss%dw(d,ni,:,m) * mesh%gauss%ww(nj,:)*mesh%gauss%rj(:,m))
                END DO
                DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
                   IF (cij(1)%ja(p) == j) THEN
                      DO d = 1, k_dim
                         cij(d)%aa(p) = cij(d)%aa(p) + x(d)
                      END DO
                      EXIT
                   ENDIF
                ENDDO
             ENDDO
          ELSE
             DO nj = 1, mesh%gauss%n_w
                j = mesh%jj(nj, m)
                DO d = 1, k_dim
                   x(d) = SUM(mesh%gauss%dw(d,nj,:,m) * mesh%gauss%ww(ni,:)*mesh%gauss%rj(:,m))
                END DO
                DO p = cij(1)%ia(i), cij(1)%ia(i+1) - 1
                   IF (cij(1)%ja(p) == j) THEN
                      DO d = 1, k_dim
                         cij(d)%aa(p) = cij(d)%aa(p) + x(d)
                      END DO
                      EXIT
                   ENDIF
                ENDDO
             ENDDO
          END IF
       ENDDO
    ENDDO
    END IF

    !===entropy viscosity matrix
    CALL st_csr(mesh%jj, resij%ia, resij%ja)
    ALLOCATE(resij%aa(SIZE(resij%ja)))
    resij%aa = 0.d0

    !===Maximum number of points in stencil
    max_nb_pt_stencil = 0
    DO i = 1, mesh%np
       max_nb_pt_stencil = MAX(max_nb_pt_stencil, mass%ia(i+1)-mass%ia(i))
    END DO


  END SUBROUTINE construct_matrices

  SUBROUTINE euler(un,unext)
    USE mesh_handling
    USE boundary_conditions
    USE fct
    !USE pardiso_solve
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(OUT) :: unext
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: rk, ulow, du
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: source
    REAL(KIND=8), DIMENSION(mesh%np)                   :: ff, deltah
    REAL(KIND=8), DIMENSION(mesh%np)                   :: rho_e_max, hmin, hmax, psi_small
    LOGICAL,      DIMENSION(mesh%np)                   :: Hdry
    REAL(KIND=8), DIMENSION(k_dim)                     :: nij, ur, ul
    REAL(KIND=8), DIMENSION(max_nb_pt_stencil)         :: ubarij, ubarji
    REAL(KIND=8), DIMENSION(k_dim+1,max_nb_pt_stencil) :: ubar
    REAL(KIND=8), DIMENSION(k_dim+1)                   :: xx
    REAL(KIND=8) :: Hstarij, Hstarji, rklocij, rklocji, ratij, ratji, &
         lambda, gammai, gammaj, sigmaij, test, oneoverh
    INTEGER :: p, p_start, p_end, i, j, k, d, ij, comp
    LOGICAL, SAVE :: once=.TRUE.
    IF (once) THEN
       isolve=-1
       isolve_m=-1
       once=.FALSE.
    END IF

    un_over_h = compute_divide_by_h(un)
    velocity = un_over_h(2:2+k_dim-1,:)

    !===Galerkin
    IF (inputs%viscosity_type=='galerkin') THEN
       dij%aa = 0.d0
       muij%aa = 0.d0
       CALL smb_2(un,rk)
       CALL divide_by_lumped(rk)
       IF (inputs%if_lumped) THEN
          unext = un+inputs%dt*rk
          !unext(1,:)  = max(un(1,:) + inputs%dt*rk(1,:),0.d0)
          !unext(2:,:) = un(2:,:)+inputs%dt*rk(2:,:)
       ELSE
          DO k = 1, inputs%syst_size
             !CALL solve_pardiso(pre_mass%aa,pre_mass%ia,pre_mass%ja,rk(k,:),ff,isolve,2)
             isolve=ABS(isolve)
             unext(k,:) = un(k,:)+inputs%dt*ff
             !IF (k==1) THEN
             !   unext(k,:)  = max(un(k,:) + inputs%dt*ff,0.d0)
             !ELSE
             !   unext(k,:) = un(k,:)+inputs%dt*ff
             !END IF
          END DO
       END IF
       RETURN
    END IF

    !===Compute first-order viscosity
    CALL compute_dij(un,.FALSE.)
    CALL compute_muij(un)
    dij%aa = MAX(dij%aa,muij%aa)
    IF (inputs%viscosity_type=='fct') THEN
       dijL%aa = dij%aa
       muijL%aa = muij%aa
    END IF

    !===Alpha viscosity
    IF (inputs%if_alpha_limit) THEN
       CALL alpha_limit(un(1,:))
    END IF
    !IF (inputs%viscosity_type=='fct') THEN
    !   dijL%aa = dij%aa
    !   muijL%aa = muij%aa
    !END IF

    !===Compute right-hand side
    IF (inputs%viscous_type=='type1') THEN
       CALL smb_1(un,rk)
    ELSE IF (inputs%viscous_type=='type2') THEN
       IF (if_roundoff_fix) THEN
          CALL smb_2_roundoff(un,rk)
       ELSE
          CALL smb_2(un,rk)
       END IF
    ELSE
       WRITE(*,*) ' BUG in Euler, viscous_type'
       STOP
    END IF
    CALL divide_by_lumped(rk)

    !===Compute First-Order solution
    IF (inputs%viscous_type=='type1')  THEN
       unext = un+inputs%dt*rk
    ELSE IF (inputs%viscous_type=='type2')  THEN
       IF (if_roundoff_fix) THEN
          unext(1,:)  = un(1,:)*(1+inputs%dt*fix_roundoff/lumped) + inputs%dt*rk(1,:)
          unext(2:,:) = un(2:,:)+inputs%dt*rk(2:,:)
       ELSE
          unext = un+inputs%dt*rk
       END IF
    ELSE
       WRITE(*,*) ' BUG in euler, wrong inputs%viscous_type'
       STOP
    END IF
    IF (inputs%viscosity_type=='viscous') THEN
       CALL check_Hmin(unext)
       RETURN
    END IF
    !write(*,*) ' ratios', maxval(limit_h/lumped),maxval(lumped/limit_h)
    !stop
    !===We assume below that we use either 'entropy_visc' or 'fct'
    IF (inputs%viscosity_type=='fct') THEN
       CALL friction(un,source)
       CALL maxmin(un(1,:),mass,hmax,hmin)
       DO i = 1, mesh%np
          deltah(i) = 0.d0
          DO p = mass%ia(i), mass%ia(i+1) - 1
             j = mass%ja(p)
             IF (i==j) CYCLE
             deltah(i) = deltah(i) + un(1,j) - un(1,i)
          END DO
          deltah(i)= ABS(deltah(i))
       END DO
       deltah = hmax-hmin
       hmin = max(hmin,0.d0)
       rho_e_max = 0.0d0
       !CALL estimate_rho_e_max(un,rho_e_max)
       DO i = 1, mesh%np
          p_start = mass%ia(i)
          p_end   = mass%ia(i+1)-1
          ubarij = 0.d0
          ubarji = 0.d0
          ij = 0
          xx = 0.d0
          test =un(1,i)
          DO p = p_start, p_end
             ij = ij + 1
             j = mass%ja(p)
             IF(i==j) THEN
                ubarij(ij)=un(1,i)
                ubarji(ij)=un(1,i)
                ubar(:,ij) =  un(1:k_dim+1,i)
             ELSE
                Hstarji = MAX(0.d0,un(1,j)+bath(j)-MAX(bath(i),bath(j)))
                Hstarij = MAX(0.d0,un(1,i)+bath(i)-MAX(bath(i),bath(j)))
                ratji = Hstarji/max(un(1,j),inputs%htiny)
                ratij = Hstarij/max(un(1,i),inputs%htiny)
                gammai = muijL%aa(p) + (dijL%aa(p)-muijL%aa(p))*Hstarij/max(un(1,i),inputs%htiny)
                gammaj = muijL%aa(p) + (dijL%aa(p)-muijL%aa(p))*Hstarji/max(un(1,j),inputs%htiny)
                DO comp = 1, k_dim+1
                   rklocij = 0.d0
                   DO k = 1, k_dim
                      rklocij  = rklocij - cij(k)%aa(p)*(velocity(k,j)*un(comp,j)-velocity(k,i)*un(comp,i))
                   END DO
                   IF (comp==1) ubarij(ij) = (rklocij + gammaj*un(1,j)+gammai*un(1,i))/max(2.d0*gammai,1.d-14)
                   IF (comp>1) THEN
                      rklocij = rklocij - cij(comp-1)%aa(p)*inputs%gravity*(un(1,j)**2-un(1,i)**2)/2
                      source(comp,i) =  source(comp,i) &
                           + inputs%gravity*cij(comp-1)%aa(p)*(-un(1,i)*(bath(j)-bath(i))+(un(1,j)-un(1,i))**2/2)
                   END IF
                   ubar(comp,ij) = rklocij/max(2.d0*dijL%aa(p),1.d-14) + (un(comp,j)+un(comp,i))/2 &
                        +((ratji-1.d0)*un(comp,j)-(ratij-1.d0)*un(comp,i)) &
                        *(dijL%aa(p)-muijL%aa(p))/max(2.d0*dijL%aa(p),1.d-14)
                END DO
                test = max(test,ubarij(ij))
                ubarji(ij) = Hstarji
                ubarij(ij) = Hstarij
                !xx = xx + inputs%dt*(2.d0*gammai*(ubar(:,ij)-un(:,i)))/lumped(i)
                xx = xx + inputs%dt*(2.d0*dijL%aa(p)*(ubar(:,ij)-un(1:k_dim+1,i)))/lumped(i)
             END IF
          END DO
          source(:,i) = source(:,i)*inputs%dt/lumped(i)
          !xx = xx + un(:,i) + source(:,i) !==xx is the low-order solution
          !WRITE(*,*) 'error xx',MAXVAL(ABS(xx-unext(:,i)))
          xx = xx + un(1:k_dim+1,i) !===low-order - source
          !write(*,*) MAXVAL(ubar(1,1:ij))-unext(1,i)
          hmax(i) = MAX(hmax(i),MAXVAL(ubar(1,1:ij))) !===Essential to get 2nd order, 1D Mannings
          hmin(i) = max(MIN(hmin(i),MINVAL(ubar(1,1:ij))),0.d0) !===Essential to get 2nd order, 1D Mannings
          !hmax(i) = MAXVAL(ubar(1,1:ij)) !===Not good
          !hmin(i) = MINVAL(ubar(1,1:ij)) !===Not good
          !write(*,*)  deltah(i), limit_h(i)
          !IF (unext(1,i) .LE. deltah(i) .OR. minval(ubarij(1:ij)).le.0) THEN
          !IF (unext(1,i) .LE. deltah(i)) THEN
          !IF (minval(ubarij(1:ij)).le.0) THEN
          !IF (unext(1,i) .LE. deltah(i)/2) THEN
          !IF (unext(1,i) .LE. limit_h(i)) THEN !===Dimensionionally correct
          IF (unext(1,i) .LE. 0.2*lumped(i)**(3-k_dim)/inputs%max_water_h) THEN !===BEST, but not dimensionally correct
             hdry(i) = .true.
          ELSE
             hdry(i) = .false.
          END IF
          !write(*,*) 'max ', unext(1,i), hmin(i)
!!$          IF (unext(1,i)> hmax(i)) THEN
!!$             write(*,*) 'max violation', xx, unext(1,i), hmax(i), test
!!$             stop
!!$          END IF
!!$
!!$          IF (unext(1,i)< hmin(i)) THEN
!!$             write(*,*) 'min violation', unext(1,i), hmin(i)
!!$             stop
!!$          END IF
          DO k = 1, ij
             !=== Kinetic energy: 0.5*(||q||^2/h)
             oneoverh = 2*ubar(1,k)/(ubar(1,k)**2+max(ubar(1,k),inputs%htiny)**2)
             rho_e_max(i) = MAX(rho_e_max(i),oneoverh*SUM(ubar(2:k_dim+1,k)**2)/2)
          END DO
          !oneoverh = 2*xx(1)/(xx(1)**2+max(xx(1),inputs%htiny)**2)
          !rho_e_max(i) = MAX(rho_e_max(i),oneoverh*SUM((xx(2:k_dim+1))**2)/2)
          !WRITE(*,*) ' rho_e_max - rhoe(unext)', rho_e_max(i)-oneoverh*SUM((xx(2:k_dim+1))**2)/2, xx(1)
       END DO
       dijL%aa=dij%aa
       muijL%aa=muij%aa
    END IF

    !===Compute entropy viscosity
    !CALL entropy_residual(un)
    CALL entropy_commutator(un)
    DO i = 1, mesh%np
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          dij%aa(p) = MIN(dij%aa(p),  1.d0*resij%aa(p))
          muij%aa(p)= MIN(muij%aa(p), 1.d0*resij%aa(p))
       END DO
    END DO
    !===If entropy viscosity only; no FCT
    IF (inputs%viscosity_type=='entropy_visc') THEN
       IF (inputs%viscous_type=='type1') THEN
          WRITE(*,*) ' Bug: entropy viscosity programmed only with type 2'
          STOP
       ELSE IF (inputs%viscous_type=='type2') THEN
          IF (if_roundoff_fix) THEN
             CALL smb_2_roundoff(un,rk)
          ELSE
             CALL smb_2(un,rk)
          END IF
       END IF
       CALL divide_by_lumped(rk)
       !===Solve and update
       IF (inputs%if_lumped) THEN
          !===Compute entropy viscosity solution
          IF (if_roundoff_fix) THEN
             unext(1,:)  = un(1,:)*(1+inputs%dt*fix_roundoff/lumped) + inputs%dt*rk(1,:)
             unext(2:,:) = un(2:,:)+inputs%dt*rk(2:,:)
          ELSE
             unext = un+inputs%dt*rk
          END IF
       ELSE
          DO k = 1, inputs%syst_size
             !CALL solve_pardiso(pre_mass%aa,pre_mass%ia,pre_mass%ja,rk(k,:),ff,isolve,2)
             isolve=ABS(isolve)
             unext(k,:) = un(k,:)+inputs%dt*ff
          END DO
       END IF
       !===If entropy viscosity + FCT
    ELSE IF (inputs%viscosity_type=='fct') THEN   !===Fct limitation using smb_2
       ulow = unext
       !WRITE(*,*) ' CHECK ulow before'
       CALL check_Hmin(ulow) !===low-order solution
       !WRITE(*,*) ' CHECK ulow AFTER'
       IF (inputs%viscous_type/='type2') THEN
          WRITE(*,*) ' BUG: FCT programmed with type2 only'
          STOP
       END IF
       IF (inputs%if_lumped) THEN
          du = 0.d0
       ELSE
          !===Compute right-hand side
          IF (inputs%viscous_type=='type1') THEN
             CALL smb_1(un,rk)
          ELSE IF (inputs%viscous_type=='type2') THEN
             CALL smb_2(un,rk)
          END IF
          CALL divide_by_lumped(rk)
          !===Compute entropy viscosity solution
          DO k = 1, inputs%syst_size
             !CALL solve_pardiso(pre_mass%aa,pre_mass%ia,pre_mass%ja,rk(k,:),ff,isolve,2)
             isolve=ABS(isolve)
             unext(k,:) = un(k,:)+inputs%dt*ff
          END DO
          du = unext-un
       END IF

       CALL compute_fctk_matrix(du,un,fctmat)
       inputs%limiter_type='avg'
       !CALL relax(un(1,:),hmin,hmax)
       !ERIC
       !CALL FCT_generic(ulow(1,:),hdry,hmax,hmin,fctmat(1),lumped,lij)
       !CALL FCT_positivity(ulow(1,:),hdry,hmax,hmin,lumped,fctmat(1),lij)
       !IF (inputs%type_test==9) THEN !===Malpasset
       !   psi_small = rho_e_max*inputs%htiny
       !   source(1,:)=0.d0
       !   CALL quadratic_limiting(ulow-source,rho_e_max,psi_small)
       !   CALL transpose_op(lij,'min')
       !END IF
       CALL apply_fctk_matrix(rk,fctmat,lij)
       !write(*,*) ' CHECK mass unext before', sum(lumped*(ulow(1,:)))
       CALL divide_by_lumped(rk)
       unext = ulow+rk
       !unext(1,:) = MAX(unext(1,:),0.d0)
       !write(*,*) ' CHECK mass after', sum(lumped*(unext(1,:)))
    END IF
    !WRITE(*,*) ' CHECK unext before'
    CALL check_Hmin(unext)
    !WRITE(*,*) ' CHECK ulow before'
    RETURN

  END SUBROUTINE euler

  SUBROUTINE compute_dij(un,if_dt)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN)  :: un
    LOGICAL :: if_dt
    INTEGER                                       :: i, p, j, d
    REAL(KIND=8)                                  :: norm_cij, lambda
    REAL(KIND=8), DIMENSION(k_dim)                :: nij
    REAL(KIND=8), DIMENSION(inputs%syst_size)     :: ur, ul

    !===Viscosity using compute_lambda
    DO i = 1, mesh%np
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          IF (i.NE.j) THEN
             DO d = 1, k_dim
                nij(d) = cij(d)%aa(p) !=== definition of cij is same as in the paper
             END DO
             norm_cij = SQRT(SUM(nij**2))
             nij=nij/norm_cij
             CALL compute_lambda_vacc(un(:,i),un(:,j),velocity(:,i),velocity(:,j),&
                  mesh_size(i),mesh_size(j),nij,lambda,if_dt)
             dij%aa(p) = norm_cij*lambda
          ELSE
             dij%aa(p) = 0.d0
          END IF
       END DO
    END DO
    CALL transpose_op(dij,'max')

    DO i = 1, mesh%np
       dij%aa(diag(i)) = -SUM(dij%aa(dij%ia(i):dij%ia(i+1)-1))
    END DO

    RETURN

  END SUBROUTINE compute_dij

  SUBROUTINE compute_muij(un)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN)  :: un
    INTEGER                                       :: i, p, j, d
    REAL(KIND=8)                                  :: norm_cij, lambda
    REAL(KIND=8), DIMENSION(k_dim)                :: nij
    REAL(KIND=8), DIMENSION(k_dim)                :: ur, ul

    !===Viscosity using speed only
    DO i = 1, mesh%np
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          IF (i.NE.j) THEN
             DO d = 1, k_dim
                nij(d) = cij(d)%aa(p) !=== definition of cij is the same as in the paper
             END DO
             ul = velocity(:,i)
             ur = velocity(:,j)
             !TESTTTTT
             lambda=MAX(MAX(-SUM(nij*ul),0.d0),MAX(SUM(nij*ur),0.d0))
             !lambda=MAX(ABS(SUM(nij*ul)),ABS(SUM(nij*ur)))
             !TESTTTTT
             muij%aa(p) = lambda
          ELSE
             muij%aa(p) = 0.d0
          END IF
       END DO
    END DO
    CALL transpose_op(muij,'max')
    DO i = 1, mesh%np
       muij%aa(diag(i)) = -SUM(muij%aa(dij%ia(i):muij%ia(i+1)-1))
    END DO

    RETURN

  END SUBROUTINE compute_muij

  SUBROUTINE transpose_op(mat,type)
    IMPLICIT NONE
    TYPE(matrice_bloc), INTENT(INOUT):: mat
    CHARACTER(LEN=3),  INTENT(IN)   :: type
    INTEGER, DIMENSION(SIZE(mat%ia)) :: iao
    INTEGER:: i, j, p, next
    IF  (type/='min' .AND. type/='max') THEN
       WRITE(*,*) ' BUG in tanspose_op'
       STOP
    END IF
    iao = mat%ia
    DO i = 1, SIZE(mat%ia)-1
       DO p = mat%ia(i), mat%ia(i+1)-1
          j = mat%ja(p)
          next = iao(j)
          iao(j) = next+1
          IF (j.LE.i) CYCLE
          IF (type=='min') THEN
             mat%aa(next) = min(mat%aa(p),mat%aa(next))
             mat%aa(p) = mat%aa(next)
          ELSE IF (type=='max') THEN
             mat%aa(next) = max(mat%aa(p),mat%aa(next))
             mat%aa(p) = mat%aa(next)
          END IF
       END DO
    END DO
  END SUBROUTINE transpose_op

  SUBROUTINE smb_1(un,rk)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(OUT) :: rk
    REAL(KIND=8), DIMENSION(inputs%syst_size,k_dim,mesh%np)        :: vv
    INTEGER :: d, i, j, k, p
    REAL(KIND=8) :: xx, Hstarij, Hstarji, ratij, ratji

    vv=flux(un)

    rk=0.d0
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          Hstarij = MAX(0.d0,un(1,i)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(1,j)+bath(j)-MAX(bath(i),bath(j)))
          ratij = Hstarij/un(1,i)
          ratji = Hstarji/un(1,j)
          IF (un(1,i).LE.inputs%htiny) THEN !Important for long-time WB
             ratij=0.d0
             ratji=0.d0
          ELSE
             ratij = Hstarij/un(1,i)
          END IF
          IF (un(1,j).LE.inputs%htiny) THEN !Important for long-time WB
             ratji=0.d0
             ratij=0.d0
          ELSE
             ratji = Hstarji/un(1,j)
          END IF

          DO k = 1, inputs%syst_size
             xx = 0.d0
             DO d = 1, k_dim
                xx = xx - cij(d)%aa(p)*(vv(k,d,j)*ratji + vv(k,d,i)*ratij)
             END DO
             rk(k,i) = rk(k,i) + xx + dij%aa(p)*(un(k,j)*ratji-un(k,i)*ratij)
          END DO
          DO k = 1, k_dim
             rk(k+1,i) = rk(k+1,i) &
                  - 0.5d0*inputs%gravity*(Hstarji**2 - Hstarij**2)*cij(k)%aa(p)
          END DO
       END DO
    END DO

    SELECT CASE(inputs%type_test)
    CASE(9,10,11,12,13,15)
       CALL friction(un,rk)
    END SELECT

    IF (inputs%if_FGN) THEN
       CALL FGN_rhs(un,rk)
    END IF

  END SUBROUTINE smb_1

  SUBROUTINE friction(un,rk)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(OUT) :: rk
    REAL(KIND=8), DIMENSION(inputs%syst_size,k_dim,mesh%np)        :: vv
    REAL(KIND=8), DIMENSION(mesh%np)  :: hloc_star, hstar, vel
    INTEGER :: d, i, j, k, p
    REAL(KIND=8) :: fric
    REAL(KIND=8), PARAMETER :: chi=1.d0

    IF (inputs%type_test==12 .OR. inputs%type_test==13) THEN
       DO k = 1, k_dim
          rk(k+1,:) = rk(k+1,:) - inputs%mannings*lumped*un(k+1,:)
       END DO
       RETURN
    END IF

    fric = inputs%gravity*inputs%mannings**2
    hstar = MAX(un(1,:),inputs%htiny)**(1.d0+inputs%eta)
    IF (k_dim==2) THEN
       vel = fric*SQRT(velocity(1,:)**2+velocity(2,:)**2)
    ELSE
       vel = fric*ABS(velocity(1,:))
    END IF
    hloc_star = chi*vel*inputs%dt

    DO i = 1, mesh%np
       DO k = 2, inputs%syst_size
          rk(k,i) = rk(k,i)  - lumped(i)*un(k,i)*vel(i)*2/(hstar(i) + MAX(hstar(i),hloc_star(i)))
       END DO
    END DO

  END SUBROUTINE friction

  SUBROUTINE smb_2_roundoff(un,rk)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(OUT) :: rk
    REAL(KIND=8), DIMENSION(inputs%syst_size,k_dim,mesh%np)        :: vv
    INTEGER :: d, i, j, k, p
    REAL(KIND=8) :: xx, hmean, Hstarij, Hstarji, ratij, ratji
    vv=flux(un)
    rk=0.d0
    fix_roundoff = 0.d0
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          Hstarij = MAX(0.d0,un(1,i)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(1,j)+bath(j)-MAX(bath(i),bath(j)))

          IF (un(1,i).LE.inputs%htiny) THEN !Important for long-time WB
             ratij=1.d0
             ratji=1.d0
          ELSE
             ratij = Hstarij/un(1,i)
          END IF

          IF (un(1,j).LE.inputs%htiny) THEN !Important for long-time WB
             ratji=1.d0
             ratji=1.d0
          ELSE
             ratji = Hstarji/un(1,j)
          END IF

          !===mass
          k=1
          xx = 0.d0
          DO d = 1, k_dim
             xx = xx - cij(d)%aa(p)*(velocity(d,j))
          END DO
          IF (i==j) THEN
             fix_roundoff(i) = fix_roundoff(i) + xx
          ELSE
             !===Fix roundoff error
             fix_roundoff(i) = fix_roundoff(i) - muij%aa(p) + (dij%aa(p)-muij%aa(p))*(-ratij)
             rk(k,i) = rk(k,i) + (xx + muij%aa(p))*un(k,j) + (dij%aa(p)-muij%aa(p))*Hstarji
          END IF
          !===rest
          DO k = 2, inputs%syst_size
             xx = 0.d0
             DO d = 1, k_dim
                xx = xx - cij(d)%aa(p)*(vv(k,d,j))
             END DO
             rk(k,i) = rk(k,i) + xx + (dij%aa(p)-muij%aa(p))*(un(k,j)*ratji-un(k,i)*ratij) &
                  + muij%aa(p)*(un(k,j)-un(k,i))
          END DO
          !===Topography contribution
          DO k = 1, k_dim
             rk(k+1,i) = rk(k+1,i) &
                  - inputs%gravity*un(1,i)*(un(1,j)+bath(j))*cij(k)%aa(p)
          END DO
       END DO
    END DO

    SELECT CASE(inputs%type_test)
    CASE(9,10,11,12,13,15)
       CALL friction(un,rk)
    END SELECT

    IF (inputs%if_FGN) THEN
       CALL FGN_rhs(un,rk)
    END IF

  END SUBROUTINE smb_2_roundoff


  SUBROUTINE apply_viscosity(un,rk)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(OUT) :: rk
    REAL(KIND=8), DIMENSION(inputs%syst_size) :: xx
    INTEGER :: i, j, k, p
    REAL(KIND=8) :: Hstarij, Hstarji, ratij, ratji

    DO i = 1, mesh%np
       xx = 0.d0
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          Hstarij = MAX(0.d0,un(1,i)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(1,j)+bath(j)-MAX(bath(i),bath(j)))
          !IF (un(1,i).LE.0.d0) THEN
          IF (un(1,i).LE.inputs%htiny) THEN
             ratij=0.d0
             ratji=0.d0
          ELSE
             ratij = Hstarij/un(1,i)
          END IF
          !IF (un(1,j).LE.0.d0) THEN
          IF (un(1,j).LE.inputs%htiny) THEN
             ratji=0.d0
             ratji=0.d0
          ELSE
             ratji = Hstarji/un(1,j)
          END IF
          DO k = 1, inputs%syst_size
             xx(k) = xx(k) + (dij%aa(p)-muij%aa(p))*(un(k,j)*ratji-un(k,i)*ratij) &
                  + muij%aa(p)*(un(k,j)-un(k,i))
          END DO
       END DO
       rk(:,i) = xx
    END DO
  END SUBROUTINE apply_viscosity

  SUBROUTINE apply_fctk_matrix(rk,fctmat,lij)
    USE mesh_handling
    IMPLICIT NONE
    TYPE(matrice_bloc), DIMENSION(inputs%syst_size),   INTENT(IN) :: fctmat
    TYPE(matrice_bloc),                       INTENT(IN) :: lij
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)             :: rk
    INTEGER :: i, k, p0, p1
    DO i = 1, mesh%np
       p0=lij%ia(i)
       p1=lij%ia(i+1)-1
       DO k = 1, inputs%syst_size
          rk(k,i) = SUM(fctmat(k)%aa(p0:p1)*lij%aa(p0:p1))
       END DO
    END DO
  END SUBROUTINE apply_fctk_matrix

  SUBROUTINE compute_fctk_matrix(du,un,fctmat)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN) :: du, un
    TYPE(matrice_bloc), DIMENSION(inputs%syst_size),   INTENT(IN) :: fctmat
    INTEGER :: i, j, k, p
    REAL(KIND=8) :: Hstarij, Hstarji, ratij, ratji

    DO i = 1, mesh%np
       DO p = fctmat(1)%ia(i), fctmat(1)%ia(i+1) - 1
          j = fctmat(1)%ja(p)
          IF (i==j) THEN
             DO k = 1, inputs%syst_size
                fctmat(k)%aa(p) = 0.d0
             END DO
             CYCLE
          END IF
          Hstarij = MAX(0.d0,un(1,i)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(1,j)+bath(j)-MAX(bath(i),bath(j)))
          !fctmat(1)%aa(p) =  -mc_minus_ml%aa(p)*(du(1,j)-du(1,i)) &
          !     + inputs%dt*(dij%aa(p)-dijL%aa(p)-(muij%aa(p)-muijL%aa(p)))*(Hstarji-Hstarij) &
          !     + inputs%dt*(muij%aa(p)-muijL%aa(p))*(un(1,j)-un(1,i))

          DO k = 1, inputs%syst_size
             fctmat(k)%aa(p) = -mc_minus_ml%aa(p)*(du(k,j)-du(k,i)) &
                  + inputs%dt*(dij%aa(p)-dijL%aa(p)-(muij%aa(p)-muijL%aa(p)))*(un_over_h(k,j)*Hstarji-un_over_h(k,i)*Hstarij) &
                  + inputs%dt*(muij%aa(p)-muijL%aa(p))*(un(k,j)-un(k,i))
          END DO

       END DO
    END DO
    !write(*,*) ' test fctmat', sum(fctmat(1)%aa)
  END SUBROUTINE compute_fctk_matrix

  SUBROUTINE alpha_limit(un)
    USE mesh_handling
    USE boundary_conditions
    USE sub_plot
    IMPLICIT NONE
    INTEGER :: i, j, p
    REAL(KIND=8), DIMENSION(mesh%np), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(mesh%np) :: alpha, alphastar
    REAL(KIND=8) :: denom, num, h_threshold
    REAL(KIND=8) :: denomstar, numstar, Hstarij, Hstarji
    REAL(KIND=8) :: amax, amaxstar
    !REAL(KIND=8), PARAMETER :: alphath=0.5d0, omalphath=1.d0/(1.d0-alphath)
    !INTEGER,      PARAMETER :: exponent=2 !Test in paper done with exponent=2.
    REAL(KIND=8), PARAMETER :: alphath=0.0d0, omalphath=1.d0/(1.d0-alphath)
    INTEGER,      PARAMETER :: exponent=2 !Test in paper done with exponent=2.

    DO i = 1, mesh%np
       !===Keep full viscosity in dry regions
       IF (un(i).LE. limit_h(i)) THEN
          !IF (un(i).LE. inputs%htiny) THEN
          alpha(i) = 1.d0
          CYCLE
       END IF
       denom =0.d0
       num  = 0.d0
       numstar = 0.d0
       denomstar = 0.d0
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          num = num + un(j) - un(i)
          denom = denom + ABS(un(j) - un(i))
          !num = num + betaij%aa(p)*(un(j) - un(i))
          !denom = denom + betaij%aa(p)*ABS(un(j) - un(i))
       END DO

       IF (ABS(num)<inputs%htiny) THEN
          alpha(i)=0.d0
       ELSE
          alpha(i) = ((ABS(num) - inputs%htiny )/ABS(denom-inputs%htiny))
       END IF
       alpha(i) = MAX(alpha(i)-alphath,0.d0)*omalphath
       alpha(i) = alpha(i)**exponent

    END DO
    IF (inputs%time + inputs%dt > inputs%Tfinal) THEN
       IF (k_dim==2)   CALL plot_scalar_field(mesh%jj, mesh%rr, alpha, 'alpha.plt')
    END IF
    !===Limit viscosity
    DO i = 1, mesh%np
       DO p = dij%ia(i), dij%ia(i+1) - 1
          j = dij%ja(p)
          IF (i.NE.j) THEN
             amax = MAX(alpha(i),alpha(j))
             dij%aa(p)  = dij%aa(p)*amax
             muij%aa(p) = muij%aa(p)*amax
          ELSE
             dij%aa(p)  = 0.d0
             muij%aa(p) = 0.d0
          END IF
       END DO
    END DO

    DO i = 1, mesh%np
       dij%aa(diag(i))  = -SUM(dij%aa(dij%ia(i):dij%ia(i+1)-1))
       muij%aa(diag(i)) = -SUM(muij%aa(dij%ia(i):dij%ia(i+1)-1))
    END DO

  END SUBROUTINE alpha_limit

  SUBROUTINE divide_by_lumped(rk)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: rk
    INTEGER :: k
    DO k = 1, inputs%syst_size
       rk(k,:) = rk(k,:)/lumped
    END DO
  END SUBROUTINE divide_by_lumped

  SUBROUTINE check_Hmin(h)
    USE boundary_conditions
    USE mesh_handling
    USE sub_plot
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: h
    INTEGER :: i
    LOGICAL, DIMENSION(mesh%np) :: check
    !check=.TRUE.   ERIC
    !check(js_D)=.FALSE.
    DO i = 1, mesh%np
       !IF (.NOT.check(i)) CYCLE
       IF (h(1,i)<0.d0) THEN
          WRITE(*,*) 'Min h<0, STOP', h(1,i)
          velocity(1,ux_js_D) = 0.d0
          IF (k_dim==2) velocity(2,uy_js_D) = 0.d0
          WRITE(*,*) 'MAXVAL(vel)', MAXVAL(ABS(velocity(1,:))), MAXVAL(ABS(velocity(k_dim,:)))
          CALL plot_scalar_field(mesh%jj, mesh%rr, h(1,:), 'h.plt')
          CALL plot_scalar_field(mesh%jj, mesh%rr, velocity(1,:), 'vx.plt')
          if (k_dim==2) CALL plot_scalar_field(mesh%jj, mesh%rr, velocity(2,:), 'vy.plt')
          STOP
       END IF
    END DO
  END SUBROUTINE check_Hmin

!!$  SUBROUTINE entropy_residual(un)
!!$    USE mesh_handling
!!$    USE pardiso_solve
!!$    USE boundary_conditions
!!$    IMPLICIT NONE
!!$    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
!!$    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: rk, unext, Entprime
!!$    REAL(KIND=8), DIMENSION(k_dim,mesh%np)  :: velnext
!!$    REAL(KIND=8), DIMENSION(mesh%np)        :: res, Ent, Entnext, maxn, minn
!!$    REAL(KIND=8), DIMENSION(mesh%np)  :: ff, rescale
!!$
!!$    INTEGER :: k, i, j, p
!!$
!!$    CALL smb_2(un,rk)
!!$    CALL divide_by_lumped(rk)
!!$    unext = un+inputs%dt*rk
!!$    velnext = compute_velocity(unext)
!!$
!!$    Ent     = inputs%gravity*un(1,:)**2
!!$    Entnext = inputs%gravity*unext(1,:)**2
!!$    DO k = 1, k_dim
!!$       Ent     = Ent     + velocity(k,:)*un(k+1,:)
!!$       Entnext = Entnext + velnext(k,:)*unext(k+1,:)
!!$    END DO
!!$    Ent     = 0.5d0*Ent
!!$    Entnext = 0.5d0*Entnext
!!$
!!$    Entprime(1,:) = inputs%gravity*un(1,:)  ! -|u|^2/2+gh
!!$    DO k = 1, k_dim
!!$       Entprime(1,:) = Entprime(1,:) -0.5d0*velocity(k,:)**2
!!$       Entprime(k+1,:) = velocity(k,:)
!!$    END DO
!!$
!!$    res = lumped*(Entnext- Ent)/inputs%dt
!!$    DO k = 1, inputs%syst_size
!!$       res = res -lumped*rk(k,:)*Entprime(k,:)
!!$    END DO
!!$
!!$    CALL maxmin(ent,dij,maxn,minn)
!!$    rescale = ABS(maxn-minn)/2 +1.d-10*inputs%gravity*max_water_h**2
!!$    res = ABS(res)/rescale
!!$
!!$    resij%aa = 0.d0
!!$    DO i = 1, mesh%np
!!$       DO p = resij%ia(i), resij%ia(i+1) - 1
!!$          j = resij%ja(p)
!!$          resij%aa(p) = max(res(i),res(j))
!!$       END DO
!!$    END DO
!!$  END SUBROUTINE entropy_residual

  SUBROUTINE entropy_commutator(un)
    USE mesh_handling
    !USE pardiso_solve
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)       :: un
    REAL(KIND=8), DIMENSION(mesh%np)                        :: scal, res, ent, rescale, minn, maxn
    REAL(KIND=8), DIMENSION(inputs%syst_size,k_dim,mesh%np) :: ff
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)       :: Entprime
    REAL(KIND=8), DIMENSION(k_dim,mesh%np)                  :: ent_flux
    REAL(KIND=8) :: xx, yy
    INTEGER :: k, i, j, p

    scal = inputs%gravity*un(1,:)**2
    !scal = scal + inputs%gravity*un(1,:)*bath !===eta+ghz ***(entropy with bathymetry)
    DO k = 1, k_dim
       scal = scal + 0.5d0*velocity(k,:)*un(k+1,:)
    END DO
    ent =  scal - 0.5d0*inputs%gravity*un(1,:)**2 !===|v|^2 h/2 + g h^2/2
    CALL maxmin(ent,dij,maxn,minn)

    DO k = 1, k_dim
       ent_flux(k,:) = velocity(k,:)*scal !===v (v|^2 h/2 + g h^2); includes vhgz if entropy with bathymetry
    END DO

    Entprime(1,:) = inputs%gravity*un(1,:) !===-|v|^2/2 + g h
    !Entprime(1,:) = Entprime(1,:) + inputs%gravity*bath !=== + g z ***(entropy with bathymetry)
    DO k = 1, k_dim
       Entprime(1,:) = Entprime(1,:) -0.5d0*velocity(k,:)**2
       Entprime(k+1,:) = velocity(k,:)  !===v
    END DO

    ff = flux(un)
    DO k = 1, k_dim
       ff(k+1,k,:) = ff(k+1,k,:) + 0.5d0*inputs%gravity*un(1,:)**2
    END DO

    res  = 0.d0
    DO i = 1, mesh%np
       xx=0.d0
       DO p = resij%ia(i), resij%ia(i+1) - 1
          j = resij%ja(p)
          !yy = inputs%gravity*un(1,i)*bath(j) !===g h z ***(entropy with bathymetry)
          DO k = 1, k_dim
             !xx = xx - cij(k)%aa(p)*Entprime(k+1,i)*yy !=== h_i v.grad(z) ***(entropy with bathymetry)
             xx = xx + cij(k)%aa(p)*(ent_flux(k,j)-SUM(Entprime(:,i)*ff(:,k,j)))
          END DO
       END DO
       res(i) = res(i) + xx
    END DO

    !rescale =MAX(ABS(maxn-minn)/2, inputs%gravity*inputs%htiny**2)
    rescale =MAX(ABS(maxn-minn)/2,inputs%epsilon_htiny*maxn)
    res = ABS(res)/rescale

    resij%aa = 0.d0
    DO i = 1, mesh%np
       DO p = resij%ia(i), resij%ia(i+1) - 1
          j = resij%ja(p)
          resij%aa(p) = max(res(i),res(j))
       END DO
    END DO
  END SUBROUTINE entropy_commutator

  SUBROUTINE compute_betaij
    USE mesh_handling
    IMPLICIT NONE
    INTEGER :: m, ni, n1, n2, i, j, i1, i2, p
    REAL(KIND=8) :: d1, d2, alpha, scal, x
    REAL(KIND=8), DIMENSION(2) :: xi, x1, x2
    betaij%aa = 0.d0
    DO m = 1, mesh%me
       DO ni = 1, 3
          n1 = MODULO(ni,3)+1
          n2 = MODULO(ni+1,3)+1
          i  = mesh%jj(ni,m)
          i1 = mesh%jj(n1,m)
          i2 = mesh%jj(n2,m)
          xi = mesh%rr(:,i)
          x1 = mesh%rr(:,i1)-xi
          x2 = mesh%rr(:,i2)-xi
          d1 = SQRT(SUM(x1**2))
          d2 = SQRT(SUM(x2**2))
          scal = SUM(x1*x2)
          alpha = ACOS(scal/(d1*d2))
          DO p = betaij%ia(i), betaij%ia(i+1) - 1
             IF (betaij%ja(p)==i1) THEN
                betaij%aa(p) = betaij%aa(p) + TAN(alpha/2)/d1
             ELSE IF (betaij%ja(p)==i2) THEN
                betaij%aa(p) = betaij%aa(p) + TAN(alpha/2)/d2
             END IF
          END DO
       END DO
    END DO

    DO i = 1, mesh%np
       x = 0.d0
       DO p = betaij%ia(i), betaij%ia(i+1) - 1
          j = betaij%ja(p)
          if ( betaij%aa(p)==0.d0 .AND. i.ne.j) THEN
             write(*,*) ' BUG', betaij%aa(p), p, i, j
          end if
          x = x + betaij%aa(p)*(mesh%rr(1,j) - mesh%rr(1,i))
       END DO
    END DO
  END SUBROUTINE compute_betaij

  SUBROUTINE maxmin(un,mat,maxn,minn)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: un
    TYPE(matrice_bloc),         INTENT(IN)  :: mat
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: maxn, minn
    REAL(KIND=8), PARAMETER :: pi=ACOS(-1.d0)
    INTEGER      :: i
    DO i = 1, SIZE(un)
       maxn(i) = MAXVAL(un(mat%ja(mat%ia(i):mat%ia(i+1)-1)))
       minn(i) = MINVAL(un(mat%ja(mat%ia(i):mat%ia(i+1)-1)))
    END DO
  END SUBROUTINE maxmin

  SUBROUTINE quadratic_limiting(ulow,cmax,psi_small)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN) :: ulow
    REAL(KIND=8), DIMENSION(mesh%np)                     :: cmax, psi_small
    REAL(KIND=8), DIMENSION(inputs%syst_size)  :: ur, Pij
    REAL(KIND=8) :: lambdai, psir, lr, a, b, c, delta, at, coeff, tm, tp
    INTEGER      :: i, j, p, k
    DO i = 1, mesh%np
       lambdai = 1.d0/(mass%ia(i+1) - 1.d0 - mass%ia(i))
       coeff = 1.d0/(lambdai*lumped(i))
       IF (psi_small(i) .LE. 0.d0) THEN
          lij%aa(mass%ia(i):mass%ia(i+1)-1)=0.d0
       END IF
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j =  mass%ja(p)
          IF (i==j) THEN
             lij%aa(p) = 0.d0
             CYCLE
          END IF
          lr = lij%aa(p) !===Use previous limiter
          DO k = 1, inputs%syst_size
             Pij(k) = fctmat(k)%aa(p)*coeff
             ur(k) = ulow(k,i) + lr*Pij(k)
          END DO
          psir = cmax(i)*ur(1) - SUM(ur(2:k_dim)**2)/2
          !IF (psir.GE.-psi_small(i)) THEN
          IF (psir>0.d0) THEN
             CYCLE
          END IF

          a = -SUM(Pij(2:k_dim+1)**2)/2 !===a is negative
          b = cmax(i)*Pij(1) - SUM(Pij(2:k_dim+1)*ulow(2:k_dim+1,i))
          c = cmax(i)*ulow(1,i) - SUM(ulow(2:k_dim+1,i)**2)/2
          delta = b**2 - 4*a*c
          at = min(a,-psi_small(i))
          IF (delta.LE.0.d0) THEN
             tp = -b/(2*at)
             tm = tp
          ELSE
             tp = (-b+sqrt(delta))/(2*at)
             tm = (-b-sqrt(delta))/(2*at)
             if (tp>tm) write(*,*) ' BUG'
          END IF
          !write(*,*) ' tm,tp', tm, tp, delta
          IF (tp>0) THEN
             lij%aa(p) = min(tp,lij%aa(p)) !===0<tp<tm
          ELSE IF (tm>0) THEN
             lij%aa(p) = min(tm,lij%aa(p)) !===tp<0<tm
          ELSE
             !lij%aa(p) = lij%aa(p) !===tp<tm<0
          END IF
          !write(*,*) lij%aa(p)
       END DO
    END DO
  END SUBROUTINE quadratic_limiting

  SUBROUTINE convex_limiting(ulow,cmin,psi_small,psi_func,psi_prime_func)
    USE boundary_conditions
    IMPLICIT NONE
    INTERFACE
       FUNCTION psi_func(u,cmin,Budget) RESULT(psi)
         USE boundary_conditions
         USE space_dim
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(inputs%syst_size) :: u
         REAL(KIND=8)                     :: cmin
         REAL(KIND=8)                     :: Budget
         REAL(KIND=8)                     :: psi
       END FUNCTION psi_func
       FUNCTION psi_prime_func(Pij,u,cmin) RESULT(psi)
         USE boundary_conditions
         USE space_dim
         IMPLICIT NONE
         REAL(KIND=8), DIMENSION(inputs%syst_size) :: u, Pij
         REAL(KIND=8)                     :: cmin
         REAL(KIND=8)                     :: psi
       END FUNCTION psi_prime_func
    END INTERFACE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN) :: ulow
    REAL(KIND=8), DIMENSION(mesh%np)                     :: cmin, psi_small
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np) :: test
    REAL(KIND=8), DIMENSION(mesh%np)  :: centrop
    REAL(KIND=8), DIMENSION(inputs%syst_size)  :: ul, ur, Pij, x
    REAL(KIND=8) :: lambdai, coeff, psir, psil, ll, lr, llold, lrold, Qplus, dQplus, Budget
    INTEGER      :: i, j, p, k, Card, it
    DO i = 1, mesh%np
       lambdai = 1.d0/(mass%ia(i+1) - 1.d0 - mass%ia(i))
       coeff = 1.d0/(lambdai*lumped(i))
       !===Budget
!!$       Qplus = 0.d0
!!$       Card  = 0
!!$       DO p = mass%ia(i), mass%ia(i+1) - 1
!!$          j = mass%ja(p)
!!$          IF (i==j) CYCLE
!!$          ul = ulow(:,i)
!!$          DO k = 1 , inputs%syst_size
!!$             Pij(k) = fctmat(k)%aa(p)*coeff
!!$             ur(k) = ulow(k,i) + lij%aa(p)*Pij(k) !===Density must be positive
!!$          END DO
!!$          dQplus = MIN(psi_func(ul,cmin(i),0.d0),psi_func(ur,cmin(i),0.d0))
!!$          IF (dQplus>0.d0) THEN
!!$             Qplus = Qplus + dQplus
!!$          ELSE
!!$             Card  = Card + 1
!!$          END IF
!!$       END DO
!!$       IF (Card.NE.0) THEN
!!$          Budget = -Qplus/Card
!!$       ELSE
!!$          Budget = -1d15*Qplus
!!$       END IF
       Budget =0.d0
       !===End Budget
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j =  mass%ja(p)
          IF (i==j) THEN
             lij%aa(p) = 0.d0
             CYCLE
          END IF
          lr = lij%aa(p) !===Use previous limiter
          DO k = 1, inputs%syst_size
             Pij(k) = fctmat(k)%aa(p)*coeff
             ur(k) = ulow(k,i) + lr*Pij(k)
          END DO
          psir = psi_func(ur,cmin(i),Budget)
          IF (psir.GE.-psi_small(i)) THEN
             CYCLE
          END IF
          ll = 0.d0
          ul = ulow(:,i)
          !psil = max(psi_func(ul,cmin(i),Budget),0.d0) !===To avoid roundoff negative
          psil = psi_func(ul,cmin(i),Budget)
          DO WHILE (ABS(psil-psir) .GT. psi_small(i))
             llold = ll
             lrold = lr
             ll = ll - psil*(lr-ll)/(psir-psil)
             lr = lr - psir/psi_prime_func(Pij,ur,cmin(i))
             if (ll< llold) then
                write(*,*) ' f1', ll , llold, psil, psir,psi_small(i)
                !stop
                ll = llold
                EXIT
             end if
             if (lr > lrold) then
                write(*,*) ' f2', lr , lrold, psil, psir,psi_small(i),psi_prime_func(Pij,ur,cmin(i))
                !write(*,*), 'cmin', cmin(i)
                !stop
                lr = lrold
                EXIT
             end if
             IF (ll.GE.lr) THEN
                ll = lr
                EXIT
             END IF
             ul = ulow(:,i) + ll*Pij
             ur = ulow(:,i) + lr*Pij
             !psil = max(psi_func(ul,cmin(i),Budget),0.d0) !===To avoid roundoff negative
             psil = psi_func(ul,cmin(i),Budget)
             psir = psi_func(ur,cmin(i),Budget)
          END DO
          IF (psir.GE.-psi_small(i)) THEN
             lij%aa(p) = lr
          ELSE
             lij%aa(p) = ll
          END IF
       END DO
       write(*,*) psil, psir
    END DO

  CONTAINS
  END SUBROUTINE convex_limiting

  FUNCTION psi_rho_e(u,emax,Budget) RESULT(psi)
    USE boundary_conditions
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size) :: u
    REAL(KIND=8)                     :: emax
    REAL(KIND=8)                     :: Budget
    REAL(KIND=8)                     :: psi
    psi = emax*u(1) - SUM(u(2:k_dim+1)**2)/(2.d0)  - Budget
  END FUNCTION psi_rho_e

  FUNCTION psi_rho_e_prime(Pij,u,emax) RESULT(psi)
    USE boundary_conditions
    USE space_dim
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size) :: u, Pij
    REAL(KIND=8)                     :: emax
    REAL(KIND=8)                     :: psi
    REAL(KIND=8) :: rho_e
    psi = emax*Pij(1) - SUM(Pij(2:k_dim+1)*u(2:k_dim+1))
  END FUNCTION psi_rho_e_prime

  SUBROUTINE estimate_rho_e_max(un,rho_e_max)
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(IN) :: un
    REAL(KIND=8), DIMENSION(mesh%np)                     :: rho_e_max
    REAL(KIND=8), DIMENSION(mesh%np)                     :: rho_e
    INTEGER :: i
    DO i = 1, mesh%np
       rho_e(i) = un(1,i)*SUM(velocity(:,i)**2)/2
    END DO
    DO i = 1, mesh%np
       rho_e_max(i)  = MAXVAL(rho_e(mass%ja(mass%ia(i):mass%ia(i+1)-1)))
    END DO
  END SUBROUTINE estimate_rho_e_max

  SUBROUTINE relax(un,minn,maxn)
    USE mesh_handling
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:)              :: un
    REAL(KIND=8), DIMENSION(:)              :: minn
    REAL(KIND=8), DIMENSION(:)              :: maxn
    REAL(KIND=8), DIMENSION(SIZE(un))       :: alpha, denom
    INTEGER      :: i, j, p, ps, pe
    REAL(KIND=8) :: x, mx, mn
    alpha = 0.d0
    denom = 1.d-14*MAX(ABS(maxn),ABS(minn))
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          IF (i==j) CYCLE
          alpha(i) = alpha(i) + un(i) - un(j)
          denom(i) = denom(i) + ABS(un(i)) + ABS(un(j))
       END DO
    END DO

    IF(inputs%limiter_type=='avg') THEN !===Average
       denom = 0.d0
       DO i = 1, SIZE(un)
          DO p = mass%ia(i), mass%ia(i+1) - 1
             j = mass%ja(p)
             IF (i==j) CYCLE
             denom(i) = denom(i) + alpha(j) + alpha(i)
          END DO
       END DO
       DO i = 1, SIZE(un)
          alpha(i) = denom(i)/(mass%ia(i+1)-mass%ia(i))/2
       END DO
       maxn = MIN(1.01*maxn,maxn + abs(alpha)/2)
       minn = MAX(0.99*minn,minn - abs(alpha)/2)
    ELSE IF(inputs%limiter_type=='minmod') THEN !===Minmod
       denom = alpha
       DO i = 1, SIZE(un)
          DO p = mass%ia(i), mass%ia(i+1) - 1
             j = mass%ja(p)
             IF (i==j) CYCLE
             IF (denom(i)*alpha(j).LE.0.d0) THEN
                denom(i) = 0.d0
             ELSE IF (ABS(denom(i)) > ABS(alpha(j))) THEN
                denom(i) = alpha(j)
             END IF
          END DO
       END DO
       alpha = denom
       maxn = MIN(1.01*maxn,maxn + abs(alpha)/2)
       minn = MAX(0.99*minn,minn - abs(alpha)/2)
    END IF
  END SUBROUTINE RELAX

!===NEW FORMULATION
    SUBROUTINE smb_2(un,rk)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(OUT) :: rk
    REAL(KIND=8), DIMENSION(inputs%syst_size,k_dim,mesh%np)        :: vv
    INTEGER :: d, i, j, k, p
    REAL(KIND=8) :: xx, hmean, Hstarij, Hstarji, ratij, ratji
    vv=flux(un)
    rk=0.d0
    fix_roundoff = 0.d0
    DO i = 1, mesh%np
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          Hstarij = MAX(0.d0,un(1,i)+bath(i)-MAX(bath(i),bath(j)))
          Hstarji = MAX(0.d0,un(1,j)+bath(j)-MAX(bath(i),bath(j)))

          DO k = 1, inputs%syst_size
             xx = 0.d0
             DO d = 1, k_dim
                xx = xx - cij(d)%aa(p)*(vv(k,d,j))
             END DO
             rk(k,i) = rk(k,i) + xx + (dij%aa(p)-muij%aa(p))*(un_over_h(k,j)*Hstarji-un_over_h(k,i)*Hstarij) &
                  + muij%aa(p)*(un(k,j)-un(k,i))
          END DO

          !===Topography contribution
          DO k = 1, k_dim
             rk(k+1,i) = rk(k+1,i) &
                  - inputs%gravity*un(1,i)*(un(1,j)+bath(j))*cij(k)%aa(p)
          END DO
       END DO
    END DO

    SELECT CASE(inputs%type_test)
    CASE(9,10,11,12,13,15)
       CALL friction(un,rk)
    END SELECT

    IF (inputs%if_FGN) THEN
       CALL FGN_rhs(un,rk)
    END IF

  END SUBROUTINE smb_2

  SUBROUTINE FGN_rhs(un,rk)
    USE mesh_handling
    USE boundary_conditions
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np)  :: un
    REAL(KIND=8), DIMENSION(inputs%syst_size,mesh%np), INTENT(OUT) :: rk
    REAL(KIND=8), DIMENSION(mesh%np)  :: lambda_bar,  eta, omega, &
         htwo_Gammap, pp, alpha, eta_over_h, heta_gammap, h, ratio
    REAL(KIND=8) :: mineta, x0
    INTEGER :: i, j, k, p

    h  = un(1,:)
    eta = un(inputs%syst_size-1,:)*compute_one_over_h(un(1,:))
    eta_over_h = un(inputs%syst_size-1,:)*compute_one_over_h(un(1,:))**2
    omega = un(inputs%syst_size,:)*compute_one_over_h(un(1,:))
    ratio =  2.d0*un(inputs%syst_size-1,:)/(eta**2+h**2+inputs%htiny)
    IF (MINVAL(eta)<-1.d-14*inputs%max_water_h) THEN
       WRITE(*,*) 'eta negative', MINVAL(eta)
       !STOP
    END IF

    alpha = inputs%lambda_bar/(3*lumped) !1D

    DO i = 1, mesh%np
       x0 =  min(2.d0,SQRT(1.d0+1.d0/(2*alpha(i)*(max(eta(i),0.d0)+inputs%htiny))))
       IF (eta(i).LE.0.d0) THEN
          pp(i) = 0.d0
          htwo_Gammap(i) = 0.d0
          heta_Gammap(i) = 0.d0
       ELSE IF (eta(i).LE.x0*h(i)) THEN
          pp(i) = -alpha(i)*inputs%gravity*eta(i)*(eta(i)**2-h(i)**2)
          htwo_Gammap(i) = 3*eta(i)**2+h(i)**2-4*un(inputs%syst_size-1,i)
          heta_Gammap(i) = 3*eta(i)**2*eta_over_h(i)+un(inputs%syst_size-1,i)-4*eta(i)**2
       ELSE
          pp(i) = -alpha(i)*inputs%gravity*(x0**2-1.d0)*un(inputs%syst_size-1,i)*h(i)
          htwo_Gammap(i) = 4*(x0-1.d0)*un(inputs%syst_size-1,i)+(1-x0**2)*h(i)**2
          heta_Gammap(i) = 4*(x0-1.d0)*eta(i)**2 + (1-x0**2)*un(inputs%syst_size-1,i)
       END IF
    END DO

    IF (k_dim==1) THEN
       lambda_bar = inputs%lambda_bar*inputs%gravity/lumped

    ELSE
       lambda_bar = inputs%lambda_bar*inputs%gravity/SQRT(lumped)
    END IF

    DO i = 1, mesh%np
       !rk(inputs%syst_size-1,i) = rk(inputs%syst_size-1,i) + lumped(i)*eta(i)*omega(i)*ratio(i)**3
       !rk(inputs%syst_size,i)   = rk(inputs%syst_size,i)   - lambda_bar(i)*lumped(i)*heta_Gammap(i)*ratio(i)
       !rk(inputs%syst_size-1,i) = rk(inputs%syst_size-1,i) + lumped(i)*un(inputs%syst_size,i)
       !rk(inputs%syst_size,i)   = rk(inputs%syst_size,i)   - lambda_bar(i)*lumped(i)*htwo_Gammap(i)
       rk(inputs%syst_size-1,i) = rk(inputs%syst_size-1,i) + lumped(i)*un(inputs%syst_size,i)*ratio(i)**3
       rk(inputs%syst_size,i)   = rk(inputs%syst_size,i)   - lambda_bar(i)*lumped(i)*htwo_Gammap(i)*ratio(i)**3
       DO p = mass%ia(i), mass%ia(i+1) - 1
          j = mass%ja(p)
          DO k = 1, k_dim
             rk(k+1,i) = rk(k+1,i) - pp(j)*cij(k)%aa(p)
          END DO
       END DO
    END DO

  END SUBROUTINE FGN_rhs

END MODULE update
