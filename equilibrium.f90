! tauss [model] -- taubal [equilibrium]
! rd    [model] -- rbar   [equilibrium]
module equilibrium
    !use universe
    use model
    implicit none ! August 27, 2016
    contains
    !subroutine search_equilibrium(exit_log1)
    subroutine search_equilibrium( guessv, momvec, obj_val, node_id, trial_id, model_msg) ! node_id == my_id
        implicit none
        real(wp), intent(out) :: obj_val
        real(wp), dimension(:), intent(in) :: guessv
        real(wp), dimension(:), intent(out) :: momvec
        integer, intent(in) :: node_id, trial_id
        integer, intent(inout) :: model_msg
        integer :: i
        !logical, intent(out) :: exit_log1
        logical :: exit_log1
        character(len=100) :: msg
        ! Used for extension 04:00pm.
        
        
        !if(my_id==0)then
        !    do i = 1, size(para)
        !        write(unit=120,fmt='(i4,e12.4,2x,a)') i, para(i), labstr(i)
        !    enddo
        !enddo
        
        model_msg = 0 ! default value set before calling the subroutine 'solve_model'.
        exit_log1 = .false. ! Initializing Exit Flag. If successful, .False.; Otherwise, .True.
        
        if(printout11)then ! .true. for running the real quantitative model
            !if(mpi_exercise_mode/=5)then
                call solve_model( guessv, momvec, node_id, trial_id, exit_log1, msg) ! defined in this module
            !else
            !    ! Make subdim dimensions be extended to ndim dimension so that we can use the usual routine to obtain the usual ndim moments instead of just a subset of the original moments.
            !    ! One rule to reach hopefully is that solve_model doesn't need to change anything. That is, its dimension remains 10 dimensions.
            !    
            !endif
        else ! .false. for testing the communication of the coarse search.
            call test_model( guessv, momvec, node_id, trial_id, exit_log1, msg) 
        endif
        
        if(exit_log1 == .true.)then ! Failure to solve the model with the given parameter setting.
            !write(*,fmt='(2/,a,a,a,/)'), ' Something wrong [', trim(adjustl(msg)), '] ' ! 4.23.2017 It needs to be turned off in production mode.
            model_msg = 1 ! Fail to solve the model.  
            momvec    = inf
            obj_val   = inf
        else
            if(printout11)then
                obj_val = objective_value(momvec,targetv) ! momvec: simulatd moments    
            else
                obj_val = test_objective_value(momvec,targetv)
            endif
        endif
    end subroutine search_equilibrium
    
    subroutine maximum_distance_vertices( mat_vertices, size_of_worst_vertices, ans_vec ) ! 7-30-2017 checked.
        real(wp), dimension(:,:), intent(in) :: mat_vertices
        real(wp), dimension(:), intent(out) :: ans_vec
        integer, intent(in) :: size_of_worst_vertices
        !real(wp) :: maximum_distance_vertices
        real(wp), dimension(:), allocatable :: centroid_vec, dist_vec
        integer :: m, n, i
        m = size(mat_vertices,dim=1) ! row
        n = size(mat_vertices,dim=2) ! col
        allocate( centroid_vec(m), dist_vec(n) )
        !! Method 1.
        !do i = 1, m
        !    centroid_vec(i) = sum(mat_vertices(i,1:n-size_of_worst_vertices))/(n-size_of_worst_vertices)        
        !enddo ! i
        !do i = 1, n
        !    dist_vec(i) = (sum((mat_vertices(:,i)-centroid_vec)**2._wp))**0.5_wp
        !enddo
        !maximum_distance_vertices = maxval(dist_vec)
        
        ! Method 2.
        do i = 1, n ! Note that it goes across the vertices including the best vertex itself.
            dist_vec(i) = sum((mat_vertices(:,i)-mat_vertices(:,1))**2._wp)**0.5_wp
        enddo ! i
        !maximum_distance_vertices = maxval(dist_vec(2:n))
        ans_vec = dist_vec
        deallocate( centroid_vec, dist_vec )
    end subroutine maximum_distance_vertices
    
    subroutine maximum_penalty_distance( vec, ans )
        implicit none
        real(wp), dimension(:), intent(in) :: vec
        real(wp), intent(out) :: ans
        real(wp) :: big, sml
        big = maxval(vec)
        sml = minval(vec)
        ans = abs(big-sml)
    end subroutine 
    
    !subroutine test2()
    !    implicit none
    !    print*, '-------------- hahaha ', obj_func_toggle
    !end subroutine test2
    
    function objective_value( mom, tar )
        real(wp) :: objective_value
        real(wp), dimension(:), intent(in) :: mom, tar
        select case(obj_func_toggle)
            case(1) ! 10 moments
                objective_value = sum(((mom-tar)/tar)**2._wp)
            case(2) ! 5 moments
                objective_value = ((mom(2)-tar(2))/tar(2))**2._wp + ((mom(3)-tar(3))/tar(3))**2._wp + ((mom(5)-tar(5))/tar(5))**2._wp + ((mom(6)-tar(6))/tar(6))**2._wp + ((mom(10)-tar(10))/tar(10))**2._wp
            case(3) ! 8 moments
                objective_value = ((mom(1)-tar(1))/tar(1))**2._wp + ((mom(4)-tar(4))/tar(4))**2._wp + ((mom(5)-tar(5))/tar(5))**2._wp + ((mom(6)-tar(6))/tar(6))**2._wp + ((mom(7)-tar(7))/tar(7))**2._wp + ((mom(8)-tar(8))/tar(8))**2._wp + ((mom(9)-tar(9))/tar(9))**2._wp + ((mom(10)-tar(10))/tar(10))**2._wp
            case(4) ! 6 moments
                objective_value = ((mom(1)-tar(1))/tar(1))**2._wp + ((mom(2)-tar(2))/tar(2))**2._wp + ((mom(5)-tar(5))/tar(5))**2._wp + ((mom(6)-tar(6))/tar(6))**2._wp + ((mom(7)-tar(7))/tar(7))**2._wp + ((mom(9)-tar(9))/tar(9))**2._wp
        end select 
        !objective_value = sum( (mom-tar)**2._wp )
    end function objective_value   
    
    function test_objective_value( mom, tar )
        real(wp) :: test_objective_value
        real(wp), dimension(:), intent(in) :: mom, tar
        real(wp), dimension(:), allocatable :: func_ray, fund_ray
        real(wp) :: temp1, temp2
        integer :: i, j, m
        select case(testfunc_idx)
            case(1) ! Ordinary norm 2 measure
                test_objective_value = sum( (mom)**2._wp )
            case(2) ! Variably Dimensioned Function
                allocate( func_ray(ndim), fund_ray(ndim) )
                do i = 1, ndim
                    func_ray(i) = mom(i)-1._wp
                    fund_ray(i) = i*func_ray(i)
                enddo ! i              
                temp1 = sum(fund_ray)
                temp2 = temp1**2._wp
                test_objective_value = sum(func_ray**2._wp) + temp1**2._wp + temp2**2._wp + 1
                deallocate( func_ray, fund_ray )
            case(3) ! only the first five dimension using the same objective function as in case 2.
                allocate( func_ray(subdim), fund_ray(subdim) )
                do i = 1, subdim
                    func_ray(i) = mom(i)-1._wp
                    fund_ray(i) = i*func_ray(i)
                enddo ! i              
                temp1 = sum(fund_ray)
                temp2 = temp1**2._wp
                test_objective_value = sum(func_ray**2._wp) + temp1**2._wp + temp2**2._wp + 1
                deallocate( func_ray, fund_ray )
            case(4) ! Extended Rosenbrock Function
                if(mod(ndim,2)/=0) write(*,*) " Extended rosenbrock function requires the size of dimensions is an even number. "
                allocate( func_ray(ndim) )
                do i = 1, ndim
                    if(mod(i,2)==0)then
                        func_ray(i) = 1._wp - mom(i-1)    
                    else
                        func_ray(i) = 10._wp*(mom(i+1)-mom(i)**2._wp)
                    endif
                enddo
                test_objective_value = sum(func_ray**2._wp) + 1._wp
                deallocate( func_ray )                
            case(5) ! Extended Powell Singular Function
                
                !write(*,*) " Function value takes non-zero only in the first eight dimensions  "
                allocate( func_ray(ndim) )
                do i = 1, ndim
                    j = mod(i,4)        
                    m = merge((i-j)/4,(i-j)/4+1,j==0)
                    if(m<=2)then
                        select case (j)
                        case (1)
                            func_ray(i) = mom(i) + 10._wp*mom(i+1) - 11._wp    
                        case (2)
                            func_ray(i) = (5._wp)**0.5_wp*(mom(i+1)-mom(i+2))
                        case (3)
                            func_ray(i) = (mom(i-1)-2._wp*mom(i) +1._wp)**2._wp
                        case (0) 
                            func_ray(i) = (10)**0.5_wp*(mom(i-3)-mom(i))**2._wp
                        end select
                    else
                        func_ray(i) = mom(i)
                    endif    
                enddo ! i 
                test_objective_value = sum(func_ray**2._wp) + 1._wp
                deallocate( func_ray )
                
        end select     
        !momvec = guessv*(sin(guessv)+0.1_wp) ! 7-21-2017, global optimum: f(x_i) = 0 for x_i = 0 for i = 1, ..., n, in [-10,10]
        !test_objective_value = sum(abs(mom*(sin(mom)+0.1_wp)))
    end function test_objective_value       
    
    subroutine test_model( guessv, momvec, node_id, trial_id, exit_log1, msg )
        implicit none
        real(wp), dimension(:), intent(in) :: guessv
        real(wp), dimension(:), intent(out) :: momvec
        integer, intent(in) :: node_id, trial_id
        logical, intent(inout) :: exit_log1
        character(len=*), intent(out) :: msg
        
        !! Random number generation
        !integer :: generator, erridx, approach, tmiddle, tstart, tend, trate, tmax
        !type(vsl_stream_state) :: river
        !real(wp) :: stoch(10)
        !call system_clock(tstart,trate,tmax)
        !generator = VSL_BRNG_MCG31
        !approach = VSL_RNG_METHOD_UNIFORM_STD
        !call system_clock(tmiddle)
        !erridx = vslnewstream(river, generator, tmiddle)
        !write(*,*) ' tmiddle ', my_id, tmiddle
        !erridx = vdrnguniform(approach, river, 10, stoch, 0._wp, 1._wp)
        !momvec = stoch
        
        !momvec = guessv**2._wp
        !write(unit=my_id+1001,fmt='(a,<ndim>f16.7)') ' guessv ', guessv ! test model
        !write(unit=my_id+1001,fmt='(a,<ndim>f16.7)') ' momvec ', momvec ! test model
        !msg = 'random number'
        
        !momvec = guessv*(sin(guessv)+0.1_wp) ! 7-21-2017, global optimum: f(x_i) = 0 for x_i = 0 for i = 1, ..., n, in [-10,10]
        
        momvec = guessv
        msg = "ok"
        
        if(printout15)then
            if( guessv(2)>3._wp .or. guessv(3)>3._wp )then
                momvec(2) = inf
                exit_log1 = .true.
            endif
        endif
        
        ! 7-21-2017 Set target_i = 0, penalty function as summation of f(x_i) and let the range to be [-10,10].
        
    endsubroutine test_model
    
    subroutine solve_model( guessv, momvec, node_id, trial_id, exit_log1, msg) ! node_id == my_id
        implicit none
        real(wp), dimension(:), intent(in) :: guessv
        real(wp), dimension(:), intent(out) :: momvec
        logical, intent(inout) :: exit_log1 ! EXIT "SOLVE MODEL" SUBROUTINE; It is initialized by its calling subroutine.
        character(len=*), intent(out) :: msg ! MESSAGE FOR EXIT THIS SUBROUTINE
        integer, intent(in) :: node_id, trial_id
        
        integer :: i, j, m, n, ti, hi, ki, yi, kpi, int1
        real(wp) :: sar, tar ! THESE DUMMY VARIABLES USED FOR MULTIPLE TESTING BLOCKS, SO DON'T UNCOMMENT MORE THAN ONE BLOCKS AT THE SAME TIME.
        real(wp), dimension(:), allocatable :: yvtemp
        integer :: colbnd1, colbnd2 ! FOR ILLUSTRATION        
        integer :: tstart, tend, trate, tmax        
        logical :: exit_bellman
        
        !exit_log1 = .false. ! Redundant. Was initialized in the evoking subroutine 'search_equilibrium' 
        
        kv1   = guessv(1)
        prtk0 = guessv(2)
        prtk1 = guessv(3)
        prtk2 = guessv(4)
        zbar  = guessv(5)
        beta  = guessv(6)
        theta = guessv(7)
        !iota  = guessv(7)
        phi1  = guessv(8)
        phi2  = guessv(9)
        phi3  = guessv(10)
        
        msg = " "
        exit_bellman = .false.
        
        !write(4000+trial_id,fmt='("initial",8x,10(e21.14,x))')  kv1, prtk0, prtk1, prtk2, zbar, beta, theta, phi1, phi2, phi3
        
        ! 4.1.2017 fnadim, fnhdim, adim, hdim all need to keep. Hard to disentangle these variables, although their functions are overlapping each other's.
        fnadim = adim ! 3.16.2017 Hardwired for preventing program complexity. 
        fnhdim = hdim ! 3.16.2017 Hardwired for preventing program complexity. 
        
        allocate( py(nmc,nmc), yvtemp(nmc), yv(0:nmc), sy(nmc), pyh(nmc,nmc), yhv(nmc), syh(nmc), survprob(14) )
        allocate( popfrac(14), kv(0:kdim-1), probtk(0:kdim-1), phi(0:kdim-1), z2(0:kdim-1), delzh(0:kdim-1), delzl(0:kdim-1) )
        allocate( efflab(14), delmeh(3), hv(hdim), av(adim), pz2(0:kdim-1,2), pka(0:kdim-1,2), rhv(fnhdim), rav(fnadim), mass_vec(adim), h_mass_vec(hdim) ) ! pt18(kdim*nmc,2*nmc) ! 4.14.2017 pz2 set the index staring from 0.
        allocate( collbdmat(1:hdim,1:(kdim-1),0:nmc,1:14), c_prf_mat(adim,0:(kdim-1),0:1,0:nmc,1:14), c_lab_mat(adim,0:(kdim-1),0:1,0:nmc,1:14) ) 
        allocate( c_grs_mat(adim,0:(kdim-1),0:1,0:nmc,1:14) ) ! before taxes business income
        
        ! a, h, k, z, y, kp, yp, op, t
        allocate( twa(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  twh(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  !vtwh(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  !vtwa(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  twk(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  !vtwk(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  twf(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  tww(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  cef(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  dcef(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &            
                  !scef(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), & 
                  cww(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  cwf(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  cwa(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &        
                  cwh(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  cwk(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  cwc(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  cww2(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  cwf2(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  cwa2(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &        
                  cwc2(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &   
                  cwh2(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  cwk2(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), & 
                  cvv(adim,hdim,0:(kdim-1),0:1,0:nmc,1:14), &
                  ced(adim,hdim,0:(kdim-1),0:1,0:nmc,1:14), &
                  sw_laborsupply(adim*hdim*1018), &
                  sw_buzcap_notuse(adim*hdim*1018), &
                  sw_labordemand(adim*hdim*1018), &
                  sw_production(adim*hdim*1018), &
                  sw_bizinvestment(adim*hdim*1018), &
                  sw_bizloan(adim*hdim*1018), &
                  sw_ini_asset(adim*hdim*1018), &
                  sw_ini_house(adim*hdim*1018), &
                  sw_nonlineartax(adim*hdim*1018), &
                  sw_worker_turned(adim*hdim*1018), &
                  sw_boss_turned(adim*hdim*1018), &
                  sw_aftertaxwealth(adim*hdim*1018), &
                  sw_taxableincome(adim*hdim*1018), &
                  sw_consumption(adim*hdim*1018), &
                  sw_socialsecurity(adim*hdim*1018), &
                  sw_worker_savtax(adim*hdim*1018), & 
                  sw_entpre_savtax(adim*hdim*1018), &
                  sw_entpre_biztax(adim*hdim*1018), &
                  sw_totinc_bx(adim*hdim*1018), &
                  sw_wealth_tax(adim*hdim*1018), &
                  !tid(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  wf(adim,hdim,0:(kdim-1),0:1,0:nmc,1:14), &
                  !ww(adim,hdim,0:(kdim-1),0:1,0:nmc,1:14), &
                  id(fnadim,hdim,0:(kdim-1),0:1,0:nmc,1:14), &  ! invariant distribution (beginning of period)
                  ppldie(fnadim,hdim,0:(kdim-1),0:1,0:nmc,1:14), &
                  id1(fnadim,hdim,0:(kdim-1),0:1,0:nmc), &
                  sid(fnadim,hdim,0:(kdim-1),0:1,0:nmc,1:14), & ! pseudo beginning of period invariant distribution
                  !sww(fnadim,hdim,0:(kdim-1),0:1,0:nmc,1:14), &
                  !afv(fnadim), & 
                  !afint(fnadim,2), & ! wint(fnadim,2), 
                  beqdist(fnadim,hdim,0:(kdim-1),0:1,0:nmc,1:14), &
                  nsav(fnadim,hdim,0:(kdim-1),0:1,0:nmc,1:14), &
                  beqdis(fnadim,hdim,0:(kdim-1),0:1,0:nmc), p_homvec(1:2,1:14), &
                  c_lab_vec(kdim-1), c_opt_vec(kdim-1) )
        
        allocate( sww(adim*hdim*1018), swf(adim*hdim*1018), swa(adim*hdim*1018), swh(adim*hdim*1018), swc(adim*hdim*1018), swk(adim*hdim*1018), sef(adim*hdim*1018), def(adim*hdim*1018) )
        allocate( sww2(adim*hdim*1018), swf2(adim*hdim*1018), swc2(adim*hdim*1018), swa2(adim*hdim*1018), swh2(adim*hdim*1018), swk2(adim*hdim*1018) )
        allocate( uwh(fnadim,fnhdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  uwk(fnadim,fnhdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  uww(fnadim,fnhdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  uwa(fnadim,fnhdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14) )
                          
        allocate( tbim(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  atwm(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), & 
                  taxm(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), & 
                  beqm(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  prom(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  ldemm(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), & 
                  lsupm(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  ent2wok(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  wok2ent(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  cspm(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  homm(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  tbi_m(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  atw_m(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  tax_m(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  beq_m(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  pro_m(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  ldem_m(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  lsup_m(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  hom_m(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  csp_m(adim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  ide(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  !nafv(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), & ! the end of period invariant distribution
                  totast_m(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  entcap_m(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  entlab_m(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  entprd_m(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  ttax_m(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  home_m(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  uwa_m(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), &
                  uwh_m(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14) )
                  ! entsize_m(fnadim,hdim,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14) )
        
        ! end period state combination ---------------------------------- 
        ! 1: local order index
        ! 2: kpx
        ! 3: ypx
        ! 4: opx
        ! 5: macro order index
        allocate(   t18vec1(6,5),   t18vec24(6,5),   t18vec56(6,5),   t18vec7(3,5) )
		allocate(    t9vec1(2,5),    t9vec24(2,5),    t9vec56(2,5),    t9vec7(1,5) )
		allocate( t1013vec1(1,5), t1013vec24(1,5), t1013vec56(2,5), t1013vec7(1,5) )  
        allocate(   t14vec1(1,5),   t14vec24(1,5),   t14vec56(2,5),   t14vec7(1,5) ) 
        ! initialized in the end period state combination subroutine
        allocate( tvector(7,3) )
        
        ! indices list for "coarse" decision ------------------------------------------
        allocate( xl14(adim*1*7,8), bl1014(hdim*7,4), bd14(1,2), bizmat(adim,0:kdim-1,0:1) )   
        allocate( xl1013(adim*1*12,8), bldum(1,1), bd1013(4,2) )
        !allocate( xl9(adim*1*nmc*13,8), bd9(4,2) ) 
        allocate( xl9(adim*1*nmc*16,8), bd9(4,2) ) 
        !allocate( xl18(adim*1*(nmc)**2*13,8), bd18(4,2) )
        allocate( xl18(adim*1*(nmc)**2*16,8), bd18(4,2) ) ! 3.29.2017
        
        allocate( term_2(adim*hdim*1018), term_3(adim*hdim*1018), term_4(adim*hdim*1018), term_5(adim*hdim*1018), term_6(adim*hdim*1018), term_7(adim*hdim*1018), term_8(adim*hdim*1018), term_9(adim*hdim*1018) ) ! 9-17-2017 

        call system_clock(tstart,trate,tmax)
        
        ! Create lists for parallelism =================================================================================================================
        call system_clock(tstart)         
        ! Index matrices and boundary list (for policy function, containing entrepenrus choose not to run a business this period) 072416, 9102016
        ! Feb 6, 2017
        call create_index_list( xl14, bd14, bl1014, '14','coarse')     
        call create_index_list( xl1013, bd1013, bldum, '10-13','coarse')     
        call create_index_list( xl9, bd9, bldum, '9','coarse')  
        call create_index_list( xl18, bd18, bldum, '18','coarse')  
        call system_clock(tend) 
        !write(*,'(a,f12.4,a)') '1: ',real(tend-tstart,wp)/real(trate,wp), ' seconds'         
        
        ! 3.30.2017 comment out
        !! indices list for "refined" decision    
        !allocate( fxl14(fnadim*hdim*7,8), fbl1014(hdim*7,4), fbd14(1,2) )   
        !allocate( fxl1013(fnadim*hdim*9,8), fbldum(1,1), fbd1013(4,2) )
        !allocate( fxl9(fnadim*hdim*nmc*10,8), fbd9(4,2) ) 
        !allocate( fxl18(fnadim*hdim*(nmc)**2*10,8), fbd18(4,2) )  
        
        call system_clock(tstart)         
        !! refined index matrices (for policy function) 072416
        !call create_index_list( fxl14, fbd14, fbl1014, '14','refined')     
        !call create_index_list( fxl1013, fbd1013, fbldum, '10-13','refined')     
        !call create_index_list( fxl9, fbd9, fbldum, '9','refined')  
        !call create_index_list( fxl18, fbd18, fbldum, '18','refined') 
        
        allocate( s1c(394*adim*hdim,7), c1s(1:adim,1:hdim,0:kdim-1,0:1,0:nmc,1:14) ) ! 3.31.2017 keep it. used for stationary distribution. 4.1.2017 correct the number as 394. Seems useless.
        
        !allocate( s2c(1018*adim*hdim,10), c2s(1:adim,1:hdim,0:kdim-1,0:1,0:nmc,0:kdim-1,0:nmc,0:2,1:14) ) ! 3.15.2017 checked.
        !call serialindices_Map2_coordinates(s2c,c2s,adim,hdim) ! COMBINATION ON COARSE GRID FOR THE UNIQUE ONE DIMENSIONAL SERIES.  
        
        allocate( s3c(1018*fnadim*fnhdim,10), c3s(1:fnadim,1:fnhdim,0:kdim-1,0:1,0:nmc,0:kdim-1,0:nmc,0:2,1:14) ) ! 
        call serialindices_Map2_coordinates(s3c,c3s,fnadim,fnhdim) ! COMBINATION ON REFINED GRID FOR THE UNIUQE ONE DIMENSIONAL SERIES. 
        
        !print*, ' test location ', c3s(10,10,1,1,1,2,3,2,5) ! ok
        call system_clock(tend) 
            
        call system_clock(tstart)         
        ! beginign of computation ===============================================================================================================================
        beta = beta**real(length,wp)
        
        ! ## Probability matrices 072316
        call set_markov_matrix( rhoy, vary, rhoyh, varyh, py, yvtemp, sy, pyh, yhv, syh )
        yv = 0._wp
        yv(1:nmc) = yvtemp 
        
        call set_survival_vector( survprob, '_survival_probability.txt' ) ! call ss(survprob,'survprob')
        survprob(1:9) = 1._wp ! In the baseline model people not retired faces no mortality risk. 10102016
        ! See Social Security Administration 2015, table 4.C6, pg. 4.50. Feb 13, 2017
        ! See the folder: E:\GoogleDrive\R_projects\research\thesis\Efficiency_units 
        
        !write(4000+trial_id,fmt='("#1",13x,2i3,8(e21.14,x))') fnadim, fnhdim, beta, yv, survprob(10:12)
        
        call set_efficiency_untis( efflab, '_efficiency_units.txt' ) ! I use the estimated 17.5 Hansen weight as model labor efficiency of the age group of 20 years ago.
        ! 22.5 Hansen weight as that of the age group of 25 years old. Feb 13, 2017
        ! See the folder: E:\GoogleDrive\R_projects\research\thesis\Efficiency_units
        
        !write(4000+trial_id,fmt='("#2",13x,10(e21.14,x))') efflab(1:10)
        
        ! #########################################################################################
        ! ## Fraction of population # no population growth is considered 072316
        popfrac(1) = 10.e10_wp ! 9-17-2017
        do i = 2, 14
            popfrac(i) = popfrac(i-1)*survprob(i-1) ! note: survprob is conditional probability of death    
        enddo
        popfrac = popfrac/sum(popfrac) ! normalization; the sum of the alive in the beginning of each age cohort    
        !write(4000+trial_id,fmt='("#3",13x,10(e21.14,x))') popfrac(5:14)
        
        ! #########################################################################################
        ! business scale vector 072316 Meh 2005 pdf-13.
        kv(0) = 0._wp
        kv(1) = kv1 ! that is k1 in the math model. to be calibrated # ! <------------1
        kv(2) = kv(1)*10._wp
        kv(3) = kv(1)*100._wp
        
        ! #########################################################################################
        ! innovation (new idea) probability, Pk(k') in the math model, to be calibrated # 072316       
        probtk(0) = prtk0 ! to be calibrated # ! <------------------------------------2
        probtk(1) = prtk1 ! to be calibrated # ! <------------------------------------3
        probtk(2) = prtk2 ! to be calibrated # ! <------------------------------------4
        probtk(3) = prtk3 ! set to zero 
        
        ! periods 1 - 8 transition matrix 072316 3.12.2017 To be calibated.
        pka(:,2) = probtk ! probability of acquiring new idea.
        pka(:,1) = 1._wp - pka(:,2) ! probability of stay put.
        !call kron(pka,py,pt18) ! redundant
        
        !write(4000+trial_id,fmt='("#4",13x,8(e21.14,x))') kv, probtk
        !write(4000+trial_id,fmt='("#5",13x,8(e21.14,x))') pka(:,2), pka(:,1)
        
        ! #########################################################################################
        ! # business technology shock levls (z2) and transition probability (phi), Meh, p.700 072316
        phi(0) = phi0 ! useless
        phi(1) = phi1 ! 0.75 (good shock: receiving advanced business project)
        phi(2) = phi2 ! 0.92 (good shock: receiving advanced business project)
        phi(3) = phi3 ! 0.97 (good shock: receiving advanced business project)
        
        do i = 0, kdim-1 ! 3.12.2017 checked. ! 4.14.2017 correct the start point to 0.
            pz2(i,1) = 1._wp - phi(i) ! bad luck
            pz2(i,2) = phi(i) ! good luck
        enddo
        
        ! level of good business shocks
        do i = 0, kdim-1
            z2(i) = zbar/phi(i) ! to be calibrated # ! <-----------------------------5   
        enddo
        
        !! #########################################################################################
        !! ### annual basis
        !delmeh(1) = 0.049 ! Meh (2005) depreication rate for good technological shock on annual basis
        !delmeh(2) = 0.059 ! Meh (2005) depreication rate for good technological shock on annual basis
        !delmeh(3) = 0.061 ! Meh (2005) depreication rate for good technological shock on annual basis
        !
        !! annual depreciation rate for different investment scales 
        !delzh = 0._wp
        !delzl = 0._wp
        !do i = 1, kdim-1
        !    delzh(i) = (delmeh(i)-PropHouseCapital*deltah)/(1._wp-PropHouseCapital)
        !    !delzl(i) = (deltak-phi(i)*delzh(i))/(1._wp-phi(i))
        !enddo
        !delzl(1:kdim-1) = 2._wp*delzh(1) ! As in Meh (2005)
        
        !!see capital_depreciation_rate.xlsx for detailed calculation and explaination
        delzh = 0._wp
        delzl = 0._wp
        
        delzl    = 0.176_wp ! Set 
        delzh(1) = 0.087_wp ! Set based on the result from Excel optimization solver Feb 13, 2017. 3.12.2017 checked.
        delzh(2) = 0.103_wp ! Set 
        delzh(3) = 0.107_wp ! Set 
        
        ! ### 5 year basis Feb 13, 2017
        delzh  = delzh*length ! converted to be on five year basis
        delzl  = delzl*length ! converted to be on five year basis. One of the determinants of reminant value of investement capital     
        deltak = deltak*length 
        deltah = deltah*length       
        
        !write(4000+trial_id,fmt='("#6",13x,12(e21.14,x))') phi, pz2(:,1), pz2(:,2)
        !write(4000+trial_id,fmt='("#7",13x,10(e21.14,x))') delzh, delzl, deltak, deltah
        
        !! ## transfer depreciation rates to 5 year basis 072316
        !do i = 1, 3
        !    delzl(i) = 1._wp - (1._wp-delzl(i))**5._wp  
        !    delzh(i) = 1._wp - (1._wp-delzh(i))**5._wp 
        !enddo        
        
        !! 7-9-2017
        !staxe = 0.42_wp
        !ptaxe = 1.4_wp
        !staxw = 0.2154_wp
        !ptaxw = 0.7646_wp
        
        ! #########################################################################################
        ! adjustment for the unit of measurement. See, ex., Nakajima (2010). 072316
        ! Done the check again on the tax parameters used by AER2009 Cagetti and De Nardi. 01-29-2017
        ! ---[pg. 107] or their Fortran code, InitialSS.f90, in lines from 351 to 356.
        ! *** See the appendix of housing_v9.pdf for the details of the outcome after scaling and variable change
        ! *** 2-4-2017
        staxebase = staxe*(45._wp/25._wp)**ptaxe ! 090616  ! Cagetti's version      
        staxwbase = staxw*(45._wp/25._wp)**ptaxw ! 090616  ! Cagetti's version
        staxbase  = stax ! 01-29-2017 this is used for Nakajima's version. Obsolete.
        
        ! #########################################################################################
        ! Used for initial guess
        AggEffLab = 0._wp
        do i = 1, 9 
            AggEffLab = AggEffLab + popfrac(i)*efflab(i)
        enddo
        AggEffLab = AggEffLab*dot_product( yv(1:nmc), sy(1:nmc) ) ! initial guess # 1   
        ! call ss(efflab,'2017efflab')
        !call ss(yv,'2017yv')
        !call ss(sy,'2017sy')
        
        
        ! 4.21.2017 comment out and moved inside the interest rate loop.
        !call set_grids_for_housing_asset(hv,9) ! ok  
        call grid(hv,hmin,hmax,2.5_wp) ! 09182016 
        call grid(rhv,hmin,hmax,2.5_wp) ! keep it, not redundant or revision needed. 10102016
        rhv = hv ! 10102016
        
        ! insert "zero" financial asset holdings
        call set_grids_for_nonhousing_asset(av,3.5_wp,'coarse') ! ok bug   
        int1 = locate(av,0._wp) 
        if(av(int1)/=0._wp)then ! 1-29-2017 Be sure to set variable "amin' to be non-positive (either zero or negative number is acceptable).
            rav = av 
            av(int1+1)=0._wp                            
            av(int1+2:adim) = rav(int1+1:adim-1)
            rav = av
        endif
        ! 3.17.2017 Actually, after this block, rav is identical to "av."
        
        rbar = 0.035_wp ! rbar hard-wired for the first run in r loop ! 9-14-2017 ! the reason why program can not replicate the result. 9-16-2017
        taubalmin = 0.0345_wp ! 9-15-2017                                         ! the reason why program can not replicate the result. 9-16-2017
        taubalmax = 0.0363_wp ! 9-15-2017                                         ! the reason why program can not replicate the result. 9-16-2017
        taubal = 0.033_wp ! 9-16-2017                                             ! the reason why program can not replicate the result. 9-16-2017
        transbeq = 0._wp     ! 9-30-2017
        !transbeq_new = 0._wp ! 9-30-2017
        
        !call evenly_distributed_since(av,10)
        !call set_grids_for_nonhousing_asset(rav,3.5_wp,'refined')
        
        ! #########################################################################################
        ! initialization of the outer loop
        iterar = 1
        epsir  = 1.0_wp
        bracketr = 1
        iteratot = 1
        iteragov8rate = 0
        
        !taubal = taubalmin ! initialization ! redundant
        ! inner loop        
        govbal2gdpmin=0.0_wp
        govbal2gdpmax=0.0_wp        
        ! outer loop --> inner loop
        if(mode6taskid==0)then
            taubalrbarmax=taubalmax ! By default, the value is 0.08
            taubalrbarmin=taubalmin ! By default, the value is 0.07
        elseif(mode6taskid==1)then
            taubalrbarmax=tauwealth+0.01_wp
            taubalrbarmin=tauwealth-0.01_wp            
        endif
        
        !write(4000+trial_id,fmt='("#8",13x,10(e21.14,x))') staxebase, staxwbase, staxbase, AggEffLab, hv(1), hv(hdim), av(1), av(adim), taubalrbarmax, taubalrbarmin
        call make_next_period_index_combination() ! 3.15.2017 checked. keep it. moved here from the inner loop. KEEP IT!! THE OPERATION CONTAINS COMBINATION INFO.
        ! 4.1.2017 revision made to introduce "end-of-period" outcomes for people are young, would-be entrepreneurs and still have the opportunity to conceive new entrepreneurial idea  
        
        if(printout1) call print_basic_model_vectors() 
        
        call system_clock(tend) 
        !write(*,'(a,f12.4,a)') '3: ',real(tend-tstart,wp)/real(trate,wp), ' seconds'         
        
        do while((epsir>epsirmin).and.(iterar<=iterarmax)) ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ [1] 
            ! The whole algorithm for stationary equilibrium has been checked. Feb 4, 2017. Cagetti's mistake has been corrected by me [has not confirmed by them].
            ! Tax scheme has been checked. Feb 4, 2017
            
            !call set_grids_for_housing_asset(hv,9) ! ok  
            
            !! 4.21.2017 moved from the block outside the interest rate loop.
            !if(iterar>1)then
            !    hmin = 0.5_wp*medwokinc
            !    hmax = 15._wp*medwokinc ! 4.21.2017 The number 12. is eye-balled based on the SCF graph, figure 5 partial estimation of age profile based on SCF data in my dissertation.
            !endif 
            !call grid(hv,hmin,hmax,2.5_wp) ! 09182016 
            !call grid(rhv,hmin,hmax,2.5_wp) ! keep it, not redundant or revision needed. 10102016
            !rhv = hv ! 10102016            
            
            ! to initialize the inner loop
            iteragov = 1
            epsigov  = 1.0_wp
            bracketgov = 1
            noneedtaubalmax = 0  
            
            if(mode6taskid==0)then
                benchrbar = rbar ! passed down the equilibrium interest rate in the last iteration of solving "benchmark" model.
            else
                if(iterar==1) rbar = benchrbar ! benchmark equilibirum interest rate as the initial price of capital.
            endif !mode6taskid
            
            ! To make the conversion to five year basis for "prices"
            rd = rbar*length ! updated. See equilibrium 534: rbarimplied are rbar are one year interest rates. "rd" and "rimplied" are five year interest rates.
            rl = (rbar+gamma1)*length
            gamma5  = rl - rd        
            
            ! "five-year" rental rate of labor efficiency. formular is checked 1-30-2017
            wage = (1._wp - alpha)*((rd+deltak)/(alpha))**(alpha/(alpha - 1._wp)) ! # 2' 08272016 See proof in housing_v9.lyx or pdf.        
            
            !write(4000+trial_id,fmt='("#9",13x,2(e21.14,x),3(18x,i3,x),(e21.14,x),2(18x,i3,x),4(e21.14,x))') epsir, epsirmin, iterar, iterarmax, iteragov, epsigov, bracketgov, noneedtaubalmax, rd, rl, gamma5, wage
            
            !benefit = merge( tauss*wage*AggEffLab/sum(popfrac(10:14)), tauss*wage*poppaysstaximplied/sum(popfrac(10:14)), iterar==1 ) ! 10122016.

            call coarse_SBE_profit_matrices(c_grs_mat,c_lab_vec,c_opt_vec) ! 09282016 3.4.2017   
            
            !write(4000+trial_id,fmt='("#10-biz",8x,3(e21.14,x))') sum(c_grs_mat), sum(c_lab_vec), sum(c_opt_vec)
            
            do while((epsigov>epsigovmin).and.(iteragov<=iteragovmax)) ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ [2]
                
                iteragov8rate = iteragov8rate + 1 !10.18.2017
                
		        if(bracketgov==1)then
		        	taubal=taubalmin ! taubal is the variable to be tried.
                    !if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 3 '
		        elseif (bracketgov==2)then
		        	taubal=taubalmax
                    !if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 4 '
		        else ! [Domain]taubal->[Image]govbal2gdp
		        	taubal=taubalmin-(taubalmax-taubalmin)*govbal2gdpmin/(govbal2gdpmax-govbal2gdpmin)
                    !if(printout5) write(unit=104,fmt='(3i3,a,2(f8.4))') iterar, iteragov, iteratot, ' 5 ', taubal, gdp
                endif ! bracketgov  
                
                
                write(unit=my_id+1001,fmt='(a,i4,x,a,i4,x,a)') 'iterar ', iterar, ', iteragov ', iteragov, '----------------------------------------------------------------------- '                
                if(mode6taskid>0.and.iterar==1.and.iteragov==1)then
                    write(my_id+1001,'((12x,"old-hmin",2x),(12x,"old-hmax",2x),(11x,"taubalmin",2x),(11x,"taubalmax",2x),(7x,"govbal2gdpmin",2x),(7x,"govbal2gdpmax",2x),(11x,"tauwealth",2x))') 
                    write(my_id+1001,'((e20.5,2x),(e20.5,2x),(e20.13,2x),(e20.13,2x),(e20.13,2x),(e20.13,2x),(e20.13,2x))') hmin, hmax, taubalmin, taubalmax, govbal2gdpmin, govbal2gdpmax, tauwealth
                else
                    write(my_id+1001,'((12x,"old-hmin",2x),(12x,"old-hmax",2x),(11x,"taubalmin",2x),(11x,"taubalmax",2x),(7x,"govbal2gdpmin",2x),(7x,"govbal2gdpmax",2x),(14x,"taubal",2x))') 
                    write(my_id+1001,'((e20.5,2x),(e20.5,2x),(e20.13,2x),(e20.13,2x),(e20.13,2x),(e20.13,2x),(e20.13,2x))') hmin, hmax, taubalmin, taubalmax, govbal2gdpmin, govbal2gdpmax, taubal
                endif
                
                !write(4000+trial_id,fmt='("#11",12x,2(e21.14,x),3(18x,i3,x),(e21.14,x),2(18x,i3,x),4(e21.14,x))') epsigov, epsigovmin, iteragov, iteragovmax, bracketgov, taubal, iterar, iteragov, taubalmin, taubalmax, govbal2gdpmin, govbal2gdpmax
                
                ! VARIABLES THAT NEED TO RENEW
                if(iterar==1.and.iteragov==1)then
                    
                    ! 1-28-2017
                    ! procedure to obtain the initial guess on GDP for getting the estimate of lump sum transfer and average income of working class households.
                    
                    if(mode6taskid==0)then
                        kndata   = ((rd+deltak)/alpha)**(1._wp/(alpha-1._wp))  
                        crplab   = AggEffLab*CorpLabFrac ! Bold assumption # 1         
                        crpcap   = kndata*crplab ! one period. Determined by the guessed AggCorpLab
                        crpprd   = crpcap**alpha*crplab**(1._wp-alpha) ! one period         
                        gdp      = crpprd/CorpOutFrac ! one period, aggregate output. Needs to be updated. # 3 
                        
                        transbeq = b2gratio*gdp ! 9-30-2017 ! <------------------------------------------------------------ should be revised 10132016.       
                        !transbeq_new = b2gratio*gdp ! 9-30-2017
                        
                        avgincw  = wage*dot_product(yv,sy)/5._wp ! wage*AggEffLab/0.73_wp ! *0.86_wp ! initial guess. assume #retiree/(#worker+#retiree+#entrepreneur) = 0.15 and #entrepreneur/(#worker+#retiree+#entrepreneur)=0.12, so population share of worker is 0.73, work force share of worker is 0.73/(0.73+0.12)=0.86
                        benefit  = tauss*wage*AggEffLab/sum(popfrac(10:14))
                    elseif(mode6taskid==1)then
                        
                        transbeq = transbeqimplied !transbeq
                        avgincw  = mean_wokinc/5._wp !avgincw
                        benefit  = benefitimplied !benefit
                        
                    endif
                    !if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 1 '
                    
                    !write(4000+trial_id,fmt='("#12-1",10x,8(e21.14,x))') kndata, crplab, crpcap, crpprd, gdp, transbeq, avgincw, benefit
                    
                elseif((iterar/=1.and.iteragov==1))then ! iteragov 5.10.2017
                    
                    ! experiment for speed up
                    ! transbeq = transbeqimplied ! 4.17.2017 `transbeqimplied` is updated in subroutine `lump_sum_transfer`.
                    transbeq = transbeqimplied ! 9-30-2017
                    !transbeq_new = transbeqimplied ! 9-30-2017
                    
                    avgincw  = mean_wokinc/5._wp ! 10102016 annual income. 4.17.2017 `mean_wokinc` is computed in subroutine `macro_statistic`. Note: `avgincw` is on annual basis.
                    ! `benefit` is updated also in subroutine `macro_statistic`. 4.17.2017 
                    !benefit  = tauss*wage*poppaysstaximplied/sum(popfrac(10:14)) ! ---- needs to be revised Feb 5, 2017 ! 4.17.2017 comment out. The exact update takes place in line 3681 of model.f90.
                    benefit  = benefitimplied ! 7-6-2017 Good. Faster in convergence. (27 mins vs 32 mins)
                    
                    !! Secant method doesn't help. It took more time to complete an expected convergence (37 mins vs the best history record 27 mins) 7-7-2017  
                    !transbeq = 0.5_wp*(transbeq + transbeqimplied)
                    !avgincw = 0.5_wp*(avgincw + mean_wokinc/5._wp)
                    !benefit = 0.5_wp*(benefit + benefitimplied) 
                    
                    ! 4.21.2017 moved from the block outside the interest rate loop.
                    !hmin = 0.5_wp*lowest_quintile_wokinc ! 10-8-2017
                    if(printout18)then
                        call grid_housing_upper_bound_adjustment(hmax, h_mass_vec) ! 10-9-2017 #1    
                    else ! tranditional method
                        hmax = 15._wp*medwokinc ! 4.21.2017 The number 12. is eye-balled based on the SCF graph, figure 5 partial estimation of age profile based on SCF data in my dissertation. Yang seems to use 17 as the multiplier.
                    endif ! printout18
                    hmin = iota*lowest_quintile_wokinc
                        
                    call grid(hv,hmin,hmax,2.5_wp) ! 09182016 
                    call grid(rhv,hmin,hmax,2.5_wp) ! keep it, not redundant or revision needed. 10102016
                    rhv = hv ! 10102016            
                    
                    !write(4000+trial_id,fmt='("#12-2",10x,5(e21.14,x))') transbeq, avgincw, benefit, hmin, hmax
                    
                    ! NOTE: Nakajima (2010) leaves the Social Security tax rate to be determined such that when the government is balancing the budget 
                    ! in each period, the model replicates the replacement ratio.
                    ! if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 2 '
                    
                    !!! **************  THE FOLLOWIGN BLOCK HAS BEEN MOVED TO LINE 1131: OUTSIDE THE GOV LOOP AND IN THE BOTTOM OF THE R LOOP AFTER NUMERICAL UPDAE OF R. ! 9-14-2017
                    !!! 8-1-2017
                    !new_amin = amin
                    !new_amax = amax
                    !call grid_boundary_inspection(new_amin,new_amax,mass_vec)
                    !do i = 1, adim-1
                    !    write(my_id+1001, '(4x,a,i3)', advance='no') "level",i
                    !enddo
                    !write(my_id+1001, '(4x,a,i3)') "level",adim
                    !write(my_id+1001, '(<adim>f12.9)') av
                    !write(my_id+1001, '(<adim>f12.9)') mass_vec
                    !write(my_id+1001, '(6x,a,x,6x,a)') "  amin", "  amax"
                    !write(my_id+1001, '(f12.7,x,f12.7)') amin, amax 
                    !write(my_id+1001, '(6x,a,x,6x,a)') " namin", " namax"
                    !write(my_id+1001, '(f12.7,x,f12.7,/)') new_amin, new_amax
                    !amin = new_amin
                    !amax = new_amax
                    !call set_grids_for_nonhousing_asset(av,3.5_wp,'coarse')
                    !rav = av ! Bug. 8-1-2017

                endif ! iteragov
                
                write(my_id+1001,'((12x,"new-hmin",2x),(12x,"new-hmax",2x))') 
                write(my_id+1001,'((e20.5,2x),(e20.5,2x))') hmin, hmax
                
                !write(*,'(2(a,i4,x),x,3(a,f8.5,x))') 'October-9-2017 iterar: ', iterar, "iteragov: ", iteragov, "hmin: ", hmin, "hmax: ", hmax, "new_hmax", new_hmax
                
                if(iterar==1.and.iteragov==1)then

                    !! 8-25-2017 version 1
                    !write(unit=my_id+1001,fmt='(a,i7.7,x,a,<ndim>f16.7)') 'No.',trial_id,'Inputs: ',guessv ! Note: trial_id falls in [1,trylen], not used directly for the list of "indexseries"
                    !write(unit=my_id+1001,fmt='(3x,a)') 'rd(5y)'
                    !write(unit=my_id+1001,fmt='(f9.4)') rd
                    !! 8-25-2017 version 2                    
                    write(my_id+1001,'("Trial no.",x,i7," ... Fresh Start ...")') trial_id 
                    do i = 1, ndim-1
                        write(my_id+1001,'((14x,"dim",i3,2x))',advance='no') i
                    enddo ! i
                    write(unit=my_id+1001,fmt='((15x,"dim",i3,2x))') ndim
                    write(my_id+1001,'(<ndim>(e20.13,2x))') guessv
                    
                    !write(unit=6,fmt='(3x,a)') 'rd(5y)'
                    !write(unit=6,fmt='(f9.4)') rd                    

                    !write(unit=my_id+1001,fmt='(2f12.6)') rbar, rbarimplied, 
                    
                    !!! 8-25-2017 version 1
                    !write(unit=my_id+1001,fmt='(2x,a,6(x,a))') 'TRANSFER(S): ', 'transbeq   ', 'avgincw    ', 'benefit    ', 'taubal     ', 'hmin       ', 'hmax       '
                    !write(unit=my_id+1001,fmt='(15x,6(f12.6))') transbeq, avgincw, benefit, taubal, hmin, hmax

                    !! 8-25-2017 version 2
                    write(unit=my_id+1001,fmt='((15x,"rd-ty",2x),(12x,"transfer",2x),(13x,"avgincw",2x),(13x,"benefit",2x),(14x,"taubal",2x),(16x,"hmin",2x),(16x,"hmax",2x))') 
                    write(unit=my_id+1001,fmt='(7(e20.13,2x))') rd, transbeq, avgincw, benefit, taubal, hmin, hmax                    
                    
                    !! 8-25-2017 comment out
                    !write(unit=my_id+1001,fmt='(2x,a,3(x,a))') 'CONVERGENCE: ', 'fundiffnow  ', 'taubalmax   ', 'taubalmin   '
                    !write(unit=my_id+1001,fmt='(15x,3(f12.6),/)') fundiffnow, taubalmax, taubalmin
                endif
                
                staxw = staxwbase*avgincw**(-ptaxw) ! that is the real denominator in the (inc/45000)**p, where 45000 is a rough estimate just for initialization
                staxe = staxebase*avgincw**(-ptaxe) ! that is the real denominator in the (inc/45000)**p, where 45000 is a rough estimate just for initialization
                stax  = staxbase*(avgincw/50000._wp)**(-ptax) ! Obsolete. Naka uses average household income in 1999. unstable in my current model ! 1-31-2017 seems redundant
                
                write(my_id+1001,'((15x,"staxw",2x),(15x,"staxe",2x))') ! isye 
                write(my_id+1001,'(2(e20.13,2x),/)') staxw, staxe ! isye
                
                ! boundary 072316
                !hmin = hmin*gdp !   
                !hmax = hmax*gdp !  
                !amax = amax*gdp !     
                      
                !call get_weight_new_grids_for_refined_grid(av,afv,wint,afint)
                
                !call set_grids_for_nonhousing_asset(afv,3.5_wp,'refined') ! ok  seems useless 10122016.   
                
                ! 3.15.2017 Being moved outside the inner loop.
                !call make_next_period_index_combination() ! 09082016 for invariant function (coded in variable.f90, ln.1364) Keep it 10012016
                
                !call bizret_mat(bizmat) ! for generating business return matrix in order to save time      ! seems redundant 090816
                
                !if(printout5) write(unit=104,fmt='(3i3,a,f8.4)') iterar, iteragov, iteratot, ' 5-1 ', taubal
                
                ! Initialization 072316 -----------------------------------------------------------------------------
                !! used in coarse policy
                tww = -99 ! tww: indicator to switch to be a labor in the beginning of the period. 1 as switches (nine states)
                twf = penalty ! the end of period distribution (nine states) in the "solving policy function" stage
                wf  = penalty ! the beginning period distribution (six states) in the "solving policy function" stage
                twa = penalty
                twh = -99 ! discrete
                twk = -99 ! discrete
                
                ! 8-24-2017 Used in convert_2d_outcome_into_series subroutine.
                sww = -99 ! indicates the begining distribution is invalid if sww=1 ! 8-24-2017
                swf = penalty
                swa = penalty
                swh = penalty
                swk = -99
                swc = penalty
                
                !! used in refined policy (xxxm); used in household decision (xxx_m)
                uww   = -99
                uwa   = penalty
                uwh   = -99
                uwk   = -99
                
                tbim  = penalty ! 0._wp
                atwm  = penalty ! 0._wp
                taxm  = penalty ! 0._wp
                beqm  = penalty ! 0._wp
                ldemm = penalty ! 0._wp
                lsupm = penalty ! 0._wp
                prom  = penalty ! 0._wp
                cspm  = penalty ! 0._wp
                wok2ent = 0._wp        
                ent2wok = 0._wp
                !nafv = penalty
                nsav = penalty

                totast_m = 0._wp
                entcap_m = 0._wp
                entlab_m = 0._wp
                entprd_m = 0._wp
                uwa_m    = 0._wp
                uwh_m    = 0._wp
                home_m   = 0._wp
                ttax_m   = 0._wp
                
                cwf  = penalty
                cwa  = penalty
                cwh  = penalty
                cww  = 1 ! default 1 means bad point. -99 indicates a good point (or normal situation).
                cwk  = -99 ! don't change it. toolbox has subroutine with hardwired constant the same as it for conditional statement.
                cwc  = penalty
                
                cef  = 0._wp
                cvv  = 0._wp ! default 1 means bad point. -99 indicates a good point. model.f90, ln.3316.
                dcef = 0._wp
                
                !write(4000+trial_id,fmt='("#13",12x,5(e21.14,x))') staxw, staxe, stax
                
                !if(num_procs==2) write(*,*) ' '
                !!!!!!! ---- kernel ----- Option I 4.8.2017 fix the wrong indexing of parallel liss (0)
                call solve_bellman_1014(exit_bellman) ! 072516 one more time 09122016 3.15-17.2017
                write(unit=my_id+1001,fmt='(a)',advance='no') ' Bellman 10-14'
                !write(*,*) merge('exit 1014','stay 1014',exit_bellman) ! 8-13-2017
                if(exit_bellman) exit_log1 = .true. ! 8-13-2017
                if(exit_bellman) exit ! exit the government loop 8-13-2017
                
                call solve_bellman_9(exit_bellman)    ! 072516 one more time 09122016 3.15-17.2017
                write(unit=my_id+1001,fmt='(a)',advance='no') ', bellman 9'
                if(exit_bellman) exit_log1 = .true. ! 8-13-2017
                if(exit_bellman) exit ! exit the government loop 8-13-2017
                
                call solve_bellman_18(exit_bellman)   ! 072516 one more time 09122016 3.15-17.2017
                write(unit=my_id+1001,fmt='(a)',advance='no') ', bellman 1-8'
                if(exit_bellman) exit_log1 = .true. ! 8-13-2017
                if(exit_bellman) exit ! exit the government loop 8-13-2017
                
                ! 4.8.2017 4:34 pm stop here.
                !! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                !! 3.17.2017 Instruction: run (1) at least once, producing solution of Bellman eqn; run (0,1,4,5) if want to check if the save-load-again-save process works; 
                !! 3.17.2017 Once the process works and the Bellman solution is obtained. Just run (4); In production model, turn off any of these subroutines (1-5). Only let (0)(1) work.
                
                !!! bookkeeping. 3.15.2017 (model.f90 and line 4212)
                call convert_2d_outcome_into_series('coarse') ! (1) `SAVING` THE OUTCOME FROM SOLVING THE BELLMAN EQUATIONS !!! Need to turn "on" in PRODUCTION mode !<-------
                
                !!!!!!! Useless block unless it is for testing the success of file IO.
                !!!!!! ---- run in development stage with the kernel block above for inspection file I/o success.
                !call convert_1d_file_into_matrix('coarse') ! (2) READ THOSE SAVED BELLMAN EQUATION SOLUTIONS
                !call convert_2d_outcome_into_series('test') ! (3) TEST WHETHER THIS FILE I/O WORKS NORMALLY.
                !!!! ---- Option II. run in development stage for time saving. only used after kernel block is executed without comment once.
                !call convert_1d_file_into_matrix('distribution-stage') ! (4) READ THOSE SAVED BELLMAN EQUATION SOLUTIONS !<-- only this one would be sufficient. comment out for development
                !call convert_2d_outcome_into_series('distribution-stage-test') ! (5)               
                
                err_dist = 1.e-12_wp
                errdist  = 1000._wp        
                inv_dist_counter = 1
                
                ! 8-24-2014 Generate report on whether the combinations on the beginning of each period are valid or not.
                call valid_beginning_period_mass(iterar,iteragov,iteratot) !! 4.2.2017 IMPORTANT! GENERATE "CVV," "(model.f90, ln.3186)", if a tested combination is valid, CVV = -99.
                if(printout7) call print_coarse_2d_brent_mat(iterar,iteragov,iteratot) ! 4.2.2017 print out cwa, cwh, cwk, cw... outcome matrices including cvv.
                ! 3.22.2017 s1c, c1s are both created here.
                ! 3.16.2017 The subroutine is used for indicating whether ALL the related end-of-period combinations of a beginning-of-period combination are valid.
                ! 3.16.2017 The subroutine has to work with subroutine make_next_period_index_combination defined above.

                call system_clock(tstart,trate,tmax)
                
                do while(errdist >= err_dist .and. inv_dist_counter<=iterainvdist) ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ [3]    
                    !dcef = cef ! 3.16.2017 ! move it inside the mass_transition.
                    !write(*,fmt='("iterainvdist:", i4)') iterainvdist ! 8-26-2017
                    !! The primary of this mass distribution convergence loop. 8-24-2017
                    call mass_transition(exit_log1, errdist, iterar, iteragov, iteratot, trial_id) ! 3.16.2017 "inv_dist_counter" is initialized in modelf.f90, ln.4725 ! 3.20.2017 Not yet only inv_dist_counter == 1. <----

                    if(inv_dist_counter==1)then
                        allocate(sef1(szperiod1),sef2(szperiod1),sef3(szperiod1) ) ! szperiod1 is set up in subroutine mass_transition. 
                        sef1 = 0._wp ! 8-24-2017
                        sef2 = 0._wp ! 8-24-2017
                        sef3 = 0._wp ! 8-24-2017
                    endif 
                    
                    if(exit_log1==.true.)then
                        if(printout3) write(unit=my_id+1001,fmt='(/,a,/)') ' Distribution error. Exit. '
                        exit ! 20132016. Should I comment out this statement?? <----------------------------------------------------------------------------
                    endif
                    
                    !write(*,fmt='(a,i4,a,i4)') 'round: ', inv_dist_counter, ', flag: ', exit_log1
                    
                    call convert_2d_distribution_into_1d() ! 3.25.2017 Don't comment out!! because Subroutine intergenerational_transfer uses it. ! The stationary distribution is stored in the matrix "sef" rather than "sef1" (defined in model.f90 and line 4679)
                    call intergenerational_transfer() ! 4.14.2017 move errdist to mass_transition.
                    
                    !call intergenerational_transfer(errdist) ! 3.23.2017 3.24.2017 done ! 4.14.2017 remove errdist, replace it in subroutine mass_transition.
                    ! call ability_transition(errdist) ! 10102016
                    !call lump_sum_transfer() ! 3.25.2017 Be moved outside the loop. New lump sum transfer. should be used only for the outer R, gov loop, I think (because it affects policy function, not the balance of government budget). Should add errdist to measure distance. 10102016
                    
                    !if(printout6.and.inv_dist_counter>=2.and.(num_procs==2.or.num_procs==1)) write(*,fmt='(a)',advance='no') ' # ' 
                    !if(printout6.and.inv_dist_counter>=2.and.(num_procs==2.or.num_procs==1)) write(*,fmt='(i4,a,f15.12,a,f15.12,x,a,i5.5,x,a,i3.3,x,a,l2)') inv_dist_counter, ' dist error ', errdist, ' sum period 1 ', sum(sef1), 'trial_id', trial_id,'iteratot',iteratot, 'exit_log1:', exit_log1
                    write(unit=my_id+1001,fmt='(a)',advance='no') ' # ' 
                    !write(unit=my_id+1001,fmt='(i4,a,f15.12,a,f15.12,x,a,i5.5,x,a,i3.3,x,a,l2)') inv_dist_counter, ' dist error ', errdist, ' sum period 1 ', sum(sef1), 'trial_id', trial_id,'iteratot',iteratot, 'exit_log1:', exit_log1 ! 7-10-2017
                    inv_dist_counter = inv_dist_counter + 1
                    iteratot = iteratot + 1
                enddo
                
                call system_clock(tend)
                !write(*,fmt='(a,f12.4,a)') ' mass dist convergence time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'                  
                write(unit=my_id+1001,fmt='(a)') ' - '
                
                if(exit_log1==.false.)then
                    call lump_sum_transfer() ! #1 3.25.2017 moved from the distribution do loop. 4.17.2017 This subroutine includes the update of `transbeq`.
                    !write(unit=6,fmt='(a,e20.13)') 'transbeqimplied: ', transbeqimplied ! 9-11-2017 Done!
                    
                    if(printout9) call print_2d_end_of_period_dist_mat(iterar,iteragov,iteratot)
                    
                    !call convert_2d_distribution_into_1d() ! #2 ! DON'T COMMENTED OUT. USED IN MACRO_STATISTICS SUBROUTINE. 9-11-2017 Done! Seems redundant. Should be removed.
                    call macro_statistics( momvec, iterar, iteragov, iteratot, exit_log1, msg,trial_id) ! #3 4.21.2017 I don't know why the medwokinc (involving allocatable array) causes error but subroutine compute_lorenz below doesn't. Maybe in the future I can use subroutine that includes the alloctable array to replace the use of alloctable array in the main program.
                    
                    !write(4000+trial_id,fmt='("#14",12x,10(e21.14,x),3(16x,i5,x))') momvec, iterar, iteragov, iteratot
                    
                    !call grid_boundary_inspection() ! moved outside 9-11-2017
                    !call convert_2d_macro_statistics_into_1d() ! SHOULD BE COMMENTED OUT IN PRODUCTION MODE
                    !if(printout7) call print_coarse_2d_brent_mat(iterar,iteragov,iteratot) ! 4.2.2017    
                    if(printout4) call printout_for_stata_lifecycle_plotting() ! SHOULD BE TURNED OFF WHEN IN PRODUCTION MODE E:\GoogleDrive\Stata\Research\dissertation\Model_lifecycle_plot.do
                endif
                !call macro_stat_and_reinitialization(iterar,iteragov,iteratot,exit_log1,msg)   
                
                deallocate( sef1, sef2, sef3 ) ! 7-5-2017. Used in 'intergenerational_transfer' defined in model.f90. 
                
                ! ##### The major code block that deals with the associated boundary parameter setting with the convergence of govbal ########
                if(exit_log1==.false.)then
                    !write(*,fmt='(3(a,i4))') ' iterar ', iterar, ' iteragov ', iteragov, ' iteratot ', iteratot
                    
                    ! ###### Note: govbal2gdp is the outcome of solving the model, see model.f90 around line 5429. ######
                    epsigov = min(abs(govbal2gdp), abs(taubalmax - taubalmin)*10.0_wp)
                    
                    !write(4000+trial_id,fmt='("#15",12x,4(e21.14,x))') epsigov, govbal2gdp, taubalmax, taubalmin
                    
                    !if(printout5) write(unit=104,fmt='(3i3,a,f8.4)') iterar, iteragov, iteratot, ' epsir=abs(fundiffnow) ', epsir
                    if(printout5) write(unit=104,fmt='("## ",3i3,6x,a,x,f8.4,/)') iterar, iteragov, iteratot, 'epsigov1:', epsigov
			        if(epsigov>epsigovmin)then
                        !write(4000+trial_id,fmt='("#16",12x,2(e21.14,x))') epsigov, epsigovmin
			        	!using bisection algorithm to update taubal
    		        	if(bracketgov==1)then
			        		if(govbal2gdp>0.0)then
           	        			noneedtaubalmax=1 ! switch: note that we obtain taubalmax, no need to look for max next round. initial value = 0.
			        			! in this case
           	        			govbal2gdpmax=govbal2gdp
           	        			taubalmax=taubalmin	        ! Shift the interview downward.
           	        			taubalmin=taubalmin-pertgov ! [Domain]
			        			! note bracketgov still ==1 ! [Image]
                                !if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 6 '
                                if(printout5) write(unit=104,fmt='(4(a,x,f8.4,x),a)') 'taubal', taubal, ' taubalmin = taubalmax - ', pertgov, ' = ', taubalmin, '  taubalmax:', taubalmax, ' bracketgov==1, obtained bal>0, save bal2gdp as bal2gdpmax, interval needs to move downward '
                                !write(4000+trial_id,fmt='("#17",12x,5(e21.14,x),(17x,i4,x))') bracketgov, govbal2gdp, govbal2gdpmax, taubalmax, taubalmin, noneedtaubalmax
			        		else
			        			govbal2gdpmin=govbal2gdp
       		        			bracketgov=2
                                !if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 7 '
                                if(printout5) write(unit=104,fmt='(2(a,x,f8.4,x),a)') 'taubal', taubal, 'govbal2gdpmin:',govbal2gdpmin,' bracketgov == 1, used taubalmin obtained govbal<0, so save bal2gdp as bal2gdpmin '
                                !write(4000+trial_id,fmt='("#18",12x,(17x,i4,x),4(e21.14,x),(17x,i4,x))') bracketgov, govbal2gdp, govbal2gdpmin, taubalmax, taubalmin, noneedtaubalmax
			        		endif
			        		
        	        		if((govbal2gdp<=0.0).and.(noneedtaubalmax==1))then ! this block will be activated when, for example, 1st round we have Max, 2nd round we happends to have Min.
			        			bracketgov=0
                                !if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 8 '
                                if(printout5) write(unit=104,fmt='(10x,a,x,f8.4,x,a)') 'taubal', taubal, ' already got bal2gdp>0 (saved as bal2gdpmax) last round, this round got bal2gdp<0 (with the input taubalmin). Start finding root '
                                !write(4000+trial_id,fmt='("#19",12x,(17x,i4,x),4(e21.14,x),(17x,i4,x))') bracketgov, govbal2gdp, govbal2gdpmin, taubalmax, taubalmin, noneedtaubalmax
        	        		endif     
			        	elseif(bracketgov==2)then
        	        		if(govbal2gdp<0.0)then ! This means we've already found MIN, found another MIN this time, and still need to find MAX.
                    			bracketgov=2 ! Commend the program to continue looking for a taubalmax (govbal>0)
                    			taubalmin=taubalmax         ! Shift the interview upward.
                    			taubalmax=taubalmax+pertgov ! [Domain]
			        			govbal2gdpmin=govbal2gdp    ! [Image]
                                !if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 9 '
                                !if(printout5) write(unit=104,fmt='(a,f8.4,a,f8.4,a)') 'taubalmin', taubalmin, 'taubalmax', taubalmax, ' bracketgov == 2 (meant to obtain bal2gdp>0), but still obtain bal2gdp<0. move interval upward '
                                if(printout5) write(unit=104,fmt='(4(a,x,f8.4,x),a)') 'taubal', taubal, 'taubalmin:', taubalmin, 'taubalmax = taubalmax + ', pertgov, ' = ', taubalmax, ' bracketgov==1, obtained bal>0, save bal2gdp as bal2gdpmax, interval needs to move downward section (1) '                                
                                !write(4000+trial_id,fmt='("#20",12x,(17x,i4,x),5(e21.14,x))') bracketgov, govbal2gdp, govbal2gdpmin, taubalmax, taubalmin, pertgov
       		        		else
			        			govbal2gdpmax=govbal2gdp	! It means we found the designated Max in this round, and should ready for tri 
        	        		    bracketgov=0 ! bracketgov=0: computed BOTH taubalmin and taubalmax 
                                !if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 10 '
                                if(printout5) write(unit=104,fmt='(4(a,x,f8.4,x),a)') 'taubal', taubal, 'bal2gdpmax', govbal2gdpmax,' bracketgov == 2 for bal2gdp>0 does catch bal2gdp>0. save bal2gdp as govbal2gdpmax '
                                !write(4000+trial_id,fmt='("#21",12x,(17x,i4,x),5(e21.14,x))') bracketgov, govbal2gdp, govbal2gdpmax, taubalmax, taubalmin, pertgov
	       	        		endif
    		        	else
        	        		! convergence criterion       
        	        		if(govbal2gdp>0.0)then
		            			taubalmax=taubal ! Positive balance means we need to lower the upper boundary of tax rate of state.
                    			govbal2gdpmax=govbal2gdp
                                !if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 11 '
                                if(printout5) write(unit=104,fmt='((x,a,x,f8.4),a,4(x,a,x,f8.4))') 'taubal', taubal, ' obtain bal2gdp>0 in contraction, renew taubalmax and bal2gdpmax', 'taubalmin',taubalmin,'taubalmax',taubalmax,'bal2gdpmin',govbal2gdpmin,'bal2gdpmax',govbal2gdpmax  
                                !write(4000+trial_id,fmt='("#22",12x,(17x,i4,x),5(e21.14,x))') bracketgov, govbal2gdp, govbal2gdpmax, taubalmax, taubalmin, pertgov
        	        		else
                    			taubalmin=taubal ! Negative balance means we need to raise the lower boundary of tax rate of state.
                    			govbal2gdpmin=govbal2gdp
                                !if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 12 '
                                if(printout5) write(unit=104,fmt='((x,a,x,f8.4),a,4(x,a,x,f8.4))') 'taubal', taubal, ' obtain bal2gdp<0 in contraction, renew taubalmin and bal2gdpmin', 'taubalmin',taubalmin,'taubalmax',taubalmax,'bal2gdpmin',govbal2gdpmin,'bal2gdpmax',govbal2gdpmax  
                                !write(4000+trial_id,fmt='("#23",12x,(17x,i4,x),5(e21.14,x))') bracketgov, govbal2gdp, govbal2gdpmin, taubalmax, taubalmin, pertgov
        	        		endif		
    		        	endif
                    endif  
                    
                    ! used for outer loop
                    rbarimplied = rimplied/5._wp ! rimplied is implied from marginal product of capital in line 3676 of model.f90
                    fundiffnow = rbar - rbarimplied ! [image of outer loop] ! note: r with bar indicates the interest rate chosen/updated each round by the convgence algorithm; r without bar is the 5-year interest rate solved for by the model.
                    if(bracketr==1.or.bracketr==21.or.bracketr==22)then       
		            	epsir=abs(fundiffnow)                              
                        if(printout5) write(unit=104,fmt='(3i3,a,f8.4)') iterar, iteragov, iteratot, ' epsir=abs(fundiffnow) ', epsir
		            else                                                ! adjust for the epsilon           
		            	epsir=min(abs(fundiffnow),abs(rbarmax-rbarmin)) ! < less strict convergence constraint.  
                        if(printout5) write(unit=104,fmt='(3i3,a,f8.4)') iterar, iteragov, iteratot, ' epsir=min(abs(fundiffnow),abs(rbarmax-rbarmin)) ', epsir
                    endif	 
		            !write(4000+trial_id,fmt='("#24",12x,4(e21.14,x),(17x,i4,x))') rbarimplied, rbar, fundiffnow, epsir, bracketr 
                    
                    ! keep track ###
                    
                    if(printout3) write(unit=my_id+1001,fmt='(" # of Tot Rnds: ",i4,", # of GovRnds: ",i5,", # of Gov8RateRnds: ",i5,2(a,f12.8),(a,i4))') iteratot, iteragov, iteragov8rate, ', taubal ', taubal, ', epsigov ', epsigov,', iteragov ', iteragov                
                    iteragov = iteragov + 1
                    !iteratot = iteratot + 1 ! 9-18-2017
                else ! exit_log1==.true. ! Macrostat issued an error message.
                    write(unit=my_id+1001,fmt='(" Macrostat exit message:",x,a)') msg 
                    write(unit=my_id+1001,fmt='(" Government Loop Update Fails!",/)')                     
                    exit
                endif ! exit_log1
            enddo ! loop on epsigov and iteragov ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ [2]
            
            if(exit_log1==.false.)then
            
                !if(printout3) write(unit=my_id+1001,fmt='(a)') ' ------------------------------------------------- '
                
                !if(printout3) write(unit=my_id+1001,fmt='((a,i4),(a,f7.4),3(a,f7.4))') '-iterar ', iterar, ' epsir ', epsir, ' input.rbar ', rbar, ' impld.rbar ', rbarimplied, ' diff.rbar ', fundiffnow            
                if(printout3) write(unit=my_id+1001,fmt='(   /,a,(11x,a), (8x,a), (8x,a),     (x,a))') 'Start: # iterar', 'epsir','new.rbar','imp.rbar','diff.rbar(fundiffnow)'
                if(printout3) write(unit=my_id+1001,fmt='(11x,i4,(f16.8),(f16.8),(f16.8),(6x,f16.8))') iterar,epsir,rbar,rbarimplied,fundiffnow
                if(printout3) write(unit=my_id+1001,fmt='(a,i5,a,i5,a,f16.8)') 'iterar: ', iterar, ', iteragov8rate: ', iteragov8rate, ', --- EPSIR: ', epsir
                
		        ! using bisection algorithm to update rbar
	            ! bisection algorithm: in the current case they try to minimize the difference of rbar and rimplied	 

                !if(printout3) write(unit=my_id+1001,fmt='(a,a,i3,3(a,f8.4))') '-old- ', 'bracketr ', bracketr, ' rbarmax ', rbarmax, ' rbarmin ', rbarmin, ' newrbar ', rbar     
                !if(printout3) write(unit=my_id+1001,fmt='(18x,4(a,f8.4))')    ' fundmax ', fundiffmax, ' fundmin ', fundiffmin, ' tbalmax ', taubalmax, ' tbalmin ', taubalmin
                
                if(printout3) write(unit=my_id+1001,fmt='(a,6(10x,a))') 'Previous: bracketr', 'rbarmax', 'rbarmin', 'fundmax', 'fundmin', 'tbalmax', 'tablmin'
                if(printout3) write(unit=my_id+1001,fmt='((13x,i5),6(x,f16.8))')   bracketr, rbarmax, rbarmin, fundiffmax, fundiffmin, taubalmax, taubalmin
                
		        if(epsir>epsirmin)then 
                    if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 15 '
                    !write(*,fmt='(a,f12.8,a,f12.8)'), ' epsir ', epsir, ' epsirmin ', epsirmin
		        	if(bracketr==1)then ! [Domain]->[Image]: rbar->fundiffnow 			
		        		if(fundiffnow>0.0_wp)then ! rbar-rimplied
                            if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 16 '
		        			! <1> mission 1: update rbar
		        			bracketr = 22 ! commend to do searching for fundiffnow<=0 next round
		        			rbarmax = rbar          ! store what have tried  
		        			fundiffmax = fundiffnow ! store what have gotten (rbar-rimplied)
                            
		        			if(abs(rbar-rbarimplied)<=radjust)then ! if not too far
		        				!rbar = rbarimplied ! for next round
                                rbar = max(rbarimplied, 0.000012_wp) ! 5.17.2017
                                if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 17 '
		        			else ! relaxation criterion
		        				rbar = rweight*rbar+(1.0_wp-rweight)*rbarimplied
                                if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 18 '
		        			endif
		        			! <2> mission 2: update the boundary for pinning down taubal
		        			!IF (dogovloop==1) THEN ! In this block, we adjust the boundary of tax rate of state.
		        				taubalmax = taubal + 0.0001_wp	 ! shift the state tax rate downward to have more private capital in the country.
		        				taubalmin = taubal - pertgov ! shift the state tax rate downward to have more private capital in the country.
		        				taubalrbarmax = taubal ! <============== used for interpolating the interval of state tax rate (used in the inner loop) based on the input "rbar" from the outer loop.
		        			!END IF
                            if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 19 '    
                            !write(4000+trial_id,fmt='("#25",12x,(17x,i4,x),7(e21.14,x))') bracketr, fundiffnow, rbarmax, radjust, rbar, taubalmax, taubalmin, taubalrbarmax
                        else ! ======================================================================================================
		        			! mission <1>: update rbar
                            if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 20 '
		        			bracketr = 21 ! commend to do searching for fundiffnow>=0 next round (fundiffnow)
		        			rbarmin = max(rbar,0.00001_wp) ! pointer of the outter loop
		        			fundiffmin = fundiffnow ! result of the outter loop (rbar-rimplied)
		        			if(abs(rbar-rbarimplied)<=radjust)then
		        				rbar = max(rbarimplied, 0.000012_wp) ! In the current result (fundiffnow>=0 in the first round), we know rimplied is relatively low. So this line set a lowe boundary to prevent a crazy negative rimplied contaminting the searching effort.
                                if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 21 '
		        			else ! if crazy, use relaxion criterion. 
		        				rbar = max(rweight*rbar+(1.0_wp-rweight)*rbarimplied, 0.000012_wp)
                                if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 22 '
		        			endif
		        			! mission <2>: adjusting the interval for pinning down taubal
		        			!IF (dogovloop == 1) THEN
		        				taubalmin = taubal - 0.0001_wp ! fix the lower boundary
		        				taubalmax = taubal + pertgov
		        				taubalrbarmin = taubal ! used for the secant method that maps the input of the otter loop (rbar) to the middle point of searching interval of the inner loop
                                if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 23 '
		        			!END IF
                            !write(4000+trial_id,fmt='("#26",12x,(17x,i4,x),7(e21.14,x))') bracketr, fundiffnow, rbarmin, radjust, rbar, taubalmax, taubalmin, taubalrbarmin
		        		endif  
		        	elseif (bracketr == 21) then
		        		! In this round (round 2), we are designated to search for fundiffnow>=0
		        		if(fundiffnow>0.0)then ! fundiffnow=rbar-rimplied
		        			! mission <1>: update rbar
                            if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 24 '
		        			bracketr=0 ! start bisection for interest rate.
		        			rbarmax=rbar          ! update/store and mission completed
		        			fundiffmax=fundiffnow ! update/store and mission completed
		        			rbar=rbarmin-(rbarmax-rbarmin)*fundiffmin/(fundiffmax-fundiffmin) ! https://en.wikipedia.org/wiki/Secant_method
		        			! mission <2>: update the tax rate of state interval
		        			!IF (dogovloop==1) THEN
		        		    taubalrbarmax = taubal ! the taubal that corresponds to rbarmax 2-4-2017 ! update with the input for the INNER loop, conditional on the case. In the current one, we have fundiffnow>=0 from the OUTTER loop, so we use the varialbe affixed with "max".
		        		    taubalinterp = taubalrbarmin + (taubalrbarmax-taubalrbarmin)* &
		        		    &(rbar-rbarmin)/(rbarmax - rbarmin) ! use the Secant method to obtain a new middle point for the state tax rate interval.
		        		    taubalmin = taubalinterp - tbalwidth*ABS(taubalrbarmax-taubalrbarmin)
		        		    taubalmax = taubalinterp + tbalwidth*ABS(taubalrbarmax-taubalrbarmin)
                            !write(4000+trial_id,fmt='("#27",12x,(17x,i4,x),7(e21.14,x))') bracketr, rbarmax, fundiffmax, rbar, taubalrbarmax, taubalinterp, taubalmin, taubalmax
		        			!END IF
                        else ! unfornunately, we didn't find the target (the rbar) which leads to fundiffnow>=0)
                            if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 25 '
		        			rbarmin = rbar ! update the lower boundary of domain
		        			fundiffmin = fundiffnow ! update the uppper boundary of image
		        			rbar = rbarmin+0.005 ! Because the intermediate search fails, it is still not the time to use the secant method. We use a coarse method to update the input rbar for the outter loop (and because we still have a relatively low interest rate rbar, so should shift UP the interval that emcompasses the real rbar, which God knows.)		    			
		        			!if(dogovloop == 1)then
		        			taubalmax = taubal + pertgov
		        			taubalmin = taubal - 0.0001
		        			taubalrbarmin = taubal
		        			!END IF
                            !write(4000+trial_id,fmt='("#28",12x,6(e21.14,x))') rbarmin, fundiffmin, rbar, taubalmax, taubalmin, taubalrbarmin
		        		endif
		        	elseif(bracketr==22)then
		        		if(fundiffnow<0.0)then !found min
                            if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 26 '
		        			bracketr=0 ! start bisection
		        			rbarmin=rbar          ! update the lower boundary of the outer loop's domain  
		        			fundiffmin=fundiffnow ! update the lower boundary of the outer loop's image
		        		    rbar = rbarmin-(rbarmax-rbarmin)*fundiffmin/(fundiffmax-fundiffmin) ! the first 
		        			!IF (dogovloop==1) THEN
		        			taubalrbarmin=taubal
		        			taubalinterp=taubalrbarmin+(taubalrbarmax-taubalrbarmin)*(rbar-rbarmin)/(rbarmax-rbarmin)
		        			taubalmin=taubalinterp-tbalwidth*ABS(taubalrbarmax-taubalrbarmin)
		        			taubalmax=taubalinterp+tbalwidth*ABS(taubalrbarmax-taubalrbarmin)		    				
		        			!END IF
                            !write(4000+trial_id,fmt='("#29",12x,(17x,i4,x),7(e21.14,x))') bracketr, rbarmin, fundiffmin, rbar, taubalrbarmin, taubalinterp, taubalmin, taubalmax
                        else ! still have a fundiffnow>-0
                            if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 27 '
		        			rbarmax=rbar
		        			fundiffmax=fundiffnow
		        			rbar=rbarmax-0.005
		        			!IF (dogovloop==1) THEN
		        				! define new bounds for taubal	
		        				! if previously found rbarmax and look for rbarmin, use previous eqm taubal as taubalmax
		        				taubalmax=taubal+0.0001
		        				taubalmin=taubal-pertgov
		        				taubalrbarmax=taubal ! change taubalrbarmin to taubalrbarmax. 2-4-2017. Correct. 5-5-2017.
		        			!END IF
                            !write(4000+trial_id,fmt='("#30",12x,6(e21.14,x))') rbarmax, fundiffmax, rbar, taubalmax, taubalmin, taubalrbarmax
		        		end if
		        	else	! if have bounds. 
		        		if(fundiffnow>0.0_wp)then
                            if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 28 '
		        			rbarmax=rbar
		        			fundiffmax=fundiffnow
		        		    rbar=rbarmin-(rbarmax-rbarmin)*fundiffmin/(fundiffmax-fundiffmin)
		        			!if(dogovloop==1)then
		        				! define new bounds for taubal	
		        				! if previously found rbarmax and look for rbarmin, use previous eqm taubal as taubalmax
		        				!taubalmin=taubal-pertgov
		        				!taubalmax=taubal+0.0001
		        			taubalrbarmax=taubal
		        			taubalinterp=taubalrbarmin+(taubalrbarmax-taubalrbarmin)*(rbar-rbarmin)/(rbarmax-rbarmin)
		        			taubalmin=taubalinterp-tbalwidth*ABS(taubalrbarmax-taubalrbarmin)
		        			taubalmax=taubalinterp+tbalwidth*ABS(taubalrbarmax-taubalrbarmin)
                            !write(4000+trial_id,fmt='("#31",12x,(17x,i4,x),7(e21.14,x))') bracketr, rbarmax, fundiffmax, rbar, taubalrbarmax, taubalinterp, taubalmin, taubalmax
		        			!endif
                        else ! fundiffnow<=0
                            if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 29 '
		        			rbarmin=rbar
		        			fundiffmin=fundiffnow
		        		    rbar=rbarmin-(rbarmax-rbarmin)*fundiffmin/(fundiffmax-fundiffmin)
		        			!if(dogovloop==1)then
		        			taubalrbarmin=taubal ! refers to the taubal that leads to a rbar-rimplied < 0.
		        			taubalinterp=taubalrbarmin+(taubalrbarmax-taubalrbarmin)*(rbar-rbarmin)/(rbarmax-rbarmin)
		        			taubalmin=taubalinterp-tbalwidth*abs(taubalrbarmax-taubalrbarmin)
		        			taubalmax=taubalinterp+tbalwidth*abs(taubalrbarmax-taubalrbarmin)
                            !write(4000+trial_id,fmt='("#32",12x,(17x,i4,x),7(e21.14,x))') bracketr, rbarmin, fundiffmin, rbar, taubalrbarmin, taubalinterp, taubalmin, taubalmax
		        			!endif
		        		endif
                    endif
                    !write(unit=102,fmt='("-----",(a,f8.4),(a,i4),3(a,f8.4))') ' epsir   ', epsir, ' iterar   ', iterar, ' rbar ', rbar, ' rimp ', rbarimplied, ' diff ', fundiffnow            
                    if(printout3) write(unit=my_id+1001,fmt='(a,(8x,a),(8x,a),(8x,a))') 'Update: bracketr', 'rbarmax', 'rbarmin', 'newrbar'
                    if(printout3) write(unit=my_id+1001,fmt='(12x,i4,3(f15.7))') bracketr, rbarmax, rbarmin, rbar 
                    if(printout3) write(unit=my_id+1001,fmt='(11x,a,3(11x,a))') 'fundmax','fundmin','tbalmax','tbalmin'
                    if(printout3) write(unit=my_id+1001,fmt='(3x,f15.7,3(3x,f15.7))') fundiffmax, fundiffmin, taubalmax, taubalmin
                endif
                if(printout5) write(unit=104,fmt='(3i3,a)') iterar, iteragov, iteratot, ' 30 '
                iterar = iterar + 1
            else ! exit_log1 == .true.
                exit ! exit the interest rate loop. 
            endif
            
            
            !! 8-1-2017
            new_amin = amin
            new_amax = amax
            call grid_boundary_inspection(new_amin, new_amax, mass_vec, exit_log1) ! 8-28-2017 revision 
            if(exit_log1==.true.) exit !10.17.2017
            write(my_id+1001,'(a)') 'Financial Asset Holdings Distribution '
            do i = 1, adim-1
                write(my_id+1001, '(9x,f11.5)', advance='no') av(i)
            enddo 
            write(my_id+1001, '(9x,f11.5)') av(adim) ! 9-15-2017
            do i = 1, adim-1
                write(my_id+1001, '(12x,a,i3)', advance='no') "level",i
            enddo
            write(my_id+1001, '(12x,a,i3)') "level",adim
            !write(my_id+1001, '(<adim>(f20.9))') av
            write(my_id+1001, '(<adim>(f20.9))') mass_vec
            
            write(my_id+1001, '(2(4x,a,x))') "old amin", "old amax"
            write(my_id+1001, '(2(f12.7,x))') amin, amax 
            write(my_id+1001, '(2(4x,a,x))') "new amin", "new amax"
            write(my_id+1001, '(2(f12.7,x),2/)') new_amin, new_amax
            
            ! housing boundary adjustment
            call grid_housing_upper_bound_adjustment(new_hmax, h_mass_vec, exit_log1) ! 10-9-2017 #2 This line is to get h_mass_vec primarily, not changing the boundary.
            if(exit_log1==.true.) exit !10.17.2017
            write(my_id+1001,'(a)') 'Current distribution of housing assets'
            do i = 1, hdim-1
                write(my_id+1001, '(9x,f11.5)', advance='no') hv(i)     
            enddo 
            write(my_id+1001, '(9x,f11.5)') hv(hdim) ! 9-15-2017
            do i = 1, hdim-1
                write(my_id+1001, '(12x,a,i3)', advance='no') "level", i
            enddo
            write(my_id+1001, '(12x,a,i3)') "level", hdim
            write(my_id+1001, '(<adim>(f20.9),/)') h_mass_vec            
            
            ! boundary update 
            amin = new_amin
            amax = new_amax
            ! grid update
            call set_grids_for_nonhousing_asset(av,3.5_wp,'coarse')
            rav = av ! Bug. 8-1-2017            
            
        end do ! do while loop on rbar ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ [1]
        
        if(.not.exit_log1.and.printout10)then
            !write(unit=127,fmt='((3x,a),3(x,a),(3x,a),(2x,a),(3x,a),(3x,a))') 'NO.', 'Rloop', 'Gloop', 'Tloop', 'rd(5y)', 'implied', 'taubal', 'errgov' 
            !write(unit=127,fmt='(6x,3(i6),4f9.4)') iterar, iteragov, iteratot, rd, rimplied, taubal, epsigov            
            write(unit=my_id+2001,fmt='(/,a,i4,a,i4,a,12f12.7)') ' node_id ', node_id, ' trial_id ', trial_id, ' guess ', guessv
            write(unit=my_id+2001,fmt='(6(2x,a))') 'TotalAsset', ' WokFinAst', ' EntFinAst', 'Captal.CRP', '  LaborCRP', '  exLabEnt'
            write(unit=my_id+2001,fmt='(6f12.6)') totast, wokfin, entfin, crpcap, crplab, entlab
            write(unit=my_id+2001,fmt='(6(2x,a))') 'Captal.SME', '  Wokhouse', '  Enthouse', 'TotalTaxes', ' OutputCRP', ' OutputEnt'
            write(unit=my_id+2001,fmt='(6f12.6)') entcap, wokhom, enthom, tottax, crpprd, entprd 
            !write(unit=127,fmt='(6(2x,a))') 'avg.entcsp', 'avg.wokcsp', 'avg.entfin', 'avg.wokfin', 'avg.enthom', 'avg.wokhom'
            !write(unit=127,fmt='(6f12.4)') mean_entcsp, mean_wokcsp, mean_entfin, mean_wokfin, mean_enthom, mean_wokhom        
            !write(unit=127,fmt='(6(2x,a))') 'avg.entinc', 'avg.wokinc', 'avg.entaxw', 'avg.wokaxw', 'GovBalance', 'govbal2Gdp'
            !write(unit=127,fmt='(6f12.4)') mean_entinc, mean_wokinc, mean_entaxw, mean_wokaxw, govbal, govbal2gdp    
            write(unit=my_id+2001,fmt='(6(2x,a))') 'e2wrat    ', 'w2erat    ', 'entrat    ', 'wokrat   ', 'retiree   ', 'medwokinc '
            write(unit=my_id+2001,fmt='(6f12.6)')  e2wrat, w2erat, entsize, woksize, (1._wp-entsize-woksize), medwokinc 
            write(unit=my_id+2001,fmt='(7(2x,a))') 'totsvt    ', 'tottax    ', 'gfrac*gdp ', 'gfrac     ', 'gdp       ', 'poor20%wok', 'taubal    '
            write(unit=my_id+2001,fmt='(7f12.6)')  totsvt, tottax, gfrac*gdp, gfrac, gdp, lowest_quintile_wokinc,taubal
            write(unit=my_id+2001,fmt='(7(2x,a))') 'dfrac     ','crpcap    ','entcap    ', 'hmin      ', 'hmax      ', 'GovBalance', 'govbal2Gdp'
            write(unit=my_id+2001,fmt='(7f12.6)') dfrac, crpcap, entcap, hmin, hmax, govbal, govbal2gdp
            write(unit=my_id+2001,fmt='(7(2x,a))') 'sumsstax  ','benefit   ', 'medwokinc ', 'taubal    ', 'taubalmax ', 'taubalmin ', 'rimplied  '
            write(unit=my_id+2001,fmt='(7f12.6)') sumsstax, benefit, medwokinc, taubal, taubalmax, taubalmin, rimplied          
            
            call compute_lorenz() 
            write(unit=my_id+2001,fmt='(3(2x,a))') 'csp_gini', 'axw_gini', 'xbi_gini'
            write(unit=my_id+2001,fmt='(3(2x,f8.4))') csp_gini, axw_gini, xbi_gini
                        
        endif ! printout10                  
        
        ! 5-9-2017 feed back added here. 
        ! such as if(exit_log1 == .false.)then
        ! .success.
        ! else
        ! .failure.
        ! endif
        
        ! initial guess, which seems to be underestimated due to the exclusion of capital income
        !! Note: avgincw is defined as the avearge non-entrepreneurial household income, which is changed over iterations.       
        
        !!!! ----------------------------- update zone endding here ----------------------------------------------------------------- !       
        !!!! See the algebra in my appendix
        !!!!staxw = staxwbase*avgincw**(-ptaxw) ! that is the real denominator in the (inc/45000)**p, where 45000 is a rough estimate just for initialization
        !!!!staxe = staxebase*avgincw**(-ptaxe) ! that is the real denominator in the (inc/45000)**p, where 45000 is a rough estimate just for initialization
        !!!!        
        !!!!! boundary 072316
        !!!!hmin = hmin*gdp 
        !!!!hmax = hmax*gdp
        !!!!amax = amax*gdp
        !!!
        !!!
        !!!! grids of housing and wealth determined by (hmin,hmax,amin,amax) 072316
        
        !! -----======######

        !! ## Keep it for future illustration
        !! TO SEE THAT ONLY WEALTH PEOPLE HAS THE ABILITY TO TAKE THE RISK OF RUNNING A BUSINESS
        !! LESS WEALTHY FAMILIES HAVE FEWER RESOURCES TO eke out their living against a NEGATIVE BUSINESS SHOCKS.
        !do i = 0, 1 ! shock
        !    do j = 0, kdim-1 ! project
        !        write(unit=008,fmt='(a,2(i3))') ' --- ', j, i
        !        tar = merge(0._wp, merge(zlow,z2(j),i==0), j==0)
        !        write(unit=008,fmt='(12(f8.2))') (bizmat(m,j,i),m=1,adim)        
        !    enddo
        !enddo        
        !
        !! ## Keep it for future illustration
        !do n = 1, 14
        !    if(n>=9)then ! see the subroutine definition of collateral_constraint_matrix
        !        colbnd1 = 0
        !        colbnd2 = 0
        !    else
        !        colbnd1 = 1 
        !        colbnd2 = nmc
        !    endif
        !    do m = colbnd1, colbnd2
        !        do j = 1, (kdim-1)
        !            do i = 1, hdim
        !                write(unit=024,fmt='(4(a,x,i3,","),"collbd=",i3)') 'h',i,'kp',j,'y',m,'t',n, collbdmat(i,j,m,n) ! (h, kp, yp, t)
        !            enddo
        !        enddo
        !    enddo
        !enddo           
        
        !! EXPERIMENT -- Good
        !allocate( ivec(118950) )
        !ivec = .false.
        !ivec = s2c(:,1)==9
        !!print*, ' count : ', count(ivec)
        !!where(s2c(:,1)==9) ivec = s2c(:,10)
        !allocate( nvec(count(ivec)) )
        !nvec = pack(s2c(:,10),s2c(:,1)==9)
        !call ssi(nvec,'nvec')
        !deallocate( nvec )
        !deallocate( ivec )

        !!# Misc.
        !write(*,fmt='(a)') '------------------------'
        !write(*,'(5x,a,4x,a,4x,a,5x,a)') '1','24','56','7'
        !write(*,'(4i6)') szt18vec1, szt18vec24, szt18vec56, szt18vec7
        !write(*,'(4i6)') szt9vec1, szt9vec24, szt9vec56, szt9vec7
        !write(*,'(4i6)') szt1013vec1, szt1013vec24, szt1013vec56, szt1013vec7
        !write(*,fmt='(a)') '------------------------'
        
        !write(*,0303) 'amax',amax,'hmin',hmin,'hmax',hmax
        !write(*,0004) 'delzh: ', delzh
        !write(*,0004) 'delzl: ', delzl		
        !write(*,0101) 'avgincw',avgincw
        !write(*,0202) 'AggEffLab',AggEffLab,'AggCorpLab',AggCorpLab        
        !write(*,0202) 'AggCorpCap',AggCorpCap,'kndata',kndata
        !write(*,0202) 'AggOut',AggOut,'transfer',transfer
        !write(*,0202) 'rd',rd,'wage',wage
        !write(*,0303) 'hmin',hmin,'hmax',hmax,'amax',amax
        !write(*,0202) 'gamma',gamma, 'rl', rl
        
        deallocate( twa, twh, twk, tww, twf, wf, id ) ! vtwh, vtwk, vtwa
        deallocate( pz2, xl14, bl1014, bd14, xl1013, bldum, bd1013, xl9, bd9, xl18, bd18 )
        !deallocate( fxl14, fbl1014, fbd14, fxl1013, fbldum, fbd1013, fxl9, fbd9, fxl18, fbd18 )
        deallocate( py, yv, yvtemp, sy, pyh, yhv, syh, survprob, popfrac, kv, probtk, phi )
        deallocate( delzh, delzl, efflab, hv, bizmat, z2, delmeh, av, collbdmat, pka, ppldie, sid ) ! pt18
        deallocate( ide, p_homvec ) ! ww, sww wint, afv, nafv, afint
        deallocate( t14vec1, t14vec24, t14vec56, t14vec7 )
        deallocate( c_prf_mat, c_lab_mat, mass_vec, h_mass_vec )
        !deallocate( ppldie_sum )
        deallocate( t18vec1, t18vec24, t18vec56, t18vec7, t9vec1, t9vec24, t9vec56, t9vec7 )
		deallocate( t1013vec1, t1013vec24, t1013vec56, t1013vec7, tvector )
        deallocate( tbim, atwm, taxm, beqm, ldemm, lsupm, prom, ent2wok, wok2ent, homm, cspm )
        deallocate( uwh, uwk, uww, uwa, beqdist, beqdis, nsav )
        deallocate( tbi_m, atw_m, tax_m, beq_m, pro_m, ldem_m, lsup_m, hom_m, csp_m )
        deallocate( totast_m, entcap_m, entlab_m, entprd_m, home_m, ttax_m, uwa_m, uwh_m, c_grs_mat )
        deallocate( s1c, c1s, cwf, cwa, cwh, cwk, cww, sww, swf, swa, swh, swk, swc, sww2, swf2, swa2, swh2, swk2 ) ! 3.31.2017 s2c and c2s is removed.
        deallocate( sw_laborsupply, sw_labordemand, sw_production, sw_bizinvestment, sw_bizloan, sw_ini_asset, sw_ini_house, sw_nonlineartax )
        deallocate( sw_worker_turned, sw_boss_turned, sw_aftertaxwealth, sw_taxableincome, sw_socialsecurity, sw_buzcap_notuse, sw_consumption ) ! csp_lorenz, xbi_lorenz
        deallocate( sw_worker_savtax, sw_entpre_savtax, sw_entpre_biztax, sw_totinc_bx, sw_wealth_tax) ! 3.27.2017
        deallocate( cww2, cwf2, cwc2, cwa2, cwh2, cwk2, rhv, rav, c_lab_vec, c_opt_vec, cwc, cef, ced, cvv, dcef, sef, def ) ! remove scef 8-24-2017 ! remeber to bring it back 10042016
        deallocate( s3c, c3s, swc2, id1 )
        deallocate( term_2, term_3, term_4, term_5, term_6, term_7, term_8, term_9 ) ! debug 9-17-2017

0004    format (a,':',4(f8.3))        
0101    format (a,':',f8.3)        
0202    format (a,':',f8.3,x,a,':',f8.3)        
0303    format (a,':',f8.3,x,a,':',f8.3,x,a,':',f8.3)        
        
    end subroutine solve_model 
    
    subroutine print_basic_model_vectors()
        implicit none
        !if(printout1)then
            call ss(hv,'hv',15,8)
            call ss(av,'av',18,12)
            !call ss(afv,'afv',18,12)
            call ss(kv,'kv',15,8)
            call ss(hv,'hv',15,8)
            call ss(survprob, 'survprob' )
            call ss(efflab,'efflab')  
            call ss(probtk,'probtk')
            call sm(pka,'pka')
            call sm(pz2,'pz2')
            call smi(t18vec1,'t18vec1')
            call smi(t18vec24,'t18vec24')
            call smi(t18vec56,'t18vec56')
            call smi(t18vec7,'t18vec7')
            call smi(t9vec1,'t9vec1')
            call smi(t9vec24,'t9vec24')
            call smi(t9vec56,'t9vec56')
            call smi(t9vec7,'t9vec7')
            call smi(t1013vec1, 't1013vec1')
            call smi(t1013vec24,'t1013vec24')
            call smi(t1013vec56,'t1013vec56')
            call smi(t1013vec7, 't1013vec7')
            call ss(delzh,'delzh1')
            call ss(delzl,'delzl1')              
            call smi(bd1013,'bd1013')
            call smi(xl1013,'xl1013')
            call smi(bd14,'bd14')
            call smi(xl9,'xl9')
            call smi(xl14,'xl14') ! end of period index combinations
            call smi(bd9,'bd9')
            call smi(xl18,'xl18')
            call smi(bd18,'bd18')
            !call smi(fbd1013,'fbd1013')
            !call smi(fxl1013,'fxl1013')
            !call smi(fbd14,'fbd14')
            !call smi(fxl9,'fxl9')
            !call smi(fbd9,'fbd9')
            !call smi(fxl18,'fxl18')
            !call smi(fbd18,'fbd18')   
            !call smi(afint,'afint')  
            !call smi(s2c,'s2c',8)          
            call smi(s3c,'s3c',8)
            call ss(rhv,'rhv')
            call ss(rav,'rav',18,12)
            call ss(popfrac,'popfrac')
            call ss(survprob,'survprob')
        !endif ! printout1            
    end subroutine print_basic_model_vectors
    
subroutine linear_combination_sobal_sequence(linear_combination_sobol, sequence_index, sobol_scaled_input, origin_input, weight_sobol, break_sobol)    
    implicit none
    real(wp), dimension(:), intent(out) :: linear_combination_sobol
    real(wp), dimension(:), intent(in) :: weight_sobol, origin_input, sobol_scaled_input
    integer, dimension(:), intent(in) :: break_sobol
    integer, intent(in) :: sequence_index
    integer :: idx
    real(wp) :: wgt
    idx = locate(real(break_sobol,wp),real(sequence_index,wp),1)
    wgt = weight_sobol(idx)
    linear_combination_sobol = wgt*origin_input + (1._wp-wgt)*sobol_scaled_input 
    !write(*,'(a,f6.3,a,i3,a,i3)') ' weight ', wgt, ' loc ', idx, ' list index ', sequence_index
    !write(*,*) ' --------------------------------------------- '
end subroutine linear_combination_sobal_sequence

subroutine linear_combination_sobol_sequence_corrected( vertex_output, sequence_index, unit_sobol_input, range_guess )
    real(wp), dimension(:), intent(out) :: vertex_output
    real(wp), dimension(:), intent(in) :: unit_sobol_input
    real(wp), dimension(:,:), intent(in) :: range_guess
    integer, intent(in) :: sequence_index
    ! scale from [0,1] to [range_guess(:,1),range_guess(:,2)]
    vertex_output = (range_guess(:,2)+range_guess(:,1))/2._wp + (unit_sobol_input-0.5_wp)*(range_guess(:,2)-range_guess(:,1))
    
end subroutine linear_combination_sobol_sequence_corrected
    
end module equilibrium
    
    
    