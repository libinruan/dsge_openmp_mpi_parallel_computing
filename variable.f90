module variable
    use toolbox
    implicit none
    character(len=30) :: labstr(145)
    real(wp) :: para(145), & ! total number of parameters in _lparameter.txt excluding boolin variables (printout1, etc).
                targetv(10), & ! target vector
                guessv(10), & ! a guess on the parameter setting mainly for mpi_exercise_mode == 0 case.
                momvec(10), & ! simulated moment vector
                startpoint(25), & ! read starting point from the result of the first stage.
                rhoy = 0.715, & ! 5 year
                vary = 0.52, & ! 5 year
                rhoyh = 0.677, & ! 5 year
                varyh = 0.37, & ! 5 year
                deltak = 0.109, & ! 1 year Nakajima non-housing capital; Meh, 0.062
                deltah = 0.017, & ! 1 year Nakajima housing capitals
                rd = 0.035, & ! 1 year
                alpha = 0.26, &
                lambda = 0.2, & ! down payment ratio. 
                theta = 0.2, & ! Share of housing consumption # 7
                sigma = 2.0, &
                nu = 0.88, & ! parameter of diminishing return to scale, to be calibrated #
                kv1 = 0.5, & ! to be calibrated # 1
                prtk0 = 0.024, & ! to be calibrated # 2
                prtk1 = 0.110, & ! to be calibrated # 3          
                prtk2 = 0.075, & ! to be calibrated # 4
                prtk3 = 0.000, & ! No innovation. Wired. Don't change it.  
                phi0 = 0.5, & ! L: my idea; seems useless        
                phi1 = 0.75, & ! Set # 8
                phi2 = 0.92, & ! Set as in Meh (2005) # 9
                phi3 = 0.97, & ! Set # 10
                beta = 0.93, & ! on annual basis. # 5
                zbar = 2.0, & ! # 6. to be calibrated # 2.374. could bring the labor employment positive in both sectors 1.8
                zlow = 0.0, & ! to be calibrated #
                zeta = 0.1, & ! house transaction cost parameter. replaced with Yang's setting.
                hmin = 1.402, & ! Yang: 1.402  
                hmax = 120.0, & ! 50 Yang: 17. // 17*aggY
                amin = 0.0, & ! 50.0, & ! Yang: 0. [Note] amin is not allowed to be positive. Otherwise the inserting-zero-into-av-vector operation in equilibrium.f90 fails.
                amax = 999.0, & ! Yang: 100. // 95*aggY
                nmax, &
                iota = 0.5, &
                dfrac = 0.27, & ! gov debt / total capital
                gfrac = 0.187, & ! gov expense / gdp OK.
                taubal = 0.033, &
                CorpLabFrac = 0.5, &
                CorpOutFrac = 0.5, & ! https://www.sba.gov/sites/default/files/rs390tot_0.pdf
                btaxw = 0.32, &   !! workers' parameters of the Gouveia Strauss effective federal individual income tax function
                ptaxw = 0.7646, & ! see Cagetti and De Nardi, 2009 Estate Taxation, p. 23  checked 9-5-2016 and 1-29-2017  
                staxw = 0.2154, & ! 
                btaxe = 0.2562, & !! entrepreneur' parameters of the G-S effective federal individual income tax function 
                ptaxe = 1.4, &    !     
                staxe = 0.42, &   !
                btax  = 0.258, & ! (a0) 
                ptax  = 0.768, & ! (a1)
                stax  = 0.031, & ! (a2) used by Kitao (0.438). Naka uses 0.031. Needs to update during global optimization.
                tausv = 0.15, & ! rate of tax on interest income
                taubp = 0.25, & ! rate of tax on business profit
                tauss = 0.074, &  ! to be calibrated # 
                PropHouseCapital = 0.5, &
                BossPropPopu = 0.12, & ! literature
                rl = 0.09, &
                penalty = -1.e6, & ! penalty for initializing matrices or performance evaluation
                epsilon = 1.e-4, & ! shift in housing service flow utility   
                wage, &        
                gamma5, &
                err_dist = 1.e-8, &
                err_svdist = 1.e-8, & ! convergence criteria.
                errdist, &
                rimplied, &
                wageimplied, &
                transbeqimplied, &
                benefitimplied, &
                entsize13, &
                woksize19, &
                gamma1 = 0.0505, &
                mu1 = 0.07, & ! maximum depreciation (1.0 in Yang's program)
                mu2 = 0.07, & ! maximum renovation (1.0 in Yang's program)
                rho1 = 0.070, & ! transaction costs of selling a house
                rho2 = 0.025, & ! transaction costs of buying a house
                b2gratio = 0.02, &
                transbeq, &
                !transbeq_new, &
                nmul = 5, &
                tauk = 0.4, &
                gsdim = 3.5, & ! the grid parameter for small region Brent optimization.	
                trsld = -1.e3, & ! expected value is reassigned as penalty level if its original value is smaller than "trsld"
                maxdist = 5.e-6, &
                tinymass = 1.e-2, &
                chunkmass = 3.e-2, &
                initau = 4._wp, &
                length = 5._wp, &
                accumass = 1.e10, &
                accupara = 1.e5, &
                errvalfc = 2.e-3, &
                momround = 1e3 ! round to the 3rd place to the right of decimal point
                
    integer ::  nmc = 3, &   
                itert = 1, &
                hdim = 10, & 
                adim = 12, &
                fnadim = 15, & ! fnadim: number of grids in the finer grid 
                fnhdim = 10, &
                kdim = 4, &
                ndim = 10, & ! number of targets
                nsbq = 10000, & ! total number of raw sobol sequence's rows
                ! length = 5, &
                iterasvmax = 10, &
                sdim   = 5, &
                nsdim = 100, &
                subdim = 5, & !10.15.2017 a subset of the original 10 dimensions
                iterainvdist = 50, &
                noamoeba = 2, &
                testfunc_idx = 1, &
                setindex = 1, &
                listlength = 10000, & ! the length of selective points after coarse search 10.15.2017
                srnumber, &
                obj_func_toggle
    
    character(len=100) :: listnumber
    
    ! ## outer loop's first task: updating the estimate of rbar ##
    real(wp) :: radjust  = 0.08 ! the maximum distance allows new rbar to take rimplied.
    real(wp) :: rweight  = 0.50 ! if rimplied is far away from rbar, interpolate with this parameter for new rbar
    real(wp) :: epsir    = 1.0  ! initialize the error on rbar (absolute distance between input and output r) 
    real(wp) :: epsirmin = 5.e-05 ! the criterion for rbar loop convergence
    integer  :: iterar   = 1
    integer  :: iteratot = 1
    integer  :: iterarmax = 40 ! maximum number of iterations on rbar
    real(wp) :: rbar ! = 0.035, this value is hard wired in equilibirum.f90 line 545.
    real(wp) :: rbarmax ! [domain]
    real(wp) :: rbarmin ! [domain]
    real(wp) :: fundiffnow ! [image]
    real(wp) :: fundiffmax ! [image]
    real(wp) :: fundiffmin ! [image]
    integer  :: bracketr ! round indicator for judging searching direction
    real(wp) :: tbalwidth  = 0.3 ! **After found the two boundary points of outer loop, with this deviation parameter, we update the interval of taubal in inner loop
    
    ! ## outer loop's second task: updating the estimate interval of taubal in inner loop ##
    real(wp) :: taubalrbarmax ! when rbar>rimplied in outer loop, we update the estimate of inner loop's interval of taubal.
    real(wp) :: taubalrbarmin ! when rbar<rimplied in outer loop, we update the estimate of inner loop's interval of taubal.
    real(wp) :: taubalinterp  ! used with the secant method; **map the outer loop guess rbar into a better inner loop guess taubal
    
    ! ## inner loop: updating taubal ##
    real(wp) :: epsigov    = 1.0
    real(wp) :: epsigovmin = 5.e-05
    integer  :: iteragov   = 1
    integer  :: iteragovmax = 40 
    integer  :: iteragov8rate
    real(wp) :: taubalmin  = 0.0345 ! [domain]
    real(wp) :: taubalmax  = 0.0363 ! [domain]
    !real(wp) :: govbal2gdp    ! [image]
    real(wp) :: govbal2gdpmax ! [image]
    real(wp) :: govbal2gdpmin ! [image]
    integer  :: noneedtaubalmax ! =1, if found a govbal2gdp>=0
    integer  :: bracketgov ! round indicator for judging searching direction
    real(wp) :: pertgov = 0.04_wp
    real(wp) :: rbarimplied
    
    ! real(wp) :: penalty = -1e+7 
    real(wp) :: inf ! defined in toolbox.f90's function 'infinity_setting'.
    
    integer ::  t, inv_dist_counter, szperiod1
    real(wp) :: new_amin, new_amax, new_hmax
    real(wp) :: staxbase, staxwbase, staxebase, kndata, avgincw, govbalance, govbal, benefit, sumsstax
    real(wp) :: AggEffLab, AggCorpLab, AggCorpCap, AggCorpOut, AggOut, AggTax, AggAst, AggInc
    real(wp) :: totast, entcap, crpcap, labsup, entlab, crplab, crpprd, entprd, totsvt, totbpt
    real(wp) :: chge2w, entsize, entinc, enttax, enthom, statax, entaxw ! ___axw: after tax wealth
    real(wp) :: chgw2e, woksize, wokinc, woktax, wokhom, fedtax, wokaxw, w2erat, e2wrat
    real(wp) :: mean_entinc, entcsp, mean_entcsp, mean_enthom, entsize1to13, mean_entaxw
    real(wp) :: mean_wokinc, wokcsp, mean_wokcsp, mean_wokhom, woksize1to13, mean_wokaxw
    real(wp) :: entfin, wokfin, mean_entfin, mean_wokfin, tottax, gdp, govbal2gdp
    real(wp) :: totefl, totexl, axw_gini, csp_gini, xbi_gini, poppaysstaximplied, medwokinc, lowest_quintile_wokinc
    real(wp) :: crp_cap_per, hug_inv_per, med_inv_per, sml_inv_per, ent_pop_per
    real(wp) :: fin2gdp_rat, hom2gdp_rat, ent_wel_per, ent_inc_per, med_wel_e2w
    real(wp) :: hug_inv_proj, med_inv_proj, sml_inv_proj, all_inv_proj
    real(wp) :: medentwel, medwokwel, all_income
    
    real(wp), dimension(:), allocatable :: yv, yhv, sy, syh, survprob, delmeh
    real(wp), dimension(:), allocatable :: popfrac, kv, probtk, phi
    real(wp), dimension(:), allocatable :: z2, delzl, delzh, efflab
    real(wp), dimension(:), allocatable :: hv, av, rhv, rav, mass_vec, h_mass_vec
    !real(wp), dimension(:), allocatable :: afv ! asset refined grid ! 4.1.2017
    real(wp), dimension(:,:), allocatable :: py, pyh, pz2, pka, pyb ! wint ! wint, aint: interpolated weight and coarse grids ! pt18
    real(wp), dimension(:,:), allocatable :: range_guess
    real(wp), dimension(:), allocatable :: term_2, term_3, term_4, term_5, term_6, term_7, term_8, term_9 ! debug 9-17-2017
    integer,  dimension(:,:), allocatable :: xl14, bl1014, bd14, xl1013, xl9, xl18, bd9, bd18 ! 4.1.2017 afint
    integer,  dimension(:,:), allocatable :: bldum, bd1013
    integer,  dimension(:,:), allocatable :: s1c, s3c ! s2c is removed 3.31.2017
    logical,  dimension(:), allocatable :: ivec
    integer,  dimension(:), allocatable :: nvec
    
    !integer,  dimension(:,:), allocatable :: fxl14, fbl1014, fbd14, fxl1013, fxl9, fxl18, fbd9, fbd18
    !integer,  dimension(:,:), allocatable :: fbldum, fbd1013          
    
    integer, dimension(:,:,:,:,:,:), allocatable :: c1s
    real(wp), dimension(:,:,:,:,:,:), allocatable :: ced
    integer, dimension(:,:,:,:,:,:,:,:,:), allocatable :: c3s ! c2s is removed 3.31.2017.
    integer,  dimension(:,:,:,:,:,:,:,:,:), allocatable :: twh, twk, tww !, vtwh, vtwk ! vtwh, vtwk, tww: indicator to switch to be a labor in the beginning of the period. 1 as switches.
    real(wp), dimension(:,:,:,:,:,:,:,:,:), allocatable :: twa, twf !, vtwa ! Trick: if you index the sole element with zero in an allocatable matrix, be sure to use colon sign. for example, 0:0.
    real(wp), dimension(:,:,:,:,:,:,:,:,:), allocatable :: cwa, cwf, cwh, cwc, cef, dcef !, scef ! 8-24-2017 Debug
    integer,  dimension(:,:,:,:,:,:,:,:,:), allocatable :: cwk, cww
    integer,  dimension(:,:,:,:,:,:), allocatable :: cvv
    
    real(wp), dimension(:,:,:,:,:,:,:,:,:), allocatable :: cwa2, cwf2, cwh2, cwc2
    integer,  dimension(:,:,:,:,:,:,:,:,:), allocatable :: cwk2, cww2    
    
    real(wp), dimension(:,:,:,:,:,:,:,:,:), allocatable :: uwa, uwh
    integer,  dimension(:,:,:,:,:,:,:,:,:), allocatable :: uwk, uww

    !real(wp), dimension(:,:,:,:,:,:,:,:,:), allocatable :: rwa, rwf, rwh
    !integer,  dimension(:,:,:,:,:,:,:,:,:), allocatable :: rwk, rww        
    
    integer, dimension(:), allocatable :: sww, swk
    real(wp), dimension(:), allocatable :: swf, swa, swh, swc, sef, def, sef1, sef2, sef3! sef1 and sef2 stationary distribution used in subrtouine ability_transition of model.f90.
    real(wp), dimension(:), allocatable :: sw_laborsupply, sw_labordemand, sw_production, sw_bizinvestment, sw_bizloan
    real(wp), dimension(:), allocatable :: sw_ini_asset, sw_ini_house, sw_nonlineartax, sw_aftertaxwealth, sw_socialsecurity
    real(wp), dimension(:), allocatable :: sw_buzcap_notuse, sw_worker_savtax, sw_entpre_savtax, sw_entpre_biztax ! 3.27.2017 add savtax and biztax.
    real(wp), dimension(:), allocatable :: sw_worker_turned, sw_boss_turned, sw_taxableincome, sw_consumption, sw_totinc_bx
    real(wp), dimension(:), allocatable :: axw_lorenz, csp_lorenz, xbi_lorenz 
    
    integer, dimension(:), allocatable :: sww2, swk2
    real(wp), dimension(:), allocatable :: swf2, swc2, swa2, swh2    

    !real(wp), dimension(:), allocatable :: labdem
    
    real(wp), dimension(:,:,:,:,:,:,:,:,:), allocatable :: tbi_m, atw_m, tax_m, beq_m, pro_m, ldem_m, lsup_m, hom_m, csp_m ! ON THE REFINED GRID
    real(wp), dimension(:,:,:,:,:,:,:,:,:), allocatable :: tbim, atwm, taxm, beqm, prom, ldemm, lsupm, homm, cspm ! ON THE COARSE GRID. tbim: taxible income, atwm: after tax wealth. production. labor demand and supply. 
    real(wp), dimension(:,:,:,:,:,:,:,:,:), allocatable :: ide, ent2wok, wok2ent ! 4.1.2017 nafv ! (nsav)-->(nafv), nsav: the beginning period refined grid; nafv: the end period refined grid
    real(wp), dimension(:,:,:,:,:,:,:,:,:), allocatable :: totast_m, entcap_m, entlab_m, entprd_m, home_m, ttax_m, uwa_m, uwh_m !, entsize_m
    real(wp), dimension(:,:), allocatable :: p_homvec, distance_vec
    
    real(wp), dimension(:,:,:,:,:,:), allocatable :: wf, id, ppldie, idtemp, sid ! ppldie: mass across states and elderly periods
    real(wp), dimension(:,:,:,:,:), allocatable :: c_prf_mat, c_lab_mat, c_grs_mat, beqdis, id1 ! invariant distribution of period one.
    real(wp), dimension(:,:,:,:,:,:), allocatable :: beqdist, nsav ! beqdis + t means it is with extra state variable: time. nsav: the begining period refined grid.
    real(wp), dimension(:,:,:), allocatable :: bizmat
    real(wp), dimension(:,:,:,:,:,:,:), allocatable :: inc1013
    real(wp), dimension(:), allocatable :: c_lab_vec, c_opt_vec
    !real(wp), dimension(:,:,:,:,:,:), allocatable :: tranmat
    integer, dimension(:,:,:,:), allocatable :: collbdmat 
    integer, dimension(:,:), allocatable :: tvector
	
    ! end period state combination
    integer, dimension(:,:), allocatable :: t18vec1, t18vec24, t18vec56, t18vec7
	integer, dimension(:,:), allocatable :: t9vec1, t9vec24, t9vec56, t9vec7
	integer, dimension(:,:), allocatable :: t1013vec1, t1013vec24, t1013vec56, t1013vec7
    integer, dimension(:,:), allocatable :: t14vec1, t14vec24, t14vec56, t14vec7
    ! range in a local block
    integer :: szt18vec1, szt18vec24, szt18vec56, szt18vec7
	integer :: szt9vec1, szt9vec24, szt9vec56, szt9vec7
	integer :: szt1013vec1, szt1013vec24, szt1013vec56, szt1013vec7
    integer :: szt14vec1, szt14vec24, szt14vec56, szt14vec7
    
    ! amoeba - general
    integer, dimension(:), allocatable :: indexseries    
    integer :: mpi_exercise_mode, trylen, sblno1, nrow, amoitrcrt, amoiniswitch, nslaves
    integer :: trial, slave, ierr, msgtype, modelmsg
    
    real(wp) :: stepsize, reinifac, amoalp, amogam, amobet, amotau, amoerrcrt, amoconvex
    real(wp) :: obj_val_1st
    real(wp), dimension(:), allocatable :: parcel, result, obj_val_vec, origin_input, selected_input, pt_input_ndim
    real(wp), dimension(:,:), allocatable :: sobolm, sobolm_scaled, mpi_sobol_scaled, mpi_sobol_mixed ! sobolm: scaled for ranges; mpi_sobol_scaled: adjusted for the starting point.
    real(wp), dimension(:,:), allocatable :: mpi_simmom_matrix, outputinput1, pointlist, pts_ndim, pts_subdim
    
    ! amoeba - break_list - mpi_exercise_mode == 0 Only.
    !real(wp), dimension(10) :: weight_list = [0.95_wp,0.9_wp,0.85_wp,0.8_wp,0.75_wp,0.7_wp,0.65_wp,0.6_wp,0.55_wp,0.5_wp] ! used for mpi_exercise_mode==1
    !integer, dimension(10) :: breaks_list = [1,101,201,301,401,501,601,701,801,901] ! used for mpi_exercise_mode==1
    real(wp), dimension(2) :: weight_list = [0._wp,0._wp] ! 8-18-2017 ! Useless now.
    integer, dimension(2) :: breaks_list = [1,50000]      ! 8-18-2017 ! Useless now.
    
    logical :: printout1, printout2, printout3, printout4, printout5, printout6, printout7, printout8, printout9, printout10, printout11, printout12 !, tausvflag
    logical :: printout13, printout14, printout15, printout16, printout17, printout18, printout19, printout20, printout21, printout22, printout23, printout24
    logical :: printout25
    logical :: receiving, status(mpi_status_size)
    character(len=50) :: node_string, trylen_string, amoeba_x_y_string, bestvertex_file
    character(:), allocatable :: solution_string, io_string, concisesolution_string
    
    ! amoeba - mpi_exercise_mode == 2
    integer :: ncol, irow, icol, AMOEBA_ID, GROUP_AMOEBA, DUMMY_WORLD, AMOEBA_WORLD, CONTRACT_ID, CONTRACT_WORLD, GROUP_CONTRACT, HEAD_ID
    integer :: HEAD_WORLD, GROUP_HEAD, sirow, sicol, slist, elist
    integer :: ROW_WORLD, ROW_ID, GROUP_ROW
    integer, dimension(:), allocatable :: id_list, refpts
    real(wp), dimension(:), allocatable :: bestvertex
    real(wp) :: bestobjval
    
contains  
    subroutine read_parameter_model( para, input_file )
        implicit none
        real(wp), dimension(:), intent(inout) :: para
        character(len=*), intent(in) :: input_file
        character(len=80) :: name_string, value_string
        integer :: n, iostat, i
        n = size( para )
        para = 0._wp
        i = 0
        open( unit=10, file=input_file, status='old', action='read', iostat=iostat )
        if( iostat==0 ) then
            do 
                read( unit=10, fmt=*, iostat=iostat ) name_string, value_string
                if( iostat/=0 ) exit
                if( scan( name_string, '!')>0 ) cycle
                select case( name_string )
                   case ('printout1') 
                        read( value_string, * ) printout1      
                   case ('printout2') 
                        read( value_string, * ) printout2    
                   case ('printout3') 
                        read( value_string, * ) printout3     
                   case ('printout4') 
                        read( value_string, * ) printout4      
                   case ('printout5') 
                        read( value_string, * ) printout5  
                   case ('printout6') 
                        read( value_string, * ) printout6        
                   case ('printout7') 
                        read( value_string, * ) printout7    
                   case ('printout8') 
                        read( value_string, * ) printout8       
                   case ('printout9') 
                        read( value_string, * ) printout9   
                   case ('printout10') 
                        read( value_string, * ) printout10   
                   case ('printout11') 
                        read( value_string, * ) printout11  
                   case ('printout12') 
                        read( value_string, * ) printout12         
                   case ('printout13') 
                        read( value_string, * ) printout13  
                   case ('printout14')                        
                        read( value_string, * ) printout14
                   case ('printout15')
                        read( value_string, * ) printout15
                   case ('printout16')     
                        read( value_string, * ) printout16
                   case ('printout17')     
                        read( value_string, * ) printout17                        
                   case ('printout18')     
                        read( value_string, * ) printout18          
                   case ('printout19')     
                        read( value_string, * ) printout19   
                   case ('printout20')     
                        read( value_string, * ) printout20    
                   case ('printout21')     
                        read( value_string, * ) printout21      
                   case ('printout22')     
                        read( value_string, * ) printout22      
                   case ('printout23')     
                        read( value_string, * ) printout23 
                   case ('printout24')     
                        read( value_string, * ) printout24     
                   case ('printout25')     
                        read( value_string, * ) printout25                        
                   case ('listnumber')
                        read(value_string, *  ) listnumber
                   case ('bestvertex_file')
                        read(value_string, * ) bestvertex_file
                   case('mpi_exercise_mode')
                       read( value_string,*) mpi_exercise_mode                         
                   !case ('tausvflag')
                   !     read( value_string, * ) tausvflag
                   case('rhoy') ! 1
                       i = i + 1 
                       read( value_string,*) rhoy 
                       labstr(i) = 'rhoy'
                       para(i) = rhoy
                   case('vary') ! 2
                       i = i + 1 
                       read( value_string,*) rhoy 
                       labstr(i) = 'vary'
                       para(i) = rhoy
                   case('rhoyh') ! 3
                       i = i + 1 
                       read( value_string,*) rhoyh 
                       labstr(i) = 'rhoyh'
                       para(i) = rhoyh
                   case('varyh') ! 4
                       i = i + 1 
                       read( value_string,*) varyh 
                       labstr(i) = 'varyh'
                       para(i) = varyh
                   case('deltak') ! 5
                       i = i + 1 
                       read( value_string,*) deltak 
                       labstr(i) = 'deltak'
                       para(i) = deltak
                   case('deltah') ! 6
                       i = i + 1 
                       read( value_string,*) deltah 
                       labstr(i) = 'deltah'
                       para(i) = deltah
                   case('rd') ! 7
                       i = i + 1 
                       read( value_string,*) rd 
                       labstr(i) = 'rd'
                       para(i) = rd
                   case('alpha') ! 8
                       i = i + 1 
                       read( value_string,*) alpha 
                       labstr(i) = 'alpha'
                       para(i) = alpha
                   case('lambda') ! 9 happened only once in the model. the collateral constraint.
                       i = i + 1 
                       read( value_string,*) lambda 
                       labstr(i) = 'lambda'
                       para(i) = lambda
                   case('theta') ! 10
                       i = i + 1 
                       read( value_string,*) theta 
                       labstr(i) = 'theta'
                       para(i) = theta
                   case('sigma') ! 11
                       i = i + 1 
                       read( value_string,*) sigma 
                       labstr(i) = 'sigma'
                       para(i) = sigma
                   case('nu') ! 12
                       i = i + 1 
                       read( value_string,*) nu 
                       labstr(i) = 'nu'
                       para(i) = nu
                   case('kv1') ! 13
                       i = i + 1 
                       read( value_string,*) kv1 
                       labstr(i) = 'kv1'
                       para(i) = kv1
                   case('prtk0') ! 14
                       i = i + 1 
                       read( value_string,*) prtk0 
                       labstr(i) = 'prtk0'
                       para(i) = prtk0
                   case('prtk1') ! 15
                       i = i + 1 
                       read( value_string,*) prtk1 
                       labstr(i) = 'prtk1'
                       para(i) = prtk1
                   case('prtk2') ! 16
                       i = i + 1 
                       read( value_string,*) prtk2 
                       labstr(i) = 'prtk2'
                       para(i) = prtk2
                   case('phi0') ! 17
                       i = i + 1 
                       read( value_string,*) phi0 
                       labstr(i) = 'phi0'
                       para(i) = phi0
                   case('phi1') ! 18
                       i = i + 1 
                       read( value_string,*) phi1 
                       labstr(i) = 'phi1'
                       para(i) = phi1
                   case('phi2') ! 19
                       i = i + 1 
                       read( value_string,*) phi2 
                       labstr(i) = 'phi2'
                       para(i) = phi2
                   case('phi3') ! 20
                       i = i + 1 
                       read( value_string,*) phi3 
                       labstr(i) = 'phi3'
                       para(i) = phi3
                   case('beta') ! 21
                       i = i + 1 
                       read( value_string,*) beta 
                       labstr(i) = 'beta'
                       para(i) = beta
                   case('zbar') ! 22
                       i = i + 1 
                       read( value_string,*) zbar 
                       labstr(i) = 'zbar'
                       para(i) = zbar
                   case('zlow') ! 23
                       i = i + 1 
                       read( value_string,*) zlow 
                       labstr(i) = 'zlow'
                       para(i) = zlow
                   case('hmin') ! 24
                       i = i + 1 
                       read( value_string,*) hmin 
                       labstr(i) = 'hmin'
                       para(i) = hmin
                   case('hmax') ! 25
                       i = i + 1 
                       read( value_string,*) hmax 
                       labstr(i) = 'hmax'
                       para(i) = hmax
                   case('amin') ! 26
                       i = i + 1 
                       read( value_string,*) amin 
                       labstr(i) = 'amin'
                       para(i) = amin
                   case('amax') ! 27
                       i = i + 1 
                       read( value_string,*) amax 
                       labstr(i) = 'amax'
                       para(i) = amax
                   case('dfrac') ! 28
                       i = i + 1 
                       read( value_string,*) dfrac 
                       labstr(i) = 'dfrac'
                       para(i) = dfrac 
                   case('gfrac') ! 29
                       i = i + 1 
                       read( value_string,*) gfrac 
                       labstr(i) = 'gfrac'
                       para(i) = gfrac
                   case('taubal') ! 30
                       i = i + 1 
                       read( value_string,*) taubal 
                       labstr(i) = 'taubal'
                       para(i) = taubal
                   case('CorpLabFrac') ! 31 
                       i = i + 1 
                       read( value_string,*) CorpLabFrac 
                       labstr(i) = 'CorpLabFrac'
                       para(i) = CorpLabFrac
                   case('CorpOutFrac') ! 32
                       i = i + 1 
                       read( value_string,*) CorpOutFrac 
                       labstr(i) = 'CorpOutFrac'
                       para(i) = CorpOutFrac
                   case('btaxw') ! 33
                       i = i + 1 
                       read( value_string,*) btaxw 
                       labstr(i) = 'btaxw'
                       para(i) = btaxw
                   case('ptaxw') ! 34
                       i = i + 1 
                       read( value_string,*) ptaxw 
                       labstr(i) = 'ptaxw'
                       para(i) = ptaxw
                   case('staxw') ! 35
                       i = i + 1 
                       read( value_string,*) staxw 
                       labstr(i) = 'staxw'
                       para(i) = staxw
                   case('btaxe') ! 36
                       i = i + 1 
                       read( value_string,*) btaxe 
                       labstr(i) = 'btaxe'
                       para(i) = btaxe
                   case('ptaxe') ! 37
                       i = i + 1 
                       read( value_string,*) ptaxe 
                       labstr(i) = 'ptaxe'
                       para(i) = ptaxe
                   case('staxe') ! 38
                       i = i + 1 
                       read( value_string,*) staxe 
                       labstr(i) = 'staxe'
                       para(i) = staxe
                   case('tauss') ! 39
                       i = i + 1 
                       read( value_string,*) tauss 
                       labstr(i) = 'tauss'
                       para(i) = tauss
                   case('PropHouseCapital') ! 40
                       i = i + 1 
                       read( value_string,*) PropHouseCapital 
                       labstr(i) = 'PropHouseCapital'
                       para(i) = PropHouseCapital
                   case('BossPropPopu') ! 41
                       i = i + 1 
                       read( value_string,*) BossPropPopu 
                       labstr(i) = 'BossProfPopu'
                       para(i) = BossPropPopu
                   case('rl') ! 42
                       i = i + 1 
                       read( value_string,*) rl 
                       labstr(i) = 'rl'
                       para(i) = rl
                   case('penalty') ! 43 
                       i = i + 1 
                       read( value_string,*) penalty 
                       labstr(i) = 'penalty'
                       para(i) = penalty
                   case('epsilon') ! 44
                       i = i + 1 
                       read( value_string,*) epsilon 
                       labstr(i) = 'epsilon'
                       para(i) = epsilon
                   case('err_dist') ! 45
                       i = i + 1 
                       read( value_string,*) err_dist 
                       labstr(i) = 'err_dist'
                       para(i) = err_dist
                   case('gamma1') ! 46
                       i = i + 1 
                       read( value_string,*) gamma1 
                       labstr(i) = 'gamma1'
                       para(i) = gamma1
                   case('mu1') ! 47
                       i = i + 1 
                       read( value_string,*) mu1 
                       labstr(i) = 'mu1'
                       para(i) = mu1
                   case('mu2') ! 48
                       i = i + 1 
                       read( value_string,*) mu2 
                       labstr(i) = 'mu2'
                       para(i) = mu2
                   case('rho1') ! 49
                       i = i + 1 
                       read( value_string,*) rho1 
                       labstr(i) = 'rho1'
                       para(i) = rho1
                   case('rho2') ! 50
                       i = i + 1 
                       read( value_string,*) rho2 
                       labstr(i) = 'rho2'
                       para(i) = rho2
                   case('nmc') ! 51
                       i = i + 1 
                       read( value_string,*) nmc 
                       labstr(i) = 'nmc'
                       para(i) = nmc
                   case('hdim') ! 52
                       i = i + 1 
                       read( value_string,*) hdim 
                       labstr(i) = 'hdim'
                       para(i) = hdim
                   case('adim') ! 53
                       i = i + 1 
                       read( value_string,*) adim 
                       labstr(i) = 'adim'
                       para(i) = adim
                   case('fnadim') ! 54
                       i = i + 1 
                       read( value_string,*) fnadim 
                       labstr(i) = 'fnadim'
                       para(i) = fnadim
                   case('kdim') ! 55
                       i = i + 1 
                       read( value_string,*) kdim 
                       labstr(i) = 'kdim'
                       para(i) = kdim
                   case('length') ! 56
                       i = i + 1 
                       read( value_string,*) length 
                       labstr(i) = 'length'
                       para(i) = length
                   case('radjust') ! 57
                       i = i + 1 
                       read( value_string,*) radjust 
                       labstr(i) = 'radjust'
                       para(i) = radjust
                   case('rweight') ! 58
                       i = i + 1 
                       read( value_string,*) rweight 
                       labstr(i) = 'rweight'
                       para(i) = rweight
                   case('epsirmin') ! 59
                       i = i + 1 
                       read( value_string,*) epsirmin 
                       labstr(i) = 'epsirmin'
                       para(i) = epsirmin
                   case('iterarmax') ! 60
                       i = i + 1 
                       read( value_string,*) iterarmax 
                       labstr(i) = 'iterarmax'
                       para(i) = iterarmax
                   case('rbar') ! 61
                       i = i + 1 
                       read( value_string,*) rbar 
                       labstr(i) = 'rbar'
                       para(i) = rbar
                   case('tbalwidth') ! 62 
                       i = i + 1 
                       read( value_string,*) tbalwidth 
                       labstr(i) = 'tbalwidth'
                       para(i) = tbalwidth
                   case('epsigovmin') ! 63
                       i = i + 1 
                       read( value_string,*) epsigovmin 
                       labstr(i) = 'epsigovmin'
                       para(i) = epsigovmin
                   case('iteragovmax') ! 64
                       i = i + 1 
                       read( value_string,*) iteragovmax 
                       labstr(i) = 'iteragovmax'
                       para(i) = iteragovmax
                   case('taubalmin') ! 65
                       i = i + 1 
                       read( value_string,*) taubalmin 
                       labstr(i) = 'taubalmin'
                       para(i) = taubalmin
                   case('taubalmax') ! 66
                       i = i + 1 
                       read( value_string,*) taubalmax 
                       labstr(i) = 'taubalmax'
                       para(i) = taubalmax
                   case('pertgov') ! 67
                       i = i + 1 
                       read( value_string,*) pertgov 
                       labstr(i) = 'pertgov'
                       para(i) = pertgov
                   case('b2gratio') ! 68
                       i = i + 1
                       read( value_string,*) b2gratio
                       labstr(i) = 'b2gratio'
                       para(i) = b2gratio
                   case('btax') ! 69
                       i = i + 1
                       read( value_string,*) btax
                       labstr(i) = 'btax'
                       para(i) = btax
                   case('ptax') ! 70
                       i = i + 1
                       read( value_string,*) ptax
                       labstr(i) = 'ptax'
                       para(i) = ptax
                   case('stax') ! 71
                       i = i + 1
                       read( value_string,*) stax
                       labstr(i) = 'stax'
                       para(i) = stax     
                   case('nmul') ! 72
                       i = i + 1
                       read( value_string,*) nmul
                       labstr(i) = 'nmul'
                       para(i) = nmul   
                   case('ndim') ! 73
                       i = i + 1
                       read( value_string,*) ndim
                       labstr(i) = 'ndim'
                       para(i) = ndim    
                   case('tauk') ! 74
                       i = i + 1
                       read( value_string,*) tauk
                       labstr(i) = 'tauk'
                       para(i) = tauk    
                   case('fnhdim') ! 75
                       i = i + 1
                       read( value_string,*) fnhdim
                       labstr(i) = 'fnhdim'
                       para(i) = fnhdim    
                   case('sdim') ! 76
                       i = i + 1
                       read( value_string,*) sdim
                       labstr(i) = 'sdim'
                       para(i) = sdim     
                   case('err_svdist') ! 77
                       i = i + 1
                       read( value_string,*) err_svdist
                       labstr(i) = 'err_svdist'
                       para(i) = err_svdist  
                   case('iterasvmax') ! 78
                       i = i + 1
                       read( value_string,*) iterasvmax
                       labstr(i) = 'iterasvmax'
                       para(i) = iterasvmax   
                   case('nsdim') ! 79
                       i = i + 1
                       read( value_string,*) nsdim
                       labstr(i) = 'nsdim'
                       para(i) = nsdim ! iterabrent    
                   case('gsdim') ! 80
                       i = i + 1
                       read( value_string,*) gsdim
                       labstr(i) = 'gsdim'
                       para(i) = gsdim ! iterabrent 
                   case('trsld') ! 81
                       i = i + 1
                       read( value_string,*) trsld
                       labstr(i) = 'trsld'
                       para(i) = trsld     
                   case('iterainvdist') ! 82
                       i = i + 1
                       read( value_string,*) iterainvdist
                       labstr(i) = 'iterainvdist'
                       para(i) = iterainvdist                        
                   case('tausv') ! 83
                       i = i + 1
                       read( value_string,*) tausv
                       labstr(i) = 'tausv'
                       para(i) = tausv
                   case('taubp') ! 84
                       i = i + 1
                       read( value_string,*) taubp
                       labstr(i) = 'taubp'
                       para(i) = taubp   
                   case('target1') ! 85
                       i = i + 1
                       read( value_string,*) targetv(1) 
                       labstr(i) = 'capital in corporate'                         
                       para(i) = targetv(1)
                   case('target2') ! 86
                       i = i + 1
                       read( value_string,*) targetv(2) 
                       labstr(i) = 'small biz projects  '
                       para(i) = targetv(2)
                   case('target3') ! 87
                       i = i + 1
                       read( value_string,*) targetv(3) 
                       labstr(i) = 'median biz projects '                                                
                       para(i) = targetv(3)
                   case('target4') ! 88
                       i = i + 1
                       read( value_string,*) targetv(4) 
                       labstr(i) = 'large biz projects  '                                                
                       para(i) = targetv(4)
                   case('target5') ! 89
                       i = i + 1
                       read( value_string,*) targetv(5) 
                       labstr(i) = 'bizmen size         '                                                
                       para(i) = targetv(5)
                   case('target6') ! 90
                       i = i + 1
                       read( value_string,*) targetv(6) 
                       labstr(i) = 'financial capital   '                                                
                       para(i) = targetv(6)
                   case('target7') ! 91
                       i = i + 1
                       read( value_string,*) targetv(7) 
                       labstr(i) = 'housing capital     '
                       para(i) = targetv(7)
                   case('target8') ! 92
                       i = i + 1
                       read( value_string,*) targetv(8) 
                       labstr(i) = 'bizmen capital      '
                       para(i) = targetv(8)
                   case('target9') ! 93
                       i = i + 1
                       read( value_string,*) targetv(9) 
                       labstr(i) = 'bizmen income       '
                       para(i) = targetv(9)
                   case('target10') ! 94
                       i = i + 1
                       read( value_string,*) targetv(10) 
                       labstr(i) = 'ratio of med assets '                       
                       para(i) = targetv(10)
                   case('lower1') ! 95
                       i = i + 1
                       read( value_string,*) range_guess(1,1) 
                       labstr(i) = 'lower size of the smallest project'                         
                       para(i) = range_guess(1,1)
                   case('upper1') ! 96
                       i = i + 1
                       read( value_string,*) range_guess(1,2) 
                       labstr(i) = 'upper size of the smallest project'                         
                       para(i) = range_guess(1,2)                       
                   case('lower2') ! 97
                       i = i + 1
                       read( value_string,*) range_guess(2,1) 
                       labstr(i) = 'lower prtk0 prob of new ideas'                         
                       para(i) = range_guess(2,1)
                   case('upper2') ! 98
                       i = i + 1
                       read( value_string,*) range_guess(2,2) 
                       labstr(i) = 'upper prtk0 prob of new ideas'                         
                       para(i) = range_guess(2,2)       
                   case('lower3') ! 99
                       i = i + 1
                       read( value_string,*) range_guess(3,1) 
                       labstr(i) = 'lower prtk1 prob of new ideas'                         
                       para(i) = range_guess(3,1)
                   case('upper3') ! 100
                       i = i + 1
                       read( value_string,*) range_guess(3,2) 
                       labstr(i) = 'upper prtk1 prob of new ideas'                         
                       para(i) = range_guess(3,2)   
                   case('lower4') ! 101
                       i = i + 1
                       read( value_string,*) range_guess(4,1) 
                       labstr(i) = 'lower prtk2 prob of new ideas'                         
                       para(i) = range_guess(4,1)
                   case('upper4') ! 102
                       i = i + 1
                       read( value_string,*) range_guess(4,2) 
                       labstr(i) = 'upper prtk2 prob of new ideas'                         
                       para(i) = range_guess(4,2)   
                   case('lower5') ! 103
                       i = i + 1
                       read( value_string,*) range_guess(5,1) 
                       labstr(i) = 'lower average tech shock'                         
                       para(i) = range_guess(5,1)
                   case('upper5') ! 104
                       i = i + 1
                       read( value_string,*) range_guess(5,2) 
                       labstr(i) = 'upper average tech shock'                         
                       para(i) = range_guess(5,2)   
                   case('lower6') ! 105
                       i = i + 1
                       read( value_string,*) range_guess(6,1) 
                       labstr(i) = 'lower discount factor'                         
                       para(i) = range_guess(6,1)
                   case('upper6') ! 106
                       i = i + 1
                       read( value_string,*) range_guess(6,2) 
                       labstr(i) = 'upper discount factor'                         
                       para(i) = range_guess(6,2)   
                   case('lower7') ! 107
                       i = i + 1
                       read( value_string,*) range_guess(7,1) 
                       labstr(i) = 'lower housing utility parameter'                         
                       para(i) = range_guess(7,1)
                   case('upper7') ! 108
                       i = i + 1
                       read( value_string,*) range_guess(7,2) 
                       labstr(i) = 'upper housing utility parameter'                         
                       para(i) = range_guess(7,2)   
                   case('lower8') ! 109
                       i = i + 1
                       read( value_string,*) range_guess(8,1) 
                       labstr(i) = 'lower phi1'                         
                       para(i) = range_guess(8,1)
                   case('upper8') ! 110
                       i = i + 1
                       read( value_string,*) range_guess(8,2) 
                       labstr(i) = 'upper phi1'                         
                       para(i) = range_guess(8,2)   
                   case('lower9') ! 111
                       i = i + 1
                       read( value_string,*) range_guess(9,1) 
                       labstr(i) = 'lower phi2'                         
                       para(i) = range_guess(9,1)
                   case('upper9') ! 112
                       i = i + 1
                       read( value_string,*) range_guess(9,2) 
                       labstr(i) = 'upper phi2'                         
                       para(i) = range_guess(9,2)   
                   case('lower10') ! 113
                       i = i + 1
                       read( value_string,*) range_guess(10,1) 
                       labstr(i) = 'lower phi3'                         
                       para(i) = range_guess(10,1)
                   case('upper10') ! 114
                       i = i + 1
                       read( value_string,*) range_guess(10,2) 
                       labstr(i) = 'upper phi3'                         
                       para(i) = range_guess(10,2) 

                   case('trylen') ! 115
                       i = i + 1
                       read( value_string,*) trylen 
                       labstr(i) = 'trylen'                         
                       para(i) = trylen    
                   case('sblno1') ! 116
                       i = i + 1
                       read( value_string,*) sblno1 
                       labstr(i) = 'sblno1'                         
                       para(i) = sblno1 
                   case('nrow') ! 117
                       i = i + 1
                       read( value_string,*) nrow 
                       labstr(i) = 'nrow'                         
                       para(i) = nrow 
                   case('stepsize') ! 118
                       i = i + 1
                       read( value_string,*) stepsize 
                       labstr(i) = 'stepsize'                         
                       para(i) = stepsize 
                   case('reinifac') ! 119
                       i = i + 1
                       read( value_string,*) reinifac 
                       labstr(i) = 'reinifac'                         
                       para(i) = reinifac 
                   case('amoalp') ! 120
                       i = i + 1
                       read( value_string,*) amoalp 
                       labstr(i) = 'amoalp'                         
                       para(i) = amoalp 
                   case('amogam') ! 121
                       i = i + 1
                       read( value_string,*) amogam 
                       labstr(i) = 'amogam'                         
                       para(i) = amogam 
                   case('amobet') ! 122
                       i = i + 1
                       read( value_string,*) amobet 
                       labstr(i) = 'amobet'                         
                       para(i) = amobet 
                   case('amotau') ! 123
                       i = i + 1
                       read( value_string,*) amotau 
                       labstr(i) = 'amotau'                         
                       para(i) = amotau 
                   case('amoerrcrt') ! 124
                       i = i + 1
                       read( value_string,*) amoerrcrt 
                       labstr(i) = 'amoerrcrt'                         
                       para(i) = amoerrcrt 
                   case('amoitrcrt') ! 125
                       i = i + 1
                       read( value_string,*) amoitrcrt  
                       labstr(i) = 'amoitrcrt'                         
                       para(i) = amoitrcrt 
                   case('amoconvex') ! 126
                       i = i + 1
                       read( value_string,*) amoconvex 
                       labstr(i) = 'amoconvex'                         
                       para(i) = amoconvex 
                   case('amoiniswitch') ! 127
                       i = i + 1
                       read( value_string,*) amoiniswitch  
                       labstr(i) = 'amoiniswitch'                         
                       para(i) = amoiniswitch  
                   case('nsbq') ! 128
                       i = i + 1
                       read( value_string,*) nsbq  
                       labstr(i) = 'nsbq'                         
                       para(i) = nsbq       
                   case('slist') ! 129
                       i = i + 1
                       read( value_string,*) slist
                       labstr(i) = 'slist'
                       para(i) = slist
                   case('elist') ! 130
                       i = i + 1
                       read( value_string,*) elist
                       labstr(i) = 'elist'
                       para(i) = elist
                   case('noamoeba') ! 131
                       i = i + 1
                       read( value_string,*) noamoeba
                       labstr(i) = 'noamoeba'
                       para(i) = noamoeba
                   case('maxdist') ! 132
                       i = i + 1
                       read( value_string,*) maxdist
                       labstr(i) = 'maxdist'
                       para(i) = maxdist                    
                   case('tinymass') ! 133
                       i = i + 1
                       read( value_string,*) tinymass
                       labstr(i) = 'tinymass'
                       para(i) = tinymass  
                   case('chunkmass') ! 134
                       i = i + 1
                       read( value_string,*) chunkmass
                       labstr(i) = 'chunkmass'
                       para(i) = chunkmass       
                   case('initau') ! 135
                       i = i + 1
                       read( value_string,*) initau
                       labstr(i) = 'initau'
                       para(i) = initau
                   case('testfunc_idx') ! 136
                       i = i + 1
                       read( value_string,*) testfunc_idx
                       labstr(i) = 'testfunc_idx'
                       para(i) = testfunc_idx          
                   case('accumass') ! 137
                       i = i + 1
                       read( value_string,*) accumass
                       labstr(i) = 'accumass'
                       para(i) = accumass  
                   case('accupara') ! 138
                       i = i + 1
                       read( value_string,*) accupara
                       labstr(i) = 'accupara'
                       para(i) = accupara
                   case('errvalfc') ! 139
                       i = i + 1
                       read( value_string,*) errvalfc
                       labstr(i) = 'errvalfc'
                       para(i) = errvalfc
                   case('iota') ! 140
                       i = i + 1
                       read( value_string,*) iota
                       labstr(i) = 'iota'
                       para(i) = iota  
                   case('momround') ! 141
                       i = i + 1
                       read( value_string,*) momround
                       labstr(i) = 'momround'
                       para(i) = momround   
                   case('setindex') ! 142
                       i = i + 1
                       read( value_string,*) setindex
                       labstr(i) = 'setindex'
                       para(i) = setindex 
                   case('listlength') ! 143
                       i = i + 1
                       read( value_string,*) listlength
                       labstr(i) = 'listlength'
                       para(i) = listlength      
                   case('subdim') ! 144
                       i = i + 1
                       read( value_string,*) subdim
                       labstr(i) = 'subdim'
                       para(i) = subdim   
                   case('obj_func_toggle') ! 145
                       i = i + 1
                       read( value_string,*) obj_func_toggle
                       labstr(i) = 'obj_func_toggle'
                       para(i) = obj_func_toggle                        
                end select
            enddo
        else
            write(*,'(a)') 'Failed to read ''para1''.'
        endif
    end subroutine read_parameter_model              
    
    subroutine set_markov_matrix( rhoy, vary, rhoyh, varyh, py, y, sy, pyh, yh, syh )
        implicit none
        real(wp), intent(in) :: rhoy, vary, rhoyh, varyh
        real(wp), dimension(:,:), intent(out) :: py, pyh
        real(wp), dimension(:), intent(out) :: y, yh, sy, syh
        integer :: n
        
        n = size(y)
        allocate( p_tau(n,n), y_tau(n), s_tau(n) )
        call tauchen( rhoy, vary**(1/2._wp), py, y, sy ) 
        deallocate( p_tau, y_tau, s_tau )
        
        n = size(yh)
        allocate( p_tau(n,n), y_tau(n), s_tau(n) )
        ! the next line makes the spread of inheritance AR(1) identical with that of labor efficiency AR(1)
        call tauchen( rhoyh, varyh**(1/2._wp), pyh, yh, syh, smin_tau*vary**(1/2._wp)/varyh**(1/2._wp) ) 
        ! use the optimal smin_tau from the AR(1) "labor" process
        ! Variable 'smin_tau' is a global variable defined in module TOOLBOXS, which passes the optimal spread parameter for 
        ! yh process, based on the information obtained in the host subroutine 'tauchen' using brent algorithm.
        ! Here is the idea: I requirs both the discretized Markov chain to have identical "level of boundary." In our case, 
        deallocate( p_tau, y_tau, s_tau )
        
        y = exp(y) 
        yh = exp(yh)
        
        !call ss(y,'y')
        !call ss(yh,'yh')
        !call sm(py,'py')
        !call sm(pyh,'pyh')
        !call ss(sy,'sy')
        !call ss(syh,'syh')  
    end subroutine set_markov_matrix        
    
    subroutine set_survival_vector( survprob, input_file ) ! V021515 
        ! conditional survival probability.
        ! examplar for array reading.
        implicit none
        real(wp), dimension(:), intent(out) :: survprob
        
        character(len=*), intent(in) :: input_file
        character(len=80) :: value_string 
        
        integer :: iostat, i
        integer :: n
        n = size(survprob)
        survprob = 1._wp
        open( unit=10, file=input_file, status='old', action='read', iostat=iostat )
        if( iostat==0 ) then
            i = 1
            do
                read( unit=10, fmt=*, iostat=iostat ) value_string
                if( iostat/=0 ) exit
                read( value_string, fmt=* ) survprob(i)
                i = i + 1
                if(i>n) exit
            enddo
        else
            write(*,'(a)') '[error] Unable to read _survival_probability.txt.'
        endif    
    end subroutine set_survival_vector     
    
    subroutine set_efficiency_untis( efflab, input_file ) ! V021515 
        implicit none
        real(wp), dimension(:), intent(out) :: efflab
        
        character(len=*), intent(in) :: input_file
        character(len=80) :: value_string 
        
        integer :: iostat, i
        integer :: n
        n = size(efflab)
        efflab = 1._wp
        open( unit=10, file=input_file, status='old', action='read', iostat=iostat )
        if( iostat==0 ) then
            i = 1
            do
                read( unit=10, fmt=*, iostat=iostat ) value_string
                if( iostat/=0 ) exit
                read( value_string, fmt=* ) efflab(i)
                i = i + 1
                if(i>n) exit
            enddo
        else
            write(*,'(a)') '[error] Unable to read _efficiency_untis.txt.'
        endif    
    end subroutine set_efficiency_untis   
    
    subroutine set_grids_for_housing_asset(hray,s)
        implicit none
        real(wp), intent(out) :: hray(:)
        integer, intent(in) :: s
        integer :: i
        do i = 1, hdim
            hray(i) = (hmin**(1._wp/s) + (hmax**(1._wp/s)-hmin**(1._wp/s))*(i-1._wp)/(hdim-1._wp))**s    
        enddo
    end subroutine 
    
    subroutine set_grids_for_nonhousing_asset(av,s,msg)
        implicit none
        real(wp), intent(out) :: av(:)
        real(wp), intent(in) :: s
        integer :: i
        character(len=*), intent(in) :: msg
        select case( msg )
            case('coarse')
                do i = 1, adim
                    av(i) = ((i-1._wp)/(adim-1._wp))**s*(amax-amin) + amin
                enddo
            case('refined')
                do i = 1, fnadim
                    av(i) = ((i-1._wp)/(fnadim-1._wp))**s*(amax-amin) + amin
                enddo
        endselect
    end subroutine set_grids_for_nonhousing_asset
    
    !subroutine evenly_distributed_since(av,n)
    !    implicit none
    !    real(wp), intent(inout) :: av(:)
    !    integer, intent(in) :: n
    !    real(wp), dimension(:), allocatable :: sv
    !    integer :: m
    !    real(wp) :: first,last
    !    m = size(av)
    !    first = av(n)
    !    last  = av(m)
    !    allocate( sv(m-n+1) )
    !    call grid(sv,first,last,1._wp)
    !    av(n:m) = sv
    !    deallocate( sv )
    !end subroutine evenly_distributed_since
    
    !subroutine asset_grid_renew(newav,upper,lower,s,msg)
    !    implicit none
    !    real(wp), intent(out) :: newav(:)
    !    real(wp), intent(in) :: s,upper,lower
    !    real(wp) :: amax, amin
    !    integer :: i 
    !    character(len=*), intent(in) :: msg
    !    amax = upper
    !    amin = lower      
    !    if(amax==amin)then
    !        newav = penalty
    !        newav(fnadim) = amax
    !    else
    !        select case( msg )
    !        case('refined')
    !            do i = 1, fnadim
    !                newav(i) = ((i-1._wp)/(fnadim-1._wp))**s*(amax-amin) + amin
    !            enddo
    !        endselect
    !    endif
    !end subroutine asset_grid_renew
    
    subroutine create_index_list(xlist,bdmat,blist,msg1,msg2)  ! 072416
        implicit none
        integer, dimension(:,:), intent(out) :: xlist, blist, bdmat
        integer, dimension(:,:), allocatable :: tarray
        !integer, dimension(:), allocatable :: tvector
        character(len=*), intent(in) :: msg1, msg2
        integer :: stp, stp_b, ax, hx, kx, zx, yx, kpx, ypx, opx, l, m, n, i, j, assdim
        
        if(msg2=='coarse')then
            assdim = adim
        elseif(msg2=='refined')then
            assdim = fnadim
        endif
        
        if(msg2=='coarse')then ! done 09102016
            select case( msg1 )
            case('14') 
            xlist = 0 !
            bdmat = 1 !
            stp = 1   ! list index
            stp_b = 1 ! boundary index
            
            allocate( tarray(7,6) ) ! no shadow options in t=14
            tarray(:,1) = [0,1,2,3,1,2,3] ! kx
            tarray(:,2) = [0,0,0,0,1,1,1] ! zx
            tarray(:,3) = 0 ! yx
            tarray(:,4) = 0 ! kpx
            tarray(:,5) = 0 ! ypx
            tarray(:,6) = 0 ! opx
            
            do l = 1, 7
                kx  = tarray(l,1) 
                zx  = tarray(l,2) 
                yx  = tarray(l,3) 
                kpx = tarray(l,4) 
                ypx = tarray(l,5) 
                opx = tarray(l,6)  
                do m = 1, 1 ! hdim !<-----------------------------------
                    
                    ! Used for value function of the beginning of period ??
                    blist(stp_b,1) = m
                    blist(stp_b,2) = kx
                    blist(stp_b,3) = zx
                    blist(stp_b,4) = yx
                    stp_b = stp_b + 1
                    
                    ! Used for value function of the end of period
                    
                    do n = 1, assdim
                        xlist(stp,1) = n 
                        xlist(stp,2) = m
                        xlist(stp,3) = kx
                        xlist(stp,4) = zx
                        xlist(stp,5) = yx
                        xlist(stp,6) = kpx
                        xlist(stp,7) = ypx
                        xlist(stp,8) = opx
                        stp = stp + 1        
                    enddo
                    
                enddo
            enddo
            bdmat(1,2) = stp - 1 ! the ending index ! only one case: opx == 0 in period 14.
            deallocate( tarray )
            
            case('10-13')
            xlist = 0
            bdmat = 1
            stp   = 1
            stp_b = 1  
            
            allocate( tarray(12,6) ) !# 1 including shadwo options t = 10-13
            tarray(:,1) = [0,1,2,3,1,2,3,1,2,3,1,2] ! kx #2
            tarray(:,2) = [0,0,0,0,1,1,1,1,1,1,1,1] ! zx
            tarray(:,3) = 0 ! yx
            tarray(:,4) = [0,0,0,0,0,0,0,1,2,3,2,3] ! kp
            tarray(:,5) = 0 ! yp
            tarray(:,6) = [0,0,0,0,0,0,0,1,1,1,2,2] ! op
            
            do l = 1, 12 !#3
                kx  = tarray(l,1) 
                zx  = tarray(l,2) 
                yx  = tarray(l,3) 
                kpx = tarray(l,4) 
                ypx = tarray(l,5) 
                opx = tarray(l,6)      
                do m = 1, 1 ! hdim
                    do n = 1, assdim
                        ! boundary matrix
                        if(l==1)then
                            if(m==1.and.n==1)then
                                bdmat(stp_b,1) = stp 
                            ! elseif(m==hdim.and.n==assdim)then
                            elseif(m==1.and.n==assdim)then
                                bdmat(stp_b,2) = stp 
                                stp_b = stp_b + 1
                            endif
                        else
                            if(l==2.or.l==8.or.l==11)then !#4
                                if(m==1.and.n==1)then
                                    bdmat(stp_b,1) = stp 
                                endif
                            elseif(l==7.or.l==10.or.l==12)then !#5
                                !if(m==hdim.and.n==assdim)then
                                if(m==1.and.n==assdim)then
                                    bdmat(stp_b,2) = stp
                                    stp_b = stp_b + 1
                                endif         
                            endif
                        endif
                        
                        xlist(stp,1) = n
                        xlist(stp,2) = m
                        xlist(stp,3) = kx
                        xlist(stp,4) = zx
                        xlist(stp,5) = yx
                        xlist(stp,6) = kpx
                        xlist(stp,7) = ypx
                        xlist(stp,8) = opx
                        stp = stp + 1        
                    enddo ! n
                enddo
            enddo
            deallocate( tarray )
            
            case('9') 
            xlist = 0
            bdmat = 1        
            stp   = 1
            stp_b = 1
            
            allocate( tarray(16,6) )
            tarray(:,1) = [0,1,2,3,1,2,3,1,2,3,0,1,2,1,2,3] ! k
            tarray(:,2) = [0,0,0,0,1,1,1,1,1,1,0,1,1,0,0,0] ! z
            tarray(:,3) = 1 ! y
            tarray(:,4) = [0,0,0,0,0,0,0,1,2,3,1,2,3,1,1,1] ! kp
            tarray(:,5) = 0 ! yp
            tarray(:,6) = [0,0,0,0,0,0,0,1,1,1,2,2,2,2,2,2] ! op
            
            do l = 1, 16
                do i = 1, nmc ! current period labor efficiency
                    kx  = tarray(l,1)
                    zx  = tarray(l,2)
                    yx  = tarray(l,3) ! set to 1
                    kpx = tarray(l,4)
                    ypx = tarray(l,5) ! set to 0
                    opx = tarray(l,6)
                    yx  = yx*i ! won't be inflated. becauase yx is renewed already.
                    
                    do m = 1, 1 ! hdim
                        do n = 1, assdim
                            if(l==1)then
                                if(i==1.and.m==1.and.n==1)then
                                    bdmat(stp_b,1) = stp 
                                elseif(i==nmc.and.m==1.and.n==assdim)then
                                    bdmat(stp_b,2) = stp 
                                    stp_b = stp_b + 1
                                endif                            
                            else ! l = 2~16
                                if(i==1.and.(l==2.or.l==8.or.l==11))then
                                    if(m==1.and.n==1)then
                                        bdmat(stp_b,1) = stp
                                    endif
                                elseif(i==nmc.and.(l==7.or.l==10.or.l==16))then ! 4.8.2017 replace l=13 with l=16.
                                    if(m==1.and.n==assdim)then
                                        bdmat(stp_b,2) = stp
                                        stp_b = stp_b + 1
                                    endif
                                endif
                            endif
                            xlist(stp,1) = n
                            xlist(stp,2) = m
                            xlist(stp,3) = kx
                            xlist(stp,4) = zx
                            xlist(stp,5) = yx
                            xlist(stp,6) = kpx
                            xlist(stp,7) = ypx
                            xlist(stp,8) = opx
                            stp = stp + 1
                        enddo ! n
                    enddo ! m
                enddo ! i
            enddo ! l
            deallocate( tarray )
            
            case('18')
            xlist = 0
            bdmat = 1        
            stp   = 1
            stp_b = 1
            !allocate( tarray(13,6) )
            allocate( tarray(16,6) ) ! 3.30.2017
            tarray(:,1) = [0,1,2,3,1,2,3,1,2,3,0,1,2,1,2,3] ! k
            tarray(:,2) = [0,0,0,0,1,1,1,1,1,1,0,1,1,0,0,0] ! z
            tarray(:,3) = 1                           ! y
            tarray(:,4) = [0,0,0,0,0,0,0,1,2,3,1,2,3,1,1,1] ! kp
            tarray(:,5) = 1                           ! yp
            tarray(:,6) = [0,0,0,0,0,0,0,1,1,1,2,2,2,2,2,2] ! op
            
            do l = 1, 16
                do i = 1, nmc ! labor now
                    do j = 1, nmc ! labor tomorrow
                        kx  = tarray(l,1)
                        zx  = tarray(l,2)
                        yx  = tarray(l,3) ! yx
                        kpx = tarray(l,4)
                        ypx = tarray(l,5) ! ypx
                        opx = tarray(l,6)
                        yx  = yx*i
                        ypx = ypx*j
                        do m = 1, 1 ! hdim. 3.30.2017 Correct set to work with 2 dimensional search in model.f90. 
                            do n = 1, assdim
                                if(l==1)then
                                    if(i==1.and.j==1.and.m==1.and.n==1)then
                                        bdmat(stp_b,1) = stp 
                                    elseif(i==nmc.and.j==nmc.and.m==1.and.n==assdim)then
                                        bdmat(stp_b,2) = stp 
                                        stp_b = stp_b + 1
                                    endif                                  
                                else ! l = 2~16
                                    if((i==1.and.j==1).and.(l==2.or.l==8.or.l==11))then
                                        if(m==1.and.n==1)then
                                            bdmat(stp_b,1) = stp    
                                        endif
                                    elseif((i==nmc.and.j==nmc).and.(l==7.or.l==10.or.l==16))then ! 4.8.2017 replace l=13 with l=16.
                                        if(m==1.and.n==assdim)then
                                            bdmat(stp_b,2) = stp
                                            stp_b = stp_b + 1
                                        endif
                                    endif    
                                endif
                                xlist(stp,1) = n
                                xlist(stp,2) = m
                                xlist(stp,3) = kx
                                xlist(stp,4) = zx
                                xlist(stp,5) = yx
                                xlist(stp,6) = kpx
                                xlist(stp,7) = ypx
                                xlist(stp,8) = opx
                                stp = stp + 1
                            enddo ! n
                        enddo ! m
                    enddo ! j
                enddo ! i
            enddo ! l
            deallocate( tarray )
            
            end select
            
        else ! msg2 == 'refined' 09102016 not yet
            select case( msg1 )
            case('14') ! last period
            xlist = 0
            bdmat = 1
            stp = 1
            stp_b = 1
            
            allocate( tarray(7,6) )
            tarray(:,1) = [0,1,2,3,1,2,3] ! kx
            tarray(:,2) = [0,0,0,0,1,1,1] ! zx
            tarray(:,3) = 0 ! yx
            
            tarray(:,4) = 0 ! kpx
            tarray(:,5) = 0 ! ypx
            tarray(:,6) = 0 ! opx
            
            do l = 1, 7
                kx  = tarray(l,1) 
                zx  = tarray(l,2) 
                yx  = tarray(l,3) 
                kpx = tarray(l,4) 
                ypx = tarray(l,5) 
                opx = tarray(l,6)  
                do m = 1, hdim
                    ! Used for value function of the beginning of period
                    blist(stp_b,1) = m
                    blist(stp_b,2) = kx
                    blist(stp_b,3) = zx
                    blist(stp_b,4) = yx
                    stp_b = stp_b + 1
                    ! Used for value function of the end of period
                    do n = 1, assdim
                        xlist(stp,1) = n 
                        xlist(stp,2) = m
                        xlist(stp,3) = kx
                        xlist(stp,4) = zx
                        xlist(stp,5) = yx
                        xlist(stp,6) = kpx
                        xlist(stp,7) = ypx
                        xlist(stp,8) = opx
                        stp = stp + 1        
                    enddo
                enddo
            enddo
            bdmat(1,2) = stp - 1 ! the ending index ! only one case: opx == 0 in period 14.
            deallocate( tarray )
            
            case('10-13')
            xlist = 0
            bdmat = 1
            stp   = 1
            stp_b = 1  
            
            allocate( tarray(9,6) ) !#1
            tarray(:,1) = [0,1,2,3,1,2,3,1,2] ! kx #2
            tarray(:,2) = [0,0,0,0,1,1,1,1,1] ! zx
            tarray(:,3) = 0 ! yx
            tarray(:,4) = [0,0,0,0,1,2,3,2,3] ! kp
            tarray(:,5) = 0 ! yp
            tarray(:,6) = [0,0,0,0,1,1,1,2,2] ! op
            
            do l = 1, 9 !#3
                kx  = tarray(l,1) 
                zx  = tarray(l,2) 
                yx  = tarray(l,3) 
                kpx = tarray(l,4) 
                ypx = tarray(l,5) 
                opx = tarray(l,6)      
                
                do m = 1, hdim
                    do n = 1, assdim
                        ! boundary matrix
                        if(l==1)then
                            if(m==1.and.n==1)then
                                bdmat(stp_b,1) = stp 
                            elseif(m==hdim.and.n==assdim)then
                                bdmat(stp_b,2) = stp 
                                stp_b = stp_b + 1
                            endif
                        else
                            if(l==2.or.l==5.or.l==8)then !#4
                                if(m==1.and.n==1)then
                                    bdmat(stp_b,1) = stp 
                                endif
                            elseif(l==4.or.l==7.or.l==9)then !#5
                                if(m==hdim.and.n==assdim)then
                                    bdmat(stp_b,2) = stp
                                    stp_b = stp_b + 1
                                endif         
                            endif
                        endif
                        
                        xlist(stp,1) = n
                        xlist(stp,2) = m
                        xlist(stp,3) = kx
                        xlist(stp,4) = zx
                        xlist(stp,5) = yx
                        xlist(stp,6) = kpx
                        xlist(stp,7) = ypx
                        xlist(stp,8) = opx
                        stp = stp + 1        
                    enddo ! n
                enddo
            enddo
            deallocate( tarray )
            
            case('9') 
            xlist = 0
            bdmat = 1        
            stp   = 1
            stp_b = 1
            
            allocate( tarray(10,6) )
            tarray(:,1) = [0,1,2,3,1,2,3,0,1,2]
            tarray(:,2) = [0,0,0,0,1,1,1,0,1,1]
            tarray(:,3) = 1
            tarray(:,4) = [0,0,0,0,1,2,3,1,2,3]
            tarray(:,5) = 0
            tarray(:,6) = [0,0,0,0,1,1,1,2,2,2] 
            
            do l = 1, 10
                do i = 1, nmc
                    kx  = tarray(l,1)
                    zx  = tarray(l,2)
                    yx  = tarray(l,3) ! set to 1
                    kpx = tarray(l,4)
                    ypx = tarray(l,5) ! set to 0
                    opx = tarray(l,6)
                    yx  = yx*i
                    
                    do m = 1, hdim
                        do n = 1, assdim
                            if(l==1)then
                                if(i==1.and.m==1.and.n==1)then
                                    bdmat(stp_b,1) = stp 
                                elseif(i==nmc.and.m==hdim.and.n==assdim)then
                                    bdmat(stp_b,2) = stp 
                                    stp_b = stp_b + 1
                                endif                            
                            else
                                if(i==1.and.(l==2.or.l==5.or.l==8))then
                                    if(m==1.and.n==1)then
                                        bdmat(stp_b,1) = stp
                                    endif
                                elseif(i==nmc.and.(l==4.or.l==7.or.l==10))then
                                    if(m==hdim.and.n==assdim)then
                                        bdmat(stp_b,2) = stp
                                        stp_b = stp_b + 1
                                    endif
                                endif
                            endif
                            xlist(stp,1) = n
                            xlist(stp,2) = m
                            xlist(stp,3) = kx
                            xlist(stp,4) = zx
                            xlist(stp,5) = yx
                            xlist(stp,6) = kpx
                            xlist(stp,7) = ypx
                            xlist(stp,8) = opx
                            stp = stp + 1
                        enddo ! n
                    enddo ! m
                enddo ! i
            enddo ! l
            deallocate( tarray )
            
            case('18')
            xlist = 0
            bdmat = 1        
            stp   = 1
            stp_b = 1
            allocate( tarray(10,6) )
            tarray(:,1) = [0,1,2,3,1,2,3,0,1,2] ! k
            tarray(:,2) = [0,0,0,0,1,1,1,0,1,1] ! z
            tarray(:,3) = 1                     ! y
            tarray(:,4) = [0,0,0,0,1,2,3,1,2,3] ! kp
            tarray(:,5) = 1                     ! yp
            tarray(:,6) = [0,0,0,0,1,1,1,2,2,2] ! op
            
            do l = 1, 10
                do i = 1, nmc ! labor now
                    do j = 1, nmc ! labor tomorrow
                        kx  = tarray(l,1)
                        zx  = tarray(l,2)
                        yx  = tarray(l,3) ! yx
                        kpx = tarray(l,4)
                        ypx = tarray(l,5) ! ypx
                        opx = tarray(l,6)
                        yx  = yx*i
                        ypx = ypx*j
                        do m = 1, hdim
                            do n = 1, assdim
                                if(l==1)then
                                    if(i==1.and.j==1.and.m==1.and.n==1)then
                                        bdmat(stp_b,1) = stp 
                                    elseif(i==nmc.and.j==nmc.and.m==hdim.and.n==assdim)then
                                        bdmat(stp_b,2) = stp 
                                        stp_b = stp_b + 1
                                    endif                                  
                                else
                                    if((i==1.and.j==1).and.(l==2.or.l==5.or.l==8))then
                                        if(m==1.and.n==1)then
                                            bdmat(stp_b,1) = stp    
                                        endif
                                    elseif((i==nmc.and.j==nmc).and.(l==4.or.l==7.or.l==10))then
                                        if(m==hdim.and.n==assdim)then
                                            bdmat(stp_b,2) = stp
                                            stp_b = stp_b + 1
                                        endif
                                    endif    
                                endif
                                xlist(stp,1) = n
                                xlist(stp,2) = m
                                xlist(stp,3) = kx
                                xlist(stp,4) = zx
                                xlist(stp,5) = yx
                                xlist(stp,6) = kpx
                                xlist(stp,7) = ypx
                                xlist(stp,8) = opx
                                stp = stp + 1
                            enddo ! n
                        enddo ! m
                    enddo ! j
                enddo ! i
            enddo ! l
            deallocate( tarray )
            
            end select            
        endif ! msg2
        
    end subroutine create_index_list   
    
    ! 3.31.2017 Here contains all the end-of-period combinations (excluding the cases where entrepreneurs hit by "good" business shock 
    subroutine serialindices_Map2_coordinates(ray,mat,sza,szh) 
        implicit none
        integer, dimension(:,:), intent(out), allocatable :: ray ! FOR OBTAINING CORRESPONDING STATE VARIABLES TO A SPECIFIC SERIAL NUMBER.
        integer, intent(in) :: sza, szh
        !integer, dimension(1:fnadim,1:fnhdim,0:kdim-1,0:1,0:nmc,0:kdim-1,0:nmc,0:2,1:14), intent(out) :: mat ! FOR OBTAINING SERIAL NUMBER.
        integer, dimension(1:sza,1:szh,0:kdim-1,0:1,0:nmc,0:kdim-1,0:nmc,0:2,1:14), intent(out) :: mat ! FOR OBTAINING SERIAL NUMBER.
        
        integer, dimension(:,:), allocatable :: vec
        integer :: idx, t, i, j, l, m, n, mi, ni, bnd, dx2
        
        mat = 0
        !allocate(ray(fnadim*fnhdim*1018,10))
        allocate(ray(sza*szh*1018,10)) ! The total number of combination is counted in variable_space_v2.xlsx.
        idx = 1
        dx2 = 1
        do t = 1, 14
            if(1<=t.and.t<=8)then
                allocate(vec(4,13))
                vec(1,:) = [0,1,2,3,1,2,3,0,1,2,1,2,3] ! k
                vec(2,:) = [0,0,0,0,1,1,1,0,1,1,0,0,0] ! z
                vec(3,:) = [0,0,0,0,1,2,3,1,2,3,1,1,1] ! kp
                vec(4,:) = [0,0,0,0,1,1,1,2,2,2,2,2,2] ! op
                bnd = 13 ! LENGTH OF THE VECTOR ABOVE
            elseif(t==9)then
                allocate(vec(4,13))
                vec(1,:) = [0,1,2,3,1,2,3,0,1,2,1,2,3] ! k
                vec(2,:) = [0,0,0,0,1,1,1,0,1,1,0,0,0] ! z
                vec(3,:) = [0,0,0,0,1,2,3,1,2,3,1,1,1] ! kp
                vec(4,:) = [0,0,0,0,1,1,1,2,2,2,2,2,2] ! op                
                bnd = 13
            elseif(10<=t.and.t<=13)then
                allocate(vec(4,9))
                vec(1,:) = [0,1,2,3,1,2,3,1,2] ! k
                vec(2,:) = [0,0,0,0,1,1,1,1,1] ! z
                vec(3,:) = [0,0,0,0,1,2,3,2,3] ! kp
                vec(4,:) = [0,0,0,0,1,1,1,2,2] ! op                
                bnd = 9
            elseif(t==14)then
                allocate(vec(4,7))
                vec(1,:) = [0,1,2,3,1,2,3] ! k
                vec(2,:) = [0,0,0,0,1,1,1] ! z
                vec(3,:) = [0,0,0,0,0,0,0] ! kp
                vec(4,:) = [0,0,0,0,0,0,0] ! op                
                bnd = 7
            endif ! t
            
            do l = 1, bnd
                do ni = 1, nmc ! yp's index: TOMORROW
                    if(t>=9)then
                        if(ni==1)then
                            n = 0 ! no labor efficiency tomorrow
                        else
                            cycle
                        endif ! ni
                    else ! t<=8
                        n = ni ! labor efficiency lottery
                    endif ! t
                    
                    do mi = 1, nmc ! y's index: TODAY
                        if(t>=10)then
                            if(mi==1)then
                                m = 0 ! for senior, today's labor efficiency is zero.
                            else
                                cycle
                            endif
                        else ! t<=9
                            m = mi 
                        endif ! t
                        
                        do j = 1, szh ! grid number of housing asset
                            do i = 1, sza ! grid number of financial asset
                                
                                ray(idx,1)  = i ! a
                                ray(idx,2)  = j ! h
                                ray(idx,3)  = vec(1,l) ! k
                                ray(idx,4)  = vec(2,l) ! z
                                ray(idx,5)  = m ! y current period's labor efficiency
                                ray(idx,6)  = vec(3,l) ! kp
                                ray(idx,7)  = n ! yp next period's labor efficiency
                                ray(idx,8)  = vec(4,l) ! op
                                ray(idx,9)  = t
                                ray(idx,10) = idx 
                                
                                mat(i,j,vec(1,l),vec(2,l),m,vec(3,l),n,vec(4,l),t) = idx
                                idx = idx + 1 
                                
                            enddo ! i
                        enddo ! j
                        
                        dx2 = dx2 + 1 ! to count the number of root combination
                        
                    enddo ! mi
                enddo ! ni
            enddo ! l
            deallocate(vec)                        
        enddo ! t
        !write(unit=120,fmt='(a,i10)'), ' The number of end period combinations: ', idx - 1
        !write(unit=120,fmt='(a,i10)'), ' The number of end period root combinations: ', dx2 - 1
    end subroutine serialindices_Map2_coordinates
      
    !subroutine set_refined_zvector(zlist,zsize,msg1)
    !    implicit none
    !    integer :: ax, hx, kx, zx, yx, i, j, l, m, n
    !    character(len=*), intent(in) :: msg1
    !    integer, intent(out) :: zsize
    !    integer, intent(out), dimension(:,:) :: zlist
    !    
    !    select case( msg1 )
    !    case('1to9')
    !    n = 1
    !    zlist = 0
    !    zsize = 0
    !    do l = 1, 7
    !        kx = tvector(l,1)
    !        zx = tvector(l,2)
    !        do yx = 1, nmc
    !            do hx = 1, hdim
    !                do ax = 1, fnadim
    !                    zlist(n,1) = ax  
    !                    zlist(n,2) = hx
    !                    zlist(n,3) = kx
    !                    zlist(n,4) = zx
    !                    zlist(n,5) = yx
    !                    n = n + 1
    !                enddo ! i
    !            enddo ! j
    !        enddo ! yx
    !    enddo ! l
    !    zsize = n - 1
    !    case('10to14')
    !    
    !    end select
    !    
    !end subroutine set_refined_zvector
    
    ! 3.15.2017 Keep it!! szt18vec1, etc. are used in subroutine valid_beginning_period_mass subroutine!!
    ! 4.1.2017 the end-of-period outcomes associated with each beginning-of-period combination is created.S
    subroutine make_next_period_index_combination() 
        implicit none
        integer :: l, apxi, opxi, kpxi, ypx, kpx, opx, t, n
        integer, dimension(:,:), allocatable :: twwmat        
        
        ! initialization of tvector GLOBALLY, although it is "not" used in the curreent subroutine
        tvector(:,1) = [0,1,2,3,1,2,3] ! current k
        tvector(:,2) = [0,0,0,0,1,1,1] ! current z
        tvector(:,3) = [0,0,0,0,0,0,0] ! current y
        
        ! initialization
        t18vec1  = -1
        t18vec24 = -1
        t18vec56 = -1
        t18vec7  = -1
        
        t9vec1  = -1 
        t9vec24 = -1
        t9vec56 = -1
        t9vec7  = -1
        
        t1013vec1  = -1
        t1013vec24 = -1
        t1013vec56 = -1
        t1013vec7  = -1
        
        t14vec1  = -1 
        t14vec24 = -1
        t14vec56 = -1
        t14vec7  = -1        
        
        ! --------------------------------------------------------- t = 1-8
        ! wage earners [t=1-8]
        l = 1
        do apxi = 1, 1 
            do opxi = 0, 1
                opx = merge(0,2,opxi==0) ! affixed with i means not directly to be used as an index
                kpx = merge(0,1,opxi==0)
                do ypx = 1, nmc ! next period's labor efficiency
                    t18vec1(l,1) = apxi ! indices of left and right brackets
                    t18vec1(l,2) = kpx
                    t18vec1(l,3) = ypx
                    t18vec1(l,4) = opx
                    l = l + 1
                enddo ! ypx 
            enddo ! opxi
        enddo ! apxi        
        szt18vec1 = l - 1
        
        ! entrepreneurs start with bad business shock [t=1-8]
        l = 1
        do apxi = 1, 1
            do opxi = 0, 1
                !opx = 0
                !kpx = 0
                opx = merge(0,2,opxi==0) ! affixed with i means not directly to be used as an index
                kpx = merge(0,1,opxi==0)            
                do ypx = 1, nmc
                    t18vec24(l,1) = apxi
                    t18vec24(l,2) = kpx
                    t18vec24(l,3) = ypx
                    t18vec24(l,4) = opx
                    l = l + 1
                enddo ! ypx
            enddo ! opxi
        enddo ! apxi
        szt18vec24 = l - 1
        
        ! small and medium scale entrepreneurs start with good business shock [t=1-8]
        l = 1
        do apxi = 1, 1
            do kpxi = 0, 1
                opx = merge(1,2,kpxi==0)
                do ypx = 1, nmc
                    t18vec56(l,1) = apxi 
                    t18vec56(l,2) = kpxi ! not index directly
                    t18vec56(l,3) = ypx
                    t18vec56(l,4) = opx
                    l = l + 1
                enddo ! ypx
            enddo ! kpxi
        enddo ! apxi
        szt18vec56 = l - 1
        
        ! large scale entrepreneurs start with good business shock [t=1-8]
        l = 1
        do apxi = 1, 1
            kpx = 3
            opx = 1
            do ypx = 1, nmc
                t18vec7(l,1) = apxi 
                t18vec7(l,2) = kpx
                t18vec7(l,3) = ypx
                t18vec7(l,4) = opx
                l = l + 1
            enddo
        enddo
        szt18vec7 = l - 1
        
        ! --------------------------------------------------------- t = 9
        ! wage earners [t=9]
        l = 1
        do apxi = 1, 1 
            do opxi = 0, 1
                opx = merge(0,2,opxi==0) ! affixed with i means not directly to be used as an index
                kpx = merge(0,1,opxi==0)
                ypx = 0
                t9vec1(l,1) = apxi
                t9vec1(l,2) = kpx
                t9vec1(l,3) = ypx
                t9vec1(l,4) = opx
                l = l + 1
            enddo ! opxi
        enddo ! apxi        
        szt9vec1 = l - 1
        
        ! entrepreneurs start with bad business shock [t=9]
        l = 1
        do apxi = 1, 1
            do opxi = 0, 1
                !opx = 0
                !kpx = 0
                ypx = 0
                opx = merge(0,2,opxi==0)
                kpx = merge(0,1,opxi==0)
                t9vec24(l,1) = apxi
                t9vec24(l,2) = kpx
                t9vec24(l,3) = ypx
                t9vec24(l,4) = opx
                l = l + 1
            enddo ! opxi
        enddo ! apxi
        szt9vec24 = l - 1  
        
        ! small and medium scale entrepreneurs start with good business shock [t=9]
        l = 1
        do apxi = 1, 1
            do kpxi = 0, 1
                opx = merge(1,2,kpxi==0)
                ypx = 0
                t9vec56(l,1) = apxi 
                t9vec56(l,2) = kpxi ! not index directly
                t9vec56(l,3) = ypx
                t9vec56(l,4) = opx
                l = l + 1
            enddo
        enddo
        szt9vec56 = l - 1   
        
        ! large scale entrepreneurs start with good business shock [t=9]
        l = 1
        do apxi = 1, 1
            kpx = 3
            opx = 1
            ypx = 0
            t9vec7(l,1) = apxi 
            t9vec7(l,2) = kpx
            t9vec7(l,3) = ypx
            t9vec7(l,4) = opx
            l = l + 1
        enddo
        szt9vec7 = l - 1   
        
        ! --------------------------------------------------------- t = 10-13
        ! retired [t=10-13]
        l = 1
        do apxi = 1, 1 
            opx = 0
            kpx = 0
            ypx = 0
            t1013vec1(l,1) = apxi
            t1013vec1(l,2) = kpx
            t1013vec1(l,3) = ypx
            t1013vec1(l,4) = opx
            l = l + 1
        enddo ! apxi        
        szt1013vec1 = l - 1    
        
        ! entrepreneurs start with bad business shock [t=10-13]
        l = 1
        do apxi = 1, 1
            opx = 0
            kpx = 0
            ypx = 0
            t1013vec24(l,1) = apxi
            t1013vec24(l,2) = kpx
            t1013vec24(l,3) = ypx
            t1013vec24(l,4) = opx
            l = l + 1
        enddo ! apxi
        szt1013vec24 = l - 1  
        
        ! small and medium scale entrepreneurs start with good business shock [t=10-13]
        l = 1
        do apxi = 1, 1
            do kpxi = 0, 1
                opx = merge(1,2,kpxi==0)
                ypx = 0
                t1013vec56(l,1) = apxi 
                t1013vec56(l,2) = kpxi ! not index directly
                t1013vec56(l,3) = ypx
                t1013vec56(l,4) = opx
                l = l + 1
            enddo
        enddo
        szt1013vec56 = l - 1        
        
        ! large scale entrepreneurs start with good business shock [t=10-13]
        l = 1
        do apxi = 1, 1
            kpx = 3
            opx = 1
            ypx = 0
            t1013vec7(l,1) = apxi 
            t1013vec7(l,2) = kpx
            t1013vec7(l,3) = ypx
            t1013vec7(l,4) = opx
            l = l + 1
        enddo
        szt1013vec7 = l - 1     
        
        ! --------------------------------------------------------- t = 14
        ! retired t=14
        l = 1
        do apxi = 1, 1
            opx = 0
            kpx = 0
            ypx = 0
            t14vec1(l,1) = apxi
            t14vec1(l,2) = kpx
            t14vec1(l,3) = ypx
            t14vec1(l,4) = opx
            l = l + 1
        enddo ! apxi
        szt14vec1 = l - 1
        
        ! entrepreneurs start with bad business shock t=14
        l = 1
        do apxi = 1, 1
            opx = 0
            kpx = 0
            ypx = 0
            t14vec24(l,1) = apxi 
            t14vec24(l,2) = kpx
            t14vec24(l,3) = ypx
            t14vec24(l,4) = opx
            l = l + 1
        enddo
        szt14vec24 = l - 1
        
        ! small and medium scale entrepreneurs start with good business shock t=14
        l = 1
        do apxi = 1, 1
            kpx = 0
            ypx = 0
            opx = 0
            t14vec56(l,1) = apxi
            t14vec56(l,2) = kpxi
            t14vec56(l,3) = ypx
            t14vec56(l,4) = opx
            l = l + 1
        enddo
        szt14vec56 = l - 1
        
        ! large scale entrepreneurs start with good business shock t=14
        l = 1
        do apxi = 1, 1
            kpx = 0
            opx = 0
            ypx = 0
            t14vec7(l,1) = apxi 
            t14vec7(l,2) = kpx
            t14vec7(l,3) = ypx
            t14vec7(l,4) = opx
            l = l + 1
        enddo
        szt14vec7 = l - 1
        
    end subroutine make_next_period_index_combination    
    
    !subroutine get_weight_new_grids_for_refined_grid(av,afv,wint,afint) ! find interval on "corase" grid
    !    implicit none
    !    real(wp), dimension(:), intent(in) :: av, afv
    !    real(wp), dimension(:,:), intent(out) :: wint
    !    integer, dimension(:,:), intent(out) :: afint
    !    real(wp) :: a, b
    !    integer :: n, i
    !    n = size(afv)
    !    do i = 1, n
    !        afint(i,1) = locate(av,afv(i))
    !        afint(i,2) = afint(i,1) + 1
    !    enddo
    !    do i = 1, n
    !        a = av(afint(i,1))
    !        b = av(afint(i,2))
    !        wint(i,1) = (b-afv(i))/(b-a) 
    !        wint(i,2) = (afv(i)-a)/(b-a)
    !    enddo
    !    
    !    if(printout3)then
    !        write(unit=055,fmt='(x,a,x,a,2(2x,a),2(4x,a),2(4x,a))') 'no','fine level','1','2','val1','val2','wgt1','wgt2'            
    !        do i = 1, n
    !            write(unit=055,fmt='(i3,3x,f8.3,2i3,2f8.3,2f8.5)') i, afv(i), afint(i,1), afint(i,2), av(afint(i,1)), av(afint(i,2)), wint(i,1), wint(i,2)
    !        enddo
    !    endif
    !    
    !end subroutine get_weight_new_grids_for_refined_grid
        
end module variable