program MPI_sandbox
    !use universe ! test on reversion. Should shown.
    use equilibrium
    implicit none
    
    integer :: tstart, tend, trate, tmax, i, j      
    logical :: exit_log1
    
    call system_clock(tstart,trate,tmax)
    call system_clock(tstart)    
    
    !call start_files_for_writing() ! open files
    
    call fmpi_init() ! USER-DEFINED SUBROUTINE IN TOOLBOX.F90 <------------------------
    call infinity_setting(inf)
    if( MY_ID == 0 .and. MPI_PROVIDED<MPI_THREAD_FUNNELED ) write(*,'(a,/)') '! [ WARNING ] The Only-Master-makes-MPI-calls setup fails.' 
    
    allocate(range_guess(ndim,2))
    call read_parameter_model(para,'_1parameter.txt')
    if(my_id==0) write(*,'(a,f20.8)') (labstr(i),para(i),i=1,128) ! works. 
    
    ! <---- here --->
    allocate(indexseries(trylen),sobolm(nsbq, ndim))
    indexseries = [(i,i=1,trylen)]
    !if(my_id==0) write(*,'(i5)') (indexseries(i),i=1,trylen) ! works. 
    !if(my_id==0) write(*,'(a,/)') ' '
    
    indexseries = indexseries + sblno1 -1 ! USED FOR OUTPUT    
    !if(my_id==0) write(*,'(i5)') (indexseries(i),i=1,trylen)
    
    write(node_string,'(i2.2)') my_id
    solution_string = trim(node_string)//'.txt'       
    open(unit=my_id+1001, file=solution_string, action='write')
    
    write(trylen_string,'(i5.5,"_",i5.5)') sblno1, sblno1+trylen-1
    io_string = 'outputinput_'//trim(trylen_string) 
    
    ! MPI BLOCK
    call get_sobol_sequence( sobolm, 0.0_wp, 1.0_wp, nsbq, ndim )  
    !call scale_sobol_original( transpose(sobolm), range_guess, sobolm_scaled )    
    !mpi_sobol_scaled = sobolm_scaled(:,sblno1:sblno1+trylen-1)    
    
    deallocate(range_guess, indexseries)
    !call search_equilibrium(exit_log1) ! <===== replace solve_model() with this one. 3.10.2017 This is the working one.
    call mpi_finalize( MPI_ERR ) !! TO BE MOVED TO THE END OF THE MAIN PROGRAM <-------  
    
    call system_clock(tend) 
    if(my_id==0) write(*,fmt='(/,a,f12.4,a,x,i3)') 'total time: ',real(tend-tstart,wp)/real(trate,wp), ' seconds', my_id
    
    !call end_files_for_writing() ! close files  
        
    !! experiment good. ===========================================
    !! 3.8.2017 Brent can handle all types of tricky probelms I faces, and the user-defined subroutine brent_localizer is good. 
    !real(wp) :: xi, xo, xmax, ymax
    !real(wp), dimension(6) :: xvec
    !xi = -5._wp
    !xo = 1._wp
    !call grid(xvec,xi,xo,1._wp)
    !call brent_localizer(ftest,xi,xo,xmax,ymax) ! func,xa,xc,xs,ys
    !print*, xmax, ymax
    !
    !real(wp), dimension(:,:), allocatable :: mat
    !real(wp), dimension(:), allocatable :: xvec, yvec, f4
    !integer :: i,j
    !logical :: msg
    !
    !allocate( mat(5,3), yvec(5), xvec(3), f4(4) )
    !call read_matrix(mat,"test_z.txt")
    !do i = 1, 5
    !    write(*,*) (mat(i,j),j=1,3)
    !    yvec(i) = 1950._wp+(i-1)*10._wp
    !enddo
    !write(*,*) (yvec(i),i=1,5)
    !do i = 1, 3
    !    xvec(i) = 10._wp+(i-1)*10._wp
    !enddo
    !write(*,*) (xvec(i), i=1,3)
    !
    !write(*,*) binterpII(yvec,xvec,mat,1979._wp,12.7_wp,penalty,msg,f4)
    !deallocate( mat, xvec, yvec )
    !
    !integer :: i
    !real(wp) :: vec(4)
    !vec = [(i,i=1,4)]
    !print*, count(vec==2)
    !! =============================================================
    
contains    
    !function ftest(x) 
    !    implicit none
    !    real(wp) :: ftest
    !    real(wp), intent(in) :: x
    !    !ftest = x**4._wp + x**3._wp - 6*x**2._wp + 4*x + 12
    !    ftest = x*sin(1._wp/x)
    !end function ftest
    
    subroutine start_files_for_writing()
        implicit none
        !open(unit=003,file="output_3buz.txt",action="write",status="replace",recl=(7*20+10)) ! dev_output   
        !write(unit=003,fmt='(a)') 'ax,hx,kx,zx,yx,kpx,ypx,t,inc,rd*a,transfer, bizret(a,k,z,dk)'
        !open(unit=004,file="output_4lbd.txt",action="write",status="replace",recl=(7*20+10)) ! dev_output 
        !write(unit=004,fmt='(a)') 'ax,hx,kx,zx,yx,kpx,ypx,t,lbd, a1, a2, a3, a4'
        open(unit=005,file="output_05_Bequest_n_Value_period14.txt",action="write",status="replace",recl=(10*20+10)) ! dev_output 
        open(unit=006,file="output_06_value_function_beginning_of_period_14.txt",action="write",status="replace") ! dev_output 
        !open(unit=007,file="output_7prerefined.txt",action="write",status="replace") ! dev_output     
        !open(unit=008,file="output_8bizmat.txt",action="write",status="replace") ! dev_output : model     
        !open(unit=009,file="output_9income.txt",action="write",status="replace")
        !open(unit=011,file="output_11test.txt",action="write",status="replace")
        !open(unit=012,file="output_12_inputs_period14.txt",action="write",status="replace")
        open(unit=013,file="output_13_inputs_coarse_value_period14.txt",action="write",status="replace")
        !open(unit=014,file="output_14income.txt",action="write",status="replace")
        !open(unit=015,file="output_15_consumption.txt",action="write",status="replace")
        !open(unit=016,file="output_16_income_period13.txt",action="write",status="replace")    
        !open(unit=017,file="output_17_corase_brent_max_indices.txt",action="write",status="replace")  
        !open(unit=018,file="output_18_corase_value_matrix_period_13.txt",action="write",status="replace")  
        !open(unit=019,file="output_19_corase_value_result_period_13.txt",action="write",status="replace")      
        !open(unit=020,file="output_20_value_function_period_14.txt",action="write",status="replace")
        !open(unit=021,file="output_21_wealth_decision.txt",action="write",status="replace")
        !open(unit=022,file="output_22_inputs_decision.txt",action="write",status="replace")
        open(unit=023,file="output_23_negative_tech_shock_impacts_on_value_function.txt",action="write",status="replace")
        open(unit=024,file="output_24_collateral.txt",status="replace",action="write")    
        !open(unit=025,file="output_25_collateral_binding.txt",status="replace",action="write")    
        open(unit=026,file="output_26_value_function_beginning_of_period.txt",status="replace",action="write")   
        open(unit=027,file="output_27_list_twa.txt",status="replace",action="write")   
        open(unit=028,file="output_28_list_twh.txt",status="replace",action="write")   
        open(unit=029,file="output_29_list_twk.txt",status="replace",action="write")   
        !open(unit=030,file="output_30_all_beginning_value_function.txt",status="replace",action="write")
        !open(unit=031,file="output_31_reconstruct_period1014_beginning_value.txt",status="replace",action="write")
        !open(unit=032,file="output_32_wf.txt",status="replace",action="write")
        !open(unit=033,file="output_33_corase_brent_max_indices.txt",status="replace",action="write")
        !open(unit=034,file="output_34_valuefunction_period8.txt",status="replace",action="write")
        !open(unit=035,file="output_35_twf_period8.txt",status="replace",action="write")
        open(unit=036,file="output_36_list_tww.txt",status="replace",action="write")
        !open(unit=037,file="output_37_read_decisions.txt",status="replace",action="write")
        !open(unit=038,file="output_38_weight.txt",status="replace",action="write")
        open(unit=039,file="output_39_coarse_twa.txt",status="replace",action="write")
        open(unit=040,file="output_40_coarse_twh.txt",status="replace",action="write")
        open(unit=041,file="output_41_coarse_twk.txt",status="replace",action="write")
        open(unit=042,file="output_42_coarse_tww.txt",status="replace",action="write")
        !open(unit=043,file="output_43_ref_ainc.txt",status="replace",action="write")
        !open(unit=044,file="output_44_ref_tax.txt",status="replace",action="write")
        !open(unit=045,file="output_45_ref_beq.txt",status="replace",action="write")
        !open(unit=046,file="output_46_ref_dist.txt",status="replace",action="write")
        !open(unit=047,file="output_47_combination.txt",status="replace",action="write")
        !open(unit=048,file="output_48_mass.txt",status="replace",action="write")
        open(unit=049,file="output_49_macro_result.txt",status="replace",action="write")
        open(unit=050,file="output_50_list_twf.txt",status="replace",action="write")
        !open(unit=051,file="output_51_debug.txt",status="replace",action="write")
        !open(unit=052,file="output_52_debug.txt",status="replace",action="write")
        !open(unit=053,file="output_53_id1.txt",status="replace",action="write")
        !open(unit=054,file="output_54_recycle.txt",status="replace",action="write")
        !open(unit=055,file="output_55_afint.txt",status="replace",action="write")    
        open(unit=056,file="output_56_cwf_value_function_end_of_period.txt",status="replace",action="write")    
        !open(unit=057,file="output_57_msg.txt",status="replace",action="write")    
        !open(unit=058,file="output_58_sid.txt",status="replace",action="write")
        open(unit=059,file="output_59_18twh.txt",status="replace",action="write")    
        open(unit=060,file="output_60_18twa.txt",status="replace",action="write")    
        open(unit=061,file="output_61_refined_Tbim.txt",status="replace",action="write")    
        open(unit=062,file="output_62_refined_Atwm.txt",status="replace",action="write")    
        open(unit=063,file="output_63_refined_Taxm.txt",status="replace",action="write")    
        open(unit=064,file="output_64_refined_Beqm.txt",status="replace",action="write")    
        open(unit=065,file="output_65_refined_Ldemm.txt",status="replace",action="write")    
        open(unit=066,file="output_66_refined_Lsupm.txt",status="replace",action="write")    
        open(unit=067,file="output_67_refined_Prom.txt",status="replace",action="write")    
        open(unit=068,file="output_68_refined_Cspm.txt",status="replace",action="write")    
        open(unit=069,file="output_69_refined_uwa.txt",status="replace",action="write")    
        open(unit=070,file="output_70_refined_uwh.txt",status="replace",action="write")    
        open(unit=071,file="output_71_refined_uwk.txt",status="replace",action="write")    
        open(unit=072,file="output_72_refined_uww.txt",status="replace",action="write")    
        open(unit=073,file="output_73_refined_wok2ent.txt",status="replace",action="write")    
        open(unit=074,file="output_74_refined_ent2wok.txt",status="replace",action="write")    
        open(unit=075,file="output_75_ide_matrix.txt",status="replace",action="write")   
        !open(unit=076,file="output_76_mass.txt",status="replace",action="write")   
        !open(unit=077,file="output_77_mass.txt",status="replace",action="write")
        !open(unit=078,file="output_78_mass.txt",status="replace",action="write")
        !open(unit=079,file="output_79_mass.txt",status="replace",action="write")
        open(unit=080,file="output_80_twwvec.txt",status="replace",action="write")
        open(unit=081,file="output_81_nsav.txt",status="replace",action="write")
        open(unit=082,file="output_82_dist_error_by_round.txt",status="replace",action="write")
        open(unit=083,file="output_83_end_dist_converge_by_period.txt",status="replace",action="write")
        open(unit=084,file="output_84_coarse_tbi_m.txt",status="replace",action="write")
        open(unit=085,file="output_85_coarse_atw_m.txt",status="replace",action="write") 
        open(unit=086,file="output_86_coarse_tax_m.txt",status="replace",action="write") 
        open(unit=087,file="output_87_coarse_beq_m.txt",status="replace",action="write") 
        open(unit=088,file="output_88_coarse_pro_m.txt",status="replace",action="write") 
        open(unit=089,file="output_89_coarse_ldem_m.txt",status="replace",action="write")
        open(unit=090,file="output_90_coarse_lsup_m.txt",status="replace",action="write")
        open(unit=092,file="output_92_coarse_csp_m.txt",status="replace",action="write")     
        open(unit=093,file="output_93_sid_sum_per_period.txt",status="replace",action="write")         
        open(unit=094,file="output_94_inv_dist_moved_components.txt",status="replace",action="write",recl=1024)         
        open(unit=095,file="output_95_beq_dist_invalid_recipients.txt",status="replace",action="write") 
        open(unit=097,file="output_97_inv_dist_current_abnormal_positive_mass.txt",status="replace",action="write") 
        !open(unit=098,file="output_98_particular_missing.txt",status="replace",action="write") 
        !open(unit=099,file="output_99_nafv.txt",status="replace",action="write") 
        open(unit=101,file="output_101_missing_housing.txt",status="replace",action="write") 
        open(unit=102,file="output_102_convergence.txt",status="replace",action="write") 
        open(unit=103,file="output_103_macro.txt",status="replace",action="write") 
        open(unit=104,file="output_104_debug.txt",status="replace",action="write") 
        open(unit=105,file="output_105_debug.txt",status="replace",action="write")
        !open(unit=106,file="output_106_debug.txt",status="replace",action="write")
        !open(unit=107,file="output_107_profitcompare.txt",status="replace",action="write")
        !open(unit=108,file="output_108_labor.txt",status="replace",action="write")
        open(unit=109,file="output_109_profit_mat.txt",status="replace",action="write")
        open(unit=110,file="output_110_detail_business_mat.txt",status="replace",action="write")
        open(unit=111,file="output_111_detail_housing_mat.txt",status="replace",action="write")
        open(unit=112,file="output_112_detail_investment_mat.txt",status="replace",action="write")
        open(unit=113,file="output_113_detail_financial_mat.txt",status="replace",action="write")
        open(unit=114,file="output_114_detail_smat.txt",status="replace",action="write")
        open(unit=115,file="output_115_coarse_cww.txt",status="replace",action="write")
        open(unit=116,file="output_116_coarse_cwf.txt",status="replace",action="write")
        open(unit=117,file="output_117_coarse_cwa.txt",status="replace",action="write")
        open(unit=118,file="output_118_coarse_cwh.txt",status="replace",action="write")
        open(unit=119,file="output_119_coarse_cwk.txt",status="replace",action="write")
        open(unit=120,file="output_120_output.txt",status="replace",action="write")        
        open(unit=121,file="output_121_stata.txt",status="replace",action="write")
        open(unit=122,file="output_122_mass.txt",status="replace",action="write")
        open(unit=123,file="output_123_coarse_cwc.txt",status="replace",action="write")
        open(unit=124,file="output_124_coarse_cef.txt",status="replace",action="write")
        open(unit=125,file="output_125_coarse_cvv.txt",status="replace",action="write")        
        open(unit=126,file="output_126_image_cef.txt",status="replace",action="write")  
        open(unit=127,file="output_127_aggregate_stats.txt",status="replace",action="write")
        open(unit=128,file="output_128_mass_transition.txt",status="replace",action="write")
        open(unit=129,file="output_129_cef.txt",status="replace",action="write")
        open(unit=130,file="output_130_sw_nonlineartax.txt",status="replace",action="write")
        open(unit=131,file="output_131_debug.txt",status="replace",action="write")
        open(unit=132,file="output_132_macro_stats.txt",status="replace",action="write")
        open(unit=133,file="output_133_debug_medinc.txt",status="replace",action="write")
    end subroutine start_files_for_writing
    
    subroutine end_files_for_writing
        implicit none
        !close(003)
        !close(004)
        close(005)
        close(006)
        !close(007)
        !close(008)
        !close(009)
        !close(011)
        !close(012)
        close(013)
        !close(014)
        !close(015)
        !close(016)
        !close(017)
        !close(018)
        !close(019)
        !close(020)
        !close(021)
        !close(022)
        close(023)
        close(024)
        !close(025)
        close(026)
        close(027)
        close(028)
        close(029)
        !close(030)
        !close(031)
        !close(032)
        !close(033)
        !close(034)
        !close(035)
        close(036)
        !close(037)
        !close(038)
        close(039)
        close(040)
        close(041)
        close(042)
        !close(043)
        !close(044)
        !close(045)
        !close(046)
        !close(047)
        !close(048)
        close(049)
        close(050)
        !close(051)
        !close(052)
        !close(053)
        !close(054)
        !close(055)
        close(056)
        !close(057)
        !close(058)
        close(059)
        close(060)
        close(061)
        close(062)
        close(063)
        close(064)
        close(065)
        close(066)
        close(067)
        close(068)
        close(069)
        close(070)
        close(071)
        close(072)
        close(073)
        close(074)
        close(075)
        !close(076)
        !close(077)
        !close(078)
        !close(079)
        close(080)
        close(081)
        close(082)
        close(083)
        close(84)
        close(85)
        close(86)
        close(87)
        close(88)
        close(89)
        close(90)
        close(92) 
        close(93)
        close(94)
        close(95)
        close(97)
        !close(98)
        !close(99)
        close(101)
        close(102)
        close(103)
        close(104)
        close(105)
        !close(106)
        !close(107)
        !close(108)
        close(109)
        close(110)
        close(111)
        close(112)
        close(113)
        close(114)
        close(115)
        close(116)
        close(117)
        close(118)
        close(119)
        close(120)   
        close(121)
        close(122)
        close(123)
        close(124)
        close(125)
        close(126)
        close(127)
        close(128)
        close(129)
        close(130)
        close(131)
        close(132)
        close(133)
    end subroutine end_files_for_writing
    
end program MPI_sandbox
    
    !integer, dimension(:,:), allocatable :: am, bm, cm
    !integer :: i, j
    !allocate( am(2,2), bm(2,2), cm(2,2) )
    !am = reshape([(i,i=2,5)],(/2,2/))
    !bm = reshape([(i,i=12,15)],(/2,2/))
    !cm = 0
    !do j = 1, 2
    !    print*, (am(j,i),i=1,2)    
    !enddo
    !do j = 1, 2
    !    print*, (bm(j,i),i=1,2)    
    !enddo    
    !where(bm>14) cm=am
    !do j = 1, 2
    !    print*, (cm(j,i),i=1,2)    
    !enddo        
    !deallocate( am, bm, cm )
    
    !real(wp) :: vec(0:3)
    !integer :: i
    !vec = [(real(i,wp),i=1,4)]
    !print*, minloc(vec)
    
    !integer :: a(3), b(3)
    !integer :: i
    !a = [(i,i=1,3)]
    !b = [(i*2,i=1,3)]
    !print*, any(a>2), any(a>4)
    
    !real(wp) :: am(2,4), bm(2,3), cm(4,12)
    !real(wp) :: dv(8), dv1(6)
    !integer :: i
    !dv = [(real(i,wp),i=1,8)]
    !dv1= [(real(i,wp),i=1,6)]
    !am = reshape(dv,(/2,4/))    
    !bm = reshape(dv1,(/2,3/))
    !call kron(am,bm,cm)
    !call sm(cm,'cm')
    
    !!! string passing
    !integer :: itemp
    !character(len=100) :: stemp
    !itemp = 1
    !call test_string(itemp,stemp)
    !print*, trim(stemp)
    !call test_string(2,stemp)
    !print*, ' if test ', stemp=='west'
    !print*, trim(stemp), len(stemp)
    
   

    ! TRICK: TO PRINT OUTCOME INSIDE INVOKED SUBROUTINES OR FUNCTIONS, DELCARE FILE-OPEN STATEMENTS IN THE MAIN PROGRAM.
    ! TRICK: THE FORMATE OF RESPECITVE OUTPUT FILE CAN RESIDE IN THE INVOKED SUBROUTINES AND FUNCTION
    ! TRICK: THE CLOSE-FILE STATEMENTS NEED TO IN THE SAME SCOPE AS THE CORRESPONDING OPEN-FILE STATMENTS.

    
    !real(wp), dimension(:), allocatable :: xray, yray
    !allocate(xray(5), yray(5))
    !xray(1) = 1.2_wp 
    !xray(2) = 0.6_wp
    !xray(3) = 0.7_wp
    !xray(4) = 1.5_wp
    !xray(5) = 1.2_wp
    !yray = 10._wp*xray
    !write(unit=*,fmt='(5f6.3)'), xray
    !write(unit=*,fmt='(5f6.3,/)'), yray
    !
    !call sort2arrayII(xray,yray,'descending')
    !
    !write(unit=*,fmt='(5f6.3)'), xray
    !write(unit=*,fmt='(5f6.3)'), yray    
    !deallocate( xray, yray )
    
    !character(len=3) :: length_para  
    
    !real(wp) :: x(2)
    !x(1) = -1._wp
    !x(1) = sqrt(x(1))
    !x(2) = 100._wp
    !print*, ' maximum is ', maxval(x)
    
    !! LOCATE SUBROUTINE
    !real(wp), dimension(:), allocatable :: tvector
    !allocate( tvector(4) )
    !tvector = [2._wp,3._wp,4._wp,5._wp]
    !print*,'the answer is ',locate(tvector,5._wp)
    !deallocate( tvector )
    
    !! RESHAPE
    !integer :: i
    !integer, dimension(:,:), allocatable :: tarray
    !integer, dimension(:), allocatable :: tvector
    !allocate( tarray(2,3), tvector(6) )
    !tvector = [(i,i=1,6)]
    !tarray = reshape(tvector,(/2,3/))
    !call smi(tarray,'tarray')
    !deallocate( tvector, tarray )
    
    !!! MINLOC 
    !!integer, dimension(:), allocatable :: testvec
    !!integer :: i
    !!allocate( testvec(10))
    !!testvec = [(i,i=1,10)]
    !!print*, testvec
    !!print*, minloc(testvec,testvec>4)
    
    !! LOCATE -99
    !integer :: i
    !real(wp), dimension(:), allocatable :: tvector
    !allocate( tvector(10) )
    !tvector = [(real(i,wp),i=1,10)]
    !where(tvector>8._wp.or.tvector<4._wp) tvector = -99._wp ! trick
    !call ss(tvector,'tvector')
    !print*, 'count -99: ', count(tvector==-99._wp) ! answer: 5
    !print*, 'minimum loc: ', minloc(tvector,tvector>-99._wp) ! the first non-negative 99 element. Answer: 4
    !print*, 'maximum loc: ', maxloc(tvector,tvector>-99._wp) ! the last non-negative 99 element. Answer: 8     
    !deallocate( tvector )  
    
    !! PRINT OUT TEST
    !integer :: i, j
    !integer, dimension(:,:,:), allocatable :: imat
    !allocate( imat(3,2,4) )
    !imat = 100
    !do j = 1, 3
    !    imat(j,:,:) = reshape([(i,i=1+j*10,8+j*10)],(/2,4/))
    !enddo
    !call smi(imat(1,:,:),'test_01')
    !deallocate( imat )
    
    !! REVERSE FOR THE INDEX
    !integer, dimension(:), allocatable :: tvector
    !allocate( tvector(4) )
    !tvector = [1,2,3,4]
    !print*, '1: ', tvector<=2
    !print*, '2: ',count(tvector<=2) ! Use this one
    !deallocate( tvector )    
    
    !! EXPERIMENT OF PACKING
    !integer, dimension(2,3,2) :: mat1
    !mat1 = reshape([(i,i=1,12)],(/2,3,2/))
    !do i = 1, 2
    !    write(*,fmt='(3i3)'), (mat1(i,j,1),j=1,3)
    !enddo
    !
    !do i = 1, 2
    !    write(*,fmt='(3i3)'), (mat1(i,j,2),j=1,3)
    !enddo    
    !vec1(1:1) = pack(mat1(1,1,:),mask=mat1(1,1,:)==1) ! NOTE: YOU NEED TO USE 1:1 TO RECEIVE ONE ELEMENT USING PACKING!!!
    !print*, ' the elemeing packed into the vector is ', pack(mat1(1,1,:),mask=mat1(1,1,:)==1)    
    
    !integer :: am(3,4)
    !am = spread((/1,2,3/),dim=2,ncopies=4)
    !do i = 1, 3
    !    write(*,fmt='(4i3)') (am(i,j),j=1,4)
    !enddo
    !print*, pack(am,am==2)
    
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
    
    !real(wp) :: vec(3)
    !vec = [(real(i,wp),i=1,3)]
    !print*, '0.9', locate(vec,0.9_wp)
    !print*, '1.0', locate(vec,1._wp)
    !print*, '2.0', locate(vec,2._wp)
    !print*, '3.0', locate(vec,3._wp)
    !print*, '3.5', locate(vec,3.5_wp)    
    
    !!! maxloc and minloc
    !integer, dimension(:), allocatable :: jvector
    !allocate( jvector(0:5) )
    !jvector = [(i,i=-3,2,1)]
    !write(*,'(6(i3),x,a,x,i3)') jvector,'minloc',minloc(jvector,jvector>0)
    !write(*,'(6(i3),x,a,x,i3)') jvector,'maxloc',maxloc(jvector,jvector>0)
    !deallocate( jvector )         
    
    !!! OPENMP EXPERIMENT 10042016
    !integer :: att, supp
    !!! OPENMP TEST 1
    !supp = 0
    !!$omp parallel default(shared) private(att) reduction(+:supp)
    !    att = 1
    !    supp = supp + &
    !        & att
    !!$omp end parallel
    !print*, ' TEST 1 ', supp
    !!! OPENMP TEST 2
    !supp = 0
    !!$omp parallel default(shared) private(att)
    !    att = 1
    !    !$omp atomic update
    !    supp = supp + &
    !        & att
    !!$omp end parallel
    !print*, ' TEST 2 ', supp    