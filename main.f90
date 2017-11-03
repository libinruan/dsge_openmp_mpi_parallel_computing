! (1) main.f90 focuses on the implementation of the Nelder-Mead simplex method on 10 dimensional parameter space.
! (2) equilibrium.f90 solves for the steady state of the economy.
! (3) model.f90 solves the Bellman equations for saving and consumption decisions across generations.
! (4) variable.f90 initializes most variables to be used in the program.
! (5) toolbox.f90 contains all the subroutines that executes evoked numerical methods.
!
! Note that AMOEBA_WORLD can be larger or smaller than CONTRACT_WORLD; however, ROW_WORLD must be the larger than either AMOEBA_WORLD or CONTRACT_WORLD.
! "nrow" = the size of ROW_WORLD, "4*noamoeba" = the size of AMOEBA_WORLD, and "ndim+1" = the size of CONTRACT_WORLD.
    
program MPI_sandbox
    !use universe ! test on reversion. Should shown.
    use equilibrium
    implicit none
    integer :: tstart, tend, trate, tmax, i, j, index_pt, index_rcv
    ! character(len=200) :: msg
    character(:), allocatable :: msg
    real(wp) :: temp(2), accu ! debug
    
    ! The mpi_exercise_mode == -2: experiment zone
    real(wp), dimension(:), allocatable :: temp1vertex
    real(wp), dimension(:,:), allocatable :: temp1matrix
    
    ! The mpi_exercise_mode == 1 case: Coarse search
    logical :: exit_log1
    
    ! The mpi_exercise_mode == -1 case: MKL experiment of random generators
    integer :: generator,erridx, approach, tmiddle
    type(vsl_stream_state) :: river
    real(wp) :: stoch(1)
    real(wp), dimension(:,:), allocatable :: mat_stage1_inputs
    
    character(:), allocatable :: stringmode6
    character(len=5) :: idmode6    
    
    call system_clock(tstart,trate,tmax)
    call system_clock(tstart)    
    
    !call start_files_for_writing() ! open files. 8-26-2017 turned on
    
    call fmpi_init() ! USER-DEFINED SUBROUTINE IN TOOLBOX.F90 <------------------------
    call infinity_setting(inf)
    
    if( MY_ID == 0 .and. MPI_PROVIDED<MPI_THREAD_FUNNELED ) write(*,'(a,/)') '! [ WARNING ] The Only-Master-makes-MPI-calls setup fails.' 
    
    allocate(range_guess(ndim,2)) ! 10.21.2017 Warning: it has to be placed prior to read parameters of the model (_parameter1.txt)
    call read_parameter_model(para,'_1parameter.txt')
    
    if(my_id==0)then ! General Operation Messages
        
        if(printout12)then
            write(*,'(a,f20.8)') (labstr(i),para(i),i=1,148) ! works. 
            write(*,*) ' '
            write(*,*) "printout17 everyone recevies bequests: ", printout17
            write(*,*) "printout18 housing upper limit extended based on consumer need: ", printout18
            write(*,*) "printout19: ", printout19
            write(*,*) "printout20: ", printout20
            write(*,*) "printout24: ", printout24
        endif
        
        if(mpi_exercise_mode==1)then
            msg = ' starting point search'
        elseif(mpi_exercise_mode==2)then
            msg = ' amoeba algorithm'
        elseif(mpi_exercise_mode==3)then
            msg = ' single node relay'
        else
            msg = ' -- not yet --'    
        endif
        
        if(mpi_exercise_mode/=0) write(*,'(2a)') 'MPI_EXERCISE_MODE: ', msg
        
        if(printout11)then
            msg = ' quantitative model'
        else
            msg = ' test model'
        endif
        write(*,'(2a)') 'PRINTOUT11: ', msg
        write(*,'(a,i4)') 'obj_func_toggle: ', obj_func_toggle
        write(*,'(a,i4)') 'mode4list: ', mode4list
        !call test2()
        
    endif

    if(mpi_exercise_mode == -2)then ! experiment zone
        
        if(my_id==0)then
            
            !real(wp), dimension(:), allocatable :: temp1vertex
            !real(wp), dimension(:,:), allocatable :: temp1matrix        
            allocate(temp1vertex(ndim),temp1matrix(ndim,subdim+1))
            
            generator = VSL_BRNG_MCG31
            approach  = VSL_RNG_METHOD_UNIFORM_STD
            call system_clock(tmiddle)
            erridx    = vslnewstream(river,generator,tmiddle)
            do i = 1, ndim         
                erridx = vdrnguniform(approach,river,1,stoch,-20._wp,20._wp)
                write(*,'(a,i3,a,(f8.2,x))') 'my_id = ', my_id, ' random number = ', stoch
                temp1vertex(i) = stoch(1)
            enddo        
            write(*,'(a,<ndim>(f8.2,x))') ' temp1vertex: ', temp1vertex
            
            call point_boundary_adjustment(temp1vertex,range_guess)
            write(*,'(a,<ndim>(f8.2,x))') ' temp1vertex: ', temp1vertex
            
            do j = 1, subdim+1
                do i = 1, ndim         
                    erridx = vdrnguniform(approach,river,1,stoch,-20._wp,20._wp)
                    !write(*,'(a,i3,a,(f8.2,x))') 'my_id = ', my_id, ' random number = ', stoch
                    temp1matrix(i,j) = stoch(1)
                enddo             
            enddo
            
            write(*,'(a)') 'Before'
            do j = 1, ndim
                write(*,'(<subdim+1>(f8.2,x))') temp1matrix(j,:) 
            enddo
            write(*,'(a)') 'After'
            call matrix_boundary_adjustment(temp1matrix,range_guess)
            do j = 1, ndim
                write(*,'(<subdim+1>(f8.2,x))') temp1matrix(j,:) 
            enddo
            
            deallocate(temp1vertex,temp1matrix)            
            
            ! write(*,*) ' locate test ', locate(real(breaks_list,wp),1.5_wp,1)
            ! write(*,*) ' locate test ', locate(real(breaks_list,wp),50._wp,1)
            ! 
            ! call read_best_point( guessv, bestvertex_file )
            ! 
            ! write(*,*) ' Filename: ', bestvertex_file
            ! write(*,*) ' Best vertex: ', guessv
            ! 
            ! write(*,'(a,x,i6)') 'Trial:', trial
            ! 
            ! do i = 1, ndim-1
            !     write(*, '((i2),2x)', advance='no') i            
            ! enddo ! i
            ! write(*, '(i2,2x)') ndim
            ! 
            ! do i = 1, ndim-1
            !     write(*, '(i2,2x)', advance='no') i
            ! enddo             
            ! write(*, '(i2,2x)') ndim
            ! write(*,'(a)') " test " 
            
            !! experiment 2, OK.
            !allocate( mat_stage1_inputs(1728,11) )
            !call read_matrix(mat_stage1_inputs,'_stage1_input.csv')
            !do i = 1, size(mat_stage1_inputs,dim=1)
            !    write(*,'(11(f8.3,x))') mat_stage1_inputs(i,:)
            !enddo ! i
            !deallocate( mat_stage1_inputs )
            
        endif
  
    elseif(mpi_exercise_mode==-1)then ! MKL experiment (random number generator)
        
        generator = VSL_BRNG_MCG31
        approach  = VSL_RNG_METHOD_UNIFORM_STD
        call system_clock(tmiddle)
        erridx    = vslnewstream(river,generator,tmiddle)
        do i = 1, 20         
            erridx = vdrnguniform(approach,river,1,stoch,0._wp,1000._wp)
            write(*,*) 'my_id = ', my_id, ' random number = ', stoch
        enddo
        
        ! accuracy experiment
        accu = 1.e8_wp
        temp = [0.1234567_wp, 0.6666666666_wp]
        print*, "temp1: ", int(temp*accu)
        print*, "temp2: ", int(temp*accu)/accu
        
        ! rounding experiment
        print*, nint(temp(1)*1e2)/1e2
        print*, nint(temp(1)*1e3)/1e3 ! 0.123
        print*, nint(temp(1)*1e4)/1e4 ! 0.1235
        print*, nint(temp(1)*1e5)/1e5 ! 0.12346
        print*, nint(temp(1)*1e6)/1e6 ! 0.123457
        print*, nint(temp(1)*1e7)/1e7 ! rounding to the 7th place to the right of decimal point
        
        print*, nint(temp*1e6)/1e6
        
    elseif(mpi_exercise_mode==0)then ! Stage 0. Building Stage with only a single node.
        
        solution_string = 'SingleNode_1001.txt'   
        concisesolution_string = 'SingleNode_2001.txt'
        open(unit=my_id+1001, file=solution_string, action='write', position='append')
        open(unit=my_id+2001, file=concisesolution_string, action='write', position='append')
        !if(printout3) write(unit=my_id+1001,fmt='(a)') ' ------------------------------------------------- '
        
        ! Directly use the parameter setting in _1parameter.txt
        guessv(1) = kv1   
        guessv(2) = prtk0 
        guessv(3) = prtk1 
        guessv(4) = prtk2 
        guessv(5) = zbar  
        guessv(6) = beta  
        guessv(7) = theta 
        !guessv(7) = iota
        guessv(8) = phi1  
        guessv(9) = phi2  
        guessv(10)= phi3     
        
        modelmsg = 0
        momvec = inf
        obj_val_1st = inf
        
        ! Don't update the parameter setting as do the other part of this program. Let's use the modified guessv array.
        call search_equilibrium( guessv, momvec, obj_val_1st, my_id, my_id, modelmsg )
        if( modelmsg == 0 )then
            write(my_id+1001, '(a,<ndim>f15.7)') 'guess  : ', guessv
            write(my_id+1001, '(a,<ndim>f15.7)') 'targetv: ', targetv
            write(my_id+1001, '(a,<ndim>f15.7)') 'moment : ', momvec
            write(my_id+1001, '(a,f15.7)') 'penalty: ', obj_val_1st 
        else ! fail to solve the model
            write(my_id+1001, '(a,a,<ndim>f15.7)') ' === Failure === ', 'guess  : ', guessv 
        endif ! modelmsg
        
        close(my_id+1001)  
        close(my_id+2001)
        
    elseif(mpi_exercise_mode==1)then ! Stage 1. Search for the optimal parameter setting on the global parameter space without the Nelder-Mead algorithm
        ! Used for all invoked nodes in the MPI implementation. 
        write(node_string,'(i3.3)') my_id
        
        if(my_id/=0)then
            solution_string = 'SlaveIntermFeedback_'//trim(node_string)//'.txt' ! 7-7-2017 Intermediate outcome of individual slave node.         
            concisesolution_string = 'SlaveMacrostat_'//trim(node_string)//'.txt' ! 7-7-2017 Macro stat of individual slave node.        
        else
            write(trylen_string,'(i5.5,"_",i5.5)') sblno1, sblno1+trylen-1
            solution_string = 'RootIntermCollect_'//trim(trylen_string)//'.txt' ! 7-7-2017 Real-time feedback from the slaves.        
            concisesolution_string = 'RootUseless_'//trim(trylen_string)//'.txt' ! 7-7-2017 Useless for the root node.        
            ! add concise solution
        endif
        open(unit=my_id+1001, file=solution_string, action='write') ! Moved here. 7-3-201
        open(unit=my_id+2001, file=concisesolution_string, action='write') ! Moved here. 7-5-201
        
        ! # 1 move inside the MPI_exercise_mode == 1!?
        allocate(indexseries(trylen),sobolm(nsbq, ndim),sobolm_scaled(ndim,nsbq),mpi_sobol_scaled(ndim,trylen))
        allocate(mpi_simmom_matrix(ndim,trylen),origin_input(ndim),mpi_sobol_mixed(trylen,ndim))
        !allocate(indexseries(trylen),sobolm(ndim,nsbq))
        indexseries = [(i,i=1,trylen)]
        !if(my_id==0) write(*,'(i5)') (indexseries(i),i=1,trylen) ! works. 
        !if(my_id==0) write(*,'(a,/)') ' '
        indexseries = indexseries + sblno1 -1 ! Shifted; USED FOR OUTPUT    
        !if(my_id==0) write(*,'(i5)') (indexseries(i),i=1,trylen) 
        
        origin_input(1) = kv1   
        origin_input(2) = prtk0 
        origin_input(3) = prtk1 
        origin_input(4) = prtk2 
        origin_input(5) = zbar  
        origin_input(6) = beta  
        origin_input(7) = theta 
        !origin_input(7) = iota
        origin_input(8) = phi1  
        origin_input(9) = phi2  
        origin_input(10)= phi3           
        
        ! Quasi-random Sobol sequence block # 2 move inside the MPI_exercise_mode == 1!?
        call get_sobol_sequence( sobolm, 0.0_wp, 1.0_wp, nsbq, ndim ) ! Generate ndim dimensional sobol sequence of length nsbq (nsbq*ndim).
        ! call scale_sobol_original( transpose(sobolm), range_guess, sobolm_scaled ) ! ndim*nsbq ! comment out 8-18-2017
        
        sobolm_scaled = transpose(sobolm) ! added 8-18-2017
        mpi_sobol_scaled = sobolm_scaled(:,sblno1:sblno1+trylen-1) ! 7-29-2017 Shifted.    
        
        if(my_id==0) call sm(sobolm,'sobolm'//trim(trylen_string)) ! checked 2017-Jul-1
        if(my_id==0) call sm(transpose(mpi_sobol_scaled),'mpi_sobol_scaled'//trim(trylen_string)) ! checked 2017-Jul-1 . Used for convex combination.        
        
        allocate( parcel(ndim), result(ndim), outputinput1(trylen,2*ndim+2), obj_val_vec(trylen) )        
        nslaves = num_procs - 1 ! Just for mpi_exercise_mode==1 case, the root processor (my_id==0) is in charge of sending out new trial and receiving the corresponding result.

        if( my_id == 0)then
            
            write(trylen_string,'(i5.5,"_",i5.5)') sblno1, sblno1+trylen-1
            
            !! 8-18-2017
            io_string = 'FinalScaledSobel_'//trim(trylen_string) 
            do i = 1, trylen
                call linear_combination_sobol_sequence_corrected( parcel, i, mpi_sobol_scaled(:,i), range_guess )
                mpi_sobol_mixed(i,:) = parcel
            enddo
            call sm(mpi_sobol_mixed,io_string)
            
            io_string = 'IOMat_'//trim(trylen_string) 
            ! [The Root, case 1] Send Initial Messages to Slave Nodes (indices ranges from 1 to nslaves)
            do i = 1, nslaves ! nslaves = num_procs-1.
                trial = i ! the index of the basic loop: 1,...,trylen; Not the indices of the adjusted Sobol sequence.
                slave = i
                
                !! comment out 8-18-2017
                !call linear_combination_sobal_sequence(parcel,trial,mpi_sobol_scaled(:,trial),origin_input,weight_list,breaks_list) ! Be sure to set up weight_list and breaks_list.    
                
                !! added 8-18-2017, then moved outside this do loop over i between 1 and nslaves
                !call linear_combination_sobol_sequence_corrected( parcel, trial, mpi_sobol_scaled(:,trial), range_guess )
                !! The formula is: vertex_output = (range_guess(:,2)+range_guess(:,1))/2._wp + (unit_sobol_input-0.5_wp)*(range_guess(:,2)-range_guess(:,1))
                !mpi_sobol_mixed(trial,:) = parcel ! bookkeeping. ! comment out 8-18-2017
                
                parcel = mpi_sobol_mixed(i,:)
                call sendjob(trial,slave,parcel) ! send from root to slave the trial parameter combination                
            enddo
            
            ! [The Root, case 1] Hear Responses from Slave Nodes
            do i = 1, trylen ! The range is correct!! 7-3-2017 ! I controls the total number of outcome expected to receive. 8-2-2017
                do
                    msgtype = 2 ! tag 1 for parameter passing; tag 2 for communicating the result between slave nodes and root. 
                    ! Non-blocking test for a message
                    call mpi_iprobe( mpi_any_source, msgtype, mpi_comm_world, &
                        & receiving, status, ierr)
                    if(receiving)then
                        ! 1. [slave] Which member node is sending the result? 
                        call mpi_recv( slave, 1, mpi_integer, mpi_any_source, &
                            & msgtype, mpi_comm_world, status, ierr)
                        ! 2. [trial] what's the returned trial index in the basic assignment loop?
                        call mpi_recv( trial, 1, mpi_integer, slave, & ! Note: trial falls in [1,trylen], rather than the shifted interval.
                            & msgtype, mpi_comm_world, status, ierr)
                        ! 3. [result] the feedback of simulated moments 
                        call mpi_recv( result, ndim, mpi_double_precision, slave, &
                            & msgtype, mpi_comm_world, status, ierr)
                        ! 4. [obj_val_1st] the value of penalty corresponds to the given trial
                        call mpi_recv( obj_val_1st, 1, mpi_double_precision, slave, &
                            & msgtype, mpi_comm_world, status, ierr)
                        
                        mpi_simmom_matrix(:,trial) = result ! Put it in column-major order to facilitate Fortran's operation; trial is the true index of the original loop, 1, ..., trylen.
                        obj_val_vec(trial) = obj_val_1st
                        
                        ! Results That Are Collected by Individual Nodes (we are now in the my_id == 0 zone)
                        if(i==1) write(my_id+1001,'(a,(12x,a),(x,a),(2x,a),(10x,"moment1"),(10x,"moment2"),(10x,"moment3"),(10x,"moment4"),(10x,"moment5"), &
                            & (10x,"moment6"),(10x,"moment7"),(10x,"moment8"),(10x,"moment9"),(9x,"moment10"),(11x,"input1"),(11x,"input2"),(11x,"input3"), &
                            & (11x,"input4"),(11x,"input5"),(11x,"input6"),(11x,"input7"),(11x,"input8"),(11x,"input9"),(10x,"input10"))') &
                            & "MyID","error", "   #trial", "   #list"
                        
                        write(my_id+1001,'(i4,(x,f16.7),(2x,i8),(2x,i8),<ndim>(x,f16.7),<ndim>(x,f16.7))') & ! 8-2-2017. Indexeries maps the basic index number to the index number for the shifted list.
                            & slave, obj_val_1st, trial, indexseries(trial), mpi_simmom_matrix(:,trial), mpi_sobol_mixed(trial,:) ! 8-2-2017. "slave" is the node returunign the computing outcome.
                        
                        ! Check to see if one more new trial is available to be assigned to the responding slave node.
                        ! Why subtract nslaves? It's because we have in the beginning scattered nslaves trials to nslaves node.
                        if(i<=trylen-nslaves)then ! 8-2-2017 The boundary is correct. There is always nslaves trials unfinished until the end of the program. 
                            !parcel = mpi_sobol_scaled(:,i+nslaves) ! Correct. 7-3-2017, but discarded becaused we use linear combination to generate parcel. 8-2-2017.
                            trial = i + nslaves
                            !! 8-19-2017 comment out
                            !call linear_combination_sobal_sequence( parcel, trial, mpi_sobol_scaled(:,trial), origin_input, weight_list, breaks_list )   
                            !mpi_sobol_mixed(trial,:) = parcel                            
                            
                            parcel = mpi_sobol_mixed(trial,:) ! 8-19-2017 added
                            call sendjob( trial, slave, parcel ) ! 8-2-2017 Now, 'trial' is different from the value of 'slave.'
                        endif
                        exit ! leave the current loop 
                    endif 
                enddo ! (unconditional)
            enddo ! i 
            
            ! [The Root, case 1] Tell All The Slave Nodes Stopping Waiting for New Trial Assignment
            do i = 1, nslaves
                trial = -1 ! the value that triggers an exit
                slave = i
                parcel = 0._wp ! a redundant place holder for new trial
                call sendjob( trial, slave, parcel )
            enddo
            
            ! [The Root, case 1] Save the overall intput and output
            outputinput1(:,1) = real(indexseries,wp) ! Index number 
            outputinput1(:,2:ndim+1) = mpi_sobol_mixed ! Input 7-9-2017 (amoeba guess)
            outputinput1(:,ndim+2:2*ndim+1) = transpose(mpi_simmom_matrix) ! Output (moment)
            outputinput1(:,2*ndim+2) = obj_val_vec ! value of objective functin 
            call sm( outputinput1, io_string, 15, 8) ! Only the root deals with the storage of complete output.
            
        else ! [The Slaves, case 1] my_id /= 0 The Block Used for Defining the Message Passing from Slave Nodes
            do
                msgtype = 1
                ! (Integer): indexi of the passed trial
                call mpi_recv( trial, 1, mpi_integer, 0, &
                    & msgtype, mpi_comm_world, status, ierr )
                if( trial == -1 ) exit
                ! (Real): the content of the passed trial; the content is determined by the root according to the value of variable 'trial'.
                call mpi_recv( parcel, ndim, mpi_double_precision, 0, &
                    & msgtype, mpi_comm_world, status, ierr )
                
                modelmsg = 0 ! 0, model is solved successfully; 1, otherwise.
                call read_parameter_model(para,'_1parameter.txt') ! 7-9-2017
                parcel = int(parcel*accupara)/accupara ! 9-24-2017
                call search_equilibrium( parcel, result, obj_val_1st, my_id, trial, modelmsg ) ! result: simulated moments. 
                
                msgtype = 2 ! tag 2 is used for sending feedback.
                ! Sending the info about the ID of the Slave Node ('my_id').
                call mpi_send( my_id, 1, mpi_integer, 0, &
                    & msgtype, mpi_comm_world, ierr )
                ! Sending the index of the corresponding trial.
                call mpi_send( trial, 1, mpi_integer, 0, &
                    & msgtype, mpi_comm_world, ierr )
                ! Sending the simulated moments.
                call mpi_send( result, ndim, mpi_double_precision, 0, &
                    & msgtype, mpi_comm_world, ierr )
                ! Sending the level of penalty to the root.
                call mpi_send( obj_val_1st, 1, mpi_double_precision, 0, &
                    & msgtype, mpi_comm_world, ierr )
            enddo                
        endif ! [The Root and Slaves, case 1]       
        
        deallocate(parcel,result,outputinput1,obj_val_vec)
        ! # 3 move inside the MPI_exercise_mode == 1!?
        deallocate(range_guess, indexseries, sobolm, sobolm_scaled, mpi_sobol_scaled) 
        deallocate(mpi_simmom_matrix,origin_input,mpi_sobol_mixed)        
        close(my_id+1001) ! 7-3-2017
        close(my_id+2001) ! 7-4-2017
        
    elseif(mpi_exercise_mode==2)then ! Amoeba algorithm
        
        if( printout16 )then ! 8-7-2017. reset a set of variables so that only one amoeba is used and only one existing best point is read. 
            slist = 1
            elist = 1
            nslaves = 1
        endif ! printout16
        
        ! 8-2-2017: again, "nrow" is the number of nodes included in each column (equivalently), which contains the member of amoeba and idle nodes that do not participate in any amoeba operations.
        ! nrow*ncol + 1 = num_procs, where the unity on the LHS represents the system node (root).
        
        allocate( refpts(elist-slist+1), selected_input(ndim), result(ndim), bestvertex(ndim) )
        
        do i = 1, elist-slist+1 ! 7-17-2017 Generate the moving list that has the sequence of starting points.
            refpts(i) = i + (slist - 1) ! 7-22-2017 Kobol sequence index shifts ahead by slist items.
        enddo
        
        !write(*,*) ' startpoint ', startpoint ! verified.
        !if(my_id==0) write(*,*) ' refpts ', refpts ! verified.            
        
        !! 7-29-2017 ### Note that "nrow" is the vertices assigned to each amoeba, whose number is not necessary equal to (ndim+1). 
        !! 7-29-2017 ### Note that variable "nrow" is probably larger than (ndim+1) so that we have a large size of vertices under the inspection of amoeba algorithm. For example, we have a parameter 
        !! space of only 10 dimension. We ask for a specification that simultaneously 4 vertices be examined in the amoeba algorithm style. In this case, 4 times 4 nodes are required for
        !! the amoeba algorithm for each starting point ( the resulting 16 nodes includes the amoeba head with rank of 0 in respective amoeba group). Suppose in this case that we only 
        !! have only one amoeba group, and then the total nodes we require for the setup of the program is 4*4 + 1 = 17. Why acquires one more node such that the total is 17? The additional
        !! node serves as the MPI root of the whole MPI implementation. In the command line, we type, in this case, "mpiexec -np 17 MPI_Sandbox."
        
        ! 7-29-2017
        ! ncol: the number of amoeba group required
        ! nrow: the number of nodes associated with an amoeba head. In the same vein, it is not necessarily equal to (ndim+1). As explained above, a single amoeba group contains idle nodes
        ! that do not participate in a contraction operation but take part in the operation of amoeba alogrithm.
        ncol = (num_procs-1)/nrow ! 7-29-2017 variable "nslaves" below equals to "ncol" defined here. 7-13-2017 The number of amoeba heads (equivalently, the number of amoebas).
        !if( abs( (num_procs-1._wp)/nrow-floor(real((num_procs-1._wp)/nrow),wp) )>5.e-10_wp .and. my_id==0 ) write(*,*) 'The total nodes for a worst vertex should be a multiple of four plus one.' ! 8-2-2017 Not necessary the case any more.
        
        ! Initialization (7-13-2017 Not redundant; To be used in creating amoeba groups and communicators.)
        irow = 0
        icol = 0
        if( my_id/=0 )then ! 7-11-2017 Except the system root, every node gets its coordinate in the physical strucutre of the system.
            irow = merge( mod(my_id,nrow), nrow, mod(my_id,nrow)/=0 )
            icol = merge( my_id/nrow+1, my_id/nrow, mod(my_id,nrow)/=0 )
        endif
        !write(6,'("my_id: ",i3," irow: ",i3," icol: ",i3)') my_id, irow, icol
        
        !! ---- AMOEBA GROUPS (collect nodes to form each amoeba groups)
        AMOEBA_ID = -10 ! initialization. The system root is NOT in any amoeba group but has a valid amoeba_id variable becasue of this defintion.
        
        !allocate( id_list(1:nrow) ) ! Bug. 7-29-2017 MPI_GATHER
        !id_list = 1
        !do i = 2, nrow
        !    id_list(i) = id_list(i-1) + 1
        !enddo
        
        allocate( id_list(1:4*noamoeba) )
        id_list = 1
        do i = 2, 4*noamoeba
            id_list(i) = id_list(i-1) + 1
        enddo
        
        AMOEBA_WORLD = MPI_COMM_NULL ! Bug. 7-22-2017 This is needed for solving the bug below. Without it, the system root is unexpectedly included in the AMOEBA_WORLD.
        ! 7-13-2017 Trick to create a communicator by rolling over the complete list of ranks with a fixed length of rolling window.
        do i = 1, ncol ! 7-13-2017 we have "ncol" amoeba groups to create.
            ! Now, we want to let every member in the group GROUP_ALL to have the same plan about the subdivision of group GROUP_ALL. 7-10-2017
            ! link https://www.rc.usf.edu/tutorials/classes/tutorial/mpi/chapter9.html#9.2.5 [step 2]. Note: Step 1 was to assigne the same member list of COMM to a new group for use. It took place in toolbox.f90.            
            !call MPI_GROUP_INCL( GROUP_ALL, nrow, id_list, GROUP_AMOEBA, mpi_err ) ! BUG 7-29-2017
            call MPI_GROUP_INCL( GROUP_ALL, 4*noamoeba, id_list, GROUP_AMOEBA, mpi_err )
            
            !if(my_id==0) print*, mpi_err, i ! <================= checked. 
            ! Then, based on the resulting subdivision, we create a new communicator for each amoeba group, DUMMY_WORLD.
            ! link https://www.rc.usf.edu/tutorials/classes/tutorial/mpi/chapter9.html#9.5.2 [step 3]            
            call MPI_COMM_CREATE( MPI_COMM_WORLD, GROUP_AMOEBA, DUMMY_WORLD, mpi_err )
            
            !if(my_id==0) print*, mpi_err, i ! <================= checked.
            if(icol .eq. i) AMOEBA_WORLD = DUMMY_WORLD ! 7-13-2017 Creating a set of new communicators in a parallel style.
            !if(irow==1.and.icol==i) write(*,*) 'mpi_err', mpi_err, ' in amoeba group # ', icol, ' My_ID is ', My_ID ! 7-13-2017 Debug. Command: mpiexec -np 23 MPI_sandbox.exe
            
            !! Trick to create the rolling windows applied across all the nodes in order            
            id_list = id_list + nrow ! 7-13-2017 Again, nrow is the total of nodes in a single amoeba.
        enddo ! i
        deallocate( id_list )
        if( amoeba_world/=mpi_comm_null ) call MPI_COMM_RANK( AMOEBA_WORLD, AMOEBA_ID, mpi_err ) ! 7-13-2017 the rank ID of nodes in the same amoeba group.
        
        write(node_string,'(i1)') icol ! The system root returns 0 for variable "icol."
        
        if( AMOEBA_ID == 0 )then ! AMOEBA_ID=0 is also the head for contraction.
                
                solution_string = trim(node_string)//'_001_amoeba_head.txt'
                !print*, solution_string ! checked
                open(unit=my_id+1001, file=solution_string, action='write') ! AMOEBA HEAD.
                
                solution_string = trim(node_string)//'_001_macrostat.txt'
                open(unit=my_id+2001, file=solution_string, action='write', position='append')
                
                solution_string = trim(node_string)//'_001_restart_list.txt'
                open(unit=my_id+3001, file=solution_string, action='write', position='append')                
                
        endif

        !write(*,'(5i6)') my_id, irow, icol, amoeba_id, merge(-1,contract_id,CONTRACT_WORLD==MPI_COMM_NULL) ! 7-13-201 verified.
                
        !! ---- 1. Initial valuation group: contract group (for each amoeba); Product: communicator CONTRACT_WORLD
        !! The system root 0 is not in any contract groups. The system root is NOT in CONTRACT_WORLD, and hence CONTRACT_WORLD/=MPI_COMM_WORLD is a False statement for the sytstem root.
        CONTRACT_ID = -10 ! initialization
        allocate(id_list(1:ndim+1)) ! 7-13-2017 Each amoeba has (ndim+1) vertices.
        do i = 1, ndim+1
            id_list(i) = i + nrow*merge(icol-1,0,my_id/=0) ! 7-11-2017 nrow is the total nodes in a single amoeba. Although the sysmtem root has the same list as the amaoeba group of index no. 1, the following code corrects this conflict naturally.
        enddo      
    
        call MPI_GROUP_INCL( GROUP_ALL, ndim+1, id_list, GROUP_CONTRACT, mpi_err )
        call MPI_COMM_CREATE( MPI_COMM_WORLD, GROUP_CONTRACT, CONTRACT_WORLD, mpi_err )
        if( CONTRACT_WORLD/=MPI_COMM_NULL )then ! 7-13-2017 Eventhough the system root has a list, this block generates useless information for it (its communicator is "NULL" still, which can be verified by the printout a few lines below.
            call MPI_COMM_RANK( CONTRACT_WORLD, CONTRACT_ID, mpi_err )    
            call MPI_COMM_SIZE( CONTRACT_WORLD, i, mpi_err )
            !write(*,'("MY_ID: ",i3," CONTRACT_ID: ",i3," SIZE: ", i3)') MY_ID, CONTRACT_ID, i ! 7-13-2017 CONFIRMED. The system root is not included in any contraction group.
        endif        
        deallocate(id_list) 
        
        ! ALLROW_WORLD 7-30-2017
        ROW_ID = -9
        allocate( id_list(1:nrow) ) 
        do i = 1, nrow
            id_list(i) = i + nrow*merge(icol-1,0,my_id/=0)    
        enddo
        call MPI_GROUP_INCL( GROUP_ALL, nrow, id_list, GROUP_ROW, mpi_err )
        call MPI_COMM_CREATE( MPI_COMM_WORLD, GROUP_ROW, ROW_WORLD, mpi_err )
        if( ROW_WORLD/=MPI_COMM_NULL )then
            call MPI_COMM_RANK( ROW_WORLD, ROW_ID, mpi_err )
        endif
        deallocate( id_list )
            
        !write(*,'(5i6)') my_id, irow, icol, amoeba_id, merge(-1,contract_id,CONTRACT_WORLD==MPI_COMM_NULL) ! 7-13-201 verified.
                
        !! ---- COMMUNICATION GROUP (among the roots of amoeba), Product: communicator HEADWORLD
        !! The system root is in the head group, whose rank is still 0 (zero).
        HEAD_ID = -10
        allocate( id_list(0:ncol) )
        id_list = 0 ! the system root.
        do i = 1, ncol
            id_list(i) = 1 + (i-1)*nrow
        enddo ! i
        call MPI_GROUP_INCL( GROUP_ALL, ncol+1, id_list, GROUP_HEAD, mpi_err )
        call MPI_COMM_CREATE( MPI_COMM_WORLD, GROUP_HEAD, HEAD_WORLD, mpi_err )
        if( HEAD_WORLD/=MPI_COMM_NULL)then
                call MPI_COMM_RANK( HEAD_WORLD, HEAD_ID, mpi_err )
                !write(*,'("MY_ID: ",i3," HEAD_ID: ",i3," ncol: ",i3)') MY_ID, HEAD_ID, ncol ! 7-17-2017 checked.
        endif 
        deallocate(id_list)

        ! Indexing within respective amoebas
        sirow = mod( AMOEBA_ID, 4 ) + 1 ! 7-17-2017 The order is [2,3,4,1]. The ID "within" a subamoeba group.
        sicol = AMOEBA_ID/4 + 1 ! 7-17-2017 The ID of a "subamoeba group" within an amoeba.
        nslaves = (num_procs-1)/nrow ! 10.14.2017. 7-29-2017 This variable equals to variable "ncol" defined above. 7-17-2017 Number of "head nodes" except for the system root; equivalently, the number of amoebas.
        
        !! 7-29-2017 checked!! 7-22-2017 **Negative number indicates the node is not in a specific group.**
        if(printout13) write(*,'("my_id ",i3," irow ",i3," icol ",i3," amoeba ",i3," contract ",i3," head ",i3, " sicol ", i3, " sirow ", i3, " row_id ", i3)') my_id, irow, icol, merge(-9,amoeba_id,AMOEBA_WORLD==MPI_COMM_NULL), merge(-9,contract_id,CONTRACT_WORLD==MPI_COMM_NULL), merge(-9,head_id,HEAD_WORLD==MPI_COMM_NULL), merge(-9,sicol,AMOEBA_WORLD==MPI_COMM_NULL), merge(-9,sirow,AMOEBA_WORLD==MPI_COMM_NULL), merge(-9,ROW_ID,ROW_WORLD==MPI_COMM_NULL) ! 7-22-2017 
        
        ! create the output filename for each amoeba member
        write(node_string,'(i1)') icol ! the index of an amoeba group.
        amoeba_x_y_string = trim(node_string)//"_"
        !write(node_string,'(i3)') AMOEBA_ID ! BUG. 7-29-2017 
        write(node_string,'(i3.3)') irow ! 7-29-2017 Neither amoeba_id nor contract_id is suitable for the purpose of creating a filename for each node contained in a single amoeba group, because the total vertices in a single amoeba may exceeds the number of a contract group (dim+1) and the number of amoeba-algorithm-associated vertices (4*noamoeba).
        amoeba_x_y_string = trim(amoeba_x_y_string)//trim(node_string)
        
        if( CONTRACT_WORLD/=MPI_COMM_NULL .and. AMOEBA_ID /= 0 )then ! 7-22-201 This excludes the system root and the amoeba heads.
            solution_string = trim(amoeba_x_y_string)//'_contract_member.txt' ! For example, 1_i_amoeba_member.txt, where i = 1,...,10 (except 0 the amoeba head).
            open(unit=my_id+1001, file=solution_string, action='write') ! AMOEBA MEMBERS.
        elseif( CONTRACT_WORLD==MPI_COMM_NULL .and. my_id/=0 )then
            solution_string = trim(amoeba_x_y_string)//'_Non-contract_member.txt'
            open(unit=my_id+1001, file=solution_string, action='write')
        endif        
        
        !Core of communication among mpi_comm_world ======================================================================================================================        
        if( my_id == 0 )then ! 7-27-2017 Checked from ln.389 to ln.482. 
            
            write(trylen_string,'(i5.5,"_",i5.5)') slist, elist
            solution_string = '0_000_Answers'//trim(trylen_string)//'.txt' ! THE SYSTEM ROOT.
            open(unit=my_id+1001, file=solution_string, action='write')
            
            ! Initialization. ! Tonight. Debug starts here.
            do i = 1, nslaves ! Across all the amoeba heads. 7-30-2017 Again, nslaves == ncol. 
                trial = refpts(i) ! 7-17-2017 Assign a trial ID number from the reference list, refpts.
                !print*, 'My_Id=', my_id, ' slave=', i, ' trial=', trial, ' nslaves: ', nslaves ! 7-17-2017 checked. 
                slave = i ! 7-17-2017 Correct. We use the rank ID in the Head_World group.

                if(printout16)then     
                    call read_best_point(selected_input,bestvertex_file)
                else ! This is for the first round after the coarse search (stage 1).
                    call read_one_point_from_short_list(trial,startpoint) ! 7-17-2017 Get a new starting point based on the input, trial.
                    selected_input = startpoint(15:24)
                endif
                
                ! The very first communication
                call sendjob(trial, slave, selected_input, opt='head_group') ! 7-17-2017 trial is a new starting point's index no. from the integrated list of favorable starting points.
                
                ! Only my_id == 0 would prints out.
                if(printout13) write(my_id+1001,'("Trial no. (refpts): ", i8, ", VERTEX ", <ndim>f12.7, " amoeba no.", i3, " (Total amoebas: ", i3,")")') trial, selected_input, slave, nslaves ! 
                if(printout13) write(*,'(/,"Trial no. (refpts): ", i8, ", VERTEX ", <ndim>f12.7, " amoeba no.", i3, " (Total amoebas: ", i3,")")') trial, selected_input, slave, nslaves ! 
                
            enddo ! i
            
            ! Listen for responses
            !do i = 1, elist-slist+1
            index_pt = nslaves ! 7-26-2017 Bug (ndim->nslaves). 7-23-2017 We have use the first ndim starting points, so we have to now start with index ndim.
            index_rcv = 0 ! 7-23-2017 The counter for book keeping the results already received.
            
            ! communication of amoeba heads.
            !do while(index_rcv>elist-slist+1) ! Bug.
            
            do
                msgtype = 2
                
                call MPI_IPROBE( MPI_ANY_SOURCE, msgtype, HEAD_WORLD, & ! 7-17-2017 Trick for the system root to scan incoming feedback from any ameoba heads.
                    & receiving, status, ierr )
                
                if( receiving )then ! check to see if the root still receives the output from the collaborated amoeba heads. 
                    
                    call MPI_RECV( slave, 1, MPI_INTEGER, MPI_ANY_SOURCE, & ! slave = the id of amoeba head.
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    call MPI_RECV( trial, 1, MPI_INTEGER, slave, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    call MPI_RECV( result, ndim, MPI_DOUBLE_PRECISION, slave, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    !print*, " ierr 1 ", ierr
                    
                    call MPI_RECV( bestvertex, ndim, MPI_DOUBLE_PRECISION, slave, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    !print*, " ierr 2 ", ierr
                    
                    call MPI_RECV( bestobjval, 1, MPI_DOUBLE_PRECISION, slave, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    !print*, " ierr 3 ", ierr
                    
                    ! bookkeeping output from each amoeba group.                    
                    write(my_id+1001,'((4x,a), &
                        & (11x,"input1"),(11x,"input2"),(11x,"input3"),(11x,"input4"),(11x,"input5"), &
                        & (11x,"input6"),(11x,"input7"),(11x,"input8"),(11x,"input9"),(10x,"input10"), &
                        & (10x,"moment1"),(10x,"moment2"),(10x,"moment3"),(10x,"moment4"),(10x,"moment5"), &
                        & (10x,"moment6"),(10x,"moment7"),(10x,"moment8"),(10x,"moment9"),(9x,"moment10"))') &
                        & "func_val"
                    
                    write(my_id+1001,'(f16.7, <ndim>(x,f16.7),<ndim>(x,f16.7))') bestobjval, bestvertex, result
                    
                    index_pt  = index_pt + 1
                    index_rcv = index_rcv + 1
                    write(my_id+1001,'(a,i3,a,i3,a,i3,/)') "-index_rcv: ", index_rcv, ",index_pt: ", index_pt, ",elist-slist+1: ", elist-slist+1
                    
                    ! Send out new starting points as long as the list does not run out.
                    if( index_pt<=elist-slist+1 )then
                        trial = refpts(index_pt)
                        call read_one_point_from_short_list(trial,startpoint)
                        selected_input = startpoint(15:24)
                        call sendjob(trial, slave, selected_input, opt='head_group')
                        if(printout13) write(my_id+1001,'(a)') '-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
                        if(printout13) write(my_id+1001,'("Trial no. (refpts): ", i8, ", VERTEX ", <ndim>f12.7, " amoeba no.", i3, " (Total amoebas: ", i3,")")') trial, selected_input, slave, nslaves ! 
                        if(printout13) write(*,'("Trial no. (refpts): ", i8, ", VERTEX ", <ndim>f12.7, " amoeba no.", i3, " (Total amoebas: ", i3,")")') trial, selected_input, slave, nslaves ! 
                    endif
                    
                endif ! receiving, if statement
                
                !write(my_id+1001,'(a,i3)') " ready to exit? ", index_rcv
                
                if( index_rcv>=elist-slist+1 ) exit ! 7-27-2017 Bug. ">"-->">=". 7-18-2017 we have collect the outcome for the full list of starting points.
                
            enddo ! indefinite do loop
            !enddo ! i
            
            ! Calling all heads of amoeba to stop subsequent evaluation. 7-18-2017
            do i = 1, nslaves
                trial = -1 ! 7-18-2017 Hard-wired for terminating the amoeba calculation.
                slave = i
                selected_input = -1._wp
                !write(my_id+1001,'(a,i3)') " stop signal ", nslaves ! 7-27-2017 checked.
                call sendjob(trial, slave, selected_input, opt='head_group') ! 7-18-2017 checked. The program terminates properly.
            enddo
            
        else ! my_id != 0
            if( HEAD_WORLD/=MPI_COMM_NULL )then ! THE AMOEBA HEAD BLOCK (communicate only with the heads of amoeba, "including" the system root.) 7-11-2017
                do 
                    msgtype = 1
                    
                    call MPI_RECV( trial, 1, MPI_INTEGER, 0, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    if( trial == -1 )then ! 7-18-2017 Let all the amoeba members know stopping any further calculation.
                        call MPI_BCAST( trial, 1, MPI_INTEGER, 0, ROW_WORLD, ierr ) ! 7-30-2017 AMOEBA_WORLD-->ROW_ID   
                        exit ! 7-18-2017 Important.
                    endif
                    
                    call MPI_RECV( selected_input, ndim, MPI_DOUBLE_PRECISION, 0, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    if(printout13)then
                        write(my_id+1001,'(a,x,i6)') '[____] Starting Point Trial:', trial
                    
                        do i = 1, ndim-1
                            write(my_id+1001, '((12x,"Input",x,i2),2x)', advance='no') i            
                        enddo ! i
                        write(my_id+1001, '((12x,"Input",x,i2),2x)') ndim
                        
                        do i = 1, ndim-1
                            write(my_id+1001, '(f20.15,2x)', advance='no') selected_input(i)   
                        enddo
                        write(my_id+1001, '(f20.15,2x)') selected_input(ndim)
                        
                        !write(my_id+1001,'(a,x,i6, ",", " vertex ",<ndim>f12.7)') 'Trial:', trial, selected_input ! 7-23-2017 checked.
                    endif
                        
                    !call MPI_BCAST( trial, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )
                    !call MPI_BCAST( selected_input, ndim, MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err ) ! 7-23-2017 maks sense.
                    call MPI_BCAST( trial, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err ) ! 7-30-2017 CONTRACT_WORLD-->ROW_WORLD
                    call MPI_BCAST( selected_input, ndim, MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err ) ! 7-30-2017 CONTRACT_WORLD-->ROW_WORLD                    
                    
                    call amoeba_algorithm( trial, selected_input, bestvertex, result, bestobjval ) ! 7-18-2017 The major part of the amoeba algorithm. Be sure to refresh the model parameter input. [Li-Pin Juan no.2]
                    
                    !! 7-19-2017 Send The local optima back to the system root from the amoeba heads.
                    msgtype = 2 ! Tag 2 is used for pass-receiving results from respective amoeba head.
                    
                    call MPI_SEND( HEAD_ID, 1, MPI_INTEGER, 0, & ! 7-18-2017 [1] Inform all its amoeba members what their amoeba ID is.
                        & msgtype, HEAD_WORLD, ierr )
                    
                    call MPI_SEND( trial, 1, MPI_INTEGER, 0, & ! [2] trial number
                        & msgtype, HEAD_WORLD, ierr ) 
                    
                    call MPI_SEND( result, ndim, MPI_DOUBLE_PRECISION, 0, & ! [3] best moments
                        & msgtype, HEAD_WORLD, ierr )
                    
                    call MPI_SEND( bestvertex, ndim, MPI_DOUBLE_PRECISION, 0, & ! [4] bets inputs
                        & msgtype, HEAD_WORLD, ierr )
                    
                    call MPI_SEND( bestobjval, 1, MPI_DOUBLE_PRECISION, 0, & ! [5] best objective vale
                        & msgtype, HEAD_WORLD, ierr )
                    
                    write(my_id+1001,'(/,a,i8.8," head_id: ",i3.3," min_val: ",f12.7)') "Outcome for trial no.", trial, head_id, bestobjval
                    write(my_id+1001,'(a,/)') '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                enddo
                
            elseif( ROW_WORLD/=MPI_COMM_NULL )then ! 7-18-2017 The amoeba members (i.e. excluding the amoeba head) in respective amoeba group.    
                
                ! Note that we prevent amoeba members except their head to communicate with the system root.
                ! By doing so, the starting point assignement is carried out only between the system root and the head of each amoeba group.
                do 
                    !call MPI_BCAST( trial, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err ) ! keeps paying attention to the announcement of the amoeba head. ! 7-29-2017 BUG   
                    call MPI_BCAST( trial, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err ) ! 7-30-2017 CONTRACT_WORLD --> ROW_WORLD
                    if( trial == -1) exit
                    !call MPI_BCAST( selected_input, ndim, MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err ) ! 7-29-2017 BUG
                    call MPI_BCAST( selected_input, ndim, MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err ) ! 7-30-2017 CONTRACT_WORLD --> ROW_WORLD
                    call amoeba_algorithm( trial, selected_input, bestvertex, result, bestobjval ) ! # Need imporvement. [Li-Pin Juan no.2]
                enddo
            endif ! HEAD_WORLD/=MPI_COMM_NULL (ROW_WORLD/=MPI_COMM_NULL )
            
        endif ! my_id == 0 or not.
            
        !! Important!!
        !if(irow==1) call read_short_list(listnumber,startpoint)
        !refpts = startpoint(15:24)
        !print*, ' debug ', ncol, num_procs, nrow ! successful command: mpiexec -np 12 (or 23) MPI_sandbox.exe (for the case of nrow == 11; 10 dimension and one more node for head.
        
        deallocate( refpts, selected_input, result, bestvertex )
        
        ! 8-2-2017 Close files created in lines 421-434.
        if( ROW_WORLD/=MPI_COMM_NULL .or. my_id == 0 )then
            close(my_id+1001)
            if( AMOEBA_ID == 0 )then
                close(my_id+2001)    
                close(my_id+3001)
            endif
        endif
        
    elseif( mpi_exercise_mode==3 )then ! Should mimic what mpi_exercise_mode==0 and 2 did. 8-7-2017.   
        
        solution_string = 'Single_node_results.txt'   
        !concisesolution_string = 'SingleNodeDetails.txt'
        !open(unit=my_id+1001, file=solution_string, action='write', position='append')
        open(unit=my_id+1001, file=solution_string, action='write', status="replace")
        !open(unit=my_id+2001, file=concisesolution_string, action='write', position='append')        
        
        i = 100
        call read_parameter_model(para, '_1parameter.txt')
        guessv(1) = kv1   
        guessv(2) = prtk0 
        guessv(3) = prtk1 
        guessv(4) = prtk2 
        guessv(5) = zbar  
        guessv(6) = beta  
        guessv(7) = theta 
        !guessv(7) = iota
        guessv(8) = phi1  
        guessv(9) = phi2  
        guessv(10)= phi3 
        modelmsg = 0
        momvec = inf
        obj_val_1st = inf
        write(msg,fmt='(i3)') i ! either 100 or 200
        write(my_id+1001,'(a,i3,a,/)') " ================================ ", i, " ================================="
        open(unit=4000+i,file="output_parameter_inspection_"//trim(msg)//".txt",action="write",status="replace")
        call search_equilibrium( guessv, momvec, obj_val_1st, 100, i, modelmsg  )
        if( modelmsg == 0 )then
            write(4000+i, '(a,<ndim>f15.7)') '1guess  : ', guessv
            write(4000+i, '(a,<ndim>f15.7)') '1targetv: ', targetv
            write(4000+i, '(a,<ndim>f15.7)') '1moment : ', momvec
            write(4000+i, '(a,f15.7)') '1penalty: ', obj_val_1st             
        else
            write(4000+i, '(a,a,<ndim>f15.7)') ' === Failure === ', 'guess  : ', guessv
        endif ! modelmsg
        close(4000+i) ! test
        
        i = 200
        call read_parameter_model(para, '_1parameter.txt')
        guessv(1) = kv1   
        guessv(2) = prtk0 
        guessv(3) = prtk1 
        guessv(4) = prtk2 
        guessv(5) = zbar  
        guessv(6) = beta  
        guessv(7) = theta 
        !guessv(7) = iota
        guessv(8) = phi1  
        guessv(9) = phi2  
        guessv(10)= phi3 
        modelmsg = 0
        momvec = inf
        obj_val_1st = inf
        write(msg,fmt='(i3)') i ! either 100 or 200
        write(my_id+1001,'(/,a,i3,a,/)') " ================================ ", i, " ================================="
        open(unit=4000+i,file="output_parameter_inspection_"//trim(msg)//".txt",action="write",status="replace")
        call search_equilibrium( guessv, momvec, obj_val_1st, 100, i, modelmsg  )
        if( modelmsg == 0 )then
            write(4000+i, '(a,<ndim>f15.7)') '2guess  : ', guessv
            write(4000+i, '(a,<ndim>f15.7)') '2targetv: ', targetv
            write(4000+i, '(a,<ndim>f15.7)') '2moment : ', momvec
            write(4000+i, '(a,f15.7)') '2penalty: ', obj_val_1st             
        else
            write(4000+i, '(a,a,<ndim>f15.7)') ' === Failure === ', 'guess  : ', guessv
        endif ! modelmsg     
        close(4000+i)
        close(my_id+1001)  
        !close(my_id+2001)    
        
    elseif(mpi_exercise_mode==4)then ! 9-30-2017
        ! 10-1-2017
        allocate( mat_stage1_inputs(mode4list,11) ) !10.21.2017 Change it!!!
        trylen = size(mat_stage1_inputs,dim=1)
        
        if(my_id==0) call read_matrix(mat_stage1_inputs,'_stage1_input.csv')
        !do i = 1, size(mat_stage1_inputs,dim=1)
        !    write(*,'(11(f8.3,x))') mat_stage1_inputs(i,:)
        !enddo ! i
        
        write(node_string,'(i3.3)') my_id
        
        if(my_id/=0)then
            solution_string = 'SlaveIntermFeedback_'//trim(node_string)//'.txt' ! 7-7-2017 Intermediate outcome of individual slave node.         
            concisesolution_string = 'SlaveMacrostat_'//trim(node_string)//'.txt' ! 7-7-2017 Macro stat of individual slave node.        
        else
            write(trylen_string,'(i5.5,"_",i5.5)') sblno1, sblno1+trylen-1
            solution_string = 'RootIntermCollect_'//trim(trylen_string)//'.txt' ! 7-7-2017 Real-time feedback from the slaves.        
            concisesolution_string = 'RootUseless_'//trim(trylen_string)//'.txt' ! 7-7-2017 Useless for the root node.        
            ! add concise solution
        endif
        
        open(unit=my_id+1001, file=solution_string, action='write') ! Moved here. 7-3-201
        open(unit=my_id+2001, file=concisesolution_string, action='write') ! Moved here. 7-5-201
        
        ! # 1 move inside the MPI_exercise_mode == 1!?
        allocate( indexseries(trylen), sobolm(nsbq, ndim), sobolm_scaled(ndim,nsbq), mpi_sobol_scaled(ndim,trylen) )
        allocate( mpi_simmom_matrix(ndim,trylen), origin_input(ndim), mpi_sobol_mixed(trylen,ndim) ) ! 10-1-2017
        !allocate(indexseries(trylen),sobolm(ndim,nsbq))
        
        if(my_id==0)then
            if(trylen<size(mat_stage1_inputs,dim=1))then ! 10-1-2017
                write(*,'(/,a,/)') "Must to set trylen bigger than the size of mat_stage1_inputs!!! "
            else
                write(*,'(/,a,i6,a,i6)') "trylen: ", trylen, " length of mat_stage1_inputs: ", size(mat_stage1_inputs,dim=1) 
            endif ! 10-1-2017
        endif ! my_id
        
        indexseries = [(i,i=1,trylen)]
        !if(my_id==0) write(*,'(i5)') (indexseries(i),i=1,trylen) ! works. 
        !if(my_id==0) write(*,'(a,/)') ' '
        
        indexseries = indexseries + sblno1 -1 ! Shifted; USED FOR OUTPUT    
        !if(my_id==0) write(*,'(i5)') (indexseries(i),i=1,trylen) 
        
        origin_input(1) = kv1   
        origin_input(2) = prtk0 
        origin_input(3) = prtk1 
        origin_input(4) = prtk2 
        origin_input(5) = zbar  
        origin_input(6) = beta  
        origin_input(7) = theta 
        !origin_input(7) = iota
        origin_input(8) = phi1  
        origin_input(9) = phi2  
        origin_input(10)= phi3           
        
        ! Quasi-random Sobol sequence block # 2 move inside the MPI_exercise_mode == 1!?
        
        call get_sobol_sequence( sobolm, 0.0_wp, 1.0_wp, nsbq, ndim ) ! Generate ndim dimensional sobol sequence of length nsbq (nsbq*ndim).
        ! call scale_sobol_original( transpose(sobolm), range_guess, sobolm_scaled ) ! ndim*nsbq ! comment out 8-18-2017
        
        sobolm_scaled = transpose(sobolm) ! added 8-18-2017
        mpi_sobol_scaled = sobolm_scaled(:,sblno1:sblno1+trylen-1) ! 7-29-2017 Shifted.    
        
        !if(my_id==0) call sm(sobolm,'sobolm'//trim(trylen_string)) ! checked 2017-Jul-1
        !if(my_id==0) call sm(transpose(mpi_sobol_scaled),'mpi_sobol_scaled'//trim(trylen_string)) ! checked 2017-Jul-1 . Used for convex combination.        
        
        allocate( parcel(ndim), result(ndim), outputinput1(trylen,2*ndim+2), obj_val_vec(trylen) )        
        nslaves = num_procs - 1 ! Just for mpi_exercise_mode==1 case, the root processor (my_id==0) is in charge of sending out new trial and receiving the corresponding result.
        
        if( my_id == 0)then
            
            write(trylen_string,'(i5.5,"_",i5.5)') sblno1, sblno1+trylen-1
            
            !! 8-18-2017
            io_string = 'FinalScaledSobel_'//trim(trylen_string) 
            do i = 1, trylen
                call linear_combination_sobol_sequence_corrected( parcel, i, mpi_sobol_scaled(:,i), range_guess )
                mpi_sobol_mixed(i,:) = parcel
            enddo
            call sm(mpi_sobol_mixed,io_string)
            
            io_string = 'IOMat_'//trim(trylen_string) 
            ! [The Root, case 1] Send Initial Messages to Slave Nodes (indices ranges from 1 to nslaves)
            do i = 1, nslaves ! nslaves = num_procs-1.
                trial = i ! the index of the basic loop: 1,...,trylen; Not the indices of the adjusted Sobol sequence.
                slave = i
                
                !! comment out 8-18-2017
                !call linear_combination_sobal_sequence(parcel,trial,mpi_sobol_scaled(:,trial),origin_input,weight_list,breaks_list) ! Be sure to set up weight_list and breaks_list.    
                
                !! added 8-18-2017, then moved outside this do loop over i between 1 and nslaves
                !call linear_combination_sobol_sequence_corrected( parcel, trial, mpi_sobol_scaled(:,trial), range_guess )
                !! The formula is: vertex_output = (range_guess(:,2)+range_guess(:,1))/2._wp + (unit_sobol_input-0.5_wp)*(range_guess(:,2)-range_guess(:,1))
                !mpi_sobol_mixed(trial,:) = parcel ! bookkeeping. ! comment out 8-18-2017
                
                !parcel = mpi_sobol_mixed(i,:)
                parcel = mat_stage1_inputs(i,2:11) ! 10-1-2017
                
                
                if(printout21)then !10.17.2017
                    call sendjob(trial,slave,parcel,srnumber=int(mat_stage1_inputs(i,1)))
                else
                    call sendjob(trial,slave,parcel) ! send from root to slave the trial parameter combination                
                endif
            enddo
            
            ! [The Root, case 1] Hear Responses from Slave Nodes
            do i = 1, trylen ! The range is correct!! 7-3-2017 ! I controls the total number of outcome expected to receive. 8-2-2017
                do
                    msgtype = 2 ! tag 1 for parameter passing; tag 2 for communicating the result between slave nodes and root. 
                    ! Non-blocking test for a message
                    call mpi_iprobe( mpi_any_source, msgtype, mpi_comm_world, &
                        & receiving, status, ierr)
                    if(receiving)then
                        ! 1. [slave] Which member node is sending the result? 
                        call mpi_recv( slave, 1, mpi_integer, mpi_any_source, &
                            & msgtype, mpi_comm_world, status, ierr)
                        ! 2. [trial] what's the returned trial index in the basic assignment loop?
                        call mpi_recv( trial, 1, mpi_integer, slave, & ! Note: trial falls in [1,trylen], rather than the shifted interval.
                            & msgtype, mpi_comm_world, status, ierr)
                        ! 3. [result] the feedback of simulated moments 
                        call mpi_recv( result, ndim, mpi_double_precision, slave, &
                            & msgtype, mpi_comm_world, status, ierr)
                        ! 4. [obj_val_1st] the value of penalty corresponds to the given trial
                        call mpi_recv( obj_val_1st, 1, mpi_double_precision, slave, &
                            & msgtype, mpi_comm_world, status, ierr)
                        
                        if(printout21)then !10.17.2017
                            call mpi_recv( srnumber, 1, mpi_integer, slave, & ! 10.17.2017
                                & msgtype, mpi_comm_world, status, ierr)    
                        endif
                        
                        mpi_simmom_matrix(:,trial) = result ! Put it in column-major order to facilitate Fortran's operation; trial is the true index of the original loop, 1, ..., trylen.
                        obj_val_vec(trial) = obj_val_1st
                        
                        ! Results That Are Collected by Individual Nodes (we are now in the my_id == 0 zone)
                        if(i==1) write(my_id+1001,'(a,(12x,a),(x,a),(2x,a),(10x,"moment1"),(10x,"moment2"),(10x,"moment3"),(10x,"moment4"),(10x,"moment5"), &
                            & (10x,"moment6"),(10x,"moment7"),(10x,"moment8"),(10x,"moment9"),(9x,"moment10"),(11x,"input1"),(11x,"input2"),(11x,"input3"), &
                            & (11x,"input4"),(11x,"input5"),(11x,"input6"),(11x,"input7"),(11x,"input8"),(11x,"input9"),(10x,"input10"))') &
                            & "MyID","error", "   #trial", "   #list"
                        
                        !write(my_id+1001,'(i4,(x,f16.7),(2x,i8),(2x,i8),<ndim>(x,f16.7),<ndim>(x,f16.7))') & ! 8-2-2017. Indexeries maps the basic index number to the index number for the shifted list.
                        !    & slave, obj_val_1st, trial, indexseries(trial), mpi_simmom_matrix(:,trial), mpi_sobol_mixed(trial,:) ! 8-2-2017. "slave" is the node returunign the computing outcome.
                        
                        if(printout21)then !10.17.2017
                            write(my_id+1001,'(i4,(x,f16.7),(2x,i8),(2x,i8),<ndim>(x,f16.7),<ndim>(x,f16.7))') & ! 10-1-2017
                                & slave, obj_val_1st, srnumber, indexseries(trial), mpi_simmom_matrix(:,trial), mat_stage1_inputs(trial,2:11)                                                 
                        else
                            write(my_id+1001,'(i4,(x,f16.7),(2x,i8),(2x,i8),<ndim>(x,f16.7),<ndim>(x,f16.7))') & ! 10-1-2017
                                & slave, obj_val_1st, trial, indexseries(trial), mpi_simmom_matrix(:,trial), mat_stage1_inputs(trial,2:11)                         
                        endif
                        
                        ! Check to see if one more new trial is available to be assigned to the responding slave node.
                        ! Why subtract nslaves? It's because we have in the beginning scattered nslaves trials to nslaves node.
                        if(i<=trylen-nslaves)then ! 8-2-2017 The boundary is correct. There is always nslaves trials unfinished until the end of the program. 
                            !parcel = mpi_sobol_scaled(:,i+nslaves) ! Correct. 7-3-2017, but discarded becaused we use linear combination to generate parcel. 8-2-2017.
                            trial = i + nslaves
                            !! 8-19-2017 comment out
                            !call linear_combination_sobal_sequence( parcel, trial, mpi_sobol_scaled(:,trial), origin_input, weight_list, breaks_list )   
                            !mpi_sobol_mixed(trial,:) = parcel                            
                            
                            !parcel = mpi_sobol_mixed(trial,:) ! 8-19-2017 added
                            parcel = mat_stage1_inputs(trial,2:11) ! 10-1-2017
                            
                            if(printout21)then !10.17.2017
                                call sendjob( trial, slave, parcel, srnumber=int(mat_stage1_inputs(trial,1)))
                            else    
                                call sendjob( trial, slave, parcel ) ! 8-2-2017 Now, 'trial' is different from the value of 'slave.'
                            endif
                        endif
                        exit ! leave the current loop 
                    endif 
                enddo ! (unconditional)
            enddo ! i 
            
            ! [The Root, case 1] Tell All The Slave Nodes Stopping Waiting for New Trial Assignment
            do i = 1, nslaves
                trial = -1 ! the value that triggers an exit
                slave = i
                parcel = 0._wp ! a redundant place holder for new trial
                if(printout21)then
                    call sendjob( trial, slave, parcel, srnumber=trial ) ! the 4th argument is redundant 10.17.2017
                else
                    call sendjob( trial, slave, parcel )
                endif
            enddo
            
            !! [The Root, case 1] Save the overall intput and output
            !outputinput1(:,1) = real(indexseries,wp) ! Index number 
            !outputinput1(:,2:ndim+1) = mpi_sobol_mixed ! Input 7-9-2017 (amoeba guess)
            !outputinput1(:,ndim+2:2*ndim+1) = transpose(mpi_simmom_matrix) ! Output (moment)
            !outputinput1(:,2*ndim+2) = obj_val_vec ! value of objective functin 
            !call sm( outputinput1, io_string, 15, 8) ! Only the root deals with the storage of complete output.
            
        else ! [The Slaves, case 1] my_id /= 0 The Block Used for Defining the Message Passing from Slave Nodes
            do
                msgtype = 1
                ! (Integer): indexi of the passed trial
                call mpi_recv( trial, 1, mpi_integer, 0, &
                    & msgtype, mpi_comm_world, status, ierr )
                if( trial == -1 ) exit
                ! (Real): the content of the passed trial; the content is determined by the root according to the value of variable 'trial'.
                call mpi_recv( parcel, ndim, mpi_double_precision, 0, &
                    & msgtype, mpi_comm_world, status, ierr )
                if(printout21)then !10.17.2017
                    call mpi_recv( srnumber, 1, mpi_integer, 0, & ! 10.17.2017
                        & msgtype, mpi_comm_world, status, ierr)    
                endif                
                
                modelmsg = 0 ! 0, model is solved successfully; 1, otherwise.
                call read_parameter_model(para,'_1parameter.txt') ! 7-9-2017
                parcel = int(parcel*accupara)/accupara ! 9-24-2017
                call search_equilibrium( parcel, result, obj_val_1st, my_id, trial, modelmsg ) ! result: simulated moments. 
                
                msgtype = 2 ! tag 2 is used for sending feedback.
                ! Sending the info about the ID of the Slave Node ('my_id').
                call mpi_send( my_id, 1, mpi_integer, 0, &
                    & msgtype, mpi_comm_world, ierr )
                ! Sending the index of the corresponding trial.
                call mpi_send( trial, 1, mpi_integer, 0, &
                    & msgtype, mpi_comm_world, ierr )
                ! Sending the simulated moments.
                call mpi_send( result, ndim, mpi_double_precision, 0, &
                    & msgtype, mpi_comm_world, ierr )
                ! Sending the level of penalty to the root.
                call mpi_send( obj_val_1st, 1, mpi_double_precision, 0, &
                    & msgtype, mpi_comm_world, ierr )
                if(printout21)then !10.17.2017
                    call mpi_send( srnumber, 1, mpi_integer, 0, & ! 10.17.2017
                        & msgtype, mpi_comm_world, ierr)    
                endif                 
                
            enddo                
        endif ! [The Root and Slaves, case 1]       
        
        deallocate(parcel,result,outputinput1,obj_val_vec)
        !! # 3 move inside the MPI_exercise_mode == 1!?
        !deallocate(range_guess, indexseries, sobolm, sobolm_scaled, mpi_sobol_scaled) 
        !deallocate(mpi_simmom_matrix,origin_input,mpi_sobol_mixed)        
        !close(my_id+1001) ! 7-3-2017
        !close(my_id+2001) ! 7-4-2017        
        
        !ncol = num_procs/nrow ! inconsistent. 10.19.2017
        
        deallocate( mat_stage1_inputs )
        
    elseif(mpi_exercise_mode==5)then !10.15.2017 Task: dimensional reduction
        
        allocate(pointlist(listlength,25), pts_ndim(subdim+1,ndim), pts_subdim(subdim+1,subdim))
        
        !do i = 1, size(pts_ndim,dim=1)
        !    write(*,'(10(f8.5,x))')  pts_ndim(i,:)
        !enddo !i
        ncol = num_procs/nrow ! Should be one anyway. Note that in this case the root (0) is a member of the single amoeba group.
        irow = my_id+1
        icol = ncol
        
        write(*,'(5(a,i3,x))'), "subdim: ", subdim, " ndim: ", ndim, " ncol: ", ncol, ' irow: ', irow, ' icol: ', icol
        
        ! The amoeba algorithm membership #1
        amoeba_id = -10
        allocate(id_list(1:4*noamoeba))
        id_list = 0
        do i = 2, 4*noamoeba
            id_list(i) = id_list(i-1) + 1
        enddo ! i
        amoeba_world = mpi_comm_null
        do i = 1, ncol
            call mpi_group_incl(group_all, 4*noamoeba, id_list, group_amoeba, mpi_err)
            call mpi_comm_create(mpi_comm_world, group_amoeba, dummy_world, mpi_err)
            if(icol.eq.i) amoeba_world = dummy_world
            id_list = id_list + nrow
        enddo ! i
        deallocate(id_list)
        if(amoeba_world/=mpi_comm_null) call mpi_comm_rank(amoeba_world, amoeba_id, mpi_err)
        if(amoeba_id==0)then
            write(node_string, '(i1)') icol
            solution_string = trim(node_string)//'id001_fid1001_intermediate.txt'
            open(unit=my_id+1001, file=solution_string, action='write')
            solution_string = trim(node_string)//'id001_fid2001_macro_stats.txt'
            open(unit=my_id+2001, file=solution_string, action='write')
            solution_string = trim(node_string)//'id001_fid3001_restart_pts.txt'
            open(unit=my_id+3001, file=solution_string, action='write')            
        endif ! amoeba_id==0
        
        ! Membership for vertices in charge of contraction #2
        contract_id = -10
        allocate(id_list(1:subdim+1))
        id_list = 0
        do i = 2, subdim+1
            id_list(i) = id_list(i-1) + 1
        enddo ! i
        call mpi_group_incl(group_all, subdim+1, id_list, group_contract, mpi_err)
        call mpi_comm_create(mpi_comm_world, group_contract, contract_world, mpi_err)
        if(contract_world/=mpi_comm_null)then
            call mpi_comm_rank(contract_world, contract_id, mpi_err)
            call mpi_comm_size(contract_world, i, mpi_err)
        endif
        !if(my_id==0) write(*,*) " size of contract world: ", i
        deallocate(id_list)
        
        ! The allrow_world #3
        row_id = -9
        allocate(id_list(1:nrow))
        id_list = 0
        do i = 2, nrow
            id_list(i) = id_list(i-1) + 1    
        enddo! i
        call mpi_group_incl(group_all, nrow, id_list, group_row, mpi_err)
        call mpi_comm_create(mpi_comm_world, group_row, row_world, mpi_err)
        if(row_world/=mpi_comm_null) call mpi_comm_rank(row_world, row_id, mpi_err)
        deallocate(id_list)
        !write(*,'(5(a,i3,x))') "myid: ", my_id, " irow: ", irow, " icol: ", icol, " amoeba_id: ", amoeba_id, " row_id: ", row_id
        
        !Note the sirow and sicol of vertices outside the amoeba_world are meanless. Be sure not to use these variables.
        sirow   = mod(amoeba_id, 4) + 1
        sicol   = amoeba_id/4 + 1
        nslaves = num_procs/nrow
        
        if(printout13) write(*,'("my_id ",i3," irow ",i3," icol ",i3," amoeba ",i3," contract ",i3," sicol ", i3, " sirow ", i3, " row_id ", i3)') &
            my_id, irow, icol, merge(-9,amoeba_id,AMOEBA_WORLD==MPI_COMM_NULL), merge(-9,contract_id,CONTRACT_WORLD==MPI_COMM_NULL), &
            merge(-9,sicol,AMOEBA_WORLD==MPI_COMM_NULL), merge(-9,sirow,AMOEBA_WORLD==MPI_COMM_NULL), merge(-9,ROW_ID,ROW_WORLD==MPI_COMM_NULL) ! 7-22-2017 
        
        write(node_string,'(i1)') icol
        amoeba_x_y_string = "icol"//trim(node_string)//"_"
        
        write(node_string,'(i1)') irow
        amoeba_x_y_string = trim(amoeba_x_y_string)//"irow"//trim(node_string)
        
        ! Only those non-root members are involved.
        if(contract_world/=mpi_comm_null .and. amoeba_id/=0)then
            solution_string = trim(amoeba_x_y_string)//'_ContractGroup.txt'
            open(unit=my_id+1001, file=solution_string, action='write')
            !write(my_id+1001,*) ' mamapapa ', amoeba_id, my_id, contract_id, row_id, contract_world/=mpi_comm_null
        elseif(contract_world==mpi_comm_null .and. amoeba_id/=0)then
            solution_string = trim(amoeba_x_y_string)//'_Non-contractGroup.txt'
            open(unit=my_id+1001, file=solution_string, action='write')
            !write(my_id+1001,*) ' mamapapa ', amoeba_id, my_id, contract_id, row_id, contract_world/=mpi_comm_null
        endif
        
        call read_matrix(pointlist,'_whole_list.csv')
        if(setindex*(subdim+1)>listlength) write(*,*) " setindex is too big to link to valid data in _whole_list.csv "
        
        pts_ndim = pointlist( 1+(setindex-1)*(subdim+1):setindex*(subdim+1), 15:24 ) ! 10 dimensions as usual. <--danger       
        
        allocate(selected_input(ndim), bestvertex(ndim), result(ndim), pt_input_ndim(ndim))
        
        !if(contract_world/=mpi_comm_null) 
        pt_input_ndim = pts_ndim(contract_id+1,:) ! pt_input_ndim is to be passed into perspective node so that each vertex has the full information of the starting point with complete 'ndim' dimensions.
        
        ! 10.20.2017 every member in the amoeba needs to be activated by calling the amoeba_algorithm. [proved]
        ! 10.15.2017 This is a required step for the current case, but essentially a redundant step for computing the result, because
        ! we want to fake that we only have one starting point (actually we have subdim+1 points to start with).
        trial = setindex !10.19.2017 fake input
        !call dimension_reduction(mpi_exercise_mode, pts_ndim, pts_subdim)
        selected_input = pt_input_ndim ! 10.15.2017 A fake best point to be fed into the algorithm so that the algorithm can be used in the mpi_exercise_mode==5 mode.
        call amoeba_algorithm(trial, selected_input, bestvertex, result, bestobjval ) !Note again that we do not use these outputs. All what we do here is to initiate the ammoeba algorithm.   
        !endif ! my_id==0
        
        close(my_id+1001)
        if(my_id==0)then
            close(my_id+2001)
            close(my_id+3001)
        endif
        deallocate(selected_input, result, bestvertex)
        deallocate(pointlist, pts_ndim, pts_subdim, pt_input_ndim)
        
    elseif(mpi_exercise_mode==6)then
        
        !print*, 'mode6taskid: ', mode6taskid
        !11.2.2017 printout final matrix (distribution, c, a, h, u)
        
        !call testsub('string1', i)
        !write(*,*) '1: ', i
        !call testsub('str2', i)
        !write(*,*) '2: ', i   
        
        call read_parameter_model(para,'_1parameter.txt') ! parameter for benchmark model
        print*, 'mode6taskid1: ', mode6taskid
        
        call read_parameter_model(para,'_1parameter_trial.txt','_1parameter_trial.txt') ! parameter for policy experiment
        print*, 'mode6taskid2: ', mode6taskid
        
        if(printout26)then
            
            print*, 'mode6taskid: ', mode6taskid
            write(idmode6,'(i3.3)') mode6taskid
            
            stringmode6 = 'sef_'//trim(idmode6)//'.txt'
            !print*, '1: ', stringmode6
            call ss(sef,stringmode6,20,8)
            stringmode6 = 'hom_'//trim(idmode6)//'.txt'
            !print*, '2: ', stringmode6
            call ss(sw_ini_house,stringmode6,20,8)
            stringmode6 = 'csp_'//trim(idmode6)//'.txt'
            !print*, '3: ', stringmode6
            call ss(sw_consumption,stringmode6,20,8)
            
        endif !mpi_exercise_mode     
        
    endif ! mpi_exercise_mode
    
    !call search_equilibrium(exit_log1) ! <===== replace solve_model() with this one. 3.10.2017 This is the working one. Obsolete, 7-3-2017.
    
    call mpi_finalize( MPI_ERR ) !! TO BE MOVED TO THE END OF THE MAIN PROGRAM <-------  
    
    call system_clock(tend) 
    if(my_id==0) write(*,fmt='(/,a,f12.4,a,x,i3)') 'total time: ',real(tend-tstart,wp)/real(trate,wp), ' seconds', my_id
    
    !call end_files_for_writing() ! close files ! 8-26-2017 
        
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
    ! NOTE: the system root never calls amoeba_algorithm, so don't use any MPI_BARRIER with MPI_COMM_WORLD; otherwise, the program would stuck where your MPI_BARRIER is placed (due to the absence of the system root).    
    subroutine amoeba_algorithm( trial, selected_input_dup, bestvertex, result, bestobjval ) ! 7-28-2017 Note: "selected_input" is a macro variable defined outside this subroutine.
        implicit none
        real(wp), dimension(:), intent(out) :: bestvertex, result
        real(wp), intent(out) :: bestobjval
        real(wp), dimension(:), intent(in) :: selected_input_dup
        real(wp), dimension(:), allocatable :: selected_input
        integer, intent(in) :: trial
        integer :: ndim, num_vertices, amo_msgtype, major_counter, good_msg_count
        integer :: shrink_counter, subworst_idx, i, j
        real(wp) :: plus_deviation, minus_deviation, true_deviation, sim_objval, best_val
        real(wp) :: tryfun, vr, ve, vcr, vco, temp_dist, temp_ans
        real(wp), dimension(:), allocatable :: sim_moments, ray_objval, dup_ray_objval, centroid, best_posi, temp_contract_vertex
        real(wp), dimension(:), allocatable :: worst_valvec, less_worse_valvec, try_vec, current_best
        real(wp), dimension(:), allocatable :: trymoms_vec, tryfun_vec, subworst_posi, distance_vec, old_objval_vec
        real(wp), dimension(:,:), allocatable :: vertex_list, mat_moments, worst_mat, try_mat, trymoms_mat
        real(wp), dimension(:,:), allocatable :: old_vertex_list, old_mat_moments
        integer, dimension(:), allocatable :: rankinglist, ray_modelmsg, old_ray_modelmsg, trymsg_vec ! Bug. ! 9-23-2017
        logical, dimension(:), allocatable :: shrink_signal_ray
        integer :: shrink_flag, mid_output_count, trymsg ! 9-23-2017
        integer :: pdim, trace_counter !10.17.2017
        character(len=4) :: str_amoeba
        
        ndim = size(bestvertex) !10.15.2017 Don't change it. Otherwise, the local variable ndim takes nothing and leads to error.
        
        if(mpi_exercise_mode/=5)then !10.17.2017 By default, we run the original amoeba algorithm done in the summer of 2017.
            pdim = ndim
            num_vertices = ndim+1
        else ! mpi_exercise_mode==5
            write(*,'(a,i3)') " subdim value: ", subdim
            pdim = subdim
            num_vertices = subdim+1
        endif
        
        !num_vertices = ndim + 1 ! 8-2-2017 Not necessary equal to the number of nodes involved in the amoeba operation.
            
        !allocate( vertex_list(ndim,ndim+1), sim_moments(ndim), mat_moments(ndim,ndim+1), ray_objval(ndim+1), ray_modelmsg(ndim+1), rankinglist(ndim+1) )
        !allocate( dup_ray_objval(ndim+1), centroid(ndim), best_posi(ndim), selected_input(ndim) )
        !allocate( worst_mat(ndim,noamoeba), worst_valvec(noamoeba), less_worse_valvec(noamoeba), subworst_posi(ndim) )
        !allocate( trymoms_vec(ndim), trymoms_mat(ndim,4*noamoeba), try_vec(ndim), try_mat(ndim,4*noamoeba), tryfun_vec(4*noamoeba), trymsg_vec(4*noamoeba) )
        !allocate( shrink_signal_ray(noamoeba), distance_vec(ndim+1), old_vertex_list(ndim,ndim+1), old_objval_vec(ndim+1), old_mat_moments(ndim,ndim+1) )
        !allocate( old_ray_modelmsg(ndim+1), current_best(ndim), temp_contract_vertex(ndim) )
        
        allocate( vertex_list(ndim,pdim+1), sim_moments(ndim), mat_moments(ndim,pdim+1), ray_objval(pdim+1), ray_modelmsg(pdim+1), rankinglist(pdim+1) )
        allocate( dup_ray_objval(pdim+1), centroid(ndim), best_posi(ndim), selected_input(ndim) )
        allocate( worst_mat(ndim,noamoeba), worst_valvec(noamoeba), less_worse_valvec(noamoeba), subworst_posi(ndim) )
        allocate( trymoms_vec(ndim), trymoms_mat(ndim,4*noamoeba), try_vec(ndim), try_mat(ndim,4*noamoeba), tryfun_vec(4*noamoeba), trymsg_vec(4*noamoeba) )
        allocate( shrink_signal_ray(noamoeba), distance_vec(pdim+1), old_vertex_list(ndim,pdim+1), old_objval_vec(pdim+1), old_mat_moments(ndim,pdim+1) )
        allocate( old_ray_modelmsg(pdim+1), current_best(ndim), temp_contract_vertex(ndim) )    

        if(my_id==0.and.printout23==.true.)then
            open(unit=my_id+5001,file='id001_fid5001_amoeba_trace.txt',action="write",status="replace")    
        endif
        
        ! Initialization
        bestvertex = 0._wp
        result = 0._wp
        bestobjval = 0._wp
        
        amo_msgtype = 0
        major_counter = 0 
        trace_counter = 0
        shrink_counter = 0
        shrink_flag = 0
        mid_output_count = 0
        temp_dist = 0._wp
        temp_ans = 0._wp
        
        selected_input = selected_input_dup !10.16.2017 Keep it. Important, although it is a fake input in the mpi_exercise_mode==5 case.
        
        do while( major_counter<=amoitrcrt )
            
            !if(amoeba_id==0) write(*,*) " major_counter ", major_counter, " trial ", trial
            
            major_counter = major_counter + 1
            
            if( amo_msgtype==0 )then ! Initialization of additional "ndim" vertex position. 7-20-2017
                
                !if( major_counter == 1 .and. printout13 )then
                !    
                !    if(printout13) write(*,'(a,i4,a,i4,a,i4)') "Info - My_id ", my_id, " Amoeba_id ", amoeba_id, " Contract_id", contract_id ! 7-21-2017 OK.
                !    call MPI_BARRIER( AMOEBA_WORLD, mpi_err )
                !    
                !    if(printout13) write(*,'(a,i3,a,<ndim>f12.7)') "Same ini. pt - amo_no ", amoeba_id, " pos ", selected_input ! 7-23-2017 checked.
                !    call MPI_BARRIER( AMOEBA_WORLD, mpi_err )
                !    
                !endif
                
                call MPI_BCAST( major_counter, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err ) ! Bug. 7-30-2017 AMOEBA_WORLD --> ROW_WORLD. [Trick] You should not put MPI_BCAST without appropraite if statement, if some vertices are not expected to take part of the broadcast operation.                 
                !call MPI_BCAST( major_counter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err ) ! Bug. Don't use it, because only members in the amoeba group are called by this subroutine.
                                
                if( CONTRACT_WORLD/=MPI_COMM_NULL )then
                    
                    if(mpi_exercise_mode/=5)then
                        if( CONTRACT_ID/=0 )then         
                            if(printout14)then ! extreme method using only the maximum coordinate.
                                ! http://bit.ly/2vpVmm0 p.276 See "why sets initau to 4._wp?"
                                true_deviation = initau*maxval(selected_input) !10.15.2017 a scalar.
                                
                                !if( contract_id==1 )then
                                !    print*, "selected_input: ", selected_input
                                !    print*, "max selected_input: ", maxval(selected_input)
                                !    print*, "true_deviation: ", true_deviation
                                !endif            
                                selected_input(CONTRACT_ID) = selected_input(CONTRACT_ID) + true_deviation
                            else
                                plus_deviation  = ( range_guess(CONTRACT_ID,2)-selected_input(CONTRACT_ID) ) !10.15.2017 a scalar
                                minus_deviation = ( range_guess(CONTRACT_ID,1)-selected_input(CONTRACT_ID) ) !10.15.2017 a scalar
                                
                                ! Pick the largger deviation for creating the initial potion of every vertex. 7-20-2017
                                true_deviation  = merge( plus_deviation, minus_deviation, abs(plus_deviation)>=abs(minus_deviation) )  
                                
                                ! 8-4-2017 https://stackoverflow.com/questions/17928010/choosing-the-initial-simplex-in-the-nelder-mead-optimization-algorithm
                                selected_input(CONTRACT_ID) = selected_input(CONTRACT_ID) + stepsize/100._wp*true_deviation ! 7-20-2017 checked. It is direction sensitive design. Don't worry.
                                selected_input = int(selected_input*accupara)/accupara ! 9-24-2017
                                
                                !if(printout13)then
                                !    write(*,'("M_id", i3, " C_id ", i3, " LB ", f8.5, " HB ", f8.5, " Final ",f12.7, " Minus ", f12.7, " Positive ", f12.7, " Devi ", f12.7)') &
                                !        & my_id, contract_id, range_guess(CONTRACT_ID,1), range_guess(CONTRACT_ID,2), selected_input(CONTRACT_ID), minus_deviation, &
                                !        & plus_deviation, true_deviation
                                !endif ! printout13           
                            endif ! printout14
                        endif ! CONTRACT_ID/=0

                        call MPI_ALLGATHER( selected_input, ndim, MPI_DOUBLE_PRECISION, vertex_list, & ! 7-21-2017 checked. 
                            & ndim, MPI_DOUBLE_PRECISION, CONTRACT_WORLD, mpi_err )                        
                    else ! mpi_exercise_mode==5
                        !pts_subdim = int(pts_subdim*accupara)/accupara !This step is done in amo_msgtype==1.
                        vertex_list = transpose(pts_ndim) !10.16.2017 pts_ndim's dimension is (subdim+1) X ndim.
                    endif ! mpi_exercise_mode==5
                    
                    if(my_id==0.and.printout23==.true.)then !#0
                        trace_counter = trace_counter + 1
                        do i = 1, pdim+1
                            write(my_id+5001,'((a,i5),2(a,i2),a,<ndim>(f10.7,x),(a,e15.8))') " trail#", trace_counter, "vertex#", i, " amotyp", amo_msgtype, "  input", vertex_list(:,i), " value ", ray_objval(i)  
                        enddo ! i
                    endif
                                 
                    if( CONTRACT_ID == 0 )then ! checked! 7-22-2017
                        write( my_id+1001, '("[amo0] major_counter:",x,i6, x,"mid_output_count:",i6)') major_counter, mid_output_count
                        write( my_id+1001, '(7x,a)') "List of initial vertices (contract group): "
                        !do j = 1, ndim                                                       
                        !    write( my_id+1001, '("Row", i3, " input ", <ndim+1>f12.7)')  j, vertex_list(j,:) 
                        !enddo
                        
                        ! 8-24-2017 Reformatting
                        write( my_id+1001,'(a,x)', advance='no' ) "vertex" 
                        do j = 1, ndim-1
                            write( my_id+1001, '((14x,"dim",x,i2),2x)', advance='no' ) j
                        enddo ! j
                        write( my_id+1001, '((14x,"dim",x,i2))') ndim
                        
                        !10.20.2017
                        if(printout24) call matrix_boundary_adjustment(vertex_list,range_guess) 
                        vertex_list = int(vertex_list*accupara)/accupara
                        
                        do j = 1, pdim+1 
                            write( my_id+1001, '(i6,x,<ndim>(f20.15,2x))' ) j, vertex_list(:,j)
                        enddo ! j
                    endif                   
                    
                endif ! CONTRACT_WORLD/=MPI_COMM_NULL
                
                if( CONTRACT_ID == 0 ) amo_msgtype = 1 ! move to eval & sort.
                call MPI_BCAST( trace_counter, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                call MPI_BCAST( amo_msgtype, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err ) ! Bug. 7-30-2017 AMOEBA_WORLD --> ROW_WORLD   
                
                !if( amoeba_world/=mpi_comm_null)then
                !    write(*,'(2i4)') my_id, amoeba_id    
                !endif
            endif ! amo_msgtype == 0 initial allocation.
            
            if( amo_msgtype == 1 )then ! Evaluation and sorting by contract group.
                write(my_id+1001, '(2/,"#### [amo1] major_counter:",x,i6, x,"mid_output_count:",i6)') major_counter, mid_output_count
                write(my_id+1001,'(a)') '[amo1] Purpose: evaluation and sorting by contract group'
                
                ! Evaluation with CONTRACT_GROUP
                if( CONTRACT_WORLD/=MPI_COMM_NULL )then ! 8-2-2017 Consistent with the statement made in line 735 (the MPI_ALLGATHER function askS the CONSTARCT_WORLD for the operation).
                    
                    if( major_counter/=1 )then ! 7-28-2017 checked. The shrinking outcome is idnetical across vertices in the contract world.
                        do i = 1, pdim+1
                            write(my_id+1001,'(" Input ", <ndim>f12.7)') vertex_list(:,i)
                        enddo                            
                    endif
                    
                    vertex_list(:,contract_id+1) = int(vertex_list(:,contract_id+1)*accupara)/accupara ! 9-24-2017
                    selected_input = vertex_list(:,CONTRACT_ID+1)
                    !write(my_id+1001,'(/,"[amo1] counter ", i3, " input ", <ndim>f12.7)') major_counter, selected_input
                    
                    ! Core computation
                    modelmsg = 0
                    call read_parameter_model(para, '_1parameter.txt')
                    call search_equilibrium( selected_input, sim_moments, sim_objval, my_id, trial, modelmsg )         
                    !if( contract_id == 1 ) modelmsg = 1 ! Used for debug. Comment out.
                    
                    !if(printout13 .and. contract_id==0 )then
                    if(printout13 )then ! debug 10.19.2017
                        write(my_id+1001, '(/,"#### [amo1] output of search_equilibrium. My_ID:", i4, ", major_count: ",i6, x, " mid_output_count: ",i6)') my_id, major_counter, mid_output_count
                        write(my_id+1001, '(" my_id ", i3, " contract_id ", i3, " sim_val ",f16.7," modelmsg ", i3)') my_id, contract_id, sim_objval, modelmsg 
                        write(my_id+1001, '(" input  ", <ndim>f12.7 )') selected_input ! 7-23-2017 check for input. Correct!
                        write(my_id+1001, '(" moment ", <ndim>f12.7 )') sim_moments
                    endif ! printout13
                    
                    !if(printout13 .and. contract_id==0 )then
                    !    write(my_id+1001, '(a)') ' haha '
                    !endif
                    
                    ! Communication after model computation ! vertex_list has been communicated across nodes before turning back to amo_msgtype=1 
                    call MPI_GATHER( sim_objval, 1, MPI_DOUBLE_PRECISION, & 
                        &            ray_objval, 1, MPI_DOUBLE_PRECISION, 0, CONTRACT_WORLD, mpi_err )
                    
                    if(printout13 .and. contract_id==0 )then
                        write(my_id+1001, '(a,3i4)') ' haha ', ndim, pdim, subdim
                    endif
                    
                    call MPI_GATHER( modelmsg, 1, MPI_INTEGER, &
                        &        ray_modelmsg, 1, MPI_INTEGER, 0, CONTRACT_WORLD, mpi_err )
                    
                    call MPI_GATHER( sim_moments, ndim, MPI_DOUBLE_PRECISION, &
                        &            mat_moments, ndim, MPI_DOUBLE_PRECISION, 0, CONTRACT_WORLD, mpi_err )           
                    
                    ! BEFORE SORTING ! 7-24-2017 Double checked with other output.
                    if( printout13 .and. contract_id==0 )then ! 7-23-2017 checked. 
                        write(my_id+1001,'(/,a)') "[amo1] (continued) gathered vertex_ist:"
                        do i = 1, pdim+1
                            write(my_id+1001,'("nonsorted-val ", f20.7, " msg ", i3, " input ", <ndim>f20.7, " mom ", <ndim>f20.7)')  ray_objval(i), ray_modelmsg(i), vertex_list(:,i), mat_moments(:,i)
                        enddo
                    endif
                    write(my_id+1001,'(a)') " "
                    
                    ! 8-6-2017
                    if(shrink_flag==1)then ! 9-23-2017 It means the program focus comes from amo_msgtype == 4.
                        
                        if( printout13 .and. contract_id==0 )then
                            write(my_id+1001,'("Comparision:")') 
                            do i = 1, pdim+1
                                write(my_id+1001,'("ray_objval: ", f16.7, x, "old_objval_vec: ", f16.7, x, " old_input ", <ndim>f12.7, " old_mom ", <ndim>f12.7)')  ray_objval(i), old_objval_vec(i), old_vertex_list(:,i), old_mat_moments(:,i)        
                            enddo ! i
                        endif ! printout13 .and. contract_id==0                        
                        
                        do i = pdim+1, 2, -1 ! Starting from the worst vertex, compare the newly-obatined value with its old one and update if it does improve. Once found one from the worst vertices imporved, stop this update procedure.
                            if(ray_objval(i)<old_objval_vec(i))then
                                do j = 1, i-1
                                    ray_objval(j)    = old_objval_vec(j)    ! 9-26-2017 store for next comparison.
                                    vertex_list(:,j) = old_vertex_list(:,j) ! 9-26-2017 store for next comparison.
                                    mat_moments(:,j) = old_mat_moments(:,j) ! 9-26-2017 store for next comparison.
                                    ray_modelmsg(j)  = old_ray_modelmsg(j)  ! 9-26-2017 store for next comparison.
                                enddo ! j
                                if( printout13 .and. contract_id==0 )then
                                    write(my_id+1001,'(/,"Shrinkage update happens from the last up to ", x, i4, x, " vertex.")') i
                                endif ! printout13 and contract_id==0
                                exit
                            endif ! ray_objval(i)<old_objval_vec(i)
                        enddo ! i
                        
                        if( printout13 .and. contract_id==0 )then
                            write(my_id+1001,'(/)') 
                            do i = 1, pdim+1
                                write(my_id+1001,'("Update, unsorted: ", f16.7, " msg ", i3, " input ", <ndim>f12.7, " mom ", <ndim>f12.7)')  ray_objval(i), ray_modelmsg(i), vertex_list(:,i), mat_moments(:,i)        
                            enddo ! i
                        endif ! printout13 .and. contract_id==0
                        
                        shrink_flag = 0 ! set it back to default value. ! 9-24-2017
                        
                    endif ! shrink_flag                    
                
                    if( contract_id==0 )then
                        ! SORTING
                        rankinglist = [(i,i=1,pdim+1)]
                        call sort_real( ray_objval, rankinglist ) ! sorted ray_objval (Penalty)
                        ray_modelmsg = ray_modelmsg( rankinglist ) ! sorted ray_modelmsg (Output logical array))
                        do j = 1, ndim
                            vertex_list(j,:) = vertex_list(j,rankinglist) ! sorted vertex_list (Trial vertex)
                            mat_moments(j,:) = mat_moments(j,rankinglist) ! sorted mat_moments (Simmulated moments))
                        enddo
                
                        if( printout13 .and. contract_id==0 )then
                            write(my_id+1001,'(/)') 
                            do i = 1, pdim+1
                                write(my_id+1001,'("Sorted: ", f16.7, " msg ", i3, " input ", <ndim>f12.7, " mom ", <ndim>f12.7)')  ray_objval(i), ray_modelmsg(i), vertex_list(:,i), mat_moments(:,i)        
                            enddo ! i
                        endif ! printout13 .and. contract_id==0                        
                        
                        good_msg_count = count( ray_modelmsg==0 ) ! 0 means the triral vertex leads to a valid output.
                        
                        ! EVALUATION
                        if( good_msg_count == pdim+1 )then ! every vertex returns a valid result. ! [REVISION NEEDED] 
                            write(my_id+1001,'(/,a)') '[amo1] No abnormal ray_modelmsg. Prepare for AMOEBA [-->amo3]'
                            amo_msgtype = 3 ! enter amoeba stage      
                            
                            ! new block 8-3-2017 to check whether the amoeba is small engough to exit the algorithm.
                            call maximum_distance_vertices( vertex_list, noamoeba, distance_vec)
                            call maximum_penalty_distance( ray_objval,temp_ans) ! 9-24-2017
                            !if(maximum_distance_vertices(vertex_list,noamoeba)>maxdist)then ! 7-30-2017 Has not fallen into the neighborhood.
                            if(maxval(distance_vec(2:pdim+1))>maxdist*sum(vertex_list(:,1)**2._wp)**0.5_wp)then
                                amo_msgtype = 3 ! Keep staying in amo_msgtype = 3.
                                
                                if(temp_ans<errvalfc) amo_msgtype = 10 ! 9-24-2017
                            else
                                amo_msgtype = 10 ! Normal stop
                            endif 
                            
                            write(my_id+1001,'(a," max|v-v1|: ", e15.7, " previous one: ", e15.7, " |v1|: ",f12.7,/)') 'Distance: ', maxval(distance_vec(2:pdim+1)), temp_dist, sum(vertex_list(:,1)**2._wp)**0.5_wp                            
                            temp_dist = maxval(distance_vec(2:pdim+1)) ! update 9-24-2017
                            
                            ! Bookkepping for restarts.
                            if( AMOEBA_ID==0 )then
                                mid_output_count = mid_output_count + 1
                                if(mid_output_count==1)then
                                    write(my_id+2001,'(i5,f24.14,<ndim>f24.14)') mid_output_count, ray_objval(1), vertex_list(:,1) 
                                    write(my_id+3001,'(i5,f24.14,<ndim>f24.14)') mid_output_count, ray_objval(1), vertex_list(:,1) 
                                    current_best = vertex_list(:,1)
                                else
                                    !temp = maxval(abs(current_best-vertex_list(:,1))) 
                                    !if( temp >1.e-8_wp )then ! There exist some changes in value of any digits.
                                        write(my_id+2001,'(i5,f24.14,<ndim>f24.14)') mid_output_count, ray_objval(1), vertex_list(:,1) 
                                        write(my_id+3001,'(i5,f24.14,<ndim>f24.14)') mid_output_count, ray_objval(1), vertex_list(:,1) 
                                        write(my_id+1001,'("mid_output_count: ",i5,x,"current lowest: ",f24.14,x,"best point: ",<ndim>f24.14)') mid_output_count, ray_objval(1), vertex_list(:,1)
                                        current_best = vertex_list(:,1)
                                    !endif
                                endif ! mid_output_count==1
                            endif
                            
                            !if(my_id==0.and.printout23==.true.)then !#2
                            !    trace_counter = trace_counter + 1
                            !    do i = 1, pdim+1
                            !        write(my_id+5001,'((a,i5),2(a,i2),a,<ndim>(f10.7,x))') " trail#", trace_counter, "vertex#", i, " amotyp", amo_msgtype, "  input", vertex_list(:,i)  
                            !    enddo ! i
                            !endif                            
                            
                        elseif( 1<=good_msg_count .and. good_msg_count<pdim+1 )then ! Simply for initial shrinkage.
                            
                            write(my_id+1001,'(/,a)') '[amo1] prepare for initial shrinkage toward the best points. [-->amo2]' ! Contract group.
                            
                            amo_msgtype = 2 ! shrinkage by contract group toward the best point with only those invalid vertices.
                            
                        else ! It seems unlikely that even the origianl input is not valid, so this case seems impossibly to happen.
                            write(my_id+1001,'(/,a)') '[amo1] prepare to Exit '
                            
                            amo_msgtype = 9 ! return and exit  
                            
                        endif ! good_msg_count
                    endif ! contract_id
                    
                endif ! CONTRACT_WORLD/=MPI_COMM_NULL
                
                call MPI_BARRIER( ROW_WORLD, mpi_err ) ! 8-2-2017 Trick: remeber to use the error message lest the MPI implementation kills itself abruptly.
                
                ! 7-28-2017 checked.
                ! 8-2-2017 Trick: a MPI_BCAST function must not be seen by nodes not included in the working amoeba world (fortunately, it doesn't matter anyway, due to the comprehension of the world ROW_WORLD.)
                if( ROW_WORLD/=MPI_COMM_NULL )then
                    call MPI_BCAST(    amo_msgtype, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                    call MPI_BCAST(  major_counter, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                    call MPI_BCAST( shrink_counter, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                    call MPI_BCAST(   ray_modelmsg, pdim+1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                    call MPI_BCAST(     ray_objval, pdim+1, MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err )
                    call MPI_BCAST(    mat_moments, ndim*(pdim+1), MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err )
                    call MPI_BCAST(    vertex_list, ndim*(pdim+1), MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err )
                endif
                
                ! if( printout13 .and. contract_id==0 )then ! 7-28-2017 That other vertices have the same info is checked. 7-23-2017 checked. ! 7-24-2017 Double checked with other output.
                !     write(my_id+1001,'(/,a)') "[amo1] After broadcast, every node has been updated."
                !     do i = 1, ndim+1
                !         write(my_id+1001,'("sorted-value  ", f20.7, " msg ", i3, " input ", <ndim>f20.7, " mom ", <ndim>f20.7)')  ray_objval(i), ray_modelmsg(i), vertex_list(:,i), mat_moments(:,i)
                !     enddo
                ! endif                
                
            endif ! amo_msgtype == 1 7-28-2017 fully checked.
            
            ! 8-3-2017 What differentiates this block from block 4? This shrinkage is performed by individual nodes in the contract world, but the shrinkage in block 4 is performed by All the nodes in the contract world.
            if( amo_msgtype == 2 )then ! Shrinkage toward the "best point" only on those invalid vertices of "contract group for the very first few runs..
                
                shrink_counter = shrink_counter + 1
                
                if( CONTRACT_WORLD/=MPI_COMM_NULL )then                    
                    
                    if(contract_id/=0) write(my_id+1001,'(/,a,a,i4,a,i6,a,i6)') '[amo2] contract shrinking toward the best ', ' shrink_counter: ', shrink_counter, ' major_count: ', major_counter, ' mid_output_count: ', mid_output_count
                    
                    do i = 1, ndim ! ["Discrimination strategy" compared to amo_msgtype == 4] 8-5-2017. Only shrink those bad vertices toward the best point. ! 7-24-2017 Double checked with other output.
                        
                        ! Just check to see if we have on the same page across vertices after transition from the case amo_msgtype=1
                        if( printout13 .and. contract_id == 0 .and. i == 1 )then
                            write(my_id+1001, '(a)', advance='no') "[amo2] Before shrinking. "
                            do j = 1, num_vertices
                                write(my_id+1001, '(i3)', advance='no') ray_modelmsg(j)
                            enddo
                            write(my_id+1001, '(a)') ' '
                            do j = 1, num_vertices
                                write(my_id+1001, '("Vertex ", i3, x,"MSG ", i3, " POS ", <num_vertices>f12.7)') j, ray_modelmsg(j), vertex_list(:,j)
                            enddo
                        endif ! if
                        
                        ! The shrinkage takes place in every node of the same amoeba. 7-23-2017 ! 7-24-2017 Double checked.
                        !if(printout13 .and. contract_id == 0 .and. i == 1 ) write(my_id+1001, '("[amo2] is going to move onto [amo1]")')
                        if( ray_modelmsg(i+1)==1 )then ! 8-5-2017 ONLY THOSE BAD VERTICES NEED TO BE IMPROVED BY SHRINKAGE. <------------------------==========SSSSSSSSSSSSS
                            if(printout13 .and. contract_id == 0) write(my_id+1001,'("SHRINKAGE at VERTEX (Contract_ID) no.", i3)') i+1                            
                            if(printout13 .and. contract_id == 0) write(my_id+1001,'(" BEST ONE: ", <ndim>f12.7)') vertex_list(:,1)
                            if(printout13 .and. contract_id == 0) write(my_id+1001,'(" BAD ONE : ", <ndim>f12.7)') vertex_list(:,i+1)
                            
                            ! RULE OF SHRINKAGE (toward the best point; the initial point obtained from the coarse search stage.)
                            
                            !vertex_list(:,i+1) = amotau*vertex_list(:,1) + (1._wp-amotau)*vertex_list(:,i+1) ! 7-23-2017 vertex_list(:,1) is the best vertex in the very first run of the computation.                              
                            vertex_list(:,i+1) = vertex_list(:,1) + amotau*(vertex_list(:,i+1)-vertex_list(:,1)) ! 8-9-2017                              
                            
                            if(printout13 .and. contract_id == 0) write(my_id+1001,'(" New vertex ", <ndim>f12.7)') vertex_list(:,i+1)
                        endif
                    enddo ! i
                    
                    !10.20.2017
                    if(printout24) call matrix_boundary_adjustment(vertex_list,range_guess) 
                    vertex_list = int(vertex_list*accupara)/accupara                    
                    
                    if(printout13 .and. contract_id == 0)then
                        write(my_id+1001,'(/,a)') "[amo2] Prepare to be re-evaluated --> [amo1]"
                        call maximum_distance_vertices( vertex_list, noamoeba, distance_vec)                        
                        write(my_id+1001,'(a," max|v-v1|: ", e15.7, " previous one: ", e15.7, " |v1|: ",f12.7,/)') 'Distance: ', maxval(distance_vec(2:pdim+1)), temp_dist, sum(vertex_list(:,1)**2._wp)**0.5_wp                            
                        temp_dist = maxval(distance_vec(2:pdim+1)) ! update 9-24-2017
                    endif
                    
                    if(my_id==0.and.printout23==.true.)then !#2
                        trace_counter = trace_counter + 1
                        do i = 1, pdim+1
                            write(my_id+5001,'((a,i5),2(a,i2),a,<ndim>(f10.7,x),(a,e15.8))') " trail#", trace_counter, "vertex#", i, " amotyp", amo_msgtype, "  input", vertex_list(:,i), " value ", ray_objval(i) 
                        enddo ! i
                    endif                    
                    
                    amo_msgtype = 1 ! Evaluation and Sorting by contract group.
                    
                endif ! CONTRACT_WORLD/=MPI_COMM_NULL     
                
                if( ROW_WORLD/=MPI_COMM_NULL )then
                    call MPI_BCAST( amo_msgtype, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err ) ! The head of amoeba_world is the same as that of contract world.
                    call MPI_BCAST( major_counter, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                    call MPI_BCAST( shrink_counter, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                    call MPI_BCAST( ray_modelmsg, pdim+1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                    call MPI_BCAST( ray_objval, pdim+1, MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err )
                    call MPI_BCAST( mat_moments, ndim*(pdim+1), MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err )
                    call MPI_BCAST( vertex_list, ndim*(pdim+1), MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err )                
                endif
                
            endif ! amo_msgtype == 2 ! checked! 7-24-2017
            
            if( amo_msgtype == 3 )then ! AMOEBA
                
                shrink_signal_ray = .false.
                centroid  = sum( vertex_list(:,1:(pdim+1)-noamoeba),dim=2 )/((pdim+1)-noamoeba)
                best_posi = vertex_list(:,1)
                best_val  = ray_objval(1) 
                
                do j = 1, noamoeba
                      subworst_idx    = num_vertices-noamoeba+j
                      worst_mat(:,j)  = vertex_list(:,subworst_idx)
                      worst_valvec(j) = ray_objval(subworst_idx)
                      less_worse_valvec(j) = ray_objval(subworst_idx-1)
                enddo
                
                if( printout13 .and. CONTRACT_ID==0 )then ! 7-28-2017 checked. 7-24-207 checked.
                    write( my_id+1001, '(/,a,i6,a,i6)') "[amo3] Amoeba preparation (matrix initialization), major_counter: ", major_counter, "mid_output_count: ", mid_output_count
                    write( my_id+1001, '(a, <ndim>f12.7)') "[amo3] centroid of the good vertices: ", centroid
                    write( my_id+1001, '(a, <ndim>f12.7)') "[amo3] current best vertex: ", best_posi 
                    write( my_id+1001, '(a, f16.7,/)') "[amo3] current minimum penalty: ", best_val                                
                    do j = 1, noamoeba
                        write( my_id+1001, '("[amo3] Worst no.",i2.2, " penalty: ", f20.7, " point: ", <ndim>f20.7,  " Less_worse_Val ", f20.7)') noamoeba-j+1, worst_valvec(j), worst_mat(:,j), less_worse_valvec(j)
                    enddo ! j
                    write(my_id+1001,'(/)') 
                endif ! printout13 .and. amoeba_id==0
                
                if( AMOEBA_WORLD/=MPI_COMM_NULL )then ! 7-24-2017 newly added
                    if( sicol<=noamoeba )then ! Bug for MPI_ALLGATHER. 7-25-2017
                        
                        subworst_posi = worst_mat(:,sicol) ! key.   
                        
                        if( sirow == 1)then
                            call amoeba_trial_vertex( centroid, subworst_posi,  'r', amoalp, amogam, amobet, try_vec ) ! reflection.   
                            if(printout24) call point_boundary_adjustment(try_vec,range_guess) !10.20.2017
                        elseif( sirow == 2 )then
                            call amoeba_trial_vertex( centroid, subworst_posi,  'e', amoalp, amogam, amobet, try_vec ) ! extension.    
                            if(printout24) call point_boundary_adjustment(try_vec,range_guess) !10.20.2017
                        elseif( sirow == 3 )then
                            call amoeba_trial_vertex( centroid, subworst_posi, 'cr', amoalp, amogam, amobet, try_vec ) ! contraction with the reflective point.    
                            if(printout24) call point_boundary_adjustment(try_vec,range_guess) !10.20.2017
                        elseif( sirow == 4 )then
                            call amoeba_trial_vertex( centroid, subworst_posi, 'co', amoalp, amogam, amobet, try_vec ) ! contraction with the original worst point.        
                            if(printout24) call point_boundary_adjustment(try_vec,range_guess) !10.20.2017
                        endif
                        
                        ! Core computation
                        modelmsg = 0
                        call read_parameter_model(para, '_1parameter.txt')
                        !call search_equilibrium( try_vec, trymoms_vec, tryfun, my_id, my_id, modelmsg ) ! 9-23-2017 comment out
                        try_vec = int(try_vec*accupara)/accupara
                        call search_equilibrium( try_vec, trymoms_vec, tryfun, my_id, my_id, trymsg ) ! 9-23-2017

                    endif
                    
                    !if( major_counter == 2) write(*,'(a,i3,a,f12.7,a,<ndim>f12.7)') " amoeba_id ", amoeba_id, " func_val ", tryfun, " try_vec ", try_vec ! BUG. 7-30-2017 
                    
                    call MPI_GATHER( try_vec, ndim, MPI_DOUBLE_PRECISION, &
                        &            try_mat, ndim, MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err )
                    call MPI_GATHER( trymoms_vec, ndim, MPI_DOUBLE_PRECISION, &
                        &            trymoms_mat, ndim, MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err )
                    call MPI_GATHER( tryfun, 1, MPI_DOUBLE_PRECISION, &
                        &        tryfun_vec, 1, MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err )
                    call MPI_GATHER( trymsg, 1, MPI_INTEGER, &
                        &        trymsg_vec, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )
                    
                    if( AMOEBA_ID == 0 .and. printout13 )then
                        write( my_id+1001, '(/,a)') "[amo3] Amoeba subgroup preliminary outcome: "   
                        i = 0
                        do j = 1, 4*noamoeba
                            
                            if(mod(j,4)==1)then ! r
                                str_amoeba = '  r '
                            elseif(mod(j,4)==2)then ! e
                                str_amoeba = '  e '
                            elseif(mod(j,4)==3)then ! cr
                                str_amoeba = ' cr '
                            else ! mod(j,4)==0 ! co 
                                str_amoeba = ' co '
                            endif ! if
                            
                            if(mod(j,4)==1)then
                                i = i + 1
                                write( my_id+1001, '("[amo3] original -- penalty: ", e15.7, " input: ", <ndim>e15.7)') worst_valvec(i), worst_mat(:,i)      
                            endif ! if 
                            write( my_id+1001, '("[amo3] no.",i2.2,x,a,x," penalty: ",e15.7," input: ", <ndim>e15.7, " moment: ",<ndim>e15.7)') j, str_amoeba, tryfun_vec(j), try_mat(:,j), trymoms_mat(:,j)
                        enddo
                    endif
                    
                    do j = 1, noamoeba
                        
                        vr  = tryfun_vec( 1+(j-1)*4 )
                        ve  = tryfun_vec( 2+(j-1)*4 )
                        vcr = tryfun_vec( 3+(j-1)*4 )
                        vco = tryfun_vec( 4+(j-1)*4 )
                        
                        subworst_idx = num_vertices-(noamoeba-j)
                        
                        if( vr<best_val )then ! CASE 1 Reflection is the improvement ove the initial best point.
                            if( ve<best_val )then ! Go extension.
                                if(contract_id==0) write(my_id+1001,'(a,i3,a)') " subamoeba#", j, ", (0)  e "
                                vertex_list(:,subworst_idx) = try_mat(:,2+(j-1)*4) 
                                ray_objval(subworst_idx)    = ve 
                                ray_modelmsg(subworst_idx)  = trymsg_vec(2+(j-1)*4) ! 9-23-2017
                                mat_moments(:,subworst_idx) = trymoms_mat(:,2+(j-1)*4) 
                                if(contract_id==0) write(my_id+1001,'(" val ",f16.7," input ",<ndim>f12.7," mom ",<ndim>f12.7)') ray_objval(subworst_idx), vertex_list(:,subworst_idx), mat_moments(:,subworst_idx) 
                            elseif( ve>=best_val )then ! Go reflection.
                                if(contract_id==0) write(my_id+1001,'(a,i3,a)') " subamoeba#", j, ", (1)  r "
                                vertex_list(:,subworst_idx) = try_mat(:,1+(j-1)*4)
                                ray_objval(subworst_idx)    = vr
                                ray_modelmsg(subworst_idx)  = trymsg_vec(1+(j-1)*4) ! 9-23-2017
                                mat_moments(:,subworst_idx) = trymoms_mat(:,1+(j-1)*4) 
                                if(contract_id==0) write(my_id+1001,'(" val ",f16.7," input ",<ndim>f12.7," mom ",<ndim>f12.7)') ray_objval(subworst_idx), vertex_list(:,subworst_idx), mat_moments(:,subworst_idx) 
                            endif
                        elseif( vr>=best_val )then
                            if( vr<less_worse_valvec(j) )then ! CASE 2 Go reflection. 
                                if(contract_id==0) write(my_id+1001,'(a,i3,a)') " subamoeba#", j, ", (2)  r "
                                vertex_list(:,subworst_idx) = try_mat(:,1+(j-1)*4)
                                ray_objval(subworst_idx)    = vr
                                ray_modelmsg(subworst_idx)  = trymsg_vec(1+(j-1)*4) ! 9-23-2017
                                mat_moments(:,subworst_idx) = trymoms_mat(:,1+(j-1)*4)
                                if(contract_id==0) write(my_id+1001,'(" val ",f16.7," input ",<ndim>f12.7," mom ",<ndim>f12.7)') ray_objval(subworst_idx), vertex_list(:,subworst_idx), mat_moments(:,subworst_idx) 
                            elseif( vr>=less_worse_valvec(j) )then ! Case 3 Go contraction.
                                if( vr<worst_valvec(j) .and. vcr<vr )then ! combination of A^R and M.
                                    if(contract_id==0) write(my_id+1001,'(a,i3,a)') " subamoeba#", j, ", (3) cr "
                                    vertex_list(:,subworst_idx) = try_mat(:,3+(j-1)*4)
                                    ray_objval(subworst_idx)    = vcr
                                    ray_modelmsg(subworst_idx)  = trymsg_vec(3+(j-1)*4) ! 9-23-2017
                                    mat_moments(:,subworst_idx) = trymoms_mat(:,3+(j-1)*4)
                                    if(contract_id==0) write(my_id+1001,'(" val ",f16.7," input ",<ndim>f12.7," mom ",<ndim>f12.7)') ray_objval(subworst_idx), vertex_list(:,subworst_idx), mat_moments(:,subworst_idx) 
                                elseif( worst_valvec(j)<=vr .and. vco<worst_valvec(j) )then ! combination of A_J and M.
                                    if(contract_id==0) write(my_id+1001,'(a,i3,a)') " subamoeba#", j, ", (4) co "
                                    vertex_list(:,subworst_idx) = try_mat(:,4+(j-1)*4)
                                    ray_objval(subworst_idx)    = vco
                                    ray_modelmsg(subworst_idx)  = trymsg_vec(4+(j-1)*4) ! 9-23-2017
                                    mat_moments(:,subworst_idx) = trymoms_mat(:,4+(j-1)*4)
                                    if(contract_id==0) write(my_id+1001,'(" val ",f16.7," input ",<ndim>f12.7," mom ",<ndim>f12.7)') ray_objval(subworst_idx), vertex_list(:,subworst_idx), mat_moments(:,subworst_idx) 
                                else ! Case 4 Shrink signal. 
                                    if( vr<worst_valvec(j) )then ! Bug exists in last version. 7-25-2017
                                        if(contract_id==0) write(my_id+1001,'(a,i3,a)') " subamoeba#", j, ", (5)  r "
                                        vertex_list(:,subworst_idx) = try_mat(:,1+(j-1)*4)
                                        ray_objval(subworst_idx)    = vr
                                        ray_modelmsg(subworst_idx)  = trymsg_vec(1+(j-1)*4) ! 9-23-2017
                                        mat_moments(:,subworst_idx) = trymoms_mat(:,1+(j-1)*4)    
                                        if(contract_id==0) write(my_id+1001,'(" val ",f16.7," input ",<ndim>f12.7," mom ",<ndim>f12.7)') ray_objval(subworst_idx), vertex_list(:,subworst_idx), mat_moments(:,subworst_idx) 
                                    else ! In this case, no update happens.
                                        if(contract_id==0) write(my_id+1001,'(a,i3,a)') " subamoeba#", j, ", (6)  o "
                                        !!! Keep the vertex intact 9-23-2017
                                        !! vertex_list(:,subworst_idx) = keep the same as in the last round. 
                                        !! ray_objval(subworst_idx)    = keep the same as in the last round.
                                        !! mat_moments(:,subworst_idx) = keep the same as in the last round.
                                    endif ! vr<worst_valvec
                                    
                                    shrink_signal_ray(j) = .true. ! 9-23-2017 raise the flag of unsuccessful update on a vertex.
                                    
                                endif
                            endif    
                        endif                            
                    enddo ! j: the do loop over all the worst vertices. 
                    
                endif ! AMEOBA_WORLD /= MPI_COMM_NULL
                
                if( contract_id == 0 .and. printout13 )then
                    write(my_id+1001,'(/,"[amo3] right after subgroup computation")') 
                    do i = 1, num_vertices
                        write(my_id+1001,'(" unsorted post-amoeba pentalty: ",f20.7," point: ", <ndim>f20.7, " moments: ", <ndim>f20.7)') ray_objval(i), vertex_list(:,i), mat_moments(:,i)
                    enddo ! i
                    
                    write(my_id+1001,'(a,<noamoeba>l3,/)') '[amo3] Shrink_single_ray, T, shrinkage needed; F, otherwise. -- ', shrink_signal_ray
                endif ! 
                
                if( count(shrink_signal_ray) == noamoeba )then ! The worst case contraction is needed in next step.
                    
                    amo_msgtype = 4 ! Shrink. toward the best
                    if(contract_id==0) write(my_id+1001,'(a,/)') '[amo3]-->[amo4] leave for shrinking... '
                    
                else ! At least one processor finds that Case 1 or 2 applies or A^C takes place in Case 3.
                    
                    ! ordering the result.
                    if(contract_id==0)then
                        rankinglist = [(i,i=1,pdim+1)]                  
                        call sort_real( ray_objval, rankinglist )
                        ray_modelmsg(:) = ray_modelmsg(rankinglist)
                        do j = 1, ndim
                            vertex_list(j,:) = vertex_list(j,rankinglist)
                            mat_moments(j,:) = mat_moments(j,rankinglist)
                        enddo ! j
                    endif
                    
                    ! Bookkepping for restarts. 8-7-2017
                    if( AMOEBA_ID==0 )then
                        mid_output_count = mid_output_count + 1
                        if( maxval(abs(current_best-vertex_list(:,1)))>1.e-8_wp )then ! There exist some changes in value of any digits.
                            write(my_id+2001,'(i5,f24.14,<ndim>f24.14)') mid_output_count, ray_objval(1), vertex_list(:,1) 
                            write(my_id+3001,'(i5,f24.14,<ndim>f24.14)') mid_output_count, ray_objval(1), vertex_list(:,1) 
                            current_best = vertex_list(:,1)
                        endif
                    endif                    
                                        
                    ! check if the exit critieria satisfies. 
                    call maximum_distance_vertices( vertex_list, noamoeba, distance_vec)
                    !if(maximum_distance_vertices(vertex_list,noamoeba)>maxdist)then ! 7-30-2017 Has not fallen into the neighborhood.
                    call maximum_penalty_distance( ray_objval,temp_ans) ! 9-24-2017
                    if(maxval(distance_vec(2:pdim+1))>maxdist*sum(vertex_list(:,1)**2._wp)**0.5_wp)then                    
                        amo_msgtype = 3 ! Keep staying in amo_msgtype = 3.
                        
                        if(my_id==0.and.printout23==.true.)then !#3
                            trace_counter = trace_counter + 1
                            do i = 1, pdim+1
                                write(my_id+5001,'((a,i5),2(a,i2),a,<ndim>(f10.7,x),(a,e15.8))') " trail#", trace_counter, "vertex#", i, " amotyp", amo_msgtype, "  input", vertex_list(:,i), " value ", ray_objval(i)  
                            enddo ! i
                        endif                          
                        
                        if(temp_ans<errvalfc) amo_msgtype = 10
                    else
                        amo_msgtype = 10 ! Natural stop.
                    endif
                    
                    !write(my_id+1001,'(a," max|v-v1|: ",f12.7," |v1|: ",f12.7)') 'Distance: ', maxval(distance_vec(2:ndim+1)), sum(vertex_list(:,1)**2._wp)**0.5_wp                
                    write(my_id+1001,'(a," max|v-v1|: ", e15.7, " previous one: ", e15.7, " |v1|: ",e15.7,/)') 'Distance: ', maxval(distance_vec(2:pdim+1)), temp_dist, sum(vertex_list(:,1)**2._wp)**0.5_wp                            
                    temp_dist = maxval(distance_vec(2:pdim+1)) ! update 9-24-2017                    
                    
                endif
                
                
                !! 7-30-2017 There are two possible scenarios that may make your MPI implmentation kills itself silendly during the run time:
                !! (1) You call MPI_BCAST without the discretion that the function should only be seen/encountered by the calling communicator; otherse should be avoided using, say, a if statement.
                !! (2) You call collective function such as MPI_GATHER, MPI_ALLGATHER without the discretion that the receiving buffer should have an approportionate size in accordance with the components of the calling communicator.
                call MPI_BCAST( amo_msgtype, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                call MPI_BCAST( major_counter, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                call MPI_BCAST( shrink_counter, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                call MPI_BCAST( ray_modelmsg, pdim+1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                call MPI_BCAST( ray_objval, pdim+1, MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err )          ! 9-23-2017 carry over.
                call MPI_BCAST( mat_moments, ndim*(pdim+1), MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err )  ! 9-23-2017 carry over.
                call MPI_BCAST( vertex_list, ndim*(pdim+1), MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err )  ! 9-23-2017 carry over.  
                
            endif ! amo_msgtype = 3 (Amoeba)
            
            ! 8-6-2017 Non-discrimination strategy: every vertex is required to shrink.
            ! 8-3-2017 What differentiates this block from block 2? This shrinkage is performed by All the nodes in the contract world, but the shrinkage in block 2 is performed RESPECTIVELY by nodes in the contract world.
            if( amo_msgtype == 4 )then ! 8-6-2017. Shrinkage strategy is NOT for the vertex initialization, compared to amo_msgtype == 2. 
                
                shrink_counter = shrink_counter + 1
                
                if( CONTRACT_WORLD/=MPI_COMM_NULL )then
                        
                    if(contract_id/=0) write(my_id+1001,'(/,a,a,i6)') '[amo4] contract shrinking toward the best', ' major_counter ', major_counter
                    
                    if( printout13 .and. contract_id == 0 )then
                        write(my_id+1001, '(a)') "[amo4] Before shrinking. "
                        do j = 1, num_vertices
                            write(my_id + 1001, '(" unsorted post-ameoba penalty: ", f16.7, " msg ", i3, " input ", <ndim>f12.7, " mom ", <ndim>f12.7)')  ray_objval(j), ray_modelmsg(j), vertex_list(:,j), mat_moments(:,j)                            
                        enddo ! j
                    endif ! if
                    
                    if(printout13 .and. contract_id == 0) write(my_id+1001,'(" best ", <ndim>f12.7, " to-be ", <ndim>f12.7)') vertex_list(:,1), vertex_list(:,i+1)
                    
                    ! Prepare for comparison in amo_msgtype = 1
                    old_objval_vec   = ray_objval   ! everyone should already have the same information due to the broadcast in amo_msgtype =3.
                    old_vertex_list  = vertex_list  ! everyone should already have the same information due to the broadcast in amo_msgtype =3.
                    old_mat_moments  = mat_moments  ! everyone should already have the same information due to the broadcast in amo_msgtype =3.
                    old_ray_modelmsg = ray_modelmsg ! everyone should already have the same information due to the broadcast in amo_msgtype =3.
                    shrink_flag      = 1
                    
                    ! overall contraction
                    !vertex_list(:,contract_id+1) = amotau*vertex_list(:,1) + (1._wp-amotau)*vertex_list(:,contract_id+1) ! 7-26-2017 
                    !temp_contract_vertex = amotau*vertex_list(:,1) + (1._wp-amotau)*vertex_list(:,contract_id+1) ! 8-9-2017
                    temp_contract_vertex = vertex_list(:,1) + amotau*(vertex_list(:,contract_id+1)-vertex_list(:,1)) ! 8-9-2017
                    
                    temp_contract_vertex = int(temp_contract_vertex*accupara)/accupara ! 9-24-2017
                    call MPI_GATHER( temp_contract_vertex, ndim, MPI_DOUBLE_PRECISION, &
                        &                     vertex_list, ndim, MPI_DOUBLE_PRECISION, 0, CONTRACT_WORLD, mpi_err )
                    
                endif
                
                if(contract_id/=0) write(my_id+1001,'(a,/)') '[amo4]-->[amo1] for post-shrinkage evaluation '
                
                amo_msgtype = 1 ! Bug. Missing. 7-26-2017 EVALUATION.
                
                !call MPI_BCAST( amo_msgtype, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                !call MPI_BCAST( major_counter, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                !call MPI_BCAST( shrink_counter, 1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                !call MPI_BCAST( ray_modelmsg, ndim+1, MPI_INTEGER, 0, ROW_WORLD, mpi_err )
                !call MPI_BCAST( ray_objval, ndim+1, MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err )
                !call MPI_BCAST( mat_moments, ndim*(ndim+1), MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err )
                call MPI_BCAST( vertex_list, ndim*(pdim+1), MPI_DOUBLE_PRECISION, 0, ROW_WORLD, mpi_err )   
                
                if( printout13 .and. contract_id == 0 )then
                    write(my_id+1001, '(/,a)') "[amo4] Trial shrinking. "
                    do j = 1, num_vertices
                        write(my_id + 1001, '("Vertex ", i3, " POS", <num_vertices>f12.7)') j, vertex_list(:,j)
                    enddo ! j
                endif ! if  
                
                if(my_id==0.and.printout23==.true.)then !#4
                    trace_counter = trace_counter + 1
                    do i = 1, pdim+1
                        write(my_id+5001,'((a,i5),2(a,i2),a,<ndim>(f10.7,x),(a,e15.8))') " trail#", trace_counter, "vertex#", i, " amotyp", amo_msgtype, "  input", vertex_list(:,i), " value ", ray_objval(i)  
                    enddo ! i
                endif                    
                
            endif ! amo_msgtype == 4       
            
            if( amo_msgtype == 9 )then ! Abnormal stop
                exit
            endif
            
            if( amo_msgtype == 10 )then ! Normal stop: successful amoeba converge.
                exit
            endif ! amo_msgtype == 10
            
            !write(my_id+1001,'(/,a," max|v-v1|: ",f12.7," |v1|: ",f12.7)') 'Distance: ', maxval(distance_vec(2:ndim+1)), sum(vertex_list(:,1)**2._wp)**0.5_wp
            
        enddo ! while ( major_counter<=amoitrcrt )
        
        if( amo_msgtype==9 )then 
            bestobjval = inf ! 0._wp
            bestvertex = 0._wp
            result     = 0._wp
        else ! including the cases such as amo_msgtype== 1 (evaluation), 3 (amoeba) or 10 (naturally exit).
            bestobjval = ray_objval(1)
            bestvertex = vertex_list(:,1)
            result     = mat_moments(:,1)
        endif
        
        deallocate( vertex_list, sim_moments, mat_moments, ray_objval, ray_modelmsg, rankinglist, dup_ray_objval )
        deallocate( centroid, best_posi, worst_mat, worst_valvec, less_worse_valvec, try_vec, trymoms_vec, subworst_posi )
        deallocate( try_mat, tryfun_vec, trymoms_mat, shrink_signal_ray, selected_input, distance_vec, old_vertex_list )
        deallocate( old_objval_vec, old_mat_moments, old_ray_modelmsg, current_best, temp_contract_vertex, trymsg_vec )
        
        close(my_id+5001)
        
    end subroutine amoeba_algorithm
    
    subroutine amoeba_trial_vertex( centroid, worst, opt, alpha, gamma, beta, amotry_vec )
        implicit none
        real(wp), dimension(:), intent(in) :: centroid, worst ! "worst" is the current worst vertex to be imporved. 7-24-2017.
        real(wp), intent(in) :: alpha, gamma, beta
        real(wp), dimension(:), intent(out) :: amotry_vec
        real(wp), dimension(:), allocatable :: reflect
        integer :: n
        character(len=*), intent(in) :: opt
        n = size( centroid )
        allocate( reflect(n) ) 
        reflect = centroid + alpha*( centroid - worst ) ! A^R_J
        select case(opt)
            case('r') ! A^R_J
                amotry_vec = reflect 
            case('e') ! A^E_J
                amotry_vec = reflect + gamma*( reflect - centroid ) ! 7-24-2017 Lee and Wiswall p.176
            case('cr') ! CONTRACTION AGAINST THE REFLECTION POINT
                amotry_vec = centroid + beta*( reflect - centroid ) ! beta*( centroid + reflect ) ! 7-24-2017 Lee and Wiswall p.176
            case('co') ! CONTRACTION AGAINST THE ORIGINAL WORST POINT
                amotry_vec = centroid + beta*( worst - centroid )   ! beta*( centroid + worst )
        endselect 
        deallocate( reflect )
    end subroutine amoeba_trial_vertex    
    
    ! The root sends trial parameter combination to slaves
    subroutine sendjob(trial,slave,parcel,opt,srnumber)
        implicit none
        integer, intent(in) :: trial, slave
        real(wp), dimension(:), intent(in) :: parcel
        integer :: msgtype, n
        character(len=*), optional :: opt
        integer, optional :: srnumber
        n = size(parcel)
        msgtype = 1 ! tag 1 is used for the communication of model parameter passing. tag 2 is used for result passing from slave.
        if(.not.present(opt))then ! original case: to members in MPI_COMM_WORLD. Used for mpi_exercise_mode==1 case.
            call mpi_send( trial, 1, mpi_integer, slave, &
                & msgtype, MPI_COMM_WORLD, ierr) ! trial variable to slave.
            call mpi_send( parcel, ndim, mpi_double_precision, slave, &
                & msgtype, MPI_COMM_WORLD, ierr) ! parcel variable(s) to slave.
        else ! if the flag string opt is shown, whatever value it is. We know it is in the mpi_exercise_mode==2 case.
            call mpi_send( trial, 1, mpi_integer, slave, &
                & msgtype, HEAD_WORLD, ierr)
            call mpi_send( parcel, ndim, mpi_double_precision, slave, &
                & msgtype, HEAD_WORLD, ierr)
        endif ! if condition of "present"
        
        if(present(srnumber))then !10.17.2017
            if(.not.present(opt))then
                call mpi_send( srnumber, 1, mpi_integer, slave, &
                    & msgtype, MPI_COMM_WORLD, ierr)            
            else
                call mpi_send( srnumber, 1, mpi_integer, slave, &
                    & msgtype, HEAD_WORLD, ierr)                            
            endif
        endif
        
    end subroutine sendjob
    
    !function ftest(x) 
    !    implicit none
    !    real(wp) :: ftest
    !    real(wp), intent(in) :: x
    !    !ftest = x**4._wp + x**3._wp - 6*x**2._wp + 4*x + 12
    !    ftest = x*sin(1._wp/x)
    !end function ftest
    
    subroutine read_one_point_from_short_list(trial,startpoint)
        implicit none
        real(wp), dimension(:), intent(out) :: startpoint
        integer, intent(in) :: trial
        integer :: n, iostat
        n = 0
        open(unit=10,file='_full_list_startingpoint.csv',status='old',action='read',iostat=iostat)
        if(iostat==0)then
            do
                n = n + 1
                read(10,*) startpoint
                if(n==trial) exit
            enddo
        else ! iostat
            write(*,*) "something wrong when calling read_one_point_from_short_list"
        endif ! iostat
        close(10)
    end subroutine read_one_point_from_short_list
    
    subroutine read_best_point( startpoint, input_file ) ! 10.15.2017 (ndim+2) columns per row.
        implicit none
        real(wp), dimension(:), intent(out) :: startpoint
        real(wp), dimension(:), allocatable :: temp_ray
        integer :: iostat
        character(len=800) :: value_string
        character(len=*), intent(in) :: input_file
        integer :: n, m
        !n = size(startpoint)
        m = 0
        allocate( temp_ray(ndim+2) )
        open( unit=10, file=input_file, status='old', action='read', iostat=iostat )
        if( iostat==0 )then
            do
                read( unit=10, fmt='(a)', iostat=iostat ) value_string
                if( iostat/=0 ) exit                
                if( scan( value_string,'!' )>0 ) cycle
                read( value_string, fmt=* ) temp_ray
                m = m + 1
            enddo ! do
        endif
        if( m==0 ) write(*,fmt='(/,a)') "something wrong with the subroutine read_best_point "
        startpoint = temp_ray(3:ndim+2) 
        deallocate( temp_ray )            
    end subroutine read_best_point
    
    subroutine read_short_list(filename, startpoint) !10.15.2017 compatible for reading file in the format of 25 columns: serial#, error, node 1, node 2, moments 1-10, inputs 1-10, group#.
        implicit none
        real(wp), dimension(:), intent(out) :: startpoint ! length of 25
        character(len=*), intent(in) :: filename
        integer :: ans1, ans2, sn, n, iostat
        
        ! catch the serial number from given filename
        ans1 = index(filename,'_',back=.true.)
        ans2 = index(filename,'.')
        read(filename(ans1+1:ans2-1),'(i)') sn
        
        n = 0
        open(unit=10,file='_full_list_startingpoint.csv',status='old',action='read',iostat=iostat) ! read line by line until reach the startpoint specifid in "filename" 
        if(iostat==0)then
            do
                n = n + 1
                read(10,*) startpoint
                if(n==sn) exit
            enddo
        else
            write(*,*) 'something wrong when opening read_short_list'
        endif
        close(10)
    end subroutine read_short_list
    
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
        open(unit=102,file="102_ConvergenceInfo.txt",status="replace",action="write") 
        open(unit=103,file="output_103_macro.txt",status="replace",action="write") 
        open(unit=104,file="output_104_debug.txt",status="replace",action="write") 
        open(unit=105,file="output_105_debug.txt",status="replace",action="write")
        !open(unit=106,file="output_106_debug.txt",status="replace",action="write")
        !open(unit=107,file="output_107_profitcompare.txt",status="replace",action="write")
        !open(unit=108,file="output_108_labor.txt",status="replace",action="write")
        !open(unit=109,file="output_109_profit_mat.txt",status="replace",action="write")
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
        open(unit=121,file="121_AssetsAgeProfile.txt",status="replace",action="write")
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
        open(unit=134,file="output_134_debug_cef.txt",status="replace",action="write")
        open(unit=135,file="output_135_debug_cef.txt",status="replace",action="write")
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
        !close(109)
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
        close(134)
        close(135)
    end subroutine end_files_for_writing
    
    subroutine dimension_reduction(mpi_exercise_mode, pts_ndim, pts_subdim)
        implicit none
        real(wp), dimension(:,:), intent(in) :: pts_ndim
        real(wp), dimension(:,:), intent(out) :: pts_subdim
        integer, dimension(:), allocatable :: index_vec
        integer, intent(in) :: mpi_exercise_mode
        integer :: i
        real(wp) :: kv1, prk0, prk1, prk2, zbar, beta, iota, phi1, phi2, phi3
        if(mpi_exercise_mode==5)then
            allocate(index_vec(subdim))
            index_vec = [1,2,3,5,6]
            do i = 1, subdim
                pts_subdim(:,i) = pts_ndim(:,index_vec(i))   
            enddo ! i
            deallocate(index_vec)
        endif ! mpi_exercise_mode
    end subroutine dimension_reduction
    
    subroutine point_boundary_adjustment(pt,ranges)
        implicit none
        real(wp), dimension(:), intent(inout) :: pt
        real(wp), dimension(:,:), intent(in) :: ranges
        integer :: m, i
        m = size(pt)
        do i = 1, m
            if(ranges(i,1)>pt(i))then
                pt(i) = ranges(i,1)            
            elseif(ranges(i,2)<pt(i))then
                pt(i) = ranges(i,2)
            endif
        enddo ! i
        return
    end subroutine point_boundary_adjustment
    
    subroutine matrix_boundary_adjustment(pts,ranges)
        implicit none
        real(wp), dimension(:,:), intent(inout) :: pts
        real(wp), dimension(:,:), intent(in) :: ranges
        integer :: m, n, i, j
        m = size(pts,dim=1) 
        n = size(pts,dim=2)
        if(m/=ndim) write(*,*) ' dimensional inconsistency in matrix_boundary_adjustment'
        do i = 1, m ! index of target
            do j = 1, n ! index of vertex
                if(ranges(i,1)>pts(i,j))then
                    pts(i,j) = ranges(i,1)
                elseif(ranges(i,2)<pts(i,j))then
                    pts(i,j) = ranges(i,2) 
                endif 
            enddo !j
        enddo !i
        return
    end subroutine matrix_boundary_adjustment
    
    subroutine testsub(i,j)
        implicit none
        character(len=*), intent(in) :: i
        integer, intent(out) :: j
        select case(i)
        case('string1')
            j = 135
            write(*,*) 'answer1=', j
        case('str2')
            j = 246
            write(*,*) 'answer2=', j
        end select
    end subroutine testsub
    
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
    
    !! Experiment for inf.
    !real(wp) :: temp1, temp2
    !temp1 = inf
    !temp2 = inf
    !print*, "answer: ", merge(" yes "," no  ",temp1==temp2) ! answer: yes   

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
    
    !! Test of maximum_distance_vertices function
    !real(wp), dimension(4,5) :: vermat
    !real(wp) :: vervec(20), temp
    !vervec = [(real(i,wp),i=1,20)]
    !vermat = reshape(vervec,(/4,5/))
    !temp = maximum_distance_vertices(vermat,2)
    !print*, "distance ", temp
    
    !! Test of maximum_distance_vertices function
    !real(wp), dimension(4,5) :: vermat
    !real(wp) :: vervec(20), temp
    !vervec = [(real(i,wp),i=1,20)]
    !vermat = reshape(vervec,(/4,5/))
    !temp = maximum_distance_vertices(vermat,2)
    !print*, merge(.true.,.false.,maximum_distance_vertices(vermat,2)>25), "distance ", temp        