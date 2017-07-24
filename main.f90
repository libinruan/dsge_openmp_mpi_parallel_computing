! (1) main.f90 focuses on the implementation of the Nelder-Mead simplex method on 10 dimensional parameter space.
! (2) equilibrium.f90 solves for the steady state of the economy.
! (3) model.f90 solves the Bellman equations for saving and consumption decisions across generations.
! (4) variable.f90 initializes most variables to be used in the program.
! (5) toolbox.f90 contains all the subroutines that executes evoked numerical methods.
    
program MPI_sandbox
    !use universe ! test on reversion. Should shown.
    use equilibrium
    implicit none
    integer :: tstart, tend, trate, tmax, i, j, index_pt, index_rcv
    character(len=200) :: msg
    
    ! The mpi_exercise_mode == 1 case: Coarse search
    logical :: exit_log1
    
    ! The mpi_exercise_mode == -1 case: MKL experiment of random generators
    integer :: generator,erridx, approach, tmiddle
    type(vsl_stream_state) :: river
    real(wp) :: stoch(1)
    
    call system_clock(tstart,trate,tmax)
    call system_clock(tstart)    
    
    !call start_files_for_writing() ! open files
    
    call fmpi_init() ! USER-DEFINED SUBROUTINE IN TOOLBOX.F90 <------------------------
    call infinity_setting(inf)
    
    if( MY_ID == 0 .and. MPI_PROVIDED<MPI_THREAD_FUNNELED ) write(*,'(a,/)') '! [ WARNING ] The Only-Master-makes-MPI-calls setup fails.' 
    
    allocate(range_guess(ndim,2))
    call read_parameter_model(para,'_1parameter.txt')
    
    if(my_id==0)then ! General Operation Messages
        
        if(printout12)then
            write(*,'(a,f20.8)') (labstr(i),para(i),i=1,131) ! works. 
            write(*,*) ' '
        endif
        
        if(mpi_exercise_mode==1)then
            msg = ' starting point search'
        elseif(mpi_exercise_mode==2)then
            msg = ' amoeba algorithm'
        else
            msg = ' -- not yet --'    
        endif
        
        if(mpi_exercise_mode/=0) write(*,'(a)') 'MPI_EXERCISE_MODE: ', msg
        
        if(printout11)then
            msg = ' quantitative model'
        else
            msg = ' test model'
        endif
        write(*,'(a)') 'PRINTOUT11: ', msg
    endif

    if(mpi_exercise_mode == -2)then ! experiment zone
        
        write(*,*) locate(real(breaks_list,wp),1.5_wp,1)
        write(*,*) locate(real(breaks_list,wp),50._wp,1)
        
    elseif(mpi_exercise_mode==-1)then ! MKL experiment (random number generator)
        
        generator = VSL_BRNG_MCG31
        approach  = VSL_RNG_METHOD_UNIFORM_STD
        call system_clock(tmiddle)
        erridx    = vslnewstream(river,generator,tmiddle)
        do i = 1, 20         
            erridx = vdrnguniform(approach,river,1,stoch,0._wp,1000._wp)
            write(*,*) 'my_id = ', my_id, ' random number = ', stoch
        enddo
        
    elseif(mpi_exercise_mode==0)then ! Stage 0. Building Stage with only a single node.
        
        solution_string = '_SingleNodeInputs.txt'   
        concisesolution_string = '_SingleNodeDetails.txt'
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
        origin_input(8) = phi1  
        origin_input(9) = phi2  
        origin_input(10)= phi3           
        
        ! Quasi-random Sobol sequence block # 2 move inside the MPI_exercise_mode == 1!?
        call get_sobol_sequence( sobolm, 0.0_wp, 1.0_wp, nsbq, ndim ) ! Generate ndim dimensional sobol sequence of length nsbq.
        call scale_sobol_original( transpose(sobolm), range_guess, sobolm_scaled ) 
        mpi_sobol_scaled = sobolm_scaled(:,sblno1:sblno1+trylen-1)    
        
        if(my_id==0) call sm(sobolm,'sobolm'//trim(trylen_string)) ! checked 2017-Jul-1
        if(my_id==0) call sm(transpose(mpi_sobol_scaled),'mpi_sobol_scaled'//trim(trylen_string)) ! checked 2017-Jul-1 . Used for convex combination.        
        
        allocate( parcel(ndim), result(ndim), outputinput1(trylen,2*ndim+2), obj_val_vec(trylen) )        
        nslaves = num_procs - 1 ! The root processor (my_id==0) is in charge of sending out new trial and receiving the corresponding result.

        if( my_id == 0)then
            write(trylen_string,'(i5.5,"_",i5.5)') sblno1, sblno1+trylen-1
            io_string = 'IOMat_'//trim(trylen_string) 
            
            ! [The Root, case 1] Send Initial Messages to Slave Nodes (indices ranges from 1 to nslaves)
            do i = 1, nslaves
                trial = i ! the index of the basic loop: 1,...,trylen; Not the indices of the adjusted Sobol sequence.
                slave = i
                call linear_combination_sobal_sequence(parcel,trial,mpi_sobol_scaled(:,trial),origin_input,weight_list,breaks_list)   
                mpi_sobol_mixed(trial,:) = parcel
                call sendjob(trial,slave,parcel) ! send from root to slave the trial parameter combination                
            enddo
            
            ! [The Root, case 1] Hear Responses from Slave Nodes
            do i = 1, trylen ! The range is correct!! 7-3-2017
                do
                    msgtype = 2 ! tag 1 for parameter passing; tag 2 for communicating the result between slave nodes and root. 
                    ! Non-blocking test for a message
                    call mpi_iprobe( mpi_any_source, msgtype, mpi_comm_world, &
                        & receiving, status, ierr)
                    if(receiving)then
                        ! which member node is sending the result?
                        call mpi_recv( slave, 1, mpi_integer, mpi_any_source, &
                            & msgtype, mpi_comm_world, status, ierr)
                        ! what's the returned trial index in the basic assignment loop?
                        call mpi_recv( trial, 1, mpi_integer, slave, & ! Note: trial falls in [1,trylen], rather than the shifted interval.
                            & msgtype, mpi_comm_world, status, ierr)
                        ! the feedback of simulated moments 
                        call mpi_recv( result, ndim, mpi_double_precision, slave, &
                            & msgtype, mpi_comm_world, status, ierr)
                        ! the value of penalty corresponds to the given trial
                        call mpi_recv( obj_val_1st, 1, mpi_double_precision, slave, &
                            & msgtype, mpi_comm_world, status, ierr)
                        
                        mpi_simmom_matrix(:,trial) = result
                        obj_val_vec(trial) = obj_val_1st
                        
                        ! Results That Are Collected by Individual Nodes (we are now in the my_id == 0 zone)
                        if(i==1) write(my_id+1001,'(a,(12x,a),(x,a),(2x,a),(10x,"moment1"),(10x,"moment2"),(10x,"moment3"),(10x,"moment4"),(10x,"moment5"), &
                            & (10x,"moment6"),(10x,"moment7"),(10x,"moment8"),(10x,"moment9"),(9x,"moment10"),(11x,"input1"),(11x,"input2"),(11x,"input3"), &
                            & (11x,"input4"),(11x,"input5"),(11x,"input6"),(11x,"input7"),(11x,"input8"),(11x,"input9"),(10x,"input10"))') &
                            & "MyID","error", "#trial", "#list"
                        
                        write(my_id+1001,'(i4,(x,f16.7),(2x,i5),(2x,i5),<ndim>(x,f16.7),<ndim>(x,f16.7))') &
                            & slave, obj_val_1st, trial, indexseries(trial), mpi_simmom_matrix(:,trial), mpi_sobol_mixed(trial,:)
                        
                        ! Check to see if one more new trial is available to be assigned to the responding slave node.
                        if(i<=trylen-nslaves)then
                            !parcel = mpi_sobol_scaled(:,i+nslaves) ! Correct. 7-3-2017 
                            trial = i + nslaves
                            call linear_combination_sobal_sequence(parcel,trial,mpi_sobol_scaled(:,trial),origin_input,weight_list,breaks_list)   
                            mpi_sobol_mixed(trial,:) = parcel                            
                            call sendjob( trial, slave, parcel )
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
        
        allocate( refpts(elist-slist+1), selected_input(ndim), result(ndim), bestvertex(ndim) )
        
        do i = 1, elist-slist+1 ! 7-17-2017 Generate the moving list that has the sequence of starting points.
            refpts(i) = i + (slist - 1) ! 7-22-2017 Kobol sequence index shifts ahead by slist items.
        enddo
        
        !write(*,*) ' startpoint ', startpoint ! verified.
        !if(my_id==0) write(*,*) ' refpts ', refpts ! verified.            
        ncol = (num_procs-1)/nrow ! 7-13-2017 The number of amoeba heads (equivalently, the number of amoeba).
        if( abs( (num_procs-1)/nrow-floor(real((num_procs-1)/nrow),wp) )>5.e-10_wp .and. my_id==0 ) write(*,*) 'The total nodes should be a multiple of four plus one'
        
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
        allocate( id_list(1:nrow) ) 
        id_list = 1
        do i = 2, nrow
            id_list(i) = id_list(i-1) + 1
        enddo
        
        AMOEBA_WORLD = MPI_COMM_NULL ! Bug. 7-22-2017 This is needed for solving the bug below. Without it, the system root is unexpectedly included in the AMOEBA_WORLD.
        ! 7-13-2017 Trick to create a communicator by rolling over the complete list of ranks with a fixed length of rolling window.
        do i = 1, ncol ! 7-13-2017 we have "ncol" amoeba groups to create.
            ! Now, we want to let every member in the group GROUP_ALL to have the same plan about the subdivision of group GROUP_ALL. 7-10-2017
            ! link https://www.rc.usf.edu/tutorials/classes/tutorial/mpi/chapter9.html#9.2.5 [step 2]. Note: Step 1 was to assigne the same member list of COMM to a new group for use. It took place in toolbox.f90.            
            call MPI_GROUP_INCL( GROUP_ALL, nrow, id_list, GROUP_AMOEBA, mpi_err )
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
        
        if(my_id/=0) call MPI_COMM_RANK( AMOEBA_WORLD, AMOEBA_ID, mpi_err ) ! 7-13-2017 the rank ID of nodes in the same amoeba group.
        
        write(node_string,'(i1)') icol ! The system root returns 0 for variable "icol."
        
        if( AMOEBA_ID == 0 )then ! AMOEBA_ID=0 is also the head for contraction.
                solution_string = trim(node_string)//'_amaoeba_head.txt'
                !print*, solution_string ! checked
                open(unit=my_id+1001, file=solution_string, action='write') ! AMOEBA HEAD.
        endif

        !write(*,'(5i6)') my_id, irow, icol, amoeba_id, merge(-1,contract_id,CONTRACT_WORLD==MPI_COMM_NULL) ! 7-13-201 verified.
                
        !! ---- 1. Initial valuation group: contract group (for every amoeba); Product: communicator CONTRACT_WORLD
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
            
        !write(*,'(5i6)') my_id, irow, icol, amoeba_id, merge(-1,contract_id,CONTRACT_WORLD==MPI_COMM_NULL) ! 7-13-201 verified.
                
        !! ---- COMMUNICATION GROUP (among the roots of amoeba), Product: communicator HEADWORLD
        !! The system root is in the head group, whose rank is still 0 (zero).
        HEAD_ID = -10
        allocate( id_list(0:ncol) )
        id_list = 0
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

        !! 7-22-2017 **Negative number indicates the node is not in a specific group.**
        if(printout13) write(*,'("my_id ",i3," irow ",i3," icol ",i3," amoeba ",i3," contract ",i3," head ",i3)') my_id, irow, icol, merge(-9,amoeba_id,AMOEBA_WORLD==MPI_COMM_NULL), merge(-9,contract_id,CONTRACT_WORLD==MPI_COMM_NULL), merge(-9,head_id,HEAD_WORLD==MPI_COMM_NULL) ! 7-22-2017 
        
        sirow = mod( AMOEBA_ID, 4 ) + 1 ! 7-17-2017 The order is [2,3,4,1]. The ID "within" a subamoeba group.
        sicol = AMOEBA_ID/4 + 1 ! 7-17-2017 The ID of a "subamoeba group" within an amoeba.
        nslaves = (num_procs-1)/nrow ! 7-17-2017 Number of "head nodes" except for the system root; equivalently, the number of amoebas.
        
        write(node_string,'(i1)') icol
        amoeba_x_y_string = trim(node_string)//"_"
        write(node_string,'(i3.3)') AMOEBA_ID
        amoeba_x_y_string = trim(amoeba_x_y_string)//node_string 
        
        !write(*,'("---my_id ",i3," irow ",i3," icol ",i3," amoeba ",i3," contract ",i3," head ",i3)') my_id, irow, icol, merge(-9,amoeba_id,AMOEBA_WORLD==MPI_COMM_NULL), merge(-9,contract_id,CONTRACT_WORLD==MPI_COMM_NULL), merge(-9,head_id,HEAD_WORLD==MPI_COMM_NULL)
        
        if( AMOEBA_WORLD/=MPI_COMM_NULL .and. AMOEBA_ID /= 0 )then ! 7-22-201 This excludes the system root and the amoeba heads.
                solution_string = trim(amoeba_x_y_string)//'_amaoeba_member.txt' ! For example, 1_i_amoeba_member.txt, where i = 1,...,10 (except 0 the amoeba head).
                open(unit=my_id+1001, file=solution_string, action='write') ! AMOEBA MEMBERS.
        endif        
        
        if( my_id == 0 )then
            write(trylen_string,'(i5.5,"_",i5.5)') slist, elist
            solution_string = 'RootIntermCollect_'//trim(trylen_string)//'.txt' ! THE SYSTEM ROOT.
            open(unit=my_id+1001, file=solution_string, action='write')
            
            ! Initialization. ! Tonight. Debug starts here.
            do i = 1, nslaves ! Across all the amoeba heads.
                trial = refpts(i) ! 7-17-2017 Assign a trial ID number from the reference list, refpts.
                !print*, 'My_Id=', my_id, ' slave=', i, ' trial=', trial, ' nslaves: ', nslaves ! 7-17-2017 checked. 
                slave = i ! 7-17-2017 Correct. We use the rank ID in the Head_World group.

                call read_one_point_from_short_list(trial,startpoint) ! 7-17-2017 Get a new starting point based on the input, trial.
                selected_input = startpoint(15:24)
                call sendjob(trial, slave, selected_input, 'head_group') ! 7-17-2017 trial is a new starting point's index no. from the integrated list of favorable starting points.
                if(printout13) write(my_id+1001,'("Step 1 assign trial ",i3," vertex_0 ", <ndim>f12.7, " to amo_group ", i3, " (total slaves ", i3,")")') trial, selected_input, slave, nslaves ! 
            enddo ! i
            
            ! Listen for responses
            !do i = 1, elist-slist+1
            index_pt = ndim ! 7-23-2017 We have use the first ndim starting points, so we have to now start with index ndim.
            index_rcv = 0 ! 7-23-2017 The counter for book keeping the results already received.
            do 
                msgtype = 2
                call MPI_IPROBE( MPI_ANY_SOURCE, msgtype, HEAD_WORLD, & ! 7-17-2017 Trick for the system root to scan incoming feedback from any ameoba heads.
                    & receiving, status, ierr )
                
                if( receiving )then
                    
                    call MPI_RECV( slave, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    call MPI_RECV( trial, 1, MPI_INTEGER, slave, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    call MPI_RECV( result, ndim, MPI_DOUBLE_PRECISION, slave, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    call MPI_RECV( bestvertex, ndim, MPI_DOUBLE_PRECISION, slave, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    call MPI_RECV( bestobjval, 1, MPI_DOUBLE_PRECISION, slave, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    ! bookkeeping output from each amoeba group.
                    write(my_id+1001,'(a,(12x,a),(x,a),(2x,a),(10x,"moment1"),(10x,"moment2"),(10x,"moment3"),(10x,"moment4"),(10x,"moment5"), &
                        & (10x,"moment6"),(10x,"moment7"),(10x,"moment8"),(10x,"moment9"),(9x,"moment10"),(11x,"input1"),(11x,"input2"),(11x,"input3"), &
                        & (11x,"input4"),(11x,"input5"),(11x,"input6"),(11x,"input7"),(11x,"input8"),(11x,"input9"),(10x,"input10"))') &
                        & "MyID","error", "#trial", "#trial"                        
                    write(my_id+1001,'(i4,(x,f16.7),(2x,i5),(2x,i5),<ndim>(x,f16.7),<ndim>(x,f16.7))') &
                        & slave, bestobjval, trial, trial, result, bestvertex
                    
                    index_pt = index_pt + 1
                    ! Send out new starting points
                    if( index_pt<=elist-slist+1 )then
                        trial = refpts(index_pt)
                        call read_one_point_from_short_list(trial,startpoint)
                        selected_input = startpoint(15:24)
                        call sendjob(trial, slave, selected_input, 'head_group')
                    endif
                    index_rcv = index_rcv + 1
                    
                endif ! receiving, if statement
                
                if( index_rcv>elist-slist+1 ) exit ! 7-18-2017 we have collect the outcome for the full list of starting points.
            
            enddo ! indefinite do loop
            !enddo ! i
            
            ! Calling all heads of amoeba to stop subsequent evaluation. 7-18-2017
            do i = 1, nslaves
                trial = -1 ! 7-18-2017 Hard-wired for terminating the amoeba calculation.
                slave = i
                selected_input = -1._wp
                call sendjob(trial, slave, selected_input, 'head_group') ! 7-18-2017 checked.
            enddo
            
        else ! my_id != 0
            if(HEAD_WORLD/=MPI_COMM_NULL)then ! THE AMOEBA HEAD BLOCK (communicate only with the heads of amoeba, "including" the system root.) 7-11-2017
                do 
                    msgtype = 1
                    
                    call MPI_RECV( trial, 1, MPI_INTEGER, 0, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    if( trial == -1 )then ! 7-18-2017 Let all the amoeba members know stopping any further calculation.
                        call MPI_BCAST( trial, 1, MPI_INTEGER, 0, AMOEBA_WORLD, ierr )    
                        exit ! 7-18-2017 Important.
                    endif
                    
                    call MPI_RECV( selected_input, ndim, MPI_DOUBLE_PRECISION, 0, &
                        & msgtype, HEAD_WORLD, status, ierr )
                    
                    if(printout13) write(my_id+1001,'(a,i3," vertex ",<ndim>f12.7)') 'Step 1 trial ', trial, selected_input ! 7-23-2017 checked.
                    
                    msgtype = 2 ! Tag 2 is used for pass-receiving results from respective amoeba head.
                    call MPI_BCAST( trial, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )
                    call MPI_BCAST( selected_input, ndim, MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err ) ! 7-23-2017 maks sense.
                    
                    call amoeba_algorithm( bestvertex, result, bestobjval ) ! 7-18-2017 The major part of the amoeba algorithm. Be sure to refresh the model parameter input. [Li-Pin Juan no.2]
                    
                    !! 7-19-2017 Send The local optima back to the system root from the amoeba heads.
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
                enddo
                
            else ! 7-18-2017 The slaves (i.e. excluding the amoeba head) in respective amoeba group.    
                
                ! Note that we prevent amoeba members except their head to communicate with the system root.
                ! By doing so, the starting point assignement is carried out only between the system root and the head of each amoeba group.
                do 
                    call MPI_BCAST( trial, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )    
                    IF( trial == -1) exit
                    call MPI_BCAST( selected_input, ndim, MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err )
                    
                    call amoeba_algorithm( bestvertex, result, bestobjval ) ! # Need imporvement. [Li-Pin Juan no.2]
                enddo
            endif ! HEAD_WORLD/=MPI_COMM_NULL
            
        endif ! my_id == 0 or not.
            
        !! Important!!
        !if(irow==1) call read_short_list(listnumber,startpoint)
        !refpts = startpoint(15:24)
        !print*, ' debug ', ncol, num_procs, nrow ! successful command: mpiexec -np 12 (or 23) MPI_sandbox.exe (for the case of nrow == 11; 10 dimension and one more node for head.
        
        deallocate( refpts, selected_input, result, bestvertex )
        
        if( AMOEBA_WORLD/=MPI_COMM_NULL .or. my_id == 0 )then
            close(my_id+1001)
        endif
        
    elseif( mpi_exercise_mode==3 )then    
        
    endif ! mpi_exercise_model
    
    !call search_equilibrium(exit_log1) ! <===== replace solve_model() with this one. 3.10.2017 This is the working one. Obsolete, 7-3-2017.
    
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
    ! NOTE: the system root never calls amoeba_algorithm, so don't use any MPI_BARRIER with MPI_COMM_WORLD; otherwise, the program would stuck where your MPI_BARRIER is placed (due to the absence of the system root).    
    subroutine amoeba_algorithm( bestvertex, result, bestobjval )
        implicit none
        real(wp), dimension(:), intent(out) :: bestvertex, result
        real(wp), intent(out) :: bestobjval
        integer :: ndim, num_vertices, amo_msgtype, major_counter, good_msg_count
        integer :: shrink_counter, subworst_idx, i, j
        real(wp) :: plus_deviation, minus_deviation, true_deviation, sim_objval, best_val
        real(wp), dimension(:), allocatable :: sim_moments, ray_objval, dup_ray_objval, centroid, best_posi
        real(wp), dimension(:), allocatable :: worst_valvec, less_worse_valvec
        real(wp), dimension(:,:), allocatable :: vertex_list, mat_moments, worst_mat
        integer, dimension(:), allocatable :: rankinglist, ray_modelmsg ! Bug.
        
        ndim = size(bestvertex)
        num_vertices = ndim + 1
        
        allocate( vertex_list(ndim,ndim+1), sim_moments(ndim), mat_moments(ndim,ndim+1), ray_objval(ndim+1), ray_modelmsg(ndim+1), rankinglist(ndim+1) )
        allocate( dup_ray_objval(ndim+1), centroid(ndim), best_posi(ndim), worst_mat(ndim,noamoeba) )
        allocate( worst_valvec(noamoeba), less_worse_valvec(noamoeba) )
        
        ! Initialization
        bestvertex = 0._wp
        result = 0._wp
        bestobjval = 0._wp
        
        amo_msgtype = 0
        major_counter = 0 
        shrink_counter = 0
        
        do while( major_counter<=amoitrcrt )
            
            if( amo_msgtype==0 )then ! Initialization of additional "ndim" vertex position. 7-20-2017
                
                major_counter = major_counter + 1 ! The result is one.
                
                if(major_counter==1 .and. printout13)then
                    if(printout13) write(*,'(a,i4,a,i4,a,i4)') "- Amoeba activity - My_id ", my_id, " Amoeba_id ", amoeba_id, " Contract_id", contract_id ! 7-21-2017 OK.
                    call MPI_BARRIER( AMOEBA_WORLD, mpi_err )
                    if(printout13) write(*,'(a,i3,a,<ndim>f12.7)') "- Same ini. pt - amo_no ", amoeba_id, " pos ", selected_input ! 7-23-2017 checked.
                    call MPI_BARRIER( AMOEBA_WORLD, mpi_err )
                endif
                
                call MPI_BCAST( major_counter, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )                
                !call MPI_BCAST( major_counter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err ) ! Bug. Don't use it, because only members in the amoeba group                 
                
                if( CONTRACT_WORLD/=MPI_COMM_NULL )then
                    
                    if( CONTRACT_ID/=0 )then ! 7-23-2017 bug exists. Not disinsection yet.
                        plus_deviation = ( range_guess(CONTRACT_ID,2)-selected_input(CONTRACT_ID) )
                        minus_deviation = ( range_guess(CONTRACT_ID,1)-selected_input(CONTRACT_ID) )
                        ! Pick the largger deviation for creating the initial potion of every vertex. 7-20-2017
                        true_deviation = merge( plus_deviation, minus_deviation, abs(plus_deviation)>=abs(minus_deviation) )                            
                        selected_input(CONTRACT_ID) = selected_input(CONTRACT_ID) + stepsize/100._wp*true_deviation ! 7-20-2017 checked. It is direction sensitive design. Don't worry.
                        
                        !write(my_id+1001,'("mid", i3, "cid ", i3, " L ", f8.5, " H ", f8.5, " final ",f12.7, " minus ", f12.7, " positive ", f12.7, " true devi ", f12.7)') &
                        if(printout13)then
                            write(*,'("M_id", i3, " C_id ", i3, " LB ", f8.5, " HB ", f8.5, " Final ",f12.7, " Minus ", f12.7, " Positive ", f12.7, " Devi ", f12.7)') &
                                & my_id, contract_id, range_guess(CONTRACT_ID,1), range_guess(CONTRACT_ID,2), selected_input(CONTRACT_ID), minus_deviation, &
                                & plus_deviation, true_deviation
                        endif ! printout13
                    endif
                       
                    call MPI_ALLGATHER( selected_input, ndim, MPI_DOUBLE_PRECISION, vertex_list, & ! 7-21-2017 checked. 
                        & ndim, MPI_DOUBLE_PRECISION, CONTRACT_WORLD, mpi_err )                        
                                       
                    if( CONTRACT_ID == 0 )then ! checked! 7-22-2017
                        write( my_id+1001, '(/,a,i4,a)') "major counter ", major_counter, " List of generated initial vertices from amoeba member: "
                        do j = 1, ndim                                                       
                            write( my_id+1001, '(a, i3, x, <ndim+1>f12.7)') " row ", j, vertex_list(j,:) 
                        enddo
                    endif
                    
                endif ! CONTRACT_WORLD/=MPI_COMM_NULL
                
                if( CONTRACT_ID == 0 ) amo_msgtype = 1 ! move to eval & sort.
                call MPI_BCAST( amo_msgtype, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )     
                
            endif ! amo_msgtype == 0 initial allocation.
            
            if( amo_msgtype == 1 )then ! Evaluation and Sorting by contract group.
                
                ! Evaluation with CONTRACT_GROUP
                if( CONTRACT_WORLD/=MPI_COMM_NULL )then
                    
                    major_counter = major_counter + 1
                    
                    selected_input = vertex_list(:,CONTRACT_ID+1)
                    write(my_id+1001,'("Round no. ", i3, " Input ", <ndim>f12.5)') major_counter, selected_input

                    ! Core computation
                    modelmsg = 0
                    call read_parameter_model(para, '_1parameter.txt')
                    call search_equilibrium( selected_input, sim_moments, sim_objval, my_id, trial, modelmsg )         
                    
                    if(printout13)then
                        write(my_id+1001,'("My_id ", i3, " contract_id ", i3)') my_id, contract_id 
                        write(my_id+1001,'(" input  ", <ndim>f12.7)') selected_input ! 7-23-2017 check for input. Correct!
                        write(my_id+1001,'(" moment ", <ndim>f12.7 )') sim_moments
                        write(my_id+1001,'(" obj ", f12.7, " msg ", i3)') sim_objval, modelmsg 
                    endif ! printout13
                    
                    ! Communication after model computation
                    call MPI_ALLGATHER( sim_objval, 1, MPI_DOUBLE_PRECISION, ray_objval, & ! Product: ray_objval
                        & 1, MPI_DOUBLE_PRECISION, CONTRACT_WORLD, mpi_err )
                    
                    call MPI_ALLGATHER( modelmsg, 1, MPI_INTEGER, ray_modelmsg, & ! Product: ray_modelmsg
                        & 1, MPI_INTEGER, CONTRACT_WORLD, mpi_err )
                    
                    call MPI_ALLGATHER( sim_moments, ndim, MPI_DOUBLE_PRECISION, mat_moments, & ! Product: mat_moments
                        & ndim, MPI_DOUBLE_PRECISION, CONTRACT_WORLD, mpi_err )           
                    
                    ! BEFORE SORTING ! 7-24-2017 Double checked with other output.
                    if( printout13 .and. contract_id==0 )then ! 7-23-2017 checked. 
                        do i = 1, ndim+1
                            write(my_id+1001,'("nonsort-val ", f12.7, " msg ", i10, " mom ", <ndim>f12.7, " input ", <ndim>f12.7)')  ray_objval(i), ray_modelmsg(i), mat_moments(:,i), vertex_list(:,i)
                        enddo
                    endif

                    ! SORTING
                    rankinglist = [(i,i=1,ndim+1)]
                    dup_ray_objval = ray_objval
                    call sort_real( ray_objval, rankinglist ) ! sorted ray_objval (Penalty)
                    ray_modelmsg  = ray_modelmsg(rankinglist)  ! sorted ray_modelmsg (Output logical array))
                    do j = 1, ndim
                        vertex_list(j,:) = vertex_list(j,rankinglist) ! sorted vertex_list (Trial vertex)
                        mat_moments(j,:) = mat_moments(j,rankinglist) ! sorted mat_moments (Simmulated moments))
                    enddo

                    good_msg_count = count( ray_modelmsg==0 ) ! 0 means the triral vertex leads to a valid output.
                    
                    if( good_msg_count == ndim+1 )then ! every vertex returns a valid result.
                        write(my_id+1001,'(/,a)') '--- amo_msgtype = 2 - enter AMOEBA'
                        amo_msgtype = 4 ! enter amoeba stage        
                    elseif( 1<good_msg_count .and. good_msg_count<ndim+1 )then
                        write(my_id+1001,'(/,a)') '--- amo_msgtype = 3 - enter SHRINKAGE'
                        amo_msgtype = 3 ! shrinkage by contract group.
                    else
                        write(my_id+1001,'(/,a)') '--- amo_msgtype = 9 - Exit '
                        amo_msgtype = 9 ! return and exit            
                    endif 
                endif ! CONTRACT_WORLD/=MPI_COMM_NULL
                
                call MPI_BCAST( amo_msgtype, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )
                call MPI_BCAST( major_counter, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )
                call MPI_BCAST( shrink_counter, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )
                call MPI_BCAST( ray_modelmsg, ndim+1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )
                call MPI_BCAST( ray_objval, ndim+1, MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err )
                call MPI_BCAST( mat_moments, ndim*(ndim+1), MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err )
                call MPI_BCAST( vertex_list, ndim*(ndim+1), MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err )
                
                if( printout13 .and. (contract_id==1 .or. contract_id==0) )then ! 7-23-2017 checked. ! 7-24-2017 Double checked with other output.
                    do i = 1, ndim+1
                        write(my_id+1001,'("sorted-VAL ", f12.7, i3, " mom ", <ndim>f12.7)')  ray_objval(i), ray_modelmsg(i), mat_moments(:,i)
                    enddo
                endif                
                
            endif ! amo_msgtype == 1
            
            if( amo_msgtype == 3 )then ! shrinkage toward the "best point" only on those invalid vertices of "contract group."
                if( CONTRACT_WORLD/=MPI_COMM_NULL )then
                    shrink_counter = shrink_counter + 1
                    major_counter = major_counter + 1
                
                    do i = 1, ndim ! Only shrink those bad vertices toward the best point. ! 7-24-2017 Double checked with other output.
                        if( printout13 .and. contract_id == 0 .and. i == 1 )then
                            do j = 1, num_vertices
                                write(my_id+1001, '(i3)', advance='no') ray_modelmsg(j)
                            enddo
                            write(my_id+1001, '(a)') ' '
                            do j = 1, ndim+1
                                write(my_id+1001, '("Bshrink.rank ", i3, x,"msg ", i3, x, <num_vertices>f12.7)') j, ray_modelmsg(j), vertex_list(:,j)
                            enddo
                        endif ! if
                        
                        ! The shrinkage takes place in every node of the same amoeba. 7-23-2017 ! 7-24-2017 Double checked.
                        if( ray_modelmsg(i+1)==1 )then
                            
                            if(printout13 .and. contract_id == 0) write(my_id+1001,'("shrink takes place at vertex ", i3)') i+1
                            if(printout13 .and. contract_id == 0) write(my_id+1001,'(" best ", <ndim>f12.7, " to-be ", <ndim>f12.7)') vertex_list(:,1), vertex_list(:,i+1)
                            
                            ! CORE OF SHRINKAGE
                            vertex_list(:,i+1) = amotau*vertex_list(:,1) + (1._wp-amotau)*vertex_list(:,i+1) ! 7-23-2017 vertex_list(:,1) is the best vertex in the very first run of the computation.                              
                            
                            if(printout13 .and. contract_id == 0) write(my_id+1001,'(i3, x, "Post-shrink ",i3," vertex ",<ndim>f12.7)') i, ray_modelmsg(i+1), vertex_list(:,i+1)
                            
                        endif
                        
                    enddo ! i
                    
                    amo_msgtype = 1 ! Evaluation and Sorting by contract group.
                    
                endif ! CONTRACT_WORLD/=MPI_COMM_NULL     
                
                call MPI_BCAST( amo_msgtype, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err ) ! The head of amoeba_world is the same as that of contract world.
                call MPI_BCAST( major_counter, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )
                call MPI_BCAST( shrink_counter, 1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )
                call MPI_BCAST( ray_modelmsg, ndim+1, MPI_INTEGER, 0, AMOEBA_WORLD, mpi_err )
                call MPI_BCAST( ray_objval, ndim+1, MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err )
                call MPI_BCAST( mat_moments, ndim*(ndim+1), MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err )
                call MPI_BCAST( vertex_list, ndim*(ndim+1), MPI_DOUBLE_PRECISION, 0, AMOEBA_WORLD, mpi_err )                
                
            endif ! amo_msgtype == 3 ! checked! 7-24-2017
            
            if( amo_msgtype == 4 )then
                major_counter = major_counter + 1
                centroid = sum( vertex_list(:,1:(ndim+1)-noamoeba),dim=2 )/((ndim+1)-noamoeba)
                best_posi = vertex_list(:,1)
                best_val  = ray_objval(1) 
                
                do j = 1, noamoeba
                      subworst_idx = num_vertices-noamoeba+j
                      worst_mat(:,j) = vertex_list(:,subworst_idx)
                      worst_valvec(j) = ray_objval(subworst_idx)
                      less_worse_valvec(j) = ray_objval(subworst_idx-1)
                enddo
                
                if( AMOEBA_ID==0 )then
                    do j = 1, noamoeba
                                    
                    enddo ! j
                endif
                
            endif ! amo_msgtype = 4
            
            if( amo_msgtype == 2 )then ! shrinkage toward the centroid on those invalid vertices of "amoeba group."
                        
            endif ! amo_msgtype == 2            
   
        enddo ! while ( major_counter<=amoitrcrt )

        deallocate( vertex_list, sim_moments, mat_moments, ray_objval, ray_modelmsg, rankinglist, dup_ray_objval )
        deallocate( centroid, best_posi, worst_mat, worst_valvec, less_worse_valvec )
        
    end subroutine amoeba_algorithm
    
    ! The root sends trial parameter combination to slaves
    subroutine sendjob(trial,slave,parcel,opt)
        implicit none
        integer, intent(in) :: trial, slave
        real(wp), dimension(:), intent(in) :: parcel
        integer :: msgtype, n
        character(len=*), optional :: opt
        n = size(parcel)
        msgtype = 1 ! tag 1 is used for the communication of model parameter passing. tag 2 is used for result passing from slave.
        if(.not.present(opt))then
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
    
    subroutine read_short_list(filename,startpoint)
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