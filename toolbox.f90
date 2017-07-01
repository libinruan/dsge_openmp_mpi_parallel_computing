include 'mkl_vsl.f90'
    
module fmpi_header
    include 'mpif.h'
end module fmpi_header

module fmpi_variable
    integer :: MPI_ERR, NUM_PROCS, MY_ID, GROUP_ALL, MPI_PROVIDED = 0 
    integer, parameter :: MY_ROOT = 0
end module fmpi_variable

subroutine fmpi_init
    use fmpi_header
    use fmpi_variable
    implicit none
    call MPI_INIT_THREAD( MPI_THREAD_FUNNELED, MPI_PROVIDED, MPI_ERR ) ! ALL MPI CALLS MADE BY MASTER THREAD; whenever your program uses threading, you should use MPI_Init_thread().
    call MPI_COMM_RANK( MPI_COMM_WORLD, MY_ID, MPI_ERR ) ! MPI_COMM_WORLD is declared in "mpif.h"
    call MPI_COMM_SIZE( MPI_COMM_WORLD, NUM_PROCS, MPI_ERR )
    call MPI_COMM_GROUP( MPI_COMM_WORLD, GROUP_ALL, MPI_ERR ) ! get the base group from MPI_comm_world
end subroutine fmpi_init

module toolbox
    use omp_lib
    use fmpi_header
    use fmpi_variable
    use mkl_vsl_type ! For sobol sequence
    use mkl_vsl      ! For sobol sequence
    use f95_precision, ONLY: WP => DP
    use, intrinsic :: ieee_arithmetic ! for detecting NaN and infinite value.     
    implicit none    
    !  the shared variables used in Tauchen block.
    real(wp), dimension(:,:), allocatable :: p_tau
    real(wp), dimension(:), allocatable :: y_tau,s_tau
    real(wp) :: sigma_tau, rho_tau, smin_tau
    integer :: n_tau    
contains        
    subroutine sobol_sequence_transformed( origin, rng )
        implicit none
        real(wp), dimension(:,:), intent(in) :: rng
        real(wp), dimension(:,:), intent(inout) :: origin
        real(wp) :: mul, add
        integer :: nrow_seq, npara, i, j
        nrow_seq = size(origin,dim=1)
        npara = size(rng,dim=1)
        do i = 1, npara
            mul = rng(i,2) - rng(i,1)
            add = rng(i,1)
            do j = 1, nrow_seq
                origin(j,i) = origin(j,i)*mul + add
                !if(i==8) print*,rng(i,2),rng(i,1),mul,add
            enddo
        enddo
    end subroutine sobol_sequence_transformed    
    
    subroutine para_range( n1, n2, nprocs, irank, ista, iend )
        implicit none
        integer :: n1, n2, nprocs, irank, iwork1, iwork2
        integer, intent(out) :: ista, iend
        iwork1 = (n2-n1+1)/nprocs
        iwork2 = mod(n2-n1+1, nprocs)
        ista = irank*iwork1 + n1 + min(irank,iwork2)
        iend = ista + iwork1 - 1
        if(iwork2 > irank) iend = iend + 1        
    end subroutine para_range
    
    subroutine make_vertices( base, rng, list, delta )
        implicit none
        real(wp), dimension(:), intent(in) :: base
        real(wp), dimension(:,:), intent(in) :: rng
        real(wp), dimension(:,:), intent(out) :: list
        real(wp) :: lower_distance, upper_distance, dist
        integer :: ndim, i
        real(wp), optional, intent(in) :: delta
        ndim = size( base )
        list(1,:) = base
        ! choose the distance to the "farrest" endpoint as the step size.
        do i = 2, ndim+1 
            list(i,:) = base
            lower_distance = abs( list(i,i-1) - rng(i-1,1) )
            upper_distance = abs( list(i,i-1) - rng(i-1,2) )
            ! dist = merge( upper_distance, lower_distance, upper_distance<=lower_distance) 
            dist = merge( delta*upper_distance, -delta*lower_distance, upper_distance>=lower_distance) 
            list(i,i-1) = list(i,i-1) + dist 
        enddo
    end subroutine make_vertices    
    
    !subroutine bind_within_range(ori,rng,dev,msg)
    !    implicit none
    !    real(wp), dimension(:,:), intent(in) :: rng
    !    real(wp), dimension(:,:), intent(inout) :: ori
    !    real(wp), dimension(:), optional, intent(in) :: dev
    !    real(wp), dimension(:,:), allocatable :: orp
    !    logical, optional :: msg
    !    integer :: m, n, k, i
    !    m = size(ori,1)
    !    n = size(ori,2)
    !    k = size(rng,1)
    !    msg = .false.
    !    allocate( orp(m,n) )
    !    orp = ori 
    !    if(n/=k) stop 'Error in ''bind_within_range'''
    !    do i = 1, n
    !        ori(:,i) = merge(rng(i,1),ori(:,i),rng(i,1)>ori(:,i))    
    !        ori(:,i) = merge(ori(:,i)-2._wp*dev(i),ori(:,i),rng(i,2)<ori(:,i)) ! Pull it backward to the range because we modify the best vertex by moving one of its coordinate "upward"   
    !    enddo
    !    if( present(msg) ) then
    !        if(sum(abs(ori-orp))>0._wp) msg = .true.    
    !    endif
    !    deallocate( orp )
    !end subroutine bind_within_range

    subroutine amoeba_try( centroid,worst,opt,alpha,gamma,beta, amotry_vec )     
        ! alpha=1, gamma=1.0, beta=0.5 ! Lee and Wiswall, page 8, A parallel Implementation of the Simplex ... routine.
        implicit none
        real(wp), dimension(:), intent(in) :: centroid, worst
        character(len=*), intent(in) :: opt
        real(wp) :: alpha, gamma, beta
        ! real(wp), dimension(:), allocatable :: amotry_vec! Bug: you cannot use allocatable on input
        real(wp), dimension(:), intent(out) :: amotry_vec
        real(wp), dimension(:), allocatable :: reflect
        integer :: n
        n = size( centroid )
        allocate( reflect(n) )
        reflect  = centroid + alpha*( centroid - worst )
        select case(opt)
            case('r') ! reflection point
                amotry_vec = reflect
            case('e') ! extension point
                amotry_vec = reflect + gamma*( reflect - centroid )
            case('cr') ! contraction point with reflection point
                amotry_vec = beta*( centroid + reflect )
            case('co') ! contraction point with original point
                amotry_vec = beta*( centroid + worst )     
        endselect    
        deallocate( reflect )
    end subroutine amoeba_try
    
    !function dis_best_worst(ray, idx)
    !    implicit none
    !    real(wp), dimension(:),intent(in) :: ray
    !    real(wp) :: dis_best_worst
    !    real(wp), dimension(:), allocatable :: temp
    !    integer, dimension(:), allocatable :: key
    !    integer, intent(out) :: idx
    !    integer :: n, i
    !    logical :: info1
    !    n = size(ray)
    !    allocate(key(n), temp(n))
    !    key = [(i,i=1,n)]
    !    temp = ray
    !    call dlasrt2('I', n, temp, key, info1)
    !    dis_best_worst = abs(temp(1)-temp(n))
    !    idx = key(1)
    !    deallocate(key, temp)        
    !end function dis_best_worst
    
    function brent(ax,bx,cx,func,xmin)
        ! look for the minimizor XMIN for the penalty function FUNC.
        ! Note that BX is between AX and CX.
        implicit none
        real(wp), intent(in) :: ax,bx,cx
        real(wp), intent(out) :: xmin
        real(wp) :: brent
        interface
          function func(x)
          use f95_precision, only: wp => dp
          implicit none
          real(wp), intent(in) :: x
          real(wp) :: func
          end function func
        end interface
        integer, parameter :: itmax=100
        real(wp), parameter :: tol=sqrt(epsilon(ax)),cgold=0.381966011250105_wp,zeps=1.0e-3_wp*epsilon(ax)
        integer :: iter
        real(wp) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
        a=min(ax,cx)
        b=max(ax,cx)
        v=bx
        w=v
        x=v
        e=0.0
        fx=func(x)
        fv=fx
        fw=fx
        do iter=1,itmax
          xm=0.5_wp*(a+b)
          tol1=tol*abs(x)+zeps
          tol2=2.0_wp*tol1
          if (abs(x-xm) <= (tol2-0.5_wp*(b-a))) then
            xmin=x
            brent=fx
            return
          end if
          if (abs(e) > tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.0_wp*(q-r)
            if (q > 0.0) p=-p
            q=abs(q)
            etemp=e
            e=d
            if (abs(p) >= abs(0.5_wp*q*etemp) .or. &
              p <= q*(a-x) .or. p >= q*(b-x)) then
              e=merge(a-x,b-x, x >= xm )
              d=cgold*e
            else
              d=p/q
              u=x+d
              if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
            end if
          else
            e=merge(a-x,b-x, x >= xm )
            d=cgold*e
          end if
          u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
          fu=func(u)
          if (fu <= fx) then
            if (u >= x) then
              a=x
            else
              b=x
            end if
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
          else
            if (u < x) then
              a=u
            else
              b=u
            end if
            if (fu <= fw .or. w == x) then
              v=w
              fv=fw
              w=u
              fw=fu
            else if (fu <= fv .or. v == x .or. v == w) then
              v=u
              fv=fu
            end if
          end if
        end do
        stop 'brent: exceed maximum iterations'
    end function brent

    subroutine tauchen(rho,sigma,p,y,s,spread) ! V021415.
    ! Approximate an univariate AR(1) process by Markov chain
    !       y(t) = rho y(t-1)+ sigma sqrt(1-rho^2) e(t), where e(t)~N(0,1)
    ! INPUTS: rho - serial correlation coefficient
    !         sigma - coefficient of variation (Li: sigma^2 is the variance of y(t). See reference Kitao, page 50)
    ! OUTPUTS: p : n-by-n matrix of Markov transition probabilities
    !          y : n-by-1 vector of symmetric and evenly-spaced Markov state space
    !          s : n-by-1 vector of stationary distribution   
    ! Note that sigma_tau, rho_tau, p_tau, y_tau, s_tau are global variables defined in the declaration area of module TOOLBOX.
        implicit none
        real(wp), intent(in) :: rho,sigma
        real(wp), intent(in), optional :: spread
        real(wp), dimension(:,:) :: p
        real(wp), dimension(:) :: y,s
        integer :: n
        real(wp) :: smin,emin
        sigma_tau=sigma
        rho_tau=rho
        n=size(y)
        n_tau=n
        if(present(spread))then
            call tauch(spread)
        else ! Note that err is a function defined below
            emin=brent(1._wp,2.5_wp,4._wp,err,smin) ! Use BRENT to look for the optimal spread SMIN for Tauchen's markove chain matrix.
            !write(*,'(2(a,f8.5))') 'brent on tauchen parameter, error: ', emin, ' smin: ', smin
            call tauch(smin)
            smin_tau = smin ! info to be kept
        endif
        p=p_tau
        y=y_tau
        s=s_tau ! statinoary distribution obtained from p_tau using the ERGODIC subroutine
    end subroutine tauchen

    subroutine tauch(m) ! V021415.
    ! following Tauchen's method to calculate the transition matrix, levels of state, and stationary distribution
    ! Note: this subroutine invokes global variables (defined in TOOLBOX): p_tau, y_tau, s_tau, rho_tau, sigma_tau 
    !       from subroutine TAUCHEN.
        implicit none
        real(wp), intent(in) :: m
        real(wp) :: ybar,dy,sr
        real(wp), dimension(size(y_tau)) :: yd
        integer :: i
        ybar = m*sigma_tau
        call grid(y_tau,-ybar,ybar,1._wp)
        dy = ybar/(n_tau-1) ! n_tau is initialized in sburoutine 'tauchen'
        sr = sigma_tau*sqrt(1-rho_tau**2)
        do i = 1,n_tau
            yd = (y_tau-rho_tau*y_tau(i)+dy)/sr ! Note that the starting state is y_tau(i)
            call vdcdfnorm(n_tau,yd,p_tau(i,:)) ! Normal distribution. row levels represent the starting states; column levels represent the ending states.
        enddo
        p_tau(:,n_tau)=1
        do i=n_tau,2,-1
            p_tau(:,i) = p_tau(:,i) - p_tau(:,i-1)    
        enddo
        call ergodic(p_tau,s_tau)
    end subroutine tauch
    
    function err(m) ! V021415.
    ! Note: this function invokes global variables (defined in TOOLBOX): p_tau, y_tau, s_tau, rho_tau, sigma_tau 
    ! ____from subroutine TAUCH which inherits these parameters from subroutine TAUCHEN.
        implicit none
        real(wp), intent(in) :: m
        real(wp) :: err, rho_, sigma_
        call tauch(m)
        call markovtest(p_tau,y_tau,s_tau,rho_,sigma_) ! generate rho_ and sigma_ based on p_tau, y_tau, and s_tau.
        err = (log(1-rho_)/log(1-rho_tau)-1)**2+(sigma_/sigma_tau-1)**2
        ! write(*,'(a,f8.5)') 'err: ', err
    end function err

    subroutine ergodic(p,s)
        use lapack95, only: geev ! Computes the eigenvalues and left and right eigenvectors of a general matrix.
        implicit none
        real(wp), dimension(:,:), intent(in) :: p
        real(wp), dimension(:), intent(out) :: s
        real(wp), dimension(size(s),size(s)) :: ip,vl
        real(wp), dimension(size(s)) :: wr,wi,dw
        integer :: m,info,uw,w1(1)
        real(wp) :: ds
        m = size(s)
        ip=p
        call geev(ip,wr,wi,vl=vl,info=info) ! vl: the eigen vector.
        dw=abs(sqrt(wr*wr+wi*wi)-1)
        w1=minloc(dw) ! looks for the smallest eigen value.
        uw=count(dw<1000*epsilon(dw)) ! epsilon(X)=E: returns the smallest number E of the same kind as X such that 1 + E > 1.
        if(uw<1) print'(a)','Warning: No unitary eigenvalue is found. Stationary distribution of Markov chain does not exist.'
        if(uw>1) print'(a)','Warning: More than one unitary eigenvalue is found. Stationary distribution of Markov chain is not unique.'
        if(uw<1 .or. uw>1) print'(a,f8.5,a,f8.5)','Using eigenvalue ',wr(w1(1)),'+i',wi(w1(1))
        s=vl(:,w1(1))/sum(vl(:,w1(1))) ! normalize the eigen vector.
        if(any(s<0)) then ! In the following code block, all that it does is to evenly redistribute the mass of negative element to those non-negative.
            print '(a)','The stationary distribution of Markvo chain has negative values. Rebalancing...'
            ds=sum(s,mask=s<0)/count(s>=0)
            where(s<0) 
                s=0
            elsewhere
                s=s+ds
            endwhere            
        endif        
    end subroutine ergodic

    subroutine markovtest(p,y,s,rho,sigma)
        use blas95, only: gemv ! Computes a matrix-vector product using a general matrix
        implicit none
        real(wp), dimension(:,:), intent(in) :: p
        real(wp), dimension(:), intent(in) :: y,s
        real(wp), intent(out) :: rho, sigma
        real(wp), dimension(size(y)) :: py
        real(wp) :: eyy, ey2, e2y
        integer :: n
        n=size(y)
        call gemv(p,y,py)
        eyy = sum(s*y*py) ! that is, the sum (expectation w.r.t. s) of component-wiese product of y at t and y at t+1
        ey2 = sum(s*y**2) ! the expectation of y^2 itself.
        e2y = sum(s*y)**2 ! the square of the expectation of y
        rho = (eyy-e2y)/(ey2-e2y) ! calculate coefficient of correlation
        sigma=sqrt(ey2-e2y) ! calculate standard deviation
    end subroutine markovtest

    subroutine ucase(string)
    ! Purpose: shift a character string to upper case on any processor
        implicit none
        character(len=*), intent(inout) :: string
        integer :: i,length
        length=len(string)
        do i=1,length
            if(lge(string(i:i),'a') .and. lle(string(i:i),'z')) then !.True. if lexically A follows B in lge(A,B)
                string(i:i)=achar(iachar(string(i:i))-32) !"iachar" gets the position in the ASCII character set; it is the inverse of "achar".
            endif
        enddo
    end subroutine ucase

    subroutine grid(x,xmin,xmax,s)
    ! Purpose: Generate grid x on [xmin,xmax] using spacing parameter s
    ! s=1       linear spacing
    ! s>1       skewed-to-the-right grid spacing with power s [most of the data are on the left side of the distribution]
    ! 0<s<1     skewed-to-the-left grid spacing with power s [most of the data are on the right side of the distribution]
    ! s<0       geometric spacing with distance changing by a factor -s^(1/(n-1))
    ! s=-1      logarithmic spacing with distances changing by a factor (xmax-xmin+1)^(1/(n-1))
    ! s=0       logarithmic spacing with distances changing by a factor (xmax/xmin)^(1/(n-1)), only if xmax,xmin>0
        implicit none
        real(wp), dimension(:), intent(out) :: x
        real(wp), intent(in) :: xmin, xmax, s
        real(wp) :: c ! growth rate of grid subintervals for logarithmic spacing
        integer :: n,i
        n=size(x)
        forall(i=1:n) x(i)=(i-1)/real(n-1,wp)
        if(s>0.0_wp) then
            x=x**s*(xmax-xmin)+xmin
        else
            if(s==-1.0_wp) then
                c=xmax-xmin+1
            else
                c=-s
            endif
            x=((xmax-xmin)/(c-1))*(c**x)-((xmax-c*xmin)/(c-1))
        endif
    end subroutine grid

    subroutine tool(funct,input,output)
        implicit none
        interface
           function funct(x)
                implicit none
                integer :: x,funct
            end function funct
        end interface
        integer :: input,output
        output = funct(input)
    end subroutine tool

    subroutine sm(x,y,leng,deci)
    ! x: matrix to be stored; y: the file name without filename extension, leng: length of a displayed number, deci: number of decimal spaces
        implicit none
        real(wp), dimension(:,:), intent(in) :: x
        character(len=*), intent(in) :: y
        integer, optional :: leng, deci
        integer :: a, b
        character(len=40) :: fname, sname
        integer :: m,n,i,j
        m = size(x,dim=1)
        n = size(x,dim=2)
        fname = y // '.txt'
        if( present(leng) )then
            a = leng
            b = deci
        else
            a = 9
            b = 4
        endif      
        write(sname,'(a,i8,a,i2.2,a,i2.2,a)') '(', n, 'f', a, '.', b, '))'
        open(unit=0,file=fname,action="write",status="replace",recl=(a*n+10)) ! CHANGE THE RECL PARAMETER FOR EXAMPEL, F15.4 THEN USE 16*N+10
        do i = 1,m
            write(unit=0,fmt=sname) (x(i,j),j=1,n)    
        enddo 
        close(0)
    end subroutine sm  
    
    subroutine sm_science(x,y,leng,deci)
    ! x: matrix to be stored; y: the file name without filename extension, leng: length of a displayed number, deci: number of decimal spaces
        implicit none
        real(wp), dimension(:,:), intent(in) :: x
        character(len=*), intent(in) :: y ! TRICK
        integer :: leng, deci
        character(len=40) :: fname, sname
        integer :: m,n,i,j
        m = size(x,dim=1)
        n = size(x,dim=2)
        fname = y // '.txt'
        write(sname,'(a,i8,a,i2.2,a,i2.2,a)') '(', n, 'ES', leng, '.', deci, '))'
        open(unit=0,file=fname,action="write",status="replace",recl=(leng*n+10)) ! CHANGE THE RECL PARAMETER FOR EXAMPEL, F15.4 THEN USE 16*N+10
        do i = 1,m
            write(unit=0,fmt=sname) (x(i,j),j=1,n)    
        enddo 
        close(0)
    end subroutine sm_science      
    
    subroutine smi(x,y,opt)
    ! x: matrix to be stored; y: the file name without filename extension
        implicit none
        integer, dimension(:,:), intent(in) :: x
        character(len=*), intent(in) :: y
        character(len=40) :: fmt1, fname
        integer, optional :: opt
        integer :: m,n,i,j
        m  = size(x,dim=1)
        n  = size(x,dim=2)
        fname = y // '.txt'
        if(.not.present(opt))then
            write(fmt1,'(a,i8,a)') '(',n,'(i6))'
        else
            write(fmt1,'(a,i8,a,i8,a)') '(',n,'(i',opt,'))'
        endif
        open(unit=0,file=fname,action="write",status="replace",recl=(7*n+10))
        do i = 1,m
            write(unit=0,fmt=fmt1) (x(i,j),j=1,n)    
        enddo 
        close(0)
    end subroutine smi      
    
    subroutine ss(x,y,leng,deci) 
    ! save sequence x in .txt format with filename y
    ! e.g. call ss(mat(:,1),'mat')
    ! or
    ! write(test,'(i6,a)') na,'mat'  !where na is an integer
    ! call ss(mat(1,:),trim(test)) 
        real(wp), dimension(:), intent(in) :: x
        integer, optional, intent(in) :: leng, deci
        integer :: a, b
        character(len=*), intent(in) :: y ! arbitrary length input string
        character(len=40) :: sname, fname     ! string for format and filename
        integer :: n, i
        n=size(x)! # of columns in x
        fname=y//'.txt'  ! concatenate
        if(present(leng).and.present(deci))then
            a = leng
            b = deci
        else
            a = 9
            b = 4
        endif
        write(sname,'(a,i2.2,a,i2.2,a)') '(f', a, '.', b, ')'
        !(above) integer to character: &
        !  setup the output format per line
        open(unit=0,file=fname,action="write",status='replace',recl=(a+10))               ! loop to save each row &
        do i=1,n
            write( unit=0,fmt=sname ) x(i)            !   in order
        enddo    
        close(0)                               ! 
    end subroutine ss

    subroutine ssi(x,y) ! save sequence x in .txt format with filename y
    !e.g. call ss(mat(:,1),'mat')
    ! or
    ! write(test,'(i6,a)') na,'mat'  !where na is an integer
    ! call ss(mat(1,:),trim(test)) 
        integer, dimension(:), intent(in) :: x
        character(len=*), intent(in) :: y ! arbitrary length input string
        character(len=40) :: fmt1,str     ! string for format and filename
        integer :: n,i
        n=size(x)! # of columns in x
        str=y//'.txt'  ! concatenate
        write(fmt1,'(a)') '(i6)' 
        !(above) integer to character: &
        !  setup the output format per line
        open(unit=0,file=str,action="write",status='replace',recl=(16+10))               ! loop to save each row &
        do i=1,n
           write(unit=0,fmt=fmt1) x(i)         !   in order
        end do
        close(0)                            
    end subroutine ssi  

    subroutine rouwenhorst(rho,sigma,p,y,s)
    ! Rouwenhorst method to approximate univariate AR(1) process by an n-state Markov chain
    !   z(t) = rho z(t-1) + sigma sqrt(1-rho^2) e(t), e(t)~N(0,1)
    ! Inputs: rho - serial correlation coefficient (i.e., first-order autocorrelation, http://www.stat.yale.edu/~jtc5/251/stochastic-processes.pdf),
    !         sigma - coefficient of standard variation of z(t)
    ! Outputs: p is an n-by-n matrix of Markov transition probabilities
    !          y is an n-by-1 vector of symmetric and evenly-spaced Markov state space
    !          s is an n-by-1 vector of stationary distribution (binomial)
    ! Reference: http://karenkopecky.net/Rouwenhorst_WP.pdf pages 6 to '9'
    !            See page 9: a. Rouwenhort's {y} is particularly apt for approximating Gaussian AR(1) processes
    !                        b. can always the unconditional moment of {z}           
    !                        c. can be applied to any stationary AR(1) process, even with high persistence. 
        implicit none
        real(wp), intent(in) :: rho, sigma
        real(wp), dimension(:,:), intent(out) :: p
        real(wp), dimension(:), intent(out) :: y,s
        real(wp) :: ybar, q
        integer :: n 
        n = size(y)
        n_tau=n
        if(size(p,dim=1)/=n .or. size(p,dim=2)/=n) then
            print '(a,i3,a,i3)','[Error] Rouwenhorst: p must be a square matrix of size ',n,' x ',n
            stop 'Program terminated by Rouwenhorst'
        endif
        if(size(s)/=n) then
            print '(a,i3)','Rouwenhorst: y and s must be vectors of the same size ',n
            stop 'Program terminated by Rouwenhorst'
        endif
        ybar = sigma*sqrt(real(n-1,wp))
        q = (1+rho)/2
        call rhmat(p,q) 
        call grid(y,-ybar,ybar,1.0_wp) 
        call binom(s) 
    !contains
    end subroutine rouwenhorst
            
    recursive subroutine rhmat(p,q)
    !rhmat: rouwenhorst matrix
        implicit none
        real(wp), dimension(:,:), intent(out) :: p
        real(wp), intent(in) :: q
        real(wp), dimension(size(p,dim=1)-1,size(p,dim=2)-1) :: p1
        integer :: h
        h = size(p,dim=1)
        if(size(p,dim=2)/=h) stop 'P must be a square matrix'
        if(h<2) stop 'P must be at least 2-by-2 matrix'
        if(h==2) then
            p = reshape((/q,1-q,1-q,q/),(/2,2/))
        else
            call rhmat(p1,q)
            !http://som.yale.edu/~geert/A%20General%20Equilibrium%20Approach%2072406.pdf, page 327
            p=0
            p(1:h-1,1:h-1)=q*p1
            p(1:h-1,2:h)=(1-q)*p1+p(1:h-1,2:h)
            p(2:h,1:h-1)=(1-q)*p1+p(2:h,1:h-1)
            p(2:h,2:h)=q*p1+p(2:h,2:h)
            p(2:h-1,:)=p(2:h-1,:)/2
        endif
    end subroutine rhmat

    subroutine binom(f)
    !binom: Binomial probability mass function with p = 1/2
    !Reference: http://som.yale.edu/~geert/A%20General%20Equilibrium%20Approach%2072406.pdf, page 328
        implicit none
        real(wp), dimension(:), intent(out) :: f
        integer :: n,k
        n=size(f)
        f(1)=2.0_wp**(1-n)
        do k=1,n_tau-1
            f(k+1)=f(k)*(n-k)/k
        enddo
    end subroutine binom

    subroutine read_matrix(mat,fname) ! double-precision floating point version
	    implicit none
	    integer :: n,j,iostat
	    real(wp), dimension(:,:), intent(inout) :: mat
	    real(wp), dimension(:), allocatable :: temp
        character(len=*) :: fname
	    character(len=40) :: fmt1
	    j=size(mat,dim=2)
	    allocate(temp(j))
	    open(unit=10,file=fname,status='old',action='read',iostat=iostat) 
	    n=0
	    if(iostat==0) then
		    do
			    read(unit=10,fmt=*,iostat=iostat) temp
			    if(iostat/=0) exit
			    n=n+1
			    mat(n,:)=temp(:)
		    enddo		
	    endif
	    deallocate(temp)
        close(10)
    end subroutine read_matrix    

    subroutine read_matrixi(mat,fname) ! floating point version
	    implicit none
	    integer :: n,j,iostat
	    integer, dimension(:,:), intent(inout) :: mat
	    integer, dimension(:), allocatable :: temp
        character(len=*) :: fname
	    character(len=40) :: fmt1
	    j=size(mat,dim=2)
	    allocate(temp(j))
	    open(unit=10,file=fname,status='old',action='read',iostat=iostat) 
	    n=0
	    if(iostat==0) then
		    do
			    read(unit=10,fmt=*,iostat=iostat) temp
                ! print*,temp
			    if(iostat/=0) exit
			    n=n+1
			    mat(n,1:j)=temp(1:j)
		    enddo		
	    endif
	    deallocate(temp)
        close(10)
    end subroutine read_matrixi

    subroutine spline(x,y,y2,yp1,ypn)
   	    use lapack95, only: pttrf, pttrs
   	    implicit none
   	    real(wp), dimension(:), intent(in) :: x,y
   	    real(wp), dimension(:), intent(out) :: y2
   	    real(wp), intent(in), optional :: yp1,ypn
   	    real(wp), dimension(size(x)) :: b
   	    real(wp), dimension(size(x)-1) :: a
   	    integer :: n,info
   	    n=size(x)
   	    if (size(y)/=n .or. size(y2)/=n) then
   	    	print *, 'spline: x,y and y2 must be of the same size'
   	    	stop 'program terminated by spline'
   	    end if
   	    a=x(2:n)-x(1:n-1)
   	    b(1)=a(1)
   	    b(n)=a(n-1)
   	    b(2:n-1)=x(3:n)-x(1:n-2)
   	    b=2.0_wp*b
   	    y2(1:n-1)=(y(2:n)-y(1:n-1))/a
   	    if (present(ypn)) then
   	    	y2(n)=ypn
   	    else
   	    	y2(n)=y2(n-1)
   	    	a(n-1)=0.0_wp
   	    end if
   	    y2(2:n)=y2(2:n)-y2(1:n-1)
   	    if (present(yp1)) then
   	    	y2(1)=y2(1)-yp1
   	    else
   	    	y2(1)=0.0_wp
   	    	a(1)=0.0_wp
   	    end if
   	    y2=6.0_wp*y2
   	    call pttrf(b,a,info)
   	    call pttrs(b,a,y2,info)
    end subroutine spline
    
    function splint(x,y,y2,xi)
   	    implicit none
   	    real(wp), dimension(:), intent(in) :: x,y,y2
   	    real(wp), intent(in) :: xi
   	    real(wp) :: splint
   	    real(wp) :: a,b,d
   	    integer :: n,i
   	    n=size(x)
   	    if (size(y)/=n .or. size(y2)/=n) then
   	    	print *, 'splint: x,y and y2 must be of the same size'
   	    	stop 'program terminated by splint'
   	    end if
   	    i=max(min(locate(x,xi),n-1),1)
   	    d=x(i+1)-x(i)
   	    if (d == 0.0) stop 'bad x input in splint'
   	    a=(x(i+1)-xi)/d
   	    b=(xi-x(i))/d
   	    splint=a*y(i)+b*y(i+1)+((a**3-a)*y2(i)+(b**3-b)*y2(i+1))*(d**2)/6.0_wp
    end function splint

    ! Locate the LOWER BOUND of the interval that encloses the point "x"    
    pure function locate(xx,x) 
    ! The code logic and the example of function "locate" below are checked. 1-29-2017
    ! vec = [1,2,3], ans(0.5)=0, ans(1)=1, ans(2)=2, ans(3)="2", ans(3.5)= 3. So when you see in this case 0 or 3, the value falls ouside the range.
   	    implicit none
   	    real(wp), dimension(:), intent(in) :: xx
   	    real(wp), intent(in) :: x
   	    integer :: locate
   	    integer :: n,il,im,iu
   	    n  = size(xx)
   	    il = 0
   	    iu = n+1
   	    do
   	    	if ( iu-il <= 1 ) exit
            im = (iu+il)/2
   	    	if ( x >= xx(im) ) then
   	    		il=im
   	    	else
   	    		iu=im
   	    	endif
   	    enddo
   	    if ( x == xx(1) ) then
   	    	locate=1
   	    elseif ( x == xx(n) ) then
   	    	locate=n-1
   	    else
   	    	locate=il
   	    endif
    end function locate

    subroutine bsearch(func,imin,imax,is,ys,message)
    ! Binary search. Return the maximizer is and the correponding value ys which it attains: ys = func(is)
    ! Input : the strictly concave function func(x), where x is the index of the discrete asset grids
    !         the lower bound and upper bound that form an interval where the binary search works on
    ! Output : the index of the maximizer: is
    !          the corresponding function value : ys
        implicit none
        integer, intent(in) :: imin,imax
        integer, intent(out) :: is
        real(wp), optional, intent(out) :: ys
        character(len=*), intent(in) :: message
        interface
            function func(x)
            use f95_precision, only: wp=>dp
            implicit none
            integer, intent(in) :: x
            real(wp) :: func
            end function func        
        end interface
        integer :: il,im,iu,i3(3),i3s(1)
        real(wp),dimension(3) :: y3
        if(imin>imax) then
            write(6,'(a,a)'), '! imin>imax @ ', message
            ys = 0._wp
            is = 0
        elseif(imin==imax) then
            is=imax
            ys=func(imax)
        else
            il=imin
            iu=imax
            do 
                if(iu-il<=2) exit
                im=(iu+il)/2
                if(func(im+1)>func(im))then
                    il=im
                else
                    iu=im+1
                endif
            enddo
            i3 = [il,il+1,iu]
            y3 = [func(il),func(il+1),func(iu)]
            i3s = maxloc(y3)
            is = i3(i3s(1))
            ys = y3(i3s(1))
        endif   
    end subroutine bsearch
    
    subroutine bsearch_err( s, i, min, max )
        implicit none
        character(len=*), intent(in) :: s
        integer, intent(in) :: i, min, max
        write(6,'(a,x,3(i4))'), s,i,min,max
    end subroutine bsearch_err
    
    subroutine brentmaxrhs(func,x,i,imax,xs,ys,numidx,vec) ! 3.7.2017 don't change it.
    ! input : func: the function where the local optimization works on
    !         x, i, imax : the asset grid vector, the index of the maximizer from the coarse grid, and the index corresponding to nonnegative consumption
    ! output : xs : the next-period asset holdings (note it is a real type variable, not a integer one.)
    !          ys : the associated value of func attained by choosing xs 
    ! optional : numidx, vec
        implicit none
        real(wp), dimension(:), intent(in) :: x
        integer, intent(in) :: i,imax
        real(wp), intent(out) :: xs, ys
        integer, intent(in), optional :: numidx, vec(:)
        logical :: exist ! 
        interface
            function func(x)
            use f95_precision, only: wp=>dp
            implicit none
            real(wp), intent(in) :: x
            real(wp) :: func
            end function func       
        end interface
        real(wp), parameter :: e=0.0001_wp
        real(wp) :: xl,xu,xi,yi,xe,ye
        
        if(i>1)then
            xl=x(i-1) ! set up the lower bound   
        endif
        
        xi=x(i) ! the coarse maxima
        yi=func(xi) 
        if(i<imax)then ! set up the upper bound
            xu=x(i+1)
        endif
        if(i>1 .and. i<imax)then ! concave case. 3.7.2016
            ys=brentmax(xl,xi,xu,func,xs) ! normal cases.     
        else ! likely to be a corner solution case.
            if(i==1)then ! the coarse maxima is located at the lower bound
                xu = x(i+1) ! fix IK's bug 1/2
                xe = xi+(xu-xi)*e
                ye = func(xe)
                if(ye<yi)then ! corner solution confirmed.
                    xs = xi
                    ys = yi
                else ! there exist a tiny curve insdie the interval.
                    ys = brentmax(xi,(xi+xu)/2,xu,func,xs)
                endif
            else
                if(i==imax)then
                    xl = x(i-1) ! fix IK's bug 2/2
                    xe = xi-(xi-xl)*e
                    ye = func(xe)
                    if(ye<yi)then ! corner solution confirmed.
                        xs = xi
                        ys = yi
                    else ! there exist a tiny curve insdie the interval.
                        ys = brentmax(xl,(xl+xi)/2,xi,func,xs)
                    endif
                else ! Li-Pin's extension for returning message. 3.7.2017 I guess this error message is useless. I think there is nothing wrong with IK's code.
                    if (present(numidx) .and. present(vec) ) then
                        inquire(file="brentmax.txt", exist=exist)
                        if (exist) then
                            open(12, file="brentmax.txt", status="old", position="append", action="write")    
                        else
                            open(12, file="brentmax.txt", status="new", action="write")
                        endif
                        write(12, '(a,<numidx>i3)') 'i<1 or i>imax', vec(1:numidx)
                        close(12)
                    else
                        write(*,'(a)') 'upper corner solution found.'
                    endif
                endif
                
            endif
        endif    
    end subroutine brentmaxrhs

    function brentmax(ax,bx,cx,func,xmax) ! 3.7.2017 don't change anything in this function.
        ! look for the minimizor XMIN for the penalty function FUNC.
        ! Note that BX is between AX and CX. 
        ! Return the maximum value the function FUNC attains.
        implicit none
        real(wp), intent(in) :: ax,bx,cx
        real(wp), intent(out) :: xmax
        real(wp) :: brentmax
        
        interface
          function func(x)
          use f95_precision, only: wp => dp
          implicit none
          real(wp), intent(in) :: x
          real(wp) :: func
          end function func
        end interface
        
        integer, parameter :: itmax=100
        real(wp), parameter :: tol=sqrt(epsilon(ax)),cgold=0.381966011250105_wp,zeps=1.0e-3_wp*epsilon(ax)
        integer :: iter
        real(wp) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
        a=min(ax,cx)
        b=max(ax,cx)
        v=bx
        w=v
        x=v
        e=0.0
        fx=-func(x) !fixed caused the original version is for minimum. My program looks for maximum.
        fv=fx
        fw=fx
        do iter=1,itmax
          xm=0.5_wp*(a+b)
          tol1=tol*abs(x)+zeps
          tol2=2.0_wp*tol1
          
          if (abs(x-xm) <= (tol2-0.5_wp*(b-a))) then ! the program found the maximum (in my application)
            xmax=x
            brentmax=-fx ! invert the minimum to be the maximum of my application. inverse the sign of the minimum so that it becomes the maximum we are looking for.
            return
          end if
          
          if (abs(e) > tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.0_wp*(q-r)
            if (q > 0.0) p=-p
            q=abs(q)
            etemp=e
            e=d
            if (abs(p) >= abs(0.5_wp*q*etemp) .or. &
              p <= q*(a-x) .or. p >= q*(b-x)) then
              e=merge(a-x,b-x, x >= xm )
              d=cgold*e
            else
              d=p/q
              u=x+d
              if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
            end if
          else
            e=merge(a-x,b-x, x >= xm )
            d=cgold*e
          end if
          
          u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
          fu=-func(u) ! correct. I adapt my maximization problem into the minimum searching algorithm tentatively. 
          if (fu <= fx) then
            if (u >= x) then
              a=x
            else
              b=x
            end if
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
          else
            if (u < x) then
              a=u
            else
              b=u
            end if
            if (fu <= fw .or. w == x) then
              v=w
              fv=fw
              w=u
              fw=fu
            else if (fu <= fv .or. v == x .or. v == w) then
              v=u
              fv=fu
            end if
          end if
        end do
        ! stop 'brent: exceed maximum iterations'
        write(*,'(a)') 'brent: exceed maximum iterations'
    end function brentmax
    
    ! Function: set up a proper interval that corresponds to a concave function which Brent can only work for. 3.7.2017
    ! Input:
    !   xa: lower bound of financial asset holding
    !   xc: upper bound of financial asset holding
    subroutine brent_localizer(func,xa,xc,xs,ys)
        implicit none
        real(wp), intent(out) :: xs, ys
        real(wp), intent(in) :: xa,xc
        integer :: m, log0(1), size
        real(wp) :: ya, yb, yc, yaa, ycc, xaa, xcc
        real(wp) :: nxb, nxa, nxc
        real(wp), dimension(:), allocatable :: xvec, yvec
        real(wp), parameter :: e=0.001_wp      
        
        interface
            function func(x)
            use f95_precision, only: wp=>dp
            implicit none
            real(wp), intent(in) :: x
            real(wp) :: func
            end function func       
        end interface
        
        size = 5 ! 5 is good FOR NOW; given a level of housing assets, the value function is concave. 09302016 
        ! 5 is good FOR NOW. we just want to have a good inivial guess for refined brent method. 
        nxa = xa 
        nxc = xc
        allocate( xvec(size), yvec(size) )
        
        call grid(xvec,nxa,nxc,1.0_wp)
        do m = 1, size
            yvec(m) = func(xvec(m))                
        enddo
        log0 = maxloc(yvec)
        if(log0(1)==1)then ! a corner solution is "likely" to be at lower bound.
            nxa = xvec(1)
            nxc = xvec(2)
            nxb = (nxa+nxc)/2._wp
        elseif(log0(1)==size)then ! a corner solution is "likely" to be at upper bound.
            nxa = xvec(log0(1)-1) 
            nxc = xvec(log0(1))
            nxb = (nxa+nxc)/2._wp
        else
            nxa = xvec(log0(1)-1) 
            nxc = xvec(log0(1)+1)
            nxb = (nxa+nxc)/2._wp
        endif
        ys = brentmax(nxa,nxb,nxc,func,xs)                              
        deallocate( xvec, yvec )   
        return
    end subroutine brent_localizer

    subroutine splie2(x1a,x2a,ya,y2a)
        ! precompute the auxiliary second-derivative table for bicubic spline
	    implicit none
	    real(wp), dimension(:), intent(in) :: x1a,x2a
	    real(wp), dimension(:,:), intent(in) :: ya
	    real(wp), dimension(:,:), intent(out) :: y2a
	    integer :: j,m,ndum
	    m=size(x1a)
	    ndum=size(x2a)
	    do j=1,m
	    	call spline(x2a,ya(j,:),y2a(j,:))
	    end do
    end subroutine splie2
      
    function splin2(x1a,x2a,ya,y2a,x1,x2)
        ! the main unit for bicubic spline interpolations
	    implicit none
	    real(wp), dimension(:), intent(in) :: x1a,x2a
	    real(wp), dimension(:,:), intent(in) :: ya,y2a
	    real(wp), intent(in) :: x1,x2
	    real(wp) :: splin2
	    integer :: j,m,ndum
	    real(wp), dimension(size(x1a)) :: yytmp,y2tmp2
	    m=size(x1a)
	    ndum=size(x2a)
	    do j=1,m
	    	yytmp(j)=splint(x2a,ya(j,:),y2a(j,:),x2)
	    end do
	    call spline(x1a,yytmp,y2tmp2)
	    splin2=splint(x1a,yytmp,y2tmp2,x1)
    end function splin2
    
    function linterp(pt,val,x)
        implicit none
        real(wp), dimension(:), intent(in) :: pt,val
        real(wp), intent(in) :: x
        integer :: i,n1,n2
        real(wp) :: d,a,b,linterp
        n1=size(pt)
        n2=size(val)
        if( n1/=n2 ) stop 'program terminated by linterp.'
        i=max(min(locate(pt,x),n1-1),1)
        d=pt(i+1)-pt(i)
        if( d==0.0 ) stop 'bad x input in linterp.'
        a=(pt(i+1)-x)/d ! lower part's weight
        b=(x-pt(i))/d ! upper part's weight
        linterp= a*val(i) + b*val(i+1)
    end function linterp
        
    function binterp(x,y,sf,xp,yp)
        implicit none
        real(wp), dimension(:), intent(in) :: x,y ! (x,y) : row and column
        real(wp), intent(in) :: xp,yp ! a point on 2D plane
        real(wp), dimension(:,:), intent(in) :: sf ! surface
        real(wp) :: ax,bx,ay,by,binterp
        integer :: n1,n2,m1,m2,ix,iy
        n1=size(x)
        n2=size(sf,dim=1)
        m1=size(y)
        m2=size(sf,dim=2)
        if( n1/=n2 .or. m1/=m2 ) stop 'program terminated by binterp.'
        ix=max(min(locate(x,xp),n1-1),1)
        iy=max(min(locate(y,yp),m1-1),1)
        ax=( x(ix+1)-xp )/( x(ix+1)-x(ix) ) ! weight for the part close to the northwest corner
        bx=( xp-x(ix) )/( x(ix+1)-x(ix) )     
        ay=( y(iy+1)-yp )/( y(iy+1)-y(iy) ) ! weight for the part close to the northwest corner
        by=( yp-y(iy) )/( y(iy+1)-y(iy) )
        binterp= sf(  ix,  iy)*ax*ay + &  ! NW
                 sf(ix+1,  iy)*bx*ay + &  ! SW
                 sf(  ix,iy+1)*ax*by + &  ! NE
                 sf(ix+1,iy+1)*bx*by      ! SE
    end function binterp
    
    ! Function: bilinear interpolation 3.11.2017
    ! Inputs:
    !   x: first dimension vector (row) Note: It is not the usual x-axis (or horizontal axis) on a two-dimensional cartesian system.
    !   y: second dimension vector (column). Note: It is not the usual y-axis (or vertical axis) on a two-dimensional cartesian system.
    !   sf: thrid dimension matrix corresponding to the combinations of (x,y)
    !   xp: first dimension point to be interpolated
    !   yp: second dimension point to be interpolated
    !   penalty: --
    ! Outputs:
    !   msg: logical type of variable
    !   f4: the values corresponding to the combinations located on the four corners.
    function binterpII(x,y,sf,xp,yp,penalty,msg,f4) ! 3.11.2017 This function is verified that produced exactly the same results as the INTERP2 subroutine in MATALB (LINK: http://www.obs.ujf-grenoble.fr/scci/logiciels/matlab61/help/techdoc/ref/interp2.html)
        implicit none
        real(wp), dimension(:), intent(in) :: x,y ! (x,y) : row and column
        real(wp), intent(in) :: xp,yp ! a point on 2D plane
        real(wp), dimension(:,:), intent(in) :: sf ! surface
        real(wp) :: ax,bx,ay,by,binterpII,fnw,fsw,fne,fse
        real(wp), intent(out) :: f4(4)
        integer :: n1,n2,m1,m2,ix,iy
        real(wp), intent(in) :: penalty
        logical, intent(inout) :: msg
        msg = .false.
        n1=size(x)
        n2=size(sf,dim=1)
        m1=size(y)
        m2=size(sf,dim=2)
        if( n1/=n2 .or. m1/=m2 ) stop 'program terminated by binterp.'
        ix=max(min(locate(x,xp),n1-1),1)
        iy=max(min(locate(y,yp),m1-1),1)
        ax=( x(ix+1)-xp )/( x(ix+1)-x(ix) ) ! weight for the part close to the northwest corner
        bx=( xp-x(ix) )/( x(ix+1)-x(ix) )     
        ay=( y(iy+1)-yp )/( y(iy+1)-y(iy) ) ! weight for the part close to the northwest corner
        by=( yp-y(iy) )/( y(iy+1)-y(iy) )
        fnw = sf(  ix,  iy) 
        fsw = sf(ix+1,  iy)
        fne = sf(  ix,iy+1)
        fse = sf(ix+1,iy+1)
        f4(1) = fnw
        f4(2) = fsw
        f4(3) = fne
        f4(4) = fse
        !if(any(f4==penalty))then
        if(any(abs(f4-penalty)<penalty/1.e5))then ! 4.23.2017 fixed the issue caused by high hmin.
            binterpII = penalty
            msg = .true.
        else            
            binterpII = fnw*ax*ay + &  ! NW
                        fsw*bx*ay + &  ! SW
                        fne*ax*by + &  ! NE
                        fse*bx*by      ! SE
            msg = .false.
        endif 
    end function binterpII    
    
    ! --------------- sort 2 arrays ---------------------
    ! Three components : sort2array, indexx_wp, arth
    
    subroutine sort2arrayII(arr,slave1,opt) ! IN ASCENDING ORDER 
    ! sort arr and in turn rearrange slave1 along the rid 
    ! of sorting the corresponding elements in arr
    ! This subroutine has to work with several 
    ! complementary subroutines: indexx_wp and arth
      implicit none
      real(wp), dimension(:), intent(inout) :: arr,slave1
      integer, dimension(size(arr)) :: index
      character(len=*), intent(in) :: opt
      if(opt=='ascending')then
        call indexx_wp(arr,index)
        arr=arr(index)
        slave1=slave1(index)
      elseif(opt=='descending')then
        arr=-arr
        call indexx_wp(arr,index)
        arr=arr(index)
        arr=-arr
        slave1=slave1(index)
      endif
    end subroutine sort2arrayII
    
    subroutine sort2array(arr,slave1) ! IN ASCENDING ORDER 
    ! sort arr and in turn rearrange slave1 along the rid 
    ! of sorting the corresponding elements in arr
    ! This subroutine has to work with several 
    ! complementary subroutines: indexx_wp and arth
      implicit none
      real(wp), dimension(:), intent(inout) :: arr,slave1
      integer, dimension(size(arr)) :: index
      call indexx_wp(arr,index)
      arr=arr(index)
      slave1=slave1(index)
    end subroutine sort2array    
    
    subroutine indexx_wp(arr,index)
	  implicit none
	  real(wp), dimension(:), intent(in) :: arr
	  integer, dimension(:), intent(out) :: index
	  integer, parameter :: nn=15, nstack=50
	  real(wp) :: a
	  integer :: n,k,i,j,indext,jstack,l,r
	  integer, dimension(nstack) :: istack
	  n=size(index)
	  index=arth(1,1,n)
	  jstack=0
	  l=1
	  r=n
	  do
	  	if (r-l < nn) then
	  		do j=l+1,r
	  			indext=index(j)
	  			a=arr(indext)
	  			do i=j-1,l,-1
	  				if (arr(index(i)) <= a) exit
	  				index(i+1)=index(i)
	  			end do
	  			index(i+1)=indext
	  		end do
	  		if (jstack == 0) return
	  		r=istack(jstack)
	  		l=istack(jstack-1)
	  		jstack=jstack-2
	  	else
	  		k=(l+r)/2
	  		call swap(index(k),index(l+1))
	  		call icomp_xchg(index(l),index(r))
	  		call icomp_xchg(index(l+1),index(r))
	  		call icomp_xchg(index(l),index(l+1))
	  		i=l+1
	  		j=r
	  		indext=index(l+1)
	  		a=arr(indext)
	  		do
	  			do
	  				i=i+1
	  				if (arr(index(i)) >= a) exit
	  			end do
	  			do
	  				j=j-1
	  				if (arr(index(j)) <= a) exit
	  			end do
	  			if (j < i) exit
	  			call swap(index(i),index(j))
	  		end do
	  		index(l+1)=index(j)
	  		index(j)=indext
	  		jstack=jstack+2
	  		if (jstack > nstack) stop 'indexx: nstack too small'
	  		if (r-i+1 >= j-l) then
	  			istack(jstack)=r
	  			istack(jstack-1)=i
	  			r=j-1
	  		else
	  			istack(jstack)=j-1
	  			istack(jstack-1)=l
	  			l=i
	  		end if
	  	end if
	  end do
    contains
	  subroutine icomp_xchg(i,j)
	    integer, intent(inout) :: i,j
	    integer :: swp
	    if (arr(j) < arr(i)) then
	    	swp=i
	    	i=j
	    	j=swp
	    end if
	  end subroutine icomp_xchg
	  
	  subroutine swap(a,b)
	    integer, intent(inout) :: a,b
	    integer :: dum
	    dum=a
	    a=b
	    b=dum
	  end subroutine swap
    end subroutine indexx_wp
    
    function arth(first,increment,n)
      integer, intent(in) :: first,increment
      integer, intent(in) :: n
      integer, dimension(n) :: arth
      integer :: k,k2
      real(wp) :: temp
      if (n > 0) arth(1)=first
      if (n <= 16) then
      	do k=2,n
      		arth(k)=arth(k-1)+increment
      	end do
      else
      	do k=2,8
      		arth(k)=arth(k-1)+increment
      	end do
      	temp=increment*8
      	k=8
      	do
      		if (k >= n) exit
      		k2=k+k
      		arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
      		temp=temp+temp
      		k=k2
      	end do
      end if
    end function arth    
    
    ! generate lorenz curve and the gini index
    subroutine lorenz(f,x,fx,gini)
    ! compute lorenz curve and gini coefficient by sorting vector x and using density f
    ! note: to link dlasrt2 add mkl_scalapack_core.lib to project's additional dependencies
    ! Verified by Stata Ineqdeco module. Accuracy to at least thousandth place. 10082016.
    	implicit none
    	real(wp), dimension(:), intent(in) :: f
    	real(wp), dimension(:), intent(inout) :: x
    	real(wp), dimension(:), intent(out) :: fx
    	real(wp), intent(out) :: gini
    	integer, dimension(size(x)) :: key
    	integer :: n,i,info
    	n=size(x)
    	if (size(f)/=n) then
    		print '(a,i3)', 'lorenz: x and f must be of the same size ',n
    	!	stop 'program terminated by lorenz'
    	end if
    	if (size(fx)/=n) then
    		print '(a,i3)', 'lorenz: x and fx must be of the same size ',n
    	!	stop 'program terminated by lorenz'
    	end if
    	key=(/ (i,i=1,n) /)
    	call dlasrt2('i',n,x,key,info)
    	fx=f(key)
    	x=x*fx
    	gini=x(1)*fx(1)
    	do i=2,n
    		x(i)=x(i)+x(i-1)
    		gini=gini+(x(i)+x(i-1))*fx(i)
    		fx(i)=fx(i)+fx(i-1)
    	end do
    	gini=1-gini/x(n)
    	x=x/x(n)
    end subroutine lorenz  

    function maxlocation(vec)
        implicit none
        real(wp), dimension(:), intent(in) :: vec
        integer :: maxlocation, ans(1)
        ans = maxloc(vec)
        maxlocation = ans(1)        
    end function maxlocation
    
    function minlocation(vec)
        implicit none
        real(wp), dimension(:), intent(in) :: vec
        integer :: minlocation, ans(1)
        ans = minloc(vec)
        minlocation = ans(1)
    end function minlocation
    
    ! 4.21.2017 checked. Perfect with the python percentile_wiki.py
    ! INPUTS, SERIES: VECTOR OF (UNSORTED) LEVELS TO BE SORTED BEFORE OBTAINING THE PERCENTILE FROM THEM. 
    !         WEIGHTS: VECTOR OF WEIGHTS CORRESPONDING INPUT VECTOR SERIES (NORMALIZATION IS NOT REQUIRED)
    !         QPROP (%): THE PERCNETILE WHOSE VALUE WE WOULD LIKE TO KNOW. MUST BE 0<=QPROP<=1.
    ! OUTPUTS, QUANT: THE VALUE BELOW WHICH THE GIVEN FFRACTION OF OBSERVATIONS IN THE GIVEN OBSERVATIONS FALL.
    !          EMSG (OPTIONAL) : IF THE TOTAL WEIGHT SUMS TO ZERO, EMSG IS 1 INDICATING ERROR HAPPENS; OTHERWISE EMSG = 0.
    ! 5.8.2015 COMPARED WITH CAGETTI AND DENARDI'S VERSION, THEY HAVE TWO BUGS IN THEIR OWN PERCENTILE SUBROUTINE 
    ! EXAMPLE: call weighted_percentile(series,weights,0.02,answer,error_message)
    ! IF YOU ARE NOT SURE WHETHER THE WEIGHT VECTOR IS NONNEGATIVE, BE SURE TO USE THE OPTIMAL VARIABLE EMSG
    ! REFERENCE: http://en.wikipedia.org/wiki/Percentile The WEIGHTED PERCENTILE METHOD
    subroutine weighted_percentile( series,weights,qprop,quant,emsg )
        implicit none
        real(wp), dimension(:), intent(in) :: series, weights
        real(wp), intent(in) :: qprop
        real(wp), intent(out) :: quant
        real(wp), dimension(:), allocatable :: ssorded, wsorded, csum, pvec, rkey
        real(wp) :: prop
        integer, dimension(:), allocatable :: key
        integer, dimension(1) :: maxi
        integer :: m, n, i
        integer, optional :: emsg
        if( present(emsg) ) emsg = 0
        m = size(series)
        prop = 100._wp*qprop
        allocate( key(m) )
        key = 1
        n = sum(key,weights>1.e-10_wp) ! GET RID OF ZERO WEIGHTS
        deallocate( key )
        !write(*,fmt='(a,i8)') "the length is: ", n
        if (n>0) then
            allocate( ssorded(n), wsorded(n), rkey(n), key(n), csum(n), pvec(n) )
            ssorded = pack(series,weights>1.e-10_wp)
            wsorded = pack(weights,weights>1.e-10_wp)
            key = [(i, i=1,n)]
            rkey = real(key,wp)
            call sort2array(ssorded,rkey)
            key = int(rkey)
            wsorded = wsorded(key)
            csum(1) = wsorded(1)
            do i = 2, n 
                csum(i) = csum(i-1)+wsorded(i)
            enddo
            do i = 1, n ! COMPUTE THE "PERCENT RANK" OF THAT VALUE IN THE ORDERED LIST OF SIZE N
                pvec(i) = 100._wp/csum(n)*(csum(i)-wsorded(i)/2._wp)   
                write(unit=133,fmt='(i5,2x,4(x,f20.15))')  i, pvec(i), csum(i), wsorded(i), ssorded(i)
            enddo
            ! THEN WE TAKE THE OBTAINED PERCENT RANK AND LOCATE THE PERCENTILE VALUE WE ARE INTERESTED
            maxi = maxloc(pvec,pvec<=prop)
            m = maxi(1)
            !write(*,fmt='(a,i8)') "the max position is: ", m
            if ( prop<=pvec(1) ) then
                quant = ssorded(1)    
            elseif ( prop>pvec(n) ) then
                quant = ssorded(n)
            else
                if(abs((pvec(m+1)-pvec(m)))<1.e-11)then
                    write(*,fmt='(a,2(f20.15))'), "##### ZERO ######", pvec(m+1), pvec(m)
                    write(*,fmt='(2f25.12)') sum(ssorded), sum(wsorded)
                    write(*,fmt='(4f25.12)') csum(m), csum(m+1), wsorded(m), wsorded(m+1) 
                endif
                quant = ssorded(m) + (prop-pvec(m))/(pvec(m+1)-pvec(m))*(ssorded(m+1)-ssorded(m))
            endif
        else
            if( present(emsg) ) emsg=1
        endif
    end subroutine weighted_percentile    
    
    ! INPUT:
    !    c : centroid
    !    r : parameter range
    !    b : best
    !  fac : the percentage points added to best
    !  nth : round the number to nth decimal place
    ! OUTPUT:
    !    cp: revised centroid
    subroutine avoid_stalemate(c,b,r,dev,fac,nth)
        implicit none
        real(wp), dimension(:), intent(in) :: c,b
        real(wp), dimension(:,:), intent(in) :: r
        real(wp), dimension(:), intent(out) :: dev
        real(wp), intent(in) :: fac
        integer, intent(in) :: nth
        real(wp), dimension(:), allocatable :: cp, bp, incr
        integer :: n
        n = size(c)
        allocate( cp(n), bp(n), incr(n) )
        cp = c
        bp = b
        incr = abs(r(:,2)-r(:,1))
        cp = nint(cp*10._wp**nth)/(10._wp**nth)
        bp = nint(bp*10._wp**nth)/(10._wp**nth)
        ! dev = merge( fac*cp, 0._wp, cp == bp )
        dev = merge( fac*incr, 0._wp, cp == bp )
        deallocate( cp, bp, incr ) 
    end subroutine avoid_stalemate
    
    subroutine avoid_stalemate2(cmat,b,dev,fac,nth)
        implicit none
        real(wp), dimension(:,:), intent(in) :: cmat
        real(wp), dimension(:), intent(in) :: b
        real(wp), dimension(:), intent(out) :: dev
        real(wp), intent(in) :: fac
        integer, intent(in) :: nth
        real(wp), dimension(:), allocatable :: bp
        real(wp), dimension(:,:), allocatable :: cmatp
        integer :: n,m
        m = size(cmat,1)
        n = size(cmat,2)
        allocate( cmatp(m,n), bp(n))
        cmatp = cmat
        bp = b
        !cmatp = nint(cmatp*10._wp**nth)/(10._wp**nth)
        !bp = nint(bp*10._wp**nth)/(10._wp**nth)
        cmatp = cmatp - spread(bp,1,m)
        dev = merge( fac*bp, 0._wp, sum(cmatp,1)==0._wp )
        deallocate( cmatp, bp ) 
    end subroutine avoid_stalemate2     
    
    ! INPUT: 
    !   ray: the serires of which we would like to know the average.
    !   idx: the latest data point in the series.
    !   len: the size of the window.
    function moving_average(ray,idx,len)    
        implicit none
        real(wp), dimension(:), intent(in) :: ray
        integer, intent(in) :: idx, len
        real(wp) :: moving_average
        integer :: n
        if( idx < len ) then
            moving_average = 9999._wp
        else
            moving_average = sum(ray(idx-len+1:idx))/len    
        endif
    end function moving_average
    
    function avg(vec)
        implicit none
        real(wp), dimension(:) :: vec
        real(wp) :: avg
        integer :: n
        n = size(vec)
        if (n==0) write(*,'("Denominator is zero when taking average.")')
        avg = sum(vec)/real(n,wp)
    end function avg
    
    subroutine difference_series(ray,diff)
        implicit none
        real(wp), dimension(:), intent(in) :: ray
        real(wp), dimension(:), intent(out) :: diff
        integer :: n
        n = size(ray)
        diff = 0._wp
        diff(2:n) = ray(1:n-1)
        diff = ray - diff
    end subroutine difference_series
    
    ! INPUT:
    !    i : represents the index of the latest data point in the ray.
    ! FUNCTION:
    !   return the absolute difference between the current step (i) and the last step (i-1).
    subroutine percentage_absolute_difference_series(ray,diff,i)
        implicit none
        real(wp), dimension(:), intent(in) :: ray
        real(wp), dimension(:), intent(inout) :: diff
        integer :: n,i
        diff = sqrt(-1._wp)
        if( i>=2 ) then
            diff(2:i) = ray(1:i-1)
            diff(2:i) = abs(ray(2:i) - diff(2:i))/diff(2:i)
        endif
        diff(1:i-1) = diff(2:i)
    end subroutine percentage_absolute_difference_series  
    
    subroutine infinity_setting(var)
        implicit none
        real(wp), intent(inout) :: var
        if( ieee_support_inf(var) )then
            var = -ieee_value(var, ieee_negative_inf)
        endif        
    end subroutine infinity_setting
    
    ! See http://math.nju.edu.cn/help/mathhpc/doc/intel/mkl/vslnotes.pdf page 32 for generating 100 3-dimensional quasi-random vectors in the (2,3)^3 hypercube using Sobol BRNG.
    ! Be sure to add the command "include 'mkl_vsl.f90'" in the beginning of module toolbox(.f90)
    ! Reference: 
    ! See how to set up a module for mkl_vsl https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/515066
    ! See "brng" parameter selection https://software.intel.com/en-us/mkl-developer-reference-c-brng-parameter-definition#TBL10-2
    ! Check the error code which should be 0 if it worked. https://goparallel.sourceforge.net/mastering-intel-mkl-support-functions/ 
    ! How to include "errcheck.inc" https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/292173
    subroutine get_sobol_sequence( sobolm, lbnd, ubnd, nsbq, ndim )
        implicit none
        integer, intent(in) :: nsbq, ndim
        integer :: brng, method, errcode
        real(wp), intent(out) :: sobolm(nsbq, ndim)
        real(wp) :: lbnd, ubnd
        type(vsl_stream_state) :: stream 
        real(wp) :: sobvec(ndim)
        integer :: i
        brng    = VSL_BRNG_SOBOL             ! The sobol sequence generator's parameters.
        method  = VSL_RNG_METHOD_UNIFORM_STD ! The sobol sequence generator's parameters.    
        errcode = vslnewstream( stream, brng, ndim )
        !if(errcode/=0) write(*,*) errcode, ' something wrong with get_sobol_sequence [1] '
        !errcode = vdrnguniform( method, stream, nsbq*ndim, sobolm, lbnd, ubnd )
        !if(errcode/=0) write(*,*) errcode, ' something wrong with get_sobol_sequence [2] '
        do i = 1, nsbq
            errcode = vdrnguniform(method, stream, ndim, sobvec, lbnd, ubnd)    
            if(errcode/=0) write(*,*) errcode, ' something wrong with get_sobol_sequence [1] '
            !write(*,*) sobvec
            !write(*,'(/,a)') ' --- ' 
            sobolm(i,:) = sobvec(:)
        enddo        
    end subroutine get_sobol_sequence
    
    real(wp) function linint(x,y,xi) ! LINEAR INTERPOLATION, SINGLE POINT
        implicit none
        real(wp), dimension(:), intent(in) :: x, y
        real(wp), intent(in) :: xi
        real(wp):: a, b, d
        integer :: n, i
        n = size(x)
        if (size(y)/=n) then
            print*, '[Error] the length of x and y are not the same.'
            stop 'prorram terminated by linint'
        end if
	    i=max(min(locate(x,xi),n-1),1) ! ALWAYS BEING ASSIGNED TO THE LOWER BOUND OF THE BRACKET
	    d=x(i+1)-x(i)
	    if (d == 0.0_wp) stop 'bad x input in linint'
	    a = (x(i+1)-xi)/d
	    b = (xi-x(i))/d
	    linint = a*y(i)+b*y(i+1)        
    end function linint
    
    ! 3.17.2017 linear interpolation.
    subroutine linwgt(grid,pt,xv,wv) ! in: real, grid, pt. out: wv: real weight, xv: integer location
        implicit none
        real(wp), dimension(:), intent(in) :: grid
        real(wp), intent(in) :: pt
        real(wp), dimension(:), intent(out) :: wv
        integer, dimension(:), intent(out) :: xv
        integer :: n
        n = size(grid)
        xv(1) = locate(grid,pt)
        if(xv(1) == 0)then
            xv(1) = 1
            xv(2) = xv(1)+1
            wv(1) = 1
            wv(2) = 0
        elseif(xv(1) == n)then
            xv(1) = n-1
            xv(2) = xv(1)+1
            wv(1) = 0
            wv(2) = 1
        else
            xv(2) = xv(1)+1
            wv(1) = (grid(xv(2))-pt)/(grid(xv(2))-grid(xv(1))) 
            wv(2) = (pt-grid(xv(1)))/(grid(xv(2))-grid(xv(1)))            
        endif
        !xv(2) = merge(n,xv(1)+1,n==xv(1)) ! adjust for outside the boundary
        !xv(1) = merge(n-1,xv(1),n==xv(2)) ! adjust for outside the boundary
        !xv(1) = merge(1,xv(1),1==xv(2))
        !xv(2) = merge(2,xv(1)+1,xv(1)==1)

    end subroutine linwgt
    
    subroutine linintv(x,y,xi,yi) ! LINEAR INTERPOLATION, VECTOR
        implicit none
        real(wp), dimension(:), intent(in) :: x, y, xi
        real(wp), dimension(:), intent(out):: yi
        real(wp):: a, b, d
        integer :: m, n, i, j
        n = size(x)
        if ( size(y)/=n ) then
            print*, '[Error] the length of x and y are not the same.'
            stop 'prorram terminated by linintv'
        end if        
        m = size(xi)
        if ( size(yi)/=m ) then
            print*, '[Error] the length of xi and yi are not the same.'
            stop 'prorram terminated by linintv'
        end if 
        do j = 1, m
		    i = max( min( locate(x,xi(j)),n-1 ), 1)
		    d = x(i+1)-x(i)
		    if (d == 0.0_wp ) then
		    	stop 'bad x input in linintv'
		    endif
		    a = ( x(i+1)-xi(j) )/d
		    b = ( xi(j)-x(i) )/d
		    yi(j) = a*y(i) +b*y(i+1)            
        enddo
    end subroutine linintv     
    
    subroutine cumsum(vec,cumsumvec) 
    ! CUMULATIVE SUM
    ! INPUT: 
    !        VEC IS AN ONE-DIMENSIONAL ARRAY
        implicit none
        real(wp), intent(in), dimension(:) :: vec
        real(wp), intent(out), dimension(:), allocatable :: cumsumvec
        integer :: m, i
        m = size(vec)
        allocate( cumsumvec(m) )
        cumsumvec = 0._wp
        cumsumvec(1) = vec(1)
        do i = 2, m
            cumsumvec(i) = cumsumvec(i-1) + vec(i)    
        enddo
    end subroutine cumsum
    
    subroutine sort_real(ray,key)
    ! SORTING A REAL ARRAY
    ! INPUT: RAY: UNSORTED REAL ARRAY
    ! OUTPUT: RAY: SORTED ARRAY; KEY: CORRESPONDING INTEGER ORDER 
        implicit none
        real(wp), dimension(:), intent(inout) :: ray
        integer, dimension(:), intent(out) :: key
        real(wp), dimension(:), allocatable :: rkey
        integer :: m, i
        m = size(ray)
        allocate( rkey(m) )
        key = [(i,i=1,m)]
        rkey = real(key,wp)
        call sort2array(ray,rkey)
        key = int(rkey)
    end subroutine sort_real
    
    function reverse_vec(vec)
        implicit none
        real(wp), dimension(:), intent(in) :: vec
        real(wp), dimension(:), allocatable :: reverse_vec
        integer :: m, i
        m = size(vec)
        allocate( reverse_vec(m) )
        do i = 1, m
            reverse_vec(i) = vec(m-i+1)    
        enddo
    end function reverse_vec
    
    ! Input must be double-precision floating point numbers. It should be 4.0_wp rather than only 4.0, for example.
    subroutine kron(am,bm,cm,pos) ! verified by comparison with Matlab. 3-16-2015
        real(wp), intent(in) :: am(:,:), bm(:,:)
        real(wp), intent(out) :: cm(:,:)
        integer :: i, j, ar, ac, br, bc, cr, cc
        real(wp), allocatable :: a1(:,:), ah(:,:), bh(:,:), temp(:)
        integer, allocatable :: idxv(:)
        character(len=*), optional :: pos
        character(len=200) :: message
        ar = size(am,dim=1)
        ac = size(am,dim=2)
        br = size(bm,dim=1)
        bc = size(bm,dim=2)
        cr = size(cm,dim=1)
        cc = size(cm,dim=2)
        if (present(pos)) then
            message = 'error in kron product for '//pos
            if( ar*br/=cr .or. ac*bc/=cc ) then
                print*, message
                stop 
            endif
        endif
        allocate( a1(ar*br,1), ah(cr,cc), bh(cr,cc), temp(ar*br), idxv(ac) )
        do j = 1, ac
            a1 = reshape( spread(am(:,j),1,br), (/ar*br,1/) )
            temp = a1(:,1)
            ah(:,1+(j-1)*bc:j*bc) = spread(temp,2,bc)
        enddo
        do j = 1, bc
            a1 = reshape( spread(bm(:,j),2,ar), (/ar*br,1/) )
            temp = a1(:,1)
            idxv = [(j+(i-1)*bc,i=1,ac)]
            bh(:,idxv) = spread( temp,2,ac )
        enddo    
        cm = ah*bh
    end subroutine kron
    
    subroutine scale_sobol_original( oldone, rng, newone )
        implicit none
        real(wp), dimension(:,:), intent(in) :: oldone, rng
        real(wp), dimension(:,:), intent(inout) :: newone
        integer :: i, npar
        real(wp) :: length, lowerend
        npar = size(oldone,dim=1)
        do i = 1, npar
            newone(i,:) = oldone(i,:)*(rng(i,2)-rng(i,1)) + rng(i,1)    
        enddo
    end subroutine scale_sobol_original        
    
    subroutine find_valid_indices(ary,penalty,idx_s,idx_e) ! 072516  
        implicit none
        real(wp), dimension(:), intent(in) :: ary
        real(wp), intent(in) :: penalty
        integer, intent(out) :: idx_s, idx_e
        integer :: n, i, idx(1)
        integer, dimension(:), allocatable :: temp, incstemp
        n = size(ary)
        allocate( temp(n), incstemp(n) )
        temp = -n
        where(ary/=penalty) temp = 0
        incstemp = [(i,i=1,n)]
        temp = temp + incstemp
        idx = minloc(temp,temp>0)
        idx_s = idx(1)
        idx = maxloc(temp,temp>0)
        idx_e = idx(1)
        deallocate( temp,incstemp)
    end subroutine find_valid_indices
    
    function zbrent(func,x1,x2,tol)
        implicit none
        real(wp), intent(in) :: x1,x2,tol
        real(wp) :: zbrent
        interface
          function func(x)
          use f95_precision, only: wp => dp ! <---- needed
          implicit none
          real(wp), intent(in) :: x
          real(wp) :: func
          end function func
        end interface
        integer, parameter :: itmax=100
        real(wp), parameter :: eps=epsilon(x1)
        integer :: iter
        real(wp) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
        a=x1
        b=x2
        fa=func(a)
        fb=func(b)
        if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
          stop 'root must be bracketed for zbrent'
        c=b
        fc=fb
        do iter=1,itmax
          if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
            c=a
            fc=fa
            d=b-a
            e=d
          end if
          if (abs(fc) < abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
          end if
          tol1=2.0_wp*eps*abs(b)+0.5_wp*tol
          xm=0.5_wp*(c-b)
          if (abs(xm) <= tol1 .or. fb == 0.0) then
            zbrent=b
            return
          end if
          if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
            s=fb/fa
            if (a == c) then
              p=2.0_wp*xm*s
              q=1.0_wp-s
            else
              q=fa/fc
              r=fb/fc
              p=s*(2.0_wp*xm*q*(q-r)-(b-a)*(r-1.0_wp))
              q=(q-1.0_wp)*(r-1.0_wp)*(s-1.0_wp)
            end if
            if (p > 0.0) q=-q
            p=abs(p)
            if (2.0_wp*p  <  min(3.0_wp*xm*q-abs(tol1*q),abs(e*q))) then
              e=d
              d=p/q
            else
              d=xm
              e=d
            end if
          else
            d=xm
            e=d
          end if
          a=b
          fa=fb
          b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
          fb=func(b)
        end do
        stop 'zbrent: exceeded maximum iterations'
        zbrent=b
    end function zbrent    
    
end module toolbox