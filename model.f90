module model 
    use variable ! test
    implicit none
    real(wp) :: inc, pzv(2) 
    real(wp), dimension(:,:,:), allocatable :: v2d
    !real(wp), dimension(:,:,:,:), allocatable :: v2dypx   
    integer :: ax, hx, kx, zx, yx, kpx, ypx, opx, apx, hpx, zpx ! rearrange later.
    real(wp) :: a, h, k, z, y, kp, yp, dk, hp, ap
    !$omp threadprivate(ax,hx,kx,zx,yx,kpx,ypx,opx,apx,hpx,zpx,a,h,k,z,y,kp,yp,dk,hp,inc,pzv,v2d,ap,net_worth)
    ! 3.11.2017 The area above is checked.
contains
    real(wp) elemental function ymin(ki) result(y) ! Ymin(k). The minimum business return before paying back borrowing interest. 2016-6-23
        integer, intent(in) :: ki
        real(wp) :: d
        if(ki==0)then ! no busniess idea
            y = 0._wp
        elseif(ki<=kdim-1)then ! if hit by bad luck, the entrepreneur hires no employer and just let the capital depreciates.
            d = delzl(ki)
            y = (1._wp-d)*kv(ki) ! see Meh's eqn 5.
        endif
    end function ymin
    
    real(wp) elemental function ftaxwok(x) result(y) ! a0: btax_, a1: ptax_, a2: stax_
        real(wp), intent(in) :: x
        y = (btaxw-btaxw*(staxw*x**ptaxw+1._wp)**(-1._wp/ptaxw))*x ! cagetti worker
        !y = (btax-btax*(stax*x**ptax+1._wp)**(-1._wp/ptax))*x ! kitao ordinary
        !y = btax*(x*(1._wp-tauss/2._wp)-((x*(1._wp-tauss/2._wp))**(-ptax)+stax)**(-1._wp/ptax))
    end function ftaxwok
    
    real(wp) elemental function ftaxent(x) result(y)
        real(wp), intent(in) :: x
        y = (btaxe-btaxe*(staxe*x**ptaxe+1._wp)**(-1._wp/ptaxe))*x ! cagetti entrepreneur
        !y = (btax-btax*(stax*x**ptax+1._wp)**(-1._wp/ptax))*x ! kitao ordinary
        !y = btax*(x*(1._wp-tauss/2._wp)-((x*(1._wp-tauss/2._wp))**(-ptax)+stax)**(-1._wp/ptax))
    end function ftaxent 
    
    real(wp) elemental function ttaxwok(x) result(y)
        real(wp), intent(in) :: x
        if(x>0._wp)then
            select case(mode6taskid)
                case(0)
                    y = length*ftaxwok(x/length) + taubal*x
                case(1)
                    y = length*ftaxwok(x/length)
                case(2)
                    y = length*ftaxwok(x/length) + taubal*x
            end select
        else
            y = 0._wp
        endif
    end function ttaxwok
    
    real(wp) elemental function ttaxent(x) result(y)
        real(wp), intent(in) :: x
        if(x>0._wp)then
            select case(mode6taskid)
                case(0)                  
                    y = length*ftaxent(x/length) + taubal*x
                case(1)
                    y = length*ftaxent(x/length)
                case(2)
                    y = length*ftaxent(x/length) + taubal*x                    
            end select
        else
            y = 0._wp
        endif
    end function ttaxent
    
    real(wp) elemental function HomTrnsCst(h,hp) result(y)
        implicit none
        real(wp), intent(in) :: h, hp
        ! y = zeta*(h-hp)**2._wp/2._wp
        y = merge(rho1*h+rho2*hp,0._wp,hp<=(1._wp-mu1)*h.or.hp>=(1._wp+mu2)*h)
    end function HomTrnsCst
    
    real(wp) elemental function u(x,h) result(y)
        implicit none
        real(wp), intent(in) :: x, h
        real(wp) :: g
        g = x**(1._wp-theta)*(h+1.e-7_wp)**(theta) ! omega: weights of non-housing in utility function !10.25.2017
        y = (g**(1._wp-sigma)-1._wp)/(1._wp-sigma) ! sigma: risk aversion coefficient
    end function u    
    
    real(wp) elemental function bu(x) result(y)
        implicit none
        real(wp), intent(in) :: x
        y = x**(1._wp-sigma)/(1._wp-sigma)
    end function bu    
    
    subroutine solve_bellman_1014(exit_bellman) ! 072516 072816 8-13-2017 [0]
        implicit none
        logical, intent(inout) :: exit_bellman ! inivialized in calling function 8-13-2017 [1]
        logical :: exit_brent                  ! 8-13-2017 [1] -------------------------------
        logical, dimension(:), allocatable :: qray ! 8-13-2017 [2]
        
        ! (1)
        integer :: tstart, tend, trate, tmax
        integer, dimension(:,:), allocatable :: xlist, blist1014, bd ! index list, and boundary
        integer  :: i, j, n, m, l, q 
        integer  :: loc1(1)
        real(wp) :: intfund, capgain, labsup, labdem, bgp ! "bgp" stands for business gross profit.
        integer  :: opxi ! 09232016
        ! (2)
        real(wp) :: svcmp 
        integer  :: stp1, dsdim ! Bug 09242016
        logical  :: log1 
        real(wp), dimension(:), allocatable :: shv, lowap, highap, smaxv, apvec
        real(wp) :: best, bestl, bestu, bdis, deb1(3)
        !logical :: cond1
        ! (3)
        character(len=2) :: str1
        write(str1,fmt='(i2)') hdim

        call system_clock(tstart,trate,tmax)
        call system_clock(tstart)

        ! Program starts ---------------------------------------------------------------------------------------
        
        allocate( xlist(adim*1*7,8), bd(1,2) )   
        bd = bd14
        xlist = xl14 ! coarse combination with entrepreneurs choosing not to run a business this period
        t = 14
        
        allocate( qray(bd(1,2)) ) ! 8-13-2017 [3]
        qray = .false.            ! 8-13-2017 
        
        !$omp parallel default(shared) private(n,capgain,intfund,stp1,log1,svcmp,shv,lowap,highap,smaxv,apvec,best,bdis,bestl,bestu,dsdim,m,loc1,labsup,labdem,bgp,exit_brent) ! [2] ! 8-13-2017 [4]
        !$omp do
        do q = bd(1,1), bd(1,2) ! income in period 14 could be used later. Proceed all the way toward the end. Nonstoping process.
            do n = 1, hdim ! experiment on each level of housing asset
                
                exit_brent = .false. ! 8-13-2017 [5]
                
                ax  = xlist(q,1)
                hx  = n ! column-major search over asset dimension for each level of housing assest holding
                kx  = xlist(q,3)
                zx  = xlist(q,4)
                yx  = xlist(q,5)
                kpx = xlist(q,6)
                ypx = xlist(q,7)
                opx = xlist(q,8)     
                
                a  = av(ax)
                h  = hv(hx)
                k  = kv(kx)
                z  = merge(0._wp, merge(zlow,z2(kx),zx==0), kx==0)
                y  = yv(yx)
                kp = kv(kpx)
                yp = yv(ypx)
                dk = merge(0._wp, merge(delzl(kx),delzh(kx),zx==0), kx==0)
                
                ! (1) Income NOTE: FOLLOWING KITAO 2008 (PG. 18), WE SEPARATE CAPITAL INCOME FROM THE REST OF INCOME.  
                if(kx==0)then
                    
                    if(printout25)then 
                        !inc = benefit + merge(transbeq,0._wp,printout17) ! 10.26.2017  
                        if(mode6taskid==0)then
                            inc = merge(benefit+merge(transbeq,0._wp,printout17), &
                                        benefit+merge(transbeq,0._wp,printout17)+merge(rd*a,0._wp,a>0._wp), &
                                        tausv>0._wp) ! taxable 3-11-2018
                            inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a, (1._wp+rd)*a, a>=0._wp), &
                                              merge(a,(1._wp+rd)*a,a>0._wp), &
                                              tausv>0._wp) & ! change in net worth 3-11-2018
                                  + (1._wp - deltah) * h - ttaxwok(inc) 
                        elseif(mode6taskid==1)then ! tausv should be set to zero in this case.
                            inc = benefit + merge(transbeq,0._wp,printout17) + (1._wp + rd) * a + (1._wp - deltah) * h 
                            inc = (1._wp - taubal) * inc
                        elseif(mode6taskid==2)then ! 3-11-2018 in this case, capital income tax is always zero.
                            net_worth = benefit + merge(transbeq,0._wp,printout17) + (1._wp + rd) * a + (1._wp - deltah) * h
                            inc = benefit+merge(transbeq,0._wp,printout17)+merge(rd*a,0._wp,a>0._wp) 
                            if(net_worth > 5._wp*exempbar)then
                                inc = inc + merge(a,(1._wp+rd)*a,a>0._wp) + (1._wp-deltah)*h - ttaxwok(inc) - 0.02_wp*(net_worth-5._wp*exempbar) - 0.01_wp*(4._wp*exempbar)   
                            elseif(net_worth>exempbar .and. net_worth<=5._wp*exempbar)then
                                inc = inc + merge(a,(1._wp+rd)*a,a>0._wp) + (1._wp-deltah)*h - ttaxwok(inc) - 0.01_wp*(net_worth - exempbar)   
                            else
                                inc = inc + merge(a,(1._wp+rd)*a,a>0._wp) + (1._wp-deltah)*h - ttaxwok(inc)
                            endif
                        else 
                            print*, "error in taxation 1014-1"
                        endif !mode6taskid
                    else !obsolete and wrong. 11.3.2017
                        inc = benefit ! 3.26.2017
                        inc = inc + merge( (1._wp+rd)*a, (1._wp+(1._wp-tausv)*rd)*a, a<0._wp) &
                              + (1._wp-deltah)*h - ttaxwok(inc) + merge(transbeq,0._wp,printout17) ! 3.25.2017
                    endif ! printout25 10.26.2017
                else ! 3.29.2017 This case can take care of the negative business shock.
                    labsup = 0._wp ! [3] ! 3.26.2017 useless in the current lifecycle stage.
                    intfund = merge(a - k, 0._wp, a > k) !! savings: sources of interest income ! 3.26.2017                             
                    bgp = c_grs_mat(ax, kx, zx, yx, t)
                    if(printout25)then
                        !inc = inc + bgp  + merge(transbeq,0._wp,printout17)    
                        if(mode6taskid==0)then
                            !inc = benefit + bgp  + merge(transbeq,0._wp,printout17)    
                            inc = bgp + merge(benefit+merge(transbeq,0._wp,printout17), &
                                              benefit+merge(transbeq,0._wp,printout17)+merge(rd*intfund,0._wp,intfund>0._wp), &
                                              tausv>0._wp)
                            inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*intfund, 0._wp, intfund>=0._wp), &
                                              merge(intfund,0._wp,intfund>0._wp), &
                                              tausv>0._wp) & ! change in net worth 3-11-2018
                                  + (1._wp - deltah) * h - ttaxent(inc)                             
                        elseif(mode6taskid==1)then
                            inc = bgp + benefit + merge(transbeq, 0._wp, printout17)                               
                            inc = inc + merge((1._wp + rd) * intfund, 0._wp, intfund>0._wp) + (1._wp-deltah) * h
                            inc = (1._wp - taubal) * inc
                        elseif(mode6taskid==2)then ! 3-1-2018 needs to be revised heavily. ! 2-25-2018 business asset is included in taxable net worth.
                            net_worth = bgp + benefit + merge(transbeq,0._wp,printout17) + merge((1._wp + rd) * intfund, 0._wp, intfund>0._wp) + (1._wp - deltah) * h
                            inc = bgp + benefit + merge(transbeq,0._wp,printout17) + merge(rd*intfund,0._wp,intfund>0._wp) 
                            if(net_worth > 5._wp*exempbar)then
                                inc = inc + merge(intfund,0._wp,intfund>0._wp) + (1._wp-deltah)*h - ttaxent(inc) - 0.02_wp*(net_worth-5._wp*exempbar) - 0.01_wp*(4._wp*exempbar)   
                            elseif(net_worth>exempbar .and. net_worth<=5._wp*exempbar)then
                                inc = inc + merge(intfund,0._wp,intfund>0._wp) + (1._wp-deltah)*h - ttaxent(inc) - 0.01_wp*(net_worth - exempbar)   
                            else
                                inc = inc + merge(intfund,0._wp,intfund>0._wp) + (1._wp-deltah)*h - ttaxent(inc)
                            endif                            
                            
                            !inc = benefit + bgp + merge(transbeq, 0._wp, printout17) ! 3-1-2018 may be wrong.
                            !net_worth = (1._wp - deltah) * h + merge((1._wp + rd) * intfund, merge(0._wp, (1._wp + rd) * a, a > 0._wp), intfund > 0._wp)
                            !if(net_worth > exempbar)then
                            !    inc = inc + (1._wp - taubal) * net_worth - ttaxent(inc)
                            !else
                            !    net_worth = (1._wp - deltah) * h + merge((1._wp + (1._wp - tausv) * rd) * intfund, merge(0._wp, (1._wp + rd) * a, a > 0._wp), intfund > 0._wp)
                            !    inc = inc + net_worth - ttaxent(inc)
                            !endif
                        else
                            print*, "error in taxation 1014-2"                            
                        endif !mode6askid
                    else
                        inc = inc + bgp ! 4.9.2017
                        ! merge((1._wp+(1._wp-tausv)*rd)*k, (1._wp+(1._wp-tausv)*rd)*k+(1._wp+rd)*a ,a>0._wp)
                        !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       (1._wp+(1._wp-tausv)*rd)*k, intfund>0._wp), & ! removed 9-12-2017
                        
                        !! 10.13.2017 comment out
                        !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       merge((1._wp+(1._wp-tausv)*rd)*k, (1._wp+(1._wp-tausv)*rd)*k+(1._wp+rd)*a ,a>0._wp), intfund>0._wp), & ! added 9-12-2017
                        !                  merge((1._wp+(1._wp-tausv)*rd)*intfund,                                                                               0._wp, intfund>0._wp), zx==0) &
                        !      + (1._wp-deltah)*h - ttaxent(inc) + merge(transbeq,0._wp,printout17) ! 4.9.2017                    
                        
                        ! 10.13.2017
                        inc = inc + merge( (1._wp+(1._wp-tausv)*rd)*intfund, merge( 0._wp, (1._wp+rd)*a, a>0._wp), intfund>0._wp) & ! 10.13.2017
                              + (1._wp-deltah)*h - ttaxent(inc) + merge(transbeq,0._wp,printout17) ! 4.9.2017                    
                    endif ! printout25
                    
                    
                    ! 3.26.2017 From left to right are the cases listed below:
                    ! (1) intfund > 0._wp & zx==0 earn rd interest from intfund (a). no borrowing (0). tax on what earns. no investment. (needs to pay back, if any)
                    ! (2) intfund <= 0._wp & zx==0 earn rd interest what borrow and own (k-a)+a. did have borrowing (k-a). tax on what earns. (needs to pay back, if any)
                    ! (3) intfund > 0._wp & zx==1 earn rd interes from intfund. tax on what earns. do investment. (needs to pay back, if any)
                    ! (4) intfund <= 0._wp & zx==1 no rd interest. no tax needed. do investment. (needs to pay back, if any)
                    ! 3.26.2017 When zx==0, bgp (c_grs_mat)  represents the negative value of borrowing cost. 
                    
                    !capgain = merge( merge( rd*intfund, 0._wp, intfund>0._wp), merge(, , intfun>0._wp), zx==1) ! ---> -r(k-a) ! 3.26.2017 seems useless
                    ! inc = benefit + merge(0._wp,capgain,tausvflag) + c_grs_mat(ax,kx,zx,yx,t) ! taxable income (business profits + interest income + SS benefits)
                    !inc = benefit + merge( 0._wp, capgain, tausvflag) ! 3.5.2017 removed the business benefits.
                    !labdem = c_lab_vec(kx) ! [4] ! 3.26.2017 useless. 
                    !inc = inc + a + (1._wp-tausv)*merge(capgain, 0._wp, tausvflag) + (1._wp-deltah)*h + (1._wp-taubp)*c_grs_mat(ax, kx, zx, yx, t) - ttaxent(inc) + transbeq ! 3.26.2017 comment out
                    
                    ! 3.25.2017 stop here. needs to change the formulation of inc, capgain.
                    !inc = inc + a + (1._wp-tausv)*merge(capgain, 0._wp, tausvflag) + (1._wp-deltah)*h + (1._wp-taubp)*c_grs_mat(ax, kx, zx, yx, t) - ttaxent(inc) + transbeq ! 3.26.2017 comment out
                    
                    !! 3.8.2017 labor income taxes. removed. and also removed "ssbtax" in the "inc" equation below.
                    !if(zx==0)then
                    !    ssbtax = 0._wp ! [5]
                    !else
                    !    ssbtax = (tauss/2._wp)*wage*labdem ! [6] ! 3.5.2017 divided tauss by 2. match up the SS tax (payroll tax) of his employees.
                    !endif
                    
                    ! 3.5.2017 remove ssbtax. avoiding being levied SS tax twice (the first time is levied in the form formulated in entrepreneurial business profit equation.
                    ! inc = inc + a + (1._wp-tausv)*merge(capgain, 0._wp, tausvflag) + (1._wp-deltah)*h + (1._wp-taubp)*c_grs_mat(ax, kx, zx, yx, t) - ttaxent(inc) + transbeq - ssbtax ! [7] 3.8.2017 remove "ssbtax"
                    ! c_grs_mat in ln. 1449.
                    
                    !if(zx==0)then
                    !    if(a-k>=0._wp)then ! No external financing. Just save in the banking account instead of investing in his business project.
                    !        !inc = inc + a + (1._wp-tausv)*merge(capgain, 0._wp, tausvflag) + (1._wp-deltah)*h + (1._wp-taubp)*c_grs_mat(ax, kx, zx, yx, t) - ttaxent(inc) + transbeq            
                    !        ! 3.25.2017 Still have capital income from saving.
                    !        inc = inc + a + (1._wp-tausv)*merge(capgain, 0._wp, tausvflag) + (1._wp-deltah)*h + c_grs_mat(ax, kx, zx, yx, t) - ttaxent(inc) + transbeq            
                    !    else ! a-k<0._wp have borrowing.
                    !        ! 3.25.2017 Do not have capital income from saving.
                    !        inc = inc + a + (1._wp-tausv)*merge(capgain, 0._wp, tausvflag) + (1._wp-deltah)*h + c_grs_mat(ax, kx, zx, yx, t) - ttaxent(inc) + transbeq         
                    !    endif
                    !else ! zx==1
                    !    if(a-k>=0._wp)then
                    !        
                    !    else ! a-k<0._wp
                    !        
                    !    endif                        
                    !endif
                    
                    ! 3.9.2017 I added + rd*merge(0._wp,a,a>=0._wp), although it seems to be redundant for cases of entrepreneurs (since their borrowing costs reflfect in their business profits).
                    ! 3.12.2017 confirmed. There is no need to add the borrowing cost for the reason stated in the line above.
                endif              
                
                ! (2) Portfolio next period
                if(inc>0._wp)then    
                    ! new zone -------------------------------------------- I
                    stp1  = 1
                    log1  = .false. ! ln. 316. 
                    svcmp = 100._wp ! Store up-to-date optimal value for convergence comparison
                    allocate( shv(nsdim), lowap(nsdim), highap(nsdim), smaxv(nsdim), apvec(nsdim) ) ! 3.10.2017 First run, "dsdim" takes nsdim; afterwards, "dsdim" takes sdim.
                    ! 3.10.2017 "nsdim" has better to be bigger than "hdim," I think.
                    do while(stp1<iterasvmax .and. log1==.false.) ! log1=.true. when only penalty can be found.
                        
                        ! ------------------- housing assset holding grid zooming-in -------------------------- 8-13-2017
                        if(stp1>1)then ! "dsdim" is safe regardless of the value of stp1. confirmed.
                            best = shv(loc1(1)) ! obtain the maximum from the previous round.
                            !if(printout5) write(unit=109,fmt='(a,2i4,f12.8)'), ' new stp old idx ', stp1, loc1(1), best
                            
                            ! It is the region to define the new search area for this round.
                            if(loc1(1)/=1 .and. loc1(1)/=dsdim)then ! Note: "dsdim" is initialized in the first round defined above (stp1==1). 
                                bdis = abs(shv(loc1(1)+1)-shv(loc1(1)-1)) ! Just expanded the search area centered at the maximizer of last round
                            else ! if it is corner solution (I think it shouldn't happen, because the algorithm always expands evenly around the maximizer after the first search except for the exception that a corner solution happens in the first round.
                                if(loc1(1)==1)then
                                    bdis = abs(shv(1)-shv(2))     
                                else ! loc1(1)==dsdim
                                    bdis = abs(shv(dsdim)-shv(dsdim-1))
                                endif
                            endif
                            
                            ! set up the exact boundary points for the new search in the current round.
                            if(loc1(1)==1)then ! In this case, we only loose up the lower boundary a little bit, as opposed to the relatively large expansion for the upper boundary point.
                                !bestl = best
                                bestl = best - 0.2_wp*bdis ! 10012016
                                bestu = best + 0.8_wp*bdis
                                if(bestl<hv(1)) bestl = hv(1) ! 10012016
                            elseif(loc1(1)==dsdim)then 
                                !bestu = best
                                bestu = best + 0.2_wp*bdis ! 10012016
                                bestl = best - 0.8_wp*bdis
                                if(bestu>hv(hdim)) bestu = hv(hdim) ! 10012016
                            else ! a normal case.
                                bestl = best - 0.8_wp*bdis
                                bestu = best + 0.8_wp*bdis
                                
                                if(bestl<hv(1)) bestl = hv(1) ! 10012016 Just in case the lower boundary poitn transpass the lower limit of our grid.
                                if(bestu>hv(hdim)) bestu = hv(hdim) ! 10012016 Similar reason to that above.
                            endif

                            dsdim = sdim ! smaller (5) 3.10.2017 zoom-in search area (using less grid points by default).
                            deallocate( shv, lowap, highap, smaxv, apvec )
                            allocate( shv(sdim), lowap(sdim), highap(sdim), smaxv(sdim), apvec(sdim) ) ! small region brent's method --- part II

                            ! call grid(shv,bestl,bestu,1._wp) 
                            !call grid(shv,bestl,bestu,1._wp) ! 3.8.2017 changed the 4th argument from 1._wp to gsdim. 3.13.2017 It may be better to keep using 1, since we like the new round covers the old optimal value and evenly spaced around it.
                            call grid(shv,bestl,bestu,gsdim) !10.26.2017
                            
                        else ! stp1 = 1 First run. So we just generate the financial asset vector that covers the full range of financial asset holdings. 3.5.2017
                            dsdim = nsdim ! In the first round, we use larger number of grid points. (20)          
                            if(hx==1)then ! given a level of housing asset, generate a grid for comparison of results from different financial asset holdings. 3.5.2017
                                call grid(shv,hv(1),hv(hdim),gsdim) ! "shv" is the vector of next period's housing asset holding.    
                            else ! hx
                                !if(cwh(ax,hx-1,kx,zx,yx,kpx,ypx,opx,t)/=penalty)then
                                !    call grid(shv,cwh(ax,hx-1,kx,zx,yx,kpx,ypx,opx,t),hv(hdim),gsdim)
                                !else
                                    call grid(shv,hv(1),hv(hdim),gsdim) ! "shv" is the vector of next period's housing asset holding. The default full range.
                                !endif ! cwh                                 
                            endif ! hx
                        endif ! stp                               
                        
                        smaxv  = penalty ! result vector for comparison
                        lowap  = -1 ! lower bound of ap
                        highap = -1 ! upper bound of ap

                        !if(printout5) write(unit=109,fmt='(a,9a)'), ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                        !if(printout5) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,kpx,ypx,opx,t
                        
                        ! ------------------------------- Brent zone -------------------------------------- 8-13-2017
                        do m = 1, dsdim
                            hp = shv(m) ! 3.6.2017 needed in brent_localizer's evoked function; pass into the thread "shv(m)" is the trial next period housing asset holding. 
                            call col_valid_fin(shv(m),kpx,ypx,t,lowap(m),deb1) ! output: lowap(m)
                            !if(printout5) write(unit=109,fmt='(a,i3,x,3f8.4)'), ' -c ', m, deb1 ! this printout can show how eager people want to borrow money.
                            call csp_valid_fin(h,shv(m),inc,highap(m))    ! 3.6.2017 actually, the output "highap(m)" is the maximum consumption rather than the maximum financial asset. 
                            
                            if(lowap(m)>=highap(m).or.highap(m)<0._wp)then ! 3.6.2017 the maximum consumption should not be zero.
                                !if(printout5) write(unit=109,fmt='(a,i4,f12.8,2f8.4,i5)'), ' -- ', m, shv(m), lowap(m), highap(m), merge(0,1,lowap(m)<highap(m))                                
                                cycle ! negative consumption in one of the refined housing assets level 
                            endif
                            
                            call brent_localizer(f14,lowap(m),highap(m),apvec(m),smaxv(m),exit_brent) ! 8-13-2017 [6]
                            !if(exit_brent) exit_bellman=.true. ! 8-13-2017 [7]
                            !if(exit_brent) write(*,*) ax, hx, kx, zx, yx, kpx, ypx, opx
                            
                            if(exit_brent) exit ! exit this loop. ! 8-13-2017 [7]
                            
                            !if(printout5) write(unit=109,fmt='(a,i4,f12.8,2f8.4,f15.8)'), ' -- ', m, shv(m), lowap(m), highap(m), smaxv(m)                            
                        enddo
                        
                        if(exit_brent==.false.)then ! 8-13-2017 [8]
                            ! Get the maximum from the set of local bests.
                            if(maxval(smaxv)/=penalty)then
                                !if(printout7) write(unit=109,fmt='(a,3f12.8)'), ' Gd ', apvec(maxloc(smaxv)), shv(maxloc(smaxv)), smaxv(maxloc(smaxv))                            
                                loc1 = maxloc(smaxv)
                                stp1 = stp1 + 1
                                if(abs(svcmp-smaxv(loc1(1)))<err_svdist) exit ! 3.8.2017 this one controls the convergence 
                                svcmp = smaxv(loc1(1)) ! update the best value.                            
                            else
                                log1 = .true. ! Exit. It indicates that the algorithm can not find a reasonalbe value other than the penalty.
                                svcmp = penalty
                            endif
                        else                                                           ! 
                            exit ! exit this housing asset holding grid searching zone ! 8-13-2017 [9]
                        endif ! exit_brent                                             !
                    enddo ! do while
                    ! end of new zone -------------------------------------------- I
                    
                    if(exit_brent==.false.)then ! 8-13-2017 [10]
                        ! save searching results or give a warning mark.
                        if(svcmp/=penalty)then
                            cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = smaxv(loc1(1))
                            cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = apvec(loc1(1))
                            cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = shv(loc1(1)) 
                            cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = inc - (apvec(loc1(1))+shv(loc1(1))) - HomTrnsCst(h,shv(loc1(1))) ! It should correspond to the consumption defined in function f14. 
                            !! CHECKED. NO LEAKS. 09302016
                            !if(printout5) write(unit=109,fmt='(a,9a)'), ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                            !if(printout5) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,kpx,ypx,opx,t                        
                            !if(cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t)<0._wp.and.printout5) write(unit=109,fmt='(a,f8.4)'), ' income ', cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t)
                            cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 0 ! die next period no matter how the current (last period) work status.
                            cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = -99 ! sign. -99 indicates normal.
                        else ! infeasible point due to collateral and budget constraints AND no downgrade place to go.
                            cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 2 
                        endif 
                    endif ! exit_brent 8-13-2017 [11]
                    
                    deallocate( shv, lowap, highap, smaxv, apvec ) 
                else ! NEGATIVE income (stops before hitting the criteria of negative consumption).
                    !if(printout7) write(unit=109,fmt='(a,10a)'), ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t',' insufficeint income '
                    !if(printout7) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,kpx,ypx,opx,t                    
                    cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 1 ! if cww takes a value other than -99,                
                endif ! income > 0   
                
                if(exit_brent) exit ! exit this n loop (over the indices of housing asset holding grid). 8-13-2017 [12]
                
            enddo ! n
            
            !if(exit_brent) write(*,*) ' q ', q, t
            if(exit_brent) qray(q)=.true. ! 8-13-2017 [13]
            
        enddo ! q
        !$omp end do
        !$omp end parallel ! 3.11.2017 There is no "psz" and "v2d" and "ap, apx, hpx, zpx" and "i, opxi" are used in this parallel block 
        ! ------- End of the end of period 14 computation.
        
        if(any(qray==.true.)) exit_bellman=.true. ! 8-13-2017 [14]
        deallocate( qray ) ! 8-13-2017 [15]
        
        if(exit_bellman==.false.)then ! 8-13-2017 [16]
            t = 14 
            ! ------- PRINT OUT 
            if(printout2)then        
                do q = bd(1,1), bd(1,2) ! income in period 14 could be used later.   
                    ax  = xlist(q,1)
                    hx  = xlist(q,2)
                    kx  = xlist(q,3)
                    zx  = xlist(q,4)
                    yx  = xlist(q,5)
                    kpx = xlist(q,6)
                    ypx = xlist(q,7)
                    opx = xlist(q,8)       
                         
                    if(ax==1.and.hx==1)then
                        write(unit=056,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                        write(unit=056,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                        write(unit=056,fmt='(3x,9x,'//str1//'i21)') (j,j=1,hdim) ! "str1" represents the size of housing asset grid.
                        write(unit=056,fmt='(3x,9x,'//str1//'f21.4)') (hv(j),j=1,hdim) 
                        
                        do i = 1, adim
                            write(unit=056,fmt='(i3,x,f8.4,'//str1//'(x,f20.8))') i,av(i),(cwf(i,j,kx,zx,yx,kpx,ypx,opx,t),j=1,hdim)    
                        enddo
                    
                    endif        
                enddo ! q
            endif ! printout2
            
            ! beginning of period 14. Ok consistent 6-17-2016 
            ! Just copy the end value function, because actually we only need the beginning period of value function for period 14
            !allocate( tvector(7,3) )
            
            ! Beginning of the period --- important     
            do n = 1, 7
                kx = tvector(n,1)
                zx = tvector(n,2)
                yx = tvector(n,3) ! all zeros, because no labor efficiency exists any more.
                do j = 1, hdim
                    do i = 1, adim
                        wf(i,j,kx,zx,yx,t) = cwf(i,j,kx,zx,yx,0,0,0,t)
                    enddo
                enddo   
                
                if(printout2)then
                    write(unit=026,fmt='(a,3(i3),a,i3)') '(kx,zx,yx) ',kx,zx,yx, ' t ',t     
                    write(unit=026,fmt='(3x,8x,'//str1//'(8x,i8))') (j,j=1,hdim)             
                    write(unit=026,fmt='(3x,5x,"999",'//str1//'(8x,f8.4))') (hv(j),j=1,hdim)             
                    do j = 1, adim
                        write(unit=026,fmt='(i3,x,f7.2,'//str1//'(x,f15.5))') j, av(j),(wf(j,i,kx,zx,yx,t),i=1,hdim)
                    enddo                                  
                endif ! printout2    
            enddo
        endif ! eixt_bellman ! 8-13-2017 [17]
            
        deallocate( xlist, bd )    
        
        ! No borrowing, no new business idea in the last peirod, so there is no probablility transition matrix. The beginning is the same as the         
        ! end of period 14 calculation        
            
        !! EXAM THE COARSE VALUE FUNCTION (WORKERS WHO HAVE LARGER HOUSES ARE ABLE TO BORROW)
        !write(unit=006,fmt='(a,2(i3))') '--- matrix (ax,hx) for beginning-of-period value function at (kx,zx) =',0,0
        !do i = 1, adim
        !    write(unit=006,fmt=3006) i,(twf(i,j,0,0,0,0,0,0,14),j=1,10)
        !enddo        
        !do m =  1, 3
        !    do n = 0, 1
        !        write(unit=006,fmt='(a,2(i3))') '--- matrix (ax,hx) for beginning-of-period value function at (kx,zx) =',m,n 
        !        write(unit=006,fmt='(4x,10(13x,i3))') (i,i=1,10)
        !        do i = 1, adim
        !            write(unit=006,fmt=3006) i,(twf(i,j,m,n,0,0,0,0,14),j=1,10)
        !        enddo
        !    enddo
        !enddo   
        
        call system_clock(tend)
        !if(printout6) write(unit=120,fmt='(a,i3,a,f12.4,a)') 'solve bellman period ',t,', time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'         
        if(printout6) write(*,fmt='(/,a,i3,a,f12.4,a)') 'solve bellman period ',t,', time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'         

        !! t = 10-13 ------------------------------------------------------------------------------------------------------
        
        if(exit_bellman==.false.)then ! 8-13-2017 [b0]
            allocate( xlist(adim*1*12,8), bd(4,2) )
            xlist = xl1013 ! 3.10.2017
            bd = bd1013
            
            allocate( qray(bd(4,2)) ) ! 8-13-2017 [b1]
            
            do t = 13, 10, -1
                qray = .false. ! 8-13-2017 [b2]
                call system_clock(tstart)
                do l = 1, 4 ! 4.22.2017
                !do l = 1, 1
                    !!$omp parallel default(shared) private(n,capgain,intfund,i,stp1,log1,svcmp,shv,lowap,highap,smaxv,apvec,best,loc1,bdis,bestl,bestu,dsdim,m,opxi,labsup,labdem,bgp,exit_brent) ! #1 8-13-2017 [b3]
                    allocate( v2d(adim,hdim,1:2) ) ! 3.11.2017 Note: "v2d" is passed into respective threads for interpolation!! 
                    !!$omp do ! #2
                    do q = bd(l,1), bd(l,2) ! 4.22.2017
                    !do q = 18, 18
                        do n = 1, hdim
                            
                            exit_brent = .false. ! 8-13-2017 [b4]
                            
                            ax  = xlist(q,1)
                            hx  = n
                            kx  = xlist(q,3)
                            zx  = xlist(q,4)
                            yx  = xlist(q,5)
                            kpx = xlist(q,6)
                            ypx = xlist(q,7)
                            opx = xlist(q,8)                          
                            
                            a  = av(ax)
                            h  = hv(hx)
                            k  = kv(kx)
                            z  = merge(0._wp, merge(zlow,z2(kx),zx==0), kx==0)
                            y  = yv(yx)
                            kp = kv(kpx)
                            yp = yv(ypx)
                            dk = merge(0._wp, merge(delzl(kx), delzh(kx), zx==0), kx==0)                          
                            
                            ! income (step 1)
                            if(kx==0)then
                                if(printout25)then
                                    ! 10.26.2017 revision
                                    !inc = benefit + merge(transbeq, 0._wp, printout17)
                                    if(mode6taskid==0)then
                                        inc = merge(benefit+merge(transbeq,0._wp,printout17), &
                                                    benefit+merge(transbeq,0._wp,printout17)+merge(rd*a,0._wp,a>0._wp), &
                                                    tausv>0._wp) ! taxable 3-11-2018
                                        inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a, (1._wp+rd)*a, a>=0._wp), &
                                                          merge(a,(1._wp+rd)*a,a>0._wp), &
                                                          tausv>0._wp) & ! change in net worth 3-11-2018
                                              + (1._wp - deltah) * h - ttaxwok(inc)                                   
                                    elseif(mode6taskid==1)then ! wealth tax (not applied on period income but beginning capital, regardless of asset types).
                                        !inc = inc + (1._wp-taubal)*merge( (1._wp+rd)*a, (1._wp+(1._wp-tausv)*rd)*a, a<0._wp) + (1._wp-taubal)*(1._wp-deltah)*h - ttaxwok(inc)  
                                        inc = benefit + merge(transbeq, 0._wp, printout17) + (1._wp + rd) * a + (1._wp - deltah) * h 
                                        inc = (1._wp - taubal) * inc
                                        !inc = inc + (1._wp - taubal) * merge( 0._wp, (1._wp + (1._wp - tausv) * rd) * a, a < 0._wp) &
                                        !    + (1._wp - taubal) * (1._wp - deltah) * h - ttaxwok(inc) + merge(0._wp, (1._wp + (1._wp - tausv) * rd) * a, a > 0._wp) 
                                    elseif(mode6taskid==2)then
                                        net_worth = benefit + merge(transbeq,0._wp,printout17) + (1._wp + rd) * a + (1._wp - deltah) * h
                                        inc = benefit+merge(transbeq,0._wp,printout17)+merge(rd*a,0._wp,a>0._wp) 
                                        if(net_worth > 5._wp*exempbar)then
                                            inc = inc + merge(a,(1._wp+rd)*a,a>0._wp) + (1._wp-deltah)*h - ttaxwok(inc) - 0.02_wp*(net_worth-5._wp*exempbar) - 0.01_wp*(4._wp*exempbar)   
                                        elseif(net_worth>exempbar .and. net_worth<=5._wp*exempbar)then
                                            inc = inc + merge(a,(1._wp+rd)*a,a>0._wp) + (1._wp-deltah)*h - ttaxwok(inc) - 0.01_wp*(net_worth - exempbar)   
                                        else
                                            inc = inc + merge(a,(1._wp+rd)*a,a>0._wp) + (1._wp-deltah)*h - ttaxwok(inc)
                                        endif                                      
                                    else
                                        print*, "error in taxation 1014-3"    
                                    endif !mode6taskid
                                else
                                    ! 3.26.2017 revision
                                    inc = benefit
                                    inc = inc + merge( (1._wp+rd)*a, (1._wp+(1._wp-tausv)*rd)*a, a<0._wp) + (1._wp-deltah)*h - ttaxwok(inc) + merge(transbeq,0._wp,printout17)
                                endif !printout25
                                
                                !capgain = merge(rd*a,0._wp,rd*a>0._wp) ! retirees can only have capital income from savings.
                                !inc = benefit + merge(0._wp,capgain,tausvflag)
                                !inc = inc + a + rd*merge(a, 0._wp, a<0._wp) + (1._wp-tausv)*merge(capgain, 0._wp, tausvflag) + (1._wp-deltah)*h - ttaxwok(inc) + transbeq  ! total wealth before the introduction of the adjustment costs of housing units
                            else
                                labsup = 0._wp ! [3] ! 3.26.2017 useless in this lifecycle stage.
                                intfund = merge(a-k, 0._wp, a>k) ! excess internal fund bears 
                                !inc = benefit
                                bgp = c_grs_mat(ax, kx, zx, yx, t)
                                
                                !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       (1._wp+(1._wp-tausv)*rd)*k, intfund>0._wp), &
                                !                  merge((1._wp+(1._wp-tausv)*rd)*intfund,                      0._wp, intfund>0._wp), zx==0) &
                                !      + merge((1._wp-taubp)*bgp, bgp, zx==1) + (1._wp-deltah)*h - ttaxent(inc) + transbeq ! 3.26.2017                            
                                
                                if(printout25)then
                                    !inc = inc + bgp + merge(transbeq, 0._wp, printout17) ! 10.26.2017      
                                    !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       (1._wp+(1._wp-tausv)*rd)*k, intfund>0._wp), & ! removed 9-12-2017  
                                    
                                    !!!10.13.2017
                                    !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       merge((1._wp+(1._wp-tausv)*rd)*k, (1._wp+(1._wp-tausv)*rd)*k+(1._wp+rd)*a ,a>0._wp), intfund>0._wp), & ! added 9-12-2017
                                    !                  merge((1._wp+(1._wp-tausv)*rd)*intfund,                                                                               0._wp, intfund>0._wp), zx==0) &
                                    !      + (1._wp-deltah)*h - ttaxent(inc) + merge(transbeq, 0._wp, printout17) ! 3.26.2017                               
                                    
                                    !10.13.2017
                                    if(mode6taskid==0)then
                                        inc = bgp + merge(benefit+merge(transbeq,0._wp,printout17), &
                                                          benefit+merge(transbeq,0._wp,printout17)+merge(rd*intfund,0._wp,intfund>0._wp), &
                                                          tausv>0._wp)
                                        inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*intfund, 0._wp, intfund>=0._wp), &
                                                          merge(intfund,0._wp,intfund>0._wp), &
                                                          tausv>0._wp) & ! change in net worth 3-11-2018
                                              + (1._wp - deltah) * h - ttaxent(inc)                           
                                    elseif(mode6taskid==1)then ! wealth tax (not applied on period income but beginning capital, regardless of asset types).
                                        inc = bgp + benefit + merge(transbeq, 0._wp, printout17) ! 10.26.2017 
                                        !inc = inc + (1._wp - taubal) * (merge(bgp, 0._wp, bgp>0._wp) + merge((1._wp + (1._wp - tausv) * rd) * intfund, 0._wp, intfund > 0._wp)) & ! 10.13.2017
                                        !      + (1._wp - taubal) * (1._wp - deltah) * h - ttaxent(inc) + merge(bgp, 0._wp, bgp<=0._wp)    
                                        inc = inc + merge((1._wp + rd) * intfund, 0._wp, intfund > 0._wp) + (1._wp - deltah) * h
                                        inc = (1._wp - taubal) * inc
                                    elseif(mode6taskid==2)then
                                        net_worth = bgp + benefit + merge(transbeq,0._wp,printout17) + merge((1._wp + rd) * intfund, 0._wp, intfund>0._wp) + (1._wp - deltah) * h
                                        inc = bgp + benefit + merge(transbeq,0._wp,printout17) + merge(rd*intfund,0._wp,intfund>0._wp) 
                                        if(net_worth > 5._wp*exempbar)then
                                            inc = inc + merge(intfund,0._wp,intfund>0._wp) + (1._wp-deltah)*h - ttaxent(inc) - 0.02_wp*(net_worth-5._wp*exempbar) - 0.01_wp*(4._wp*exempbar)   
                                        elseif(net_worth>exempbar .and. net_worth<=5._wp*exempbar)then
                                            inc = inc + merge(intfund,0._wp,intfund>0._wp) + (1._wp-deltah)*h - ttaxent(inc) - 0.01_wp*(net_worth - exempbar)   
                                        else
                                            inc = inc + merge(intfund,0._wp,intfund>0._wp) + (1._wp-deltah)*h - ttaxent(inc)
                                        endif  
                                    else
                                        print*, "error in taxation 1014-4"    
                                    endif
                                else
                                    inc = inc + bgp ! 4.10.2017
                                    !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       (1._wp+(1._wp-tausv)*rd)*k, intfund>0._wp), & ! removed 9-12-2017  
                                    
                                    !!!10.13.2017
                                    !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       merge((1._wp+(1._wp-tausv)*rd)*k, (1._wp+(1._wp-tausv)*rd)*k+(1._wp+rd)*a ,a>0._wp), intfund>0._wp), & ! added 9-12-2017
                                    !                  merge((1._wp+(1._wp-tausv)*rd)*intfund,                                                                               0._wp, intfund>0._wp), zx==0) &
                                    !      + (1._wp-deltah)*h - ttaxent(inc) + merge(transbeq, 0._wp, printout17) ! 3.26.2017                               
                                    
                                    !10.13.2017
                                    inc = inc + merge( (1._wp+(1._wp-tausv)*rd)*intfund, merge( 0._wp, (1._wp+rd)*a, a>0._wp), intfund>0._wp) & ! 10.13.2017
                                          + (1._wp-deltah)*h - ttaxent(inc) + merge(transbeq, 0._wp, printout17) ! 3.26.2017                                                               
                                endif !printout25
                                    
                                !inc = inc + a + (1._wp-tausv)*merge(capgain,0._wp,tausvflag) + (1._wp-deltah)*h + (1._wp-taubp)*c_grs_mat(ax, kx, zx, yx, t) - ttaxent(inc) + transbeq
                                !capgain = merge(rd*intfund,0._wp,rd*intfund>0._wp)
                                ! inc = benefit + merge(0._wp,capgain,tausvflag) + c_grs_mat(ax,kx,zx,yx,t) ! taxable income
                                
                                !labdem = c_lab_vec(kx) ! [4] ! 0204pm stop at here. ! 3.26.2017 seems useless.
                                !! 3.10.2017 Social security tax has been taxed in business profits, so I removed ssbtax below when calculating the disposable income.
                                !if(zx==0)then
                                !    ssbtax = 0._wp ! [5]   
                                !else
                                !    ssbtax = tauss*wage*labdem ! [6]
                                !endif
                                
                                !inc = inc + a + (1._wp-tauk)*merge(capgain,0._wp,tausvflag) + (1._wp-deltah)*h - ttaxent(inc) + transbeq - ssbtax ! [7]                         
                                ! inc = inc + a + rd*merge(a,0._wp,a<0._wp.and.zx==0) + (1._wp-tausv)*merge(capgain,0._wp,tausvflag) + (1._wp-deltah)*h + (1._wp-taubp)*c_grs_mat(ax, kx, zx, yx, t) - ttaxent(inc) + transbeq
                                ! 3.12.2017 entrepreneur's disposable income shouldn't include + rd*merge(a,0._wp,a<0._wp.and.zx==0)
                                ! For one, there is no possible for him to run a zero scale project (kx/=0 anyway in this block), and for two, the corresponding situation is that
                                ! in this condition block he either earns riskless interest rates (which is a cased covered in the disposable income equation above) or he borrows and pay the interests
                                ! through his business profits equation. Therefore, there is no need to consider the term, + rd*merge(a,0._wp,a<0._wp.and.zx==0), in any case.
                            endif  
                            
                            if(inc>0._wp)then
                                
                                if(opx/=0)then ! Entrepreneurial case. The agent will receive idiosyncratic business shock.
                                    pzv(1) = pz2(kpx,1) ! probability of technology shock is determined by the size of next-period project (kpx).
                                    pzv(2) = pz2(kpx,2) ! probability of technology shock is determined by the size of next-period project (kpx).        
                                    do i = 1, 2 ! There are two types of business shock. Note that this variable is not the state index of business shocks, which ranges from 0 to 1, rather than 1 to 2.
                                        v2d(:,:,i) = wf(:,:,kpx,i-1,0,t+1)    
                                    enddo
                                else ! For workers, there is no business shock, so the value function reduces to a single level in the business shock dimension.
                                    v2d(:,:,1) = wf(:,:,kpx,0,0,t+1)    
                                endif ! opx
                                
                                stp1  = 1
                                log1  = .false.
                                svcmp = 100._wp ! used for convergence
                                allocate( shv(nsdim), lowap(nsdim), highap(nsdim), smaxv(nsdim), apvec(nsdim) ) ! first run
                                do while(stp1<iterasvmax .and. log1==.false.)
                                    
                                    ! ACCURACY BLOCK --- START
                                    if(stp1>1)then
                                        best = shv(loc1(1))
                                        !if(printout7) write(unit=109,fmt='(a,2i4,f12.8)'), ' new stp old idx ', stp1, loc1(1), best
                                        if(loc1(1)/=1.and.loc1(1)/=dsdim)then
                                            bdis = abs(shv(loc1(1)+1)-shv(loc1(1)-1))
                                        else
                                            if(loc1(1)==1)then
                                                bdis = abs(shv(1)-shv(2))       
                                            else ! loc1(1)==dsdim
                                                bdis = abs(shv(dsdim)-shv(dsdim-1))
                                            endif
                                        endif
                                        
                                        !if(loc1(1)==1)then
                                        !    bestl = best
                                        !    bestu = best + 0.8_wp*bdis
                                        !elseif(loc1(1)==dsdim)then
                                        !    bestu = best
                                        !    bestl = best - 0.8_wp*bdis
                                        !else
                                        !    bestl = best - 0.8_wp*bdis
                                        !    bestu = best + 0.8_wp*bdis
                                        !endif

                                        if(loc1(1)==1)then
                                            !bestl = best
                                            bestl = best - 0.2_wp*bdis ! 10012016
                                            bestu = best + 0.8_wp*bdis
                                            if(bestl<hv(1)) bestl = hv(1) ! 10012016
                                        elseif(loc1(1)==dsdim)then
                                            !bestu = best
                                            bestu = best + 0.2_wp*bdis ! 10012016
                                            bestl = best - 0.8_wp*bdis
                                            if(bestu>hv(hdim)) bestu = hv(hdim) ! 10012016
                                        else
                                            bestl = best - 0.8_wp*bdis
                                            bestu = best + 0.8_wp*bdis
                                            if(bestl<hv(1)) bestl = hv(1) ! 10012016
                                            if(bestu>hv(hdim)) bestu = hv(hdim) ! 10012016
                                        endif                                    
                                        
                                        dsdim = sdim ! zoom-in search area (using less grid points by default).
                                        deallocate( shv, lowap, highap, smaxv, apvec )
                                        allocate( shv(sdim), lowap(sdim), highap(sdim), smaxv(sdim), apvec(sdim) ) ! small region brent's method --- part II
                                        !if(printout7) write(unit=109,fmt='(a,3f8.4)'), ' boundary ', bestl, bestu, bdis
                                        !call grid(shv,bestl,bestu,1._wp) ! 3.13.2017 keep using 1._wp rather than gsdim.
                                        call grid(shv,bestl,bestu,gsdim) ! 3.13.2017 keep using 1._wp rather than gsdim.
                                    else ! stp1 = 1
                                        dsdim = nsdim ! bigger grid              
                                        if(hx==1)then
                                            call grid(shv,hv(1),hv(hdim),gsdim)    
                                        else ! hx
                                            !if(cwh(ax,hx-1,kx,zx,yx,kpx,ypx,opx,t)/=penalty)then
                                            !    call grid(shv,cwh(ax,hx-1,kx,zx,yx,kpx,ypx,opx,t),hv(hdim),gsdim)
                                            !else
                                                call grid(shv,hv(1),hv(hdim),gsdim) ! the default full range
                                            !endif ! cwh                                 
                                        endif ! hx
                                    endif ! stp
                                    ! ACCURACY BLOCK --- END                               
                                    
                                    smaxv = penalty
                                    lowap  = -1
                                    highap = -1
                                    
                                    ! 4.23.2017 # 1 & 2
                                    if(printout7) write(unit=109,fmt='(a,9a)'),  ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                                    if(printout7) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax, hx, kx, zx, yx, kpx, ypx, opx, t
                                    
                                    do m = 1, dsdim
                                        if(printout7) write(unit=109,fmt='(i4)') m ! 4.23.2017 # 3
                                        
                                        hp = shv(m) ! 3.12.2017 Pass into the thread. Note that hp keeps being refined throughout the process.
                                        call col_valid_fin(shv(m),kpx,ypx,t,lowap(m),deb1)
                                        !if(printout7) write(unit=109,fmt='(a,i3,x,3f8.4)'), ' -c ', m, deb1 ! this printout can show how eager people want to borrow money.
                                        call csp_valid_fin(h,shv(m),inc,highap(m))      
                                        
                                        if(lowap(m)>=highap(m).or.highap(m)<0._wp)then
                                            !if(printout7) write(unit=109,fmt='(a,i4,f12.8,2f8.4,i5)'), ' -A ', m, shv(m), lowap(m), highap(m), merge(0,1,lowap(m)<highap(m))
                                            cycle ! negative consumption in one of the refined housing assets level 
                                        endif
                                        call brent_localizer(f1013,lowap(m),highap(m),apvec(m),smaxv(m),exit_brent) ! <-------------------------------------------- ! 8-13-2017 [b5]
                                        !if(printout7) write(unit=109,fmt='(a,i4,f12.8,2f8.4,e13.4)'), ' -B ', m, shv(m), lowap(m), highap(m), smaxv(m)
                                        
                                        if(exit_brent) exit ! exit this loop over m. ! 8-13-2017 [b6]
                                        
                                        if(printout7) write(unit=109,fmt='(f15.4)') smaxv(m) ! 4.23.2017 # 4
                                        !if(smaxv(m)<penalty/1.e5) smaxv(m)=penalty ! 4.23.2017 <---------------------------
                                    enddo
                                    
                                    if(exit_brent==.false.)then ! 8-13-2017 [b7]
                                        ! get the maximum from the set of local bests.
                                        if(maxval(smaxv)/=penalty)then
                                            !if(printout7) write(unit=109,fmt='(a,3f12.8)'), ' Gd ', apvec(maxloc(smaxv)), shv(maxloc(smaxv)), smaxv(maxloc(smaxv))
                                            loc1 = maxloc(smaxv)
                                            stp1 = stp1 + 1
                                            if(abs(svcmp-smaxv(loc1(1)))<err_svdist) exit
                                            svcmp = smaxv(loc1(1))                            
                                        else
                                            log1 = .true.
                                            svcmp= penalty
                                            !if(printout7) write(unit=109,fmt='(a)'), ' Ba '
                                        endif
                                    else                                                           !
                                        exit ! exit this housing asset holding grid searching zone ! 8-13-2017 [b8]    
                                    endif ! exit_brent                                             !
                                    
                                enddo ! do while      
                                
                                if(exit_brent==.false.)then ! 8-13-2017 [b9]
                                    ! save searching results or give a warning mark.
                                    if(svcmp/=penalty)then ! SOME CASES NEED TO DO COMPARISON 
                                        cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = smaxv(loc1(1))
                                        cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = apvec(loc1(1))
                                        cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = shv(loc1(1)) 
                                        cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = -99 ! -99 indicates a normal situation occurs.
                                        if(opx==0)then ! Once retired, retired forever.
                                            cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 0 
                                            cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = inc - (apvec(loc1(1))+shv(loc1(1))) - HomTrnsCst(h,shv(loc1(1))) 
                                            !! CHECKED. NO LEAKS. 09302016 OK!
                                            !if(cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t)<0._wp)then
                                            !if(printout5) write(unit=109,fmt='(a,9a,x,a)'), ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t', 'A'
                                            !if(printout5) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,kpx,ypx,opx,t                        
                                            !if(printout5) write(unit=109,fmt='(a,f15.4)'), ' income ', cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t)                                                                        
                                            !endif
                                        else ! opx/=0 3.30.2017 Cases where further comparison needs to make. 3.30.2017 Stop here. 11:00 am.
                                            cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = kpx ! The benchmark for comparison is the new business idea.
                                            cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = inc - (apvec(loc1(1))+shv(loc1(1))) - HomTrnsCst(h,shv(loc1(1)))
                                            if(opx==1)then ! 3.11.2017 Engage in the existing production technology.
                                                !if(printout5) write(unit=109,fmt='(a,9a,x,a)'), ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t', ' B'
                                                !if(printout5) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,kpx,ypx,opx,t                                          
                                                !if(printout5) write(unit=109,fmt='(2f15.4,"----------------",i3)'), cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t), cwf(ax,hx,kx,zx,yx,0,0,0,t), merge(1,0,cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t)>cwf(ax,hx,kx,zx,yx,0,0,0,t))
                                                if( cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) < cwf(ax,hx,kx,zx,yx,0,0,0,t))then ! 3.11.2017 Choose to be retired next period. 3.29.2017 If there is a tie, I assume the agent prefer to being an entrepreneur.
                                                    cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwf(ax,hx,kx,zx,yx,0,0,0,t) 
                                                    cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwa(ax,hx,kx,zx,yx,0,0,0,t)
                                                    cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwh(ax,hx,kx,zx,yx,0,0,0,t)
                                                    cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cww(ax,hx,kx,zx,yx,0,0,0,t)
                                                    cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwk(ax,hx,kx,zx,yx,0,0,0,t)
                                                    cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwc(ax,hx,kx,zx,yx,0,0,0,t)
                                                    !! CHECKED. NO LEAKS. 09302016 OK.
                                                    !if(cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t)<0._wp)then
                                                    !    if(printout5) write(unit=109,fmt='(a,9a,x,a)'), ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t', ' B'
                                                    !    if(printout5) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,kpx,ypx,opx,t   
                                                    !    if(printout5) write(unit=109,fmt='(2f15.4)'), cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t), cwf(ax,hx,kx,zx,yx,0,0,0,t) 
                                                    !    if(printout5) write(unit=109,fmt='(a,f15.4)'), ' income ', cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t)                                                                        
                                                    !endif                                       
                                                endif
                                            else ! opx==2 ! choose to engage in a smaller investment project.
                                                opxi = merge(opx-1,opx-2,kx/=0) ! Compared opx==2 (advanced project) with opxi (the default project in this period).
                                                if( cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) < cwf(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t))then
                                                    cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwf(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t) 
                                                    cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwa(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)
                                                    cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwh(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)
                                                    cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cww(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t) 
                                                    cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwk(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)                                            
                                                    cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwc(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t) 
                                                    !! CHECKED. NO LEAKS. 09302016
                                                    !if(cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t)<0._wp)then
                                                    !if(printout5) write(unit=109,fmt='(a,9a,x,a)'), ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t', ' C'
                                                    !if(printout5) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,kpx,ypx,opx,t                        
                                                    !if(printout5) write(unit=109,fmt='(a,9i4)'), ' source ', ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t     
                                                    !if(printout5) write(unit=109,fmt='(a,f15.4)'), ' source ', cwf(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)                                           
                                                    !if(printout5) write(unit=109,fmt='(a,f15.4)'), ' income ', cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t)                                                                                                           
                                                endif
                                            endif ! opx==1
                                        endif ! opx/=0
                                    else ! infeasible point due to collateral and budget constraints. DOWNGRADE WITHOUT ANY COMPARISON
                                        if(opx==0)then ! no where to downgrade their expenses or costs
                                            cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 2 ! warning singal is issued. That is, a bad state combination is found.
                                        elseif(opx==1)then ! may instead choose to be a worker
                                            cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwf(ax,hx,kx,zx,yx,0,0,0,t) 
                                            cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwa(ax,hx,kx,zx,yx,0,0,0,t)
                                            cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwh(ax,hx,kx,zx,yx,0,0,0,t)
                                            cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cww(ax,hx,kx,zx,yx,0,0,0,t)   
                                            cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwk(ax,hx,kx,zx,yx,0,0,0,t)                                    
                                            cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwc(ax,hx,kx,zx,yx,0,0,0,t)
                                            !! CHECKED. NO LEAKS. 09302016
                                            !if(cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t)<0._wp)then
                                            !if(printout5) write(unit=109,fmt='(a,9a,x,a)'), ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t', ' D'
                                            !if(printout5) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,kpx,ypx,opx,t   
                                            !if(printout5) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,0,0,0,t 
                                            !if(printout5) write(unit=109,fmt='(a,f15.4)'), ' source ', cwc(ax,hx,kx,zx,yx,0,0,0,t)
                                            !if(printout5) write(unit=109,fmt='(a,f15.4)'), ' income ', cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t)                                                                        
                                            !endif                                     
                                        elseif(opx==2)then ! may instead choose to engage in an investment project of smaller scale.
                                            opxi = merge(opx-1,opx-2,kx/=0)
                                            cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwf(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t) 
                                            cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwa(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)
                                            cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwh(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)
                                            cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cww(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)   
                                            cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwk(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)                                                                                
                                            cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwc(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)   
                                            !! CHECKED. NO LEAKS. 09302016
                                            !if(cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t)<0._wp)then
                                            !if(printout5) write(unit=109,fmt='(a,9a,x,a)'), ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t', ' E'
                                            !if(printout5) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,kpx,ypx,opx,t    
                                            !if(printout5) write(unit=109,fmt='(a,9i4)'), ' source ', ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t
                                            !if(printout5) write(unit=109,fmt='(a,f15.4)'), ' source ', cwc(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)                                       
                                            !if(printout5) write(unit=109,fmt='(a,f15.4)'), ' income ', cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t)                                                                        
                                            !endif                                       
                                        endif    
                                    endif 
                                endif ! exit_brent 8-13-2017 [b10]
                                
                                deallocate( shv, lowap, highap, smaxv, apvec )                             
                                
                            else ! inc<0._wp ! bug fixed 09232016
                                ! 3.11.2017 Note: there is no need to allocate matrices: shv, lowap, highap, smaxv, and apvec!!
                                cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 3    
                            endif ! inc
                            
                            if(exit_brent) exit ! exit this n loop (over the indices of housing asset hodling grid). 8-13-2017 [b11]
                            
                        enddo ! n
                        
                        if(exit_brent) qray(q)=.true. ! 8-13-2017 [b12]
                        
                    enddo ! q
                    !!$omp end do ! #3
                    deallocate( v2d )
                    !!$omp end parallel ! #4
                    
                    if(any(qray==.true.)) exit_bellman=.true. ! 8-13-2017 [b13]
                    if(exit_bellman) exit ! exit the loop over l ! 8-13-2017 [b14]
                enddo ! l
                
                if(exit_bellman==.false.)then ! 8-13-2017 [b15]
                    if(printout2)then 
                        do l = 1, 4                
                            do q = bd(l,1), bd(l,2) ! income in period 14 could be used later.
                                ax  = xlist(q,1)
                                hx  = xlist(q,2)
                                kx  = xlist(q,3)
                                zx  = xlist(q,4)
                                yx  = xlist(q,5)
                                kpx = xlist(q,6)
                                ypx = xlist(q,7)
                                opx = xlist(q,8)  
                                         
                                if(ax==1.and.hx==1)then ! printout only once
                                    write(unit=056,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                                    write(unit=056,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                                    write(unit=056,fmt='(3x,9x,'//str1//'i21)') (j,j=1,hdim) 
                                    write(unit=056,fmt='(3x,9x,'//str1//'f21.4)') (hv(j),j=1,hdim) 
                                    do i = 1, adim
                                        write(unit=056,fmt='(i3,x,f8.4,'//str1//'(x,f20.8))') i,av(i),(cwf(i,j,kx,zx,yx,kpx,ypx,opx,t),j=1,hdim)    
                                    enddo
                                endif 
                                
                            enddo ! q            
                        enddo ! l
                    endif ! printout2           
                    
                    ! 3.14.2017 The uncertainty comes solely from the innovation of business ideas.
                    do n = 1, 7 ! Beginning of the period. 3.12.2017 Grouping of combination at the "beginning" of period. 
                        kx = tvector(n,1)
                        zx = tvector(n,2)
                        yx = tvector(n,3)                
                        if(n<=4)then
                            where(cww(:,:,kx,zx,yx,0,0,0,t)==-99) & ! the condition is needed, so that we dont' need to screen out expectation of utils of unrealistic states
                                & wf(:,:,kx,zx,yx,t) = cwf(:,:,kx,zx,yx,0,0,0,t)
                        elseif(5<=n.and.n<=6)then
                            where(cww(:,:,kx,zx,yx,kx,0,1,t)==-99.and.cww(:,:,kx,zx,yx,kx+1,0,2,t)==-99) & ! two possibilities
                                & wf(:,:,kx,zx,yx,t) = (1._wp-probtk(kx))*cwf(:,:,kx,zx,yx,kx,0,1,t) + probtk(kx)*cwf(:,:,kx,zx,yx,kx+1,0,2,t) ! 3.12.2017 "probtk" is the probability to get new idea.
                        elseif(n==7)then
                            where(cww(:,:,kx,zx,yx,kx,0,1,t)==-99) &
                                & wf(:,:,kx,zx,yx,t) = cwf(:,:,kx,zx,yx,kx,0,1,t)
                        endif
                        
                        if(printout2)then
                            write(unit=026,fmt='(a,3(i3),a,i3)') '(kx,zx,yx) ',kx,zx,yx, ' t ',t     
                            write(unit=026,fmt='(3x,8x,'//str1//'(8x,i8))') (j,j=1,hdim)             
                            write(unit=026,fmt='(3x,5x,"999",'//str1//'(8x,f8.4))') (hv(j),j=1,hdim)             
                            do j = 1, adim
                                write(unit=026,fmt='(i3,x,f7.2,'//str1//'(x,f15.5))') j, av(j),(wf(j,i,kx,zx,yx,t),i=1,hdim)
                            enddo                                  
                        endif ! printout2
                    enddo ! n
                endif ! exit_bellman ! 8-13-2017 [b16]
                    
                call system_clock(tend)
                !if(printout6) write(unit=120,fmt='(a,i3,a,f12.4,a)') 'solve bellman period ',t,', time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'  
                if(printout6) write(*,fmt='(a,i3,a,f12.4,a)') 'solve bellman period ',t,', time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'  
                
                if(exit_bellman) exit ! exit the loop over t. 8-13-2017 [b17]
                
            enddo ! t
            deallocate( qray, xlist, bd ) ! 8-13-2017 [b18]
        endif ! exit_bellman 8-13-2017 [b19]
    end subroutine solve_bellman_1014   
    
    function f14(x) ! coarse value function for financial holdings. 09232016
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: f14, csp
        csp = inc - (x+hp) - HomTrnsCst(h,hp)
        if(csp<=0._wp.or.(x+hp)<0._wp)then ! 3.9.2017 Basically, the second condition is tested in subroutine "col_valid_fin." 3.11.2017 allow agents to leave zero bequest (replace x+hp<=0 with x+hp<0).
            f14 = penalty
            return
        endif        
        f14 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp)
    end function f14   
    
    function f1013(x) ! Refer to line 331. ! Use only linear interpolation. 072516. 072816.
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: f1013, csp, vp(2)
        integer :: n
        logical :: msg
        real(wp) :: vec(4)

        csp = inc - (x+hp) - HomTrnsCst(h,hp)   
        if(csp<=0._wp.or.(x+hp)<0._wp)then ! 3.11.2017 replace x+hp<=0 with x+hp<0.
            f1013 = penalty
            return
        endif
        if(opx==0)then
            vp(1) = binterpII(av,hv,v2d(:,:,1),x,hp,penalty,msg,vec)
            
            !! This block should be commented out because as long as the income is positive, we should let them interpolate between penalty and adjacent non-penalty.
            !if(msg==.false.)then
            !    f1013 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp) + survprob(t)*beta*vp(1)
            !else
            !    f1013 = penalty
            !endif
            
            f1013 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp) + survprob(t)*beta*vp(1) ! 3.11.2017 confirmed. bequest = x + hp
        else ! opx/=0             
            do n = 1, 2
                vp(n) = binterpII(av,hv,v2d(:,:,n),x,hp,penalty,msg,vec)
            enddo ! n
            !if(msg==.false..and.all(vp/=penalty))then
            !    f1013 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp) + survprob(t)*beta*dot_product(pzv,vp) ! dot_product(pzv,dum_wf)               
            !else
            !    f1013 = penalty
            !endif 
            f1013 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp) + survprob(t)*beta*dot_product(pzv,vp) ! dot_product(pzv,dum_wf)               
        endif
    end function f1013    
    
    subroutine col_valid_fin(lhpx,kpx,ypx,t,lowap,vec) ! t: current period 09232016. 10132016. 3.6.2017. 5.10.2017 used in solve_bellman_1014().
    ! lhpx: the level we conduct experiments playing with
        implicit none
        real(wp), intent(in) :: lhpx ! 3.6.2017 actually lhpx (its evoked name is "shv(m)") becomes useless in this subroutine after I replace lhpx (the trial next)) with h (the current) in the defition of "house".
        integer, intent(in) :: kpx, ypx, t
        real(wp) :: labinc, bizret, invscl, sswel, house
        real(wp) :: v1, v2, r2
        real(wp), intent(out) :: lowap
        real(wp), optional, intent(inout) :: vec(3)
        
        if(t<=8)then
            labinc = (1._wp-tauss/2._wp)*wage*efflab(t+1)*yv(ypx) ! next period's labor income 3.5.2017 added tauss/2.
            sswel  = 0._wp
        else
            labinc = 0._wp
            sswel  = benefit
        endif
        
        bizret = merge(0._wp, (1._wp-delzl(kpx))*kv(kpx), kpx==0) ! the worst business return
        invscl = merge(0._wp, kv(kpx), kpx==0) ! next period investment scale
        ! house  = (1._wp-lambda)*lhpx ! 3.5.2017 this line shows that "lhpx" should be the known current level of housing asset holding. <<<<============= not yet.
        house = (1._wp-lambda)*h ! 3.6.2017 replace lhpx with h (the CURRENT period's level of housing asset holdings).
        !v1  = invscl - (bizret+labinc+sswel)/(1._wp+rl) - house ! borrowing constraints (silos) 3.5.2017 change rd to rl as used in Meh. <<<<========
        v1  = invscl - (bizret)/(1._wp+rl) - house ! 3.4.2018
        !v2  = -(1._wp-deltah)*lhpx/(1._wp+rd) ! to rule out negative bequests, Yang 2009, pdf. 6. ! comment out 09282016
        v2  = -lhpx ! 09252016 Naka positive wealth. 10132016. 3.5.2017 agree that THROUGHOUT THE LIFECYCLE, total asset has to be nonnegative.
        lowap  = merge(v1, v2, v1>v2) ! lowap = max(v1,v2)
        lowap  = merge(av(1), lowap, av(1)>lowap) ! no less than the first grid point of financial asset holdings.
        
        ! debug zone
        if(present(vec))then
            vec(1) = v1
            vec(2) = v2
            vec(3) = lowap 
        endif
    end subroutine col_valid_fin
    
    subroutine csp_valid_fin(h,lhpx,inc,highap) ! 09232016. 3.6.2017 this should be the maximum of consumption rather than that of the financial asset holding.
    ! 3.6.2017 highap actually stands for the maximum consumption.
    ! lhpx: the trial next period's housing asset. lhpx=shv(m).
        implicit none
        real(wp), intent(in) :: h, lhpx, inc
        real(wp), intent(out) :: highap
        highap = inc - lhpx - HomTrnsCst(h,lhpx)
        highap = merge(av(adim),highap,highap>av(adim)) ! 3.6.2017 deal with available consumption exceeds the maximum finanacial asset hodling in the grid we specified.
    end subroutine csp_valid_fin
    
    subroutine solve_bellman_9(exit_bellman) ! 8-13-2017 [0]
        implicit none
        logical, intent(inout) :: exit_bellman ! 8-13-2017 [1]
        logical :: exit_brent                  ! 8-13-2017 [1]
        logical, dimension(:), allocatable :: qray ! 8-13-2017 [2]
        
        integer :: n, i, j, q, m, l
        integer :: tstart, tend, trate, tmax
        integer, dimension(:,:), allocatable :: xlist, bd  
        real(wp) :: intfund, capgain, labsup, labdem, ssbtax, bgp ! [1]
        integer :: opxi, loc1(1)
        
        real(wp) :: svcmp 
        integer  :: stp1, dsdim
        logical  :: log1 
        real(wp), dimension(:), allocatable :: shv, lowap, highap, smaxv, apvec, pyv, vpe1, vpe2     
        real(wp) :: best, bestl, bestu, bdis
        
        character(len=2) :: str1
        write(str1,fmt='(i2)') hdim
    
        ! main 
        call system_clock(tstart,trate,tmax)
        call system_clock(tstart)      
        
        allocate( xlist(adim*1*nmc*16,8), bd(4,2), vpe1(nmc), vpe2(nmc), pyv(nmc) ) ! 3.29.2017
        xlist = xl9
        bd = bd9
        t = 9
        
        allocate( qray(bd(4,2)) ) ! 8-13-2017 [3]
        qray = .false.            ! 8-13-2017
        
        do l = 1, 4
            !$omp parallel default(shared) private(n,capgain,intfund,i,stp1,log1,svcmp,shv,lowap,highap,smaxv,apvec,best,bdis,bestl,bestu,dsdim,m,loc1,opxi,labsup,labdem,ssbtax,bgp,exit_brent) ! [2] ! 8-13-2017 [4]
            allocate( v2d(adim,hdim,1:2) )
            !$omp do
            do q = bd(l,1), bd(l,2)
                do n = 1, hdim
                    
                    exit_brent = .false. ! 8-13-2017 [5]
                    
                    ax  = xlist(q,1)
                    hx  = n
                    kx  = xlist(q,3)
                    zx  = xlist(q,4)
                    yx  = xlist(q,5)
                    kpx = xlist(q,6)
                    ypx = xlist(q,7)
                    opx = xlist(q,8)                          
                    
                    a  = av(ax)
                    h  = hv(hx)
                    k  = kv(kx)
                    z  = merge(0._wp, merge(zlow,z2(kx),zx==0), kx==0)
                    y  = yv(yx)
                    kp = kv(kpx)
                    yp = yv(ypx)
                    dk = merge(0._wp, merge(delzl(kx),delzh(kx),zx==0), kx==0)  
                    
                    ! income (step 1)
                    labsup = efflab(t)*yv(yx)
                    if(kx==0)then
                        if(printout25)then

                            if(mode6taskid==0)then
                                ssbtax = tauss/2._wp*wage*labsup 
                                inc = merge(merge(transbeq,0._wp,printout17), &
                                            merge(transbeq,0._wp,printout17)+merge(rd*a,0._wp,a>0._wp), &
                                            tausv>0._wp) & ! 3-11-2018
                                      + wage*labsup - ssbtax !10.26.2017 a part of pre-tax income is not taxable by the current U.S. law
                                inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a, (1._wp+rd)*a, a>=0._wp), &
                                                  merge(a,(1._wp+rd)*a,a>0._wp), &
                                                  tausv>0._wp) &
                                      + (1._wp-deltah)*h - ttaxwok(inc)
                            elseif(mode6taskid==1)then ! wealth tax (not applied on period income but beginning capital, regardless of asset types).
                                ssbtax = tauss/2._wp*wage*labsup 
                                inc = merge(transbeq,0._wp,printout17) + wage*labsup - ssbtax !10.26.2017 a part of pre-tax income is not taxable by the current U.S. law                                   
                                inc = inc + (1._wp - deltah) * h + (1._wp + rd) * a
                                inc = (1._wp - taubal) * inc
                            elseif(mode6taskid==2)then
                                net_worth = merge(transbeq,0._wp,printout17) + (1._wp + rd) * a + (1._wp - deltah) * h + wage*labsup - ssbtax                                
                                ssbtax = tauss/2._wp*wage*labsup 
                                inc = merge(transbeq,0._wp,printout17) + merge(rd*a,0._wp,a>0._wp) + wage*labsup - ssbtax !10.26.2017 a part of pre-tax income is not taxable by the current U.S. law                                
                                if(net_worth > 5._wp*exempbar)then
                                    inc = inc + merge(a,(1._wp+rd)*a,a>0._wp) + (1._wp-deltah)*h - ttaxwok(inc) - 0.02_wp*(net_worth-5._wp*exempbar) - 0.01_wp*(4._wp*exempbar)       
                                elseif(net_worth>exempbar .and. net_worth<=5._wp*exempbar)then
                                    inc = inc + merge(a,(1._wp+rd)*a,a>0._wp) + (1._wp-deltah)*h - ttaxwok(inc) - 0.01_wp*(net_worth - exempbar)   
                                else
                                    inc = inc + merge(a,(1._wp+rd)*a,a>0._wp) + (1._wp-deltah)*h - ttaxwok(inc)
                                endif                                
                            else
                                print*, "error in taxation 9-5"
                            endif 
                        else
                            inc = wage*labsup
                            ssbtax = tauss/2._wp*wage*labsup 
                            inc = inc + merge( (1._wp+rd)*a, (1._wp+(1._wp-tausv)*rd)*a, a<0._wp) + (1._wp-deltah)*h - ttaxwok(inc) + merge(transbeq,0._wp,printout17) - ssbtax ! 3.26.2017 revision                            
                        endif !printout25
                        
                    else
                        intfund = merge(a-k, 0._wp, a>k) ! excess internal fund bears 
                        !inc = wage*labsup
                        bgp = c_grs_mat(ax, kx, zx, yx, t)
                        ssbtax = tauss/2._wp*wage*labsup 
                        
                        if(printout25)then
                        
                            !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       (1._wp+(1._wp-tausv)*rd)*k, intfund>0._wp), &
                            !                  merge((1._wp+(1._wp-tausv)*rd)*intfund,                      0._wp, intfund>0._wp), zx==0) &
                            !      + merge((1._wp-taubp)*bgp, bgp, zx==1) + (1._wp-deltah)*h - ttaxent(inc) + transbeq - ssbtax ! 3.26.2017                    

                            !inc = inc + bgp + merge(transbeq,0._wp,printout17) - ssbtax ! 4.10.2017
                            !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       (1._wp+(1._wp-tausv)*rd)*k, intfund>0._wp), & ! removed 9-12-2017
                            
                            !!!10.13.2017 comment out
                            !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       merge((1._wp+(1._wp-tausv)*rd)*k, (1._wp+(1._wp-tausv)*rd)*k+(1._wp+rd)*a ,a>0._wp), intfund>0._wp), & ! added 9-12-2017
                            !                  merge((1._wp+(1._wp-tausv)*rd)*intfund,                                                                               0._wp, intfund>0._wp), zx==0) &
                            !      + (1._wp-deltah)*h - ttaxent(inc) + merge(transbeq,0._wp,printout17) - ssbtax ! 4.10.2017
                            if(mode6taskid==0)then
                                inc = bgp + merge(merge(transbeq,0._wp,printout17), &
                                                  merge(transbeq,0._wp,printout17)+merge(rd*intfund,0._wp,intfund>0._wp), &
                                                  tausv>0._wp) + wage*labsup - ssbtax
                                inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*intfund, 0._wp, intfund>=0._wp), &
                                                  merge(intfund,0._wp,intfund>0._wp), &
                                                  tausv>0._wp) & ! change in net worth 3-11-2018
                                      + (1._wp - deltah) * h - ttaxent(inc)                                  
                                
                                !inc = wage*labsup + bgp + merge(transbeq,0._wp,printout17) - ssbtax ! 4.10.2017
                                !inc = inc + merge( (1._wp+(1._wp-tausv)*rd)*intfund, merge( 0._wp, (1._wp+rd)*a, a>0._wp), intfund>0._wp) & ! 10.13.2017
                                !      + (1._wp-deltah)*h - ttaxent(inc)                             
                            elseif(mode6taskid==1)then ! wealth tax (not applied on period income but beginning capital, regardless of asset types).
                                inc = bgp + merge(transbeq, 0._wp, printout17) + wage * labsup - ssbtax ! 4.10.2017
                                !inc = inc + (1._wp - taubal) * (merge(bgp, 0._wp, bgp > 0._wp) + merge((1._wp + (1._wp - tausv) * rd) * intfund, 0._wp, intfund>0._wp)) & ! 10.13.2017
                                !      + (1._wp - taubal) * (1._wp - deltah) * h - ttaxent(inc) + merge(bgp, 0._wp, bgp <= 0._wp)   
                                inc = inc + merge((1._wp + rd) * intfund, 0._wp, intfund>0._wp) + (1._wp - deltah) * h
                                inc = (1._wp - taubal) * inc
                            elseif(mode6taskid==2)then    
                                net_worth = bgp + merge(transbeq,0._wp,printout17) + merge((1._wp + rd) * intfund, 0._wp, intfund>0._wp) + (1._wp - deltah) * h + wage*labsup - ssbtax
                                inc = bgp + merge(transbeq,0._wp,printout17) + merge(rd*intfund,0._wp,intfund>0._wp) + wage*labsup - ssbtax
                                if(net_worth > 5._wp*exempbar)then
                                    inc = inc + merge(intfund,0._wp,intfund>0._wp) + (1._wp-deltah)*h - ttaxent(inc) - 0.02_wp*(net_worth-5._wp*exempbar) - 0.01_wp*(4._wp*exempbar)   
                                elseif(net_worth>exempbar .and. net_worth<=5._wp*exempbar)then
                                    inc = inc + merge(intfund,0._wp,intfund>0._wp) + (1._wp-deltah)*h - ttaxent(inc) - 0.01_wp*(net_worth - exempbar)   
                                else
                                    inc = inc + merge(intfund,0._wp,intfund>0._wp) + (1._wp-deltah)*h - ttaxent(inc)
                                endif 
                            else
                                print*, "error in taxation 9-6"    
                            endif !mode6taskid
                            
                        else
                            
                            !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       (1._wp+(1._wp-tausv)*rd)*k, intfund>0._wp), &
                            !                  merge((1._wp+(1._wp-tausv)*rd)*intfund,                      0._wp, intfund>0._wp), zx==0) &
                            !      + merge((1._wp-taubp)*bgp, bgp, zx==1) + (1._wp-deltah)*h - ttaxent(inc) + transbeq - ssbtax ! 3.26.2017                    

                            inc = inc + bgp ! 4.10.2017
                            !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       (1._wp+(1._wp-tausv)*rd)*k, intfund>0._wp), & ! removed 9-12-2017
                            
                            !!!10.13.2017 comment out
                            !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       merge((1._wp+(1._wp-tausv)*rd)*k, (1._wp+(1._wp-tausv)*rd)*k+(1._wp+rd)*a ,a>0._wp), intfund>0._wp), & ! added 9-12-2017
                            !                  merge((1._wp+(1._wp-tausv)*rd)*intfund,                                                                               0._wp, intfund>0._wp), zx==0) &
                            !      + (1._wp-deltah)*h - ttaxent(inc) + merge(transbeq,0._wp,printout17) - ssbtax ! 4.10.2017
                            
                            inc = inc + merge( (1._wp+(1._wp-tausv)*rd)*intfund, merge( 0._wp, (1._wp+rd)*a, a>0._wp), intfund>0._wp) & ! 10.13.2017
                                  + (1._wp-deltah)*h - ttaxent(inc) + merge(transbeq,0._wp,printout17) - ssbtax ! 4.10.2017
                        
                        endif
                        
                        !capgain = merge(rd*intfund, 0._wp, rd*intfund>0._wp)
                        ! inc = merge(wage*labsup,0._wp,zx==0) + merge(0._wp,capgain,tausvflag) + c_grs_mat(ax,kx,zx,yx,t) ! taxable income. Note: entrepreneurs always rent his labor efficency out.
                        !inc = merge(0._wp,capgain,tausvflag) + wage*labsup ! 3.12.2017 
                        
                        !if(zx==0)then
                        !    ssbtax = tauss*wage*labsup     
                        !else
                        !    labdem = c_lab_vec(kx) ! [4]
                        !    ssbtax = merge(tauss*wage*(labdem-labsup)+2._wp*tauss*wage*labsup,2*tauss*wage*labdem,labdem>labsup)
                        !endif
                        
                        !ssbtax = tauss/2._wp*wage*labsup ! 3.12.2017 correct.
                        !inc = inc + a + rd*merge(0._wp, a, a>=0._wp) + (1._wp-tausv)*merge(capgain,0._wp,tausvflag) + (1._wp-deltah)*h + (1._wp-taubp)*c_grs_mat(ax, kx, zx, yx, t) - ttaxent(inc) + transbeq - ssbtax                    
                        !inc = inc + a + (1._wp-tausv)*merge(capgain,0._wp,tausvflag) + (1._wp-deltah)*h + (1._wp-taubp)*c_grs_mat(ax, kx, zx, yx, t) - ttaxent(inc) + transbeq - ssbtax ! 3.12.2017 remove the duplicated borrowing cost
                        ! 3.13.2017 In this case, the agent engages in a business project anyway, so we don't need to consider the case where he borrows without running a business project.
                        ! 3.13.2017 In other words, he is either in the case that he has redundant capital to save in the bank or he needs to borrow to run his business (k>0) subject to higher than r_d interest rate.
                    endif   
                    
                    if(inc>0._wp)then
                        if(opx/=0)then
                            pzv(1) = pz2(kpx,1) ! probability of technology shock is determined by the size of next-period project (kpx).
                            pzv(2) = pz2(kpx,2) ! probability of technology shock is determined by the size of next-period project (kpx).    
                            do i = 1, 2
                                v2d(:,:,i) = wf(:,:,kpx,i-1,0,t+1) ! 3.31.2017 checked.    
                            enddo
                        else
                            v2d(:,:,1) = wf(:,:,kpx,0,0,t+1)    
                        endif ! opx
                        
                        stp1  = 1
                        log1  = .false.
                        svcmp = 100._wp ! used for convergence
                        allocate( shv(nsdim), lowap(nsdim), highap(nsdim), smaxv(nsdim), apvec(nsdim) ) ! first run
                        do while(stp1<iterasvmax .and. log1==.false.)

                            ! ACCURACY BLOCK --- START
                            if(stp1>1)then
                                best = shv(loc1(1))
                                !if(printout7) write(unit=109,fmt='(a,2i4,f12.8)'), ' new stp old idx ', stp1, loc1(1), best
                                if(loc1(1)/=1.and.loc1(1)/=dsdim)then
                                    bdis = abs(shv(loc1(1)+1)-shv(loc1(1)-1))
                                else
                                    if(loc1(1)==1)then
                                        bdis = abs(shv(1)-shv(2))       
                                    else ! loc1(1)==dsdim
                                        bdis = abs(shv(dsdim)-shv(dsdim-1))
                                    endif
                                endif
                                
                                !if(loc1(1)==1)then
                                !    bestl = best
                                !    bestu = best + 0.8_wp*bdis
                                !elseif(loc1(1)==dsdim)then
                                !    bestu = best
                                !    bestl = best - 0.8_wp*bdis
                                !else
                                !    bestl = best - 0.8_wp*bdis
                                !    bestu = best + 0.8_wp*bdis
                                !endif

                                if(loc1(1)==1)then
                                    !bestl = best
                                    bestl = best - 0.2_wp*bdis ! 10012016
                                    bestu = best + 0.8_wp*bdis
                                    if(bestl<hv(1)) bestl = hv(1) ! 10012016
                                elseif(loc1(1)==dsdim)then
                                    !bestu = best
                                    bestu = best + 0.2_wp*bdis ! 10012016
                                    bestl = best - 0.8_wp*bdis
                                    if(bestu>hv(hdim)) bestu = hv(hdim) ! 10012016
                                else
                                    bestl = best - 0.8_wp*bdis
                                    bestu = best + 0.8_wp*bdis
                                    if(bestl<hv(1)) bestl = hv(1) ! 10012016
                                    if(bestu>hv(hdim)) bestu = hv(hdim) ! 10012016
                                endif                                
                                
                                dsdim = sdim ! zoom-in search area (using less grid points by default).
                                deallocate( shv, lowap, highap, smaxv, apvec )
                                allocate( shv(sdim), lowap(sdim), highap(sdim), smaxv(sdim), apvec(sdim) ) ! small region brent's method --- part II
                                !if(printout7) write(unit=109,fmt='(a,3f8.4)'), ' boundary ', bestl, bestu, bdis
                                !call grid(shv,bestl,bestu,1._wp) ! good 3.13.2017 keep using 1._wp with the same reason stated in periods 10-14.
                                call grid(shv,bestl,bestu,gsdim)
                            else ! stp1 = 1
                                dsdim = nsdim ! bigger grid              
                                if(hx==1)then
                                    call grid(shv,hv(1),hv(hdim),gsdim)    
                                else ! hx
                                    !if(cwh(ax,hx-1,kx,zx,yx,kpx,ypx,opx,t)/=penalty)then ! comment out 10092016
                                    !    call grid(shv,cwh(ax,hx-1,kx,zx,yx,kpx,ypx,opx,t),hv(hdim),gsdim)
                                    !else
                                        call grid(shv,hv(1),hv(hdim),gsdim) ! the default full range
                                    !endif ! cwh                                 
                                endif ! hx
                            endif ! stp
                            ! ACCURACY BLOCK --- END                             
                            
                            !if(stp1>1)then
                            !    deallocate( shv, lowap, highap, smaxv, apvec )
                            !    allocate( shv(sdim), lowap(sdim), highap(sdim), smaxv(sdim), apvec(sdim) ) ! small region brent's method --- part II
                            !endif
                            !! smart move 09212016 important for smoothness.
                            !if(hx>1)then
                            !    if(cwh(ax,hx-1,kx,zx,yx,kpx,ypx,opx,t)/=penalty)then
                            !        call grid(shv,cwh(ax,hx-1,kx,zx,yx,kpx,ypx,opx,t),hv(hdim),gsdim)
                            !    else
                            !        call grid(shv,hv(1),hv(hdim),gsdim) ! the default full range
                            !    endif
                            !else
                            !    call grid(shv,hv(1),hv(hdim),gsdim)
                            !endif                                
                            
                            smaxv = penalty
                            lowap  = -1
                            highap = -1

                            !if(printout7) write(unit=109,fmt='(a,9a)'), ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                            !if(printout7) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,kpx,ypx,opx,t
                            
                            !debcnt1 = 0
                            do m = 1, dsdim
                                hp = shv(m) ! pass into the thread
                                call col_valid_fin(shv(m),kpx,ypx,t,lowap(m))
                                call csp_valid_fin(h,shv(m),inc,highap(m))      
                                !if(printout7) write(unit=109,fmt='(a,2f8.4,i5)'), ' -- ', lowap(m), highap(m), merge(0,1,lowap(m)<highap(m))
                                if(lowap(m)>=highap(m).or.highap(m)<0._wp)then
                                    !if(printout7) write(unit=109,fmt='(i3)'), debcnt1
                                    cycle ! negative consumption in one of the refined housing assets level 
                                endif
                                call brent_localizer(f9,lowap(m),highap(m),apvec(m),smaxv(m),exit_brent) ! <-------------------------------------------- ! 8-13-2017 [6]
                                !if(printout7) write(unit=109,fmt='(a,2f8.4)'), '    ', apvec(m), smaxv(m)
                                !if(printout7.and.m==sdim) write(unit=109,fmt='(a)'), '    '
                                
                                if(exit_brent) exit ! exit this loop over m. ! 8-13-2017 [7]
                                
                            enddo
                            
                            if(exit_brent==.false.)then ! 8-13-2017 [8]
                                ! get the maximum from the set of local bests.
                                if(maxval(smaxv)/=penalty)then
                                    loc1 = maxloc(smaxv)
                                    stp1 = stp1 + 1
                                    if(abs(svcmp-smaxv(loc1(1)))<err_svdist) exit
                                    svcmp = smaxv(loc1(1))                            
                                else
                                    log1 = .true.
                                    svcmp= penalty
                                    !if(printout7.and.debcnt1==sdim) write(unit=109,fmt='(a)'), ' all penalty  ' ! correct 09222016
                                    !if(printout7.and.debcnt1/=sdim) write(unit=109,fmt='(a,10i3)'), ' something wrong in smaxv ', ax,hx,kx,zx,yx,kpx,ypx,opx,t, debcnt1 ! no leaks 09222016
                                endif
                            else                                                           ! 
                                exit ! exit this housing asset holding grid searching zone ! 8-13-2017 [9]
                            endif ! exit_brent                                             !
                            
                        enddo ! do while                              
                        
                        if(exit_brent==.false.)then ! 8-13-2017 [10]
                            ! save searching results or give a warning mark.
                            if(svcmp/=penalty)then ! SOME CASES NEED TO DO COMPARISON 
                                cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = smaxv(loc1(1))
                                cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = apvec(loc1(1))
                                cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = shv(loc1(1))    
                                cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = -99
                                if(opx==0)then
                                    cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 0
                                    cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = inc - (apvec(loc1(1))+shv(loc1(1))) - HomTrnsCst(h,shv(loc1(1)))                                
                                else ! opx/=0
                                    cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = kpx
                                    cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = inc - (apvec(loc1(1))+shv(loc1(1))) - HomTrnsCst(h,shv(loc1(1)))                                
                                    if(opx==1)then ! just choose to be a worker instead
                                        if( cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) < cwf(ax,hx,kx,zx,yx,0,0,0,t))then ! 3.12.2017 Entrepreneurs > Workers > Retirees. 3.31.2017 correct.
                                            cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwf(ax,hx,kx,zx,yx,0,0,0,t) 
                                            cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwa(ax,hx,kx,zx,yx,0,0,0,t)
                                            cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwh(ax,hx,kx,zx,yx,0,0,0,t)
                                            cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwk(ax,hx,kx,zx,yx,0,0,0,t)
                                            cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cww(ax,hx,kx,zx,yx,0,0,0,t)
                                            cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwc(ax,hx,kx,zx,yx,0,0,0,t)                                        
                                        endif
                                    else ! opx==2 ! choose to engage in a smaller investment project.
                                        !opxi = merge(opx-1,opx-2,kx/=0)
                                        ! 3.30.2017 only need to pay attention to kp and op.
                                        opxi = merge(opx-1,opx-2,zx==1) ! 3.30.2017
                                        if( cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) < cwf(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t))then
                                            cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwf(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t) 
                                            cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwa(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)
                                            cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwh(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)
                                            cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwk(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)                                            
                                            cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cww(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t) 
                                            cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwc(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)                                         
                                        endif
                                    endif ! opx==1
                                endif ! opx/=0
                            else ! infeasible point due to collateral and budget constraints. DOWNGRADE WITHOUT ANY COMPARISON
                                if(opx==0)then ! no where to downgrade their expenses or costs
                                    cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 1 
                                elseif(opx==1)then ! may instead choose to be a worker
                                    cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwf(ax,hx,kx,zx,yx,0,0,0,t) 
                                    cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwa(ax,hx,kx,zx,yx,0,0,0,t)
                                    cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwh(ax,hx,kx,zx,yx,0,0,0,t)
                                    cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwk(ax,hx,kx,zx,yx,0,0,0,t)                                    
                                    cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cww(ax,hx,kx,zx,yx,0,0,0,t)   
                                    cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwc(ax,hx,kx,zx,yx,0,0,0,t)                                
                                elseif(opx==2)then ! may instead choose to engage in an investment project of smaller scale.
                                    opxi = merge(opx-1,opx-2,zx==1) ! 3.31.2017 
                                    cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwf(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t) 
                                    cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwa(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)
                                    cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwh(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)
                                    cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwk(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)                                                                                
                                    cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cww(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)   
                                    cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwc(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)                                
                                endif    
                            endif 
                        endif ! exit_brent 8-13-2017 [11]
                        
                        deallocate( shv, lowap, highap, smaxv, apvec )     
                        
                    else
                        cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 1                           
                    endif ! inc>0
                    
                    if(exit_brent) exit ! exit this n loop (over the indices of housing asset holding grid).  8-13-2017 [12]
                    
                enddo ! n
                
                if(exit_brent) qray(q)=.true. ! 8-13-2017 [13]
                
            enddo ! q
            !$omp end do
            deallocate( v2d )
            !$omp end parallel
            
            if(any(qray==.true.)) exit_bellman=.true. ! 8-13-2017 [14]
            if(exit_bellman) exit ! exit the loop over l ! 8-13-2017 [15]
        enddo ! l
        
        if(exit_bellman==.false.)then ! 8-13-2017 [16]
            ! Review done! 3.13.2017 afternoon 4:25 pm
            if(printout2)then            
                do l = 1, 4                
                    do q = bd(l,1), bd(l,2) ! income in period 14 could be used later.
                        ax  = xlist(q,1)
                        hx  = xlist(q,2)
                        kx  = xlist(q,3)
                        zx  = xlist(q,4)
                        yx  = xlist(q,5)
                        kpx = xlist(q,6)
                        ypx = xlist(q,7)
                        opx = xlist(q,8)  
                                 
                        if(ax==1.and.hx==1)then
                            write(unit=056,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                            write(unit=056,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                            write(unit=056,fmt='(3x,9x,'//str1//'i21)') (j,j=1,hdim) 
                            write(unit=056,fmt='(3x,9x,'//str1//'f21.4)') (hv(j),j=1,hdim) 
                            do i = 1, adim
                                write(unit=056,fmt='(i3,x,f8.4,'//str1//'(x,f20.8))') i,av(i),(cwf(i,j,kx,zx,yx,kpx,ypx,opx,t),j=1,hdim)    
                            enddo
                        endif 
                        
                    enddo ! q            
                enddo ! l
            endif ! printout2  
            
            ! 3.14.2017 Whether we have a valide beginning-of-period combination depends on whether the corresponding end-of-period combinations are both valid.
            ! 3.14.2017 The uncertainty comes solely from the innovation of business ideas.<-------
            ! 3.31.2017 Stop here 12:23 pm. revise 2<=n.and.n<=4
            do n = 1, 7
                kx = tvector(n,1)
                zx = tvector(n,2)
                do yx = 1, nmc
                    if(n==1)then ! 3.14.2017 Workers have the chance to conceive a business idea. 3.31.2017 Taking expectation.
                        where(cww(:,:,kx,zx,yx,0,0,0,t)==-99.and.cww(:,:,kx,zx,yx,1,0,2,t)==-99) & ! workers might come across an entrepreneurial idea. 3.14.2017
                            & wf(:,:,kx,zx,yx,t) = (1._wp-probtk(kx))*cwf(:,:,kx,zx,yx,0,0,0,t) + probtk(kx)*cwf(:,:,kx,zx,yx,1,0,2,t)    
                    elseif(2<=n.and.n<=4)then ! 3.14.2017 These are the cases where entrepreneurs receive bad busniess shock in the beginning of period, which forces them to exit entrepreneurship.
                        where(cww(:,:,kx,zx,yx,0,0,0,t)==-99.and.cww(:,:,kx,zx,yx,1,0,2,t)==-99) & ! 3.31.2017 entrepreneurs receive bad shocks and are "doomed" to exit entrepreneurship. 3.14.2017
                            & wf(:,:,kx,zx,yx,t) = (1._wp-probtk(0))*cwf(:,:,kx,zx,yx,0,0,0,t) + probtk(0)*cwf(:,:,kx,zx,yx,1,0,2,t) ! 3.31.2017 correct the error that probtk's index is not 0.
                            !& wf(:,:,kx,zx,yx,t) = 1._wp ! (1._wp-probtk(0)) ! *cwf(:,:,kx,zx,yx,0,0,0,t) + probtk(0)*cwf(:,:,kx,zx,yx,1,0,2,t) ! 3.31.2017 correct the error that probtk's index is not 0.
                    elseif(5<=n.and.n<=6)then
                        where(cww(:,:,kx,zx,yx,kx,0,1,t)==-99.and.cww(:,:,kx,zx,yx,kx+1,0,2,t)==-99) &
                            & wf(:,:,kx,zx,yx,t) = (1._wp-probtk(kx))*cwf(:,:,kx,zx,yx,kx,0,1,t) + probtk(kx)*cwf(:,:,kx,zx,yx,kx+1,0,2,t)
                    elseif(n==7)then
                        where(cww(:,:,kx,zx,yx,kx,0,1,t)==-99) &
                            & wf(:,:,kx,zx,yx,t) = cwf(:,:,kx,zx,yx,kx,0,1,t) ! 3.14.2017 The entrepreneur who runs the largest project has no more room to access more to advanced bussiness project.
                    endif
                    
                    if(printout2)then
                        write(unit=026,fmt='(a,3(i3),a,i3)') '(kx,zx,yx) ',kx,zx,yx, ' t ',t     
                        write(unit=026,fmt='(3x,8x,'//str1//'(8x,i8))') (j,j=1,hdim)             
                        write(unit=026,fmt='(3x,5x,"999",'//str1//'(8x,f8.4))') (hv(j),j=1,hdim)             
                        do j = 1, adim
                            write(unit=026,fmt='(i3,x,f7.2,'//str1//'(x,f15.5))') j, av(j),(wf(j,i,kx,zx,yx,t),i=1,hdim)
                        enddo                                  
                    endif ! printout2                 
                    
                enddo ! yx          
            enddo ! n 
        endif ! exit_bellman 8-13-2017 [17]
        
        call system_clock(tend)
        !if(printout6) write(unit=120,fmt='(a,i3,a,f12.4,a)') 'solve bellman period ',t,', time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'          
        if(printout6) write(*,fmt='(a,i3,a,f12.4,a)') 'solve bellman period ',t,', time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'          
        deallocate( qray, xlist, bd ) ! [18]
 
    end subroutine solve_bellman_9
    
    function f9(x)
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: f9, csp, vp(2), vec(4)
        integer :: n
        logical :: msg

        csp = inc - (x+hp) - HomTrnsCst(h,hp)
        
        if(csp<=0._wp .or. (x+hp)<0._wp)then ! 3.13.2017 replace x+hp<=0 with x+hp<0
            f9 = penalty
            return
        endif        
        
        if(opx==0)then
            vp(1) = binterpII(av,hv,v2d(:,:,1),x,hp,penalty,msg,vec)
            f9 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp) + survprob(t)*beta*vp(1)
            !if(msg==.false.)then
            !    f9 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp) + survprob(t)*beta*vp(1)
            !else
            !    f9 = penalty
            !endif            
        else ! opx/=0
            do n = 1, 2
                vp(n) = binterpII(av,hv,v2d(:,:,n),x,hp,penalty,msg,vec)    
            enddo ! n
            f9 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp) + survprob(t)*beta*dot_product(pzv,vp) ! dot_product(pzv,dum_wf)               
            !if(msg==.false..and.all(vp/=penalty))then
            !    f9 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp) + survprob(t)*beta*dot_product(pzv,vp) ! dot_product(pzv,dum_wf)               
            !else
            !    f9 = penalty
            !endif             
        endif
    end function f9
    
    !subroutine coarse_SBE_profit_matrices(c_prf_mat,c_grs_mat)
    ! 3.4.2017 checked. 
    ! output: 
    !   c_grs_mat: before-tax business profit
    !   c_lab_vec: labor efficiency demand
    !   c_opt_vec: production
    ! 3.26.2017 revised. let this subroutine return negative profit (i.e., loss) so that the borrowing cost can be easily used in the disposable income computation.
    subroutine coarse_SBE_profit_matrices(c_grs_mat, c_lab_vec, c_opt_vec) 
        implicit none
        real(wp), intent(out), dimension(adim,0:(kdim-1),0:1,0:nmc,1:14) :: c_grs_mat ! c_prf_mat, 
        real(wp), dimension(:), allocatable :: nv
        real(wp), dimension(:), intent(out) :: c_lab_vec, c_opt_vec
        integer :: n, t
        real(wp) :: product, labcost, bizincm, extfund, tblincm, r, dk, zv
        character(len=5) :: str1, str2
        write(str1,fmt='(i5)') ndim
        write(str2,fmt='(i5)') adim        
        allocate( nv(ndim) ) ! , c_lab_vec(kdim-1), c_opt_vec(kdim-1) )
        !c_prf_mat = 0._wp
        c_grs_mat = 0._wp ! 9-14-2017
        c_lab_vec = 0._wp ! 9-14-2017
        c_opt_vec = 0._wp ! 9-14-2017
        
        ! 3.26.2017 production and labor employment obtained here should be only used when business shock is postiive.
        if(printout28)then
            ! Do nothing
        else
            do n = 1, kdim-1
                ! c_lab_vec(n) = (wage/(z2(n)*nu*(1._wp-alpha))*kv(n)**(-alpha*nu))**(1._wp/((1._wp-alpha)*nu-1._wp))            
                c_lab_vec(n) = ((1._wp+tauss/2._wp)*wage/(z2(n)*nu*(1._wp-alpha))*kv(n)**(-alpha*nu))**(1._wp/((1._wp-alpha)*nu-1._wp)) ! condiser social security tax 3.3.2017             
                c_opt_vec(n) = z2(n)*(kv(n)**(alpha)*c_lab_vec(n)**(1._wp-alpha))**nu ! production. should only be used when business shock is positive. 
            enddo
        endif ! printout28
        
        do t = 1, 14
            do yx = 0, nmc ! the labor efficiency of the entrepreneur
                if(t>=10)then ! old
                    if(yx/=0) cycle ! Entrepreneurs of retirement age have NO labor efficiency.   
                else ! t<10 young
                    if(yx==0) cycle ! Young entrepreneurs have non-zero labor efficiency shock.
                endif
                do zx = 1, 0, -1 ! business shock
                    do kx = 1, (kdim-1) ! so we ignore cases that kx==0. 09282016
                        if(zx==1.and.kx==0) cycle ! Stupid. There is no such case, because kx-loop starts with 1. 3.3.2017
                        do ax = 1, adim
                            
                            !r  = merge( rd+(kv(kx)-av(ax))/kv(kx)*gamma5, rd, kv(kx)>av(ax)) ! individual borrowing costs or saving interset rate.
                            !r  = merge( merge(rd+(kv(kx)-av(ax))/kv(kx)*gamma5, rd+gamma5, av(ax)>0._wp), rd, kv(kx)>av(ax)) ! 10.13.2017
                            r  = merge(rd+(kv(kx)-av(ax))/kv(kx)*gamma5, rd, kv(kx)>av(ax)) ! 3-11-2018
                            dk = merge( delzl(kx), delzh(kx), zx==0) ! depreciation
                            zv = merge( zlow, z2(kx), zx==0) ! level of business shock                            
                            
                            if(zx==0)then ! entrepreneurs face bad business shock
                                
                                product = 0._wp ! no production, for all business projects.   
                                labcost = 0._wp ! 10.13.2017 correct. Refer to pg. 696, Meh 2005, the firm's optimization problem.                        
                                
                            else ! zx/=0 entrepreneurs face good business shock.
                                product = c_opt_vec(kx) ! zv*(kv(kx)**(alpha)*c_lab_vec(kx)**(1._wp-alpha))**nu
                                labcost = (1._wp+tauss/2._wp)*wage*c_lab_vec(kx)
                                
                            endif ! zx
                                
                            ! it covers the case where net financial position is negative. 3.4.2017
                            !extfund = merge(kv(kx)-av(ax), 0._wp, kv(kx)-av(ax)>0._wp) ! 3.25.2017 irrelevant to the level of zx.
                            !extfund = merge( merge(kv(kx)-av(ax), kv(kx), av(ax)>0._wp), 0._wp, kv(kx)-av(ax)>0._wp) ! 3-2-2018
                            extfund = merge(kv(kx)-av(ax), 0._wp, kv(kx)-av(ax)>0._wp) ! 3-11-2018
                            
                            ! Checked with Kitao, 2008, pdf. 5; Meh, 2005, pdf. 9. Kitao considers capital interest into the profit; Meh doesn't. That's why the discrepancy arises.
                            ! Anyway, I don't add the capital gain (saving interest) obtained from own funds. <--------- Important!! 09282016
                            ! Details can be found in housing_v9.pdf 3.3.2017
                            
                            !c_grs_mat(ax,kx,zx,yx,t) = product - dk*kv(kx) - r*extfund - labcost ! "taxable" business income !<------------------------------- (1) 09252016 correct!
                            !c_grs_mat(ax,kx,zx,yx,t) = product + (1._wp-dk)*kv(kx) - (1._wp+r)*extfund - labcost ! 3.4.2017, 3.25.2017 when hit by a negative business shock, this business porfit becomases the tota capital costs the entrepreneur needs to pay.
                            
                            if(zx==1)then
                                c_grs_mat(ax,kx,zx,yx,t) = product + (1._wp-dk)*kv(kx) - (1._wp+r)*extfund - labcost ! 3.26.2017 probably taking negative values reflect the loss due to borrowing cost and bad business shock
                            else ! zx==0 no production, no labcost, temparary assumption: business capital can be idle. 
                                !c_grs_mat(ax,kx,zx,yx,t) = product + (1._wp-dk)*kv(kx) - (1._wp+r)*extfund - labcost
                                c_grs_mat(ax,kx,zx,yx,t) = (1._wp-dk)*kv(kx) - (1._wp+r)*extfund ! 10.13.2017 remnant capital after depreciation and principle (as well as interest) needed to be paid back to the financial intermediary.
                            endif ! zx

                            bizincm = c_grs_mat(ax,kx,zx,yx,t) ! used in the print-out block
                            
                            !c_prf_mat(ax,kx,zx,yx,t) = product + (1._wp-dk)*kv(kx) - (1._wp+r)*(kv(kx)-av(ax)) - labcost - ttaxent(bizincm)
                            ! Excess own fund accruses saving interests, not counted as business income
                            ! Just used for business decision making, not the real entrepreneurial disposable income for household decision making.
                            
                            !c_prf_mat(ax,kx,zx,yx,t) = product + (1._wp-dk)*kv(kx) - (1._wp+r)*extfund - labcost - ttaxent(c_grs_mat(ax,kx,zx,yx,t)) ! uselss but needs to keep it
                            
                        enddo ! ax
                    enddo ! kx
                enddo ! zx
            enddo ! yx
        enddo ! t
        
        deallocate( nv ) ! c_lab_vec, c_opt_vec ) ! 3.5.2017 obsolete.
    end subroutine coarse_SBE_profit_matrices
    
    subroutine asset_array_for_brentmaxrhs(as,ae,bs,be,rs,re,smeg) ! 09102016 ok
    ! a* are end points of the array of domain of last period value function (those used for expectation)
    ! b* are end points of the array of domain of coarse value function just obtained from the proceeding block
        implicit none
        integer, intent(in) :: as, ae, bs, be ! the end indices of two arrays. 
        integer, intent(out) :: rs, re
        character(len=*), intent(out) :: smeg

        if((as<=bs.and.be<=ae).or.(bs<=as.and.ae<=be))then ! one is enclosed by the other.
            if(as<=bs.and.be<=ae)then
                rs = bs
                re = be
            else
                rs = as
                re = ae
            endif
            if(rs/=re)then
                smeg = 'normal'        
            else
                smeg = 'single'
            endif
        elseif(bs<as.or.ae<be)then
            if(bs<as)then
                rs = as
                re = be
            else
                rs = bs
                re = ae
            endif
            if(rs<re)then
                smeg = 'normal'        
            elseif(rs==re)then
                smeg = 'single'
            else
                smeg = 'abnormal'
            endif
        endif 
        
        !else
        !    smeg = 'abnormal'     
        !endif          
    end subroutine asset_array_for_brentmaxrhs
    
    subroutine multiple_arrays_for_brentmaxrhs(as,ae,bs,be,rs,re,smeg) ! 072516 09102016
        implicit none
        integer, intent(in) :: bs, be
        integer, dimension(:), intent(in) :: as, ae
        integer, intent(out) :: rs, re
        integer :: aas, aae
        character(len=*), intent(out) :: smeg
        aas = maxval(as)
        aae = minval(ae)
        call asset_array_for_brentmaxrhs(aas,aae,bs,be,rs,re,smeg)   
    end subroutine multiple_arrays_for_brentmaxrhs
    
    subroutine collbd_for_brentmaxrhs(apx,afx,alx,cafx,calx,newstr,newend,smeg,collbd,answer) ! 072516
        implicit none
        integer, intent(in) :: afx(2), alx(2), cafx, calx, newstr, newend,collbd
        integer, intent(inout) :: apx
        integer, intent(out) :: answer(2)
        character(len=*), intent(inout) :: smeg  
        if(smeg=='normal')then
            if(collbd>=newstr.and.collbd<newend)then
                answer(1)=collbd
                answer(2)=newend
            elseif(collbd==newend)then
                answer(1)=collbd
                answer(2)=collbd
                smeg='single'      
            elseif(collbd<newstr)then
                answer(1)=newstr
                answer(2)=newend                
            else
                answer(1)=-2
                answer(2)=-2
                smeg='abnormal'
            endif
        elseif(smeg=='single')then
            if(newstr<collbd)then
                answer(1)=-3
                answer(2)=-3
                smeg='abnormal'
            else
                answer(1)=newstr
                answer(2)=newend
            endif
        else ! smeg='abnormal'
            answer = -1
        endif
        
        if(smeg/='abnormal')then ! 09112016 added
            if(apx<answer(1))then
                apx = answer(1)    
            endif
        endif
        
        ! if(apx<collbd.or.apx<newstr)then smeg='abnormal' ! 09112016 removed
        
    end subroutine collbd_for_brentmaxrhs
    
    subroutine solve_bellman_18(exit_bellman) ! 8-13-2017 [0]
        implicit none
        logical, intent(inout) :: exit_bellman     ! [1]
        logical :: exit_brent                      !
        logical, dimension(:), allocatable :: qray ! [2]
        
        integer :: tstart, tend, trate, tmax
        integer :: n, i, j, m, q, l
        integer, dimension(:,:), allocatable :: xlist, bd
        real(wp) :: capgain, intfund, labsup, labdem, ssbtax, bgp
        
        real(wp) :: svcmp 
        integer  :: stp1, loc1(1)
        logical  :: log1 
        integer :: opxi
        !real(wp), dimension(:,:), allocatable :: pbym ! probability transition matrix of business-and-labor pairs
        real(wp), dimension(:), allocatable :: shv, lowap, highap, smaxv, apvec
        real(wp) :: pyv(nmc), vpe1(nmc), vpe2(nmc)  
        real(wp) :: best, bdis, bestl, bestu, dsdim
        character(len=2) :: str1
        write(str1,fmt='(i2)') hdim
        ! debug delcaration zone
        
        call system_clock(tstart,trate,tmax)

        allocate( xlist(adim*1*(nmc)**2*16,8), bd(4,2) )  ! 3.30.2017
        xlist = xl18
        bd = bd18
        
        allocate( qray(bd(4,2)) ) ! [3-1]
        
        do t = 8, 1, -1
            qray = .false.            ! [3-2]
            call system_clock(tstart)
            do l = 1, 4
                !$omp parallel default(shared) private(n,capgain,intfund,i,stp1,log1,svcmp,shv,lowap,highap,smaxv,apvec,best,bdis,bestl,bestu,dsdim,m,loc1,opxi,labsup,labdem,ssbtax,bgp,exit_brent) ! [4]
                allocate( v2d(adim,hdim,1:2) )
                !$omp do
                do q = bd(l,1), bd(l,2)
                    do n = 1, hdim
                        
                        exit_brent = .false. ! [5]
                        
                        ax  = xlist(q,1)
                        hx  = n
                        kx  = xlist(q,3)
                        zx  = xlist(q,4)
                        yx  = xlist(q,5)
                        kpx = xlist(q,6)
                        ypx = xlist(q,7)
                        opx = xlist(q,8) ! actually, it corresponds to the level of kpx.
                        
                        a = av(ax)
                        h = hv(hx)
                        k = kv(kx)
                        z = merge(0._wp, merge(zlow,z2(kx),zx==0), kx==0)
                        y = yv(yx)
                        kp = kv(kpx)
                        yp = yv(ypx)
                        dk = merge(0._wp, merge(delzl(kx),delzh(kx),zx==0), kx==0)   
                        
                        labsup = efflab(t)*yv(yx)
                        if(kx==0)then
                            if(printout25)then

                                if(mode6taskid==0)then
                                    ssbtax = tauss/2._wp*wage*labsup
                                    inc = merge(merge(0._wp, transbeq, printout17==.False. .and. t/=1), &
                                                merge(0._wp, transbeq, printout17==.False. .and. t/=1)+merge(rd*a,0._wp,a>0._wp), &
                                                tausv>0._wp) &
                                          + wage*labsup - ssbtax ! 10.26.2017 revision                                     
                                    inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a, (1._wp+rd)*a, a>=0._wp), &
                                                      merge(a,(1._wp+rd)*a,a>0._wp), &
                                                      tausv>0._wp) & 
                                          + (1._wp-deltah)*h - ttaxwok(inc)
                                elseif(mode6taskid==1)then ! wealth tax (not applied on period income but beginning capital, regardless of asset types).
                                    ssbtax = tauss/2._wp*wage*labsup
                                    inc = merge(0._wp, transbeq, printout17==.False. .and. t/=1) + wage*labsup - ssbtax ! 10.26.2017 revision                                     
                                    !inc = inc + (1._wp - taubal) * merge( 0._wp, (1._wp + (1._wp - tausv) * rd) * a, a < 0._wp) &
                                    !    + (1._wp - taubal) * (1._wp - deltah) * h - ttaxwok(inc) + merge(0._wp, (1._wp + (1._wp - tausv) * rd) * a, a > 0._wp)
                                    inc = inc + (1._wp - deltah) * h + (1._wp + rd) * a
                                    inc = (1._wp - taubal) * inc
                                elseif(mode6taskid==2)then
                                    ssbtax = tauss/2._wp*wage*labsup
                                    net_worth = merge(0._wp, transbeq, printout17==.False. .and. t/=1) + (1._wp + rd) * a + (1._wp - deltah) * h + wage*labsup - ssbtax
                                    inc = merge(0._wp, transbeq, printout17==.False. .and. t/=1) + merge(rd*a,0._wp,a>0._wp) + wage*labsup - ssbtax                                   
                                    if(net_worth > 5._wp*exempbar)then
                                        inc = inc + merge(a,(1._wp+rd)*a,a>0._wp) + (1._wp-deltah)*h - ttaxwok(inc) - 0.02_wp*(net_worth-5._wp*exempbar) - 0.01_wp*(4._wp*exempbar)       
                                    elseif(net_worth>exempbar .and. net_worth<=5._wp*exempbar)then
                                        inc = inc + merge(a,(1._wp+rd)*a,a>0._wp) + (1._wp-deltah)*h - ttaxwok(inc) - 0.01_wp*(net_worth - exempbar)
                                    else    
                                        inc = inc + merge(a,(1._wp+rd)*a,a>0._wp) + (1._wp-deltah)*h - ttaxwok(inc)
                                    endif
                                else
                                    print*, "error in taxation 18-7"
                                endif
                            else
                                inc = wage*labsup
                                ssbtax = tauss/2._wp*wage*labsup
                                inc = inc + merge( (1._wp+rd)*a, (1._wp+(1._wp-tausv)*rd)*a, a<0._wp) + (1._wp-deltah)*h - ttaxwok(inc) + merge(0._wp, transbeq, printout17==.False. .and. t/=1) - ssbtax ! 3.26.2017 revision
                            endif
                                
                            !capgain = merge(rd*a,0._wp,rd*a>0._wp)
                            !inc = merge(0._wp,capgain,tausvflag) + wage*labsup
                            !inc = inc + a + rd*merge(0._wp, a, a>=0._wp) + (1._wp-tausv)*merge(capgain,0._wp,tausvflag) + (1._wp-deltah)*h - ttaxwok(inc) + transbeq - ssbtax 
                        else ! entrepreneurs
                            intfund = merge(a-k,0._wp,a>k) ! excess internal fund bears 
                            !inc = wage*labsup
                            bgp = c_grs_mat(ax, kx, zx, yx, t)
                            ssbtax = tauss/2._wp*wage*labsup ! 3.14.2017 correct.  
                            if(printout25)then
                                !inc = inc + bgp + merge(0._wp, transbeq, printout17==.False. .and. t/=1) - ssbtax ! 4.10.2017 
                                !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       (1._wp+(1._wp-tausv)*rd)*k, intfund>0._wp), & ! removed 9-12-2017
                                
                                !!!10.13.2017
                                !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       merge((1._wp+(1._wp-tausv)*rd)*k, (1._wp+(1._wp-tausv)*rd)*k+(1._wp+rd)*a ,a>0._wp), intfund>0._wp), & ! added 9-12-2017
                                !                  merge((1._wp+(1._wp-tausv)*rd)*intfund,                                                                               0._wp, intfund>0._wp), zx==0) &
                                !      + (1._wp-deltah)*h - ttaxent(inc) + merge(0._wp, transbeq, printout17==.False. .and. t/=1) - ssbtax ! 4.10.2017                            
                                if(mode6taskid==0)then
                                    inc = bgp + merge(merge(transbeq,0._wp,printout17), &
                                                      merge(transbeq,0._wp,printout17)+merge(rd*intfund,0._wp,intfund>0._wp), &
                                                      tausv>0._wp) + wage*labsup - ssbtax
                                    inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*intfund, 0._wp, intfund>=0._wp), &
                                                      merge(intfund,0._wp,intfund>0._wp), &
                                                      tausv>0._wp) & ! change in net worth 3-11-2018
                                          + (1._wp - deltah) * h - ttaxent(inc)                                      
                                    
                                    !inc = bgp + merge(0._wp, transbeq, printout17==.False. .and. t/=1) + wage*labsup - ssbtax ! 4.10.2017 
                                    !inc = inc + merge( (1._wp+(1._wp-tausv)*rd)*intfund, merge( 0._wp, (1._wp+rd)*a, a>0._wp), intfund>0._wp) & ! 10.13.2017 
                                    !      + (1._wp-deltah)*h - ttaxent(inc)                                 
                                elseif(mode6taskid==1)then ! wealth tax (not applied on period income but beginning capital, regardless of asset types).
                                    inc = bgp + merge(0._wp, transbeq, printout17==.False. .and. t/=1) + wage * labsup - ssbtax ! 4.10.2017 
                                    !inc = inc + (1._wp - taubal) * (merge(bgp, 0._wp, bgp > 0._wp) + merge((1._wp + (1._wp - tausv) * rd) * intfund, 0._wp, intfund > 0._wp)) & ! 10.13.2017 
                                    !      + (1._wp - taubal) * (1._wp - deltah) * h - ttaxent(inc) + merge(bgp, 0._wp, bgp <= 0._wp)                                      
                                    inc = inc + merge( (1._wp+rd)*intfund, 0._wp, intfund>0._wp) + (1._wp - deltah) * h
                                    inc = (1._wp - taubal) * inc
                                elseif(mode6taskid==2)then
                                    net_worth = bgp + merge(transbeq,0._wp,printout17) + merge((1._wp + rd) * intfund, 0._wp, intfund>0._wp) + (1._wp - deltah) * h + wage*labsup - ssbtax
                                    inc = bgp + merge(transbeq,0._wp,printout17) + merge(rd*intfund,0._wp,intfund>0._wp) + wage*labsup - ssbtax 
                                    if(net_worth > 5._wp*exempbar)then
                                        inc = inc + merge(intfund,0._wp,intfund>0._wp) + (1._wp-deltah)*h - ttaxent(inc) - 0.02_wp*(net_worth-5._wp*exempbar) - 0.01_wp*(4._wp*exempbar)   
                                    elseif(net_worth>exempbar .and. net_worth<=5._wp*exempbar)then
                                        inc = inc + merge(intfund,0._wp,intfund>0._wp) + (1._wp-deltah)*h - ttaxent(inc) - 0.01_wp*(net_worth - exempbar)   
                                    else
                                        inc = inc + merge(intfund,0._wp,intfund>0._wp) + (1._wp-deltah)*h - ttaxent(inc)
                                    endif                                     
                                    
                                    !inc = wage*labsup + bgp + merge(0._wp, transbeq, printout17==.False. .and. t/=1) - ssbtax ! 2-24-2018
                                    !net_worth = (1._wp - deltah) * h + merge((1._wp + rd) * intfund, merge(0._wp, (1._wp + rd) * a, a > 0._wp), intfund > 0._wp)
                                    !if(net_worth > exempbar)then
                                    !    inc = inc + (1._wp - taubal) * net_worth - ttaxent(inc)
                                    !else
                                    !    net_worth = (1._wp - deltah) * h + merge((1._wp + (1._wp - tausv) * rd) * intfund, merge(0._wp, (1._wp + rd) * a, a > 0._wp), intfund > 0._wp)
                                    !    inc = inc + net_worth - ttaxent(inc)
                                    !endif                                     
                                else
                                    print*, "error in taxation 18-8"    
                                endif !mode6taskid
                            else
                                inc = inc + bgp ! 4.10.2017
                                !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       (1._wp+(1._wp-tausv)*rd)*k, intfund>0._wp), & ! removed 9-12-2017
                                
                                !!!10.13.2017
                                !inc = inc + merge(merge((1._wp+(1._wp-tausv)*rd)*a,       merge((1._wp+(1._wp-tausv)*rd)*k, (1._wp+(1._wp-tausv)*rd)*k+(1._wp+rd)*a ,a>0._wp), intfund>0._wp), & ! added 9-12-2017
                                !                  merge((1._wp+(1._wp-tausv)*rd)*intfund,                                                                               0._wp, intfund>0._wp), zx==0) &
                                !      + (1._wp-deltah)*h - ttaxent(inc) + merge(0._wp, transbeq, printout17==.False. .and. t/=1) - ssbtax ! 4.10.2017                            
                                
                                inc = inc + merge( (1._wp+(1._wp-tausv)*rd)*intfund, merge( 0._wp, (1._wp+rd)*a, a>0._wp), intfund>0._wp) & ! 10.13.2017 
                                      + (1._wp-deltah)*h - ttaxent(inc) + merge(0._wp, transbeq, printout17==.False. .and. t/=1) - ssbtax ! 4.10.2017                            
                            endif
                                
                            !inc = merge(wage*labsup,0._wp,zx==0) + merge(0._wp,capgain,tausvflag) + c_grs_mat(ax,kx,zx,yx,t) ! taxable income
                            !inc = merge(0._wp, capgain, tausvflag) + wage*labsup ! 3.14.2017 
                            
                            !! 3.14.2017 comment out.
                            !if(zx==0)then
                            !    ssbtax = tauss*wage*labsup    
                            !else
                            !    labdem = c_lab_vec(kx)
                            !    ssbtax = merge(tauss*wage*(labdem-labsup)+2._wp*tauss*wage*labsup,2*tauss*wage*labdem,labdem>labsup)
                            !endif
                            
                            !inc = inc + a + (1._wp-tausv)*merge(capgain,0._wp,tausvflag) + (1._wp-deltah)*h + (1._wp-taubp)*c_grs_mat(ax, kx, zx, yx, t) - ttaxent(inc) + transbeq  - ssbtax ! 3.14.2017 added the business after-tax profit.                   
                        endif 
                        
                        if(inc>0._wp)then
                            if(opx/=0)then ! opx is "given." not chosen. it corresponds to the level of kpx actually. !!!10.26.2017
                                pzv(1) = pz2(kpx,1) ! probability of technology shock is determined by the size of next-period's business project.
                                pzv(2) = pz2(kpx,2) ! probability of technology shock is determined by the size of next-period's business project.    
                                do i = 1, 2
                                    v2d(:,:,i) = wf(:,:,kpx,i-1,ypx,t+1)
                                enddo                                
                            else
                                v2d(:,:,1) = wf(:,:,kpx,0,ypx,t+1)
                            endif ! opx
                            
                            stp1 = 1
                            log1 = .false.
                            svcmp= 100._wp
                            allocate( shv(nsdim), lowap(nsdim), highap(nsdim), smaxv(nsdim), apvec(nsdim) ) ! first run                                  
                            do while(stp1<iterasvmax.and.log1==.false.)
                                
                                ! ACCURACY BLOCK --- START
                                if(stp1>1)then
                                    best = shv(loc1(1))
                                    !if(printout7) write(unit=109,fmt='(a,2i4,f12.8)'), ' new stp old idx ', stp1, loc1(1), best
                                    if(loc1(1)/=1.and.loc1(1)/=dsdim)then
                                        bdis = abs(shv(loc1(1)+1)-shv(loc1(1)-1))
                                    else
                                        if(loc1(1)==1)then
                                            bdis = abs(shv(1)-shv(2))       
                                        else ! loc1(1)==dsdim
                                            bdis = abs(shv(dsdim)-shv(dsdim-1))
                                        endif
                                    endif
                                    
                                    !if(loc1(1)==1)then
                                    !    bestl = best
                                    !    bestu = best + 0.8_wp*bdis
                                    !elseif(loc1(1)==dsdim)then
                                    !    bestu = best
                                    !    bestl = best - 0.8_wp*bdis
                                    !else
                                    !    bestl = best - 0.8_wp*bdis
                                    !    bestu = best + 0.8_wp*bdis
                                    !endif
                                    
                                    if(loc1(1)==1)then
                                        !bestl = best
                                        bestl = best - 0.2_wp*bdis ! 10012016
                                        bestu = best + 0.8_wp*bdis
                                        if(bestl<hv(1)) bestl = hv(1) ! 10012016
                                    elseif(loc1(1)==dsdim)then
                                        !bestu = best
                                        bestu = best + 0.2_wp*bdis ! 10012016
                                        bestl = best - 0.8_wp*bdis
                                        if(bestu>hv(hdim)) bestu = hv(hdim) ! 10012016
                                    else
                                        bestl = best - 0.8_wp*bdis
                                        bestu = best + 0.8_wp*bdis
                                        if(bestl<hv(1)) bestl = hv(1) ! 10012016
                                        if(bestu>hv(hdim)) bestu = hv(hdim) ! 10012016
                                    endif                                    

                                    dsdim = sdim ! bigger grid
                                    deallocate( shv, lowap, highap, smaxv, apvec )
                                    allocate( shv(sdim), lowap(sdim), highap(sdim), smaxv(sdim), apvec(sdim) ) ! small region brent's method --- part II
                                    !if(printout7) write(unit=109,fmt='(a,3f8.4)'), ' boundary ', bestl, bestu, bdis
                                    !call grid(shv,bestl,bestu,1._wp) ! good 3.14.2017 keep using 1._wp with the same reason stated in periods 10-14.
                                    call grid(shv,bestl,bestu,gsdim)
                                else ! stp1 = 1
                                    dsdim = nsdim ! smaller grid              
                                    if(hx==1)then
                                        call grid(shv,hv(1),hv(hdim),gsdim)    
                                    else ! hx
                                        !if(cwh(ax,hx-1,kx,zx,yx,kpx,ypx,opx,t)/=penalty)then ! comment out 10092016.
                                        !    call grid(shv,cwh(ax,hx-1,kx,zx,yx,kpx,ypx,opx,t),hv(hdim),gsdim)
                                        !else
                                            call grid(shv,hv(1),hv(hdim),gsdim) ! the default full range
                                        !endif ! cwh                                 
                                    endif ! hx
                                endif ! stp
                                ! ACCURACY BLOCK --- END                                
                                
                                smaxv = penalty
                                lowap  = -1
                                highap = -1

                                !if(printout7) write(unit=109,fmt='(a,9a)'), ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                                !if(printout7) write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,kpx,ypx,opx,t
                                
                                !debcnt1 = 0
                                do m = 1, dsdim
                                    hp = shv(m) ! pass into the thread
                                    call col_valid_fin(shv(m),kpx,ypx,t,lowap(m))
                                    call csp_valid_fin(h,shv(m),inc,highap(m))      
                                    !if(printout7) write(unit=109,fmt='(a,2f8.4,i5)'), ' -- ', lowap(m), highap(m), merge(0,1,lowap(m)<highap(m))
                                    if(lowap(m)>=highap(m).or.highap(m)<0._wp)then
                                        !if(printout7) write(unit=109,fmt='(i3)'), debcnt1
                                        cycle ! negative consumption in one of the refined housing assets level 
                                    endif
                                    call brent_localizer(f18,lowap(m),highap(m),apvec(m),smaxv(m),exit_brent) ! <-------------------------------------------- [6]
                                    !if(printout7) write(unit=109,fmt='(a,2f8.4)'), '    ', apvec(m), smaxv(m)
                                    !if(printout7.and.m==sdim) write(unit=109,fmt='(a)'), '    '
                                    
                                    if(exit_brent) exit ! exit this loop over m. ! [7]
                                    
                                enddo
                                
                                if(exit_brent==.false.)then ! [8]
                                    ! get the maximum from the set of local bests.
                                    if(maxval(smaxv)/=penalty)then
                                        loc1 = maxloc(smaxv)
                                        stp1 = stp1 + 1
                                        if(abs(svcmp-smaxv(loc1(1)))<err_svdist) exit
                                        svcmp = smaxv(loc1(1))                            
                                    else
                                        log1 = .true.
                                        svcmp= penalty
                                        !if(printout7.and.debcnt1==sdim) write(unit=109,fmt='(a)'), ' all penalty  ' ! correct 09222016
                                        !if(printout7.and.debcnt1/=sdim) write(unit=109,fmt='(a,10i3)'), ' something wrong in smaxv ', ax,hx,kx,zx,yx,kpx,ypx,opx,t, debcnt1 ! no leaks 09222016
                                    endif
                                else                                                           !
                                    exit ! exit this housing asset housing grid searching zone ! [9]    
                                endif ! exit_brent                                             !
                            enddo ! do while                              
                            
                            if(exit_brent==.false.)then ! [10]
                                ! save searching results or give a warning mark.
                                if(svcmp/=penalty)then ! SOME CASES NEED TO DO COMPARISON 
                                    cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = smaxv(loc1(1))
                                    cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = apvec(loc1(1))
                                    cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = shv(loc1(1))  
                                    cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = -99
                                    if(opx==0)then
                                        cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 0
                                        cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = inc - (apvec(loc1(1))+shv(loc1(1))) - HomTrnsCst(h,shv(loc1(1)))                                
                                    else ! opx/=0
                                        cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = kpx
                                        cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = inc - (apvec(loc1(1))+shv(loc1(1))) - HomTrnsCst(h,shv(loc1(1)))                                
                                        if(opx==1)then ! just choose to be a worker instead
                                            if( cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) < cwf(ax,hx,kx,zx,yx,0,ypx,0,t))then ! 3.14.2017 Entrepreneurs > Workers > Retirees.
                                                cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwf(ax,hx,kx,zx,yx,0,ypx,0,t) 
                                                cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwa(ax,hx,kx,zx,yx,0,ypx,0,t)
                                                cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwh(ax,hx,kx,zx,yx,0,ypx,0,t)
                                                cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwk(ax,hx,kx,zx,yx,0,ypx,0,t)
                                                cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cww(ax,hx,kx,zx,yx,0,ypx,0,t)
                                                cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwc(ax,hx,kx,zx,yx,0,ypx,0,t)
                                            endif
                                        else ! opx==2 ! choose to engage in a smaller investment project.
                                            opxi = merge(opx-1,opx-2,zx==1) ! 3.31.2017
                                            if( cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) < cwf(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t) )then
                                                cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwf(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t) 
                                                cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwa(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)
                                                cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwh(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)
                                                cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwk(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)                                            
                                                cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cww(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t) 
                                                cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwc(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t) 
                                            endif
                                        endif ! opx==1
                                    endif ! opx/=0
                                else ! infeasible point due to collateral and budget constraints. DOWNGRADE WITHOUT ANY COMPARISON
                                    if(opx==0)then ! no where to downgrade their expenses or costs
                                        cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 1 ! -99 indicates a normal case.
                                    elseif(opx==1)then ! may instead choose to be a worker
                                        cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwf(ax,hx,kx,zx,yx,0,ypx,0,t) 
                                        cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwa(ax,hx,kx,zx,yx,0,ypx,0,t)
                                        cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwh(ax,hx,kx,zx,yx,0,ypx,0,t)
                                        cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cww(ax,hx,kx,zx,yx,0,ypx,0,t)   
                                        cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwk(ax,hx,kx,zx,yx,0,ypx,0,t)                                    
                                        cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwc(ax,hx,kx,zx,yx,0,ypx,0,t)   
                                    elseif(opx==2)then ! may instead choose to engage in an investment project of smaller scale.
                                        opxi = merge(opx-1,opx-2,zx==1) ! 3.31.2017
                                        cwf(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwf(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t) 
                                        cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwa(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)
                                        cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwh(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)
                                        cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cww(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)   
                                        cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwk(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)                                                                                
                                        cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cwc(ax,hx,kx,zx,yx,kpx-1,ypx,opxi,t)   
                                    endif    
                                endif 
                            endif ! exit_brent [11]
                        
                            deallocate( shv, lowap, highap, smaxv, apvec )  
                            
                        else ! inc<0._wp
                            cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 1                        
                        endif ! inc
                        
                        if(exit_brent) exit ! exit this n loop (over the indices of housing asset holding grid). 8-13-2017 [12]
                        
                    enddo ! n
                    
                    if(exit_brent) qray(q)=.true. ! [13]
                    
                enddo ! q  
                !$omp end do
                deallocate( v2d )
                !$omp end parallel
                
                if(any(qray==.true.)) exit_bellman=.true. ! [14]
                if(exit_bellman) exit ! [15]
                
            enddo ! l
            
            if(exit_bellman==.false.)then ! [16]
                if(printout2)then            
                    do l = 1, 4                
                        do q = bd(l,1), bd(l,2) ! income in period 14 could be used later.
                            ax  = xlist(q,1)
                            hx  = xlist(q,2)
                            kx  = xlist(q,3)
                            zx  = xlist(q,4)
                            yx  = xlist(q,5)
                            kpx = xlist(q,6)
                            ypx = xlist(q,7)
                            opx = xlist(q,8)  
                                     
                            if(ax==1.and.hx==1)then
                                write(unit=056,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                                write(unit=056,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                                write(unit=056,fmt='(3x,9x,'//str1//'i21)') (j,j=1,hdim) 
                                write(unit=056,fmt='(3x,9x,'//str1//'f21.4)') (hv(j),j=1,hdim) 
                                do i = 1, adim
                                    write(unit=056,fmt='(i3,x,f8.4,'//str1//'(x,f20.8))') i,av(i),(cwf(i,j,kx,zx,yx,kpx,ypx,opx,t),j=1,hdim)    
                                enddo
                            endif 
                            
                        enddo ! q            
                    enddo ! l
                endif ! printout2  
                

                do n = 1, 7
                    !$omp parallel default(shared) private(i,j,vpe1,vpe2,pyv)
                    kx = tvector(n,1)
                    zx = tvector(n,2)
                    do yx = 1, nmc            
                        !$omp do
                        do q = 1, hdim*adim  
                            call acc_index_begin_mass_per_period(q,i,j)
                            pyv = py(yx,:)
                            if(n==1)then
                                    vpe1 = cwf(i,j,kx,zx,yx,0,1:nmc,0,t) ! still a worker
                                    vpe2 = cwf(i,j,kx,zx,yx,1,1:nmc,2,t) ! turn to be an entrepreneur
                                    if(all(cww(i,j,kx,zx,yx,0,1:nmc,0,t)==-99).and.all(cww(i,j,kx,zx,yx,1,1:nmc,2,t)==-99)) &
                                        & wf(i,j,kx,zx,yx,t) = (1._wp-probtk(kx))*dot_product(pyv,vpe1) + probtk(kx)*dot_product(pyv,vpe2)
                            elseif(2<=n.and.n<=4)then ! 3.31.2017                    
                                    vpe1 = cwf(i,j,kx,zx,yx,0,1:nmc,0,t)
                                    vpe2 = cwf(i,j,kx,zx,yx,1,1:nmc,2,t)
                                    if(all(cww(i,j,kx,zx,yx,0,1:nmc,0,t)==-99).and.all(cww(i,j,kx,zx,yx,1,1:nmc,2,t)==-99)) &
                                        & wf(i,j,kx,zx,yx,t) = (1._wp-probtk(0))*dot_product(pyv,vpe1) + probtk(0)*dot_product(pyv,vpe2) ! 3.31.2017 correct the error that probtk's index is not 0.
                            elseif(5<=n.and.n<=6)then                     
                                    vpe1 = cwf(i,j,kx,zx,yx,kx,1:nmc,1,t)
                                    vpe2 = cwf(i,j,kx,zx,yx,kx+1,1:nmc,2,t)
                                    if(all(cww(i,j,kx,zx,yx,kx,1:nmc,1,t)==-99).and.all(cww(i,j,kx,zx,yx,kx+1,1:nmc,2,t)==-99)) &
                                        & wf(i,j,kx,zx,yx,t) = (1._wp-probtk(kx))*dot_product(pyv,vpe1) + probtk(kx)*dot_product(pyv,vpe2)             
                            elseif(n==7)then                     
                                    vpe1 = cwf(i,j,kx,zx,yx,kx,1:nmc,1,t)
                                    if(all(cww(i,j,kx,zx,yx,kx,1:nmc,1,t)==-99)) &
                                        & wf(i,j,kx,zx,yx,t) = dot_product(pyv,vpe1)           
                            endif
                        enddo ! q
                        !$omp end do
                    enddo ! yx
                    !$omp end parallel
                enddo ! n
                    
                do n = 1, 7
                    kx = tvector(n,1)
                    zx = tvector(n,2)
                    do yx = 1, nmc                    
                        if(printout2)then
                            write(unit=026,fmt='(a,3(i3),a,i3)') '(kx,zx,yx) ',kx,zx,yx, ' t ',t     
                            write(unit=026,fmt='(3x,8x,'//str1//'(8x,i8))') (j,j=1,hdim)             
                            write(unit=026,fmt='(3x,5x,"999",'//str1//'(8x,f8.4))') (hv(j),j=1,hdim)             
                            do j = 1, adim
                                write(unit=026,fmt='(i3,x,f7.2,'//str1//'(x,f15.5))') j, av(j),(wf(j,i,kx,zx,yx,t),i=1,hdim)
                            enddo                                  
                        endif ! printout2                      
                    enddo ! yx              
                enddo ! n
            endif ! exit_bellman [17]
            
            call system_clock(tend)
            !if(printout6) write(unit=120,fmt='(a,i3,a,f12.4,a)') 'solve bellman period ',t,', time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'     
            if(printout6) write(*,fmt='(a,i3,a,f12.4,a)') 'solve bellman period ',t,', time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'  
            
            if(exit_bellman) exit ! exit the loop over t. [18]
            
        enddo ! t
        deallocate( qray, xlist, bd ) ! [19]    
        
    end subroutine solve_bellman_18
    
    function f18(x)
        implicit none
        real(wp), intent(in) :: x
        real(wp) :: f18, csp, vp(2), vec(4)
        integer :: n
        logical :: msg
        
        csp = inc - (x+hp) - HomTrnsCst(h,hp)
        
        if(csp<=0._wp.or.(x+hp)<0._wp)then ! 3.14.2017 replace x+hp<=0 with x+hp<0
            f18 = penalty
            return
        endif         
        
        if(opx==0)then

            vp(1) = binterpII(av,hv,v2d(:,:,1),x,hp,penalty,msg,vec)    
            !if(msg==.false.)then
            !    f18 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp) + survprob(t)*beta*vp(1)
            !else
            !    f18 = penalty
            !endif      
            f18 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp) + survprob(t)*beta*vp(1)
        else ! opx/=0
            do n = 1, 2
                vp(n) = binterpII(av,hv,v2d(:,:,n),x,hp,penalty,msg,vec)
            enddo ! n
            
            !if(msg==.false..and.all(vp/=penalty))then
            !    f18 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp) + survprob(t)*beta*dot_product(pzv,vp) ! dot_product(pzv,dum_wf)               
            !else
            !    f18 = penalty
            !endif       
            f18 = u(csp,h) + (1._wp-survprob(t))*bu(x+hp) + survprob(t)*beta*dot_product(pzv,vp) ! dot_product(pzv,dum_wf)     
        endif
    end function f18 ! 3.14.2017 11:00 pm stop here.   
    
    subroutine acc_index_begin_mass_per_period(q,i,j)
        implicit none
        integer, intent(in) :: q
        integer, intent(out) :: i, j
        ! =IF(MOD(B2,$F$2)=0,$F$2,MOD(B2,$F$2))
        i = merge(adim,mod(q,adim),mod(q,adim)==0)
        ! =(B2-C2)/$F$2+1
        j = (q-i)/adim+1
    end subroutine acc_index_begin_mass_per_period
    
    subroutine acc_index_begin_mass_per_period_fnadim(q,i,j)
        implicit none
        integer, intent(in) :: q
        integer, intent(out) :: i, j
        i = merge(fnadim,mod(q,fnadim),mod(q,fnadim)==0)
        ! =(B2-C2)/$F$2+1
        j = (q-i)/fnadim+1        
    end subroutine acc_index_begin_mass_per_period_fnadim
        
    subroutine print_coarse_2d_brent_mat(iterar,iteragov,iteratot)
        implicit none
        integer, intent(in) :: iterar, iteragov, iteratot
        integer, dimension(:,:), allocatable :: tarray
        integer :: i, j, l, m, n, t, upbnd
        character(len=2) :: str1
        write(str1,fmt='(i2)') hdim
        
        !vtwa = twa
        !where(vtwa==penalty) vtwa=av(1)
        !vtwh = twh
        !where(tww==1) vtwh=1
        !vtwk = twk
        !where(tww==1) vtwk=0            
        
        do t = 14, 1, -1 ! <---------------------
            if(t<=9)then
                upbnd = 13
                allocate( tarray(13,6) )
                
                tarray(:,1) = [0,1,2,3,1,2,3,0,1,2,1,2,3] ! k
                tarray(:,2) = [0,0,0,0,1,1,1,0,1,1,0,0,0] ! z
                tarray(:,3) = 1                     ! y ! 4.9.2017 it doesn't matter what value it takes.
                tarray(:,4) = [0,0,0,0,1,2,3,1,2,3,1,1,1] ! kp
                tarray(:,5) = 1                     ! yp ! 4.9.2017 it doesn't matter what value it takes.
                tarray(:,6) = [0,0,0,0,1,1,1,2,2,2,2,2,2] ! op   
                
                !tarray(:,1) = [0,1,2,3,1,2,3,1,2,3,0,1,2,1,2,3]
                !tarray(:,2) = [0,0,0,0,1,1,1,1,1,1,0,1,1,0,0,0]
                !tarray(:,3) = 1 ! changes below ! 4.9.2017 it doesn't matter what value it takes.
                !tarray(:,4) = [0,0,0,0,0,0,0,1,2,3,1,2,3,1,1,1]
                !tarray(:,5) = 0 ! changes below ! 4.9.2017 it doesn't matter what value it takes.
                !tarray(:,6) = [0,0,0,0,0,0,0,1,1,1,2,2,2,2,2,2]    
                
            elseif(t>=10.and.t<=13)then
                upbnd = 9 ! 12 
                allocate( tarray(9,6) ) ! 4.9.2017 replace 12.
                
                tarray(:,1) = [0,1,2,3,1,2,3,1,2] ! k
                tarray(:,2) = [0,0,0,0,1,1,1,1,1] ! z
                tarray(:,3) = 1                   ! y
                tarray(:,4) = [0,0,0,0,1,2,3,2,3] ! kp
                tarray(:,5) = 1                   ! yp
                tarray(:,6) = [0,0,0,0,1,1,1,2,2] ! op   
                    
                !tarray(:,1) = [0,1,2,3,1,2,3,1,2,3,1,2] ! kx #2
                !tarray(:,2) = [0,0,0,0,1,1,1,1,1,1,1,1] ! zx
                !tarray(:,3) = 0 ! yx changes below ! 4.9.2017 it doesn't matter what value it takes.
                !tarray(:,4) = [0,0,0,0,0,0,0,1,2,3,2,3] ! kp
                !tarray(:,5) = 0 ! yp changes below ! 4.9.2017 it doesn't matter what value it takes.
                !tarray(:,6) = [0,0,0,0,0,0,0,1,1,1,2,2] ! op                
            
            elseif(t==14)then
                upbnd = 7
                allocate( tarray(7,6) )
                tarray(:,1) = [0,1,2,3,1,2,3] ! k
                tarray(:,2) = [0,0,0,0,1,1,1] ! z
                tarray(:,3) = 1               ! y changes below ! 4.9.2017 it doesn't matter what value it takes.
                tarray(:,4) = [0,0,0,0,0,0,0] ! kp
                tarray(:,5) = 1               ! yp changes below ! 4.9.2017 it doesn't matter what value it takes.
                tarray(:,6) = [0,0,0,0,0,0,0] ! op                  
            endif
            
            do l = 1,upbnd
                kx  = tarray(l,1)
                zx  = tarray(l,2)
                yx  = tarray(l,3) ! yx changes below
                kpx = tarray(l,4)
                ypx = tarray(l,5) ! ypx changes below
                opx = tarray(l,6)   
                
                do i = 1, nmc ! yx
                    do j = 1, nmc ! ypx
                        if(t==9)then
                            if(j==1)then ! 4.9.217 correct
                                yx  = i
                                ypx = 0 ! 4.9.2017 correct only happen once for each i when j==1
                            else ! j/=1
                                cycle
                            endif
                        elseif(t>=10)then
                            if(i==1.and.j==1)then
                                yx  = 0 ! 4.9.2017 correct no labor efficiency this period.
                                ypx = 0 ! 4.9.2017 correct no labor efficiency this period.
                            else
                                cycle    
                            endif
                        else ! t<=8
                            yx  = i ! 4.9.2017 cycles in both the current and next periods
                            ypx = j ! 4.9.2017 cycles in both the current and next periods  
                        endif
                        
                        write(unit=115,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                        write(unit=115,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                        write(unit=115,fmt='(3x,7x,'//str1//'i13)') (n,n=1,hdim) 
                        write(unit=115,fmt='(3x,i7,'//str1//'f13.2)') -999, (hv(n),n=1,hdim)    
                        
                        write(unit=116,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                        write(unit=116,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                        write(unit=116,fmt='(3x,7x,'//str1//'i13)') (n,n=1,hdim) 
                        write(unit=116,fmt='(3x,i7,'//str1//'f13.3)') -999, (hv(n),n=1,hdim) 
                        
                        write(unit=117,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                        write(unit=117,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                        write(unit=117,fmt='(3x,7x,'//str1//'i13)') (n,n=1,hdim) 
                        write(unit=117,fmt='(3x,i7,'//str1//'f13.3)') -999, (hv(n),n=1,hdim) 
                        
                        write(unit=118,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                        write(unit=118,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                        write(unit=118,fmt='(3x,7x,'//str1//'i13)') (n,n=1,hdim) 
                        write(unit=118,fmt='(3x,i7,'//str1//'f13.3)') -999, (hv(n),n=1,hdim) 
                        
                        write(unit=119,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                        write(unit=119,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                        write(unit=119,fmt='(3x,7x,'//str1//'i13)') (n,n=1,hdim) 
                        write(unit=119,fmt='(3x,i7,'//str1//'f13.2)') -999, (hv(n),n=1,hdim)                          
                        
                        write(unit=123,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                        write(unit=123,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                        write(unit=123,fmt='(3x,7x,'//str1//'i13)') (n,n=1,hdim) 
                        write(unit=123,fmt='(3x,i7,'//str1//'f13.2)') -999, (hv(n),n=1,hdim)      
                        
                        write(unit=124,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                        write(unit=124,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                        write(unit=124,fmt='(3x,7x,'//str1//'i13)') (n,n=1,hdim) 
                        write(unit=124,fmt='(3x,i7,'//str1//'f13.2)') -999, (hv(n),n=1,hdim)                          
                        
                        do m = 1, adim
                            
                            write(unit=115,fmt='(i3,x,f6.2,'//str1//'(x,i12  ))')   m,av(m),(cww(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim)  
                            write(unit=116,fmt='(i3,x,f6.2,'//str1//'(x,f12.3))')   m,av(m),(cwf(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim)  
                            write(unit=117,fmt='(i3,x,f6.2,'//str1//'(x,f12.3))')   m,av(m),(cwa(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim)  
                            write(unit=118,fmt='(i3,x,f6.2,'//str1//'(x,f12.3))')   m,av(m),(cwh(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim)  
                            write(unit=119,fmt='(i3,x,f6.2,'//str1//'(x,i12  ))')   m,av(m),(cwk(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim)                          
                            write(unit=123,fmt='(i3,x,f6.2,'//str1//'(x,f12.3  ))') m,av(m),(cwc(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim) 
                            write(unit=124,fmt='(i3,x,f6.2,'//str1//'(x,f12.9  ))') m,av(m),(cef(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim) 
                            
                            !write(unit=115,fmt='(f12.5,'//str1//'(x,i12))')   av(m), (cww(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim)
                            !write(unit=116,fmt='(f12.5,'//str1//'(x,es12.5))')   av(m), (cwf(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim)   
                            !write(unit=117,fmt='(f12.5,'//str1//'(x,es12.5))')   av(m), (cwa(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim)
                            !write(unit=118,fmt='(f12.5,'//str1//'(x,es12.5))')   av(m), (cwh(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim)
                            !write(unit=119,fmt='(f12.5,'//str1//'(x,i12))')   av(m), (cwk(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim)
                        enddo ! m
                    enddo ! j
                enddo ! i
            enddo ! l
            deallocate( tarray )
            
            ! cvv -------------------------
            if(1<=t.and.t<=9)then
                do n = 1, 7
                    kx = tvector(n,1)
                    zx = tvector(n,2)
                    do m = 1, nmc
                        yx = m
                        !if(printout1)then
                        write(unit=125,fmt='(a,3(i3),a,i3)') '(kx,zx,yx) ',kx,zx,yx, ' t ',t     
                        write(unit=125,fmt='(3x,'//str1//'(x,i4))') (j,j=1,hdim)             
                        do j = 1, adim
                            write(unit=125,fmt='(i3,'//str1//'(x,i4))') j, (cvv(j,i,kx,zx,yx,t),i=1,hdim)
                        enddo                                   
                        !endif ! printout1                           
                    enddo ! m        
                enddo ! n                 
            else
                do n = 1, 7
                    kx = tvector(n,1)
                    zx = tvector(n,2)
                    yx = 0
                    
                    write(unit=125,fmt='(a,3(i3),a,i3)') '(kx,zx,yx) ',kx,zx,yx, ' t ',t     
                    write(unit=125,fmt='(3x,'//str1//'(x,i4))') (j,j=1,hdim)             
                    do j = 1, adim
                        write(unit=125,fmt='(i3,'//str1//'(x,i4))') j, (cvv(j,i,kx,zx,yx,t),i=1,hdim) 
                    enddo                                   
                enddo                  
            endif ! wf condition            
        enddo ! t
    end subroutine print_coarse_2d_brent_mat
    
    subroutine print_2d_end_of_period_dist_mat(iterar,iteragov,iteratot)
        implicit none
        integer, intent(in) :: iterar, iteragov, iteratot
        integer, dimension(:,:), allocatable :: tarray
        integer :: i, j, l, m, n, t, upbnd
        character(len=2) :: str1
        write(str1,fmt='(i2)') hdim
        
        !vtwa = twa
        !where(vtwa==penalty) vtwa=av(1)
        !vtwh = twh
        !where(tww==1) vtwh=1
        !vtwk = twk
        !where(tww==1) vtwk=0            
        
        do t = 14, 1, -1 ! <---------------------
            if(t<=9)then
                upbnd = 13
                allocate( tarray(13,6) )
                
                tarray(:,1) = [0,1,2,3,1,2,3,0,1,2,1,2,3] ! k
                tarray(:,2) = [0,0,0,0,1,1,1,0,1,1,0,0,0] ! z
                tarray(:,3) = 1                     ! y ! 4.9.2017 it doesn't matter what value it takes.
                tarray(:,4) = [0,0,0,0,1,2,3,1,2,3,1,1,1] ! kp
                tarray(:,5) = 1                     ! yp ! 4.9.2017 it doesn't matter what value it takes.
                tarray(:,6) = [0,0,0,0,1,1,1,2,2,2,2,2,2] ! op   
                
                !tarray(:,1) = [0,1,2,3,1,2,3,1,2,3,0,1,2,1,2,3]
                !tarray(:,2) = [0,0,0,0,1,1,1,1,1,1,0,1,1,0,0,0]
                !tarray(:,3) = 1 ! changes below ! 4.9.2017 it doesn't matter what value it takes.
                !tarray(:,4) = [0,0,0,0,0,0,0,1,2,3,1,2,3,1,1,1]
                !tarray(:,5) = 0 ! changes below ! 4.9.2017 it doesn't matter what value it takes.
                !tarray(:,6) = [0,0,0,0,0,0,0,1,1,1,2,2,2,2,2,2]    
                
            elseif(t>=10.and.t<=13)then
                upbnd = 9 ! 12 
                allocate( tarray(9,6) ) ! 4.9.2017 replace 12.
                
                tarray(:,1) = [0,1,2,3,1,2,3,1,2] ! k
                tarray(:,2) = [0,0,0,0,1,1,1,1,1] ! z
                tarray(:,3) = 1                   ! y
                tarray(:,4) = [0,0,0,0,1,2,3,2,3] ! kp
                tarray(:,5) = 1                   ! yp
                tarray(:,6) = [0,0,0,0,1,1,1,2,2] ! op   
                    
                !tarray(:,1) = [0,1,2,3,1,2,3,1,2,3,1,2] ! kx #2
                !tarray(:,2) = [0,0,0,0,1,1,1,1,1,1,1,1] ! zx
                !tarray(:,3) = 0 ! yx changes below ! 4.9.2017 it doesn't matter what value it takes.
                !tarray(:,4) = [0,0,0,0,0,0,0,1,2,3,2,3] ! kp
                !tarray(:,5) = 0 ! yp changes below ! 4.9.2017 it doesn't matter what value it takes.
                !tarray(:,6) = [0,0,0,0,0,0,0,1,1,1,2,2] ! op                
            
            elseif(t==14)then
                upbnd = 7
                allocate( tarray(7,6) )
                tarray(:,1) = [0,1,2,3,1,2,3] ! k
                tarray(:,2) = [0,0,0,0,1,1,1] ! z
                tarray(:,3) = 1               ! y changes below ! 4.9.2017 it doesn't matter what value it takes.
                tarray(:,4) = [0,0,0,0,0,0,0] ! kp
                tarray(:,5) = 1               ! yp changes below ! 4.9.2017 it doesn't matter what value it takes.
                tarray(:,6) = [0,0,0,0,0,0,0] ! op                  
            endif
            
            do l = 1,upbnd
                kx  = tarray(l,1)
                zx  = tarray(l,2)
                yx  = tarray(l,3) ! yx changes below
                kpx = tarray(l,4)
                ypx = tarray(l,5) ! ypx changes below
                opx = tarray(l,6)   
                
                do i = 1, nmc ! yx
                    do j = 1, nmc ! ypx
                        if(t==9)then
                            if(j==1)then ! 4.9.217 correct
                                yx  = i
                                ypx = 0 ! 4.9.2017 correct only happen once for each i when j==1
                            else ! j/=1
                                cycle
                            endif
                        elseif(t>=10)then
                            if(i==1.and.j==1)then
                                yx  = 0 ! 4.9.2017 correct no labor efficiency this period.
                                ypx = 0 ! 4.9.2017 correct no labor efficiency this period.
                            else
                                cycle    
                            endif
                        else ! t<=8
                            yx  = i ! 4.9.2017 cycles in both the current and next periods
                            ypx = j ! 4.9.2017 cycles in both the current and next periods  
                        endif
                        
                        write(unit=129,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                        write(unit=129,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                        write(unit=129,fmt='(3x,7x,'//str1//'i19)') (n,n=1,hdim) 
                        write(unit=129,fmt='(3x,i7,'//str1//'f19.2)') -999, (hv(n),n=1,hdim)    
                        
                        write(unit=130,fmt='(7a)') '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                        write(unit=130,fmt='(7i4)') kx,zx,yx,kpx,ypx,opx,t
                        write(unit=130,fmt='(3x,7x,'//str1//'i19)') (n,n=1,hdim) 
                        write(unit=130,fmt='(3x,i7,'//str1//'f19.2)') -999, (hv(n),n=1,hdim)                                        
                        
                        do m = 1, adim
                            write(unit=129,fmt='(i3,x,f6.2,'//str1//'(x,f18.12  ))')  m,av(m),(cef(m,n,kx,zx,yx,kpx,ypx,opx,t),n=1,hdim)  
                            write(unit=130,fmt='(i3,x,f6.2,'//str1//'(x,f18.12  ))')  m,av(m),(sw_nonlineartax(c3s(m,n,kx,zx,yx,kpx,ypx,opx,t)),n=1,hdim)  
                        enddo ! m
                    enddo ! j
                enddo ! i
            enddo ! l
            deallocate( tarray )
            
            !! cvv ------------------------- start of period 
            !if(1<=t.and.t<=9)then
            !    do n = 1, 7
            !        kx = tvector(n,1)
            !        zx = tvector(n,2)
            !        do m = 1, nmc
            !            yx = m
            !            !if(printout1)then
            !            write(unit=125,fmt='(a,3(i3),a,i3)') '(kx,zx,yx) ',kx,zx,yx, ' t ',t     
            !            write(unit=125,fmt='(3x,'//str1//'(x,i4))') (j,j=1,hdim)             
            !            do j = 1, adim
            !                write(unit=125,fmt='(i3,'//str1//'(x,i4))') j, (cvv(j,i,kx,zx,yx,t),i=1,hdim)
            !            enddo                                   
            !            !endif ! printout1                           
            !        enddo ! m        
            !    enddo ! n                 
            !else
            !    do n = 1, 7
            !        kx = tvector(n,1)
            !        zx = tvector(n,2)
            !        yx = 0
            !        
            !        write(unit=125,fmt='(a,3(i3),a,i3)') '(kx,zx,yx) ',kx,zx,yx, ' t ',t     
            !        write(unit=125,fmt='(3x,'//str1//'(x,i4))') (j,j=1,hdim)             
            !        do j = 1, adim
            !            write(unit=125,fmt='(i3,'//str1//'(x,i4))') j, (cvv(j,i,kx,zx,yx,t),i=1,hdim) 
            !        enddo                                   
            !    enddo                  
            !endif ! wf condition            
            
        enddo ! t
    end subroutine print_2d_end_of_period_dist_mat    
    
    ! 3.16.2017 check to see whether all the subsequent outcome related to a specific combination are valid.
    subroutine valid_beginning_period_mass(iterar, iteragov, iteratot) ! 8-26-2017 should begin here...
        implicit none
        integer, intent(in) :: iterar, iteragov, iteratot
        integer, dimension(:), allocatable :: twwvec
        integer, dimension(:,:), allocatable :: twwmat    
        integer :: t, n, yxi, l, nupbnd, m
        character(len=3) :: str1, str2
        integer :: i ! printout for debug 8-26-2017
        
        !real(wp), dimension(:), allocatable :: nafvvec
        !integer :: t, l, n, m, nupbnd, ans_h,ans_k,yxi, i
        !real(wp) :: ans_a
        !integer :: pt(2)
        !real(wp) :: wt(2), nafvvec_max, upper, lower
        !character(len=3) :: str1, str2
        m = 1 ! 3.17.2017
        write(str1,fmt='(i3)') hdim
        !call smi(tvector,'tvector')
        do t = 1, 14
            do n = 1, 7
                kx = tvector(n,1)
                zx = tvector(n,2)
                
                if(1<=t.and.t<=8)then
                    if(n==1)then ! 3.16.2017 szt...vec.. have all been checked!
                        allocate( twwvec(szt18vec1), twwmat(szt18vec1,4) ) ! , nafvvec(szt18vec1) )  
                        nupbnd = szt18vec1
                        twwvec = -99
                        twwmat = t18vec1
                    elseif(2<=n.and.n<=4)then
                        allocate( twwvec(szt18vec24), twwmat(szt18vec24,4) ) ! , nafvvec(szt18vec24) )    
                        nupbnd = szt18vec24
                        twwvec = -99
                        twwmat = t18vec24                       
                    elseif(5<=n.and.n<=6)then
                        allocate( twwvec(szt18vec56), twwmat(szt18vec56,4) ) ! , nafvvec(szt18vec56) )    
                        nupbnd = szt18vec56
                        twwvec = -99
                        twwmat = t18vec56                        
                    elseif(n==7)then
                        allocate( twwvec(szt18vec7), twwmat(szt18vec7,4) ) ! , nafvvec(szt18vec7) )    
                        nupbnd = szt18vec7
                        twwvec = -99
                        twwmat = t18vec7  
                    endif ! n condition
                elseif(t==9)then
                    if(n==1)then
                        allocate( twwvec(szt9vec1), twwmat(szt9vec1,4) )  ! , nafvvec(szt9vec1) )  
                        nupbnd = szt9vec1
                        twwvec = -99
                        twwmat = t9vec1
                    elseif(2<=n.and.n<=4)then
                        allocate( twwvec(szt9vec24), twwmat(szt9vec24,4) ) ! nafvvec(szt9vec24) )    
                        nupbnd = szt9vec24
                        twwvec = -99
                        twwmat = t9vec24                       
                    elseif(5<=n.and.n<=6)then
                        allocate( twwvec(szt9vec56), twwmat(szt9vec56,4) ) ! , nafvvec(szt9vec56) )    
                        nupbnd = szt9vec56
                        twwvec = -99
                        twwmat = t9vec56                        
                    elseif(n==7)then
                        allocate( twwvec(szt9vec7), twwmat(szt9vec7,4) ) ! , nafvvec(szt9vec7) )    
                        nupbnd = szt9vec7
                        twwvec = -99
                        twwmat = t9vec7                       
                    endif ! n condition                    
                elseif(10<=t.and.t<=13)then
                    if(n==1)then
                        allocate( twwvec(szt1013vec1), twwmat(szt1013vec1,4) ) ! , nafvvec(szt1013vec1) )  
                        nupbnd = szt1013vec1
                        twwvec = -99
                        twwmat = t1013vec1
                    elseif(2<=n.and.n<=4)then
                        allocate( twwvec(szt1013vec24), twwmat(szt1013vec24,4) ) ! , nafvvec(szt1013vec24) )    
                        nupbnd = szt1013vec24
                        twwvec = -99
                        twwmat = t1013vec24                       
                    elseif(5<=n.and.n<=6)then
                        allocate( twwvec(szt1013vec56), twwmat(szt1013vec56,4) ) ! , nafvvec(szt1013vec56) )    
                        nupbnd = szt1013vec56
                        twwvec = -99
                        twwmat = t1013vec56                        
                    elseif(n==7)then
                        allocate( twwvec(szt1013vec7), twwmat(szt1013vec7,4) ) ! , nafvvec(szt1013vec7) )    
                        nupbnd = szt1013vec7 
                        twwvec = -99
                        twwmat = t1013vec7                       
                    endif ! n condition 
                elseif(t==14)then    
                    allocate( twwvec(szt1013vec1), twwmat(szt1013vec1,4) ) ! , nafvvec(szt1013vec1) )
                        nupbnd = szt1013vec1
                        twwvec = -99
                        twwmat = t1013vec1                    
                endif                     
                
                do yxi = 1, nmc ! 3.16.2017 The current period index of variable "yx"
                    !! 3.16.2017 comment out. new block is shown right below this block.
                    !if(10<=t.and.yxi/=1)then 
                    !    cycle    
                    !elseif(10<=t.and.yxi==1)then
                    !    yx = 0
                    !else
                    !    yx = yxi
                    !endif
                    
                    ! The index of current labor efficiency
                    ! I use this trick to reduce the total levels of labor efficiency to only one for the case that agents are in retirement.
                    
                    if(10>t)then
                        yx = yxi    
                    elseif(t>=10)then
                        if(yxi/=1)then ! we only need to count retirenent once (and here we use yxi==1 to indicate the retirement case) 3.16.2017
                            cycle
                        else ! retired.
                            yx = 0
                        endif
                    endif
                    
                    do hx = 1, hdim
                        do ax = 1, adim
                            !write(unit=109,fmt='(a,9a)'), ' ---- ', '  ax','  hx','  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                            !write(unit=109,fmt='(a,9i4)'), ' ---- ', ax,hx,kx,zx,yx,kpx,ypx,opx,t                         
                            !if(ax==1.and.hx==1) write(unit=127,fmt='(9i4,a,i4)') ax,hx,kx,zx,yx,kpx,ypx,opx,t, ' -- ', yxi
                            
                            ! 4.1.2017 Seems useless.
                            s1c(m,1) = ax 
                            s1c(m,2) = hx
                            s1c(m,3) = kx
                            s1c(m,4) = zx
                            s1c(m,5) = yx
                            s1c(m,6) = t
                            s1c(m,7) = m
                            c1s(ax, hx, kx, zx, yx, t) = m
                            !write(unit = 127,fmt='(8(i8," "))') s1c(m,:), c1s(ax, hx, kx, zx, yx, t)
                            m = m + 1
                            
                            ! step 1.
                            do l = 1, nupbnd
                                kpx = merge(merge(twwmat(l,2)+kx,twwmat(l,2),5<=n.and.n<=6), 0, t<14) ! 3.16.2017 The 2nd variable in twwmat represents "kpx."
                                ypx = twwmat(l,3) ! ypx
                                opx = twwmat(l,4) ! opx                               
                                twwvec(l) = cww( ax, hx, kx, zx, yx, kpx, ypx, opx, t) ! cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t)   
                            enddo
                            write(str2,fmt='(i3)') nupbnd
                            !write(unit=109,fmt='(a,'//str2//'i4)'), ' twwvec ', twwvec   
                            
                            ! step 2. The essence of this subroutine
                            if(all(twwvec==-99)) cvv(ax,hx,kx,zx,yx,t) = -99 ! 4.1.2017 
                        
                        enddo
                    enddo
     
                    !! ## Block for bginning period mass ! 8-26-2017 match between files generated at different time 
                    !if(printout2)then
                    !    write(6,'(a)') ' Validity in the beginning of periods '
                    !    write(unit=081,fmt='((4x,"iterar",2x),(2x,"iteragov",2x),(2x,"iteratot",2x),(8x,"kx",2x),(8x,"zx",2x),(8x,"yx",2x),(9x,"t",2x))')  
                    !    write(unit=081,fmt='(7(i10,2x))') iterar, iteragov, iteratot, kx, zx, yx, t
                    !    write(unit=081,fmt='(" no.",'//str1//'(x,i3))') (i,i=1,hdim)                     
                    !    do m = 1, fnadim
                    !        !write(unit=081,fmt='(i4,'//str1//'(x,f8.4))')   m, (nsav(m,i,kx,zx,yx,t),i=1,hdim)        
                    !        write(unit=081,fmt='(i4,'//str1//'(x,i3))')   m, (cvv(m,i,kx,zx,yx,t),i=1,hdim)        
                    !    enddo ! m
                    !endif ! printout2
                    
                enddo ! yxi
                deallocate( twwvec, twwmat ) ! nafvvec
            enddo ! n loop
        enddo ! t 1-14  
        
        !write(*,*) ' I am here. Done! '
    end subroutine valid_beginning_period_mass 
       
    ! 3.15.2017 checked. verified how "s2c" is generated.
    subroutine convert_2d_outcome_into_series(opt) ! sza, szh must be consistent with the option tag opt.
        implicit none
        character(len=*), intent(in) :: opt
        integer :: m, i, a, h, k, z, y, kp, yp, op, t, idx
        select case(opt)
            case('coarse')    
                m = size(s3c(:,1))
                do i = 1, m
                    a   = s3c(i,1)
                    h   = s3c(i,2)
                    k   = s3c(i,3)
                    z   = s3c(i,4)
                    y   = s3c(i,5)
                    kp  = s3c(i,6)
                    yp  = s3c(i,7)
                    op  = s3c(i,8)
                    t   = s3c(i,9)
                    idx = s3c(i,10)
                    sww(idx) = cww(a,h,k,z,y,kp,yp,op,t)
                    swf(idx) = cwf(a,h,k,z,y,kp,yp,op,t)
                    swa(idx) = cwa(a,h,k,z,y,kp,yp,op,t)
                    swh(idx) = cwh(a,h,k,z,y,kp,yp,op,t)
                    swk(idx) = cwk(a,h,k,z,y,kp,yp,op,t)
                    swc(idx) = cwc(a,h,k,z,y,kp,yp,op,t)
                enddo
                
                ! 8-25-2017 The accuracy is high to the right of decimal point at least 15 digits.
                if(printout2)then
                    call ssi(sww,'sww')
                    call ss(swf,'swf',20,8)
                    call ss(swa,'swa',20,8)
                    call ss(swh,'swh',20,8)
                    call ssi(swk,'swk')
                    call ss(swc,'swc',20,8)
                    !write(6,'(a)') ' finish printing out sww, swf, swa, swh, swk and swc '
                endif ! 8-25-2017 stop here.
                
            case('test') ! used to see if the backup works normally. 3.15.2017
                m = size(s3c(:,1))
                do i = 1, m
                    a   = s3c(i,1)
                    h   = s3c(i,2)
                    k   = s3c(i,3)
                    z   = s3c(i,4)
                    y   = s3c(i,5)
                    kp  = s3c(i,6)
                    yp  = s3c(i,7)
                    op  = s3c(i,8)
                    t   = s3c(i,9)
                    idx = s3c(i,10)
                    sww2(idx) = cww2(a,h,k,z,y,kp,yp,op,t)
                    swa2(idx) = cwa2(a,h,k,z,y,kp,yp,op,t)
                    swf2(idx) = cwf2(a,h,k,z,y,kp,yp,op,t)
                    swh2(idx) = cwh2(a,h,k,z,y,kp,yp,op,t)
                    swk2(idx) = cwk2(a,h,k,z,y,kp,yp,op,t)
                    swc2(idx) = cwc2(a,h,k,z,y,kp,yp,op,t)
                enddo
                !call ss(swa2,'swa2',20,7)
                !call ss(swf2,'swf2',20,7)
                !call ss(swh2,'swh2',20,7)
                !call ss(swc2,'swc2',20,7)
                !call ssi(sww2,'sww2')
                !call ssi(swk2,'swk2')
            case('distribution-stage-test')
                m = size(s3c(:,1))
                do i = 1, m
                    a   = s3c(i,1)
                    h   = s3c(i,2)
                    k   = s3c(i,3)
                    z   = s3c(i,4)
                    y   = s3c(i,5)
                    kp  = s3c(i,6)
                    yp  = s3c(i,7)
                    op  = s3c(i,8)
                    t   = s3c(i,9)
                    idx = s3c(i,10)
                    sww2(idx) = cww(a,h,k,z,y,kp,yp,op,t)
                    swa2(idx) = cwa(a,h,k,z,y,kp,yp,op,t)
                    swf2(idx) = cwf(a,h,k,z,y,kp,yp,op,t)
                    swh2(idx) = cwh(a,h,k,z,y,kp,yp,op,t)
                    swk2(idx) = cwk(a,h,k,z,y,kp,yp,op,t)
                    swc2(idx) = cwc(a,h,k,z,y,kp,yp,op,t)
                enddo
                !call  ss(swa2,'swa2',20,7)
                !call ssi(sww2,'sww2')
                !call  ss(swf2,'swf2',20,7)
                !call  ss(swh2,'swh2',20,7)
                !call  ss(swc2,'swc2',20,7)
                !call ssi(swk2,'swk2')                
        end select
    end subroutine convert_2d_outcome_into_series
    
    subroutine convert_1d_file_into_matrix(opt)
        implicit none
        character(len=*), intent(in) :: opt   
        integer :: sza, szh
        select case(opt)
            case('coarse')
                sza = adim
                szh = hdim
                call read_series2mat(cwa2,'swa.txt',sza,szh)
                call read_series2mat(cwf2,'swf.txt',sza,szh)
                call read_series2mat(cwh2,'swh.txt',sza,szh)
                call read_series2mat(cwc2,'swc.txt',sza,szh)
                call read_series2mat_i(cww2,'sww.txt',sza,szh)
                call read_series2mat_i(cwk2,'swk.txt',sza,szh)
                
            case('distribution-stage')
                sza = adim
                szh = hdim
                ! FOR TESTING WETHER THIS FUNCTION WORKS NORMALLY, COMMEND OUT THIS BLOCK, AND REMOVE EXCLAMATION OF THE ABOVE BLOCK.
                call read_series2mat(cwa,'swa.txt',sza,szh)
                call read_series2mat(cwf,'swf.txt',sza,szh)
                call read_series2mat(cwh,'swh.txt',sza,szh)
                call read_series2mat(cwc,'swc.txt',sza,szh)
                call read_series2mat_i(cww,'sww.txt',sza,szh)
                call read_series2mat_i(cwk,'swk.txt',sza,szh)   
                
            end select
    contains
        subroutine read_series2mat(mat,fname,sza,szh)
            implicit none
            integer :: sza, szh
            real(wp), dimension(sza,szh,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), intent(out) :: mat
            character(len=*), intent(in) :: fname
            integer :: iostat, m, a, h, k, z, y, kp, yp, op, t, idx
            real(wp) :: val
            m = 0
            open(unit=17,file=fname,status='old',action='read',iostat=iostat) 
                if(iostat==0)then
                    do
                        read(unit=17,fmt=*,iostat=iostat) val
                        !write(*,fmt='(a,i10)') ' I am here ', iostat
                        if(iostat/=0) exit
                        m = m + 1
                        a   = s3c(m,1)
                        h   = s3c(m,2)
                        k   = s3c(m,3)
                        z   = s3c(m,4)
                        y   = s3c(m,5)
                        kp  = s3c(m,6)
                        yp  = s3c(m,7)
                        op  = s3c(m,8)
                        t   = s3c(m,9)
                        mat(a,h,k,z,y,kp,yp,op,t) = val
                        !write(*,fmt='(f15.7)') val
                    enddo
                    if(m/=sza*szh*1018) write(*,fmt='(a)') ' something wrong when reading serial files '
                else
                    write(*,fmt='(a)') ' something wrong with read_series2mat '
                endif
            close(10)
        end subroutine read_series2mat 
        
        subroutine read_series2mat_i(mat,fname,sza,szh)
            implicit none
            integer :: sza, szh
            integer, dimension(sza,szh,0:(kdim-1),0:1,0:nmc,0:(kdim-1),0:nmc,0:2,1:14), intent(out) :: mat
            character(len=*), intent(in) :: fname
            integer :: iostat, m, a, h, k, z, y, kp, yp, op, t, idx
            integer:: val
            m = 0
            open(unit=17,file=fname,status='old',action='read',iostat=iostat) 
                if(iostat==0)then
                    do
                        read(unit=17,fmt=*,iostat=iostat) val
                        !write(*,fmt='(a,i10)') ' I am here ', iostat
                        if(iostat/=0) exit
                        m = m + 1
                        a   = s3c(m,1)
                        h   = s3c(m,2)
                        k   = s3c(m,3)
                        z   = s3c(m,4)
                        y   = s3c(m,5)
                        kp  = s3c(m,6)
                        yp  = s3c(m,7)
                        op  = s3c(m,8)
                        t   = s3c(m,9)
                        mat(a,h,k,z,y,kp,yp,op,t) = val
                        !write(*,fmt='(f15.7)') val
                    enddo
                    if(m/=sza*szh*1018) write(*,fmt='(a)') ' something wrong when reading serial files '
                else
                    write(*,fmt='(a)') ' something wrong with read_series2mat '
                endif
            close(10)
        end subroutine read_series2mat_i         
    end subroutine convert_1d_file_into_matrix
    
    ! crnmat(.,1): coordinate x (a); crnmat(.,2): coordinate y (h); crnmat(.,3): signal; crnmat(.,4): weight ! 4.14.2017
    subroutine create_four_corners(ap,hp,rav,rhv,crnmat,vmat,indicator)
        implicit none
        real(wp), intent(in) :: ap ,hp
        real(wp), dimension(:), intent(in) :: rav, rhv
        integer,  dimension(:,:), intent(in) :: vmat
        real(wp), dimension(4,4), intent(out) :: crnmat
        integer :: adx, hdx, i
        integer, intent(out) :: indicator
        real(wp) :: ax, bx, ay, by, sum1
        indicator = 0
        adx = locate(rav,ap)
        hdx = locate(rhv,hp)
        !print*, adx,hdx
        if(adx+1<=adim.and.hdx+1<=hdim)then
            crnmat(1,1) = real(adx  ,wp) ! NW 
            crnmat(2,1) = real(adx+1,wp) ! SW
            crnmat(3,1) = real(adx  ,wp) ! NE
            crnmat(4,1) = real(adx+1,wp) ! SE
            crnmat(1,2) = real(hdx  ,wp) ! NW
            crnmat(2,2) = real(hdx  ,wp) ! SW
            crnmat(3,2) = real(hdx+1,wp) ! NE
            crnmat(4,2) = real(hdx+1,wp) ! SE
            do i = 1, 4
                crnmat(i,3) = vmat(int(crnmat(i,1)),int(crnmat(i,2))) ! either 1 or -99    
            enddo
            
            ax = (rav(adx+1)-ap)/(rav(adx+1)-rav(adx))
            bx = (ap-rav(adx))/(rav(adx+1)-rav(adx))
            ay = (rhv(hdx+1)-hp)/(rhv(hdx+1)-rhv(hdx))
            by = (hp-rhv(hdx))/(rhv(hdx+1)-rhv(hdx))
            
            crnmat(1,4) = merge(ax*ay,0._wp,crnmat(1,3)==-99._wp) ! <----- the integer to real number conversion may be problematic.
            crnmat(2,4) = merge(bx*ay,0._wp,crnmat(2,3)==-99._wp) ! <----- the integer to real number conversion may be problematic.
            crnmat(3,4) = merge(ax*by,0._wp,crnmat(3,3)==-99._wp) ! <----- the integer to real number conversion may be problematic.
            crnmat(4,4) = merge(bx*by,0._wp,crnmat(4,3)==-99._wp) ! <----- the integer to real number conversion may be problematic.
            sum1 = sum(crnmat(:,4))
            crnmat(:,4) = crnmat(:,4)/sum1      
            
        else
            indicator = 1
        endif
            
    end subroutine create_four_corners
    
    subroutine adjust_corners(hp,crnmat,vmat,erridx)
        implicit none
        real(wp), intent(inout) :: crnmat(4,4), hp
        integer,  dimension(:,:), intent(in) :: vmat
        integer :: i
        integer, intent(out) :: erridx
        real(wp) :: tp1, tp2, wgt1, wgt2
        logical :: flag1
        flag1 = .true.
        erridx = 0
        do while(flag1)
            crnmat(1,1) = real(int(crnmat(1,1)) + 1, wp)
            crnmat(2,1) = real(int(crnmat(2,1)) + 1, wp)
            crnmat(3,1) = real(int(crnmat(3,1)) + 1, wp)
            crnmat(4,1) = real(int(crnmat(4,1)) + 1, wp)
            if((int(crnmat(4,1)) + 1)>adim)then ! 10132016
                erridx = 1
                exit
            endif
            tp1 = vmat(crnmat(2,1),crnmat(2,2)) ! SW
            tp2 = vmat(crnmat(4,1),crnmat(4,2)) ! SE
            if(tp1==-99._wp.and.tp2==-99._wp)then
                crnmat(:,3) = 1._wp ! reset (flag)
                crnmat(:,4) = 0._wp ! reset (weight)
                crnmat(2,3) = -99._wp ! signal for valid point
                crnmat(4,3) = -99._wp
                wgt1 = (rhv(crnmat(4,2))-hp)/(rhv(crnmat(4,2))-rhv(crnmat(2,2)))
                wgt2 = (hp-rhv(crnmat(2,2)))/(rhv(crnmat(4,2))-rhv(crnmat(2,2)))
                crnmat(2,4) = wgt1 ! new weight
                crnmat(4,4) = wgt2 ! new weight
                flag1 = .false.
            elseif(tp1==-99._wp)then
                crnmat(:,3) = 1._wp
                crnmat(:,4) = 0._wp
                crnmat(2,3) = -99._wp
                crnmat(2,4) = 1._wp
                flag1 = .false.    
            elseif(tp2==-99._wp)then
                crnmat(:,3) = 1._wp
                crnmat(:,4) = 0._wp
                crnmat(4,3) = -99._wp
                crnmat(4,4) = 1._wp 
                flag1 = .false.
            endif
        enddo
    end subroutine adjust_corners
    
    subroutine convert_2d_distribution_into_1d()
        implicit none
        integer :: m, i, a, h, k, z, y, kp, yp, op, t, idx
        !integer :: tstart, tend, trate, tmax
        m = size(s3c(:,1)) ! the length of the list s2c.
        !call system_clock(tstart,trate,tmax)
        !call system_clock(tstart)        
        ! 3.25.2017 parallel remains faster than sequential one. On average, 
        !$omp parallel do default(shared) private(a,h,k,z,y,kp,yp,op,t,idx)
        do i = 1, m
            a   = s3c(i,1)
            h   = s3c(i,2)
            k   = s3c(i,3)
            z   = s3c(i,4)
            y   = s3c(i,5)
            kp  = s3c(i,6)
            yp  = s3c(i,7)
            op  = s3c(i,8)
            t   = s3c(i,9)
            idx = s3c(i,10)
            sef(idx) = cef(a,h,k,z,y,kp,yp,op,t) ! 4.16.2017 `cef` is reset in subroutine mass_transition.
        enddo        
        !$omp end parallel do
        !call system_clock(tend)
        !write(*,fmt='(a,f12.4,a,f12.9)') ' 2d-1d time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds', sum(sef(1:5000))          
        !call ss(sef,'sef',25,12)
        !call ss_science(sef,'sef_e',20,13) ! verified. 9-11-2017 
    end subroutine convert_2d_distribution_into_1d
    
    subroutine convert_2d_macro_statistics_into_1d()
        implicit none
        !integer :: m, i, a, h, k, z, y, kp, yp, op, t, idx
        !m = size(s2c(:,1))
        !!$omp parallel do default(shared) private(a,h,k,z,y,kp,yp,op,t,idx)
        !do i = 1, m
        !    a   = s2c(i,1)
        !    h   = s2c(i,2)
        !    k   = s2c(i,3)
        !    z   = s2c(i,4)
        !    y   = s2c(i,5)
        !    kp  = s2c(i,6)
        !    yp  = s2c(i,7)
        !    op  = s2c(i,8)
        !    t   = s2c(i,9)
        !    idx = s2c(i,10)
        !    sef(idx) = cef(a,h,k,z,y,kp,yp,op,t)
        !enddo        
        !!$omp end parallel do
        call ss(sw_taxableincome,'sw_taxableincome',25,12)
        call ss(sw_aftertaxwealth,'sw_aftertaxwealth',25,12)
        call ss(swc,'sw_consumption',25,12)
        
    end subroutine convert_2d_macro_statistics_into_1d    
    
    subroutine mass_transition(exit_log1,error,idxr,idxg,idxt,trial_id) ! 3.16.2017 4:08 pm stop here. 3.17.2017 3:13 pm begin here again.
        implicit none
        
        integer :: sz, sz2, q, n, tdx, kd, kdold, kpxn, kpxni, zpxn, ypxn, i, opxn, idx, axi, xv(2), j, sumerr, erridx
        real(wp) :: sum0, sum1, sum2, sum3, sum4, crnmat(4,4), beqsum, wv(2), biztrs, pz2prob
        logical :: flag1
        integer, intent(in) :: idxr, idxg, idxt, trial_id
        integer :: tstart, tend, trate, tmax
        logical, intent(inout) :: exit_log1
        real(wp), intent(out) :: error
        integer :: sn1
        !real(wp) :: sum_accurate
        ! debug remember to remove the variables in the following block
        real(wp) :: sum5
        integer :: di, locmax(9)
        character(len=3) :: str1, str_idxr, str_idxg, str_idxt
        
        call system_clock(tstart,trate,tmax) 
        exit_log1 = .false.
        sz  = size(s3c(:,1))
        cef = 0._wp ! matrix initialization. 4.16.2017 Correct. Needs initialization here.
        sumerr = 0 ! error variable initialization. An integer.
        error = 100000._wp
        
        if(inv_dist_counter==1)then ! 8-24-2017 For the first round (inv_dist_counter==1), we assume "Uniform Distribution" in the beginning of the first period, and in the block below calculate the corresponding mass distirbution in the end of the first period.
            
            allocate(ivec(sz))
            ivec = .false.
            ivec = s3c(:,9)==1 ! logical. length of period one.
            
            szperiod1 = count(ivec) !<----------------- important. Length of period one. macro variable.
            allocate(nvec(szperiod1))
            nvec = pack(s3c(:,10),ivec) ! 3.18.2017 shrink full vector to an appropriate subset. Trick: the size of the destination ("nvec") of packing must be explicit and exact.

            ! 8-24-2017 We only allow people with {financial asset} = transbeq and {non-financial asset} = the minimum housing services have non-zero mass.
            call linwgt(rav,transbeq,xv,wv) ! 3.17.2017 Is it redundant? No, it is used in the block right below. 11:13 am 3.17.2017 Tempararily stop here.
            
            ! 3.18.2017 rav = av, checked.
            
            !! Debug 8-24-2017
            !sum_accurate = 0._wp ! 8-24-2017
            !scef = 0._wp ! 8-24-2017
            
            ! 8-24-2017 NOTE!!! Don't use openmp here, because subsequent revision includes summation (of mass in a more accurate way) without taking care of the side effect of parallelization.
            !!$omp parallel do default(shared) private(n,t,biztrs) ! 4.12.2017 Too much overhead involoves if we use parallelism. Turn off openmp here.
            do q = 1, szperiod1
                n   = nvec(q) ! nvec is the index list.
                ax  = s3c(n,1) ! a
                hx  = s3c(n,2) ! h
                kx  = s3c(n,3) ! k
                zx  = s3c(n,4) ! z
                yx  = s3c(n,5) ! y
                kpx = s3c(n,6) ! kp
                ypx = s3c(n,7) ! yp
                opx = s3c(n,8) ! op
                t   = s3c(n,9) ! t
                
                ! 3.17.2017 Initialization of mass distribution (In the current version, I assume people enter the model without any endowment and live on the minimum housing service.)
                ! 3.17.2017 I assume at the end of the first period, people with their own endowed human capital and the lump sum bequest choose their the end-of-period decision of the first peirod.
                if((ax==xv(1).or.ax==xv(2)).and.(hx==1))then ! 3.17.2017 revised this block for distribution based on linear interpolation
                    ! if(cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t)/=-99)then ! 3.18.2017 removed
                    if(cvv(ax,hx,kx,zx,yx,t)/=-99)then ! 3.18.2017 added
                        !write(unit=127,fmt='(9i4,2f8.4)') ax,hx,kx,zx,yx,kpx,ypx,opx,t,biztrs,cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)
                        cycle
                    else
                        if(ax==xv(1))then ! 4.10.2017 stop here. 
                            biztrs = merge( merge(pka(0,1),pka(0,2),opx==0), merge(pka(kx,1),pka(kx,2),opx==1), zx==0) ! 4.1.2017 Consistent with variable_space_v3.xlsx 4.11.2017 correctly following the lottery probability.
                            !!$omp critical
                            cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = wv(1)*py(yx,ypx)*biztrs
                            !!$omp end critical
                            !write(unit=127,fmt='(9i4,2f8.4)') ax,hx,kx,zx,yx,kpx,ypx,opx,t,biztrs,cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)
                            
                            !! Debug ! 8-24-2017
                            !scef(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t) ! 8-24-2017
                            !sum_accurate = sum_accurate + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t) ! 8-24-2017
                        else ! ax==xv(2)
                            biztrs = merge( merge(pka(0,1),pka(0,2),opx==0), merge(pka(kx,1),pka(kx,2),opx==1), zx==0) ! 4.1.2017 consistent with variable_space_v3.xlsx 4.11.2017 correctly following the lottery probability.
                            !!$omp critical
                            cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = wv(2)*py(yx,ypx)*biztrs   
                            !!$omp end critical
                            !write(unit=127,fmt='(9i4,2f8.4)') ax,hx,kx,zx,yx,kpx,ypx,opx,t,biztrs,cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)
                            
                            !! Debug ! 8-24-2017 
                            !scef(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t) ! 8-24-2017
                            !sum_accurate = sum_accurate + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t) ! 8-24-2017
                        endif
                    endif
                    !cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = 1._wp ! comment out this line for distribution mass according to the linear interpolation result.
                else
                    cycle    
                endif
            enddo
            !!$omp end parallel do
            
            call system_clock(tend) 
            
            sum1 = sum(cef(:,:,:,:,:,:,:,:,1))
            
            cef(:,:,:,:,:,:,:,:,1) = cef(:,:,:,:,:,:,:,:,1)/sum1*popfrac(1)
            
            !! 8-24-2017 Checked. Conclusion: sum1 equals to sum_accurate, the later uses traditional way for summation; the former uses the intrisic funtion "sum".
            !write(*,'(" Mama, I am here. ")') 
            !write(*,'(3(a,f20.15,2x))') 'sum1  =', sum1, 'agg_sum1  =', sum(cef(:,:,:,:,:,:,:,:,1)), 'pop(1)=', popfrac(1) ! 4.1.2017
            !write(*,'(2(a,f20.15,2x),/)') 'sumacc=', sum_accurate, 'agg_sumacc=', sum(scef(:,:,:,:,:,:,:,:,1)) ! 8-24-2017
            !if(printout6) write(*,fmt='(a,i3,a,f12.4,a)') 'initial distribution ',t,', time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'     
            
            !! 8-26-2017 to print out the mass sequence in the similar way as convert_2d_distribution_into_1d
            !if(printout2)then
            !    write(*,'(a)') ' initial mass distribution (uniform distribution) '
            !    do i = 1, size(s3c(:,1))
            !        write(unit=134,fmt='("no.",i10,2x,"mass=",e20.15)') s3c(i,10), cef(s3c(i,1),s3c(i,2),s3c(i,3),s3c(i,4),s3c(i,5),s3c(i,6),s3c(i,7),s3c(i,8),s3c(i,9))       
            !    enddo ! i
            !endif ! printout2
            
            deallocate( ivec, nvec )
            
        else
            ! 8-24-2017 sef1 is obtained in subroutine intergenerational transfer.
            sef1 = sef1/sum(sef1)*popfrac(1) ! NORMALIZATION 10042016 (it is initialized in distribution loop of equilibrium.f90) ! 4.14.2017 END-of-1st period. correct.
            
            ! ADJUSTMENT FOR INVALID POINT WITH POSITIVE MASS BY MOVING THE MASS TO NEXT VALID POINT IN ITS NEIGHBORHOOD. 10102016
            cef = 0._wp ! 4.22.2017 Initialization 8-24-2017
            ! 4.12.2017 sequence converted to matrix form. Lower part is the minor adjustment for invalid combination with mass (seems redundant, because intergenerational_transition screen out invalid transition.)
            do q = 1, szperiod1
                !n   = nvec(q)
                ax  = s3c(q,1) ! using `q` rather than the `n` in last block can save our trouble in allocating `nvec` and `ivec`
                hx  = s3c(q,2)
                kx  = s3c(q,3)
                zx  = s3c(q,4)
                yx  = s3c(q,5)
                kpx = s3c(q,6)
                ypx = s3c(q,7)
                opx = s3c(q,8)
                t   = s3c(q,9)
                idx = s3c(q,10) ! 3.17.2017 replace the line below with this one.
                !idx = c3s(ax,hx,kx,zx,yx,kpx,ypx,opx,t)
                
                ! 4.14.2017 This if-statement seems redundant, because again in subroutien intergenerational transition invalid points of 1st period with positive mass has been screened out and mass has been reallocated.
                if(cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t)/=-99.and.sef1(idx)>0._wp)then ! 4.14.2017 If invalid points of 1st period has mass, shift it to the next financial asset level.
                    flag1 = .true.
                    axi   = ax
                    do while(flag1)
                        axi = axi + 1
                        if(axi>adim)then
                            print*, ' something wrong in model line 3284 '    
                        else
                            if(cww(axi,hx,kx,zx,yx,kpx,ypx,opx,t)==-99)then
                                cef(axi,hx,kx,zx,yx,kpx,ypx,opx,t) = sef1(idx)
                                flag1 = .false. ! 4.1.2017 stop message sent out.
                            endif
                        endif
                    enddo ! flag1
                else
                    cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t) = sef1(idx) ! 4.14.2017 This is the conversion from sequence to structure we need.
                endif
            enddo ! q   
            !write(*,*) ' 1st end-of-period mass: mama ', sum(cef(:,:,:,:,:,:,:,:,1)), popfrac(1) ! 4.1.2017
        endif ! inv_dist_counter == 1
        
        allocate(ivec(sz))
        ! MASS TRANSITION BLOCK--------------------------------------------------- BEGINNING
        do tdx = 1, 13 ! current period
        !do tdx = 1, 1
            ivec = .false. ! initialization
            ivec = s3c(:,9)==tdx ! get true or false sequence. 
            allocate( nvec(count(ivec==.true.)))
            nvec = pack(s3c(:,10),ivec) ! 3.17.2017 packing indices of combination series
            sum2 = 0._wp ! both of which are "shared."
            sum3 = 0._wp ! both of which are "shared."
            sum5 = 0._wp ! debug
            sz2  = count(ivec==.true.) ! All the indices in the current period.
            
            !$omp parallel default(shared) private(n,t,sn1,sum4,zpxn,crnmat,i,kpxn,opxn,ypxn,kpxni,kd,kdold,erridx,pz2prob) reduction(+:sumerr) !## 1 ! 4.12.2017
            !$omp do !## 2
            do q = 1, sz2 ! 3.17.2017 "sz2" is the number of all combinations (regardless of its validity).
                ! the combination of state vairables denotes the mass to be moved down toward the next period. 3.17.2017
                n   = nvec(q) ! it is required to get the index list.
                ax  = s3c(n,1)
                hx  = s3c(n,2)
                kx  = s3c(n,3)
                zx  = s3c(n,4)
                yx  = s3c(n,5)
                kpx = s3c(n,6)
                ypx = s3c(n,7)
                opx = s3c(n,8)
                t   = s3c(n,9)
                sn1 = s3c(n,10) ! index 4.12.2017 It seems redundant, because `n` has the same definition.
                
                !! OPENMP EXPERIMENTS. WORKS GREAT. 10042016
                !!$omp atomic update
                !sum2 = sum2 + 1
                !!$omp critical
                !write(unit=122,fmt='(i6)') sn1
                !!$omp end critical
                
                !if(cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t)/=-99.and.cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)>0._wp)then
                !    write(unit=121,fmt='(9a)') '  ax', '  hx', '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                !    write(unit=121,fmt='(9i4)') ax,hx,kx,zx,yx,kpx,ypx,opx,t  
                !endif
                
                !if(t==2)then
                !    !write(unit=121,fmt='(9a)') '  ax', '  hx', '  kx','  zx','  yx',' kpx',' ypx',' opx','   t'
                !    write(unit=121,fmt='(9i4)') ax,hx,kx,zx,yx,kpx,ypx,opx,t  
                !endif                
                
                !$omp atomic update !## 3
                sum2 = sum2 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t) ! 3.17.2017 add up mass regardless of whether the combination is valid or not.
                
                !write(unit=122,fmt='(9a)') '  ax', '  hx', '  kx', '  zx', '  yx', ' kpx', ' ypx', ' opx', '   t'
                !write(unit=122,fmt='(9i4,f20.15)') ax,hx,kx,zx,yx,kpx,ypx,opx,t,sum2                  
                
                ! 3.20.2017 These seem `useless` to impose these two conditions. Just keep it.
                if(cww(ax,hx,kx,zx,yx,kpx,ypx,opx,t)/=-99)then
                    if(cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)>0._wp) write(unit=131,fmt='(a,x,9i4)'), 'cww/=-99 but cef>0: ', ax, hx, kx, zx, yx, kpx, ypx, opx, t
                    cycle ! 3.17.2017 skip if the combination is not valid.
                endif
                
                !if(cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)==0._wp)then
                if(abs(cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t))<1.e-10_wp)then
                    cycle ! 3.17.2017 skip if the combination has no mass at all.
                endif
                
                !sum5 = sum5 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t) ! debug
                
                ap = cwa(ax,hx,kx,zx,yx,kpx,ypx,opx,t) ! beginning of next period's financial ! the target needs to be interpolated.
                hp = cwh(ax,hx,kx,zx,yx,kpx,ypx,opx,t) ! beginning of next period's housing   ! the target needs to be interpolated.
                kd = cwk(ax,hx,kx,zx,yx,kpx,ypx,opx,t) ! beginning of next period's career    ! the target needs to be interpolated.        
                
                !write(unit=131,fmt='(9(x,a))') ' ax',' hx',' kx',' zx',' yx','kpx','ypx','opx','  t'
                !write(unit=131,fmt='(9i4)') ax, hx, kx, zx, yx, kpx, ypx, opx, t
                
                ! crnmat(.,1): coordinate x (a); crnmat(.,2): coordinate y (h); crnmat(.,3): signal; crnmat(.,4): weight
                if(t>=1.and.t<=7)then
                    sum4 = 0._wp
                    di   = 1 ! debug
                    do zpxn = 0, 1 ! 4.12.2017 zpxn == 0 cases results in the same situation a worker faces.
                        if(zpxn==0)then ! 4.12.2017 Uncertainty generator for the current period.
                            
                            call create_four_corners(ap,hp,rav,rhv,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx) ! interpolation: output is `crnmat`; `cvv` is inoput for signalling bad points
                            if(erridx==1)then ! 4.12.2017 It's very unlikely to trigger this condition.
                                sumerr = sumerr + erridx
                                !print*, 'warning -------------------------- 1'
                                cycle
                            endif
                            ! 4.12.2017
                            if(all(crnmat(:,3)==1._wp)) call adjust_corners(hp,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx) 
                            if(erridx==1)then
                                sumerr = sumerr + erridx
                                !print*, 'warning -------------------------- 2'
                                cycle
                            endif       
                            
                            do i = 1, 4 ! loop over each interpolated point
                                
                                if(crnmat(i,3)/=-99._wp)then
                                    !write(unit=131,fmt='("<",i4,">")') di
                                    !di = di + 1
                                    cycle
                                endif
                                
                                do kpxn = 0, 1 ! 4.12.2017 Uncertainty generator. either keep being a worker or turning to be an entrepreneur. <---- stop here.
                                    !kpxni = merge(1,2,kpxn==0)
                                    opxn = merge(2,0,kpxn==1)
                                    do ypxn = 1, nmc ! Uncertainty generator.
                                        ! 3.19.2017 nature happens: kpxn (opxn) and ypxn (both are controlled in this do-loops); innate: ypx and 
                                        !cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) + &
                                        !    & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pka(kd,kpxn+1)*py(ypx,ypxn)*survprob(t) ! three uncertainties: ideas, labor efficiency and survival rate.  
                                        
                                        pz2prob = merge(1._wp,pz2(kd,zpxn+1),kd==0) ! 4.12.2017 There is no business shock to wage earners. Correct.
                                        !$omp atomic update !## 4
                                        cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) + &
                                            & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2prob*pka(0,kpxn+1)*py(ypx,ypxn)*survprob(t) ! 4.12.2017 replace `kd` with `0` because people hit by bad luck they can only be a worker or run the smallest business project next period.                                         
                                        sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pka(0,kpxn+1)*py(ypx,ypxn)*survprob(t)*pz2prob   
                                        
                                        !write(unit=131,fmt='("<",i4,">",2e15.8)') di, cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pka(0,kpxn+1)*py(ypx,ypxn)*survprob(t)*pz2prob, sum4   
                                        !di = di + 1                                        
                                        
                                        !if(cww(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1)/=-99.and.crnmat(i,4)>0._wp) &
                                        !    & write(unit=121,fmt='(9i4,a,i3,f8.4,f8.4,4f8.4,4f8.4)'), int(crnmat(i,1)),int(crnmat(i,2)),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1,' something wrong ', i, crnmat(i,3), crnmat(i,4), crnmat(:,3), crnmat(:,4)
                                    enddo ! ypxn       
                                enddo ! kpxn
                            enddo ! i
                            
                        else ! zpxn==1 ! 4.12.2017 Stop here 11:23 am. copy ln.2811.
                            
                            if(kd/=0)then ! 4.14.2017 actually this step is redundant.
                                
                                call create_four_corners(ap,hp,rav,rhv,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)
                                if(erridx==1)then
                                    sumerr = sumerr + erridx
                                    !sum5 = sum5 + 1._wp
                                    !print*, 'warning -------------------------- 3'
                                    cycle
                                endif
                                if(all(crnmat(:,3)==1._wp))then
                                    call adjust_corners(hp,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)
                                    !write(unit=127,fmt='(a,i6,9i4,4f12.5," zpxn ",i4)') 'norm ',sn1,ax,hx,kx,zx,yx,kpx,ypx,opx,t,crnmat(:,3),zpxn
                                endif  
                                if(erridx==1)then
                                    sumerr = sumerr + erridx
                                    !sum5 = sum5 + 1._wp
                                    !write(unit=127,fmt='(a,i6,9i4,4f12.5," zpxn ",i4)') 'norm ',sn1,ax,hx,kx,zx,yx,kpx,ypx,opx,t,crnmat(:,3),zpxn
                                    !print*, 'warning -------------------------- 4'
                                    cycle
                                endif    
                                
                                do i = 1, 4
                                    if(crnmat(i,3)/=-99._wp) cycle
                                    !if(zpxn==1)then ! good shocks (I)
                                    
                                    if(kd/=kdim-1)then ! kd == 1 or 2
                                        do kpxn = 0, 1 ! evolution of business ideas (II)
                                            kpxni = merge(1,2,kpxn==0) ! indicates the column number corresponds to the probability of new or existing business idea.
                                            opxn  = merge(1,2,kpxn==0) ! indicates whether stay put or engage in advanced project.
                                            do ypxn = 1, nmc ! (III)
                                                !$omp atomic update !## 5
                                                cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn+kd,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn+kd,ypxn,opxn,t+1) + &
                                                    & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*py(ypx,ypxn)*survprob(t)  
                                                sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*py(ypx,ypxn)*survprob(t)  
                                                
                                                !write(unit=131,fmt='("<",i4,">",2e15.8)') di, cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*py(ypx,ypxn)*survprob(t), sum4   
                                                !di = di + 1                                                                                        
                                                
                                            enddo
                                        enddo ! kpxn
                                    else ! kd==kdim-1 (that is, the largest business project)
                                        do kpxn = 3, 3
                                            kpxni = 1 ! first column (II)
                                            opxn  = 1 ! stay put
                                            do ypxn = 1, nmc ! the evolution of next perk=od's labor efficiency (III)
                                                !$omp atomic update !## 6
                                                cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) + &
                                                    & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*py(ypx,ypxn)*survprob(t)      
                                                sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*py(ypx,ypxn)*survprob(t)      
                                                
                                                !write(unit=131,fmt='("<",i4,">",2e15.8)') di, cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*py(ypx,ypxn)*survprob(t), sum4   
                                                !di = di + 1                                                                                        
                                                
                                            enddo
                                        enddo ! kpxn                                    
                                    endif ! kd
                                enddo ! i
                            endif ! kd/=0
                            
                        endif ! zpxn                         
                    enddo ! zpxn ! 4.12.2017    
                    
                    !$omp atomic update !## 7
                    sum3 = sum3 + sum4 ! sum3 is shared so we uses atomic.
                    !write(unit=122,fmt='(36x,f20.15)') sum3 
                    !if(abs(sum4-cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t))<1.e-10)then
                    !    write(unit=131,fmt='(2e15.8)') sum4, cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)
                    !else
                    !    write(unit=131,fmt='(2e15.8,a)') sum4, cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t), " --------------- -----------------------------------------------"
                    !endif
                    
                    ! 3.19.2017 Stop here 4:30 pm.
                elseif(t==8)then ! 4.12.2017
                    sum4 = 0._wp
                    do zpxn = 0, 1
                        if(zpxn==0)then
                            call create_four_corners(ap,hp,rav,rhv,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)
                            if(erridx==1)then
                                sumerr = sumerr + erridx
                                !print*, 'warning -------------------------- 5'
                                cycle
                            endif
                            if(all(crnmat(:,3)==1._wp)) call adjust_corners(hp,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)
                            if(erridx==1)then
                                sumerr = sumerr + erridx
                                !print*, 'warning -------------------------- 6'
                                cycle
                            endif  
                            do i = 1, 4
                                if(crnmat(i,3)/=-99._wp) cycle
                                do kpxn = 0, 1 ! control the column number of future business idea.
                                    opxn = merge(2,0,kpxn==1) ! Note: in the current block (the worker case), it is either 2 or zero.
                                    ypxn = 0 ! no more labor efficiency realized at the end of period 9.
                                    pz2prob = merge(1._wp,pz2(kd,zpxn+1),kd==0) ! 4.12.2017
                                    !$omp atomic update !## 8
                                    cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) + &
                                        & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2prob*pka(0,kpxn+1)*survprob(t) ! no bability of labor efficiency; no probability of current business shock.
                                    sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pka(0,kpxn+1)*survprob(t)*pz2prob ! no business shock, no labor efficiency.
                                enddo ! kpxn
                            enddo ! i    
                        else ! zpxn==1
                            if(kd/=0)then
                                call create_four_corners(ap,hp,rav,rhv,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)
                                if(erridx==1)then
                                    sumerr = sumerr + erridx
                                    !print*, 'warning -------------------------- 7'
                                    !sum5 = sum5 + 1._wp
                                    cycle
                                endif
                                if(all(crnmat(:,3)==1._wp))then
                                    call adjust_corners(hp,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)
                                    !write(unit=127,fmt='(a,i6,9i4,4f12.5," zpxn ",i4)') 'norm ',sn1,ax,hx,kx,zx,yx,kpx,ypx,opx,t,crnmat(:,3),zpxn
                                endif  
                                if(erridx==1)then
                                    sumerr = sumerr + erridx
                                    !print*, 'warning -------------------------- 8'
                                    !sum5 = sum5 + 1._wp
                                    !write(unit=127,fmt='(a,i6,9i4,4f12.5," zpxn ",i4)') 'norm ',sn1,ax,hx,kx,zx,yx,kpx,ypx,opx,t,crnmat(:,3),zpxn
                                    cycle
                                endif
                                
                                do i = 1, 4
                                    if(crnmat(i,3)/=-99._wp) cycle
                                    if(kd/=kdim-1)then
                                        do kpxn = 0, 1 ! evolution of business ideas (II)
                                            kpxni = merge(1,2,kpxn==0)
                                            opxn  = merge(1,2,kpxn==0)
                                            ypxn = 0 ! no labor efficiency arrived at the end of the 9-th period.
                                            !$omp atomic update !## 9
                                            cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn+kd,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn+kd,ypxn,opxn,t+1) + &
                                                & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*survprob(t) ! no labor efficiency.  
                                            sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*survprob(t)                                                  
                                        enddo ! kpxn
                                    else ! kd==kdim-1 The largest business project.
                                        kpxn = 3
                                        kpxni = 1
                                        opxn  = 1
                                        ypxn  = 0 ! no more labor efficiency realized since the end of the 9-th period.
                                        !$omp atomic update !## 10
                                        cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) + &
                                            & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*survprob(t) ! no labor efficiency shock      
                                        sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*survprob(t) ! no labor efficiency shock.                                  
                                    endif ! kd
                                enddo ! i
                            endif ! kd/=0 or not    
                        endif ! zpxn = 1 or not
                    enddo ! zpxn
                    
                elseif(t>=9.and.t<=12)then
                    sum4 = 0._wp
                    do zpxn = 0, 1
                        if(zpxn==0)then
                            call create_four_corners(ap,hp,rav,rhv,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)
                            if(erridx==1)then
                                sumerr = sumerr + erridx
                                !print*, 'warning -------------------------- 9'
                                cycle
                            endif
                            if(all(crnmat(:,3)==1._wp)) call adjust_corners(hp,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)
                            if(erridx==1)then
                                sumerr = sumerr + erridx
                                !print*, 'warning -------------------------- 10'
                                cycle
                            endif  
                            do i = 1, 4
                                if(crnmat(i,3)/=-99._wp) cycle
                                kpxn = 0 ! once the agent is hit by a bad luck, retiree next period
                                opxn = 0
                                ypxn = 0
                                pz2prob = merge(1._wp,pz2(kd,zpxn+1),kd==0) ! 4.12.2017
                                !$omp atomic update !## 11
                                cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) + &
                                    & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2prob*survprob(t)   
                                sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*survprob(t)*pz2prob
                            enddo ! i                           
                        else ! zpxn==1
                            if(kd/=0)then
                                call create_four_corners(ap,hp,rav,rhv,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)
                                if(erridx==1)then
                                    sumerr = sumerr + erridx
                                    !print*, 'warning -------------------------- 11'
                                    cycle
                                endif
                                if(all(crnmat(:,3)==1._wp)) call adjust_corners(hp,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)   
                                if(erridx==1)then
                                    sumerr = sumerr + erridx
                                    !print*, 'warning -------------------------- 12'
                                    cycle
                                endif  
                                
                                do i = 1, 4
                                    if(crnmat(i,3)/=-99._wp) cycle
                                    if(kd/=kdim-1)then ! not-the-largest business projects
                                        do kpxn = 0, 1 ! kpxn: for advancing the index of a project
                                            kpxni = merge(1,2,kpxn==0) ! kpxni: for matrix index
                                            opxn  = merge(1,2,kpxn==0)
                                            ypxn  = 0
                                            !$omp atomic update !## 12
                                            cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn+kd,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn+kd,ypxn,opxn,t+1) + &
                                                & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*survprob(t) ! 3.20.2017 only business shock and new business idea would realize. 
                                            sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*survprob(t)                                                  
                                        enddo ! kpxn
                                    else ! kd==kdim-1
                                        kpxn = 3
                                        kpxni = 1
                                        opxn  = 1
                                        ypxn  = 0
                                        !$omp atomic update !## 13
                                        cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) + & ! 3.20.2017 stay put
                                            & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*survprob(t)      
                                        sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*pka(kd,kpxni)*survprob(t)                                              
                                    endif ! kd
                                enddo ! i    
                            endif ! kd
                        endif ! zpxn
                    enddo ! zpxn
                    
                    !$omp atomic update !## 14
                    sum3 = sum3 + sum4
                elseif(t==13)then
                    sum4 = 0._wp
                    if(kd==0)then ! retirees. OK 4.14.2017
                        zpxn = 0 ! no business shock any more for retrees.
                        call create_four_corners(ap,hp,rav,rhv,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)    
                        if(erridx==1)then
                            sumerr = sumerr + erridx
                            !print*, 'warning -------------------------- 13'
                            cycle
                        endif
                        if(all(crnmat(:,3)==1._wp)) call adjust_corners(hp,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)
                        if(erridx==1)then
                            sumerr = sumerr + erridx
                            !print*, 'warning -------------------------- 14'
                            cycle
                        endif                             
                        do i = 1, 4
                            if(crnmat(i,3)/=-99._wp) cycle
                            kpxn = 0
                            opxn = 0
                            ypxn = 0
                            !$omp atomic update !## 15
                            cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) + &
                                & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*survprob(t) ! no business shock. dead for sure.   
                            sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*survprob(t)
                        enddo ! i                           
                    else ! kd/=0 Now we are at period 13.
                        do zpxn = 0, 1
                            !kdold = kd
                            !if(zpxn==0) kd=0 ! even in the final period, receiving negative business shock forces entrepreneurs to exit entrepreneurship.
                            call create_four_corners(ap,hp,rav,rhv,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)
                            if(erridx==1)then
                                sumerr = sumerr + erridx
                                !print*, 'warning -------------------------- 15'
                                cycle
                            endif
                            if(all(crnmat(:,3)==1._wp)) call adjust_corners(hp,crnmat,cvv(:,:,kd,zpxn,ypx,t+1),erridx)                                
                            if(erridx==1)then
                                sumerr = sumerr + erridx
                                !print*, 'warning -------------------------- 16'
                                cycle
                            endif
                            
                            do i = 1, 4
                                if(crnmat(i,3)/=-99._wp) cycle
                                
                                ! 4.14.2017 Actually the following condition statements all converge to the same calculation because of the introduction of the model's terminal condition at this point of time.
                                if(zpxn==1)then
                                    if(kd/=kdim-1)then
                                        kpxn = 0
                                        opxn = 0
                                        ypxn = 0
                                        
                                        !$omp atomic update !## 16
                                        cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) + &
                                            & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*survprob(t) ! no improvement as in the largest project. 
                                        sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*survprob(t)
                                    else ! kd==kdim-1
                                        kpxn = 0 ! die for sure
                                        opxn = 0
                                        ypxn = 0
                                        !$omp atomic update !## 17
                                        cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),kd,zpxn,ypx,kpxn,ypxn,opxn,t+1) + &
                                            & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*survprob(t)      
                                        sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*survprob(t)
                                    endif ! kd
                                else ! zpxn==0
                                    kpxn  = 0 ! die for sure
                                    opxn  = 0
                                    ypxn  = 0
                                    !$omp atomic update !## 18
                                    cef(crnmat(i,1),crnmat(i,2),0,zpxn,ypx,kpxn,ypxn,opxn,t+1) = cef(crnmat(i,1),crnmat(i,2),0,zpxn,ypx,kpxn,ypxn,opxn,t+1) + &
                                        & cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*survprob(t)   
                                    sum4 = sum4 + cef(ax,hx,kx,zx,yx,kpx,ypx,opx,t)*crnmat(i,4)*pz2(kd,zpxn+1)*survprob(t) ! 6:45 pm. 3.20.2017 kd is correct. Don't change it.  
                                endif ! zpxn==0
                            enddo ! i
                        enddo ! zpsn                         
                    endif ! kd
                    !$omp atomic update !## 19
                    sum3 = sum3 + sum4                    
                endif ! t      
            enddo ! q
            !$omp end do !## 20
            !$omp end parallel !## 21           
                        
            if(sumerr>0) exit_log1 = .true. ! 10132016. ! 3.27.2017 It is highly unlikely to have an error issued. But if it does, then it'd better to stop the program.
            
            deallocate( nvec )  
            if(printout8)then
                write(unit=128,fmt='(a,i3,(a,e25.15))'),   ' t+1 ', tdx+1, ' Expected arrival (sum2) ', sum2*survprob(tdx)
                write(unit=128,fmt='(a,i3,(a,e25.15))'),   ' t+1 ', tdx+1, ' Total mass transferred  ', sum3
                write(unit=128,fmt='(a,i3,(a,e25.15))'),   ' t+1 ', tdx+1, ' Result of transferred   ', sum(cef(:,:,:,:,:,:,:,:,tdx+1))
                write(unit=128,fmt='(a,i3,(a,e25.15),/)'), ' t+1 ', tdx+1, ' Mass target             ', popfrac(tdx+1)
            endif
        enddo ! tdx    
        
        deallocate( ivec )
        !call print_2period_cef() ! image of distribution
        
        if(inv_dist_counter>1)then
            error = maxval(dcef-cef)    
        endif ! inv_dist_counter
        
        !! 8-26-2017 to print out the mass sequence in the similar way as convert_2d_distribution_into_1d. 
        !! 8-27-2017 should start here. locmax should be multi-dimensiona becuase of dim of dcef and cef.
        !if(printout2)then
        !    locmax = maxloc(dcef-cef) 
        !    write(*,'(a,f20.13,a,9i10)') ' final mass distribution, error = ', error, ' maxloc ', locmax
        !    write(str1,fmt='(i3)') inv_dist_counter+300
        !    open(unit=300+inv_dist_counter,file="output_"//str1//"_cef.txt",action="write",status="replace")
        !    do i = 1, size(s3c(:,1))
        !        write(unit=300+inv_dist_counter,fmt='(i10,x,9(i2,x),2(e25.13,x))') &
        !            & s3c(i,10), s3c(i,1),  s3c(i,2), s3c(i,3), s3c(i,4), s3c(i,5), s3c(i,6), s3c(i,7), s3c(i,8), s3c(i,9), cef(s3c(i,1),s3c(i,2),s3c(i,3),s3c(i,4),s3c(i,5),s3c(i,6),s3c(i,7),s3c(i,8),s3c(i,9)), abs(dcef(s3c(i,1),s3c(i,2),s3c(i,3),s3c(i,4),s3c(i,5),s3c(i,6),s3c(i,7),s3c(i,8),s3c(i,9))-cef(s3c(i,1),s3c(i,2),s3c(i,3),s3c(i,4),s3c(i,5),s3c(i,6),s3c(i,7),s3c(i,8),s3c(i,9)))       
        !    enddo ! i
        !    close(300+inv_dist_counter)
        !endif ! printout2     

        ! 9-18-2017
        if(printout2)then
            !locmax = maxloc(dcef-cef) 
            !write(*,'(a,f20.13,a,9i10)') ' final mass distribution, error = ', error, ' maxloc ', locmax
            !write(str1,fmt='(i3)') inv_dist_counter+300
            write(str1, fmt='(i3)') trial_id+300
            write(str_idxr,fmt='(i3.3)') idxr
            write(str_idxg,fmt='(i3.3)') idxg
            write(str_idxt,fmt='(i3.3)') idxt
            open(unit=5000,file="output_"//str1//"_"//str_idxr//"-"//str_idxg//"-"//str_idxt//"_cef.txt",action="write",status="replace")
            do i = 1, size(s3c(:,1))
                write(unit=5000,fmt='(i10,x,9(i2,x),(e25.13,x))') &
                    & s3c(i,10), s3c(i,1),  s3c(i,2), s3c(i,3), s3c(i,4), s3c(i,5), s3c(i,6), s3c(i,7), s3c(i,8), s3c(i,9), cef(s3c(i,1),s3c(i,2),s3c(i,3),s3c(i,4),s3c(i,5),s3c(i,6),s3c(i,7),s3c(i,8),s3c(i,9))
            enddo ! i
            close(5000)
        endif ! printout2             
        
        dcef = cef ! 4.14.2017 for comparison in next round.
    end subroutine mass_transition  
    
    !subroutine intergenerational_transfer(error)
    subroutine intergenerational_transfer()
        implicit none
        integer :: sz, xv(2), q, idx, zxi, yxi, i, j, l, n, kpxi, opxi
        real(wp) :: wv(2), sum0, sum1, pyhprob, pyprob, sum2
        !real(wp), intent(out) :: error ! 4.14.2017 comment out
        integer :: tstart, tend, trate, tmax
        !real(wp) :: dsum1 ! debug
        
        !call system_clock(tstart,trate,tmax)
        !call system_clock(tstart)
        sz = size(s3c(:,1)) ! all the agents in the economy including the entrant. The length of the full list.
        ! sef1-sef3 are all of size szperiod1.
        sef1 = 0._wp ! mass recevied by valid 1st period combinations
        sef2 = sef(1:szperiod1) ! 4.1.2017 obtain in subroutine convert_2d_distribution_into_1d.
        sum0 = 0._wp ! mass that would have transferred to invalid 1st period combinations
        sef3 = 0._wp
        sum1 = 0._wp
        call linwgt(rav,transbeq,xv,wv)
        
        !!! 3.24.2017 use parallel computing the time is only half of the time used by sequential version (0.34-0.205 vs 0.669 seconds).
        !$omp parallel do default(shared) private(zxi,yxi,t,i,pyhprob,j,l,pyprob,n,kpxi,opxi,idx) ! #1#
        !dsum1 = 0._wp
        do q = 1, sz  
            ax  = s3c(q,1)
            hx  = s3c(q,2)
            kx  = s3c(q,3) ! The business project upon the parent's death.
            zxi = s3c(q,4) ! The business shock upon the parent's death.
            yxi = s3c(q,5) 
            kpx = s3c(q,6) ! useless for intergenerational transition
            ypx = s3c(q,7) ! useless for intergenerational transition
            opx = s3c(q,8) ! same as above
            t   = s3c(q,9)
            
            !sum1 = sum1 + cef(ax,hx,kx,zxi,yxi,kpx,ypx,opx,t)
            !dsum1 = dsum1 + cef(ax,hx,kx,zxi,yxi,kpx,ypx,opx,t) ! debug
            do i = 1, nmc ! child's labor efficiency at the beginning of period one.
                
                if(t>=10)then
                    pyhprob = syh(i)    
                else
                    pyhprob = pyh(yxi,i)
                endif                
                
                do j = 1, 2 ! linear interpolation for transbequets on financial asset dimensiona.
                    if(cvv(xv(j),1,kx,zxi,i,1)==-99)then ! the child inherits the parant's career attributes (kx and zxi) upon death.
                        do l = 1, nmc ! child's labor efficiency at the end of peirod one.
                            pyprob = py(i,l) ! the transition probability of child's labor efficiency. <-----
                            do n = 1, 2 ! child's business idea innovation (go advance or not) <-----
                                
                                if(zxi==0)then ! bad shock
                                    !if(kx==0)then ! 4.14.2017 comment out BUG!
                                    kpxi = merge(0,1,n==1)
                                    opxi = merge(0,2,n==1)                                        
                                    !else ! kx/=0 ! 4.14.2017 comment out
                                    !    kpxi = 0
                                    !    opxi = 0
                                    !endif ! kx
                                else ! zxi /= 0 good shock
                                    if(kx/=kdim-1)then
                                        kpxi = kx+(n-1) ! merge(kx+(n-1),kx+(n-1),n==1)
                                        opxi = n
                                    else ! kx==kdim-1
                                        if(n==2)then
                                            cycle ! 4.12.2017 correct.
                                        else ! n == 1
                                            kpxi = 3
                                            opxi = 1
                                        endif ! n
                                    endif
                                endif ! zxi
                                idx = c3s(xv(j),1,kx,zxi,i,kpxi,l,opxi,1) ! the addresss at the end of period one for storeing the transfered mass 
                                !$omp atomic update ! #2#
                                sef1(idx) = sef1(idx) + wv(j)*cef(ax,hx,kx,zxi,yxi,kpx,ypx,opx,t)*(1._wp-survprob(t))*pyhprob*pyprob*pka(kx,n)    
                                
                            enddo ! n
                        enddo ! l
                    else ! cvv "bad combination"
                        !$omp atomic update !#3#
                        sum0 = sum0 + wv(j)*cef(ax,hx,kx,zxi,yxi,kpx,ypx,opx,t)*(1._wp-survprob(t))*pyhprob ! Relative to the above eqa. for sef1, two items are absent: pyprob and pka.
                        
                        !$omp critical ! #4$
                        !write(unit=128,fmt='(a,i4,a,8i4)') ' round: ', inv_dist_counter, ' bad combination: ', ax,hx,kx,zxi,yxi,kpx,ypx,opx
                        !$omp end critical ! #5#
                    endif ! cvv(.)
                enddo ! j
            enddo ! i
        enddo ! q
        !$omp end parallel do ! #6#
        
        !call system_clock(tend)
        !write(*,fmt='(a,f12.4,a)') ' time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'         
        !print*, ' valid transfer ', sum(sef1), ' total transfer ', sum(sef1)+sum0
        
        if(sum0>0._wp)then ! 4.1.2017 adjusted for bad combination.
            !print*, ' I am here -------------------------- ', sum0 ! `sum0` is supposed to be small. ~1.e-6_wp
            sef3 = sef1 ! let sef3 be the backup of the newly obtained period one's distribution.
            sum2 = sum(sef3(1:szperiod1))
            do q = 1, szperiod1
                !if(sef1(q)/=0._wp)then
                if(sef3(q)>0._wp)then
                    sef3(q) = sef3(q)/sum2
                    sef1(q) = sef1(q) + sum0*sef3(q)
                endif
            enddo ! q
        endif ! sum0
        
        !print*, ' ------------- dsum1: ', dsum1 ! In a nomral case, `dsum1` is supposed to be one.
        ! note: return "sef1"... 4.1.2017 
        
        !print*, ' valid transfer plus equal addition ', sum(sef1(1:szperiod1)), ' sum of equal addition ', sum0, ' sum allocation ratio (should be one) ', sum(sef3)
        
        !error = maxval(abs(sef1-sef2)) ! sef1: new end of 1st period mass; sef2: old one. ! 4.14.2017 comment out
    end subroutine intergenerational_transfer
    
    ! 4.17.2017 accidental bequests.
    subroutine lump_sum_transfer()
    ! Yang 2009, accidental bequests are distributed evenly to new agents at age 25, which endogenously determines bt. However, in the model formulation, the bequests are evenly distributed to all the surviving agents using a lump-sum transfer.
    ! Nakajima 2010, The government imposes a 100% estate tax rate on `accidental bequests` and distributed all the proceeds equally to all the surviving agents using a lump-sump transfer.
    ! Conesa et al. 2009 There are no annuity markets and therefore a fraction of households leaves unintended bequests, denoted by Tr_t, that are redistributed in a lump-sum manner across individuals currently alive.
    ! 3.24.2017 (Fri) checked.
        implicit none
        integer :: sz, q, yxi, zxi
        real(wp) :: sum2
        integer :: tstart, tend, trate, tmax
        sz = size(s3c(:,1))
        !if(inv_dist_counter==1)then
        !    !allocate(sef1(szperiod1),sef2(szperiod1))    !<------------------
        !endif
        
        call system_clock(tstart,trate,tmax)
        call system_clock(tstart)        
        sum2 = 0._wp
        ! 3.24.2017 Parallel computing is good. Don't comment out it. 0.016 vs 0.038
        !$omp parallel do default(shared) private(zxi,yxi,t) reduction(+:sum2)
        do q = 1, sz
            ax  = s3c(q,1)
            hx  = s3c(q,2)
            kx  = s3c(q,3)
            zxi = s3c(q,4)
            yxi = s3c(q,5)
            kpx = s3c(q,6)
            ypx = s3c(q,7)
            opx = s3c(q,8)
            t   = s3c(q,9)
            !if(cww(ax,hx,kx,zxi,yxi,kpx,ypx,opx,t)/=-99) cycle
            !if(cef(ax,hx,kx,zxi,yxi,kpx,ypx,opx,t)==0._wp) cycle 
            if(abs(cef(ax,hx,kx,zxi,yxi,kpx,ypx,opx,t))<1.e-10_wp) cycle 
            sum2 = sum2 + cef(ax,hx,kx,zxi,yxi,kpx,ypx,opx,t)*(1._wp-survprob(t))*((1._wp+rd)*cwa(ax,hx,kx,zxi,yxi,kpx,ypx,opx,t)+(1._wp-deltah)*cwh(ax,hx,kx,zxi,yxi,kpx,ypx,opx,t))
        enddo ! q
        !$omp end parallel do
        call system_clock(tend)
        
        !transbeqimplied = sum2 ! five year per capita. 3.24.2017 4.17.2017 This quantity is on five-year baiss. Correct!
        transbeqimplied = merge(sum2, sum2/popfrac(1) , printout17) ! 9-30-2017
        !write(*,fmt='(a,f12.4,a,f12.4)') ' lump sum transfer time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds. transbeqimplied: ', transbeqimplied          
    end subroutine lump_sum_transfer
    
    subroutine macro_statistics(momvec,idxr,idxg,idxt,exit_log1,msg,trial_id)
        implicit none
        integer :: sz, q, idx, lsz, i, j
        real(wp) :: capgain, intfund, bgp, massvec(3) ! 9-12-2017 capgain is useless. bgp refers to business gross profit. ! 9-17-2017 inc should be removed.
        logical, dimension(:), allocatable :: lvece, lvecw, tvec
        real(wp), dimension(:), allocatable :: svec, rvec
        integer, intent(in) :: idxr, idxg, idxt,trial_id ! (iterar, iteragov, iteratot)
        !real(wp) :: totast1, wokfin1, entfin1, wokhom1, enthom1, entcap1, crpcap1, totefl1, entlab1, crplab1, tottax1, entprd1, crpprd1, gdp1, totsvt1, govbal1, govbal2gdp1, hug_inv_proj1, med_inv_proj1, sml_inv_proj1, all_inv_proj1, hug_inv_per1, med_inv_per1, sml_inv_per1, ent_wel_per1, woksize1, chgw2e1, w2erat1, &
        !    entsize1, chge2w1, e2wrat1, medwokinc1, lowest_quintile_wokinc1, medwokwel1, medentwel1, entcsp1, wokcsp1, entinc1, wokinc1, all_income1, ent_inc_per1
        integer :: szwok, szent
        integer :: tstart, tend, trate, tmax
        character(len=4) :: str1
        character(len=3) :: str2
        ! real(wp) :: r_govbal
        ! output
        real(wp), dimension(:), intent(out) :: momvec
        logical, intent(inout) :: exit_log1
        character(len=*), intent(out) :: msg
        
        call system_clock(tstart,trate,tmax)
        call system_clock(tstart)        
        
        ! 9-18-2017
        sef = int(sef*accumass)/accumass
        
        sz = size(s3c(:,1))
        allocate( svec(sz), tvec(sz) ) ! svec: real, tvec: logical.
        
        tvec = .false. ! 9-12-2017
        svec = 0._wp   ! 9-12-2017
        term_2 = 0._wp ! 9-17-2017
        term_3 = 0._wp
        term_4 = 0._wp
        term_5 = 0._wp
        term_6 = 0._wp
        term_7 = 0._wp
        term_8 = 0._wp
        term_9 = 0._wp
        
        ! macro statistics ! 4.16.2017 moved from equilibrium.f90 This movement solved a NaN-in-tottax problem.
        sw_laborsupply   = 0._wp
        sw_labordemand   = 0._wp
        sw_production    = 0._wp
        sw_bizinvestment = 0._wp
        sw_bizloan       = 0._wp
        sw_ini_asset     = 0._wp
        sw_ini_house     = 0._wp
        sw_nonlineartax  = 0._wp
        sw_socialsecurity= 0._wp
        sw_worker_turned = 0._wp
        sw_boss_turned   = 0._wp
        sw_buzcap_notuse = 0._wp
        sw_aftertaxwealth= 0._wp
        sw_consumption   = 0._wp
        sw_worker_savtax = 0._wp
        sw_entpre_savtax = 0._wp
        sw_entpre_biztax = 0._wp  
        sw_taxableincome = 0._wp
        sw_totinc_bx     = 0._wp
        sw_wealth_tax    = 0._wp
        sw_net_worth     = 0._wp
        
        ! 3.27.2017 0.31 vs 0.7 parallel still outperform sequential.
        !$omp parallel do default(shared) private(t,idx,capgain,intfund,bgp) 
        do q = 1, sz
            ax  = s3c(q,1)
            hx  = s3c(q,2)
            kx  = s3c(q,3)
            zx  = s3c(q,4)
            yx  = s3c(q,5)
            kpx = s3c(q,6)
            ypx = s3c(q,7)
            opx = s3c(q,8)
            t   = s3c(q,9)  ! 9-12-2017 not in the list of threadprivate.
            idx = s3c(q,10) ! 9-12-2017 not in the list of threadprivate.
            
            if(sef(idx)<=0._wp) cycle ! "sef" (and "cef") is obtained from subroutine convert_2d_distribution_into_1d. 
            !atw = 0._wp
            sw_ini_asset(idx) = rav(ax) ! 4.15.2015 Sum over the beginning-of-period financial asset holdings. Correct. 
            sw_ini_house(idx) = rhv(hx) ! 4.15.2017 Sum over the beginning-of-period housing asset holdings. Correct. 
            
            a = rav(ax)
            h = rhv(hx)
            
            !! =======================  Business behviors. Supply and demand of labor. =============== ! 4.15.2017 Stop here. 4:00 pm.
            if(kx==0)then ! salaried worker (or retiree)
                sw_production(idx) = 0._wp ! production in coarse_SBE_profit_matrices
                sw_bizinvestment(idx) = 0._wp ! k
                sw_bizloan(idx) = 0._wp ! (k-a)
                sw_buzcap_notuse(idx) = 0._wp ! 10.13.2017
                
                if(t>=10)then ! retired
                    sw_laborsupply(idx) = 0._wp ! labor supply to the whole economy
                    sw_labordemand(idx) = 0._wp ! SME sector's labor demand
                else ! working
                    sw_laborsupply(idx) = efflab(t)*yv(yx) ! labor supply to the whole economy
                    sw_labordemand(idx) = 0._wp ! SME sector's labor demand
                endif ! t
            else ! entrepreneurs
                if(t>=10)then ! old entrepreneurs
                    sw_laborsupply(idx) = 0._wp 
                    if(zx==0)then ! choose to be a worker.
                        sw_labordemand(idx) = 0._wp 
                        sw_production(idx) = 0._wp
                        sw_bizinvestment(idx) = 0._wp ! 9-13-2017 kv(kx) ! IDLE CAPITAL <------ Bug 10122016 kept in the entrepreneurial sector, but not put into production
                        sw_buzcap_notuse(idx) = kv(kx) ! 3.25.2017 ! 4.15.2017 consistent with the formulation of disposable income. 
                        !sw_bizloan(idx) = merge( kv(kx)-sw_ini_asset(idx), 0._wp, kv(kx)>sw_ini_asset(idx))
                        sw_bizloan(idx) = merge( merge(kv(kx)-a, kv(kx), a>0._wp), 0._wp, kv(kx)>a) !10.13.2017
                    else ! good shocks
                        sw_labordemand(idx) = c_lab_vec(kx)
                        sw_production(idx) = c_opt_vec(kx)
                        sw_bizinvestment(idx) = kv(kx)
                        sw_buzcap_notuse(idx) = 0._wp ! 3.25.2017
                        !sw_bizloan(idx) = merge( kv(kx)-sw_ini_asset(idx), 0._wp, kv(kx)>sw_ini_asset(idx))
                        sw_bizloan(idx) = merge( merge(kv(kx)-a, kv(kx), a>0._wp), 0._wp, kv(kx)>a) !10.13.2017
                    endif ! zx 
                else ! young entrepreneurs (t<=9)
                    if(zx==0)then
                        sw_laborsupply(idx) = efflab(t)*yv(yx) ! 0._wp ! 10122016 Supplied to the other enterprise or corporate sector.  
                        sw_labordemand(idx) = 0._wp
                        sw_production(idx) = 0._wp
                        sw_bizinvestment(idx) = 0._wp ! 9-13-2017 removed kv(kx) ! idle capital being locked in the entrepreneurial sector.
                        sw_buzcap_notuse(idx) = kv(kx) ! 3.25.2017
                        !sw_bizloan(idx) = merge( kv(kx)-sw_ini_asset(idx), 0._wp, kv(kx)>sw_ini_asset(idx)) ! 3.27.2017 sw_ini_asset(idx) rather than sw_ini_asset(ax)     
                        sw_bizloan(idx) = merge( merge(kv(kx)-a, kv(kx), a>0._wp), 0._wp, kv(kx)>a) !10.13.2017
                    else ! good shocks zx>0._wp
                        !! 3.25.2017 comment out
                        !if(c_lab_vec(kx)>efflab(t)*yv(yx))then
                        !    sw_laborsupply(idx) = efflab(t)*yv(yx) ! own labor efficiency 
                        !    sw_labordemand(idx) = c_lab_vec(kx) ! labor demand (including own labor supply) 20132016                 
                        !else ! weak labor demand
                        !    sw_laborsupply(idx) = c_lab_vec(kx) ! own labor efficiency is not used up. 10122016. <--------------
                        !    sw_labordemand(idx) = c_lab_vec(kx) ! labor demand 
                        !endif ! c_lab_vec>efflab(t)*yv(yx)
                        
                        ! 3.25.2017 revision made.
                        sw_laborsupply(idx) = efflab(t)*yv(yx) ! 0._wp ! 10122016 Supplied to the other enterprise or corporate sector.  
                        sw_labordemand(idx) = c_lab_vec(kx)
                        sw_production(idx)  = c_opt_vec(kx)                 
                        sw_bizinvestment(idx) = kv(kx) ! idle capital being locked in the entrepreneurial sector.
                        sw_buzcap_notuse(idx) = 0._wp ! 3.25.2017
                        !sw_bizloan(idx) = merge( kv(kx)-sw_ini_asset(idx), 0._wp, kv(kx)>sw_ini_asset(idx))    
                        sw_bizloan(idx) = merge( merge(kv(kx)-a, kv(kx), a>0._wp), 0._wp, kv(kx)>a) !10.13.2017
                    endif ! zx
                endif ! t
            endif ! kx == 0 or not (worker/retiree or entrepreneur)
            flon(idx) = sw_bizloan(idx)
            
            ! ===================================Income section=============================================== 3.26.2017 4:51 pm stop here. 4.15.2017 stop here. 9:20 pm.
            if(kx==0)then ! worker or retirees.
                !capgain = merge(rd*a, 0._wp, rd*a>0._wp)   
                !inc = merge(benefit,0._wp,t>=10) + merge(0._wp,capgain,tausvflag) + wage*sw_laborsupply(idx)
                
                if(printout25)then

                    if(mode6taskid==0)then
                        sw_socialsecurity(idx) = tauss/2._wp*wage*sw_laborsupply(idx)
                        sw_taxableincome(idx)  = merge(benefit, 0._wp, t>=10) + merge(0._wp, transbeq, printout17==.False. .and. t/=1) &
                                                 + merge(0._wp, merge(rd*a,0._wp,a>0._wp), tausv>0._wp) &
                                                 + wage*sw_laborsupply(idx) - sw_socialsecurity(idx)
                        sw_nonlineartax(idx)   = ttaxwok(sw_taxableincome(idx)) ! 9-17-2017
                        sw_worker_savtax(idx)  = merge(tausv*rd*a, 0._wp, a>0._wp) ! 10.27.2017 negative interest rate is acceptible.
                        sw_totinc_bx(idx)      = sw_taxableincome(idx) + merge(merge(rd*a,0._wp,a>0._wp), 0._wp, tausv>0._wp) !10.13.2017
                        sw_net_worth(idx)      = sw_totinc_bx(idx) + (1._wp-deltah)*h + merge(a,(1._wp+rd)*a,a>0._wp)
                        fnon(idx) = sw_nonlineartax(idx)
                        fsax(idx) = sw_worker_savtax(idx)                         
                        sw_wealth_tax(idx)     = 0._wp !(1._wp+rd)*a + (1._wp-deltah)*h ! useless in mode6taskid==0 case.
                        sw_aftertaxwealth(idx) = sw_net_worth(idx) - sw_nonlineartax(idx) - sw_worker_savtax(idx)
                        !sw_aftertaxwealth(idx) = sw_taxableincome(idx) + merge((1._wp+rd)*a, (1._wp+(1._wp-tausv)*rd)*a, a<0._wp) &
                        !                       + (1._wp-deltah)*h - ttaxwok(sw_taxableincome(idx))
                        fwtx(idx) = sw_wealth_tax(idx)
                    elseif(mode6taskid==1)then ! wealth tax (not applied on period income but beginning capital, regardless of asset types).
                        sw_socialsecurity(idx) = tauss/2._wp*wage*sw_laborsupply(idx)
                        sw_taxableincome(idx)  = merge(benefit, 0._wp, t>=10) + merge(0._wp, transbeq, printout17==.False. .and. t/=1) &
                                                 + merge(rd*a, 0._wp, a>0._wp) & 
                                                 + wage*sw_laborsupply(idx) - sw_socialsecurity(idx)
                        sw_nonlineartax(idx)   = 0._wp ! 3-1-2018 ttaxwok(sw_taxableincome(idx)) ! 9-17-2017
                        sw_worker_savtax(idx)  = 0._wp ! merge(tausv*rd*a, 0._wp, a>0._wp) ! 10.27.2017 negative interest rate is acceptible.
                        sw_totinc_bx(idx)      = sw_taxableincome(idx) ! + merge(rd*a, 0._wp, a>0._wp) !10.13.2017
                        sw_net_worth(idx)      = sw_totinc_bx(idx) + (1._wp-deltah)*h + merge(a, (1._wp+rd)*a, a>0._wp) ! 3-1-2018 a comprehensive amount of income except for the SS contribution.
                        fnon(idx) = sw_nonlineartax(idx) ! 0._wp actually. 3-1-2018
                        fsax(idx) = sw_worker_savtax(idx)                          
                        sw_wealth_tax(idx)     = taubal * sw_net_worth(idx) ! 3-3-2018
                        sw_aftertaxwealth(idx) = merge((1._wp - taubal) * sw_net_worth(idx), sw_net_worth(idx), sw_wealth_tax(idx) > 0._wp) ! 3-1-2018
                        sw_wealth_tax(idx)     = merge(sw_wealth_tax(idx), 0._wp, sw_wealth_tax(idx) > 0._wp)
                        fwtx(idx) = sw_wealth_tax(idx)
                    elseif(mode6taskid==2)then
                        sw_socialsecurity(idx) = tauss/2._wp*wage*sw_laborsupply(idx)
                        sw_taxableincome(idx)  = merge(benefit, 0._wp, t>=10) + merge(0._wp, transbeq, printout17==.False. .and. t/=1) &
                                                 + merge(rd*a,0._wp,a>0._wp) &
                                                 + wage*sw_laborsupply(idx) - sw_socialsecurity(idx)
                        sw_nonlineartax(idx)   = ttaxwok(sw_taxableincome(idx)) ! 9-17-2017
                        sw_worker_savtax(idx)  = 0._wp ! 10.27.2017 negative interest rate is acceptible.
                        sw_totinc_bx(idx)      = sw_taxableincome(idx) !10.13.2017
                        sw_net_worth(idx)      = sw_totinc_bx(idx) + (1._wp-deltah)*h + merge(a, (1._wp+rd)*a, a>0._wp) ! 3-1-2018 a comprehensive amount of income except for the SS contribution.
                        fnon(idx) = sw_nonlineartax(idx)
                        fsax(idx) = sw_worker_savtax(idx)                         
                        net_worth = sw_net_worth(idx)
                        if(net_worth > 5._wp*exempbar)then
                            sw_wealth_tax(idx) = 0.02_wp*(net_worth-5._wp*exempbar) + 0.01_wp*(4._wp*exempbar)
                            sw_aftertaxwealth(idx) = sw_net_worth(idx) - sw_nonlineartax(idx) - sw_wealth_tax(idx)        
                        elseif(net_worth>exempbar .and. net_worth<=5._wp*exempbar)then
                            sw_wealth_tax(idx) = 0.01_wp*(net_worth - exempbar)
                            sw_aftertaxwealth(idx) = sw_net_worth(idx) - sw_nonlineartax(idx) - sw_wealth_tax(idx)
                        else    
                            sw_wealth_tax(idx) = 0._wp
                            sw_aftertaxwealth(idx) = sw_net_worth(idx) - sw_nonlineartax(idx) - sw_wealth_tax(idx)
                        endif                       
                        fwtx(idx) = sw_wealth_tax(idx)
                    else
                        print*, "error in macro statistics-1"    
                    endif !mode6taskid
                    
                else
                    
                    sw_taxableincome(idx) = merge(benefit, 0._wp, t>=10) + wage*sw_laborsupply(idx) ! 3.27.2017 ! 9-12-2017 Interest income will be taxed below.
                    sw_nonlineartax(idx)   = ttaxwok(sw_taxableincome(idx)) ! 9-17-2017
                    sw_socialsecurity(idx) = tauss/2._wp*wage*sw_laborsupply(idx) ! 3.25.2017        
                    sw_worker_savtax(idx)  = merge(tausv*rd*a, 0._wp, a>0._wp) ! 10.27.2017 negative interest rate is acceptible.
                    sw_totinc_bx(idx)      = sw_taxableincome(idx) + merge(rd*a, 0._wp, a>0._wp) !10.13.2017
                    sw_aftertaxwealth(idx) = sw_taxableincome(idx) + merge( (1._wp+rd)*a, (1._wp+(1._wp-tausv)*rd)*a, a<0._wp) &
                                           & + (1._wp-deltah)*h - ttaxwok(sw_taxableincome(idx)) + merge(0._wp, transbeq, printout17==.False. .and. t/=1) - sw_socialsecurity(idx) ! 3.27.2017 ! 9-12-2017 at this line, inc ~= wealth before tax.
                endif ! printout25
                
                sw_worker_turned(idx) = merge( 1._wp, 0._wp, s3c(idx,3)==0.and.swk(idx)==1) ! 10.13.2017 useless !4.16.2017, 3.25.2017 swk=cwk. If cwk==1, the career "NEXT" peirod is to work for a wage.
            else ! kx>0 Boss
                intfund = merge(a-kv(kx), 0._wp, a-kv(kx)>0._wp) ! 3.27.2017 sw_ini_asset(idx) rather than sw_ini_asset(ax). ! sw_bizloan(idx) ! Bug intfund wrong!! 10122016 
                bgp = c_grs_mat(ax, kx, zx, yx, t) 
                if(printout25)then                  
                    if(mode6taskid==0)then
                        sw_socialsecurity(idx) = tauss/2._wp*wage*sw_laborsupply(idx)
                        sw_taxableincome(idx)  = bgp + merge(benefit, 0._wp, t>=10) + merge(0._wp, transbeq, printout17==.False. .and. t/=1) &
                                                 + merge(0._wp, merge(rd*intfund,0._wp,intfund>0._wp), tausv>0._wp) &
                                                 + wage * sw_laborsupply(idx) - sw_socialsecurity(idx) ! 10.26.2017 revision                    
                        sw_nonlineartax(idx)   = ttaxent(sw_taxableincome(idx)) ! 9-17-2017
                        sw_entpre_savtax(idx)  = merge( tausv*rd*intfund, 0._wp, intfund>0._wp) !10.13.2017
                        sw_totinc_bx(idx)      = sw_taxableincome(idx) + merge(merge(rd*intfund,0._wp,intfund>0._wp), 0._wp, tausv>0._wp)
                        sw_net_worth(idx)      = sw_totinc_bx(idx) + (1._wp-deltah)*h + merge(intfund, 0._wp, intfund>0._wp) ! 3-1-2017 add bgp. Useless.                
                        fnon(idx) = sw_nonlineartax(idx)
                        fsax(idx) = sw_entpre_savtax(idx)                          
                        sw_wealth_tax(idx)     = 0._wp ! 3-1-2018 Useless here. used only for mode6taskid=1. Here, this variable is for saving the information of household net worth.
                        sw_aftertaxwealth(idx) = sw_net_worth(idx) - sw_nonlineartax(idx) - sw_entpre_savtax(idx) ! 3-1-2018
                        fwtx(idx) = sw_wealth_tax(idx)                        
                    elseif(mode6taskid==1)then ! wealth tax case. tausv should be set to zero, so we can ignore sw_XXX_savtax.
                        sw_socialsecurity(idx) = tauss/2._wp*wage*sw_laborsupply(idx)
                        sw_taxableincome(idx)  = bgp + merge(benefit, 0._wp, t>=10) + merge(0._wp, transbeq, printout17==.False. .and. t/=1) &
                                                 + merge(rd * intfund, 0._wp, intfund > 0._wp) & ! 10.26.2017 revision                    
                                                 + wage * sw_laborsupply(idx) - sw_socialsecurity(idx)
                        sw_nonlineartax(idx)   = 0._wp ! ttaxent(sw_taxableincome(idx)) ! 3-1-2018 Useless in mode6taskid==1. 9-17-2017
                        sw_entpre_savtax(idx)  = 0._wp ! 3-1-2018 merge( tausv*rd*intfund, 0._wp, intfund>0._wp) !10.13.2017
                        sw_totinc_bx(idx)      = sw_taxableincome(idx) !10.13.2017, 3-1-2018
                        sw_net_worth(idx)      = sw_totinc_bx(idx) + merge(intfund, 0._wp, intfund > 0._wp) + (1._wp-deltah)*h  ! 3-1-2018 needs to be revised for the tax-on-toprich case in the future.
                        fnon(idx) = sw_nonlineartax(idx)
                        fsax(idx) = sw_entpre_savtax(idx)                          
                        sw_wealth_tax(idx)     = taubal * sw_net_worth(idx) ! variable used in the conditional statement below 3-3-2018
                        sw_aftertaxwealth(idx) = merge((1._wp - taubal) * sw_net_worth(idx), sw_net_worth(idx), sw_wealth_tax(idx) > 0._wp) ! 3-3-2018      
                        sw_wealth_tax(idx)     = merge(sw_wealth_tax(idx), 0._wp, sw_wealth_tax(idx) > 0._wp)    
                        fwtx(idx) = sw_wealth_tax(idx)
                    elseif(mode6taskid==2)then   
                        sw_socialsecurity(idx) = tauss/2._wp*wage*sw_laborsupply(idx)
                        sw_taxableincome(idx)  = bgp + merge(benefit, 0._wp, t>=10) + merge(0._wp, transbeq, printout17==.False. .and. t/=1) &
                                                 + merge(rd * intfund, 0._wp, intfund > 0._wp) & ! 10.26.2017 revision                    
                                                 + wage * sw_laborsupply(idx) - sw_socialsecurity(idx)
                        sw_nonlineartax(idx)   = ttaxent(sw_taxableincome(idx)) ! 9-17-2017
                        sw_entpre_savtax(idx)  = 0._wp
                        sw_totinc_bx(idx)      = sw_taxableincome(idx) !10.13.2017
                        sw_net_worth(idx)      = sw_totinc_bx(idx) + merge(intfund, 0._wp, intfund > 0._wp) + (1._wp-deltah)*h 
                        fnon(idx) = sw_nonlineartax(idx)
                        fsax(idx) = sw_entpre_savtax(idx)  
                        net_worth = sw_net_worth(idx)
                        if(net_worth > 5._wp*exempbar)then
                            sw_wealth_tax(idx) = 0.02_wp*(net_worth-5._wp*exempbar) + 0.01_wp*(4._wp*exempbar)
                            sw_aftertaxwealth(idx) = sw_net_worth(idx) - sw_nonlineartax(idx) - sw_wealth_tax(idx)        
                        elseif(net_worth>exempbar .and. net_worth<=5._wp*exempbar)then
                            sw_wealth_tax(idx) = 0.01_wp*(net_worth - exempbar)
                            sw_aftertaxwealth(idx) = sw_net_worth(idx) - sw_nonlineartax(idx) - sw_wealth_tax(idx)
                        else    
                            sw_wealth_tax(idx) = 0._wp
                            sw_aftertaxwealth(idx) = sw_net_worth(idx) - sw_nonlineartax(idx) - sw_wealth_tax(idx)
                        endif  
                        fwtx(idx) = sw_wealth_tax(idx)
                    else
                        print*, "error in macro statistics-2"
                    endif
                else
                    sw_taxableincome(idx)  = merge(benefit, 0._wp, t>=10) + wage*sw_laborsupply(idx) + bgp ! 4.16.2017 add bgp. 
                    sw_nonlineartax(idx)   = ttaxent(sw_taxableincome(idx)) ! 9-17-2017
                    sw_socialsecurity(idx) = tauss/2._wp*wage*sw_laborsupply(idx)            
                    sw_entpre_savtax(idx)  = merge( tausv*rd*intfund, 0._wp, intfund>0._wp) !10.13.2017
                    sw_totinc_bx(idx)      = sw_taxableincome(idx) + merge( rd*intfund, 0._wp, intfund>0._wp) !10.13.2017
                    sw_aftertaxwealth(idx)  = sw_taxableincome(idx) +  merge( (1._wp+(1._wp-tausv)*rd)*intfund, merge( 0._wp, (1._wp+rd)*a, a>0._wp), intfund>0._wp) & !10.13.2017
                                            & + (1._wp-deltah)*h - ttaxent(sw_taxableincome(idx)) + merge(0._wp, transbeq, printout17==.False. .and. t/=1) - sw_socialsecurity(idx) ! 4.16.2017 revision
                endif
                sw_boss_turned(idx) = merge( 1._wp, 0._wp, (s3c(idx,3)/=0.and.swk(idx)==0)) ! 3.27.2017 it is only a portion of the whole population exiting entrepreneurship.
                
                !! 3.25.2017 comment out.
                !if(zx==0)then
                !    sw_socialsecurity(idx) = tauss*wage*sw_laborsupply(idx) ! applicable to old entrepreneurs.
                !else
                !    sw_socialsecurity(idx) = merge( tauss*wage*(c_lab_vec(kx)-sw_laborsupply(idx))+2*tauss*wage*sw_laborsupply(idx), 2*tauss*wage*sw_laborsupply(idx), c_lab_vec(kx)>sw_laborsupply(idx)) ! revised in the future
                !endif
                !inc = inc + a + (1._wp-tausv)*merge(capgain,0._wp,tausvflag) + (1._wp-deltah)*h + (1._wp-taubp)*c_grs_mat(ax,kx,zx,yx,t) - sw_nonlineartax(idx) + transbeq - sw_socialsecurity(idx) ! 3.25.2017 replace tauk with tausv
                       
            endif ! kx
            sw_consumption(idx) = cwc(ax,hx,kx,zx,yx,kpx,ypx,opx,t) ! 3.25.2017 no interpolation required because I don't use refine grid any more.
            
            ! if(sw_nonlineartax(idx)/=sw_nonlineartax(idx)) write(unit=128,fmt='(9(a,x,f8.4),(a,x,f8.4))') 'ax',ax,'hx',hx,'kx',kx,'zx',zx,'yx',yx,'kpx',kpx,'ypx',ypx,'opx',opx,'t',t,' consumpt ', sw_nonlineartax(idx)
            ! write(unit=128,fmt='(i8,9(a,x,i4,", "),2(a,x,f8.4))') idx, 'ax',ax,'hx',hx,'kx',kx,'zx',zx,'yx',yx,'kpx',kpx,'ypx',ypx,'opx',opx,'t',t,' consumpt ', sw_nonlineartax(idx), ' mass ', sef(idx)
        enddo ! q
        !$omp end parallel do
        
        !call system_clock(tend)
        !!if(printout6) write(unit=120,fmt='(a,i3,a,f12.4,a)') 'solve bellman period ',t,', time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'         
        !if(printout6) write(*,fmt='(a,f12.4,a)') 'solve macro statistic time: ', real(tend-tstart,wp)/real(trate,wp), ' seconds'    
        
        !tvec = (abs(sw_buzcap_notuse-0._wp)<1.e-10_wp).and.s3c(:,3)>0.and.s3c(:,4)==1 ! 5.9.2017 Note: s3c(:,3): k; s3c(:,4): z. ! 9-13-2017 this is entrepreneur's version uses sw_buzcap_notuse<0
        !tvec = s3c(:,3)>0.and.s3c(:,4)==1 ! 9-13-2017 Entrepreneurs running a business. !10.13.2017 comment out
        tvec = s3c(:,3)>0
        
        ! 4.16.2017 12:05pm stop here.
        !totast = dot_product(sef,sw_ini_asset) ! 10.18.2017 It's net financial asset holdings (assets net of liabilities).
        totast = sum(sef*sw_ini_asset, sef>=0._wp) ! 10.26.2017
        svec   = sef*sw_ini_asset
        !wokfin = sum(svec,s3c(:,3)==0.and.s3c(:,9)<=9) ! retirees are included into the working class household hereafter.
        if(printout22)then ! 10.18.2017
            wokfin = sum(svec, tvec==.false. .and. sw_ini_asset>0._wp .and. sef>=0._wp)  
            entfin = sum(svec, tvec==.true.  .and. sw_ini_asset>0._wp .and. sef>=0._wp)   
        else
            wokfin = sum(svec, tvec==.false. .and. sef>=0._wp)
            entfin = sum(svec, tvec==.true.  .and. sef>=0._wp)             
        endif 
        svec   = sef*sw_ini_house
        !wokhom = sum(svec,s3c(:,3)==0.and.s3c(:,9)<=9)
        wokhom = sum(svec, tvec==.false. .and. sef>0._wp) ! 3.27.2017 added s3c condition. 7-2-2017 Including retired people.
        enthom = sum(svec, tvec==.true.  .and. sef>0._wp) ! 3.27.2017 added s3c condition.
        
        svec   = sef*(sw_bizinvestment + sw_buzcap_notuse) ! 9-13-2017 definition of bizinvestment defined above is updated. only z>0 people have positive bizinvestment.
        entcap = sum(svec, tvec==.true. .and. sef>0._wp ) ! 3.27.2017 Only sum up those figures not coming from idle capital. ! 8-14-2017 ok.
        if(printout22)then
            crpcap = (wokfin+entfin)/(1._wp+dfrac)-entcap
        else
            crpcap = totast/(1._wp+dfrac)-entcap ! Cagetti and De Nardi, TAUCfinalSS.f90, line 1273. a constant fraction of total capital. ! <========================<<
        endif
            
        ! 3.27.2017 following cagetti and de nardi's 2009 AER paper. See their InitialSS.f90 in AER_2009_codes.zip (unzipped as MS200405121codes). Their model
        ! also consider government debt accoutns for a fraction of total capital.
        ! totast = entcap + crpcap + dfrac*(entcap + crpcap) 
        ! 8-18-2017 So we know the functional relationship: crpcap = f(totast, entcap, dfrac).
        
        svec   = sef*sw_laborsupply
        ! totefl = sum(svec,s3c(:,3)==0) ! Total efficiency unit of labor is provided by workers. Note that entrepreneurs' efficiency unit stay only in the entrepreneurial sector.
        totefl = sum(svec, sef>0._wp) ! the total efficiency units provided in the economy. 10122016.
        
        svec   = sef*sw_labordemand ! So it is easy to obtain the efficiency unit of labor put into the non-entrepreneurial sector by subtracting external labor demand of entreprises from the total efficiency units provided by workers.
        entlab = sum(svec,tvec==.true. .and. sef>0._wp) ! 3.27.2017 screen out the contribution that comes from entrepreneurs whose business is not in operation ! total external efficiency units of labor demanded by entreprises.
        crplab = totefl - entlab 
        
        tottax = sum(sef*sw_nonlineartax, sef>0._wp) ! 9-13-2017 it needs revision, because it doesn't contain social security taxes and capital gain tax from savings interest.
        
        entprd = sum(sef*sw_production, sef>0._wp) ! 9-12-2017 Actually it is the production of entrepreneurial sector.
        crpprd = wage*crplab + (rd+deltak)*crpcap ! 3.27.2017 Cobb-Douglas production function's property (TAUCfinal.f90 line 1281). Note that wage, rd, and deltak are all on five years basis.
        ! 3.28.2017 Any tax rates do not show up in the equation above, as formualted in Imorhoroglu (1995) a life-cycle analysis of social security.
        
        gdp = crpprd + entprd ! 5.10.2017 As in  Fern-Villa and Dirk Krueger, GDP doens't include housing service in their definition.
        ! 8-18-2017 The definition of GDP is computed in the same way as in Cagetti and De Nardi AER_2009 paper.

        ! Cagetti: gdp = [wage*totlcorp+(rbar+delt)*totkcorp] + inck.
        ! Cagetti: `inck` (gross entr. output) includes wages paid, interest owed and depr. They collectively equal to the first term of my enrepreneur's profit (that is, the entrepreneurial production)
        
        ! 3.27.2017 corporate tax and capital income tax needs to be added into the tax revenue. Stop here. 4:09 pm <---------------
        !sw_worker_savtax, sw_entpre_savtax, sw_entpre_biztax ! There is no business taxes right?!
        
        totsvt = sum(sef*sw_worker_savtax, sef>0._wp .and. sw_worker_savtax>0._wp) + sum(sef*sw_entpre_savtax, sef>0._wp .and. sw_entpre_savtax>0._wp)
        !totbpt = sum(sef*sw_entpre_biztax) ! 4.16.2017 comment out
        !govbal = (totsvt + totbpt + tottax) - gfrac*gdp - rd*dfrac*(crpcap+entcap) ! 3.27.2017 excluding SS benefit, because it runs in the way like pay as you go.

        if(mode6taskid==0)then    
            govbal = (totsvt + tottax) - gfrac*gdp - rd*dfrac*(crpcap+entcap) ! 4.16.2017 TAUCfinal.f90 line 1280 + 1445.
            tottaxrev = totsvt + tottax ! used for mpi_exercise_mode=6
            !tauwealth = totsvt/sum(sef*sw_wealth_tax, sw_wealth_tax>0._wp .and. sef>0._wp) ! don't be confused with the name sw_wealth_tax. In the case of mode6taskid=0, it is the household net worth.
        elseif(mode6taskid==1)then
            govbal = sum(sef*sw_wealth_tax, sef>0._wp .and. sw_wealth_tax>0._wp) - tottaxrev ! 2-24-2018 tausv is set to zero, so we don't need sw_XXX_savtax here (and totsvt).
        elseif(mode6taskid==2)then    
            ! needs to be revised 2-24-2018
            govbal = (sum(sef*sw_wealth_tax, sef>0._wp .and. sw_wealth_tax>0._wp) + tottax + totsvt) - tottaxrev
        else
            print*, "error in govbal"
        endif !mode6taskid
        
        ! r_govbal = nint(govbal*10e14_wp)/10e14_wp
        !write(4000+trial_id,fmt='("macro-0",8x,11(e21.14,x))') totsvt, tottax, gfrac, gdp, rd, dfrac, crpcap, entcap, (totsvt + tottax), gfrac*gdp, rd*dfrac*(crpcap+entcap)
        !write(4000+trial_id,fmt='("macro-1",8x,9(e21.14,x))') totefl, entlab, crplab, tottax, entprd, crpprd, gdp, govbal !, r_govbal
        
        ! 3.27.2017 Note that Cagetti and De Nardi regard "SS transfer" as an item of government expenses and therefore formulated in the equation above.
        ! 3.27.2017 But in my model, Social security system runs in the way like pay as you go and therefore it is self-financed from the socail security tax.
        ! 3.27.2017 Therefore, in my formulation of variable "govbal" I don't include the revenue from socail security tax.
        govbal2gdp = govbal/gdp
        
        hug_inv_proj = sum(sef,tvec==.true..and.s3c(:,3)==3 .and. sef>=0._wp)
        med_inv_proj = sum(sef,tvec==.true..and.s3c(:,3)==2 .and. sef>=0._wp)
        sml_inv_proj = sum(sef,tvec==.true..and.s3c(:,3)==1 .and. sef>=0._wp)   
        
        hug_inv_proj = merge(hug_inv_proj, 0._wp, hug_inv_proj>=0._wp)
        med_inv_proj = merge(med_inv_proj, 0._wp, med_inv_proj>=0._wp)
        sml_inv_proj = merge(sml_inv_proj, 0._wp, sml_inv_proj>=0._wp)
        
        all_inv_proj = hug_inv_proj + med_inv_proj + sml_inv_proj
        
        if(all_inv_proj>=0._wp)then
            hug_inv_per  = hug_inv_proj/all_inv_proj
            med_inv_per  = med_inv_proj/all_inv_proj
            sml_inv_per  = sml_inv_proj/all_inv_proj
        else
            exit_log1 = .true.
            msg = ' non-positive entrepreneur mass '
        endif
        
        ent_wel_per  = (entfin+enthom)/(wokfin+entfin+wokhom+enthom)
        
        deallocate( svec, tvec )      
          
        !if(exit_log1==.false.)then ! There is no warning flag raised. 3.29.2017 stop here. ! 4.22.2017 comment out so that I know the program has something wrong.
            
        ! 3.28.2017 group size at the "END" of period
        ! 3.29.2017 From today onward, the worker cohort should include retired (because it is rare for entrepreneurs to retire in their retirement ages, so most of retirees are just pre-working-household household head.
        ! 3.29.2017 but it should EXCLUDE the non-in-labor-force people. and only look at periods 1 to 9.
        
        ! 5.10.2017 cross sectional size (depends on agent's decision)
        ! correct 5.9.2017
        ! All workers (including entrepreneur-turned workers). Excluding retirees.
        !woksize = sum(sef,s3c(:,3)==0.and.s3c(:,9)<=9) + sum(sef,s3c(:,3)/=0.and.s3c(:,4)==0.and.s3c(:,9)<=9) ! 4.17.2017 self-chosen workers + entrepreneur-turned workers. No retirees are included. 
        !chgw2e  = sum(sef,s3c(:,3)==0.and.swk==1.and.s3c(:,9)<=9) + sum(sef,s3c(:,3)/=0.and.s3c(:,4)==0.and.swk==1.and.s3c(:,9)<=9) ! a subset of the above 4.17.2017 `swk` is generated in the subroutine `convert_2d_outcome_into_series.` with `coarse` flag.
        woksize = sum(sef, s3c(:,3)==0.and.s3c(:,9)<=9.and.sef>=0._wp) !10.14.2017
        chgw2e  = sum(sef, s3c(:,3)==0.and.s3c(:,9)<=9.and.swk==1.and.sef>=0._wp) !10.14.2017
        w2erat  = chgw2e/woksize
        
        !write(4000+trial_id,fmt='("macro-2",8x,8(e21.14,x))') govbal2gdp, hug_inv_per, med_inv_per, sml_inv_per, ent_wel_per, woksize, chgw2e, w2erat
        
        ! All entrepreneurs
        !entsize = sum(sef,s3c(:,3)>0.and.s3c(:,4)==1) ! 4.17.2017 original ent size including "bad" luck entrepreneurs. 
        !! Entrepreneurs who are hit by bad business shock
        !chge2w  = sum(sef,s3c(:,3)>0.and.s3c(:,4)==1.and.swk==0) ! 4.17.2017 a subset of the above ! 5.10.2017 remove s3c(:,9)<=9. Correct.
        entsize = sum(sef, s3c(:,3)>0.and.sef>=0._wp) !10.14.2017
        chge2w  = sum(sef, s3c(:,3)>0.and.swk==0.and.sef>=0._wp) !10.14.2017
        e2wrat  = chge2w/entsize
        
        ! 4.17.2017 aggregate stats --------------------------------------
        ! 5.9.2017 lvece and lvecw both are correct.
        allocate(lvece(sz),lvecw(sz))
        lvece = .false.
        lvecw = .false.
        
        !! 5.10.2017 cross sectional size (not depends on agent's decision). Do not consider retirees.
        !lvece = s3c(:,3)>0.and.s3c(:,4)==1.and.abs(sw_buzcap_notuse-0._wp)<1.e-7_wp
        !lvecw = (lvece==.false..and.s3c(:,9)<=9) ! 5.10.2017 remove this line: (s3c(:,3)==0.and.s3c(:,9)<=9).or. and s3c(:,3)/=0.and.s3c(:,4)==0
        
        !lvece = s3c(:,3)>0.and.s3c(:,4)==1 ! 9-15-2017
        !lvecw = (s3c(:,3)==0.and.s3c(:,9)<=9).or.(s3c(:,3)>0.and.s3c(:,4)==0.and.s3c(:,9)<=9) ! 9-15-2017
        
        !10.14.2017
        lvece = s3c(:,3)>0 ! 10.14.2017 including back luck entrepreneurs.
        lvecw = s3c(:,3)==0.and.s3c(:,9)<=9
        
        szwok = count(lvecw)
        allocate(svec(szwok),rvec(szwok))
        svec = 0._wp
        rvec = 0._wp
        ! Method 1 causes stack overflow in laptop. Ok on cluster.
        !svec = pack(sef,lvecw==.true.)
        !rvec = pack(sw_taxableincome,lvecw==.true.)

        ! Method 2 
        j = 0
        do i = 1, sz
            if(lvecw(i)==.true.)then
                j = j + 1    
                svec(j) = merge(sef(i), 0._wp, sef(i)>=0._wp)
                !rvec(j) = sw_taxableincome(i) !10.13.207 comment out
                rvec(j) = sw_totinc_bx(i) !10.13.2017
            endif
        enddo
        
        call weighted_percentile(rvec,svec,0.5_wp,medwokinc)
        call weighted_percentile(rvec,svec,0.2_wp,lowest_quintile_wokinc) ! 4.23.2017
        
        j = 0
        do i = 1, sz
            if(lvecw(i)==.true.)then
                j = j + 1
                svec(j) = merge(sef(i), 0._wp, sef(i)>=0._wp)
                rvec(j) = sw_ini_asset(i) + sw_ini_house(i)
            endif
        enddo
        call weighted_percentile(rvec,svec,0.5_wp,medwokwel)
        deallocate(svec,rvec)
        
        ! Extended.
        szent = count(lvece)
        allocate(svec(szent),rvec(szent))
        svec = 0._wp
        rvec = 0._wp
        j = 0 
        do i = 1, sz
            if(lvece(i)==.true.)then
                j = j + 1
                svec(j) = merge(sef(i), 0._wp, sef(i)>=0._wp)
                rvec(j) = sw_ini_asset(i) + sw_ini_house(i)
            endif
        enddo
        call weighted_percentile(rvec,svec,0.5_wp,medentwel)    
        !call weighted_percentile(rvec,svec,suprich_prop,medallinc)    
        
        deallocate(svec,rvec)
            
        szwok = count(sef>0._wp)
        allocate(svec(szwok), rvec(szwok))
        j = 0
        do i = 1, sz
            if(sef(i)>0._wp)then
                j = j + 1
                svec(j) = sef(i)
                rvec(j) = sw_net_worth(i)
            endif    
        enddo
        call weighted_percentile(rvec,svec,0.5_wp,medallinc)    
        deallocate(svec, rvec)
        
        med_wel_e2w = medentwel/medwokwel
        !write(4000+trial_id,fmt='("macro-3",8x,8(e21.14,x))') entsize, chge2w, e2wrat, medwokinc, lowest_quintile_wokinc, medwokwel, medentwel, med_wel_e2w
        
        entcsp = sum(sef*sw_consumption, lvece.and.sw_consumption>0._wp.and.sef>=0._wp) ! 4.17.2017 It excludes bad luck entrepreneurs.
        wokcsp = sum(sef*sw_consumption, lvecw.and.sw_consumption>0._wp.and.sef>=0._wp) ! 4.17.2017 It includes entrepreneur-turned workers and retirees.
        
        entinc = sum(sef*sw_totinc_bx, lvece.and.sw_totinc_bx>0._wp.and.sef>=0._wp) !10.13.2017 10.18.2017
        wokinc = sum(sef*sw_totinc_bx, lvecw.and.sw_totinc_bx>0._wp.and.sef>=0._wp) !10.13.2017 10.18.2017
        
        !all_income  = entinc + wokinc ! sum(sef*sw_taxableincome) !10.13.2017 comment out
        all_income  = sum(sef*sw_totinc_bx,sw_totinc_bx>0._wp.and.sef>=0._wp) !10.13.2017
        ent_inc_per = entinc/all_income
        
        entaxw = sum(sef*sw_aftertaxwealth, lvece.and.sef>=0._wp) ! <--- 9-15-2017
        wokaxw = sum(sef*sw_aftertaxwealth, lvecw.and.sef>=0._wp) ! <--- 9-15-2017
        
        ! 4.17.2017 stats per capita -------------------------------------
        entsize = sum(sef, lvece.and.sef>=0._wp)
        woksize = sum(sef, lvecw.and.sef>=0._wp)
        
        mean_entfin = entfin/entsize
        mean_wokfin = wokfin/woksize
        mean_enthom = enthom/entsize
        mean_wokhom = wokhom/woksize

        mean_entinc = entinc/entsize 
        mean_wokinc = wokinc/woksize ! 11.1.2017 Important!!! `mean_wokinc` is updated and will be used in equilibrium.f90's line 413. 
        mean_entcsp = entcsp/entsize
        mean_wokcsp = wokcsp/woksize
        mean_entaxw = entaxw/entsize
        mean_wokaxw = wokaxw/woksize
        
        !write(4000+trial_id,fmt='("macro-4",8x,8(e21.14,x),2(11x,i10,x))') wokinc, all_income, ent_inc_per, entaxw, wokaxw, entsize, woksize, sum(sw_aftertaxwealth), count(lvece), count(lvecw)
        
        ! 5.10.2017 Update zone 
        ! nakajima 2010, pdf p.34 ------------------------------------------------------------------------
        rimplied = alpha*(crpcap/crplab)**(alpha-1._wp)-deltak ! Note: further adjustment needed to divide itself by 5 year in equilibrium.f90.             
        sumsstax = sum(sef*sw_laborsupply,sef>=0._wp)*wage*tauss ! 4.17.2017 tauss/2._wp*wage*sw_laborsupply(idx)
        ! 4.17.2017 stop here 4:09 pm.
        !poppaysstaximplied = sumsstax/(tauss*wage) ! 3.25.2017 May need to be revised!! the base that determines the total amount of social security benefits. ! 4.17.2017 comment out.
        !benefit  = sumsstax/sum(popfrac(10:14)) ! 4.17.2017 update benefit.
        benefitimplied  = sumsstax/sum(popfrac(10:14)) ! 7-6-2017 On Mesabi, conducting update after a longer interval is faster in terms of convergence. So I use beenfitimplied and update only 
        
        deallocate(lvece,lvecw)
        
        !! 10.13.2017
        !if(printout19)then
        !    write(my_id+1001,'(/,a)') ' ======================================= original '
        !    write(my_id+1001,'(8(4x,a,x))') 'totast', 'wokfin', 'entfin', 'wokhom', 'enthom', 'entcap', 'crpcap', '   gdp' 
        !    write(my_id+1001,'(8(e10.3,x))') totast, wokfin, entfin, wokhom, enthom, entcap, crpcap, gdp
        !    write(my_id+1001,'(a)') ' ======================================= new      '
        !    write(my_id+1001,'(8(4x,a,x))') 'totast', 'wokfin', 'entfin', 'wokhom', 'enthom', 'entcap', 'crpcap', '   gdp' 
        !    write(my_id+1001,'(8(e10.3,x),/)') totast, sum(sef*sw_ini_asset,ttvec==.false.), sum(sef*sw_ini_asset,ttvec==.true.), sum(sef*sw_ini_house,ttvec==.false.), sum(sef*sw_ini_house,ttvec==.true.), &
        !        sum(sef*(sw_buzcap_notuse+sw_bizinvestment),ttvec==.true.), totast/(1._wp+dfrac)-sum(sef*(sw_buzcap_notuse+sw_bizinvestment),ttvec==.true.), & 
        !        (sum(sef*(sw_buzcap_notuse+sw_bizinvestment),ttvec==.true.) + totast/(1._wp+dfrac)-sum(sef*(sw_buzcap_notuse+sw_bizinvestment))) ! 10.13.2017
        !endif         
        
        ! 4.16.2017 stop here. 8:30 pm.
        ! negative value 3.27.2017 stop here. (1) ---------------------------------------------------------------------------------------------------
        if(crpcap<0._wp)then 
            exit_log1 = .true.
            msg = ' crpcap<0 '
        endif
        if(crplab<0._wp)then ! 3.27.2017 Cagetti and De Nardi have the same (the only one) statement too. (2)
            exit_log1 = .true.
            msg = ' crplab<0 '
        endif
        
        if(entlab<0._wp)then ! (3)
            exit_log1 = .true.
            msg = ' entlab<0 '
        endif
        if(entcap<0._wp)then ! (4)
            exit_log1 = .true.
            msg = ' entcap<0 '
        endif 
        
        !! Not a number (NaN) cases usually caused by division by zero, or overflow errors. ---------------------------------------------------------
        if(crpcap/=crpcap)then ! (1)
            exit_log1 = .true.
            msg = ' crpcap NaN '
        endif
        if(crplab/=crplab)then ! (2)
            exit_log1 = .true.
            msg = ' crplab NaN '
        endif
        if(entlab/=entlab)then ! (3)
            exit_log1 = .true.
            msg = ' entlab NaN '
        endif
        if(entcap/=entcap)then ! (4)
            exit_log1 = .true.
            msg = ' entcap NaN '
        endif        
        
        !if( sml_inv_per<0._wp .or. med_inv_per<0._wp .or. hug_inv_per<0._wp )then ! (5) 10.18.2017
        !    exit_log1 = .true.
        !    msg = ' negative mass '
        !endif
        
        if(exit_log1==.false.)then ! the steady state is obtained successfully.

            momvec(1)  = entcap/(crpcap+entcap) ! crpcap/(crpcap+entcap) ! 8-18-2017 continue using this expression.
            momvec(2)  = sml_inv_per
            momvec(3)  = med_inv_per
            momvec(4)  = hug_inv_per
            momvec(5)  = entsize !10.14.2017 Only entrepreneurs experiencing good business shocks.
            !momvec(6)  = (crpcap+entcap)/gdp !wokfin + entfin ~= (crpcap+entcap) + <government debt> = (crpcap+entcap)*(1+dfrac). 8-19-201 Correct. See Cagetti and DeNardi's 2009 AER code, line 2197 in InitialSS.f90.
            momvec(6)  = (wokfin + entfin)/gdp !10.18.2017
            momvec(7)  = (enthom+wokhom)/gdp !10.14.2017, 8-19-2017
            momvec(8)  = (entfin+enthom)/(entfin+enthom+wokfin+wokhom) ! total asset helds by entrepreneurs, regardless of business shocks 7-2-2017, including the amount belonging to retired people.
            momvec(9)  = ent_inc_per
            momvec(10) = medentwel/medwokwel !ok. median net worth.
            
            momvec = nint(momvec*momround)/momround ! 10-11-2017
            
            if(printout6.and.(num_procs==2.or.num_procs==1)) write(*,fmt='(a,3e15.7)') 'Distribution: ', sml_inv_per, med_inv_per, hug_inv_per
            
            ! Inspect whether all elements of momvec are available.
            do i = 1, size(momvec)
                if(momvec(i)/=momvec(i))then
                    exit_log1 = .true.
                    msg = 'at least one momvec is not available'
                endif
            enddo
            
        else ! Divergence happends. Bad input. No sensible result is obtained.
            momvec = penalty 
            write(my_id+1001,'(/,a,a,/)') '====== Exit message ====== ', msg
        endif 
        
        if(printout19)then          
            write(my_id+1001,'(/,a)') ' ======================================= macrostatistics for each government loop '
            write(my_id+1001,'(8(4x,a,x))') 'netast', 'wokfin', 'entfin', 'wokhom', 'enthom', 'entcap', 'crpcap', '   gdp' 
            write(my_id+1001,'(8(e10.3,x))') totast, wokfin, entfin, wokhom, enthom, entcap, crpcap, gdp
            !write(my_id+1001,'(a)') ' ======================================= new      '
            !write(my_id+1001,'(8(4x,a,x))') 'totast', 'wokfin', 'entfin', 'wokhom', 'enthom', 'entcap', 'crpcap', '   gdp' 
            !write(my_id+1001,'(8(e10.3,x))') totast, wokfin1, entfin1, wokhom1, enthom1, entcap1, crpcap1, gdp1
            !write(my_id+1001,'(a)') ' --------------------------------------- original '
            write(my_id+1001,'(5(4x,a,x),2(x,a,x),2(3x,a,x))') 'hugPrj', 'medPrj', 'smlPrj', 'w2eRat', 'e2wRat', 'entIncRat', 'medWelRat', 'entSize', 'wokSize' 
            write(my_id+1001,'(9(e10.3,x))') hug_inv_per, med_inv_per, sml_inv_per, w2erat, e2wrat, ent_inc_per, med_wel_e2w, entsize, woksize
            
            if(mode6taskid==0)then
                write(my_id+1001,'(5(2x,a,x),(x,a,x))') ' entCapR', 'asst2Gdp', 'home2Gdp', 'netAsstR', 'tottaxrv'
                write(my_id+1001,'(6(e10.3,x))') momvec(1), momvec(6), momvec(7), momvec(8), tottaxrev
            elseif(mode6taskid==1)then
                write(my_id+1001,'(5(2x,a,x),(x,a,x))') ' entCapR', 'asst2Gdp', 'home2Gdp', 'netAsstR', 'tottaxrv', 'txRwealth'
                write(my_id+1001,'(6(e10.3,x))') momvec(1), momvec(6), momvec(7), momvec(8), tottaxrev, taubal                
            elseif(mode6taskid==2)then
                write(my_id+1001,'(5(2x,a,x),(x,a,x))') ' entCapR', 'asst2Gdp', 'home2Gdp', 'netAsstR', 'tottaxrv', 'txRwealth'
                write(my_id+1001,'(6(e10.3,x))') momvec(1), momvec(6), momvec(7), momvec(8), tottaxrev, taubal                                
            endif
            
            write(my_id+1001,'(2(8x,a,x),(6x,a,x),(2x,a,x),2(3x,a,x))') 'rd', 'rl', 'wage', 'transbeq', 'avgincw', 'benefit'
            write(my_id+1001,'(6(e10.3,x))') rd, rl, wage, transbeq, avgincw, benefit
            write(my_id+1001,'(a,/)') ' ==================================== End of macrostatistics for each government loop '
            
            fsef = sef
            fhom = sw_ini_house
            fcsp = sw_consumption
            finc = sw_totinc_bx
            fast = sw_ini_asset
            fbuz = sw_bizinvestment
            faxw = sw_aftertaxwealth
            fnwr = sw_net_worth
            
        endif ! printout19        
        
        !write(4000+trial_id,fmt='("macro-5",8x,8(e21.14,x))') rimplied, sumsstax, benefitimplied, momvec(1:5)
        
        !write(4000+trial_id,fmt='("macro-6",8x,5(e21.14,x))') momvec(6:10)
        
        !!if(iteratot==1)then
        !    write(str1,fmt='(i4)') trial_id+5000
        !    write(str2,fmt='(i3.3)') idxt
        !    open(unit=5000+trial_id,file="output_"//str1//"-"//str2//"_taxelements.txt",action="write",status="replace")
        !    do i = 1, size(s3c(:,1))
        !        write(unit=5000+trial_id,fmt='(9e25.18)') sw_aftertaxwealth(i), sef(i), sw_aftertaxwealth(i)*sef(i)
        !    enddo ! i
        !    close(5000+trial_id)
        !!endif ! iteratot
        
    end subroutine macro_statistics  
    
    !!! 8-1-2017 version 1.
    !subroutine grid_boundary_inspection( amin, amax, mass_vec )
    !    implicit none
    !    real(wp), intent(inout) :: amin, amax
    !    integer :: sz_list, size_av, i
    !    integer :: valid_lower_limit, valid_upper_limit
    !    real(wp) :: mass_lower_limit, mass_upper_limit, temp, damin, damax
    !    real(wp), dimension(:), intent(out) :: mass_vec
    !    logical, dimension(:), allocatable :: tvec
    !    sz_list = size(s3c(:,1))
    !    size_av = size(av)
    !    
    !    damin = amin
    !    damax = amax
    !    
    !    allocate( tvec(sz_list) )
    !    
    !    ! To aggregate the corresponding mass for a specific level of financial asset holdings.
    !    do i = 1, size_av
    !        tvec = .false.
    !        tvec = s3c(:,1)==i
    !        mass_vec(i) = sum(sef,tvec)
    !    enddo
    !    
    !    valid_lower_limit = 1
    !    valid_upper_limit = size_av
    !    mass_lower_limit = mass_vec(valid_lower_limit)
    !    mass_upper_limit = mass_vec(valid_upper_limit)
    !    
    !    ! Lower bound. Stage 1. Shrink from the bottom of the array up until the mass of the new bottom grid is larger than a predetermined threshold (tinymass).
    !    do i = 1, size_av-1 ! 8-1-2017 should start here.
    !        if( mass_vec(i)<tinymass .and. mass_vec(i+1)>=tinymass )then ! the case that the grid may need to shrink.
    !            valid_lower_limit = i+1 ! update.
    !            mass_lower_limit = mass_vec(i+1)
    !            exit ! Exit the loop as long as we find the "first" grid with mass larger than a threshold and being preceeded by a grid with mass smaller than the threshold.
    !        endif
    !    enddo 
    !    amin = av( valid_lower_limit )
    !    
    !    ! Lower bound. Stage 2. Extend from the bottom of the array if the current bottom grid's mass is larger than a threshold (chunkmass).
    !    if( mass_lower_limit>chunkmass )then
    !        temp = 0.5_wp*abs( av(valid_lower_limit+1)-av(valid_lower_limit) )            
    !        amin = av(valid_lower_limit) - temp
    !    endif
    !
    !    ! Upper bound, Stage 1. Shrink from the top of the array down until the mass of the new top grid is larger than a predetermined threshold (tinymass).
    !    do i = size_av, 2, -1
    !        if( mass_vec(i)<tinymass .and. mass_vec(i-1)>=tinymass )then
    !            valid_upper_limit = i-1
    !            mass_upper_limit = mass_vec(i-1)
    !            exit
    !        endif
    !    enddo ! i
    !    amax = av( valid_upper_limit )
    !    
    !    ! Upper bound , Stage 2. Extend from the top of the array if the current top grid's mass is larer than a threshold (chunkmass).
    !    if( mass_upper_limit>chunkmass )then
    !        temp = 0.5_wp*abs( av(valid_upper_limit)-av(valid_upper_limit-1) )
    !        amax = av(valid_upper_limit) + temp
    !    endif 
    !    
    !    deallocate( tvec )
    !end subroutine grid_boundary_inspection
    
    subroutine grid_boundary_inspection( amin, amax, mass_vec, flg )
        implicit none
        real(wp), intent(inout) :: amin, amax
        real(wp), dimension(:), intent(out) :: mass_vec
        real(wp) :: temp
        integer :: size_list, size_av, i, size_idxvec, idx_max, idx_min
        integer, dimension(:), allocatable :: intvec, idxvec
        logical, dimension(:), allocatable :: logvec
        logical, optional, intent(inout) :: flg
        
        !write(*,'(a)') ' enter grid inspection '
        
        size_list = size(s3c(:,1)) 
        size_av = size(av)
        
        allocate( logvec(size_list), intvec(size_list) )
        intvec = [(i,i=1,size_list)]
        
        ! vector of a repeated group with members indexed in increasing order, [1,2,...,size_av] 8-28-2017
        intvec = mod(intvec,size_av)
        where(intvec==0) intvec=size_av
        
        logvec = .false.
        logvec = sef>tinymass
        size_idxvec = count(logvec)
        allocate( idxvec(size_idxvec) )
        idxvec = pack(intvec,logvec)
        ! search for the maximum index of financial asset holding

        idx_max = maxval(idxvec)
        ! search for the minimum index of financial asset holding
        idx_min = maxval(-idxvec)
        idx_min = -idx_min
        
        if( (idx_min<1 .or. idx_max<1) .or. (idx_min>adim .or. idx_max>adim) )then  !10.17.2017
            flg = .true.
            return
        endif
        
        ! return (to be updated)
        amax = av(idx_max)
        amin = av(idx_min)
        
        !! calculate distribution of financial asset holdings
        !mass_vec = 0._wp
        do i = 1, size_av
            mass_vec(i) = sum(sef,s3c(:,1)==i)            
        enddo ! i
        
        ! if having fat tail, expand the corresponding limit (by raising it if it is a upper bound, and vice versa.)
        if(mass_vec(idx_max)>chunkmass)then !10.17.2017
            if(idx_max>idx_min)then
                temp = 0.5_wp*abs( av(idx_max)-av(idx_max-1) )
            else
                temp = 0.5_wp*abs( av(idx_max) )
            endif
            amax = amax + temp ! increase up a little bit the upper bound
        endif
        
        if(mass_vec(idx_min)>chunkmass)then !10.17.2017
            if(idx_max>idx_min)then
                temp = 0.5_wp*abs( av(idx_min)-av(idx_min+1) )
            else
                temp = 0.5_wp*abs( av(idx_min) )
            endif
            amin = amin - temp ! lower down a little bit the lower bound
        endif
        
        deallocate( idxvec, logvec, intvec )
        
        !write(*,'(a)') ' exit grid inspection '
    end subroutine grid_boundary_inspection
    
    subroutine grid_housing_upper_bound_adjustment(new_hmax, h_mass_vec, flg)
        implicit none
        real(wp), intent(out) :: new_hmax
        real(wp), dimension(:), intent(out) :: h_mass_vec
        real(wp) :: temp
        integer :: i, size_list, size_hv, idx_max, size_idxvec
        integer, dimension(:), allocatable :: intvec, idxvec
        logical, dimension(:), allocatable :: logvec
        logical, optional, intent(inout) :: flg
        
        size_list = size(s3c(:,1))
        size_hv = size(hv)
        allocate(intvec(size_list),logvec(size_list))
        intvec = [(i,i=1,size_list)]
        intvec = mod(intvec,size_hv)
        where(intvec==0) intvec = size_hv
        
        logvec = .false.
        logvec = sef>tinymass
        size_idxvec = count(logvec)
        allocate( idxvec(size_idxvec) )
        idxvec = pack(intvec, logvec)
        idx_max = maxval(idxvec)
        
        if( idx_max<1 .or. idx_max>hdim )then
            flg = .true.
            return
        endif
        
        new_hmax = hv(idx_max)
        do i = 1, size_hv
            h_mass_vec(i) = sum(sef, s3c(:,2)==i)
        enddo
        if(h_mass_vec(idx_max)>chunkmass)then
            if(idx_max>1)then
                temp = 0.5_wp*abs( hv(idx_max)-hv(idx_max-1))
            else !idx_max==1
                temp = 0.5_wp*abs( hv(idx_max) )
            endif
            new_hmax = new_hmax + temp
        endif
        deallocate(intvec,logvec,idxvec)
    end subroutine grid_housing_upper_bound_adjustment
    
    subroutine read_series2series(vec,fname)
        implicit none
        integer :: iostat, m, n
        real(wp), dimension(:), intent(out) :: vec
        character(len=*), intent(in) :: fname
        real(wp) :: val
        m = 0
        n = size(vec)
        open(unit=107,file=fname,status='old',action='read',iostat=iostat)
            if(iostat==0)then
                do 
                    read(unit=107,fmt=*,iostat=iostat) val    
                    if(iostat/=0) exit
                    m = m + 1
                    vec(m) = val
                enddo
            else
                write(*,fmt='(a)') ' something wrong with read_series2series '
            endif
            if(m/=n) write(*,fmt='(a)') ' inconsistency in read_series2series '
        close(107)
    end subroutine read_series2series
    
    subroutine compute_lorenz()
        implicit none
        integer :: n, m
        logical, dimension(:), allocatable :: lvec
        real(wp), dimension(:), allocatable :: nvec, dvec
        
        !call read_series2series(sef,'sef.txt')
        !call read_series2series(sw_aftertaxwealth,'sw_aftertaxwealth.txt')
        !call read_series2series(sw_consumption,'sw_consumption.txt')
        m = size(sef)
        n = count(sef>0._wp)
        allocate( lvec(m), nvec(n), dvec(n) )
        lvec = sef>0._wp
        !dvec = pack(sef,sef>0._wp) ! numerical level
        !nvec = pack(sw_aftertaxwealth,sef>0._wp)
        dvec = pack(sef,lvec)
        allocate(axw_lorenz(n),csp_lorenz(n), xbi_lorenz(n))
        nvec = pack(sw_aftertaxwealth,lvec)
        call lorenz(dvec,nvec,axw_lorenz,axw_gini)
        nvec = pack(sw_consumption,lvec)
        call lorenz(dvec,nvec,csp_lorenz,csp_gini) ! Note: lorenz(f,x,fx,gini)
        !nvec = pack(sw_taxableincome,lvec)
        nvec = pack(sw_totinc_bx,lvec)
        call lorenz(dvec,nvec,xbi_lorenz,xbi_gini)
                
        !print*, ' axw_gini ', axw_gini
        !print*, ' csp_gini ', csp_gini
        !print*, ' xbi_gini ', xbi_gini
        !call ss(axw_lorenz,'axw_lorenz_y',20,12)
        !call ss(nvec,'axw_lorenz_x',20,12)
        deallocate( lvec, dvec, nvec, axw_lorenz, csp_lorenz, xbi_lorenz )
        
    end subroutine compute_lorenz
    
    ! 4.17.2017 Save for Stata age-profile plot with 1) mass 2) financial asset and 3) non-financial asset for each combination.
    subroutine printout_for_stata_lifecycle_plotting()
        implicit none
        integer :: q, idx
        do q = 1, size(s3c(:,1))
            ax  = s3c(q,1)
            hx  = s3c(q,2)
            kx  = s3c(q,3)
            zx  = s3c(q,4)
            yx  = s3c(q,5)
            kpx = s3c(q,6)
            ypx = s3c(q,7)
            opx = s3c(q,8)
            t   = s3c(q,9)
            idx = s3c(q,10)
            write(unit=121,fmt='(10(i6,",",x),f25.15,",",x,f25.15,",",x,f25.15)') ax,hx,kx,zx,yx,kpx,ypx,opx,t,idx,sef(idx),sw_ini_asset(idx),sw_ini_house(idx)
        enddo ! q
    end subroutine printout_for_stata_lifecycle_plotting
        
end module model