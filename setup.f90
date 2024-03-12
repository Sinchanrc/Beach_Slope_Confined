module setup

    use particle
    implicit none

    public
    
    contains

    subroutine init()

        use initialize

        implicit none

        integer :: i,j,k,m

        ! rhomax=1108.373!1090.0_dp!(1.0_dp+atwood)*rhomin/(1.0_dp-atwood)
        rhomax=rho*rel_den

        fpy=floor(real(wc,dp)/(2*real(pr,dp)))+1
        fpx=floor(real((wl),dp)/(2*real(pr,dp)))!+1
        prrealx= ((real(wl,dp)/(fpx))/2.0)!*(fpx-1)/fpx
        ! prrealy= ((wc/(fpy-1))/2.0)*(fpy-1)/fpy
        ! prrealx= ((real(wl,dp)/(fpx-1))/2.0_dp)
        prrealy= ((real(wc,dp)/(fpy-1))/2.0_dp)
        ! fpx=fpx-1
        fmass=((wc*wl*1.0_dp)*rho)/(fpx*fpy)
        domain_shift=4*prrealx


        spy=floor(real(soilh,dp)/(2*real(prrealy,dp)))+4
        spx=floor(real(soill,dp)/(2*real(prrealx,dp)))
        solidx=((soill/(spx))/2.0)!*(spx-1)/spx
        solidy=((real(soilh,dp)/(spy-1))/2.0_dp)

        L=L!+2*prrealx
        H=H+2*prrealy

        bny=floor((real(H,dp)/(2*real(br,dp))))!+1
        bnx=floor((real(L,dp)/(2*real(br,dp))))
        brrealx= ((real(L,dp)/(bnx))/2.0_dp)!*(bnx-1)/real(bnx,dp)
        brrealy= ((real(H,dp)/(bny))/2.0_dp)!*(bny-1)/real(bny,dp)
        ! bny=bny-1
        ! bnx=bnx-1
        ! h1=(max(prrealx,prrealy)*hfac)/(2.0_dp)   !4.8 *sqrt(por)
        h1=(max(prrealx,prrealy)*hfac)/(2.0_dp*sqrt(por)) 
        totc=max((floor(real(wc,dp)/(incr*h1)))*(floor(wl/(incr*h1))),1)
        fplist=ceiling(real((fpx*fpy),dp)/(totc)*1.0_dp) 
        fplistmax=ceiling(1.20_dp*real(fplist,dp))
        cellx= ceiling(real(L,dp)/(incr*h1))+4
        celly= ceiling(real(H,dp)/(incr*h1))+3     
        dl=(min(2*prrealx,2*prrealy))
        dl1=(min(prrealx,prrealy)+min(brrealx,brrealy))
        bny=floor((real(lhs_btm,dp)/(2*real(brrealy,dp))))
        domain_shift=5*prrealx


        bnx=bnx+2*(bl)
        bny=bny+2*(bl)

        ! allocate(blist(bny,bl),dpcell(celly,cellx),flist(fpy,fpx+5))

        xl=(2*bl-1)*brrealx+brrealx 
        yl=(2*bl-1)*brrealy+solidy+set_ht
        xu=xu+(2*bl-1)*brrealx+brrealx
        yu=yu+(2*bl-1)*brrealy+solidy+set_ht
        ! line_grad=(yu-yl)/(xu-xl)
        line_grad=(line_grad*22)/(180.0_dp*7)

        fpy=floor(real(wc,dp)/(2*real(prrealy,dp)/sqrt(por)))!+1
        fpx=floor(real((wl),dp)/(2*real(prrealx,dp)/sqrt(por)))

        allocate(blist(bny,bl),dpcell(celly,cellx),flist(fpy,fpx+5))

        xrcutoff=((brrealx)*((2*bl)-1))+(fpx-1)*2*prrealx/sqrt(por)+3*prrealx/sqrt(por)+domain_shift-0.001_dp
        xlcutoff=((brrealx)*((2*bl)-1))+brrealx+domain_shift!+0.001_dp
        ytcutoff=((brrealy)*((2*bl)-1))+(fpy-1)*2*prrealy/sqrt(por)-2*prrealy/sqrt(por)

        icount=0
        count=0
        lam=0.0010_dp*(h1**2)
        sx=2
        sy=2
        ex=cellx-1
        ey=celly-1
        con_fac=(rel_den-1.0_dp)**(-1) 
    

        !$omp parallel default(shared)

        !Allocating the cells containing the particles
        !$omp do private(i,j) schedule(runtime) collapse(2) 
            do j=1,cellx
                do i=1,celly
                allocate(dpcell(i,j)%plist(fplistmax),dpcell(i,j)%ftemp(fplistmax), &
                dpcell(i,j)%cellid(2),dpcell(i,j)%porlist(fplistmax))

                ! if ((dpcell(i,j)%ybot<=ytcutoff).and.(dpcell(i,j)%xleft<=xlcutoff)) then
                !     allocate(dpcell(i,j)%exitlist(ceiling(0.25*fplistmax)))
                !     dpcell(i,j)%exitbuff=.true.
                ! end if

                do k=1,fplistmax 
                    dpcell(i,j)%ftemp(k)%part=>dpcell(i,j)%plist(k)
                end do


                end do
            end do
        !$omp end do

        !Determining the spatial limts of the cells
        !$omp do private(i,j) schedule(runtime) collapse(2) 
            do j=1,cellx
                do i=1,celly
                dpcell(i,j)%cellid(1)=j
                dpcell(i,j)%cellid(2)=i
                
                dpcell(i,j)%ybot=(-i+celly-1)*(incr*h1)
                dpcell(i,j)%ytop=(-i+celly)*(incr*h1)
                dpcell(i,j)%xleft=(j-2)*(incr*h1)
                dpcell(i,j)%xright=(j-1)*(incr*h1)

                end do
            end do
        !$omp end do
        
        !$omp end parallel

        xlcutoff=((brrealx)*((2*bl)-1))+brrealx
        ytcutoff=((brrealy)*((2*bl)-1))+(fpy-1)*2*prrealy/sqrt(por)+prrealy/sqrt(por)

        ! co=10*sqrt(2*abs(g)*wc)

        ! iparm=0
        ! iparm(1)=0
        ! iparm(2)=3
        ! iparm(4)=61
        ! iparm(10)=13
        ! iparm(11)=1
        ! iparm(13)=1
        ! iparm(27)=1

        allocate(reserve%tank(500))
        reserve%si=reserve_par
        

        return
    end subroutine init
    
end module setup