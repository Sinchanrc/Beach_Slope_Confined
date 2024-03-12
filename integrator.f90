module integrator

    use initialize
    use particle
    use domain
    use kernel

    implicit none

    contains

    pure function vis(oden) result(dmu)

        real(dp),intent(in) :: oden
        real(dp),parameter :: tem=20.0_dp
        real(dp) :: dmu,Acn,Bcn,Sp

        Sp=oden-1000.0_dp
        Acn=1.474_dp*(1e-3)+1.5_dp*(1e-5)*tem-3.927_dp*(1e-8)*tem**2
        Bcn=1.073_dp*(1e-5)-8.5_dp*(1e-8)*tem+2.230_dp*(1e-10)*tem**2
        dmu=exp(-0.00379418_dp+(0.604129_dp/(139.18_dp+tem)))*(1e-3)* &
            (1.0_dp + Acn*Sp + Bcn*Sp**2)



    end function


    subroutine cellshiftalt
        implicit none  
        
        integer :: i,j,k

        ! Preparing particle transfers to surrounding cells    
        !$omp do private(i,j,k) schedule (runtime) collapse(2)
            do j=sx,ex
                do i=sy,ey
                dpcell(i,j)%temfct=0
                dpcell(i,j)%tn%count=0
                dpcell(i,j)%ts%count=0
                dpcell(i,j)%te%count=0
                dpcell(i,j)%tw%count=0
                dpcell(i,j)%tne%count=0
                dpcell(i,j)%tnw%count=0
                dpcell(i,j)%tse%count=0
                dpcell(i,j)%tsw%count=0
                dpcell(i,j)%elist=0


                if (dpcell(i,j)%ptot/=0) then

                    do k=1,dpcell(i,j)%ptot ! must be carried out serially
                    if (dpcell(i,j)%plist(k)%tid==3) then
                    associate(x1=>dpcell(i,j)%plist(k)%x,y1=>dpcell(i,j)%plist(k)%y)
                        if((x1>=dpcell(i,j)%xleft).and.(x1<dpcell(i,j)%xright)) then
                        if((y1>=dpcell(i,j)%ybot).and.(y1<dpcell(i,j)%ytop))then
                            if ((x1>xlcutoff).or.(dpcell(i,j)%plist(k)%buffer).or. &
                            ((x1<=xlcutoff).and.(y1<=ytcutoff))) then
                            dpcell(i,j)%temfct=dpcell(i,j)%temfct+1
                            dpcell(i,j)%ftemp(dpcell(i,j)%temfct)%part=>dpcell(i,j)%plist(k)
                            ! dpcell(i,j)%ftemp(dpcell(i,j)%temfct)%dead=merge(.false.,.true.,x1<=domain_ex)
                            ! dpcell(i,j)%ftemp(dpcell(i,j)%temfct)%buffer=merge(.false.,.true.,x1<=domain_en)
                            ! else
                            !     dpcell(i,j)%elist=dpcell(i,j)%elist+1
                            !     dpcell(i,j)%exitlist(dpcell(i,j)%elist)=dpcell(i,j)%plist(k)%pid
                                
                            end if
                        elseif(y1<dpcell(i,j)%ybot) then
                            if ((x1>xlcutoff).or.(dpcell(i,j)%plist(k)%buffer).or. &
                            ((x1<=xlcutoff).and.(y1<=ytcutoff))) then
                            dpcell(i,j)%ts%count=dpcell(i,j)%ts%count+1
                            dpcell(i,j)%ts%list(dpcell(i,j)%ts%count)=dpcell(i,j)%plist(k)
                            ! dpcell(i,j)%ts%list(dpcell(i,j)%ts%count)%dead=merge(.false.,.true.,x1<=domain_ex)
                            ! dpcell(i,j)%ts%list(dpcell(i,j)%ts%count)%buffer=merge(.false.,.true.,x1<=domain_en)
                            ! else
                            !     dpcell(i,j)%elist=dpcell(i,j)%elist+1
                            !     dpcell(i,j)%exitlist(dpcell(i,j)%elist)=dpcell(i,j)%plist(k)%pid
                            end if
                        elseif(y1>=dpcell(i,j)%ytop) then
                            if ((x1>xlcutoff).or.(dpcell(i,j)%plist(k)%buffer).or. &
                            ((x1<=xlcutoff).and.(y1<=ytcutoff))) then
                            dpcell(i,j)%tn%count=dpcell(i,j)%tn%count+1
                            dpcell(i,j)%tn%list(dpcell(i,j)%tn%count)=dpcell(i,j)%plist(k)
                            ! dpcell(i,j)%tn%list(dpcell(i,j)%tn%count)%dead=merge(.false.,.true.,x1<=domain_ex)
                            ! dpcell(i,j)%tn%list(dpcell(i,j)%tn%count)%buffer=merge(.false.,.true.,x1<=domain_en)
                            ! else
                            !     dpcell(i,j)%elist=dpcell(i,j)%elist+1
                            !     dpcell(i,j)%exitlist(dpcell(i,j)%elist)=dpcell(i,j)%plist(k)%pid
                            end if
                        end if
                        elseif(x1>=dpcell(i,j)%xright) then
                        if(y1>=dpcell(i,j)%ytop)then
                            if ((x1>xlcutoff).or.(dpcell(i,j)%plist(k)%buffer).or. &
                            ((x1<=xlcutoff).and.(y1<=ytcutoff))) then
                            dpcell(i,j)%tne%count=dpcell(i,j)%tne%count+1
                            dpcell(i,j)%tne%list(dpcell(i,j)%tne%count)=dpcell(i,j)%plist(k)
                            ! dpcell(i,j)%tne%list(dpcell(i,j)%tne%count)%dead=merge(.false.,.true.,x1<=domain_ex)
                            ! dpcell(i,j)%tne%list(dpcell(i,j)%tne%count)%buffer=merge(.false.,.true.,x1<=domain_en)
                            ! else
                            !     dpcell(i,j)%elist=dpcell(i,j)%elist+1
                            !     dpcell(i,j)%exitlist(dpcell(i,j)%elist)=dpcell(i,j)%plist(k)%pid
                            end if
                        elseif((y1>=dpcell(i,j)%ybot).and.(y1<dpcell(i,j)%ytop)) then
                            if ((x1>xlcutoff).or.(dpcell(i,j)%plist(k)%buffer).or. &
                            ((x1<=xlcutoff).and.(y1<=ytcutoff))) then
                            dpcell(i,j)%te%count=dpcell(i,j)%te%count+1
                            dpcell(i,j)%te%list(dpcell(i,j)%te%count)=dpcell(i,j)%plist(k)
                            ! dpcell(i,j)%te%list(dpcell(i,j)%te%count)%dead=merge(.false.,.true.,x1<=domain_ex)
                            ! dpcell(i,j)%te%list(dpcell(i,j)%te%count)%buffer=merge(.false.,.true.,x1<=domain_en)
                            ! else
                            !     dpcell(i,j)%elist=dpcell(i,j)%elist+1
                            !     dpcell(i,j)%exitlist(dpcell(i,j)%elist)=dpcell(i,j)%plist(k)%pid
                            end if
                        elseif(y1<dpcell(i,j)%ybot) then
                            if ((x1>xlcutoff).or.(dpcell(i,j)%plist(k)%buffer).or. &
                            ((x1<=xlcutoff).and.(y1<=ytcutoff))) then
                            dpcell(i,j)%tse%count=dpcell(i,j)%tse%count+1
                            dpcell(i,j)%tse%list(dpcell(i,j)%tse%count)=dpcell(i,j)%plist(k)
                            ! dpcell(i,j)%tse%list(dpcell(i,j)%tse%count)%dead=merge(.false.,.true.,x1<=domain_ex)
                            ! dpcell(i,j)%tse%list(dpcell(i,j)%tse%count)%buffer=merge(.false.,.true.,x1<=domain_en)
                            ! else
                            !     dpcell(i,j)%elist=dpcell(i,j)%elist+1
                            !     dpcell(i,j)%exitlist(dpcell(i,j)%elist)=dpcell(i,j)%plist(k)%pid
                            end if
                        end if
                        elseif (x1<dpcell(i,j)%xleft) then
                        if (y1>=dpcell(i,j)%ytop) then
                            if ((x1>xlcutoff).or.(dpcell(i,j)%plist(k)%buffer).or. &
                            ((x1<=xlcutoff).and.(y1<=ytcutoff))) then
                            dpcell(i,j)%tnw%count=dpcell(i,j)%tnw%count+1
                            dpcell(i,j)%tnw%list(dpcell(i,j)%tnw%count)=dpcell(i,j)%plist(k)
                            ! dpcell(i,j)%tnw%list(dpcell(i,j)%tnw%count)%dead=merge(.false.,.true.,x1<=domain_ex)
                            ! dpcell(i,j)%tnw%list(dpcell(i,j)%tnw%count)%buffer=merge(.false.,.true.,x1<=domain_en)
                            ! else
                            !     dpcell(i,j)%elist=dpcell(i,j)%elist+1
                            !     dpcell(i,j)%exitlist(dpcell(i,j)%elist)=dpcell(i,j)%plist(k)%pid
                            end if
                        elseif((y1>=dpcell(i,j)%ybot).and.(y1<dpcell(i,j)%ytop)) then
                            if ((x1>xlcutoff).or.(dpcell(i,j)%plist(k)%buffer).or. &
                            ((x1<=xlcutoff).and.(y1<=ytcutoff))) then
                            dpcell(i,j)%tw%count=dpcell(i,j)%tw%count+1
                            dpcell(i,j)%tw%list(dpcell(i,j)%tw%count)=dpcell(i,j)%plist(k)
                            ! dpcell(i,j)%tw%list(dpcell(i,j)%tw%count)%dead=merge(.false.,.true.,x1<=domain_ex)
                            ! dpcell(i,j)%tw%list(dpcell(i,j)%tw%count)%buffer=merge(.false.,.true.,x1<=domain_en)
                            ! else
                            !     dpcell(i,j)%elist=dpcell(i,j)%elist+1
                            !     dpcell(i,j)%exitlist(dpcell(i,j)%elist)=dpcell(i,j)%plist(k)%pid
                            end if
                        elseif (y1<dpcell(i,j)%ybot) then
                            if ((x1>xlcutoff).or.(dpcell(i,j)%plist(k)%buffer).or. &
                            ((x1<=xlcutoff).and.(y1<=ytcutoff))) then
                            dpcell(i,j)%tsw%count=dpcell(i,j)%tsw%count+1
                            dpcell(i,j)%tsw%list(dpcell(i,j)%tsw%count)=dpcell(i,j)%plist(k)
                            ! dpcell(i,j)%tsw%list(dpcell(i,j)%tsw%count)%dead=merge(.false.,.true.,x1<=domain_ex)
                            ! dpcell(i,j)%tsw%list(dpcell(i,j)%tsw%count)%buffer=merge(.false.,.true.,x1<=domain_en)
                            ! else
                            !     dpcell(i,j)%elist=dpcell(i,j)%elist+1
                            !     dpcell(i,j)%exitlist(dpcell(i,j)%elist)=dpcell(i,j)%plist(k)%pid
                            end if
                        end if
                        end if
                    end associate
                    end if
                    end do

                end if
                end do
            end do
        !$omp end do

        ! Particle transfers    
        !$omp do private(i,j,k) schedule (runtime) collapse(2)
            do j=sx,ex
                do i=sy,ey
                ! From here everything must be done sequentially

                dpcell(i,j)%ptot=dpcell(i,j)%btot+dpcell(i,j)%temfct

                if (dpcell(i,j)%temfct/=0) then
                    do k=1,dpcell(i,j)%temfct
                    dpcell(i,j)%plist(k+dpcell(i,j)%btot)=dpcell(i,j)%ftemp(k)%part
                    end do
                end if

                if ((dpcell(i+1,j)%tn%count/=0)) then
                    do k=1,dpcell(i+1,j)%tn%count
                    dpcell(i,j)%plist(k+dpcell(i,j)%ptot)=dpcell(i+1,j)%tn%list(k)
                    end do
                end if
                dpcell(i,j)%ptot=dpcell(i,j)%ptot+(dpcell(i+1,j)%tn%count)

                if (dpcell(i-1,j)%ts%count/=0) then
                    do k=1,dpcell(i-1,j)%ts%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i-1,j)%ts%list(k)
                    end do
                end if

                if (dpcell(i,j+1)%tw%count/=0) then
                    do k=1,dpcell(i,j+1)%tw%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i,j+1)%tw%list(k)
                    end do
                end if

                if (dpcell(i,j-1)%te%count/=0) then
                    do k=1,dpcell(i,j-1)%te%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i,j-1)%te%list(k)
                    end do
                end if

                if (dpcell(i-1,j-1)%tse%count/=0) then
                    do k=1,dpcell(i-1,j-1)%tse%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i-1,j-1)%tse%list(k)
                    end do
                end if

                if (dpcell(i+1,j-1)%tne%count/=0) then
                    do k=1,dpcell(i+1,j-1)%tne%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i+1,j-1)%tne%list(k)
                    end do
                end if

                if (dpcell(i-1,j+1)%tsw%count/=0) then
                    do k=1,dpcell(i-1,j+1)%tsw%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i-1,j+1)%tsw%list(k)
                    end do
                end if

                if (dpcell(i+1,j+1)%tnw%count/=0) then
                    do k=1,dpcell(i+1,j+1)%tnw%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i+1,j+1)%tnw%list(k)
                    end do
                end if

                end do
            end do
        !$omp end do
        
    end subroutine cellshiftalt

    subroutine cellshift
        implicit none 
        
        integer :: i,j,k
        
        ! Preparing particle transfers to surrounding cells    
        !$omp do private(i,j,k) schedule (runtime) collapse(2)
            do j=sx,ex
                do i=sy,ey
                dpcell(i,j)%temfct=0
                dpcell(i,j)%tn%count=0
                dpcell(i,j)%ts%count=0
                dpcell(i,j)%te%count=0
                dpcell(i,j)%tw%count=0
                dpcell(i,j)%tne%count=0
                dpcell(i,j)%tnw%count=0
                dpcell(i,j)%tse%count=0
                dpcell(i,j)%tsw%count=0

                if (dpcell(i,j)%ptot/=0) then

                    do k=1,dpcell(i,j)%ptot ! must be carried out serially
                    if (dpcell(i,j)%plist(k)%tid==3) then
                    associate(x1=>dpcell(i,j)%plist(k)%x,y1=>dpcell(i,j)%plist(k)%y)
                        if((x1>=dpcell(i,j)%xleft).and.(x1<dpcell(i,j)%xright)) then
                        if((y1>=dpcell(i,j)%ybot).and.(y1<dpcell(i,j)%ytop))then
                            dpcell(i,j)%temfct=dpcell(i,j)%temfct+1
                            dpcell(i,j)%ftemp(dpcell(i,j)%temfct)%part=>dpcell(i,j)%plist(k)
                        elseif(y1<dpcell(i,j)%ybot) then
                            dpcell(i,j)%ts%count=dpcell(i,j)%ts%count+1
                            dpcell(i,j)%ts%list(dpcell(i,j)%ts%count)=dpcell(i,j)%plist(k)
                        elseif(y1>=dpcell(i,j)%ytop) then
                            dpcell(i,j)%tn%count=dpcell(i,j)%tn%count+1
                            dpcell(i,j)%tn%list(dpcell(i,j)%tn%count)=dpcell(i,j)%plist(k)
                        end if
                        elseif(x1>=dpcell(i,j)%xright) then
                        if(y1>=dpcell(i,j)%ytop)then
                            dpcell(i,j)%tne%count=dpcell(i,j)%tne%count+1
                            dpcell(i,j)%tne%list(dpcell(i,j)%tne%count)=dpcell(i,j)%plist(k)
                        elseif((y1>=dpcell(i,j)%ybot).and.(y1<dpcell(i,j)%ytop)) then
                            dpcell(i,j)%te%count=dpcell(i,j)%te%count+1
                            dpcell(i,j)%te%list(dpcell(i,j)%te%count)=dpcell(i,j)%plist(k)
                        elseif(y1<dpcell(i,j)%ybot) then
                            dpcell(i,j)%tse%count=dpcell(i,j)%tse%count+1
                            dpcell(i,j)%tse%list(dpcell(i,j)%tse%count)=dpcell(i,j)%plist(k)
                        end if
                        elseif (x1<dpcell(i,j)%xleft) then
                        if (y1>=dpcell(i,j)%ytop) then
                            dpcell(i,j)%tnw%count=dpcell(i,j)%tnw%count+1
                            dpcell(i,j)%tnw%list(dpcell(i,j)%tnw%count)=dpcell(i,j)%plist(k)
                        elseif((y1>=dpcell(i,j)%ybot).and.(y1<dpcell(i,j)%ytop)) then
                            dpcell(i,j)%tw%count=dpcell(i,j)%tw%count+1
                            dpcell(i,j)%tw%list(dpcell(i,j)%tw%count)=dpcell(i,j)%plist(k)
                        elseif (y1<dpcell(i,j)%ybot) then
                            dpcell(i,j)%tsw%count=dpcell(i,j)%tsw%count+1
                            dpcell(i,j)%tsw%list(dpcell(i,j)%tsw%count)=dpcell(i,j)%plist(k)
                        end if
                        end if
                    end associate
                    end if
                    end do

                end if
                end do
            end do
        !$omp end do

        ! Particle transfers    
        !$omp do private(i,j,k) schedule (runtime) collapse(2)
            do j=sx,ex
                do i=sy,ey
                ! From here everything must be done sequentially

                dpcell(i,j)%ptot=dpcell(i,j)%btot+dpcell(i,j)%temfct

                if (dpcell(i,j)%temfct/=0) then
                    do k=1,dpcell(i,j)%temfct
                    dpcell(i,j)%plist(k+dpcell(i,j)%btot)=dpcell(i,j)%ftemp(k)%part
                    end do
                end if

                if ((dpcell(i+1,j)%tn%count/=0)) then
                    do k=1,dpcell(i+1,j)%tn%count
                    dpcell(i,j)%plist(k+dpcell(i,j)%ptot)=dpcell(i+1,j)%tn%list(k)
                    end do
                end if
                dpcell(i,j)%ptot=dpcell(i,j)%ptot+(dpcell(i+1,j)%tn%count)

                if (dpcell(i-1,j)%ts%count/=0) then
                    do k=1,dpcell(i-1,j)%ts%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i-1,j)%ts%list(k)
                    end do
                end if

                if (dpcell(i,j+1)%tw%count/=0) then
                    do k=1,dpcell(i,j+1)%tw%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i,j+1)%tw%list(k)
                    end do
                end if

                if (dpcell(i,j-1)%te%count/=0) then
                    do k=1,dpcell(i,j-1)%te%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i,j-1)%te%list(k)
                    end do
                end if

                if (dpcell(i-1,j-1)%tse%count/=0) then
                    do k=1,dpcell(i-1,j-1)%tse%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i-1,j-1)%tse%list(k)
                    end do
                end if

                if (dpcell(i+1,j-1)%tne%count/=0) then
                    do k=1,dpcell(i+1,j-1)%tne%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i+1,j-1)%tne%list(k)
                    end do
                end if

                if (dpcell(i-1,j+1)%tsw%count/=0) then
                    do k=1,dpcell(i-1,j+1)%tsw%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i-1,j+1)%tsw%list(k)
                    end do
                end if

                if (dpcell(i+1,j+1)%tnw%count/=0) then
                    do k=1,dpcell(i+1,j+1)%tnw%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i+1,j+1)%tnw%list(k)
                    end do
                end if

                end do
            end do
        !$omp end do
        
    end subroutine cellshift

    subroutine boun_vel()
        implicit none

        logical :: pcounter

        real(dp) :: t1

        integer :: i,j,k,m

        !$omp do schedule (runtime) private(m,t1,k,i,j,pcounter) collapse(2)
        do j=sx,ex
            do i=sy,ey 

            do k=1,dpcell(i,j)%ptot

            if (dpcell(i,j)%plist(k)%tid<=2) then
            pcounter=.false.
            t1=0.0_dp
            dpcell(i,j)%plist(k)%vx=0.0_dp
            dpcell(i,j)%plist(k)%vy=0.0_dp
            dpcell(i,j)%plist(k)%vxs=0.0_dp
            dpcell(i,j)%plist(k)%vys=0.0_dp

            do m=1,dpcell(i,j)%list(k)%count
                associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                    y=>dpcell(i,j)%list(k)%pnh(m)%ppart)

                if (x%tid==3) then

                pcounter=.true.

                dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx- &
                    x%mass*Wab(dpcell(i,j)%list(k)%dist(m),h1)* &
                    x%vx/x%density

                dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy- &
                    x%mass*Wab(dpcell(i,j)%list(k)%dist(m),h1)* &
                    x%vy/x%density

                t1=t1+x%mass*Wab(dpcell(i,j)%list(k)%dist(m),h1)/&
                        x%density

                end if

                end associate

            end do
            if (t1>=1e-6) then
            dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx/t1
            ! dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vx
            dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy/t1            
            ! dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vy

            end if

            end if

            end do

            end do
        end do
        !$omp end do
        
    end subroutine boun_vel

    subroutine int_vel
        implicit none

        real(dp) :: t1,t2,pvol,p_dist,vart,ps=1.0_dp

        integer :: i,j,k,m

        vart=0.0_dp

        !$omp do schedule (runtime) private(m,t1,t2,k,i,j) collapse(2)
        do j=sx,ex
            do i=sy,ey 

            do k=1,dpcell(i,j)%ptot

            t1=0.0_dp
            t2=0.0_dp
            if(dpcell(i,j)%plist(k)%free) then
            dpcell(i,j)%pplist(k)%varts=0.0_dp

            do m=1,dpcell(i,j)%list(k)%count
                associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                    y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                    z=>dpcell(i,j)%list(k)%klt)


                t1=t1+(x%vx/y%porosity-dpcell(i,j)%plist(k)%vx/dpcell(i,j)%pplist(k)%porosity)* &
                (dpcell(i,j)%pplist(k)%coff(1)*z(1,m)+dpcell(i,j)%pplist(k)%coff(2)* &
                    z(2,m))*(x%mass/x%density)

                t2=t2+(x%vy/y%porosity-dpcell(i,j)%plist(k)%vy/dpcell(i,j)%pplist(k)%porosity)* &
                (dpcell(i,j)%pplist(k)%coff(3)*z(1,m)+dpcell(i,j)%pplist(k)%coff(4)* &
                    z(2,m))*(x%mass/x%density)

                end associate
            end do

            dpcell(i,j)%pplist(k)%varts=ps*(dl**2)*sqrt(t1**2+t2**2)*&
            sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)/ &
                ((umax+1e-6)*dpcell(i,j)%pplist(k)%porosity)


            end if

            end do

            end do
        end do
        !$omp end do

        ! Non-pressure velocity calculation for fluid particles
        !$omp do schedule (runtime) private(m,t1,t2,k,i,j,pvol,p_dist,vart) collapse(2)
        do j=sx,ex
            do i=sy,ey            
            
            do k=1,dpcell(i,j)%ptot

            t1=0.0_dp
            t2=0.0_dp
            if((dpcell(i,j)%plist(k)%tid==3).and. &
                (.not.(dpcell(i,j)%plist(k)%buffer))) then
                dpcell(i,j)%plist(k)%vxs=0.0_dp
                dpcell(i,j)%plist(k)%vys=0.0_dp
                if (dpcell(i,j)%list(k)%count/=0) then

                do m=1,dpcell(i,j)%list(k)%count
                    associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                        y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                        z=>dpcell(i,j)%list(k)%klt)

                        if (dpcell(i,j)%plist(k)%free) then
                            vart=dpcell(i,j)%pplist(k)%varts
                            elseif ((.not.(dpcell(i,j)%plist(k)%free)).and. &
                            (x%free).and.(.not.(dpcell(i,j)%plist(k)%buffer))) then
                            vart=Wab(dpcell(i,j)%list(k)%dist(m),h1)*y%varts &
                                    /Wo(h1)
            
                            else
                            vart=0.0_dp
                        end if
            

                    pvol=x%mass/x%density
                    p_dist=(dpcell(i,j)%list(k)%dist(m)**2+lam)**(-1)

                        t1=(2*pvol*(dpcell(i,j)%plist(k)%vx- &
                        x%vx)*(dpcell(i,j)%plist(k)%x- &
                        x%x)* &
                        (dpcell(i,j)%pplist(k)%coff(1)*z(1,m)+dpcell(i,j)%pplist(k)%coff(2)* &
                        z(2,m)))*p_dist


                        t2=(2*pvol*(dpcell(i,j)%plist(k)%vx- &
                        x%vx)*(dpcell(i,j)%plist(k)%y- &
                        x%y)* &
                        (dpcell(i,j)%pplist(k)%coff(3)*z(1,m)+dpcell(i,j)%pplist(k)%coff(4)* &
                        z(2,m)))*p_dist

                        dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs+ &
                        (t1+t2)*((2*(mu*dpcell(i,j)%pplist(k)%porosity/dpcell(i,j)%plist(k)%density)* &
                        (mu*y%porosity/x%density)/ &
                        ((mu*dpcell(i,j)%pplist(k)%porosity/dpcell(i,j)%plist(k)%density)&
                        +(mu*y%porosity/x%density)))+vart) !+vart

                        if (x%tid>2) then

                        if ((dpcell(i,j)%pplist(k)%nut<1e-6).and.(y%nut<1e-6)) then
                            
                        dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs+ &
                        (t1+t2)*(dpcell(i,j)%pplist(k)%nut+y%nut)*0.50_dp
                        else 
                            dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs+ &
                        (t1+t2)*((2*dpcell(i,j)%pplist(k)%nut*y%nut)/ &
                        (dpcell(i,j)%pplist(k)%nut+y%nut))

                        end if

                        ! else

                        !     dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs+ &
                        ! (t1+t2)*(dpcell(i,j)%pplist(k)%nut)

                        end if


                        t1=(2*pvol*(dpcell(i,j)%plist(k)%vy- &
                        x%vy)*(dpcell(i,j)%plist(k)%x- &
                        x%x)* &
                        (dpcell(i,j)%pplist(k)%coff(1)*z(1,m)+dpcell(i,j)%pplist(k)%coff(2)* &
                        z(2,m)))*p_dist


                        t2=(2*pvol*(dpcell(i,j)%plist(k)%vy- &
                        x%vy)*(dpcell(i,j)%plist(k)%y- &
                        x%y)* &
                        (dpcell(i,j)%pplist(k)%coff(3)*z(1,m)+dpcell(i,j)%pplist(k)%coff(4)* &
                        z(2,m)))*p_dist

                        dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys+ &
                        (t1+t2)*((2*(mu*dpcell(i,j)%pplist(k)%porosity/dpcell(i,j)%plist(k)%density)* &
                        (mu*y%porosity/x%density)/ &
                        ((mu*dpcell(i,j)%pplist(k)%porosity/dpcell(i,j)%plist(k)%density)&
                        +(mu*y%porosity/x%density)))+vart) !+vart

                        if (x%tid>2) then

                        if ((dpcell(i,j)%pplist(k)%nut<1e-6).and.(y%nut<1e-6)) then
                            
                            dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys+ &
                            (t1+t2)*(dpcell(i,j)%pplist(k)%nut+y%nut)*0.50_dp
                        else 
                                dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys+ &
                            (t1+t2)*((2*dpcell(i,j)%pplist(k)%nut*y%nut)/ &
                            (dpcell(i,j)%pplist(k)%nut+y%nut))
        
                        end if

                        ! else

                        !     dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys+ &
                        !     (t1+t2)*(dpcell(i,j)%pplist(k)%nut)

                        end if

                        ! Turbulence Calculations

                        if (dpcell(i,j)%plist(k)%tid==3) then

                        t1=(pvol*(y%tke- &
                        dpcell(i,j)%pplist(k)%tke)* &
                        (dpcell(i,j)%pplist(k)%coff(1)*z(1,m)+dpcell(i,j)%pplist(k)%coff(2)* &
                        z(2,m)))

                        ! dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs+(t1+t2)*dpcell(i,j)%pplist(k)%porosity
                        dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs- &
                        (t1)*dpcell(i,j)%pplist(k)%porosity*2.0_dp/3


                        t1=(pvol*(y%tke- &
                        dpcell(i,j)%pplist(k)%tke)* &
                        (dpcell(i,j)%pplist(k)%coff(3)*z(1,m)+dpcell(i,j)%pplist(k)%coff(4)* &
                        z(2,m)))

                        ! dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys+(t1+t2)*dpcell(i,j)%pplist(k)%porosity
                        dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys-&
                        (t1)*dpcell(i,j)%pplist(k)%porosity*2.0_dp/3

                        end if

                        ! end if

                    end associate
                end do
                end if
                    ! dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vx+dt* &
                    ! dpcell(i,j)%plist(k)%vxs+dt*dpcell(i,j)%pplist(k)%resistx*dpcell(i,j)%pplist(k)%porosity

                    ! dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vy+dt* &
                    ! (dpcell(i,j)%plist(k)%vys+g*dpcell(i,j)%pplist(k)%porosity)+ &
                    ! dt*dpcell(i,j)%pplist(k)%resisty*dpcell(i,j)%pplist(k)%porosity  
                    
                    dpcell(i,j)%plist(k)%vxs=(dpcell(i,j)%plist(k)%vx+(dt* &
                    dpcell(i,j)%plist(k)%vxs))/(1.0_dp-dt*dpcell(i,j)%pplist(k)%resistx*dpcell(i,j)%pplist(k)%porosity)

                    dpcell(i,j)%plist(k)%vys=(dpcell(i,j)%plist(k)%vy+dt* &
                    (dpcell(i,j)%plist(k)%vys))/(1.0_dp-dt*dpcell(i,j)%pplist(k)%resisty*dpcell(i,j)%pplist(k)%porosity)

                    dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys+g*dpcell(i,j)%pplist(k)%porosity*dt

            end if

            if ((dpcell(i,j)%plist(k)%buffer)) then

                dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vx
                ! dpcell(i,j)%plist(k)%vys=g*dpcell(i,j)%pplist(k)%porosity*dt
            end if

                end do
            end do
        end do
        !$omp end do
        
    end subroutine int_vel

    subroutine comp_vel
        implicit none

        integer :: i,j,k,m
        
        real(dp) :: t1,t2

        ! New velocity calculations for particles 
        !$omp do private(m,t1,i,k,j,t2) schedule (runtime) collapse(2)
        do j=sx,ex
            do i=sy,ey
                if(dpcell(i,j)%ptot/=0) then
                
                do k=1,dpcell(i,j)%ptot
                    if ((dpcell(i,j)%plist(k)%tid==3).and. &
                    (.not.(dpcell(i,j)%plist(k)%buffer))) then
                    !New fluid particle velocities
                    dpcell(i,j)%plist(k)%vx=0.0_dp
                    dpcell(i,j)%plist(k)%vy=0.0_dp
                    if (dpcell(i,j)%list(k)%count/=0) then

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                            y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                            z=>dpcell(i,j)%list(k)%klt)

                        t2=2*dpcell(i,j)%pplist(k)%porosity*y%porosity/ &
                        (y%porosity+dpcell(i,j)%pplist(k)%porosity)

                        ! t1=(dpcell(i,j)%pplist(k)%porosity**2)*(x%mass/x%density)* &
                        !     (x%pressure-dpcell(i,j)%plist(k)%pressure)* &
                        !     (dpcell(i,j)%pplist(k)%coff(1)*z(1,m)+dpcell(i,j)%pplist(k)%coff(2)* &
                        ! z(2,m))/dpcell(i,j)%plist(k)%density

                        ! t1=((x%mass/x%density)*(dpcell(i,j)%pplist(k)%porosity**2)*&
                        ! (x%pressure+dpcell(i,j)%plist(k)%pressure)* &
                        ! (z(1,m)))/dpcell(i,j)%plist(k)%density

                        t1=2*((x%mass/x%density)*(dpcell(i,j)%pplist(k)%porosity**2)*&
                        ((x%pressure*dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)+&
                        (dpcell(i,j)%plist(k)%pressure*x%density/y%porosity))* &
                        (dpcell(i,j)%pplist(k)%coff(1)*z(1,m)+dpcell(i,j)%pplist(k)%coff(2)* &
                        z(2,m)))/((dpcell(i,j)%plist(k)%density)* &
                        ((dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)+ &
                        (x%density/y%porosity)))


                        dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx-t1*dt

                        ! t1=(dpcell(i,j)%pplist(k)%porosity**2)*(x%mass/x%density)* &
                        !     (x%pressure-dpcell(i,j)%plist(k)%pressure)* &
                        !     (dpcell(i,j)%pplist(k)%coff(3)*z(1,m)+dpcell(i,j)%pplist(k)%coff(4)* &
                        ! z(2,m))/dpcell(i,j)%plist(k)%density


                        ! t1=((x%mass/x%density)*(dpcell(i,j)%pplist(k)%porosity**2)*&
                        ! (x%pressure+dpcell(i,j)%plist(k)%pressure)* &
                        ! (z(2,m)))/dpcell(i,j)%plist(k)%density

                        t1=2*((x%mass/x%density)*(dpcell(i,j)%pplist(k)%porosity**2)*&
                        ((x%pressure*dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)+&
                        (dpcell(i,j)%plist(k)%pressure*x%density/y%porosity))* &
                        (dpcell(i,j)%pplist(k)%coff(3)*z(1,m)+dpcell(i,j)%pplist(k)%coff(4)* &
                        z(2,m)))/((dpcell(i,j)%plist(k)%density)* &
                        ((dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)+ &
                        (x%density/y%porosity)))


                        dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy-t1*dt

                        end associate
                    end do
                    end if

                    dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vxs + &
                    dpcell(i,j)%plist(k)%vx

                    dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vys + &
                    dpcell(i,j)%plist(k)%vy

                    dpcell(i,j)%pplist(k)%vel=sqrt(dpcell(i,j)%plist(k)%vx**2+ &
                    dpcell(i,j)%plist(k)%vy**2)/dpcell(i,j)%pplist(k)%porosity

                    end if

                    if ((dpcell(i,j)%plist(k)%buffer)) then
                        dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vxs
                        dpcell(i,j)%plist(k)%vy=0.0_dp
                    end if

                end do

                end if
            end do
        end do
        !$omp end do        
        
    end subroutine comp_vel

    subroutine comp_pos
        
        implicit none   
        
        integer :: i,j,k

        ! New position calculations for fluid particles
        !$omp do schedule (runtime) private(i,k,j) collapse(2)   
        do j=sx,ex
            do i=sy,ey
                do k=1,dpcell(i,j)%ptot

                    if (dpcell(i,j)%plist(k)%tid==3) then

                    dpcell(i,j)%plist(k)%x=(dpcell(i,j)%plist(k)%x+ &
                    dpcell(i,j)%plist(k)%xs+dt*dpcell(i,j)%plist(k)%vx/dpcell(i,j)%pplist(k)%porosity)/2.0_dp

                    dpcell(i,j)%plist(k)%y=(dpcell(i,j)%plist(k)%y + &
                    dpcell(i,j)%plist(k)%ys+dt*dpcell(i,j)%plist(k)%vy/dpcell(i,j)%pplist(k)%porosity)/2.0_dp

                    if (dpcell(i,j)%plist(k)%buffer) then

                            if((dpcell(i,j)%plist(k)%x>(8*prrealx+brrealx+0.001_dp)) &
                            .and.(dpcell(i,j)%plist(k)%x<xrcutoff)) then

                                dpcell(i,j)%plist(k)%buffer=.false.

                            end if

                    end if

                    end if

                end do
            end do
        end do
        !$omp end do        
        
    end subroutine comp_pos

    subroutine timestep

        use,intrinsic :: ieee_arithmetic

        implicit none

        integer :: i,j
        real(dp),parameter :: Dm=1e-9
        
        ! Finding max velocity
        !$omp single
        umax=0.0_dp
        numax=0.0_dp
        !$omp end single

        !$omp do schedule(runtime) private(i,j) collapse(2)
            do j=sx,ex 
                do i=sy,ey
                    if (dpcell(i,j)%ptot/=0) then
                    dpcell(i,j)%maxvel=maxval(abs(dpcell(i,j)%pplist(1:dpcell(i,j)%ptot)%vel))
                    dpcell(i,j)%maxeddy=maxval(abs(dpcell(i,j)%pplist(1:dpcell(i,j)%ptot)%nut))
                    end if
                    ! umax=max(umax,dpcell(i,j)%maxvel)
                    ! numax=max(numax,dpcell(i,j)%maxeddy)
                end do
            end do
            
        !$omp end do

        !$omp do schedule(runtime) private(i,j) reduction(max:umax,numax) collapse(2)
            do j=sx,ex 
                do i=sy,ey

                    umax=max(umax,dpcell(i,j)%maxvel)
                    numax=max(numax,dpcell(i,j)%maxeddy)

                end do
            end do
            
        !$omp end do

        !$omp single
        if (t<=1.0_dp) then
        dt=min(sig1*(dl)/umax,sig1*((dl)**2)/((mu/rho)+numax),0.0010_dp)
        elseif ((t>1.0_dp).and.(t<=2.50_dp)) then
        dt=real(min(real(sig1*(dl)/umax),real(sig1*((dl)**2)*rho/mu),0.002))
        ! elseif ((t>2.0_dp).and.(t<=3.0_dp)) then
        ! dt=real(min(real(sig1*(dl)/umax),real(sig1*((dl)**2)*rho/mu),0.002))
        elseif ((t>3.0_dp).and.(t<=4.0_dp)) then
        dt=real(min(real(sig1*(dl)/umax),real(sig1*((dl)**2)*rho/mu),0.004))
        else
        dt=min(sig1*(dl)/umax,sig1*((dl)**2)/((mu/rho)+numax),0.0050_dp)
        end if
        if((.not.(ieee_is_finite(numax))).or.(ieee_is_nan(numax))) then
        numax=0.0_dp
        end if
        dtsol=sig1*(dl**2)/(Dm+numax/tschmidt)
        ! dtsol=dt/20.0_dp
        solsteps=1
        if (dtsol<dt) then
        solsteps=ceiling(dt/dtsol)
        dtsol=dt/solsteps        
        end if

        dtsol=merge(dtsol,dt,dtsol<dt)

        !$omp end single 
        
    end subroutine timestep
    

end module integrator
