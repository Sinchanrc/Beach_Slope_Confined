module gradient

    use particle
    use domain
    use kernel
    use functions
    use initialize

    implicit none
    
    contains

    subroutine gradcorr(p1,cell1,n,cell2,cell3,cell4,smln) !Kernel gradient correction
        implicit none
        type(particles),intent(inout) :: p1
        type(cell),intent(inout) :: cell1
        type(cell),intent(in) :: cell2,cell3,cell4
        type(cell) :: cellarr(2,2)
        real(dp) ::res2,radial,res1,res3,res4
        integer :: counter1,c3,c4
        real(dp),intent(in) :: smln
        integer,intent(in) :: n
        real(dp),dimension(4) :: res

        res=0.0_dp
        res2=0.0_dp
        radial=0.0_dp
        res1=0.0_dp
        res3=0.0_dp
        res4=0.0_dp
        cellarr(1,1)=cell1
        cellarr(2,1)=cell2
        cellarr(1,2)=cell3
        cellarr(2,2)=cell4

            if (cell1%list(n)%count/=0) then
            do counter1=1,cell1%list(n)%count
                associate(x1=>cell1%list(n)%interlist(1,counter1), &
                y1=>cell1%list(n)%interlist(2,counter1),p2=>cell1%list(n)%interlist(3,counter1))
                if (dpcell(y1,x1)%plist(p2)%tid/=4) then
                do c3=1,2
                    do c4=1,2
                    if ((cellarr(c4,c3)%cellid(1)==x1).and.(cellarr(c4,c3)%cellid(2)==y1)) then
                        radial=cell1%list(n)%dist(counter1)

                        res2=Wabx(cellarr(c4,c3)%plist(p2),p1,radial,smln)
                        res4=Waby(cellarr(c4,c3)%plist(p2),p1,radial,smln)

                        res1=(cellarr(c4,c3)%plist(p2)%mass)*(cellarr(c4,c3)%plist(p2)%x &
                        -p1%x)/(cellarr(c4,c3)%plist(p2)%density)
                        res3=(cellarr(c4,c3)%plist(p2)%mass)*(cellarr(c4,c3)%plist(p2)%y &
                        -p1%y)/(cellarr(c4,c3)%plist(p2)%density)


                        res(1)=res(1)+res2*res1
                        res(2)=res(2)+res4*res1
                        res(3)=res(3)+res2*res3
                        res(4)=res(4)+res4*res3
                    end if
                    end do
                end do
                end if
                end associate         
            end do
            end if

            call invertmat2D(res)

            cell1%pplist(n)%coff=res

        return
    end subroutine gradcorr

    subroutine compcorr(pm1,pm2)

        implicit none

        integer,intent(in) :: pm1,pm2
        integer :: i,j,k,m

        !$omp do schedule(runtime) private(i,k,j) collapse(2)
        do j=sx,ex
            do i=sy,ey
            if(dpcell(i,j)%ptot/=0) then

                do k=1,dpcell(i,j)%ptot

                    ! call combinedcorr(dpcell(i,j)%plist(k),dpcell(i,j),k,h1)

                ! if (dpcell(i,j)%plist(k)%tid==pm1) then

                if ((dpcell(i,j)%list(k)%count/=0)) then ! &.and.(dpcell(i,j)%pplist(k)%gradvx>=0.80_dp)

                    if(dpcell(i,j)%plist(k)%x<=(dpcell(i,j)%xleft+ &
                        abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                        dpcell(i,j)%plist(k)%y<=(dpcell(i,j)%ybot+ &
                        abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
        
                        if (pm2==1) then
                        call gradcorr(dpcell(i,j)%plist(k),dpcell(i,j),k,dpcell(i+1,j),&
                        dpcell(i,j-1),dpcell(i+1,j-1),h1)
                        ! else
                        ! call sgradcorrmg(dpcell(i,j)%plist(k),dpcell(i,j),k,dpcell(i+1,j),&
                        ! dpcell(i,j-1),dpcell(i+1,j-1),h1)
                        end if
        
                    else if(dpcell(i,j)%plist(k)%x>(dpcell(i,j)%xleft+ &
                        abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                        dpcell(i,j)%plist(k)%y<=(dpcell(i,j)%ybot+ &
                        abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
        
                        if (pm2==1) then
                        call gradcorr(dpcell(i,j)%plist(k),dpcell(i,j),k,dpcell(i+1,j),&
                        dpcell(i,j+1),dpcell(i+1,j+1),h1)
                        ! else
                        ! call sgradcorrmg(dpcell(i,j)%plist(k),dpcell(i,j),k,dpcell(i+1,j),&
                        ! dpcell(i,j+1),dpcell(i+1,j+1),h1)
                        end if
        
                    else if(dpcell(i,j)%plist(k)%x>(dpcell(i,j)%xleft+ &
                        abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                        dpcell(i,j)%plist(k)%y>(dpcell(i,j)%ybot+ &
                        abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then

                        if (pm2==1) then
                        call gradcorr(dpcell(i,j)%plist(k),dpcell(i,j),k,dpcell(i-1,j),&
                        dpcell(i-1,j+1),dpcell(i,j+1),h1)
                        ! else 
                        ! call sgradcorrmg(dpcell(i,j)%plist(k),dpcell(i,j),k,dpcell(i-1,j),&
                        ! dpcell(i-1,j+1),dpcell(i,j+1),h1)
                        end if
        
                    else if(dpcell(i,j)%plist(k)%x<=(dpcell(i,j)%xleft+ &
                        abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                        dpcell(i,j)%plist(k)%y>(dpcell(i,j)%ybot+ &
                        abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
        
                        if (pm2==1) then
                        call gradcorr(dpcell(i,j)%plist(k),dpcell(i,j),k,dpcell(i-1,j-1),&
                        dpcell(i,j-1),dpcell(i-1,j),h1)
                        ! else 
                        ! call sgradcorrmg(dpcell(i,j)%plist(k),dpcell(i,j),k,dpcell(i-1,j-1),&
                        ! dpcell(i,j-1),dpcell(i-1,j),h1)
                        end if
        
                    else if(dpcell(i,j)%plist(k)%x==(dpcell(i,j)%xleft+ &
                        abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                        dpcell(i,j)%plist(k)%y==(dpcell(i,j)%ybot+ &
                        abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
        
                        if (pm2==1) then
                        call gradcorr(dpcell(i,j)%plist(k),dpcell(i,j),k,dpcell(i-1,j-1),&
                        dpcell(i,j-1),dpcell(i-1,j),h1)
                        ! else 
                        ! call sgradcorrmg(dpcell(i,j)%plist(k),dpcell(i,j),k,dpcell(i-1,j-1),&
                        ! dpcell(i,j-1),dpcell(i-1,j),h1)

                    end if

                else

                    dpcell(i,j)%pplist(k)%coff=(/1.0_dp,0.0_dp,0.0_dp,1.0_dp/)
                    

                end if
    
                end if

                ! end if

                end do

            end if             
            end do
        end do
        !$omp end do


        

    
    end subroutine compcorr

end module gradient