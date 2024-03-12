module pprobe

    use interactions
    use particle
    use initialize
    use kernel

    implicit none

    contains

    subroutine probesetup(loc)

        real(dp),intent(in),dimension(:,:) :: loc

        do i=1,pbno
            probe(i)%pos(1)=loc(1,i)
            probe(i)%pos(2)=loc(2,i)

        end do

        do j=2,cellx-1
        do i=2,celly-1
        do m=1,pbno
            if ((probe(m)%pos(1)>=dpcell(i,j)%xleft) .and. &
                (probe(m)%pos(1)<dpcell(i,j)%xright).and. &
                (probe(m)%pos(2)>=dpcell(i,j)%ybot).and. &
                (probe(m)%pos(2)<dpcell(i,j)%ytop)) then
                    probe(m)%ci=i
                    probe(m)%cj=j

            end if
        end do
        end do
        end do
        
        
    end subroutine probesetup

    subroutine psearch(list,n,cell1,cell2,cell3,cell4 &
        ,cell5,cell6,cell7,cell8,cell9,smln) 
        implicit none
        type(verlet),intent(inout) :: list
        type(cell),intent(in) :: cell1,cell2,cell3,cell4,cell5,cell6,cell7,cell8,cell9
        type(cell) :: blockset(9)
        real(dp),intent(in) :: smln
        integer,intent(in) :: n
        integer::counter1,counter2
        real(dp) :: radial
    

        blockset(1)=cell1
        blockset(2)=cell2
        blockset(3)=cell3
        blockset(4)=cell4
        blockset(5)=cell5
        blockset(6)=cell6
        blockset(7)=cell7
        blockset(8)=cell8
        blockset(9)=cell9


        list%count=0

            do counter1=1,9
            if (blockset(counter1)%ptot/=0) then
                do counter2=1,blockset(counter1)%ptot
                radial=sqrt((probe(n)%pos(1)-blockset(counter1)%plist(counter2)%x)**2+ &
                        (probe(n)%pos(2)-blockset(counter1)%plist(counter2)%y)**2 )
                if((radial<(2*blen*smln)).and.((blockset(counter1)%plist(counter2)%tid>2)))then
                    
                    list%count=list%count+1
                    list%interlist(1,list%count)=blockset(counter1)%cellid(1)
                    list%interlist(2,list%count)=blockset(counter1)%cellid(2)
                    list%interlist(3,list%count)=counter2
                    list%dist(list%count)=radial
                end if
                end do
            end if
            end do



    end subroutine psearch

    subroutine probevalue

        do k=1,pbno

            associate(i=>probe(k)%ci,j=>probe(k)%cj)
                call psearch(probe(k)%part,k,dpcell(i,j),dpcell(i,j-1), &
                dpcell(i,j+1),dpcell(i+1,j),dpcell(i-1,j),dpcell(i+1,j+1),dpcell(i+1,j-1)&
                ,dpcell(i-1,j-1),dpcell(i-1,j+1),h1)
    
            end associate
    
        end do 

        do k=1,pbno

        probe(k)%pressure=0.0_dp
        if (probe(k)%part%count/=0) then
            term1=0.0_dp
            ! probe(k)%pressure=0._dp
            do m=1,probe(k)%part%count
                associate(x=>probe(k)%part%interlist(1,m), &
                    y=>probe(k)%part%interlist(2,m), &
                    pp=>probe(k)%part%interlist(3,m))

                    term1=term1+dpcell(y,x)%plist(pp)%mass* &
                    Wab(probe(k)%part%dist(m),h1)/dpcell(y,x)%plist(pp)%density

                    probe(k)%pressure=probe(k)%pressure+dpcell(y,x)%plist(pp)%pressure* &
                                    dpcell(y,x)%plist(pp)%mass*Wab(probe(k)%part%dist(m),h1)&
                                    /dpcell(y,x)%plist(pp)%density          

                end associate
            end do

            probe(k)%pressure=probe(k)%pressure/term1



        end if
        end do

            probedata(1,iter)=t*sqrt(abs(g)/wc)
            probedata(2,iter)=probe(1)%pressure/(rho*abs(g)*wc)
        
    end subroutine probevalue

end module pprobe
