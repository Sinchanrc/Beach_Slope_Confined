module memory

    use initialize

    implicit none
    

    contains

    subroutine memalloc()

        implicit none

        integer :: i,j,k,m

        !Allocating arrays in the cells
        !$omp parallel do private(i,j,k) schedule (runtime) collapse(2) default(shared)
            do j=sx,ex
            do i=sy,ey

            if (dpcell(i,j)%btot/=0) then
                allocate(dpcell(i,j)%list(fplistmax),dpcell(i,j)%pplist(fplistmax)) 
                ! allocate(dpcell(i,j)%list2(dpcell(i,j)%btot))
                do k=1,fplistmax
                    allocate(dpcell(i,j)%list(k)%nh(fplist), &
                    dpcell(i,j)%list(k)%pnh(fplist),&
                    dpcell(i,j)%list(k)%dist(fplist),dpcell(i,j)%list(k)%klt(2,fplist))
                    ! dpcell(i,j)%list(k)%interlist=0
                    dpcell(i,j)%list(k)%dist=0.0_dp

                    allocate(dpcell(i,j)%pplist(k)%coff(4))

                    do m=1,fplist 
                        dpcell(i,j)%list(k)%nh(m)%part=>dpcell(i,j)%plist(m)
                        dpcell(i,j)%list(k)%pnh(m)%ppart=>dpcell(i,j)%pplist(m)
                    end do

                end do

                ! do k=1,dpcell(i,j)%btot 

                !     allocate(dpcell(i,j)%list2(k)%interlist(3,fplist), &
                !     dpcell(i,j)%list2(k)%dist(fplist))
                !     dpcell(i,j)%list2(k)%interlist=0
                !     dpcell(i,j)%list2(k)%dist=0.0_dp

                ! end do


            else
                allocate(dpcell(i,j)%list(fplistmax),dpcell(i,j)%pplist(fplistmax))  
                do k=1,fplistmax
                    allocate(dpcell(i,j)%list(k)%nh(fplist),dpcell(i,j)%list(k)%pnh(fplist),&
                    dpcell(i,j)%list(k)%dist(fplist),dpcell(i,j)%list(k)%klt(2,fplist))
                    ! dpcell(i,j)%list(k)%interlist=0
                    dpcell(i,j)%list(k)%dist=0.0_dp
                    allocate(dpcell(i,j)%pplist(k)%coff(4))

                    do m=1,fplist 
                        dpcell(i,j)%list(k)%nh(m)%part=>dpcell(i,j)%plist(m)
                        dpcell(i,j)%list(k)%pnh(m)%ppart=>dpcell(i,j)%pplist(m)
                    end do

                end do  

            end if

                allocate(dpcell(i,j)%tn%list(ceiling(0.25*fplist)),&
                dpcell(i,j)%ts%list(ceiling(0.25*fplist)),&
                dpcell(i,j)%te%list(ceiling(0.25*fplist)), &
                dpcell(i,j)%tw%list(ceiling(0.25*fplist)),dpcell(i,j)%tne%list(ceiling(0.25*fplist)),&
                dpcell(i,j)%tnw%list(ceiling(0.25*fplist)), &
                dpcell(i,j)%tse%list(ceiling(0.25*fplist)),dpcell(i,j)%tsw%list(ceiling(0.25*fplist)))
            end do
            end do
        !$omp end parallel do

    end subroutine memalloc

    subroutine matrixid()
        implicit none

        integer :: i,j,k,mo

        count=0

        do j=sx,ex 
            do i=sy,ey
            if (dpcell(i,j)%ptot/=0) then
                do cout=1,dpcell(i,j)%ptot

                    ! if ((dpcell(i,j)%plist(cout)%tid==3)) then
                    count=count+1
                    dpcell(i,j)%plist(cout)%pid=count
                    dpcell(i,j)%plist(cout)%matid=dpcell(i,j)%plist(cout)%pid

                end do
            end if
            end do
        end do        

        finmax=count!+reserve_par-50

        do i=1,reserve_par
            count=count+1
            reserve%tank(i)=count
        end do

        allocate(frow(finmax+1),fvec(finmax),fsol(finmax),fguess(finmax))
        ! allocate(fval(finmax*ceiling(fac2*fplistmax)),fcol(finmax*ceiling(fac2*fplistmax)))
        allocate(fmatrix(finmax),perm(finmax),pguess(finmax))

        fguess=1000.0_dp
        fsol=0.0_dp
        pguess=0.0_dp

        do j=sx,ex 
            do i=sy,ey
            if (dpcell(i,j)%ptot/=0) then
                do cout=1,dpcell(i,j)%ptot

                    ! if ((dpcell(i,j)%plist(cout)%tid/=4)) then

                    fguess(dpcell(i,j)%plist(cout)%matid)=dpcell(i,j)%plist(cout)%pressure

                    ! end if


                end do
            end if
            end do
        end do

        ! ! Allocating fmatrix

            do i=1,finmax

            allocate(fmatrix(i)%val(ceiling(fac2*fplistmax)))
            allocate(fmatrix(i)%col(ceiling(fac2*fplistmax)))

            end do

            allocate(dpar(128),tmp(finmax*(2*150+1)+(150*(150+9))/2+1))
            fsol=fguess
            
        
    end subroutine matrixid


    
end module memory