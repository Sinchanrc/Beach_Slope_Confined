module output

    use initialize
    
    implicit none

    public
    
    contains

    subroutine print_fixbd()

        implicit none

        integer :: i,j,k,m

        !Boundary particle postion
        write(result,'("boundary_",i0,".txt")')(modifier+iter)
        open(11,file='boundary/'//result,status='replace')
        do j=(sx),(ex)
        do i=(sy),(ey)
            if (dpcell(i,j)%btot/=0) then
            do cout=1,dpcell(i,j)%btot
            if (dpcell(i,j)%plist(cout)%tid==1) then
            write(11,'(F10.3,1X,F10.3,1X,F10.3)')dpcell(i,j)%plist(cout)%x,&
            dpcell(i,j)%plist(cout)%y,dpcell(i,j)%plist(cout)%pressure
            end if
            end do
            end if
        end do
        end do
        close(11)

        return
    end subroutine print_fixbd

    subroutine print_ghostbd()

        implicit none

        integer :: i,j,k,m

        !Ghost Boundary particle postion
        write(result,'("ghost_",i0,".txt")')(modifier+iter)
        open(12,file='ghost/'//result,status='replace')
        do j=(sx),(ex)
        do i=(sy),(ey)
            if (dpcell(i,j)%btot/=0) then
            do cout=1,dpcell(i,j)%btot
            if (dpcell(i,j)%plist(cout)%tid==2) then
            write(12,'(F10.3,1X,F10.3,1X,F10.3)')dpcell(i,j)%plist(cout)%x,&
            dpcell(i,j)%plist(cout)%y,dpcell(i,j)%plist(cout)%pressure
            end if
            end do
            end if
        end do
        end do
        close(12)

        return
    end subroutine print_ghostbd

    ! subroutine print_interpolation()

    !     implicit none

    !     integer :: i,j,k,m

    !     !Ghost Boundary particle postion
    !     ! write(result,'("ghost_",i0,".txt")')(modifier+iter)
    !     open(12,file='interpolation.txt',status='replace')
    !     do j=(sx),(ex)
    !     do i=(sy),(ey)
    !         if (dpcell(i,j)%btot/=0) then
    !         do cout=1,dpcell(i,j)%btot
    !         if (dpcell(i,j)%plist(cout)%tid<=2) then
    !         write(12,'(F10.3,1X,F10.3)')dpcell(i,j)%plist(cout)%x+dpcell(i,j)%plist(cout)%posshift(1),&
    !         dpcell(i,j)%plist(cout)%y+dpcell(i,j)%plist(cout)%posshift(2)
    !         end if
    !         end do
    !         end if
    !     end do
    !     end do
    !     close(12)

    !     return
    ! end subroutine print_interpolation

    subroutine print_fluid()

        implicit none

        integer :: i,j,k,m

        write(result,'("fluid_",i0,".txt")')(modifier+iter)
        open(13,file='fluid/'//result,status='replace',action='write')
            do j=(sx),(ex)
            do i=(sy),(ey)
                if (dpcell(i,j)%ptot/=0) then
                do cout=1,dpcell(i,j)%ptot
                if ((dpcell(i,j)%plist(cout)%tid==3)) then
                write(13,'(F10.3,1X,F10.3,1X,F10.3,1X,ES10.3,1X,ES10.3,1X,F10.3)')dpcell(i,j)%plist(cout)%x,&
                dpcell(i,j)%plist(cout)%y,dpcell(i,j)%plist(cout)%pressure,dpcell(i,j)%plist(cout)%vx &
                ,dpcell(i,j)%plist(cout)%vy,dpcell(i,j)%pplist(cout)%porosity
                end if
                end do
                end if
            end do
            end do
        close(13)
        
    end subroutine print_fluid

    subroutine print_porous()

        implicit none

        integer :: i,j,k,m

        ! write(result,'("porous_",i0,".txt")')(modifier+iter)
        open(14,file='porous.txt',status='replace')
            do j=(sx),(ex)
            do i=(sy),(ey)
                if (dpcell(i,j)%porct/=0) then
                do cout=1,dpcell(i,j)%porct
                ! if (dpcell(i,j)%plist(cout)%tid==4) then
                write(14,'(F10.3,1X,F10.3)')dpcell(i,j)%porlist(cout)%x,&
                dpcell(i,j)%porlist(cout)%y
                ! end if
                end do
                end if
            end do
            end do
        close(14)
        
    end subroutine print_porous

    subroutine print_free()

        implicit none

        integer :: i,j,k,m
    
            ! write(result,'("free_",i0,".txt")')(modifier+iter)
            ! open(11,file='freesurf/'//result,status='replace')
            open(11,file="probe.txt",status='replace')
                ! do j=(sx),(ex)
                ! do i=(sy),(ey)
                    ! if (dpcell(i,j)%ptot/=0) then
                    do cout=1,160
                    ! if ((dpcell(i,j)%plist(cout)%tid>2).and.(.not.(dpcell(i,j)%plist(cout)%mirror))) then
                    ! write(11,*)fmatrix(dpcell(i,j)%plist(cout)%matid)%val(:),fvec(dpcell(i,j)%plist(cout)%matid)
                    write(11,'(F10.3,1X,F10.3)')probedata(1,cout),probedata(2,cout)
                    ! end if
                    end do
                    ! end if
                ! end do
                ! end do
            close(11)
            
    end subroutine print_free

    subroutine print_probe()

        implicit none

        integer :: i,j,k,m

        ! write(result,'("free_",i0,".txt")')(modifier+iter)
        ! open(11,file='freesurf/'//result,status='replace')
        open(11,file="probeloc.txt")
            ! do j=(sx),(ex)
            ! do i=(sy),(ey)
                ! if (dpcell(i,j)%ptot/=0) then
                ! do cout=1,dpcell(i,j)%ptot
                ! if ((dpcell(i,j)%plist(cout)%tid>2).and.(.not.(dpcell(i,j)%plist(cout)%mirror))) then
                ! write(11,*)fmatrix(dpcell(i,j)%plist(cout)%matid)%val(:),fvec(dpcell(i,j)%plist(cout)%matid)
                write(11,'(F10.3,1X,F10.3)')probe(1)%pos(1),probe(1)%pos(2)
                ! end if
                ! end do
                ! end if
            ! end do
            ! end do
        close(11)
        
    end subroutine print_probe

    subroutine combined()
        implicit none

        integer :: i,j,k,m
        character(len=70):: fluid,buffer3

        write(fluid,'("fluid_",i0,".txt")')(modifier+iter)
        write(buffer3,'("buffer_",i0,".txt")')(modifier+iter)
        open(11,file='fluid/'//fluid,status='replace',action='write')
        open(14,file='buffer/'//buffer3,status='replace',action='write')

        do j=sx,ex
            do i=sy,ey
                do cout=1,dpcell(i,j)%ptot
                if ((dpcell(i,j)%plist(cout)%tid==3).and.(.not.(dpcell(i,j)%plist(cout)%buffer))) then
                write(11,'(F10.3,1X,F10.3,1X,F10.3,1X,ES10.3,1X,ES10.3,1X,F10.5)')dpcell(i,j)%plist(cout)%x,&
                dpcell(i,j)%plist(cout)%y,dpcell(i,j)%plist(cout)%pressure,dpcell(i,j)%plist(cout)%vx &
                ,dpcell(i,j)%plist(cout)%vy,(dpcell(i,j)%plist(cout)%con+0.50_dp)

                elseif ((dpcell(i,j)%plist(cout)%tid==3).and.(dpcell(i,j)%plist(cout)%buffer)) then
                write(14,'(F10.3,1X,F10.3,1X,F10.3,1X,ES10.3,1X,ES10.3,1X,F10.5)')dpcell(i,j)%plist(cout)%x,&
                dpcell(i,j)%plist(cout)%y,dpcell(i,j)%plist(cout)%pressure,dpcell(i,j)%plist(cout)%vx &
                ,dpcell(i,j)%plist(cout)%vy,(dpcell(i,j)%plist(cout)%con+0.50_dp)

                end if
                end do

            end do
        end do

        close(11)
        close(14)
    
        
    end subroutine combined
    
end module output