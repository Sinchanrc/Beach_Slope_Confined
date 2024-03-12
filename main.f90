program dam_break

    use initialize
    use setup
    use boundary
    use solver
    use output
    use internal
    use memory
    use search
    use isph
    use integrator
    use turbulence
    use part_shift
    use gradient
    use scalar
    use porous

    implicit none

    integer :: j1,i1

    call init
    call gen_boundary
    call gen_fluid    
    call memalloc

    call combined
    call print_porous
    call print_fixbd
    call print_ghostbd
    call matrixid
    iter=iter+1

    !$omp parallel default(shared)
    call effpor
    !$omp end parallel

    do while(iter<51)

        told=t
        t=t+dt

        !$omp parallel default(shared)

        call projection 
        call cellshift
        call resetid
        call neighbour
        call effpor
        call eddyvis
        call int_vel
        !$omp end parallel
        call ppesolve
    
        !$omp parallel default(shared)
        call comp_vel        
        call comp_pos
        call cellshiftalt
        call neighbour
        call freesurf
        call effpor
        call opt2_shift
        call timestep
        !$omp end parallel
        iter=iter+1
    end do

    iter=1

    do j1=sx,ex 
        do i1=sy,ey
        if (dpcell(i1,j1)%ptot/=0) then
            do cout=1,dpcell(i1,j1)%ptot

                if ((dpcell(i1,j1)%plist(cout)%tid==3).and. &
                (.not.(dpcell(i1,j1)%plist(cout)%buffer))) then

                    dpcell(i1,j1)%plist(cout)%vx=0.0_dp!entry_vel
                    dpcell(i1,j1)%plist(cout)%vy=0.0_dp

                    ! if(((dpcell(i1,j1)%plist(cout)%y-yl-prrealy)-line_grad* &
                    ! (dpcell(i1,j1)%plist(cout)%x-xl))>0.0) then
                    !     dpcell(i1,j1)%plist(cout)%vx=entry_vel
                    ! end if

                end if


            end do
        end if
        end do
    end do

    iter1=1
    iter2=1

    do while(iter<8001)

        told=t
        t=t+dt

        if(((told)<iter1*ins_1).and. &
        (t)>=iter1*ins_1) then

        call remove_buffer1
        call insert_buffer1

        iter1=iter1+1

        end if

        if(((told)<iter2*ins_2).and. &
            (t)>=iter2*ins_2) then

            call remove_buffer2
            call insert_buffer2(1.0_dp)

            iter2=iter2+1

        end if

        !$omp parallel default(shared)

        call projection 
        call cellshift
        call resetid
        call neighbour
        call effpor
        call eddyvis
        call int_vel
        !$omp end parallel
        call ppesolve
    
        !$omp parallel default(shared)
        call comp_vel        
        call comp_pos
        call cellshiftalt
        call neighbour
        call freesurf
        call effpor
        call opt2_shift
        call timestep
        !$omp end parallel
        ! call combined
        iter=iter+1
    end do

    do j1=sx,ex 
        do i1=sy,ey
        if (dpcell(i1,j1)%ptot/=0) then
            do cout=1,dpcell(i1,j1)%ptot

                ! if ((dpcell(i1,j1)%plist(cout)%tid==3).and. &
                ! (.not.(dpcell(i1,j1)%plist(cout)%buffer))) then

                !     ! if(((dpcell(i1,j1)%plist(cout)%y-yl)-line_grad* &
                !     ! (dpcell(i1,j1)%plist(cout)%x-xl))>0.0) then

                !     dpcell(i1,j1)%plist(cout)%vx=0.0_dp!entry_vel
                !     dpcell(i1,j1)%plist(cout)%vy=0.0_dp

                !     ! end if

                ! end if

                if ((dpcell(i1,j1)%plist(cout)%tid==3)) then

                    dpcell(i1,j1)%plist(cout)%con=-0.50_dp

                    if(((dpcell(i1,j1)%plist(cout)%y-yl)-line_grad* &
                    (dpcell(i1,j1)%plist(cout)%x-xl))>0.0) then

                        ! if ((dpcell(i1,j1)%pplist(cout)%porosity>0.85_dp)) then

                        dpcell(i1,j1)%plist(cout)%mass=dpcell(i1,j1)%plist(cout)%mass*rel_den
                        dpcell(i1,j1)%plist(cout)%oden=dpcell(i1,j1)%plist(cout)%oden*rel_den
                        dpcell(i1,j1)%plist(cout)%density=dpcell(i1,j1)%plist(cout)%density*rel_den
                        dpcell(i1,j1)%plist(cout)%con=0.50_dp

                        ! end if

                    end if


                end if

            end do
        end if
        end do
    end do

    time_shift=t 
    iter=1
    iter1=1
    iter2=1
    t=0.0_dp

    do while(t<time)

        told=t
        t=t+dt
        
        ! if (dtsol>=dt) then
        ! !$omp parallel default(shared)
        ! call scalart
        ! call scalarupdate(dt)
        ! !$omp end parallel
        ! else
        do j1=1,solsteps
        !$omp parallel default(shared)
        call scalart
        call scalarupdate(dtsol)
        !$omp end parallel
        end do
        ! end if

        if(((told+time_shift)<iter1*ins_1).and. &
            (t+time_shift)>=iter1*ins_1) then

            call remove_buffer1
            call insert_buffer1

            iter1=iter1+1

        end if

        if(((told+time_shift)<iter2*ins_2).and. &
            (t+time_shift)>=iter2*ins_2) then

            call remove_buffer2
            call insert_buffer2(rel_den)

            iter2=iter2+1

        end if

        !$omp parallel default(shared)

        call projection 
        call cellshift
        call resetid
        call neighbour
        call effpor
        call eddyvis
        call int_vel
        !$omp end parallel
        call ppesolve
    
        !$omp parallel default(shared)
        call comp_vel        
        call comp_pos
        call cellshiftalt
        call neighbour
        call freesurf
        call effpor
        call opt2_shift
        call densityupdate
        call timestep
        call eddyvis 
        !$omp end parallel

        if ((told<iter*displaytime).and. &
        (t>=iter*displaytime)) then

        call combined
        iter=iter+1
        end if



    end do

    call combined
    
end program dam_break

