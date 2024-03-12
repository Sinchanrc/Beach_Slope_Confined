module turbulence

    use particle
    use initialize
    use domain
    use kernel

    implicit none

    contains

    subroutine eddyvis

        implicit none
        real(dp) :: t1,t2,t3,t4
        real(dp),parameter ::cv=0.080_dp,ce=1.0_dp,cs=0.150_dp
        integer :: i,j,k,m
        
        !Mean strain and SGS tensor calculation for fluid particles
        !$omp do schedule (runtime) private(m,i,k,j,t1,t2,t3,t4) collapse(2)
        do j=sx,ex
            do i=sy,ey
                
                do k=(dpcell(i,j)%btot+1),dpcell(i,j)%ptot

                    t1=0.0_dp 
                    t2=0.0_dp 
                    t3=0.0_dp
                
                ! if ((dpcell(i,j)%plist(k)%tid==3).and.(dpcell(i,j)%pplist(k)%gradvx>=0.70_dp)) then
                    ! dpcell(i,j)%pplist(k)%meanstrn=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count ! Calculating mean strain tensor
                        associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                        y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                        z=>dpcell(i,j)%list(k)%klt(:,m))

                        t1=t1+ &
                            (x%mass*(x%vx/y%porosity- &
                            dpcell(i,j)%plist(k)%vx/dpcell(i,j)%pplist(k)%porosity)* &
                            (dpcell(i,j)%pplist(k)%coff(1)*z(1) &
                            +dpcell(i,j)%pplist(k)%coff(2)* &
                            z(2))) &
                            /x%density

                        t2=t2+ &
                            (x%mass*(x%vy/y%porosity- &
                            dpcell(i,j)%plist(k)%vy/dpcell(i,j)%pplist(k)%porosity)* &
                            (dpcell(i,j)%pplist(k)%coff(3)*z(1) &
                            +dpcell(i,j)%pplist(k)%coff(4)* &
                            z(2))) &
                            /x%density

                        t3=t3+ &
                        (0.50_dp*x%mass*((x%vx/y%porosity- &
                        dpcell(i,j)%plist(k)%vx/dpcell(i,j)%pplist(k)%porosity)* &
                        (dpcell(i,j)%pplist(k)%coff(3)*z(1)+dpcell(i,j)%pplist(k)%coff(4)* &
                        z(2))+ &
                        (x%vy/y%porosity- &
                        dpcell(i,j)%plist(k)%vy/dpcell(i,j)%pplist(k)%porosity)* &
                        (dpcell(i,j)%pplist(k)%coff(1)*z(1)+dpcell(i,j)%pplist(k)%coff(2)* &
                        z(2)))) &
                        /x%density

                        end associate
                    end do

                    ! dpcell(i,j)%pplist(k)%meanstrn(1)=t1
                    ! dpcell(i,j)%pplist(k)%meanstrn(4)=t2
                    ! dpcell(i,j)%pplist(k)%meanstrn(2)=t3
                    
                    ! dpcell(i,j)%pplist(k)%meanstrn(3)=t3


                    ! Calculating local strain

                    ! dpcell(i,j)%pplist(k)%S=sqrt(2*(dpcell(i,j)%pplist(k)%meanstrn(1)**2+ &
                    ! dpcell(i,j)%pplist(k)%meanstrn(4)**2))

                    t4=sqrt(2*(t1**2+ &
                    t2**2)+(2*t3)**2)

                    ! Calculating Eddy viscosity
                    dpcell(i,j)%pplist(k)%nut=((cs*dl)**2)*t4

                    ! Calculating tke
                    dpcell(i,j)%pplist(k)%tke=(cv/ce)*(dl*t4)**2

                    ! Calculating SGS stress tensor
                    ! dpcell(i,j)%pplist(k)%tau(1)=2*dpcell(i,j)%pplist(k)%meanstrn(1)*dpcell(i,j)%pplist(k)%nut &
                    ! -2*dpcell(i,j)%pplist(k)%tke/3.0_dp
                    
                    ! dpcell(i,j)%pplist(k)%tau(2)=2*dpcell(i,j)%pplist(k)%meanstrn(2)*dpcell(i,j)%pplist(k)%nut

                    ! dpcell(i,j)%pplist(k)%tau(3)=2*dpcell(i,j)%pplist(k)%meanstrn(3)*dpcell(i,j)%pplist(k)%nut

                    ! dpcell(i,j)%pplist(k)%tau(4)=2*dpcell(i,j)%pplist(k)%meanstrn(4)*dpcell(i,j)%pplist(k)%nut &
                    ! -2*dpcell(i,j)%pplist(k)%tke/3.0_dp
                ! else
                !     dpcell(i,j)%pplist(k)%tke=0.0_dp
                !     dpcell(i,j)%pplist(k)%nut=0.0_dp
                ! end if
                end do

            end do
        end do
        !$omp end do
        
        
    end subroutine eddyvis
    

end module turbulence
