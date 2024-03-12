module mls

    use kernel
    use initialize
    use functions
    use domain
    use particle

    implicit none

    integer :: i,j,k,m
    
    contains

    !Moving least Squares Kernel value 
    subroutine Wmls(a,b,c,smln,err)
        implicit none 
        integer,intent(in) :: a,b,c
        real(dp) ,intent(in) :: smln
        real(dp) :: array(3,3)
        integer :: is
        logical,intent(out) :: err

        array=0.0_dp

        if (dpcell(b,a)%list(c)%count/=0) then
            do is=1,dpcell(b,a)%list(c)%count
                associate(x=>dpcell(b,a)%list(c)%interlist(1,is),&
                    y=>dpcell(b,a)%list(c)%interlist(2,is),&
                    pp=>dpcell(b,a)%list(c)%interlist(3,is) )

                array(1,1)=array(1,1)+Wab(dpcell(b,a)%list(c)%dist(is),blen*smln)*&
                (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                array(2,1)=array(2,1)+Wab(dpcell(b,a)%list(c)%dist(is),blen*smln)*&
                (dpcell(b,a)%plist(c)%x+dpcell(b,a)%plist(c)%posshift(1)-dpcell(y,x)%plist(pp)%x)*&
                (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                array(3,1)=array(3,1)+Wab(dpcell(b,a)%list(c)%dist(is),blen*smln)*&
                (dpcell(b,a)%plist(c)%y+dpcell(b,a)%plist(c)%posshift(2)-dpcell(y,x)%plist(pp)%y)*&
                (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                array(1,2)=array(1,2)+Wab(dpcell(b,a)%list(c)%dist(is),blen*smln)*&
                (dpcell(b,a)%plist(c)%x+dpcell(b,a)%plist(c)%posshift(1)-dpcell(y,x)%plist(pp)%x)*&
                (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                array(2,2)=array(2,2)+Wab(dpcell(b,a)%list(c)%dist(is),blen*smln)*&
                ((dpcell(b,a)%plist(c)%x+dpcell(b,a)%plist(c)%posshift(1)-dpcell(y,x)%plist(pp)%x)**2)&
                *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                array(3,2)=array(3,2)+Wab(dpcell(b,a)%list(c)%dist(is),blen*smln)*&
                (dpcell(b,a)%plist(c)%x+dpcell(b,a)%plist(c)%posshift(1)-dpcell(y,x)%plist(pp)%x)*&
                (dpcell(b,a)%plist(c)%y+dpcell(b,a)%plist(c)%posshift(2)-dpcell(y,x)%plist(pp)%y)*&
                (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                array(1,3)=array(1,3)+Wab(dpcell(b,a)%list(c)%dist(is),blen*smln)*&
                (dpcell(b,a)%plist(c)%y+dpcell(b,a)%plist(c)%posshift(2)-dpcell(y,x)%plist(pp)%y)*&
                (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                array(2,3)=array(2,3)+Wab(dpcell(b,a)%list(c)%dist(is),blen*smln)*&
                (dpcell(b,a)%plist(c)%x+dpcell(b,a)%plist(c)%posshift(1)-dpcell(y,x)%plist(pp)%x)*&
                (dpcell(b,a)%plist(c)%y+dpcell(b,a)%plist(c)%posshift(2)-dpcell(y,x)%plist(pp)%y)*&
                (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                array(3,3)=array(3,3)+Wab(dpcell(b,a)%list(c)%dist(is),blen*smln)*&
                ((dpcell(b,a)%plist(c)%y+dpcell(b,a)%plist(c)%posshift(2)-dpcell(y,x)%plist(pp)%y)**2)&
                *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)
                end associate
            end do

                call invertmat3D(array,err)


                if (.not.(err)) then
                dpcell(b,a)%pplist(c)%bmls(1)=array(1,1)
                dpcell(b,a)%pplist(c)%bmls(2)=array(2,1)
                dpcell(b,a)%pplist(c)%bmls(3)=array(3,1)
                end if
        end if


        
    end subroutine

    subroutine comp_ghost

        implicit none

        ! logical :: error
        real(dp) :: t1,t2

        !$omp do schedule(runtime) private(m,i,k,t1) collapse(2)
        do j=sx,ex
            do i=sy,ey
                if(dpcell(i,j)%ptot/=0) then
                    
                    do k=1,dpcell(i,j)%ptot

                    if (dpcell(i,j)%plist(k)%tid<=2) then

                        ! dpcell(i,j)%plist(k)%pressure=0.0d0
                        dpcell(i,j)%plist(k)%vx=0.0_dp
                        dpcell(i,j)%plist(k)%vy=0.0_dp
                        t1=0.0_dp
                        if (dpcell(i,j)%list(k)%count/=0) then
                            ! call Wmls(j,i,k,h1,error)

                            ! if (.not.(error)) then
                            ! associate(num1=>dpcell(i,j)%list2(k)%count)
                            ! do m=1,num1
                            !     associate(x=>dpcell(i,j)%list2(k)%interlist(1,m), &
                            !     y=>dpcell(i,j)%list2(k)%interlist(2,m), &
                            !     pp=>dpcell(i,j)%list2(k)%interlist(3,m))

                            !     t1=Wab(dpcell(i,j)%list2(k)%dist(m),blen*h1)*&
                            !     (dpcell(i,j)%pplist(k)%bmls(1)+dpcell(i,j)%pplist(k)%bmls(2)* &
                            !     (dpcell(i,j)%plist(k)%x&
                            !     +dpcell(i,j)%plist(k)%posshift(1)-dpcell(y,x)%plist(pp)%x)+ &
                            !     dpcell(i,j)%pplist(k)%bmls(3)*(dpcell(i,j)%plist(k)%y&
                            !     +dpcell(i,j)%plist(k)%posshift(2)-dpcell(y,x)%plist(pp)%y))*&
                            !     (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                            !     dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx- &
                            !     (dpcell(y,x)%plist(pp)%vx)*t1

                            !     dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy- &
                            !     (dpcell(y,x)%plist(pp)%vy)*t1

                            !     end associate
                            ! end do
                            ! end associate
                            ! ! dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vx
                            ! ! dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vy

                            ! else

                            do m=1,dpcell(i,j)%list2(k)%count                           
                            associate(x=>dpcell(i,j)%list2(k)%interlist(1,m), &
                                y=>dpcell(i,j)%list2(k)%interlist(2,m), &
                                pp=>dpcell(i,j)%list2(k)%interlist(3,m))

                                ! dpcell(i,j)%list(k)%error=.true.

                                dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx- &
                                Wab(dpcell(i,j)%list2(k)%dist(m),h1)*dpcell(y,x)%plist(pp)%vx* &
                                dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density

                                dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy- &
                                Wab(dpcell(i,j)%list2(k)%dist(m),h1)*dpcell(y,x)%plist(pp)%vy* &
                                dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density

                                t1=t1+dpcell(y,x)%plist(pp)%mass* &
                                Wab(dpcell(i,j)%list2(k)%dist(m),h1)/dpcell(y,x)%plist(pp)%density

                            end associate
                            end do

                
                            dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx/t1                            
                            dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy/t1
                            
                            ! if (dpcell(i,j)%plist(k)%tid==1) then
                                ! dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vx!0.0_dp
                                ! dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vy!0.0_dp
                            ! end if
                            ! end if

                        end if

                    end if

                    end do                    
                end if           
            end do
        end do
        !$omp end do

        ! !$omp parallel do schedule(runtime) private(t2,i,k) collapse(2) default(shared)
        ! do j=sx,ex
        !     do i=sy,ey
        !         if(dpcell(i,j)%btot/=0) then
                    
        !             do k=1,dpcell(i,j)%btot
        !                 if ((dpcell(i,j)%plist(k)%tid==2).and.(dpcell(i,j)%list(k)%count/=0)) then

        !                     if ((dpcell(i,j)%plist(k)%xnorm==1).and.(dpcell(i,j)%plist(k)%ynorm==0)) then
        !                         t2=0.0d0

        !                         elseif ((dpcell(i,j)%plist(k)%ynorm==1).and.(dpcell(i,j)%plist(k)%xnorm==0)) then
        !                         t2=-g*rho*dpcell(i,j)%plist(k)%posshift(2)

        !                         elseif ((dpcell(i,j)%plist(k)%ynorm==1).and.(dpcell(i,j)%plist(k)%xnorm==1)) then
        !                         t2=-g*rho*dpcell(i,j)%plist(k)%posshift(2)
        !                     end if

        !                     dpcell(i,j)%plist(k)%pressure=dpcell(i,j)%plist(k)%pressure+t2

        !                 end if
        !             end do

        !         end if           
        !     end do
        ! end do
        ! !$omp end parallel do
            
    end subroutine comp_ghost

    subroutine ghost_en

        implicit none

        real(dp) :: t1

        !$omp do schedule(runtime) private(m,i,k,t1) collapse(2)
        do j=sx,ex
            do i=sy,ey
                if(dpcell(i,j)%ptot/=0) then
                    
                    do k=1,dpcell(i,j)%ptot

                    if (dpcell(i,j)%plist(k)%tid<=2) then

                        ! dpcell(i,j)%plist(k)%pressure=0.0d0
                        dpcell(i,j)%pplist(k)%tke=0.0_dp
                        dpcell(i,j)%pplist(k)%nut=0.0_dp
                        t1=0.0_dp
                        if (dpcell(i,j)%list(k)%count/=0) then

                            do m=1,dpcell(i,j)%list2(k)%count                           
                            associate(x=>dpcell(i,j)%list2(k)%interlist(1,m), &
                                y=>dpcell(i,j)%list2(k)%interlist(2,m), &
                                pp=>dpcell(i,j)%list2(k)%interlist(3,m))

                                dpcell(i,j)%pplist(k)%tke=dpcell(i,j)%pplist(k)%tke- &
                                Wab(dpcell(i,j)%list2(k)%dist(m),h1)*dpcell(y,x)%pplist(pp)%tke* &
                                dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density

                                dpcell(i,j)%pplist(k)%nut=dpcell(i,j)%pplist(k)%nut- &
                                Wab(dpcell(i,j)%list2(k)%dist(m),h1)*dpcell(y,x)%pplist(pp)%nut* &
                                dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density


                                t1=t1+dpcell(y,x)%plist(pp)%mass* &
                                Wab(dpcell(i,j)%list2(k)%dist(m),h1)/dpcell(y,x)%plist(pp)%density

                            end associate
                            end do

                
                            dpcell(i,j)%pplist(k)%tke=dpcell(i,j)%pplist(k)%tke/t1 
                            dpcell(i,j)%pplist(k)%nut=dpcell(i,j)%pplist(k)%nut/t1                            
                            
                        end if

                    end if

                    end do                    
                end if           
            end do
        end do
        !$omp end do

    end subroutine ghost_en
    
end module mls