module internal

    use initialize
    use particle

    implicit none

    contains

    subroutine gen_fluid()

        implicit none

        integer :: i,j,k,m,bfy1
        integer :: step=0

        real(dp) :: bounlen,bounlen2

        bfy1=fpy

        bounlen=fpy*2*prrealy
        bounlen2=open_lhs*2*prrealy

        ins_1=(4*prrealx/sqrt(por))/abs(entry_vel/0.381_dp)
        ins_2=4*prrealx/((bounlen/bounlen2)*abs(entry_vel))

        rv_buf_r=((brrealx)*((2*bl)-1))+(fpx-1)*2*prrealx/sqrt(por)+prrealx/sqrt(por)+domain_shift&
                    +3*2*prrealx/sqrt(por)+0.5*prrealx/sqrt(por)

        rv_buf_l=2*prrealx+brrealx+0.5*prrealx

        !$omp parallel default(shared)
        !Setting up fluid particle positions and pressure
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j =1,fpx
            do i=fpy,1,-1
                flist(i,j)%y=((brrealy)*((2*bl)-1))+(fpy-i)*2*prrealy/sqrt(por)+prrealy/sqrt(por)

                flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx/sqrt(por)+prrealx/sqrt(por)+domain_shift

            flist(i,j)%vx=entry_vel!*3.0_dp*(((flist(i,j)%y-(brrealy)*(2*bl-1))/bounlen) &
                            ! -0.5_dp*((flist(i,j)%y-(brrealy)*(2*bl-1))/bounlen)**2)
            flist(i,j)%vy=0.0_dp

                flist(i,j)%pressure=(-flist(i,j)%y+((brrealy*distfac)*((2*bl)-1))+wc+prrealy/sqrt(por))*rho*abs(g)


            end do
            end do
        !$omp end do

        !Setting particle id for fluid particles
        !$omp single 
            count=0
            finmin=count+1
            do j =1,fpx
                do i=fpy,1,-1

                    count=count+1
                    flist(i,j)%pid=count
                    flist(i,j)%tid=3


                if(((flist(i,j)%y-yl)-line_grad*(flist(i,j)%x-xl))>0.0) then
                    flist(i,j)%tid=0
                end if
                
                end do
            end do
        
        !$omp end single

        !Distributing fluid particles to cells
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j=2,cellx-1
                do i=2,celly-1
                do l1=1,fpx
                    do k=1,fpy
                    if ((flist(k,l1)%x>=dpcell(i,j)%xleft) .and. &
                        (flist(k,l1)%x<dpcell(i,j)%xright).and. &
                        (flist(k,l1)%y>=dpcell(i,j)%ybot).and. &
                        (flist(k,l1)%y<dpcell(i,j)%ytop).and.(flist(k,l1)%tid/=0))then
                        dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)=flist(k,l1)

                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%mass=fmass
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%density=rho

                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%oden=dpcell(i,j)%plist&
                        (dpcell(i,j)%ptot)%density
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%ovol=fmass/rho

                    end if
                    end do
                end do
                end do
            end do
        !$omp end do

        !Buffer on RHS

        !$omp single

            deallocate(flist)
            allocate(flist(bfy1,5),buffer1(bfy1,2))
            k=1
            
        !$omp end single

            !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j =1,5
                do i=bfy1,1,-1

                flist(i,j)%y=((brrealy)*((2*bl)-1))+(bfy1-i)*2*prrealy/sqrt(por)+prrealy/sqrt(por)
                flist(i,j)%x=((brrealx)*((2*bl)-1))+(fpx-1)*2*prrealx/sqrt(por)+prrealx/sqrt(por)+domain_shift&
                    +j*2*prrealx/sqrt(por)

                flist(i,j)%buffer=.true.
                flist(i,j)%vx=entry_vel!*3.0_dp*(((flist(i,j)%y-(brrealy)*(2*bl-1))/bounlen) &
                            ! -0.5_dp*((flist(i,j)%y-(brrealy)*(2*bl-1))/bounlen)**2)
                flist(i,j)%vy=0.0_dp
                flist(i,j)%pressure=(-flist(i,j)%y+((brrealy*distfac)*((2*bl)-1))+wc+prrealy/sqrt(por))*rho*abs(g)

                if(j>3) then

                buffer1(i,(j-3))%y=((brrealy)*((2*bl)-1))+(bfy1-i)*2*prrealy/sqrt(por)+prrealy/sqrt(por)
                buffer1(i,(j-3))%x=((brrealx)*((2*bl)-1))+(fpx-1)*2*prrealx/sqrt(por)+prrealx/sqrt(por)+domain_shift&
                    +j*2*prrealx/sqrt(por)

                buffer1(i,(j-3))%buffer=.true.                

                buffer1(i,(j-3))%vx=entry_vel!*3.0_dp*(((flist(i,j)%y-(brrealy)*(2*bl-1))/bounlen) &
                            ! -0.5_dp*((flist(i,j)%y-(brrealy)*(2*bl-1))/bounlen)**2)
                buffer1(i,(j-3))%vy=0.0_dp
        
                buffer1(i,(j-3))%pressure=(-flist(i,j)%y+((brrealy*distfac)*((2*bl)-1))+wc+prrealy/sqrt(por))*rho*abs(g)
                k=k+1


                end if
                end do
            end do
            !$omp end do

        !$omp single 
            do j =1,5
                do i=bfy1,1,-1

                    count=count+1
                    flist(i,j)%pid=count
                    flist(i,j)%tid=3
                end do
            end do
        
        !$omp end single

        !Distributing RHS buffer to cells
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j=2,cellx-1
                do i=2,celly-1
                do l1=1,5
                    do k=1,bfy1
                    if ((flist(k,l1)%x>=dpcell(i,j)%xleft) .and. &
                        (flist(k,l1)%x<dpcell(i,j)%xright).and. &
                        (flist(k,l1)%y>=dpcell(i,j)%ybot).and. &
                        (flist(k,l1)%y<dpcell(i,j)%ytop))then
                        dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)=flist(k,l1)

                        dpcell(i,j)%entrybuff=.true.

                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%mass=fmass
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%density=rho

                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%oden=dpcell(i,j)%plist&
                        (dpcell(i,j)%ptot)%density
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%ovol=fmass/rho

                    end if
                    end do
                end do
                end do
            end do
        !$omp end do

        !Counting number of cells on RHS
        !$omp single

            entrycounter1=0

            do j=2,cellx-1
                do i=2,celly-1

                    if (dpcell(i,j)%entrybuff) then

                        entrycounter1=entrycounter1+1

                    end if
                end do
            end do

        !$omp end single 

        !$omp single

            allocate(entrycell1(entrycounter1))
            step=1

            ! Pointing to cells containing entry points on RHS
            do j=2,cellx-1
                do i=2,celly-1

                    if (dpcell(i,j)%entrybuff) then

                        entrycell1(step)%bcell=>dpcell(i,j)
                        step=step+1

                    end if

                end do
            end do
        
            deallocate(flist)
            fpy=floor(real(coastal_ht,dp)/(2*real(prrealy,dp)))+1
            fpx=floor(real(2.0_dp,dp)/(2*real(prrealx,dp)))!+1
            allocate(flist(fpy,fpx))
        !$omp end single

        !!!!!!!!!!!!!!!!!!! Setting up porous media(non-deform)!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Setting up fluid particle positions and pressure
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j =1,fpx
                do i=fpy,1,-1
    
    
                flist(i,j)%y=((brrealy)*((2*bl)-1))+(fpy-i)*2*prrealy+prrealy+set_ht
                ! if (mod(i,2)==0) then
                ! flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx+2*prrealx
                ! else
                flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx+prrealx+domain_shift
                ! end if
                flist(i,j)%vx=0.0_dp!-2.5_dp/(3600*24*por)
                flist(i,j)%vy=0.0_dp
                if (i==1) then
                    flist(i,j)%pressure=0.0_dp
    
                else
                    flist(i,j)%pressure=(-flist(i,j)%y+((brrealy*distfac)*((2*bl)-1))+coastal_ht &
                    +prrealy+set_ht)*rho*abs(g)
                
                end if

                end do
                end do
        !$omp end do

        !Setting particle id for fluid particles
        !$omp single 
            do j =1,fpx
                do i=fpy,1,-1

                        count=count+1
                        flist(i,j)%pid=count
                        flist(i,j)%tid=3


                        if(((flist(i,j)%y-yl+prrealy)-line_grad*(flist(i,j)%x-xl))<0.0) then
                            flist(i,j)%tid=0

                        else
                            flist(i,j)%y=flist(i,j)%y+2*dl

                        end if


                ! end if
                
                end do
            end do
        
        !$omp end single

        !Distributing fluid particles to cells
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j=2,cellx-1
                do i=2,celly-1
                do l1=1,fpx
                    do k=1,fpy

                    if ((flist(k,l1)%x>=dpcell(i,j)%xleft) .and. &
                        (flist(k,l1)%x<dpcell(i,j)%xright).and. &
                        (flist(k,l1)%y>=dpcell(i,j)%ybot).and. &
                        (flist(k,l1)%y<dpcell(i,j)%ytop).and.(flist(k,l1)%tid/=0))then
                        dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)=flist(k,l1)
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%mass=fmass*1.0
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%density=rho*1.0
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%oden=dpcell(i,j)%plist&
                        (dpcell(i,j)%ptot)%density

                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%ovol=fmass/rho

                    end if


                    end do
                end do
                end do
            end do
        !$omp end do

        !Buffer on LHS
        !$omp single

            deallocate(flist)
            allocate(flist(open_lhs,5),buffer2(open_lhs,2))
            k=1

        !$omp end single

        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j =1,5
                do i=open_lhs,1,-1
                    flist(i,j)%y=2*bny*brrealy-brrealy+prrealy+(open_lhs-i)*2*prrealy
                    flist(i,j)%x=(j-1)*2*prrealx+brrealx
                flist(i,j)%vx=0.0_dp!(-bounlen/bounlen2)*entry_vel!*3.0_dp*(((flist(i,j)%y-(brrealy)*(2*bny-1))/bounlen2) &
                                ! -0.5_dp*((flist(i,j)%y-(brrealy)*(2*bny-1))/bounlen2)**2)
                flist(i,j)%vy=0.0_dp

                    flist(i,j)%buffer=.true.
                    flist(i,j)%pressure=(-flist(i,j)%y+((brrealy*distfac)*((2*bl)-1))+wc+prrealy)*rho*abs(g)

                    if (j<=2) then

                    buffer2(i,j)%y=2*bny*brrealy-brrealy+prrealy+(open_lhs-i)*2*prrealy
                    buffer2(i,j)%x=(j-1)*2*prrealx+brrealx

                buffer2(i,j)%vx=0.0_dp!(-bounlen/bounlen2)*entry_vel!*3.0_dp*(((flist(i,j)%y-(brrealy)*(2*bny-1))/bounlen2) &
                                ! -0.5_dp*((flist(i,j)%y-(brrealy)*(2*bny-1))/bounlen2)**2)
                buffer2(i,j)%vy=0.0_dp

                    buffer2(i,j)%buffer=.true.
                    buffer2(i,j)%pressure=(-flist(i,j)%y+((brrealy*distfac)*((2*bl)-1))+wc+prrealy)*rho*abs(g)*rel_den
                    k=k+1

                    end if


                end do
                end do
            !$omp end do

        !$omp single 
            do j =1,5
                do i=1,open_lhs

                        count=count+1
                        flist(i,j)%pid=count
                        flist(i,j)%tid=3

                
                end do
            end do
        
        !$omp end single

        !Distributing LHS buffer
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j=2,cellx-1
                do i=2,celly-1
                do l1=1,5
                    do k=1,open_lhs

                    if ((flist(k,l1)%x>=dpcell(i,j)%xleft) .and. &
                        (flist(k,l1)%x<dpcell(i,j)%xright).and. &
                        (flist(k,l1)%y>=dpcell(i,j)%ybot).and. &
                        (flist(k,l1)%y<dpcell(i,j)%ytop).and.(flist(k,l1)%tid/=0))then
                        dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)=flist(k,l1)
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%mass=fmass*1.0
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%density=rho*1.0
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%oden=dpcell(i,j)%plist&
                        (dpcell(i,j)%ptot)%density
                        dpcell(i,j)%entrybuff=.true.

                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%ovol=fmass/rho

                    end if


                    end do
                end do
                end do
            end do
        !$omp end do

        !$omp single
        
            deallocate(flist)

        !$omp end single

        !Counting number of cells on LHS
        !$omp single

            entrycounter2=0

            do j=2,cellx/2
                do i=2,celly-1

                    if (dpcell(i,j)%entrybuff) then

                        entrycounter2=entrycounter2+1

                    end if
                end do
            end do

        !$omp end single 

        !$omp single

            allocate(entrycell2(entrycounter2))
            step=1

            ! Pointing to cells containing entry points on RHS
            do j=2,cellx/2
                do i=2,celly-1

                    if (dpcell(i,j)%entrybuff) then

                        entrycell2(step)%bcell=>dpcell(i,j)
                        step=step+1

                    end if

                end do
            end do

        !$omp end single

        !$omp end parallel
    
    end subroutine gen_fluid

    subroutine insert_buffer1()

        implicit none

        integer :: i,j,k

        !$omp parallel do schedule(runtime) private(j) default(shared)
            do j=1,entrycounter1
                do i=1,size(buffer1,2)
                    do k=1,size(buffer1,1)
                    if ((buffer1(k,i)%x>=entrycell1(j)%bcell%xleft) .and. &
                        (buffer1(k,i)%x<entrycell1(j)%bcell%xright).and. &
                        (buffer1(k,i)%y>=entrycell1(j)%bcell%ybot).and. &
                        (buffer1(k,i)%y<entrycell1(j)%bcell%ytop))then
                        entrycell1(j)%bcell%ptot=entrycell1(j)%bcell%ptot+1
                        entrycell1(j)%bcell%plist(entrycell1(j)%bcell%ptot)=buffer1(k,i)
                        entrycell1(j)%bcell%plist(entrycell1(j)%bcell%ptot)%tid=3

                        entrycell1(j)%bcell%plist(entrycell1(j)%bcell%ptot)%mass=fmass
                        entrycell1(j)%bcell%plist(entrycell1(j)%bcell%ptot)%density=rho

                        entrycell1(j)%bcell%plist(entrycell1(j)%bcell%ptot)%oden=entrycell1(j)%bcell &
                        %plist(entrycell1(j)%bcell%ptot)%density
                        entrycell1(j)%bcell%plist(entrycell1(j)%bcell%ptot)%ovol=fmass/rho
                        entrycell1(j)%bcell%plist(entrycell1(j)%bcell%ptot)%con=-0.50_dp

                    end if
                    end do
                end do
            end do
        !$omp end parallel do
        
    end subroutine insert_buffer1

    subroutine insert_buffer2(relden)

        implicit none

        integer :: i,j,k
        real(dp),intent(in) :: relden

        !$omp parallel do schedule(runtime) private(j) default(shared)
            do j=1,entrycounter2
                do i=1,size(buffer2,2)
                    do k=1,size(buffer2,1)
                    if ((buffer2(k,i)%x>=entrycell2(j)%bcell%xleft) .and. &
                        (buffer2(k,i)%x<entrycell2(j)%bcell%xright).and. &
                        (buffer2(k,i)%y>=entrycell2(j)%bcell%ybot).and. &
                        (buffer2(k,i)%y<entrycell2(j)%bcell%ytop))then
                        entrycell2(j)%bcell%ptot=entrycell2(j)%bcell%ptot+1
                        entrycell2(j)%bcell%plist(entrycell2(j)%bcell%ptot)=buffer2(k,i)
                        entrycell2(j)%bcell%plist(entrycell2(j)%bcell%ptot)%tid=3

                        entrycell2(j)%bcell%plist(entrycell2(j)%bcell%ptot)%mass=fmass*relden
                        entrycell2(j)%bcell%plist(entrycell2(j)%bcell%ptot)%density=rho*relden

                        entrycell2(j)%bcell%plist(entrycell2(j)%bcell%ptot)%oden=entrycell2(j)%bcell &
                        %plist(entrycell2(j)%bcell%ptot)%density
                        entrycell2(j)%bcell%plist(entrycell2(j)%bcell%ptot)%ovol=fmass/rho
                        entrycell2(j)%bcell%plist(entrycell2(j)%bcell%ptot)%con=0.50_dp

                    end if
                    end do
                end do
            end do
        !$omp end parallel do
        
    end subroutine insert_buffer2

    subroutine remove_buffer1()

        implicit none

        integer :: i,j,k

        !$omp parallel do default(shared) schedule(static)
        do i=1,entrycounter1

            entrycell1(i)%bcell%temfct=0

            do j=1,entrycell1(i)%bcell%ptot

                if ((.not.(entrycell1(i)%bcell%plist(j)%buffer)).or. &
                (entrycell1(i)%bcell%plist(j)%buffer.and.(entrycell1(i)%bcell%plist(j)%x<rv_buf_r))) then

                    entrycell1(i)%bcell%temfct=entrycell1(i)%bcell%temfct+1
                    entrycell1(i)%bcell%ftemp(entrycell1(i)%bcell%temfct)%part=>entrycell1(i)%bcell%plist(j)

                end if

            end do

            entrycell1(i)%bcell%ptot=entrycell1(i)%bcell%temfct

            do j=1,entrycell1(i)%bcell%temfct

                entrycell1(i)%bcell%plist(j)=entrycell1(i)%bcell%ftemp(j)%part

            end do

        end do
        !$omp end parallel do
                
    end subroutine remove_buffer1

    subroutine remove_buffer2()

        implicit none

        integer :: i,j,k

        !$omp parallel do default(shared) schedule(static)
        do i=1,entrycounter2

            entrycell2(i)%bcell%temfct=0

            do j=1,entrycell2(i)%bcell%ptot

                if ((.not.(entrycell2(i)%bcell%plist(j)%buffer)).or. &
                (entrycell2(i)%bcell%plist(j)%buffer.and.(entrycell2(i)%bcell%plist(j)%x>rv_buf_l))) then

                    entrycell2(i)%bcell%temfct=entrycell2(i)%bcell%temfct+1
                    entrycell2(i)%bcell%ftemp(entrycell2(i)%bcell%temfct)%part=>entrycell2(i)%bcell%plist(j)

                end if

            end do

            entrycell2(i)%bcell%ptot=entrycell2(i)%bcell%temfct

            do j=1,entrycell2(i)%bcell%temfct

                entrycell2(i)%bcell%plist(j)=entrycell2(i)%bcell%ftemp(j)%part

            end do

        end do
        !$omp end parallel do
                
    end subroutine remove_buffer2
    
end module internal