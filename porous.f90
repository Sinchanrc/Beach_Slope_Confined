module porous

    use particle
    use kernel
    use domain
    use initialize

    implicit none

    contains

    pure function vis(oden) result(dmu)

        real(dp),intent(in) :: oden
        real(dp),parameter :: tem=20.0_dp
        real(dp) :: dmu,Acn,Bcn,Sp

        Sp=oden-1000.0_dp
        ! Acn=1.474_dp*(1e-3)+1.5_dp*(1e-5)*tem-3.927_dp*(1e-8)*tem**2
        ! Bcn=1.073_dp*(1e-5)-8.5_dp*(1e-8)*tem+2.230_dp*(1e-10)*tem**2
        Acn=0.001652_dp
        Bcn=0.0000083_dp
        ! dmu=exp(-0.00379418_dp+(0.604129_dp/(139.18_dp+tem)))*(1e-3)* &
        !     (1.0_dp + Acn*Sp + Bcn*Sp**2)

        dmu=(1.002_dp+Acn*Sp+Bcn*Sp**2)*0.001_dp



    end function

    subroutine porsearch(list,n,cell1,cell2,cell3,cell4,smln,solidfrac)
        use functions
        implicit none
        type(verlet),intent(inout) :: list
        type(cell),intent(inout) :: cell1,cell2,cell3,cell4
        ! type(cell) :: blockset(4)
        real(dp),intent(in) :: smln
        integer,intent(in) :: n
        integer ::counter1,counter2,intct
        real(dp) :: radial,redfac
        real(dp),intent(out) :: solidfrac


        solidfrac=0.0_dp
        redfac=1.0_dp !0.9830_dp

            do counter2=1,cell1%porct

                radial=dist(cell1%plist(n),cell1%porlist(counter2))

                if((radial<(2*smln)))then

                    solidfrac=solidfrac+(1.0_dp-por*redfac)*Wab(radial,h1)* &
                    cell1%porlist(counter2)%mass/cell1%porlist(counter2)%density


                end if

            end do

            do counter2=1,cell2%porct

                radial=dist(cell1%plist(n),cell2%porlist(counter2))

                if((radial<(2*smln)))then

                    solidfrac=solidfrac+(1.0_dp-por*redfac)*Wab(radial,h1)* &
                    cell2%porlist(counter2)%mass/cell2%porlist(counter2)%density



                end if

            end do

            do counter2=1,cell3%porct

                radial=dist(cell1%plist(n),cell3%porlist(counter2))

                if((radial<(2*smln)))then

                    solidfrac=solidfrac+(1.0_dp-por*redfac)*Wab(radial,h1)* &
                    cell3%porlist(counter2)%mass/cell3%porlist(counter2)%density


                end if

            end do

            do counter2=1,cell4%porct

                radial=dist(cell1%plist(n),cell4%porlist(counter2))

                if((radial<(2*smln)))then

                    solidfrac=solidfrac+(1.0_dp-por*redfac)*Wab(radial,h1)* &
                    cell4%porlist(counter2)%mass/cell4%porlist(counter2)%density


                end if

            end do

    end subroutine porsearch

    subroutine effpor
    
        implicit none

        real(dp) :: solidfrac,lamk=650.0_dp,d50=0.001020_dp,Fch,kper !299.71_dp

        type(cell),pointer :: cell1,cell2,cell3,cell4

        integer :: i,j,k,m

        !$omp do schedule (runtime) &
        !$omp private(m,k,i,j,solidfrac,Fch,kper,cell1,cell2,cell3,cell4) collapse(2)
        do j=sx,ex
            do i=sy,ey                

            do k=1,dpcell(i,j)%ptot

                solidfrac=0.0_dp
                cell1=>dpcell(i,j)
                dpcell(i,j)%pplist(k)%resistx=0.0_dp
                dpcell(i,j)%pplist(k)%resisty=0.0_dp

                if (dpcell(i,j)%plist(k)%tid==3) then

                dpcell(i,j)%pplist(k)%porosity=1.0_dp

                if(dpcell(i,j)%plist(k)%x<=(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y<=(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then

                    cell2=>dpcell(i,j-1)
                    cell3=>dpcell(i+1,j)
                    cell4=>dpcell(i+1,j-1)
    

                else if(dpcell(i,j)%plist(k)%x>(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y<=(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    cell2=>dpcell(i,j+1)
                    cell3=>dpcell(i+1,j)
                    cell4=>dpcell(i+1,j+1)

                else if(dpcell(i,j)%plist(k)%x>(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y>(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then

                    cell2=>dpcell(i-1,j+1)
                    cell3=>dpcell(i-1,j)
                    cell4=>dpcell(i,j+1)

                else if(dpcell(i,j)%plist(k)%x<=(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y>(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then

                    cell2=>dpcell(i-1,j)
                    cell3=>dpcell(i-1,j-1)
                    cell4=>dpcell(i,j-1)

                else if(dpcell(i,j)%plist(k)%x==(dpcell(i,j)%xleft+ &
                    abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                    dpcell(i,j)%plist(k)%y==(dpcell(i,j)%ybot+ &
                    abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then

                    cell2=>dpcell(i-1,j)
                    cell3=>dpcell(i-1,j-1)
                    cell4=>dpcell(i,j-1)

                end if

                call porsearch(dpcell(i,j)%list(k),k,cell1,cell2,cell3,cell4,h1,solidfrac)

                if (solidfrac>(1.0_dp-por)) then

                        solidfrac=1.0_dp-por

                end if

                dpcell(i,j)%pplist(k)%porosity=1.0_dp-solidfrac
                dpcell(i,j)%plist(k)%density=dpcell(i,j)%pplist(k)%porosity* &
                                            dpcell(i,j)%plist(k)%oden

                Fch=1.750_dp/(sqrt(150*dpcell(i,j)%pplist(k)%porosity**3))

                ! dpcell(i,j)%pplist(k)%resistx=-((mu/dpcell(i,j)%plist(k)%oden)*dpcell(i,j)%plist(k)%vx*lamk* &
                ! ((1.0_dp-dpcell(i,j)%pplist(k)%porosity)**2)/((dpcell(i,j)%pplist(k)%porosity**3)*(d50**2))) !&
                ! ! -((Fch*dpcell(i,j)%plist(k)%vx*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)* &
                ! ! sqrt(lamk)*(1.0_dp-dpcell(i,j)%pplist(k)%porosity))/(sqrt(dpcell(i,j)%pplist(k)%porosity**3)*d50))

                ! dpcell(i,j)%pplist(k)%resisty=-((mu/dpcell(i,j)%plist(k)%oden)*dpcell(i,j)%plist(k)%vy*lamk* &
                ! ((1.0_dp-dpcell(i,j)%pplist(k)%porosity)**2)/((dpcell(i,j)%pplist(k)%porosity**3)*(d50**2))) !&
                ! ! -((Fch*dpcell(i,j)%plist(k)%vy*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)* &
                ! ! sqrt(lamk)*(1.0_dp-dpcell(i,j)%pplist(k)%porosity))/(sqrt(dpcell(i,j)%pplist(k)%porosity**3)*d50))

                dpcell(i,j)%pplist(k)%resistx=-((vis(dpcell(i,j)%plist(k)%oden)/dpcell(i,j)%plist(k)%oden)*lamk* &
                ((1.0_dp-dpcell(i,j)%pplist(k)%porosity)**2)/((dpcell(i,j)%pplist(k)%porosity**3)*(d50**2))) &
                -((Fch*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)* &
                sqrt(lamk)*(1.0_dp-dpcell(i,j)%pplist(k)%porosity))/(sqrt(dpcell(i,j)%pplist(k)%porosity**3)*d50))

                dpcell(i,j)%pplist(k)%resisty=-((vis(dpcell(i,j)%plist(k)%oden)/dpcell(i,j)%plist(k)%oden)*lamk* &
                ((1.0_dp-dpcell(i,j)%pplist(k)%porosity)**2)/((dpcell(i,j)%pplist(k)%porosity**3)*(d50**2))) &
                -((Fch*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)* &
                sqrt(lamk)*(1.0_dp-dpcell(i,j)%pplist(k)%porosity))/(sqrt(dpcell(i,j)%pplist(k)%porosity**3)*d50))



                end if

            end do

            end do
        end do
        !$omp end do
        
    end subroutine effpor
    

end module porous
