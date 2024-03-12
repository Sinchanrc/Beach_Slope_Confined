module wcsph

    use initialize
    use domain
    use kernel


    
    implicit none

    contains

    subroutine presstoden

        !$omp parallel do schedule(runtime) default(shared) private(i,k) collapse(2)
        do j=sx,ex
        do i=sy,ey
            if(dpcell(i,j)%ptot/=0) then
                
                do k=1,dpcell(i,j)%ptot

                        dpcell(i,j)%plist(k)%density=rho*(((dpcell(i,j)%plist(k)%pressure*t_gam/(rho*(co**2))) &
                        +1)**(1/t_gam))


                end do
                
            end if           
        end do
        end do
        !$omp end parallel do
        
    end subroutine presstoden

    subroutine dentopress
        
        !$omp parallel do schedule (runtime) default(shared) private(i,k) collapse(2)     
        do j=sx,ex
            do i=sy,ey
                if(dpcell(i,j)%ptot/=0) then

                do k=1,dpcell(i,j)%ptot

                    if (dpcell(i,j)%plist(k)%tid>2) then

                    dpcell(i,j)%plist(k)%pressure=(rho*(co**2))* &
                                        (((dpcell(i,j)%plist(k)%density/rho)**t_gam)-1)/t_gam

                    if(dpcell(i,j)%plist(k)%pressure<0.0d0) then
                        dpcell(i,j)%plist(k)%pressure=0.0d0
                        dpcell(i,j)%plist(k)%density=rho
                    end if

                    end if

                end do

                end if
            end do
        end do
        !$omp end parallel do
        
    end subroutine dentopress

    subroutine calmasscons
        
        !$omp parallel do schedule (runtime) private(term1,term2,m,i,k) default(shared) collapse(2)   
        do j=sx,ex
            do i=sy,ey
                if(dpcell(i,j)%ptot/=0) then

                do k=1,dpcell(i,j)%ptot
                if (dpcell(i,j)%plist(k)%tid>2) then
                    term1=0.0d0
                    term2=0.0d0

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))    

                            term2=(rho*((((((brrealy*distfac)*((2*bl)-1))+wc-dpcell(y,x)%plist(pp)%y)*rho*abs(g)* &
                            t_gam/(rho*(co**2)))+1)**(1/t_gam)))- &
                            (rho*((((((brrealy*distfac)*((2*bl)-1))+wc-dpcell(i,j)%plist(k)%y)*rho*abs(g)* &
                            t_gam/(rho*(co**2)))+1)**(1/t_gam)))


                            term1=term1-(dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%vx- &
                            dpcell(i,j)%plist(k)%vx)* &
                            (dpcell(i,j)%plist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                            dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%plist(k)%coff(2)* &
                            Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))) &
                            *(dpcell(i,j)%plist(k)%density/dpcell(y,x)%plist(pp)%density) &
                            +delt*2*h1*co*(dpcell(y,x)%plist(pp)%density-dpcell(i,j)%plist(k)%density-term2)* &
                            (dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)* &
                            (dpcell(i,j)%plist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                            dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%plist(k)%coff(2)* &
                            Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))* &
                            (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)/ &
                            (dpcell(i,j)%list(k)%dist(m)**2+lam)

                            term2=term2-(dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%vy- &
                            dpcell(i,j)%plist(k)%vy)* &
                            (dpcell(i,j)%plist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                            dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%plist(k)%coff(4)* &
                            Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))) &
                            *(dpcell(i,j)%plist(k)%density/dpcell(y,x)%plist(pp)%density) &
                            +delt*2*h1*co*(dpcell(y,x)%plist(pp)%density-dpcell(i,j)%plist(k)%density-term2)* &
                            (dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)* &
                            (dpcell(i,j)%plist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                            dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%plist(k)%coff(4)* &
                            Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))* &
                            (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)/ &
                            (dpcell(i,j)%list(k)%dist(m)**2+lam)
    
                        end associate
                    end do

                    dpcell(i,j)%plist(k)%tden=dpcell(i,j)%plist(k)%density+(term1+term2)*dt

                end if
                end do

                end if
            end do
        end do
        !$omp end parallel do
        
    end subroutine calmasscons

    subroutine updateden
        
        !$omp parallel do schedule(runtime) default(shared) private(i,k) collapse(2)
        do j=sx,ex
        do i=sy,ey
            if(dpcell(i,j)%ptot/=0) then
                
                do k=1,dpcell(i,j)%ptot

                if (dpcell(i,j)%plist(k)%tid>2) then

                        dpcell(i,j)%plist(k)%density=dpcell(i,j)%plist(k)%tden

                end if


                end do
                
            end if           
        end do
        end do
        !$omp end parallel do

    end subroutine updateden
    

end module wcsph
