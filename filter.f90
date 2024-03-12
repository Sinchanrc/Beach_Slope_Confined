module filter

    use particle
    use initialize
    use kernel
    use domain

    implicit none
    
    contains

    subroutine shep_fil_int(pm1)

        implicit none

        integer,intent(in) ::pm1
        
        !Density for internal particles
        !$omp parallel do schedule(runtime) private(k,m,term1,i) default(shared) collapse(2)
        do j=sx,ex
            do i=sy,ey
            if (dpcell(i,j)%ptot/=0) then
                
                do k=1,dpcell(i,j)%ptot
                if ((dpcell(i,j)%plist(k)%tid==pm1).or.(pm1==0)) then
                    if(dpcell(i,j)%list(k)%count/=0) then
                        dpcell(i,j)%plist(k)%density=0.0_dp
                        do m=1,dpcell(i,j)%list(k)%count                           
                            associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                                y=>dpcell(i,j)%list(k)%interlist(2,m), &
                                pp=>dpcell(i,j)%list(k)%interlist(3,m))
                            
                                dpcell(i,j)%plist(k)%density=dpcell(i,j)%plist(k)%density+ &
                                Wab(dpcell(i,j)%list(k)%dist(m),h1)*dpcell(y,x)%plist(pp)%mass

                            end associate
                        end do
                    end if

                    dpcell(i,j)%plist(k)%density=dpcell(i,j)%plist(k)%density+&
                    dpcell(i,j)%plist(k)%mass*Wo(h1)

                end if
                end do
                
            end if
            end do
        end do
        !$omp end parallel do

        !Shepherd filtering for density for fluid particles
        !$omp parallel do schedule(runtime) private(k,m,term1,i) default(shared) collapse(2)
        do j=sx,ex
            do i=sy,ey
            if (dpcell(i,j)%ptot/=0) then
                
                do k=1,dpcell(i,j)%ptot 
                if ((dpcell(i,j)%plist(k)%tid==pm1).or.(pm1==0))then                   
                    if(dpcell(i,j)%list(k)%count/=0) then                        
                        term1=0.0_dp
                        do m=1,dpcell(i,j)%list(k)%count                           
                            associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                                y=>dpcell(i,j)%list(k)%interlist(2,m), &
                                pp=>dpcell(i,j)%list(k)%interlist(3,m))
                            
                                term1=term1+dpcell(y,x)%plist(pp)%mass* &
                                Wab(dpcell(i,j)%list(k)%dist(m),h1)/dpcell(y,x)%plist(pp)%density

                            end associate
                        end do
                    end if

                    term1=term1+dpcell(i,j)%plist(k)%mass*Wo(h1)/dpcell(i,j)%plist(k)%density
                    dpcell(i,j)%pplist(k)%tden=dpcell(i,j)%plist(k)%density/term1

                end if
                end do
                
            end if
            end do
        end do
        !$omp end parallel do

        !$omp parallel do schedule(runtime) private(i,k) default(shared) collapse(2)
        do j=sx,ex
                do i=sy,ey
                    if (dpcell(i,j)%ptot/=0) then
                        
                        do k=1,dpcell(i,j)%ptot
                        if ((dpcell(i,j)%plist(k)%tid==pm1).or.(pm1==0)) then                    
                            dpcell(i,j)%plist(k)%density=dpcell(i,j)%pplist(k)%tden
                        end if
                        end do
                        
                    end if
                end do
        end do 
        !$omp end parallel do       
            
    end subroutine shep_fil_int

end module filter
