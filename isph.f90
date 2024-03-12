module isph

    use particle
    use initialize
    use domain
    use kernel
    use functions
    use solver
    use output


    implicit none

    contains

    subroutine projection
        implicit none

        integer :: i,j,k,m
        
        ! Intermediate fluid position calc
        !$omp do schedule(runtime) collapse(2) private(i,j,k)
            do j=sx,ex
                do i=sy,ey
                if (dpcell(i,j)%ptot/=0) then
                
                    do k=1,dpcell(i,j)%ptot

                    if (dpcell(i,j)%plist(k)%tid==3) then
                    dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%x
                    dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+(dt* &
                    dpcell(i,j)%plist(k)%vx/dpcell(i,j)%pplist(k)%porosity)
                    dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%y
                    dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+(dt* &
                    dpcell(i,j)%plist(k)%vy/dpcell(i,j)%pplist(k)%porosity)
                    end if
                    end do

                end if
                end do
            end do
        !$omp end do        
        

        
    end subroutine projection
    
    subroutine freesurf
        implicit none 

        real(dp) :: heff,t1,t2
        integer :: i,j,k,m


        ! Free surface(x%tid==3) then
        !$omp do private(m,i,k,j,heff) schedule (runtime) collapse(2)    
        do j=sx,ex
            do i=sy,ey
                do k=1,dpcell(i,j)%ptot
                ! Free Surface Identification
                if ((dpcell(i,j)%plist(k)%tid==3).and. &
                (.not.(dpcell(i,j)%plist(k)%buffer))) then

                heff=h1
                if(dpcell(i,j)%pplist(k)%gradvx<=(lamfs*maxdivr)) then
                
                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%nh(m)%part)
    
                        x%vicinity=.true.

    
                        end associate
                    end do

                end if

                end if

                end do
            end do
        end do
        !$omp end do
        
    end subroutine freesurf

    subroutine resetid
        implicit none 

        integer :: mat_count,t1,i,j

        !$omp single

        mat_count=0

        do j=1,cellx 
            do i=1,celly
            if (dpcell(i,j)%ptot/=0) then
                do t1=1,dpcell(i,j)%ptot

                    mat_count=mat_count+1
                    dpcell(i,j)%plist(t1)%matid=mat_count

                end do
            end if
            end do
        end do

        !$omp end single
    
        
    end subroutine resetid

    subroutine ppesolve

        implicit none

        integer :: i,j,k,m,dia_pt,ins_pt

        real(dp) :: lamp,lamp2,t1,t2,rho_i,rho_j,pvol,p_dist,rho_ij,dia_term

        !$omp parallel do schedule(runtime) default(shared)
        do i=1,finmax 
            fmatrix(i)%sz=1
            fmatrix(i)%val(:)=0.0_dp
            fmatrix(i)%val(1)=1.0_dp
            fmatrix(i)%col(:)=i
            fvec(i)=0.0_dp
        end do
        !$omp end parallel do

        
        
        ! Preparing the coff matrix for fluid part in CSR
        !$omp parallel do schedule(runtime) default(shared) &
        !$omp private(m,t1,t2,k,i,j,lamp,lamp2,rho_i,rho_j,pvol,p_dist,rho_ij,dia_pt,ins_pt,dia_term) &
        !$omp collapse(2)
        do j=sx,ex
            do i=sy,ey
            ! if (dpcell(i,j)%ptot/=0) then

                

                do k=1,dpcell(i,j)%ptot

                dia_pt=1
                ins_pt=0
                dia_term=0.0_dp
                        
                ! if (dpcell(i,j)%plist(k)%tid/=4) then
                associate(pos=>dpcell(i,j)%plist(k)%matid,num=>dpcell(i,j)%list(k)%count)
                    t1=0.0_dp
                    t2=0.0_dp
                    fmatrix(pos)%sz=num+1
                    ! fmatrix(pos)%val(:)=0.0_dp
                    ! fmatrix(pos)%col(:)=pos
                    fvec(pos)=0.0_dp
                    num2=dpcell(i,j)%list(k)%count
                    rho_i=dpcell(i,j)%plist(k)%density*dpcell(i,j)%pplist(k)%porosity**(-1)
                    if((num/=0)) then   !.or.(dpcell(i,j)%plist(k)%free/=1)
                    do m=1,num2
                        associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                            y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                            z=>dpcell(i,j)%list(k)%klt)

                        if (pos>x%matid) then

                            dia_pt=m+1
                            ins_pt=m

                        else

                            ins_pt=m+1

                        end if

                            ! if (x%tid/=4) then
                        

                            if (dpcell(i,j)%plist(k)%free) then

                                lamp=(1.0_dp+x%oden/(rho*rel_den))**(-1)
                                ! lamp=0.50_dp

                            ! elseif((pll<dpcell(i,j)%pplist(k)%gradvx).and.(dpcell(i,j)%pplist(k)%gradvx<pul)) then

                            ! lamp=((1.0_dp+x%density/ &
                            ! (y%porosity*rho*1.0))**(-1)) &
                            ! *(1.0_dp-cos(22.0_dp*(dpcell(i,j)%pplist(k)%gradvx-pll)/ &
                            !                 ((pul-pll)*7.0_dp)))
                            else 

                                lamp =1.0_dp

                            end if

                            ! rho_i=dpcell(i,j)%plist(k)%density*dpcell(i,j)%pplist(k)%porosity**(-2)
                            rho_j=x%density*y%porosity**(-1)
                            pvol=x%mass*x%density**(-1)
                            p_dist=((dist(x,dpcell(i,j)%plist(k))**2)+lam)**(-1)

                            rho_ij=2*x%density*dpcell(i,j)%plist(k)%density* &
                                (((y%porosity)*(dpcell(i,j)%pplist(k)%porosity))**(-1))* &
                                ((x%density*y%porosity**(-1))+ &
                                (dpcell(i,j)%plist(k)%density*dpcell(i,j)%pplist(k)%porosity**(-1)))**(-1)

                            if (dpcell(i,j)%plist(k)%tid==3) then

                                t1=2*pvol*(&
                                (dpcell(i,j)%pplist(k)%coff(1)*z(1,m)+dpcell(i,j)%pplist(k)%coff(2)* &
                                z(2,m))) &
                                *(dpcell(i,j)%plist(k)%x-x%x)*(p_dist)

                                t2=2*pvol*(&
                                (dpcell(i,j)%pplist(k)%coff(3)*z(1,m)+dpcell(i,j)%pplist(k)%coff(4)* &
                                z(2,m))) &
                                *(dpcell(i,j)%plist(k)%y-x%y)*(p_dist)
                                


                                ! if (.not.(x%buffer)) then

                                ! fmatrix(pos)%val(m)=(-(t1+t2))*lamp
                                ! fmatrix(pos)%col(m)=x%matid

                                fmatrix(pos)%val(ins_pt)=(-(t1+t2))*lamp
                                fmatrix(pos)%col(ins_pt)=x%matid

                                ! else 
                                !     fvec(pos)=fvec(pos)+(t1+t2)*lamp*x%pressure

                                ! end if

                                ! fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(t1+t2)
                                dia_term=dia_term+(t1+t2)


                                if (x%tid==3) then


                                t1=(pvol*(x%vxs*rho_ij/y%porosity- &
                                dpcell(i,j)%plist(k)%vxs*rho_ij/dpcell(i,j)%pplist(k)%porosity)* &
                                (dpcell(i,j)%pplist(k)%coff(1)*z(1,m)+dpcell(i,j)%pplist(k)%coff(2)* &
                                z(2,m)))


                                t2=(pvol*(x%vys*rho_ij/y%porosity- &
                                dpcell(i,j)%plist(k)%vys*rho_ij/dpcell(i,j)%pplist(k)%porosity)* &
                                (dpcell(i,j)%pplist(k)%coff(3)*z(1,m)+dpcell(i,j)%pplist(k)%coff(4)* &
                                z(2,m)))

                                else 

                                    t1=(pvol*(x%vxs*rho_i/dpcell(i,j)%pplist(k)%porosity- &
                                    dpcell(i,j)%plist(k)%vxs*rho_i/dpcell(i,j)%pplist(k)%porosity)* &
                                    (dpcell(i,j)%pplist(k)%coff(1)*z(1,m)+dpcell(i,j)%pplist(k)%coff(2)* &
                                    z(2,m)))


                                    t2=(pvol*(x%vys*rho_i/dpcell(i,j)%pplist(k)%porosity- &
                                    dpcell(i,j)%plist(k)%vys*rho_i/dpcell(i,j)%pplist(k)%porosity)* &
                                    (dpcell(i,j)%pplist(k)%coff(3)*z(1,m)+dpcell(i,j)%pplist(k)%coff(4)* &
                                    z(2,m)))

                                end if

                                fvec(pos)=fvec(pos)+(t1+t2)*lamp/real(dt,dp) !lamp



                            elseif (dpcell(i,j)%plist(k)%tid==2) then


                                t1=2*pvol*(&
                                z(1,m)) &
                                *(dpcell(i,j)%plist(k)%x-x%x)*(p_dist)

                                t2=2*pvol*(&
                                z(2,m)) &
                                *(dpcell(i,j)%plist(k)%y-x%y)*(p_dist)
                                
                                if (x%tid<=2) then

                                    ! fmatrix(pos)%val(m)=-(t1+t2)!*y%lamp
                                    ! fmatrix(pos)%col(m)=x%matid
                                    ! fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(t1+t2)!*y%lamp
    
                                    fmatrix(pos)%val(ins_pt)=-(t1+t2)!*y%lamp
                                    fmatrix(pos)%col(ins_pt)=x%matid
                                    dia_term=dia_term+(t1+t2)

                                end if


                                ! fmatrix(pos)%val(1)=-1.0_dp!(t1+t2)
                                ! fmatrix(pos)%col(1)=dpcell(i,j)%plist(k)%wall(1)!x%matid
                                ! fmatrix(pos)%val(num+1)=1.0_dp!fmatrix(pos)%val(num+1)+(t1+t2)

                                ! end if

                            elseif (dpcell(i,j)%plist(k)%tid==1) then

                                if (x%tid==3) then

                                    if (x%free) then

                                        lamp2=(1.0_dp+dpcell(i,j)%plist(k)%oden/(rho*rel_den))**(-1)
                                        ! lamp2=0.50_dp
        
                                    ! elseif((pll<y%gradvx).and.(y%gradvx<pul)) then
        
                                    ! lamp2=((1.0_dp+dpcell(i,j)%plist(k)%density/(rho*1.0))**(-1)) &
                                    ! *(1.0_dp-cos(22.0_dp*(y%gradvx-pll)/ &
                                    !                 ((pul-pll)*7.0_dp)))
                                    else 
        
                                        lamp2 =1.0_dp
        
                                    end if

                                    if (x%buffer) then

                                        lamp2=1.0_dp

                                    end if

                                t1=2*dpcell(i,j)%plist(k)%mass*(&
                                (y%coff(1)*z(1,m)+y%coff(2)* &
                                z(2,m))) &
                                *(dpcell(i,j)%plist(k)%x-x%x)*(p_dist)/ &
                                ((dpcell(i,j)%plist(k)%density))

                                t2=2*dpcell(i,j)%plist(k)%mass*(&
                                (y%coff(3)*z(1,m)+y%coff(4)* &
                                z(2,m))) &
                                *(dpcell(i,j)%plist(k)%y-x%y)*(p_dist)/ &
                                ((dpcell(i,j)%plist(k)%density))

                                ! t1=2*pvol*(&
                                ! Wabx(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(dpcell(i,j)%plist(k)%x-x%x)*(p_dist)

                                ! t2=2*pvol*(&
                                ! Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(dpcell(i,j)%plist(k)%y-x%y)*(p_dist)

                                ! if (.not.(x%buffer)) then

                                ! fmatrix(pos)%val(m)=-(t1+t2)*lamp2!*y%lamp
                                ! fmatrix(pos)%col(m)=x%matid

                                fmatrix(pos)%val(ins_pt)=-(t1+t2)*lamp2!*y%lamp
                                fmatrix(pos)%col(ins_pt)=x%matid
                                ! else 
                                !     fvec(pos)=fvec(pos)+(t1+t2)*lamp2*x%pressure
                                ! end if

                                ! fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(t1+t2)!*y%lamp
                                dia_term=dia_term+(t1+t2)

                                    t1=(dpcell(i,j)%plist(k)%mass*(x%vxs*rho_j/y%porosity- &
                                    dpcell(i,j)%plist(k)%vxs*rho_j/dpcell(i,j)%pplist(k)%porosity)* &
                                    (y%coff(1)*z(1,m)+y%coff(2)* &
                                    z(2,m)))/(dpcell(i,j)%plist(k)%density)


                                    t2=(dpcell(i,j)%plist(k)%mass*(x%vys*rho_j/y%porosity- &
                                    dpcell(i,j)%plist(k)%vys*rho_j/dpcell(i,j)%pplist(k)%porosity)* &
                                    (y%coff(3)*z(1,m)+y%coff(4)* &
                                    z(2,m)))/(dpcell(i,j)%plist(k)%density)

                                    ! if (.not.(x%buffer)) then

                                    fvec(pos)=fvec(pos)+(t1+t2)*lamp2/real(dt,dp)  !lamp

                                    ! end if

                                
                                elseif (x%tid<=2) then

                                    ! t1=2*x%mass*(&
                                    ! Wabx(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                    ! *(dpcell(i,j)%plist(k)%x-x%x)/ &
                                    ! (((dist(x,dpcell(i,j)%plist(k))**2)+lam)* &
                                    ! (x%density*dpcell(i,j)%plist(k)%density))

                                    ! t2=2*x%mass*(&
                                    ! Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                    ! *(dpcell(i,j)%plist(k)%y-x%y)/ &
                                    ! (((dist(x,dpcell(i,j)%plist(k))**2)+lam)* &
                                    ! (x%density*dpcell(i,j)%plist(k)%density))

                                    t1=2*pvol*(&
                                    z(1,m)) &
                                    *(dpcell(i,j)%plist(k)%x-x%x)*(p_dist)

                                    t2=2*pvol*(&
                                    z(2,m)) &
                                    *(dpcell(i,j)%plist(k)%y-x%y)*(p_dist)

                                    ! fmatrix(pos)%val(m)=-(t1+t2)!*y%lamp
                                    ! fmatrix(pos)%col(m)=x%matid
                                    ! fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(t1+t2)!*y%lamp
                                    
                                    fmatrix(pos)%val(ins_pt)=-(t1+t2)!*y%lamp
                                    fmatrix(pos)%col(ins_pt)=x%matid
                                    dia_term=dia_term+(t1+t2)

                                end if
                            
                            end if

                        ! end if


                        end associate
                    end do
                    end if

                    ! fmatrix(pos)%col(num+1)=pos

                    fmatrix(pos)%val(dia_pt)=dia_term
                    fmatrix(pos)%col(dia_pt)=pos



                    ! if ((num==0)) then  !.or.( dpcell(i,j)%plist(k)%free==1)
                    ! fmatrix(pos)%val(dia_pt)=1.0_dp
                    ! end if

                    t1=fmatrix(pos)%val(dia_pt)

                    ! if (dpcell(i,j)%plist(k)%tid/=2) then

                    fvec(pos)=fvec(pos)/real(t1,dp)   ! Jacobi Preconditioning
                    
                    do m=1,ceiling(fac2*fplistmax)
                    fmatrix(pos)%val(m)=fmatrix(pos)%val(m)/real(t1,dp) !Jacobi Preconditioning
                    end do 

                    ! end if


                end associate

                ! end if

                end do

            ! end if            
            end do
        end do
        !$omp end parallel do

        ! call format(fmatrix,fval,frow,fcol,finmax)
        ! call bicgstab(tl,fguess,finmax,fval,frow,fcol,fvec,fsol)
        call fgmres

        ! call pardisosolver(fval,frow,fcol,fvec,fsol)
        
        ! Assigning pressures
        !$omp parallel do schedule(runtime) default(shared) private(i,k,j) collapse(2)       
            do j=sx,ex
                do i=sy,ey
                    if (dpcell(i,j)%ptot/=0) then

                        do k=1,dpcell(i,j)%ptot

                        ! if (.not.(dpcell(i,j)%plist(k)%buffer)) then
                        
                        dpcell(i,j)%plist(k)%pressure=fsol(dpcell(i,j)%plist(k)%matid)

                        ! end if

                        end do

                    end if
                end do
            end do
        !$omp end parallel do 

    end subroutine ppesolve

end module isph



