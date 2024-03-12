module matrix

    use solver
    use initialize
    use gradient
    use domain
    use kernel
    use functions
    use particle
    use interactions
    use search
    use setup
    use memory
    use mls


    implicit none
    
    contains

    subroutine initbounsol()

        implicit none

        integer :: sign

        count=0
        do j=sx,ex 
            do i=sy,ey
            if (dpcell(i,j)%btot/=0) then
                do cout=1,dpcell(i,j)%btot

                    if ((dpcell(i,j)%plist(cout)%tid<=2)) then
                    count=count+1
                    dpcell(i,j)%plist(cout)%matid=count
                    end if

                end do
            end if
            end do
        end do

        binmax=count

        allocate(bval(binmax*fplistmax),bcol(binmax*fplistmax),brow(binmax+1), &
        bvec(binmax),bmatrix(binmax),bsol(binmax),bguess(binmax))

        bguess=30000.0d0
        bsol=30000.0d0

        do k=1,binmax
            
            allocate(bmatrix(k)%val(fplistmax),bmatrix(k)%col(fplistmax))

        end do

        call compcorr(1,1)

        !$omp parallel do schedule(dynamic) private(i,j,k,m,term1,term2,sign) &
        !$omp default(shared) collapse(2) ! preparing the coff matrix for boundary part in CSR
            do j=sx,ex
                do i=sy,ey
                    if (dpcell(i,j)%btot/=0) then
                    
                    do k=1,dpcell(i,j)%btot
                        if (dpcell(i,j)%plist(k)%tid<=2) then
                        associate(pos=>dpcell(i,j)%plist(k)%matid,num=>dpcell(i,j)%list(k)%count)

                        bmatrix(pos)%sz=num+1
                        bmatrix(pos)%val(:)=0.0d0
                        bmatrix(pos)%col(:)=1
                        bvec(pos)=0.0d0
                        num2=dpcell(i,j)%list(k)%count
                            do m=1,num2
                            associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                                y=>dpcell(i,j)%list(k)%interlist(2,m), &
                                pp=>dpcell(i,j)%list(k)%interlist(3,m))

                                term1=8*dpcell(y,x)%plist(pp)%mass*(dpcell(i,j)%plist(k)%coff(1)* &
                                Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)+ &
                                dpcell(i,j)%plist(k)%coff(2)* &
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))* &
                                (dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/(dpcell(i,j)%list(k)%dist(m))**2

                                term2=8*dpcell(y,x)%plist(pp)%mass*(dpcell(i,j)%plist(k)%coff(3)* &
                                Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)+ &
                                dpcell(i,j)%plist(k)%coff(4)* &
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))* &
                                (dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/(dpcell(i,j)%list(k)%dist(m))**2


                                ! if (dpcell(y,x)%plist(pp)%tid<=2) then
                                ! bmatrix(pos)%val(m)=-(term1+term2)/(dpcell(y,x)%plist(pp)%density+ &
                                ! dpcell(i,j)%plist(k)%density)**2
                                ! bmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid
                                ! else
                                ! bvec(pos)=bvec(pos)+dpcell(y,x)%plist(pp)%pressure*(term1+term2)/&
                                ! (dpcell(y,x)%plist(pp)%density+dpcell(i,j)%plist(k)%density)**2
                                ! end if
                                ! bmatrix(pos)%val(num+1)=bmatrix(pos)%val(num+1)+(term1+term2)/&
                                ! (dpcell(y,x)%plist(pp)%density+dpcell(i,j)%plist(k)%density)**2
                                

                            !     if (dpcell(i,j)%plist(k)%ynorm==1) then

                            !     if ((dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)>0.0d0) then
                            !     sign=1
                            !     else
                            !     sign=-1
                            !     end if

                            !     if (.not.(dpcell(i,j)%plist(k)%mirror)) then
                            !         sign=0
                            !     end if

                            !     bvec(pos)=bvec(pos)+sign*g*abs((dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y))*rho &
                            !             *(term1+term2)/(dpcell(y,x)%plist(pp)%density+dpcell(i,j)%plist(k)%density)**2

                            !     end if


                            end associate
                            end do
                        bmatrix(pos)%col(num+1)=pos

                        term1=bmatrix(pos)%val(num+1)
                        if (term1<1e-6) then
                            term1=1.0d0
                        end if
                        bvec(pos)=bvec(pos)/term1   ! Jacobi Preconditioning
                        do m=1,fplistmax!num+1
                        bmatrix(pos)%val(m)=bmatrix(pos)%val(m)/term1 !Jacobi Preconditioning
                        end do                 
                        end associate
                        end if                
                    end do

                    end if
                end do
            end do
        !$omp end parallel do

        ! call format(bmatrix,bval,brow,bcol,binmax)

        ! call bicgstab(tl,bguess,binmax,bval,brow,bcol,bvec,bsol)

        ! Assigning boundary pressures
        !$omp parallel do schedule(dynamic) private(i,j,k) default(shared) collapse(2)  
            do j=sx,ex
            do i=sy,ey
                if (dpcell(i,j)%btot/=0) then
                do k=1,dpcell(i,j)%btot
                    if (dpcell(i,j)%plist(k)%tid<=2) then
                    dpcell(i,j)%plist(k)%pressure=bsol(dpcell(i,j)%plist(k)%matid)
                    end if
                end do
                end if
            end do
            end do
        !$omp end parallel do   

        deallocate(bval,bcol,brow,bvec,bmatrix,bsol,bguess)

        ! end if

    end subroutine initbounsol
    
end module matrix