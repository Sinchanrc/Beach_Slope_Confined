!------------------------------------------------------------------------------
! IITKGP, PhD
!------------------------------------------------------------------------------
!
! MODULE:  Functions
!
!> @author
!> Sinchan Roy Chowdhury}
!
! DESCRIPTION: 
!> Functions to calculate particle distance and matrix inversions
!
! REVISION HISTORY:
! 
! 
!------------------------------------------------------------------------------

module functions

    use particle
    use domain

    implicit none
    
    contains

    !Calculating the distance between 2 particles
    function dist(p1,p2) result(res2) 
        implicit none
        real(dp) :: res2
        type(particles),intent(in) :: p1,p2

        res2=0.0_dp

        associate(x1=>p1%x,x2=>p2%x,y1=>p1%y, &
                    y2=>p2%y)

            res2=sqrt(((x1-x2)**2)+((y1-y2)**2))
            
        end associate

    end function dist

    !Calculating the distance between a interpolation particle and internal particle
    ! function bpdist(p1,p2) result(res2) 
    !     implicit none
    !     real(dp) :: res2
    !     type(particles),intent(in) :: p1,p2

    !     res2=0.0_dp

    !     associate(x1=>p1%x,x2=>p2%x,y1=>p1%y, &
    !                 y2=>p2%y)

    !         res2=sqrt(((x1+p1%posshift(1)-x2)**2)+((y1+p1%posshift(2)-y2)**2))
            
    !     end associate

    ! end function

    ! Invert 2D matrix
    pure subroutine invertmat2D(array) 

        use,intrinsic :: ieee_arithmetic
        implicit none
        real(dp),dimension(4),intent(inout) :: array
        real(dp):: arr(2,2),det
            det=0.0_dp
            arr(:,:)=0.0_dp
            arr(1,1)=array(4)
            arr(1,2)=-array(2)
            arr(2,1)=-array(3)
            arr(2,2)=array(1)
            det= array(1)*array(4)-array(2)*array(3)
            if ((abs(det)>=1e-5).and.(.not.(ieee_is_nan(det))) &
            .and.(ieee_is_finite(det))) then
            arr(:,:)=arr(:,:)/real(det,dp)
            array(1)=arr(1,1)
            array(2)=arr(1,2)
            array(3)=arr(2,1)
            array(4)=arr(2,2)
            else

            array(1)=1.0_dp
            array(2)=0.0_dp
            array(3)=0.0_dp
            array(4)=1.0_dp
            end if

            ! array(1)=1.0_dp
            ! array(2)=0.0_dp
            ! array(3)=0.0_dp
            ! array(4)=1.0_dp


    end subroutine

    !Invert a 3D matrix
    subroutine matinv3(A,err) 
        implicit none
        !! Performs a direct calculation of the inverse of a 3×3 matrix.
        real(dp) , intent(inout) :: A(3,3)   !! Matrix
        real(dp)              :: B(3,3)   !! Inverse matrix
        real(dp)           :: detinv
        logical,intent(out) :: err

        ! Calculate the inverse determinant of the matrix
        detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

        if (detinv<1e-6) then
            err=.true.
            return
        else
            err=.false.
        end if

        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

        A(:,:)=B(:,:)

    end subroutine

    subroutine invertmat3D(array,err) ! Invert 3D matrix

        use,intrinsic :: ieee_arithmetic
        implicit none
        real(dp),dimension(3,3),intent(inout) :: array
        real(dp):: arr(3,3),det
        logical,intent(out) :: err

            det= array(1,1)*(array(2,2)*array(3,3)-array(3,2)*array(2,3)) &
            -array(1,2)*(array(2,1)*array(3,3)-array(3,1)*array(2,3)) &
            +array(1,3)*(array(2,1)*array(3,2)-array(3,1)*array(2,2))

            if ((det<1e-5).or.(ieee_is_nan(det)).or.(.not.(ieee_is_finite(det)))) then
            err=.true.
            array=0._dp
            return
            else
            err=.false.
            end if

            arr(1,1)=((array(2,2)*array(3,3)-array(3,2)*array(2,3)))
            arr(2,1)=-((array(1,2)*array(3,3)-array(3,2)*array(1,3)))
            arr(3,1)=((array(1,2)*array(2,3)-array(2,2)*array(1,3)))
            arr(1,2)=-((array(2,1)*array(3,3)-array(3,1)*array(2,3)))
            arr(2,2)=((array(1,1)*array(3,3)-array(3,1)*array(1,3)))
            arr(3,2)=-((array(1,1)*array(2,3)-array(2,1)*array(1,3)))
            arr(1,3)=((array(2,1)*array(3,2)-array(3,1)*array(2,2)))
            arr(2,3)=-((array(1,1)*array(3,2)-array(3,1)*array(1,2)))
            arr(3,3)=((array(2,2)*array(1,1)-array(2,1)*array(1,2)))

            array=arr/det


    end subroutine

    subroutine push_part(inputpt,dcell,mass,den,vel )
        use initialize
        implicit none

        type(buffer),intent(in) :: inputpt
        type(cell),intent(inout) :: dcell
        real(dp),intent(in) :: mass,den,vel 

        dcell%ptot=dcell%ptot+1
        dcell%plist(dcell%ptot)%mass=mass
        dcell%plist(dcell%ptot)%density=den
        dcell%plist(dcell%ptot)%x=inputpt%x
        dcell%plist(dcell%ptot)%y=inputpt%y
        dcell%plist(dcell%ptot)%pressure=0.0_dp
        dcell%plist(dcell%ptot)%vx=vel
        dcell%plist(dcell%ptot)%vy=0.0_dp
        dcell%plist(dcell%ptot)%domain=.false.
        dcell%plist(dcell%ptot)%range=.true.
        dcell%plist(dcell%ptot)%buffer=.true.
        dcell%plist(dcell%ptot)%tid=3
        dcell%plist(dcell%ptot)%pid=reserve%tank(reserve%si)
        reserve%si=reserve%si-1
        
    end subroutine push_part

    subroutine pull_part(dcell)
        use initialize

        implicit none

        integer :: i

        type(cell),intent(in) :: dcell

        do i=1,dcell%elist
        reserve%si=reserve%si+1
        reserve%tank(reserve%si)=dcell%exitlist(i)
        end do
        
    
        
    end subroutine pull_part
    
end module functions