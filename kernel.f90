!------------------------------------------------------------------------------
! IITKGP, PhD
!------------------------------------------------------------------------------
!
! MODULE: Kernel and Kernel Derivatives
!
!> @author
!> Sinchan Roy Chowdhury}
!
! DESCRIPTION: 
!> Functions to calculate kernel values and their derivatives for particle pairs
!  and Moving least squares values for ghost boundary
!
! REVISION HISTORY:
! 
! 
!------------------------------------------------------------------------------




module kernel

    use particle
    use domain
    use initialize

    implicit none

    real(dp),parameter :: alp=0.1_dp
    
    contains
    
    function Wab(radial,smln) result(res2)
        implicit none
        real(dp) :: res2,t1,t2
        real(dp),intent(in) :: smln
        real(dp),intent(in) :: radial 

        res2=0.0_dp

        ! res2=(7.0_dp/(4*3.140_dp*(smln**2)))*(1+(2*radial/smln))* & !Wendland c2
        !         (1.0_dp-(radial/(2*smln)))**4

        ! res2=(9.0_dp/(4*3.140_dp*(smln**2)))*(1+(3*radial/smln)+ &
        ! (35.0_dp/12)*((radial/smln)**2))* & !Wendland c4
        ! (1.0_dp-(radial/(2*smln)))**6

        ! res2=(2/(3.14*smln**2))*((3*(radial/(4*smln))**2)- &
        ! (3*radial/(4*smln))+3/4.0d0)          !Quadratic

        !Hybrid Kernel 

        t1=merge((4.0_dp-(6.0*(radial**2)/smln**2)+ &
        (3.0*(radial**3)/smln**3)),(2.0_dp-radial/smln)**3, &
        (radial<smln)) !Cubic spline

        t2=merge(((radial/smln)**3-(6.0*radial/smln)+6.0_dp), &
        ((2.0_dp-radial/smln)**3),(radial<smln)) !Hyperbolic

        res2=alp*t1*(10._dp/(28.0*3.14*smln**2)) + &
            (1.0_dp-alp)*t2*(1.0_dp/(3.0*3.14*smln**2))

        

    end function Wab

    function Wo(smln) result(res2)
        implicit none
        real(dp) :: res2,t1,t2
        real(dp),intent(in) ::smln 

        res2=0.0_dp

        ! res2=(7.0_dp/(4.0_dp*3.140_dp*(smln**2))) !Wendland c2
        ! res2=(9.0_dp/(4.0_dp*3.140_dp*(smln**2))) !Wendland c4

        ! res2=3/(2*3.14*smln**2) !Quadratic

        !Hybrid

        t1=(4.0_dp) !Cubic spline

        t2=(6.0_dp) !Hyperbolic

        res2=alp*t1*(10._dp/(28.0*3.14*smln**2)) + &
            (1.0_dp-alp)*t2*(1.0_dp/(3.0*3.14*smln**2))




    end function Wo

    function Wabx(p1,p2,radial,smln) result(res)
        implicit none
        type(particles),intent(in) :: p1,p2
        real(dp) :: res,t1,t2
        real(dp),intent(in) :: smln
        real(dp) ,intent(in):: radial

        res=0.0_dp
        ! res=(35.0_dp*(p1%x-p2%x)*((1.0_dp-radial/(2*smln)))**3)/ & !Wendland c2
        !     (4.0_dp*3.140_dp*smln**4)

        ! res=((p1%x-p2%x)*((1.0_dp-radial/(2*smln)))**5)* & !Wendland c4
        ! ((35.0_dp*radial/(3*(smln**3)))+(14.0_dp/(3*(smln**2))))* &
        ! (9.0_dp/(4.0_dp*3.140_dp*(smln**2)))

        ! res=(2/(3.14*smln**2))*((3/(8*smln**2))-(3/(4*smln*radial))) &
        !         *(p2%x-p1%x)                         !Quadratic


        !Hybrid Kernel 

        t1=merge(((9.0*radial/smln**3)-(12.0_dp/smln**2))*(p2%x-p1%x),&
        ((3.0/(radial*smln))*((2.0_dp-radial/smln)**2))*(p1%x-p2%x), &
        (radial<smln)) !Cubic spline

        t2=merge(((3.0*radial/smln**3)-6.0_dp/(radial*smln))*(p2%x-p1%x), &
        ((3.0/(radial*smln))*((2.0_dp-radial/smln)**2))*(p1%x-p2%x),& 
        (radial<smln)) !Hyperbolic

        res=alp*t1*(10._dp/(28.0*3.14*smln**2)) + &
            (1.0_dp-alp)*t2*(1.0_dp/(3.0*3.14*smln**2))


    end function Wabx

    function Waby(p1,p2,radial,smln) result(res)
        implicit none
        type(particles),intent(in) :: p1,p2
        real(dp) :: res,t1,t2
        real(dp),intent(in) :: smln
        real(dp) ,intent(in):: radial

        res=0.0_dp
        ! res=(35.0_dp*(p1%y-p2%y)*((1.0_dp-radial/(2*smln)))**3)/ & !Wendland
        !     (4.0_dp*3.140_dp*smln**4)

        ! res=((p1%y-p2%y)*((1.0_dp-radial/(2*smln)))**5)* & !Wendland c4
        ! ((35.0_dp*radial/(3*(smln**3)))+(14.0_dp/(3*(smln**2))))* &
        ! (9.0_dp/(4.0_dp*3.140_dp*(smln**2)))

        ! res=(2/(3.14*smln**2))*((3/(8*smln**2))-(3/(4*smln*radial))) &
        ! *(p2%y-p1%y)                         !Quadratic


        !Hybrid Kernel 

            t1=merge(((9.0*radial/smln**3)-(12.0_dp/smln**2))*(p2%y-p1%y),&
            ((3.0/(radial*smln))*((2.0_dp-radial/smln)**2))*(p1%y-p2%y), &
            (radial<smln)) !Cubic spline
    
            t2=merge(((3.0*radial/smln**3)-6.0_dp/(radial*smln))*(p2%y-p1%y), &
            ((3.0/(radial*smln))*((2.0_dp-radial/smln)**2))*(p1%y-p2%y),& 
            (radial<smln)) !Hyperbolic
    
            res=alp*t1*(10._dp/(28.0*3.14*smln**2)) + &
                (1.0_dp-alp)*t2*(1.0_dp/(3.0*3.14*smln**2))


    end function Waby

    function hWab(radial,smln) result(res2)
        implicit none
        real(dp) :: res2,t1,t2
        real(dp),intent(in) :: smln
        real(dp),intent(in) :: radial 

        res2=0.0_dp

        t1=merge((4.0_dp-(6.0*(radial**2)/smln**2)+ &
        (3.0*(radial**3)/smln**3)),(2.0_dp-radial/smln)**3, &
        (radial<smln)) !Cubic spline

        t2=merge(((radial/smln)**3-(6.0*radial/smln)+6.0_dp), &
        ((2.0_dp-radial/smln)**3),(radial<smln)) !Hyperbolic

        res2=alp*t1*(10._dp/(28.0*3.14*smln**2)) + &
            (1.0_dp-alp)*t2*(1.0_dp/(3.0*3.14*smln**2))

        

    end function hWab

    function hWo(smln) result(res2)
        implicit none
        real(dp) :: res2,t1,t2
        real(dp),intent(in) ::smln 

        res2=0.0_dp
        t1=(4.0_dp) !Cubic spline

        t2=(6.0_dp) !Hyperbolic

        res2=alp*t1*(10._dp/(28.0*3.14*smln**2)) + &
            (1.0_dp-alp)*t2*(1.0_dp/(3.0*3.14*smln**2))




    end function hWo

    function hWabx(p1,p2,radial,smln) result(res)
        implicit none
        type(particles),intent(in) :: p1,p2
        real(dp) :: res,t1,t2
        real(dp),intent(in) :: smln
        real(dp) ,intent(in):: radial

        res=0.0_dp

        t1=merge(((9.0*radial/smln**3)-(12.0_dp/smln**2))*(p2%x-p1%x),&
        ((3.0/(radial*smln))*((2.0_dp-radial/smln)**2))*(p1%x-p2%x), &
        (radial<smln)) !Cubic spline

        t2=merge(((3.0*radial/smln**3)-6.0_dp/(radial*smln))*(p2%x-p1%x), &
        ((3.0/(radial*smln))*((2.0_dp-radial/smln)**2))*(p1%x-p2%x),& 
        (radial<smln)) !Hyperbolic

        res=alp*t1*(10._dp/(28.0*3.14*smln**2)) + &
            (1.0_dp-alp)*t2*(1.0_dp/(3.0*3.14*smln**2))


    end function hWabx

    function hWaby(p1,p2,radial,smln) result(res)
        implicit none
        type(particles),intent(in) :: p1,p2
        real(dp) :: res,t1,t2
        real(dp),intent(in) :: smln
        real(dp) ,intent(in):: radial

        res=0.0_dp

        t1=merge(((9.0*radial/smln**3)-(12.0_dp/smln**2))*(p2%y-p1%y),&
        ((3.0/(radial*smln))*((2.0_dp-radial/smln)**2))*(p1%y-p2%y), &
        (radial<smln)) !Cubic spline

        t2=merge(((3.0*radial/smln**3)-6.0_dp/(radial*smln))*(p2%y-p1%y), &
        ((3.0/(radial*smln))*((2.0_dp-radial/smln)**2))*(p1%y-p2%y),& 
        (radial<smln)) !Hyperbolic

        res=alp*t1*(10._dp/(28.0*3.14*smln**2)) + &
            (1.0_dp-alp)*t2*(1.0_dp/(3.0*3.14*smln**2))


    end function hWaby

end module kernel