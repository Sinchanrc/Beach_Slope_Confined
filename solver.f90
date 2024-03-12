module solver

    use particle 
    implicit none
    
    contains

    function vecvec(l1s,a,b) result(res) ! Vector-vector multiplication
        implicit none 
        real(dp),intent(in) :: a(:),b(:)
        integer,intent(in) :: l1s
        integer :: ks
        real(dp) :: res
            ks=0
            res=0.0_dp

            !$omp parallel do simd reduction(+:res) schedule(static) default(shared)
            do ks=1,l1s
    
                res=real(res,dp)+real(a(ks),dp)*real(b(ks),dp)
    
            end do
            !$omp end parallel do simd

    end function

    function kahanvecvec(l1s,a,b) result(res) ! Vector-vector multiplication
        implicit none 
        real(dp),intent(in) :: a(:),b(:)
        integer,intent(in) :: l1s
        integer :: ks
        real(dp) :: res,c,y,t
            ks=0
            res=0.0_dp
            c=0.0_dp
            y=0.0_dp
            t=0.0_dp

            !$omp parallel do simd reduction(+:res,c) schedule(static) default(shared) private(y)
            do ks=1,l1s
    
                y=real(a(ks),dp)*real(b(ks),dp)-c
                t=res+y
                c=(t-res)-y
                res=t
    
            end do
            !$omp end parallel do simd
            res=res-c

    end function

    subroutine spmv(ls,val,row,col,b,res) ! Sparse Matrix-vector multiplication in CSR format
        implicit none
        integer,intent(in) :: ls
        integer,intent(in) :: row(ls+1),col(:)
        integer:: is,js
        real(dp),intent(in) :: val(:),b(ls)
        real(dp),intent(out) :: res(ls)
        real(dp) :: comp

        is=0
        js=0
        res(1:ls)=0.0_dp

        !$omp parallel do private(js,is,comp) default(shared) schedule(static)
            do is=1,ls
                comp=0.0_dp
            !$omp simd reduction(+:comp)
            do js=row(is),(row(is+1)-1)                
                if(col(js)/=0) then
                ! res(is)=real(res(is),dp)+real(val(js),dp)*real(b(col(js)),dp)
                comp=comp+real(val(js),dp)*real(b(col(js)),dp)
                end if
            end do
            !$omp end simd
            res(is)=comp
            end do
        !$omp end parallel do


    end subroutine

    subroutine spmv2(ls,b,res) ! Sparse Matrix-vector multiplication in CSR format
        use initialize
        implicit none
        integer,intent(in) :: ls
        integer:: is,js
        real(dp),intent(in) :: b(ls)
        real(dp),intent(out) :: res(ls)
        real(dp) :: comp,sol(ls)

        is=0
        js=0
        sol(1:ls)=0.0_dp
        
        !$omp parallel default(shared)

        !$omp do private(js,is) schedule(static)
            do is=1,ls
                ! comp=0.0_dp
                ! res(is)=0.0_dp
            ! !$omp simd reduction(+:comp)
            do js=1,fmatrix(is)%sz               

                sol(is)=sol(is)+fmatrix(is)%val(js)*b(fmatrix(is)%col(js))

            end do
            ! !$omp end simd
            ! res(is)=comp
            end do
        !$omp end do

        !$omp do private(is) schedule(static)
        do is=1,ls

            res(is)=sol(is)

        end do

        !$omp end do

        !$omp end parallel



    end subroutine

    subroutine bicgstab(tls,guess,ls,val,row,col,bvecs,sl) ! BiCGStab iterative solver
        implicit none
        integer,intent(in) :: ls
        real(dp),intent(in) :: val(:),bvecs(ls)   
        integer,intent(in) :: row(ls+1),col(:)
        real(dp) :: rs(ls),p(ls),res(ls),alpha,mm1(ls),&
                            mm2(ls),si(ls),wi,betas,prevr,reschk(ls)
        real(dp),intent(out) ::  sl(ls)
        real(dp),intent(inout) :: guess(ls)
        real(dp),intent(in) :: tls
        integer :: tcal=1,ks,n
        real(dp) :: error,gerror,itls


        itls=1e-2
    
        error=0.0_dp
        gerror=0.0_dp
        rs(1:ls)=0.0_dp
        p(1:ls)=0.0_dp
        alpha=0.0_dp
        mm1(1:ls)=0.0_dp
        mm2(1:ls)=0.0_dp
        si(1:ls)=0.0_dp
        wi=0.0_dp
        betas=0.0_dp
        prevr=0.0_dp
        sl(1:ls)=guess(1:ls)
        n=0
        tcal=1
        res(1:ls)=0.0_dp
        ks=0

        ! call spmv(ls,val,row,col,guess,res)
        call spmv2(ls,guess,res)
        !$omp parallel do simd default(shared) schedule(static)
            do ks=1,ls
            rs(ks)=real(bvecs(ks),dp)-real(res(ks),dp)
            p(ks)=real(rs(ks),dp)
            res(ks)=real(rs(ks),dp)
            end do
        !$omp end parallel do simd

        ! !$omp parallel do default(shared) reduction(+:gerror) schedule(static)
        !     do ks=1,ls 
        !         gerror=gerror+(rs(ls))**2
        !         ! gerror=gerror+bvecs(ls)**2
        !     end do
        ! !$omp end parallel do 

        ! gerror=sqrt(gerror)


        do while(tcal==1)
            
            n=n+1
            ! call spmv(ls,val,row,col,p,mm1)
            call spmv2(ls,p,mm1)
            prevr=real(vecvec(ls,res,rs),dp)
            alpha=(prevr/real(vecvec(ls,mm1,res),dp))

            !$omp parallel do simd default(shared) schedule(static)
            do ks=1,ls
                si(ks)=real(rs(ks),dp)-real(alpha*mm1(ks),dp)
            end do
            !$omp end parallel do simd

            ! call spmv(ls,val,row,col,si,mm2)
            call spmv2(ls,si,mm2)
            wi=(real(vecvec(ls,mm2,si),dp)/real(vecvec(ls,mm2,mm2),dp))

            !$omp parallel do simd default(shared) schedule(static)
            do ks=1,ls
                sl(ks)=real(sl(ks),dp)+real(alpha*p(ks),dp)+real(wi*si(ks),dp)
                rs(ks)=real(si(ks),dp)-real(wi*mm2(ks),dp)
            end do
            !$omp end parallel do simd
            
            betas=(real(vecvec(ls,rs,res),dp)/real(prevr,dp))*(real(alpha,dp)/real(wi,dp))

            !$omp parallel do simd default(shared) schedule(static)
            do ks=1,ls 
                p(ks)=real(rs(ks),dp)+real(betas*(p(ks)-wi*mm1(ks)),dp)
            end do
            !$omp end parallel do simd

            ! call spmv2(ls,sl,reschk)

            ! error=0.0_dp

            ! !$omp parallel do default(shared) reduction(+:error) schedule(runtime)
            ! do ks=1,ls 
            !     error=error+(bvecs(ls)-reschk(ls))**2
            ! end do
            ! !$omp end parallel do
            ! error=sqrt(error)

            ! if(real((error/gerror),dp)>itls) then
            !     tcal=1
            !     guess=sl
            ! else
            !     tcal=0
            ! end if

            do ks=1,ls

            if(abs(real((sl(ks)-guess(ks)),dp)/real((guess(ks)+1e-5),dp))>itls) then
                tcal=1
                guess=sl
                exit
            else
                tcal=0
            end if
            end do 

        end do

        !$omp parallel do simd default(shared) schedule(static)
        do ks=1,ls
            guess(ks)=sl(ks)
        end do
        !$omp end parallel do simd
    
        ! write(*,*)"Tol=",error/gerror
    
    end subroutine bicgstab

    subroutine format(array,val,row,col,len) ! Creating standard CSR format for solver
        use initialize
        implicit none
        type(matrixsetup),intent(in) :: array(:)
        real(dp),intent(out) :: val(:)
        integer,intent(out) :: row(:)
        integer,intent(out) ::col(:)
        integer,intent(in) :: len
        integer :: is,js

        nnz=0

        row(1)=1
        val(:)=0.0_dp
        do is=1,len
            row(is+1)=array(is)%sz+row(is)
            nnz=nnz+array(is)%sz
            do js=1,array(is)%sz
            val(row(is)-1+js)=array(is)%val(js)
            col(row(is)-1+js)=array(is)%col(js)
            end do
        end do

    end subroutine   


    subroutine pformat(array,val,row,col,len) ! Creating standard CSR format for solver
        use initialize
        implicit none
        type(matrixsetup),intent(in) :: array(:)
        real(dp),intent(out) :: val(:)
        integer,intent(out) :: row(:)
        integer,intent(out) ::col(:)
        integer,intent(in) :: len
        integer :: is,js

        row(1)=1
        val(:)=0.0_dp
        !!$omp parallel do default(shared) schedule(static) private(js,is)
        do is=1,len
            row(is+1)=ceiling(fac2*fplistmax)+(is-1)*ceiling(fac2*fplistmax)+1
            do js=1,ceiling(fac2*fplistmax)
            val(row(is)-1+js)=array(is)%val(js)
            col(row(is)-1+js)=array(is)%col(js)
            end do
        end do
        !!$omp end parallel do

    end subroutine     

    subroutine fgmres()

        use mkl_rci
        use initialize
        ! use blas95

        implicit none

        integer :: ipar(128),rci_request,gmresct=0 
        

        ! allocate(dpar(128),tmp(finmax*(2*150+1)+(150*(150+9))/2+1))

        call dfgmres_init(finmax,fsol,fvec,rci_request,ipar,dpar,tmp)

            ipar(9)=1
            ipar(10)=0
            ipar(12)=1
            dpar(1)=0.001



        11   call dfgmres(finmax,fsol,fvec,rci_request,ipar,dpar,tmp)

            if (rci_request==0) then

                call dfgmres_get(finmax,fsol,fvec,rci_request,ipar,dpar,tmp,gmresct)

            elseif (rci_request==1) then

                ! call mkl_dcsrgemv('N',finmax, fval,frow, fcol, tmp(ipar(22)), tmp(ipar(23)))
                ! call spmv (finmax,fval,frow,fcol,tmp(ipar(22)),tmp(ipar(23)))
                call spmv2(finmax,tmp(ipar(22)),tmp(ipar(23)))

                goto 11


            end if

            ! deallocate(dpar,tmp)
        
    end subroutine

    subroutine pardisosolver(fval1,frow1,fcol1,fvec1,fsol1)

        use omp_lib
        use mkl_pardiso
        use initialize

        implicit none


        real(dp),intent(inout) :: fval1(finmax*ceiling(fac2*fplistmax)),fvec1(finmax)
        real(dp),intent(inout) :: fsol1(finmax)
        integer,intent(in) :: frow1(finmax+1),fcol1(finmax*ceiling(fac2*fplistmax))

        real(dp),allocatable ::ar(:)
        integer :: i,perm1(finmax),idum(1)
        integer,allocatable ::ja(:),ia(:)
        real(dp) :: ddum(1)

        call pardisoinit (pt, mtype, iparm)

        ! pt(:)%DUMMY=0
        ! iparm=0

        iparm(1) = 1 ! no solver default
        iparm(2) = 0 ! fill-in reordering from METIS
        iparm(3) = OMP_GET_NUM_THREADS() ! numbers of processors
        ! iparm(4) = 31 ! no iterative-direct algorithm
        iparm(5) = 0 ! no user fill-in reducing permutation
        iparm(6) = 0 ! =0 solution on the first n components of x
        iparm(7) = 0 ! not in use
        iparm(8) = 2 ! numbers of iterative refinement steps
        iparm(9) = 0 ! not in use
        iparm(10) = 13 ! perturb the pivot elements with 1E-13
        iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
        iparm(12) = 0 ! not in use
        iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
        iparm(25)=0
        iparm(27)=1

        ! iparm(:)=0


        ! iparm(1)=1
        ! iparm(2)=2
        ! iparm(3)=OMP_GET_MAX_THREADS()
        ! iparm(8)=4
        ! iparm(10)=13
        ! iparm(11)=1
        ! iparm(13)=1
        ! ! iparm(15)=45000
        ! ! iparm(16)=45000
        ! ! iparm(17)=90000
        ! iparm(25)=1
        ! iparm(27)=1

        allocate(ar(nnz),ja(nnz),ia(finmax+1))

        perm1(:)=0

        do i=1,nnz 

            ar(i)=fval1(i)
            ja(i)=fcol1(i)

        end do

        do i=1,(finmax+1)

            ia(i)=frow1(i)

        end do

        fsol1=0.0_dp

        ! phase = 11

        ! call pardiso(pt,maxfct,mnum,mtype,phase,finmax,ar,ia,ja,idum,1,iparm,msglvl,ddum,ddum,err1)

        ! phase = 22

        ! call pardiso(pt,maxfct,mnum,mtype,phase,finmax,ar,ia,ja,idum,1,iparm,msglvl,ddum,ddum,err1)

        phase = 13

        call pardiso(pt,maxfct,mnum,mtype,phase,finmax,ar,ia,ja,perm1,1,iparm,msglvl,fvec1,fsol1,err1)

        phase=-1

        call pardiso (pt, maxfct, mnum, mtype, phase, finmax, ddum, idum,idum, idum, nrhs, iparm, msglvl, ddum, ddum, err1)

        ! deallocate(ar,ja,ia)

        

    end subroutine

    
end module solver