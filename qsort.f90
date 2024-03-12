module qsort1
use initialize
use particle
implicit none


!call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)


contains

recursive subroutine QSort(A,nA)

    ! DUMMY ARGUMENTS
    integer, intent(in) :: nA
    type (group), dimension(nA), intent(in out) :: A

    ! LOCAL VARIABLES
    integer :: left, right
    real :: random
    integer :: pivot
    type (group) :: temp
    integer :: marker

    if (nA > 1) then

        call random_number(random)
        pivot = A(int(random*real(nA-1))+1)%value   ! random pivor (not best performance, but avoids worst-case)
        left = 0
        right = nA + 1

        do while (left < right)
            right = right - 1
            do while (A(right)%value > pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(left)%value < pivot)
                left = left + 1
            end do
            if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
            end if
        end do

        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
        call QSort(A(:marker-1),marker-1)
        call QSort(A(marker:),nA-marker+1)

    end if
end subroutine QSort


subroutine pformat(array,val,row,col,len) 
    use initialize
    use particle
    implicit none
    type(matrixsetup),intent(inout) :: array(:)
    double precision,intent(out) :: val(:)
    integer,intent(out) :: row(:)
    integer,intent(out) ::col(:)
    integer,intent(in) :: len
    integer :: is,js,rancol,ks,tcal=0
    type(group) :: Alist(len,fplistmax)
    real :: random,ranint

        do is=1,len
        do js =array(is)%sz+1,fplistmax
            do while(tcal==0)
            if (array(is)%col(array(is)%sz)<=finmax/2) then
                call random_number(ranint)
                rancol=finmax/2 + FLOOR(((finmax)+1-finmax/2)*ranint)
                do ks=1,js
                    if(array(is)%col(ks)==rancol) then
                    tcal=0
                    exit
                    else
                    tcal=1
                    array(is)%col(js)=rancol
                    end if
                end do
            else
                call random_number(ranint)
                rancol=1 + FLOOR(((finmax/2))*ranint)
                do ks=1,js
                    if(array(is)%col(ks)==rancol) then
                    tcal=0
                    exit
                    else
                    tcal=1
                    array(is)%col(js)=rancol
                    end if
                end do
            end if
            end do
        end do
        end do

        do is=1,len
        do js =1,fplistmax
            Alist(is,js)%order=js
            Alist(is,js)%value=array(is)%col(js)
        end do
        call QSort(Alist(is,:),fplistmax)
        end do



        row(1)=1
        val=0.0d0
        do is=1,len
            row(is+1)=fplistmax+row(is)
            do js=1,fplistmax
            val(row(is)-1+js)=array(is)%val(Alist(is,js)%order)
            col(row(is)-1+js)=array(is)%col(Alist(is,js)%order)
            end do
        end do
    
end subroutine pformat

end module qsort1