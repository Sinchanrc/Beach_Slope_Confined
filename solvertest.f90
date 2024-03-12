program test
      use particle
      use solver

      implicit none

! #include "lisf.h"
        ! integer, parameter :: dp = selected_real_kind(15, 100) 

        real(dp),dimension(5) ::b=(/8.0_dp,6.0_dp,4.0_dp,9.0_dp,8.0_dp/),guess,sol
        real(dp),dimension(13) :: val
        integer,dimension(6) :: row=(/1,4,6,8,11,14/)
        integer,dimension(13)  ::  col=(/1,3,4,2,4,1,3,2,3,4,1,2,5/)
        real(dp) :: tls=0.010_dp,zer=0.0d0
    
        val=1.0_dp
        guess=1.0_dp

        call bicgstab(tls,guess,5,val,row,col,b,sol)
        write(*,*)sol!,4*atan(1.0_dp)
!         LIS_MATRIX :: A
!         LIS_SCALAR,allocatable :: b(:),x(:)
!         LIS_SOLVER :: solver
!         LIS_INTEGER,allocatable :: ptr(:),index1(:)
!         LIS_SCALAR,allocatable :: value(:)
!         LIS_INTEGER :: ierr,nnz=13,n=5
!         integer ::i
!         real,allocatable :: b1(:),ptr1(:),index2(:)
!         call lis_initialize(ierr)

!         allocate(ptr(0:n),index1(0:12),value(0:12),b1(5),b(0:4),x(0:4),&
!                   index2(13),ptr1(6))

!         index2=(/0,2,3,1,3,0,2,1,2,3,0,1,4/)
!         value=1.0
!         ptr1=(/0,3,5,7,10,13/)
!         b1=(/8,6,4,9,8/)
!         do i=0,12
!             index1(i)=index2(i+1)
!         end do
!         do i=0,n
!             ptr(i)=ptr1(i+1)
!         end do
!         do i=0,n-1
!             b(i)=b1(i+1)
!         end do
!         call lis_matrix_create(0,A,ierr)
!         call CHKERR(ierr)
!         call lis_matrix_set_size(A,0,n,ierr)
!         call CHKERR(ierr)

!         call lis_matrix_set_csr(nnz,ptr,index1,value,A,ierr)
!         call CHKERR(ierr)
!         call lis_matrix_assemble(A,ierr)
!       !   call CHKERR(ierr)
!       !   call lis_vector_create(0,b,ierr)
!       !   call CHKERR(ierr)
!       !   call lis_vector_set_size(b,0,n,ierr)
!       !   call lis_vector_create(0,x,ierr)
!       !   call CHKERR(ierr)
!       !   call lis_vector_set_size(x,0,n,ierr)
!       !   do i=1,n
!       !       CALL lis_vector_set_value(LIS_INS_VALUE,i,b1(i),b,ierr)
!       !   end do
!         call lis_solver_create(solver,ierr)
!         call CHKERR(ierr)
!         call lis_solver_set_option('-i bicg -p jacobi',solver,ierr)
!         call CHKERR(ierr)
!         call lis_solver_set_option('-tol 1.0e-12',solver,ierr)
!         call CHKERR(ierr)
!         call lis_solve(A,b,x,solver,ierr)
!         call CHKERR(ierr)
!       ! !   call lis_array_print(x)
!       !   write(*,*)x


        
!         call lis_finalize(ierr)

!       stop
      end program
