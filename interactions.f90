!------------------------------------------------------------------------------
! IITKGP, PhD
!------------------------------------------------------------------------------
!
! MODULE: Particle interactions
!
!> @author
!> Sinchan Roy Chowdhury}
!
! DESCRIPTION: 
!> Datatype stores the cell id and the location of interacting
!  particle in the particle list of that cell.
! 
!
! REVISION HISTORY:
! 
! 
!------------------------------------------------------------------------------

module interactions

    use particle

    implicit none

    type verlet

        type(exchange),allocatable :: nh(:)
        type(exchageprop),allocatable :: pnh(:)
        real(dp),allocatable :: dist(:),klt(:,:) ! Radial distance of interaction particle
        integer,allocatable :: interlist(:,:) ! index of cell and location number in list of particle
        integer :: count=0 ! Total 
        logical :: error=.false.

    end type verlet

    type pprob
        real(dp) :: pressure=0._dp,pos(2)=0._dp
        integer :: ci,cj
        type(verlet) :: part
    end type pprob

    
end module interactions