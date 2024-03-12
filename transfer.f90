!------------------------------------------------------------------------------
! IITKGP, PhD
!------------------------------------------------------------------------------
!
! MODULE: Transfer datatype for moving particles inbetween cells
!
!> @author
!> Sinchan Roy Chowdhury}
!
! DESCRIPTION: 
!> Particle type definition.tid=1> fixed,2>ghost,3>interpolation,4>internal
!
! REVISION HISTORY:
! 
! 
!------------------------------------------------------------------------------

module transfer

    use particle

    implicit none

    type trfvec

        type(particles),allocatable,dimension(:) :: list !Particles to be transferred
        integer :: count=0 !No of particles to be transferred


    end type trfvec
    
end module transfer