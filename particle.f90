!------------------------------------------------------------------------------
! IITKGP, PhD
!------------------------------------------------------------------------------
!
! MODULE: Particle
!
!> @author
!> Sinchan Roy Chowdhury}
!
! DESCRIPTION: 
!> Particle type definition.tid=1> fixed,2>ghost,3>internal
!
! REVISION HISTORY:
! 
! 
!------------------------------------------------------------------------------


module particle

    use, intrinsic :: iso_fortran_env, dp=>real64


    implicit none

    public

    ! integer, parameter,public :: dp = selected_real_kind(15) 


    type,public :: particles

        
        
        real(dp),public  :: mass=0.0_dp,density=0.0_dp,pressure=0.0_dp, &
                            vx=0.0_dp,vy=0.0_dp,vxs=0.0_dp,vys=0.0_dp,x=0.0_dp, &
                            y=0.0_dp,xs=0.0_dp,ys=0.0_dp,con=0.0_dp,ovol=0.0_dp,oden=0.0_dp
                            
        !Variables depending on dimension
        integer,public :: tid=0,pid=0,matid=0,bounlvl=0
        ! real(dp),allocatable,public :: meanstrn(:),tau(:),shift(:),bmls(:),coff(:),coffmg(:)
        logical,public :: mirror=.false.,xnorm=.false.,ynorm=.false.,vicinity=.false.,&
                            free=.false.,buffer=.false.,domain=.true.,range=.true.
        ! integer,allocatable :: wall(:)
        ! real(dp),allocatable,public :: posshift(:)

    end type particles

    type,public:: matrixsetup
        integer,public :: sz ! Number of elements
        real(dp),allocatable,dimension(:),public :: val !Storing the A matrix vales for a particle
        integer,allocatable,dimension(:),public :: col ! Corresponding column numbers
    end type matrixsetup

    type group
        integer :: order    ! original order of unsorted data
        integer :: value       ! values to be sorted by
    end type group

    type,public :: properties 
        real(dp),allocatable,public :: meanstrn(:),tau(:),coff(:),bmls(:)
        real(dp),public :: cdiff=0.0_dp,gradvx=0.0_dp,tden=0.0_dp,vel=0.0_dp,&
                            nut=0.0_dp,tke=0.0_dp,porosity=1.0_dp,resistx=0.0_dp, &
                            resisty=0.0_dp,varts=0.0_dp
                
        

        logical ::inpore=.false.
    end type properties

        type buffer

        real(dp) :: x=0.0_dp,y=0.0_dp
        integer :: id=0

    end type

    type reservoir

        integer,allocatable ::tank(:)
        integer:: si=0

    end type

    type exchange

        type(particles),pointer :: part

    end type

    type exchageprop

        type(properties),pointer :: ppart

    end type
        
    
    
end module particle