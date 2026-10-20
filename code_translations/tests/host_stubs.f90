!=======================================================================
! Minimal stand-ins for the host-code modules that the original
! al_nw_smearing_correction.f uses. Only for driver_original.f90.
!=======================================================================
module params
    implicit none
    integer, parameter :: idp = 8          ! Real kind
    integer :: MyProcNum = 1               ! Only processor 1 writes output
end module params

module ns_contr
    use params, only: idp
    implicit none
    real(kind=idp), dimension(1) :: timestep0 = 0.1d0   ! Time step per grid level
    integer :: grlev = 1                                ! Grid level
    real(kind=idp) :: density = 1.225d0                 ! Air density [kg/m^3]
    real(kind=idp) :: time = 0.d0                       ! Simulation time [s]
end module ns_contr

module TurbineData
    implicit none
    integer, parameter :: aclnamelen = 4
    character(len=aclnamelen) :: ACLproject = 'test'    ! Output file prefix
end module TurbineData

! CPU timer used for the timing column of the diagnostics files
function rtc()
    implicit none
    real*8 :: rtc
    call cpu_time(rtc)
end function rtc
