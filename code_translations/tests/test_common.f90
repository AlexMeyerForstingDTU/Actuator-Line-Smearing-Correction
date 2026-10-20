!=======================================================================
! Test cases shared by driver_original.f90 and driver_fortran.f90.
! driver_cpp.cpp and driver_python.py repeat the same definitions.
!
! Test blade (Note_review_of_PhD_comments.md, Sect. 5): NREL 5MW-like,
! 64 uniform trailing points, cylinder polars inboard of 10 m and a
! thin-airfoil-like polar outboard, 8 m/s at 9.2 rpm, frozen "CFD"
! inflow with induction a = 0.3.
!=======================================================================
module test_common
    implicit none

    integer, parameter :: ns = 63, nv = ns + 1, len2d = 181, ntab = 17
    real*8, parameter :: pi = 4.d0*atan(1.d0), d2r = pi/180.d0, r2d = 180.d0/pi
    real*8, parameter :: wind = 8.d0, a_ind = 0.3d0, rho = 1.225d0
    real*8, parameter :: omega = 9.2d0*2.d0*pi/60.d0

    ! Chord and twist of the blade at the table radii
    real*8, parameter, dimension(ntab) :: rtab = [2.8667d0, 5.6d0, 8.3333d0, 11.75d0, &
        15.85d0, 19.95d0, 24.05d0, 28.15d0, 32.25d0, 36.35d0, 40.45d0, 44.55d0, 48.65d0, &
        52.75d0, 56.1667d0, 58.9d0, 61.6333d0]
    real*8, parameter, dimension(ntab) :: ttab = [13.308d0, 13.308d0, 13.308d0, 13.308d0, &
        11.480d0, 10.162d0, 9.011d0, 7.795d0, 6.544d0, 5.361d0, 4.188d0, 3.125d0, 2.319d0, &
        1.526d0, 0.863d0, 0.370d0, 0.106d0]
    real*8, parameter, dimension(ntab) :: ctab = [3.542d0, 3.854d0, 4.167d0, 4.557d0, &
        4.652d0, 4.458d0, 4.249d0, 4.007d0, 3.748d0, 3.502d0, 3.256d0, 3.010d0, 2.764d0, &
        2.518d0, 2.313d0, 2.086d0, 1.419d0]

    ! One test case: frozen inflow marched in time
    type test_case
        real*8 :: eps = 10.d0              ! Smearing length [m]
        real*8 :: cone = 0.d0              ! Cone angle [deg]
        real*8 :: dt = 0.1d0               ! Time step [s]
        integer :: n_steps = 300
        logical :: l1point = .false.       ! Ang Li correction
        logical :: bent = .false.          ! Pre-bent, swept blade with out-of-plane offsets
        logical :: root_at_zero = .false.  ! First trailing point at r = 0
        real*8 :: max_resid = 1.d-6        ! Convergence limit [m/s]
        integer :: pitch_step_at = 0       ! From this step on: pitch = 2 deg (0 = never)
        integer :: n_rotors = 1, n_blades = 1
        real*8 :: pitch_22 = 0.d0          ! Pitch of rotor 2, blade 2 [deg]
    end type test_case

    ! Geometry, polars and inflow of the test blade
    type test_blade
        real*8, dimension(nv) :: rv
        real*8, dimension(ns) :: rc, chord, twist, fbend, sweep, dsweep
        real*8, dimension(ns, 3) :: cyl_cp, vcyl
        real*8, dimension(nv, 3) :: cyl_vtp
        real*8, dimension(len2d, ns) :: alpha_2d, cl_2d, cd_2d
    end type test_blade

contains

    ! Parameters of a named test case
    function get_case(name) result(c)
        implicit none

        ! --- Input  ---
        character(len=*), intent(in) :: name

        ! --- Output ---
        type(test_case) :: c

        select case (trim(name))
        case ('eps10_dt001')
            c%dt = 0.01d0
            c%n_steps = 3000
        case ('eps10')
        case ('cone25')
            c%cone = 2.5d0
        case ('eps25')
            c%eps = 2.5d0
        case ('nocorr')
            c%eps = 0.01d0
        case ('bent')
            c%cone = 2.5d0
            c%bent = .true.
            c%l1point = .true.
        case ('pitchstep')
            c%n_steps = 400
            c%pitch_step_at = 301
        case ('noconv')
            c%n_steps = 2
            c%max_resid = 0.d0
        case ('multi')
            c%n_steps = 100
            c%n_rotors = 2
            c%n_blades = 2
            c%pitch_22 = 1.d0
        case ('rootzero')
            c%n_steps = 1
            c%root_at_zero = .true.
        case default
            stop 'Error: unknown test case.'
        end select
    end function get_case

    ! Piecewise linear table lookup, constant beyond the ends
    function table_interp(x, xt, yt) result(y)
        implicit none

        ! --- Input  ---
        real*8, intent(in) :: x
        real*8, intent(in), dimension(:) :: xt, yt

        ! --- Output ---
        real*8 :: y

        ! --- Dummy variables ---
        integer :: k, n

        n = size(xt)
        if (x <= xt(1)) then
            y = yt(1)
        else if (x >= xt(n)) then
            y = yt(n)
        else
            k = 2
            do while (x > xt(k))
                k = k + 1
            end do
            y = yt(k-1) + (x - xt(k-1))/(xt(k) - xt(k-1))*(yt(k) - yt(k-1))
        end if
    end function table_interp

    ! Geometry, polars and frozen cylindrical inflow of the test blade
    subroutine build_blade(c, b)
        implicit none

        ! --- Input  ---
        type(test_case), intent(in) :: c

        ! --- Output ---
        type(test_blade), intent(out) :: b

        ! --- Dummy variables ---
        integer :: k, l, m
        real*8 :: r0, t

        ! 1. Trailing points (uniform) and control points in between
        r0 = 1.5d0
        if (c%root_at_zero) r0 = 0.d0
        do k = 1, nv
            b%rv(k) = r0 + (63.d0 - r0)*dble(k-1)/dble(nv-1)
        end do
        b%rc = 0.5d0*(b%rv(1:ns) + b%rv(2:nv))
        do l = 1, ns
            b%chord(l) = table_interp(b%rc(l), rtab, ctab)
            b%twist(l) = table_interp(b%rc(l), rtab, ttab)
        end do

        ! 2. Blade shape: straight, or pre-bent and swept with out-of-plane offsets
        b%fbend = 0.d0
        b%sweep = 0.d0
        b%dsweep = 0.d0
        b%cyl_cp = 0.d0
        b%cyl_vtp = 0.d0
        b%cyl_cp(:, 1) = b%rc
        b%cyl_vtp(:, 1) = b%rv
        if (c%bent) then
            do l = 1, ns
                t = b%rc(l)/63.d0
                b%fbend(l) = 4.d0*t*t
                b%sweep(l) = 3.d0*t*t*t
                b%dsweep(l) = 0.5d0*b%sweep(l)
                b%cyl_cp(l, 3) = -1.5d0*t*t
            end do
            do k = 1, nv
                t = b%rv(k)/63.d0
                b%cyl_vtp(k, 3) = -1.5d0*t*t
            end do
        end if

        ! 3. Polars: cylinder inboard of 10 m, smooth thin-airfoil-like outboard
        do l = 1, ns
            do m = 1, len2d
                b%alpha_2d(m, l) = -90.d0 + dble(m-1)
                if (b%rc(l) < 10.d0) then
                    b%cl_2d(m, l) = 0.d0
                    b%cd_2d(m, l) = 0.4d0
                else
                    b%cl_2d(m, l) = 1.4d0*tanh(2.d0*pi*(b%alpha_2d(m, l) + 3.d0)*d2r/1.4d0)
                    b%cd_2d(m, l) = 0.01d0
                end if
            end do
        end do

        ! 4. Frozen "CFD" inflow in the cylindrical frame: uniform a = 0.3
        b%vcyl = 0.d0
        b%vcyl(:, 3) = wind*(1.d0 - a_ind)
    end subroutine build_blade

    ! Rotor power and thrust from the section forces of one blade
    subroutine power_thrust(b, Fstn, factor, P, T)
        implicit none

        ! --- Input  ---
        type(test_blade), intent(in) :: b
        real*8, intent(in), dimension(ns, 3) :: Fstn
        real*8, intent(in) :: factor          ! Number of blades the forces stand for

        ! --- Output ---
        real*8, intent(out) :: P, T

        ! --- Dummy variables ---
        integer :: l

        P = 0.d0
        T = 0.d0
        do l = 1, ns
            P = P + factor*Fstn(l, 2)*omega*b%rc(l)*(b%rv(l+1) - b%rv(l))
            T = T + factor*Fstn(l, 3)*(b%rv(l+1) - b%rv(l))
        end do
    end subroutine power_thrust

    ! Write one named result line: name v1 v2 ...
    subroutine write_values(unit, name, values)
        implicit none

        ! --- Input  ---
        integer, intent(in) :: unit
        character(len=*), intent(in) :: name
        real*8, intent(in), dimension(:) :: values

        write(unit, '(a, *(1x, es24.16e3))') name, values
    end subroutine write_values

end module test_common
