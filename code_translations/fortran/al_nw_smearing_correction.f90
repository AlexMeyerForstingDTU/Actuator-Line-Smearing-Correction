!=======================================================================
! Fast smearing correction for the actuator line (free-form Fortran)
!
! Translated from al_nw_smearing_correction.f. Method: Meyer Forsting
! et al. (2019, 2020), with the near-wake model of Pirrung et al. (2016,
! 2017a, 2017b) and the 2020 corrigenda. Equation numbers refer to these
! papers.
!
! Usage:
!   call init_nw_smearing_correction(n_rotors, n_blades, n_sections, project)
!   then every time step, for each rotor i and blade j:
!   call nw_smearing_correction(i, j, time, dt, density, ..., Fstn)
!=======================================================================
module al_nw_smearing_correction
    implicit none

    real*8, parameter :: pi = 4.d0*atan(1.d0)
    real*8, parameter :: d2r = pi/180.d0, r2d = 180.d0/pi

    ! Fit constants of the mixing factor k_Phi, Pirrung et al. (2017b).
    ! Row i gives a_hr,i. P: h/r > 0, Eq. (14). N: h/r < 0, Eq. (13).
    real*8, parameter, dimension(4, 5) :: nw_P = reshape([ &
        -1.64636754184988d0, 8.14821474772595d0, -12.1784861715161d0, 5.02652654794305d0, 21.7713111885100d0, &
        -0.499006672527007d0, 6.08465268675599d0, -15.1712005446736d0, 14.8254135311781d0, -2.42318681495122d0, &
        3.90835606394205d0, -18.7662330478906d0, 39.1243282362892d0, -29.4870086106427d0, -60.2947309276912d0, &
        -1.60622843516794d0, 7.42952711674470d0, -15.8594777074799d0, 11.6870163960092d0, -195.060865844227d0], &
        [4, 5], order=[2, 1])
    real*8, parameter, dimension(4, 5) :: nw_N = reshape([ &
        1.01933206359188d0, -0.135668151345576d0, 0.395524705910666d0, 0.0801822997046969d0, 44.8347503471217d0, &
        12.9874512870379d0, 49.9999999999916d0, 0.00235345218088083d0, 11.3116074964299d0, 3935.34323353354d0, &
        -0.690159136848837d0, 101.238775997879d0, -0.00154467589177718d0, 3.99520376143011d0, 0.394541495293350d0, &
        -0.269253409815015d0, 49.9999999999843d0, -0.00247975579295973d0, 0.403642518874381d0, 1.16610462658971d0], &
        [4, 5], order=[2, 1])

    ! Settings
    logical :: kaya = .true.            ! Switch: correction on/off
    real*8 :: nw_max_resid = 1.d-6      ! Convergence limit of the correction velocity [m/s]
    real*8 :: x_aoa = 0.5d0             ! Chordwise position of the AoA in the CFD [-]

    ! Number of control points (CP, sections) and trailing points per blade
    integer :: n_s = 0, n_v = 0

    ! History and geometry constants of one rotor blade
    type blade_state
        logical :: refresh = .true.                     ! Recompute the geometry constants
        ! First and last CP felt by each vortex (n_v, 2)
        integer, allocatable, dimension(:, :) :: cp_loop_ind
        ! History of the last time step: (n_s), (n_s, 3) and (n_v, n_s)
        real*8, allocatable, dimension(:) :: W_history, vn_history, vt_history, Gamma_history
        real*8, allocatable, dimension(:, :) :: dyn_comp_history
        real*8, allocatable, dimension(:, :) :: X_history, Y_history
        ! Geometry constants (n_v, n_s): root-corrected Phi_C, Phi of a straight vortex,
        ! amplitudes, fit coefficients a_hr,i (n_v, n_s, 4) and smearing factors
        real*8, allocatable, dimension(:, :) :: phi, phi_s, dx, dy
        real*8, allocatable, dimension(:, :, :) :: a
        real*8, allocatable, dimension(:, :) :: smear_fac
        integer :: unit = 0                             ! Diagnostics file, 0 = none
    end type blade_state

    type(blade_state), allocatable, dimension(:, :) :: blade        ! (n_rotors, n_blades)

contains

    !-------------------------------------------------------------------
    ! Allocate the state of all rotor blades and open the output files
    ! <project>.nwRR.blBB ('' = no output)
    !-------------------------------------------------------------------
    subroutine init_nw_smearing_correction(n_rotors, n_blades, n_sections, project)
        implicit none

        ! --- Input  ---
        integer, intent(in) :: n_rotors, n_blades, n_sections
        character(len=*), intent(in) :: project

        ! --- Dummy variables ---
        integer :: i, j
        character(len=16) :: suffix

        n_s = n_sections
        n_v = n_sections + 1
        allocate(blade(n_rotors, n_blades))
        do i = 1, n_rotors
            do j = 1, n_blades
                allocate(blade(i, j)%cp_loop_ind(n_v, 2))
                blade(i, j)%cp_loop_ind(:, 1) = 1
                blade(i, j)%cp_loop_ind(:, 2) = n_s
                allocate(blade(i, j)%W_history(n_s), blade(i, j)%vn_history(n_s))
                allocate(blade(i, j)%vt_history(n_s), blade(i, j)%Gamma_history(n_s))
                allocate(blade(i, j)%dyn_comp_history(n_s, 3))
                allocate(blade(i, j)%X_history(n_v, n_s), blade(i, j)%Y_history(n_v, n_s))
                allocate(blade(i, j)%phi(n_v, n_s), blade(i, j)%phi_s(n_v, n_s))
                allocate(blade(i, j)%dx(n_v, n_s), blade(i, j)%dy(n_v, n_s))
                allocate(blade(i, j)%a(n_v, n_s, 4), blade(i, j)%smear_fac(n_v, n_s))
                blade(i, j)%W_history = 0.d0
                blade(i, j)%vn_history = 0.d0
                blade(i, j)%vt_history = 0.d0
                blade(i, j)%Gamma_history = 0.d0
                blade(i, j)%dyn_comp_history = 0.d0
                blade(i, j)%X_history = 0.d0
                blade(i, j)%Y_history = 0.d0
                if (len_trim(project) > 0) then
                    write(suffix, '(a, i2.2, a, i2.2)') '.nw', i, '.bl', j
                    open(newunit=blade(i, j)%unit, file=trim(project)//trim(suffix), &
                         status='replace')
                end if
            end do
        end do
    end subroutine init_nw_smearing_correction

    !-------------------------------------------------------------------
    ! Corrected blade forces of rotor i, blade j for one time step.
    ! The smeared AL forces leave part of the induction out of the CFD
    ! velocities; the near-wake model supplies that part and corrects
    ! the angle of attack. Positions (r, theta, z) are in the cylindrical
    ! rotor frame; velocities in section axes (s, t, n) or in the
    ! cylindrical frame.
    !-------------------------------------------------------------------
    subroutine nw_smearing_correction(i, j, time, dt, density, cyl_cp, cyl_vtp, vstn_cfd, &
                                      vcyl_cfd, eps_cp, omega, chord, twist, pitch, cl_2d, cd_2d, &
                                      alpha_2d, cone, fbend, sweep, dsweep, l1point, Fstn)
        implicit none

        ! --- Input  ---
        integer, intent(in) :: i, j                         ! Rotor and blade number
        real*8, intent(in) :: time, dt, density             ! [s], [s], [kg/m^3]
        real*8, intent(in), dimension(n_s, 3) :: cyl_cp     ! CP positions
        real*8, intent(in), dimension(n_v, 3) :: cyl_vtp    ! Trailing-point positions
        real*8, intent(in), dimension(n_s, 3) :: vstn_cfd   ! CFD velocity at the CPs (s, t, n)
        real*8, intent(in), dimension(n_s, 3) :: vcyl_cfd   ! CFD velocity at the CPs (cylindrical)
        real*8, intent(in) :: eps_cp                        ! Smearing length [m]
        real*8, intent(in) :: omega                         ! Rotor speed [rad/s]
        real*8, intent(in), dimension(n_s) :: chord, twist  ! Chord [m] and twist [deg] at the CPs
        real*8, intent(in) :: pitch                         ! Blade pitch [deg] (+ve to feather)
        real*8, intent(in), dimension(:, :) :: cl_2d, cd_2d, alpha_2d   ! Polars (len_2d, n_s) [deg]
        real*8, intent(in) :: cone                          ! Cone [deg] (+ve upstream)
        real*8, intent(in), dimension(n_s) :: fbend, sweep, dsweep  ! Bend, sweep, delta sweep [deg]
        logical, intent(in) :: l1point                      ! Ang Li correction on/off

        ! --- Output ---
        real*8, intent(out), dimension(n_s, 3) :: Fstn      ! Forces per unit length (s, t, n) [N/m]

        ! --- Dummy variables ---
        logical :: recompute_relax, converged
        integer :: sub_iter, k, l
        real*8 :: h, h_r, delta, core_cutoff, W_k
        real*8 :: relax, relax_safety, resmax, cpu_t1, cpu_t2
        real*8, dimension(n_s) :: r_cp, W, vn, vt, W_iter, Wn_iter, Wt_iter, section_res
        real*8, dimension(n_s) :: scone, v_rot, vn_rot, vt_rot, un, ut, vrel_cp
        real*8, dimension(n_s) :: phi_cp, phi_cp_cfd, phi_cp_cyl, phi_cp_14, thetadot
        real*8, dimension(n_s) :: alpha_cp, alpha_cp_cfd, alpha_cp_14, alpha_cp_34
        real*8, dimension(n_s) :: Gamma_steady, Gamma, cl, cd, lift, drag
        real*8, dimension(n_v) :: r_vtp, phi_vtp, vrel_vtp, dGamma, dbeta
        real*8, dimension(n_s, 3) :: vstn, dyn_comp
        real*8, dimension(n_v, n_s) :: phi_star, X, Y

        if (.not. kaya) return
        r_cp = cyl_cp(:, 1)
        r_vtp = cyl_vtp(:, 1)
        ! Important: r = 0 at the first trailing point is a singularity
        if (r_vtp(1) < 0.1d0) then
            write(*, '(a)') '!!!! WARNING smearing correction switched off !!!!'
            write(*, '(a)') 'First vortex trailing point is zero, please modify'
            kaya = .false.
            return
        end if
        call cpu_time(cpu_t1)

        ! b is short for blade(i, j)
        associate(b => blade(i, j))

        ! 1. Preamble: start from the previous time step
        relax_safety = 0.2d0
        recompute_relax = .true.
        W = b%W_history
        vn = b%vn_history
        vt = b%vt_history
        scone = sin((cone + fbend)*d2r)
        v_rot = -(omega*r_cp)                       ! Rotational speed in section axes
        vn_rot = -v_rot*sin(-pitch*d2r)
        vt_rot = v_rot*cos(-pitch*d2r)*cos(dsweep*d2r)

        ! 2. Geometry constants: first call, or after a refresh
        if (b%refresh) then
            core_cutoff = 1.82d0*eps_cp             ! 99th percentile of the Gaussian
            b%cp_loop_ind(:, 1) = 1
            b%cp_loop_ind(:, 2) = n_s
            do k = 1, n_v
                do l = 1, n_s
                    h = r_vtp(k) - r_cp(l)
                    h_r = h/r_vtp(k)
                    call nw_hr_components(h, h_r, b%phi(k, l), b%phi_s(k, l), b%dx(k, l), &
                                          b%dy(k, l), b%a(k, l, :))
                    b%smear_fac(k, l) = smearing_factor(r_vtp(k), h_r, 0.d0, 0.d0, eps_cp)
                    ! CPs felt by vortex k: inboard and outboard ones within the cut-off
                    if (h > 0.d0) then
                        if (abs(h) > core_cutoff) b%cp_loop_ind(k, 1) = l + 1
                    else
                        if (abs(h) < core_cutoff) b%cp_loop_ind(k, 2) = l
                    end if
                end do
            end do
            b%refresh = .false.
        end if

        ! 3. Helix angles from the CFD only (MF2019 step 2), fixed during the sub-iteration
        phi_cp_cfd = atan2(vcyl_cfd(:, 3), -(v_rot + vcyl_cfd(:, 2)))
        do l = 1, n_s
            if (phi_cp_cfd(l) > pi/2.d0) phi_cp_cfd(l) = pi - phi_cp_cfd(l)
        end do
        phi_vtp(1) = phi_cp_cfd(1)
        phi_vtp(n_v) = phi_cp_cfd(n_s)
        do k = 2, n_v - 1
            phi_vtp(k) = nw_interp1d(r_cp, phi_cp_cfd, r_vtp(k))
        end do
        do k = 1, n_v
            do l = b%cp_loop_ind(k, 1), b%cp_loop_ind(k, 2)
                h = r_vtp(k) - r_cp(l)
                h_r = h/r_vtp(k)
                phi_star(k, l) = nw_phi_components(h_r, b%a(k, l, :), phi_vtp(k), b%phi(k, l), &
                                                   b%phi_s(k, l))
                delta = -(cyl_vtp(k, 3) - cyl_cp(l, 3))     ! Normal to the disc, +ve upstream
                b%smear_fac(k, l) = smearing_factor(r_vtp(k), h_r, delta, phi_vtp(k), eps_cp)
            end do
        end do

        ! 4. Sub-iterate until the correction velocities converge
        sub_iter = 0
        converged = .false.
        do while (.not. converged)
            ! Velocities and angles with the current correction
            vstn = nw2vstn(vt, vn, cone, pitch, sweep, fbend)
            un = vn_rot + vstn_cfd(:, 3) + vstn(:, 3)        ! Normal velocity
            ut = vt_rot + vstn_cfd(:, 2) + vstn(:, 2)        ! Tangential velocity (< 0)
            phi_cp = atan2(un, -ut)
            alpha_cp = phi_cp*r2d - twist
            vrel_cp = sqrt(un**2 + ut**2)
            vrel_vtp(1) = vrel_cp(1)
            vrel_vtp(n_v) = vrel_cp(n_s)
            do k = 2, n_v - 1
                vrel_vtp(k) = nw_interp1d(r_cp, vrel_cp, r_vtp(k))
            end do
            call nw_circulation(dt, vrel_cp, alpha_cp, b%Gamma_history, b%dyn_comp_history, &
                                alpha_2d, cl_2d, chord, Gamma_steady, Gamma, dGamma, dyn_comp)
            dbeta = (vrel_vtp*dt)/r_vtp                      ! Advance in beta*, 2017b Eq. (3)

            ! Induction of the vorticity shed this time step, smeared part only
            W_iter = 0.d0
            Wn_iter = 0.d0
            Wt_iter = 0.d0
            X = 0.d0
            Y = 0.d0
            do k = 1, n_v
                do l = b%cp_loop_ind(k, 1), b%cp_loop_ind(k, 2)
                    call nw_xy_components(phi_star(k, l), b%dx(k, l), b%dy(k, l), dbeta(k), &
                                          dGamma(k), b%X_history(k, l), b%Y_history(k, l), &
                                          X(k, l), Y(k, l))
                    W_k = b%smear_fac(k, l)*(X(k, l) + Y(k, l))
                    W_iter(l) = W_iter(l) + W_k
                    ! Normal and tangential parts with each vortex's helix angle (MF2019 Eq. 19)
                    Wn_iter(l) = Wn_iter(l) + W_k*cos(phi_vtp(k))
                    Wt_iter(l) = Wt_iter(l) + W_k*sin(phi_vtp(k))
                end do
            end do

            if (recompute_relax) then
                relax = nw_calc_relax(vrel_cp, alpha_cp, un, -ut, chord, dt, b%dx, b%dy, dbeta, &
                                      b%phi, relax_safety)
                recompute_relax = .false.
            end if

            ! Relax and measure the change
            W_iter = relax*W + (1.d0 - relax)*W_iter
            section_res = sqrt((W_iter - W)**2)
            phi_cp_cyl = atan2(vcyl_cfd(:, 3) + vn, -(v_rot + vcyl_cfd(:, 2) + vt))
            ! Tip vortex within 45 deg of the rotor plane: axial residual, else tangential
            if (phi_vtp(n_v) < pi/4.d0) then
                resmax = maxval(section_res*cos(phi_cp_cyl))
            else
                resmax = maxval(section_res*sin(phi_cp_cyl))
            end if

            sub_iter = sub_iter + 1
            converged = (resmax < nw_max_resid) .or. (sub_iter == 500)
            if (sub_iter == 100) recompute_relax = .true.
            if (sub_iter == 500) then
                if (relax_safety < 0.8d0) then
                    relax_safety = relax_safety + 0.1d0
                    relax = min(relax + 0.1d0*(1.d0 - relax), 0.99d0)
                    recompute_relax = .true.
                end if
                ! No convergence: keep the near-wake induction of the last time step
                W = b%W_history
                vn = b%vn_history
                vt = b%vt_history
                dyn_comp = b%dyn_comp_history
                Gamma_steady = b%Gamma_history
                X = b%X_history
                Y = b%Y_history
            else
                W = W_iter
                vn = relax*vn + (1.d0 - relax)*Wn_iter
                vt = relax*vt + (1.d0 - relax)*Wt_iter
            end if
        end do

        ! 5. Blade forces with the corrected inflow
        vstn = nw2vstn(vt, vn, cone, pitch, sweep, fbend)
        un = vn_rot + vstn_cfd(:, 3) + vstn(:, 3)
        ut = vt_rot + vstn_cfd(:, 2) + vstn(:, 2)
        phi_cp = atan2(un, -ut)
        alpha_cp = phi_cp*r2d - twist
        alpha_cp_cfd = phi_cp_cfd*r2d - twist
        vrel_cp = sqrt(un**2 + ut**2)
        ! Ang Li correction for coned or bent blades (Li et al., 2022): lift and drag
        ! at the 3/4-chord AoA, force direction at the 1/4-chord inflow angle
        if (l1point) then
            thetadot = -omega*scone*r2d
            alpha_cp_34 = alpha_cp + (0.75d0 - x_aoa)*chord/vrel_cp*thetadot
            alpha_cp_14 = alpha_cp + (0.25d0 - x_aoa)*chord/vrel_cp*thetadot
        else
            alpha_cp_34 = alpha_cp
            alpha_cp_14 = alpha_cp
        end if
        phi_cp_14 = (alpha_cp_14 + twist)*d2r
        do l = 1, n_s
            cl(l) = nw_interp1d(alpha_2d(:, l), cl_2d(:, l), alpha_cp_34(l))
            cd(l) = nw_interp1d(alpha_2d(:, l), cd_2d(:, l), alpha_cp_34(l))
        end do
        lift = 0.5d0*density*vrel_cp**2*chord*cl
        drag = 0.5d0*density*vrel_cp**2*chord*cd
        Fstn(:, 1) = 0.d0
        Fstn(:, 2) = lift*sin(phi_cp_14) - drag*cos(phi_cp_14)
        Fstn(:, 3) = lift*cos(phi_cp_14) + drag*sin(phi_cp_14)

        ! 6. Diagnostics line and history for the next time step
        if (b%unit /= 0) then
            call cpu_time(cpu_t2)
            write(b%unit, '(es15.7, i5, *(es15.7))') time, sub_iter, resmax, &
                alpha_cp_cfd, alpha_cp, phi_vtp, vt, vn, dGamma, cpu_t2 - cpu_t1
        end if
        b%W_history = W
        b%dyn_comp_history = dyn_comp
        b%Gamma_history = Gamma_steady
        b%vn_history = vn
        b%vt_history = vt
        b%X_history = X
        b%Y_history = Y

        end associate
    end subroutine nw_smearing_correction

    !-------------------------------------------------------------------
    ! Smearing factor of a vortex at its trailing point (beta = 0).
    ! Meyer Forsting et al. (2019) Eq. (21)-(22) with beta = 0, as in
    ! Meyer Forsting et al. (2020) Sect. 2.3.
    !-------------------------------------------------------------------
    function smearing_factor(r, h_r, delta, phi, epsilon) result(smear_fac)
        implicit none

        ! --- Input  ---
        real*8, intent(in) :: r, h_r        ! Trailing-point radius [m], h/r [-]
        real*8, intent(in) :: delta         ! Distance normal to the disc (+ve upstream) [m]
        real*8, intent(in) :: phi           ! Helix angle [rad]
        real*8, intent(in) :: epsilon       ! Smearing length [m]

        ! --- Output ---
        real*8 :: smear_fac

        smear_fac = exp(-((r*h_r)**2 + (delta*cos(phi))**2)/epsilon**2)
    end function smearing_factor

    !-------------------------------------------------------------------
    ! Constants of one vortex-section pair, h = r_vtp - r_cp, h_r = h/r_vtp.
    ! Pirrung et al. (2016) and (2017b).
    !-------------------------------------------------------------------
    subroutine nw_hr_components(h, h_r, phi, phi_s, dx, dy, a)
        implicit none

        ! --- Input  ---
        real*8, intent(in) :: h, h_r

        ! --- Output ---
        real*8, intent(out) :: phi, phi_s       ! Root-corrected Phi_C, Phi of a straight vortex
        real*8, intent(out) :: dx, dy           ! Amplitudes without Gamma and Phi
        real*8, intent(out), dimension(4) :: a  ! Fit coefficients a_hr,i of the mixing factor

        ! --- Dummy variables ---
        real*8 :: dw, root_corr

        ! Phi of Wang and Coton, 2016 Eq. (4)
        if (h_r > 0.d0) then
            phi = pi/4.d0*abs((1.d0 + 0.5d0*h_r)*log(1.d0 - h_r))
        else
            phi = log(1.d0 - h_r)/(1.5d0 + log(1.d0 - 0.5d0*h_r))
        end if

        ! Amplitudes, 2016 Eq. (25). Important: D_Y is negative.
        dw = 1.d0/(4.d0*pi*abs(h)*h_r)
        dx = 1.359d0*dw
        dy = -0.359d0/4.d0*dw

        ! Root correction, 2016 Eq. (27)-(29). Important: applied once, through Phi (Phi_C, Eq. 31).
        root_corr = (dx*(1.d0 - exp(-pi/(2.d0*phi))) + dy*(1.d0 - exp(-2.d0*pi/phi)))/(dx + dy)
        phi = root_corr*phi

        ! Phi of a straight vortex, 2017b Eq. (7)
        phi_s = 0.788d0*abs(h_r)

        ! Fit coefficients of the mixing factor, 2017b Eq. (10) and (12)
        if (h_r < 0.d0) then
            a = nw_N(:, 1) + nw_N(:, 2)*exp(nw_N(:, 3)*h_r) + nw_N(:, 4)*exp(nw_N(:, 5)*h_r) &
                - nw_N(:, 4) - nw_N(:, 2)
        else
            a = nw_P(:, 1) + nw_P(:, 2)*h_r + nw_P(:, 3)*h_r**2 + nw_P(:, 4)*h_r**3
        end if
    end subroutine nw_hr_components

    !-------------------------------------------------------------------
    ! Phi* for the helix angle phi: weighted mix of straight and circular
    ! Phi. Pirrung et al. (2017b) Eq. (8)-(11), with the corrigendum's Phi_C.
    !-------------------------------------------------------------------
    function nw_phi_components(h_r, a, phi, phi_c, phi_s) result(phi_star)
        implicit none

        ! --- Input  ---
        real*8, intent(in) :: h_r
        real*8, intent(in), dimension(4) :: a   ! Fit coefficients a_hr,i
        real*8, intent(in) :: phi               ! Helix angle [rad]
        real*8, intent(in) :: phi_c, phi_s      ! Root-corrected circular Phi, straight Phi

        ! --- Output ---
        real*8 :: phi_star

        ! --- Dummy variables ---
        real*8 :: k_phi

        if (h_r < 0.d0) then
            k_phi = a(1) + a(2)*exp(a(3)*(pi/2.d0 - phi)) + a(4)*exp(-8.d0*(pi/2.d0 - phi)) &
                    - a(2) - a(4)
        else
            k_phi = a(1) + a(2)*phi + a(3)*phi**2 + a(4)*phi**3
        end if
        ! The fit is only valid below 89 degrees: keep the weight in [0, 1]
        k_phi = max(min(k_phi, 1.d0), 0.d0)
        phi_star = k_phi*phi_s + (1.d0 - k_phi)*phi_c
    end function nw_phi_components

    !-------------------------------------------------------------------
    ! Advance the slow (X) and fast (Y) indicial functions by one time
    ! step: new element plus the decayed history.
    ! Pirrung et al. (2016) Eq. (25)-(26) with Phi* of (2017b).
    !-------------------------------------------------------------------
    subroutine nw_xy_components(phi_star, dx, dy, dbeta, dGamma, X_history, Y_history, X, Y)
        implicit none

        ! --- Input  ---
        real*8, intent(in) :: phi_star, dx, dy, dbeta, dGamma
        real*8, intent(in) :: X_history, Y_history

        ! --- Output ---
        real*8, intent(out) :: X, Y

        ! --- Dummy variables ---
        real*8 :: X_exp, Y_exp

        X_exp = exp(-dbeta/phi_star)
        Y_exp = exp(-4.d0*dbeta/phi_star)
        X = dGamma*dx*phi_star*(1.d0 - X_exp) + X_history*X_exp
        Y = dGamma*dy*phi_star*(1.d0 - Y_exp) + Y_history*Y_exp
    end subroutine nw_xy_components

    !-------------------------------------------------------------------
    ! Bound circulation, its time-filtered value and the trailed strengths
    !-------------------------------------------------------------------
    subroutine nw_circulation(dt, vrel, alpha, Gamma_history, dyn_comp_history, alpha_2d, cl_2d, &
                              chord, Gamma_steady, Gamma, dGamma, dyn_comp)
        implicit none

        ! --- Input  ---
        real*8, intent(in) :: dt                                    ! Time step [s]
        real*8, intent(in), dimension(:) :: vrel, alpha             ! (n_s), alpha [deg]
        real*8, intent(in), dimension(:) :: Gamma_history           ! Gamma_steady of the last step
        real*8, intent(in), dimension(:, :) :: dyn_comp_history     ! (n_s, 3)
        real*8, intent(in), dimension(:, :) :: alpha_2d, cl_2d      ! Polars (len_2d, n_s)
        real*8, intent(in), dimension(:) :: chord

        ! --- Output ---
        real*8, intent(out), dimension(size(vrel)) :: Gamma_steady, Gamma
        real*8, intent(out), dimension(size(vrel)+1) :: dGamma
        real*8, intent(out), dimension(size(vrel), 3) :: dyn_comp

        ! --- Dummy variables ---
        integer :: l, n
        real*8, dimension(size(vrel)) :: cl, t

        ! Quasi-steady circulation from the static polar
        n = size(vrel)
        do l = 1, n
            cl(l) = nw_interp1d(alpha_2d(:, l), cl_2d(:, l), alpha(l))
        end do
        Gamma_steady = 0.5d0*vrel*chord*cl

        ! Time filter: Madsen and Gaunaa (2004), Pirrung et al. (2017a) Eq. (15)-(19).
        ! The constants are A_j/2 and b_j; A3 = 0.2625 (paper: 0.2656) so that the A_j sum to 1.
        t = dt/(chord/(2.d0*vrel))
        dyn_comp(:, 1) = dyn_comp_history(:, 1)*exp(-0.3064d0*t) &
                         + 0.27735d0*(Gamma_steady + Gamma_history)*(1.d0 - exp(-0.3064d0*t))
        dyn_comp(:, 2) = dyn_comp_history(:, 2)*exp(-0.0439d0*t) &
                         + 0.0914d0*(Gamma_steady + Gamma_history)*(1.d0 - exp(-0.0439d0*t))
        dyn_comp(:, 3) = dyn_comp_history(:, 3)*exp(-3.227d0*t) &
                         + 0.13125d0*(Gamma_steady + Gamma_history)*(1.d0 - exp(-3.227d0*t))
        Gamma = dyn_comp(:, 1) + dyn_comp(:, 2) + dyn_comp(:, 3)

        ! Trailed strengths: the end vortices take all Gamma of the end sections
        dGamma(1) = Gamma(1)
        dGamma(2:n) = Gamma(2:n) - Gamma(1:n-1)
        dGamma(n+1) = -Gamma(n)
    end subroutine nw_circulation

    !-------------------------------------------------------------------
    ! Relaxation factor of the sub-iteration, Pirrung et al. (2017a)
    ! Eq. (27)-(31) with the corrected Eq. (28).
    !-------------------------------------------------------------------
    function nw_calc_relax(vrel, alpha, vn, vt, chord, dt, dx, dy, dbeta, phi, safety) result(relax)
        implicit none

        ! --- Input  ---
        real*8, intent(in), dimension(:) :: vrel, alpha     ! (n_s), alpha [deg]
        real*8, intent(in), dimension(:) :: vn, vt          ! Normal and tangential speed (vt > 0)
        real*8, intent(in), dimension(:) :: chord
        real*8, intent(in) :: dt
        real*8, intent(in), dimension(:, :) :: dx, dy       ! (n_v, n_s)
        real*8, intent(in), dimension(:) :: dbeta           ! (n_v)
        real*8, intent(in), dimension(:, :) :: phi          ! Root-corrected Phi_C (n_v, n_s)
        real*8, intent(in) :: safety                        ! Safety margin

        ! --- Output ---
        real*8 :: relax

        ! --- Dummy variables ---
        integer :: l
        real*8 :: tau, d, A1, A2, gradw

        relax = 0.d0
        do l = 1, size(vrel)
            ! Gamma_dyn/Gamma_QS: Eq. (15) averaged over the time step
            tau = (2.d0*vrel(l))/chord(l)
            d = 1 - 0.5547d0/(tau*dt*0.3064d0)*(1 - exp(-tau*dt*0.3064d0)) &
                  - 0.1828d0/(tau*dt*0.0439d0)*(1 - exp(-tau*dt*0.0439d0)) &
                  - 0.2625d0/(tau*dt*3.227d0)*(1 - exp(-tau*dt*3.227d0))
            ! B_v of the inner (l) and outer (l+1) vortex of section l, Eq. (28)
            A1 = dx(l, l)*phi(l, l)*(1.d0 - exp(-dbeta(l)/phi(l, l))) &
                 + dy(l, l)*phi(l, l)*(1.d0 - exp(-4.d0*dbeta(l)/phi(l, l)))
            A2 = dx(l+1, l)*phi(l+1, l)*(1.d0 - exp(-dbeta(l+1)/phi(l+1, l))) &
                 + dy(l+1, l)*phi(l+1, l)*(1.d0 - exp(-4.d0*dbeta(l+1)/phi(l+1, l)))
            ! Gradient, Eq. (27): alpha converted to radians
            gradw = d*pi*chord(l)*(A1 - A2)*((alpha(l)*pi/180.d0*vn(l))/vrel(l) &
                    + vrel(l)/(abs(vt(l))*((vn(l)/abs(vt(l)))**2 + 1.d0)))
            relax = max(relax, -(1.d0 + gradw)/(1.d0 - gradw))      ! Eq. (31)
        end do
        relax = relax + safety*(1.d0 - relax)
        relax = min(relax, 0.99d0)
    end function nw_calc_relax

    !-------------------------------------------------------------------
    ! Linear interpolation in the ascending table x, y at xin,
    ! linear extrapolation outside the table
    !-------------------------------------------------------------------
    function nw_interp1d(x, y, xin) result(yin)
        implicit none

        ! --- Input  ---
        real*8, intent(in), dimension(:) :: x, y
        real*8, intent(in) :: xin

        ! --- Output ---
        real*8 :: yin

        ! --- Dummy variables ---
        integer :: i, n

        n = size(x)
        if (xin < x(1)) then
            yin = y(1) + (xin - x(1))/(x(2) - x(1))*(y(2) - y(1))
        else if (xin > x(n)) then
            yin = y(n) + (xin - x(n))/(x(n) - x(n-1))*(y(n) - y(n-1))
        else
            i = 2
            do while (x(i) < xin)
                i = i + 1
            end do
            yin = y(i-1) + (xin - x(i-1))/(x(i) - x(i-1))*(y(i) - y(i-1))
        end if
    end function nw_interp1d

    !-------------------------------------------------------------------
    ! Correction velocities (cylindrical vt, vn) in section axes (n, 3)
    !-------------------------------------------------------------------
    function nw2vstn(vt, vn, cone, pitch, sweep, fbend) result(vstn)
        implicit none

        ! --- Input  ---
        real*8, intent(in), dimension(:) :: vt, vn
        real*8, intent(in) :: cone, pitch
        real*8, intent(in), dimension(:) :: sweep, fbend

        ! --- Output ---
        real*8, dimension(size(vt), 3) :: vstn

        ! --- Dummy variables ---
        integer :: l

        do l = 1, size(vt)
            vstn(l, :) = vxyz_vstn([0.d0, vt(l), vn(l)], 0.d0, 0.d0, 0.d0, cone, pitch, &
                                   sweep(l), fbend(l))
        end do
    end function nw2vstn

    !-------------------------------------------------------------------
    ! Rotate a vector from the Cartesian or cylindrical frame to section
    ! axes. Order: yaw, tilt, azimuth, cone, pitch, sweep, bend [deg].
    !-------------------------------------------------------------------
    function vxyz_vstn(a, yaw, tilt, azimut, cone, pitch, sweep, fbend) result(v)
        implicit none

        ! --- Input  ---
        real*8, intent(in), dimension(3) :: a
        real*8, intent(in) :: yaw, tilt, azimut, cone, pitch, sweep, fbend

        ! --- Output ---
        real*8, dimension(3) :: v

        v = rot_y(a, yaw)
        v = rot_x(v, tilt)
        v = rot_z(v, -azimut)
        v = rot_y(v, cone)
        v = rot_x(v, -pitch)
        v = rot_z(v, -sweep)
        v = rot_y(v, fbend)
    end function vxyz_vstn

    ! Rotation about the first axis (tilt, pitch), angle in degrees
    function rot_x(a, angle) result(v)
        implicit none
        real*8, intent(in), dimension(3) :: a
        real*8, intent(in) :: angle
        real*8, dimension(3) :: v
        real*8 :: c, s

        c = cos(angle*pi/180.d0)
        s = sin(angle*pi/180.d0)
        v = [a(1), c*a(2) + s*a(3), -s*a(2) + c*a(3)]
    end function rot_x

    ! Rotation about the second axis (yaw, cone, bend), angle in degrees
    function rot_y(a, angle) result(v)
        implicit none
        real*8, intent(in), dimension(3) :: a
        real*8, intent(in) :: angle
        real*8, dimension(3) :: v
        real*8 :: c, s

        c = cos(angle*pi/180.d0)
        s = sin(angle*pi/180.d0)
        v = [c*a(1) - s*a(3), a(2), s*a(1) + c*a(3)]
    end function rot_y

    ! Rotation about the third axis (azimuth, sweep), angle in degrees
    function rot_z(a, angle) result(v)
        implicit none
        real*8, intent(in), dimension(3) :: a
        real*8, intent(in) :: angle
        real*8, dimension(3) :: v
        real*8 :: c, s

        c = cos(angle*pi/180.d0)
        s = sin(angle*pi/180.d0)
        v = [c*a(1) - s*a(2), s*a(1) + c*a(2), a(3)]
    end function rot_z

end module al_nw_smearing_correction
