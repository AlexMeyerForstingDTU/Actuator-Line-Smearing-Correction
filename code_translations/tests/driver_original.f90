!=======================================================================
! Test driver for the original fixed-form al_nw_smearing_correction.f
!
! Usage: driver_original <case>
! Runs one test case in the current directory and writes results.txt.
! The correction itself writes the diagnostics files test.nwRR.blBB.
!=======================================================================
program driver_original
    use test_common
    use ns_contr, only: timestep0, grlev, density, time
    use nw_smearing_correction_module, only: kaya, nw_max_resid, nw_vn_history, nw_vt_history
    implicit none

    character(len=32) :: case_name
    integer :: results_unit

    call get_command_argument(1, case_name)
    open(newunit=results_unit, file='results.txt', status='replace')
    if (trim(case_name) == 'units') then
        call run_units(results_unit)
    else
        call run_case(get_case(trim(case_name)), results_unit)
    end if
    close(results_unit)

contains

    ! Unit tests of the individual routines
    subroutine run_units(u)
        implicit none

        ! --- Input  ---
        integer, intent(in) :: u

        ! --- Dummy variables ---
        real*8, parameter, dimension(7) :: hr_list = &
            [0.1d0, 0.3d0, 0.5d0, 0.7d0, 0.8d0, -0.5d0, -3.d0]
        real*8, parameter, dimension(4) :: phi_hr = [0.1d0, 0.3d0, -0.5d0, -3.d0]
        real*8, parameter, dimension(4) :: phi_deg = [5.d0, 30.d0, 10.d0, 60.d0]
        type(test_blade) :: b
        integer :: n, k, l
        real*8, parameter, dimension(2) :: relax_dt = [0.01d0, 0.1d0]
        real*8 :: phi, phi_s, dx, dy, X, Y, X_hist, Y_hist, ratio
        real*8, dimension(2) :: relax_fac
        real*8, dimension(3) :: f, v
        real*8, dimension(4) :: a, phi_star
        real*8, dimension(56) :: hr_values
        real*8, dimension(15) :: repeat_values
        real*8, dimension(7) :: yin
        real*8, dimension(1, 3) :: v1
        real*8, dimension(2, 3) :: v2
        real*8, dimension(ns) :: un, ut, vrel, alpha, G_steady, G, G_hist, vt_pos
        real*8, dimension(nv) :: dG, vrel_vtp, dbeta
        real*8, dimension(ns, 3) :: dyn, dyn_hist
        real*8, dimension(nv, ns) :: phi_pair, phi_s_pair, dx_pair, dy_pair
        real*8, dimension(nv, ns, 4) :: a_pair

        ! The original needs its module constants (fit coefficients, n_s)
        kaya = .true.
        call init_nw_smearing_correction(1, 1, ns)

        ! 1. Smearing factor, with and without an out-of-plane offset
        call smearing_factor(60.d0, 0.1d0, 0.d0, 0.d0, 10.d0, f(1))
        call smearing_factor(60.d0, 0.1d0, 2.d0, 0.3d0, 10.d0, f(2))
        call smearing_factor(10.d0, -0.5d0, 1.d0, 0.5d0, 2.5d0, f(3))
        call write_values(u, 'smearing_factor', f)

        ! 2. Geometry constants, r_vtp = 10 m: phi, phi_s, dx, dy, a(1:4) per h/r
        do n = 1, 7
            call nw_hr_components(10.d0*hr_list(n), hr_list(n), phi, phi_s, dx, dy, a)
            hr_values(8*n-7:8*n) = [phi, phi_s, dx, dy, a]
        end do
        call write_values(u, 'hr_components', hr_values)

        ! 3. Phi* for four (h/r, helix angle) pairs
        do n = 1, 4
            call nw_hr_components(10.d0*phi_hr(n), phi_hr(n), phi, phi_s, dx, dy, a)
            call nw_phi_components(phi_hr(n), a, phi_deg(n)*d2r, phi, phi_s, phi_star(n))
        end do
        call write_values(u, 'phi_components', phi_star)

        ! 4. Indicial functions: one step, then the steady state of 2000 steps
        call nw_XY_components(0.09d0, 0.0123d0, -0.0008125d0, 0.05d0, 2.d0, 0.3d0, -0.02d0, X, Y)
        call write_values(u, 'xy_single', [X, Y])
        X = 0.d0
        Y = 0.d0
        do n = 1, 2000
            X_hist = X
            Y_hist = Y
            call nw_XY_components(0.09d0, 0.0123d0, -0.0008125d0, 0.05d0, 2.d0, X_hist, Y_hist, &
                                  X, Y)
        end do
        call write_values(u, 'xy_steady', [X, Y])

        ! 5. Interpolation, inside and beyond both ends
        call nw_interp1d([0.d0, 1.d0, 2.d0, 3.d0], [0.d0, 10.d0, 15.d0, 30.d0], 4, &
                         [-1.d0, 0.d0, 0.5d0, 2.d0, 2.5d0, 3.d0, 4.d0], yin, 7)
        call write_values(u, 'interp1d', yin)

        ! 6. nw2vstn called repeatedly on a coned blade (note Sect. 4.1 and 6.1)
        do n = 1, 15
            call nw2vstn([0.d0], [-1.d0], 1, 2.5d0, 0.d0, [0.d0], [0.d0], v1)
            repeat_values(n) = v1(1, 3)
        end do
        call write_values(u, 'nw2vstn_repeat', repeat_values)
        call write_values(u, 'nw2vstn_repeat_last', v1(1, :))

        ! 7. nw2vstn and the full rotation chain with all angles non-zero
        call nw2vstn([1.5d0, -0.3d0], [-0.7d0, 2.d0], 2, 2.5d0, 3.d0, [4.d0, -2.d0], &
                     [-2.d0, 5.d0], v2)
        call write_values(u, 'nw2vstn_general', [v2(1, :), v2(2, :)])
        v = [1.d0, -2.d0, 0.5d0]
        call Vxyz_Vstn(v, 1, 10.d0, 5.d0, 30.d0, 2.5d0, 3.d0, 4.d0, -2.d0)
        call write_values(u, 'vxyz_vstn', v)

        ! 8. Circulation of the test blade without correction (cone 0, pitch 0)
        call build_blade(get_case('eps10'), b)
        un = wind*(1.d0 - a_ind)
        ut = -(omega*b%rc)
        vrel = sqrt(un**2 + ut**2)
        alpha = atan2(un, -ut)*r2d - b%twist
        G_hist = 0.d0
        dyn_hist = 0.d0
        call nw_circulation(0.1d0, vrel, alpha, G_hist, dyn_hist, b%alpha_2d, b%cl_2d, len2d, &
                            b%chord, G_steady, G, dG, dyn)
        call write_values(u, 'circulation_gamma_steady', G_steady)
        call write_values(u, 'circulation_gamma', G)
        call write_values(u, 'circulation_dgamma', dG)
        ! Steady state of the time filter: Gamma -> Gamma_steady when sum(A_j) = 1
        G_hist = 0.d0
        dyn_hist = 0.d0
        do n = 1, 300
            call nw_circulation(1.d0, vrel, alpha, G_hist, dyn_hist, b%alpha_2d, b%cl_2d, len2d, &
                                b%chord, G_steady, G, dG, dyn)
            G_hist = G_steady
            dyn_hist = dyn
        end do
        ratio = 0.d0
        do l = 1, ns
            if (abs(G_steady(l)) > 1.d-3) ratio = max(ratio, abs(G(l)/G_steady(l) - 1.d0))
        end do
        call write_values(u, 'circulation_ratio', [ratio])

        ! 9. Relaxation factor of the test blade at two time steps
        do k = 1, nv
            do l = 1, ns
                call nw_hr_components(b%rv(k) - b%rc(l), (b%rv(k) - b%rc(l))/b%rv(k), &
                    phi_pair(k, l), phi_s_pair(k, l), dx_pair(k, l), dy_pair(k, l), a_pair(k, l, :))
            end do
        end do
        call nw_interp1d(b%rc, vrel, ns, b%rv(2:nv-1), vrel_vtp(2:nv-1), ns-1)
        vrel_vtp(1) = vrel(1)
        vrel_vtp(nv) = vrel(ns)
        vt_pos = -ut
        do n = 1, 2
            dbeta = (vrel_vtp*relax_dt(n))/b%rv
            call nw_calc_relax(vrel, alpha, un, vt_pos, b%chord, relax_dt(n), dx_pair, dy_pair, &
                               dbeta, phi_pair, 0.2d0, relax_fac(n))
        end do
        call write_values(u, 'calc_relax', relax_fac)
    end subroutine run_units

    ! CFD velocities in section axes for a given cone and pitch
    function cfd_vstn(b, cone, pitch) result(vstn)
        implicit none

        ! --- Input  ---
        type(test_blade), intent(in) :: b
        real*8, intent(in) :: cone, pitch

        ! --- Output ---
        real*8, dimension(ns, 3) :: vstn

        ! --- Dummy variables ---
        integer :: l
        real*8, dimension(3) :: tmp

        do l = 1, ns
            tmp = b%vcyl(l, :)
            call Vxyz_Vstn(tmp, 1, 0.d0, 0.d0, 0.d0, cone, pitch, b%sweep(l), b%fbend(l))
            vstn(l, :) = tmp
        end do
    end function cfd_vstn

    ! March one test case in time and write its results
    subroutine run_case(c, u)
        implicit none

        ! --- Input  ---
        type(test_case), intent(in) :: c
        integer, intent(in) :: u

        ! --- Dummy variables ---
        type(test_blade) :: b
        integer :: it, i, j, n_series
        real*8 :: pitch, p, P_tot, T_tot
        real*8, dimension(ns, 3) :: vstn
        real*8, dimension(ns, 3, 2, 2) :: F
        real*8, dimension(110) :: P_series, T_series
        real*8, dimension(4) :: P_blade, T_blade

        ! 1. Inputs and host values
        call build_blade(c, b)
        kaya = .true.
        nw_max_resid = c%max_resid
        timestep0(1) = c%dt
        grlev = 1
        density = rho
        call init_nw_smearing_correction(c%n_rotors, c%n_blades, ns)

        ! 2. March in time
        n_series = 0
        do it = 1, c%n_steps
            time = it*c%dt
            pitch = 0.d0
            if (c%pitch_step_at > 0 .and. it >= c%pitch_step_at) pitch = 2.d0
            do i = 1, c%n_rotors
                do j = 1, c%n_blades
                    p = pitch
                    if (i == 2 .and. j == 2) p = c%pitch_22
                    vstn = cfd_vstn(b, c%cone, p)
                    call nw_smearing_correction(i, j, len2d, b%cyl_cp, b%cyl_vtp, vstn, b%vcyl, &
                        c%eps, omega, b%chord, b%twist, p, b%cl_2d, b%cd_2d, b%alpha_2d, &
                        c%cone, b%fbend, b%sweep, b%dsweep, c%l1point, F(:, :, i, j))
                end do
            end do
            if (c%pitch_step_at > 0 .and. it > c%n_steps - 110) then
                n_series = n_series + 1
                call power_thrust(b, F(:, :, 1, 1), 3.d0, P_series(n_series), T_series(n_series))
            end if
        end do

        ! 3. Results
        if (c%root_at_zero) then
            call write_values(u, 'kaya_after', [merge(1.d0, 0.d0, kaya)])
        else if (c%n_rotors > 1) then
            do i = 1, 2
                do j = 1, 2
                    call power_thrust(b, F(:, :, i, j), 1.d0, P_blade(2*i+j-2), T_blade(2*i+j-2))
                end do
            end do
            call write_values(u, 'P_blade', P_blade)
            call write_values(u, 'T_blade', T_blade)
        else
            call power_thrust(b, F(:, :, 1, 1), 3.d0, P_tot, T_tot)
            call write_values(u, 'P', [P_tot])
            call write_values(u, 'T', [T_tot])
            call write_values(u, 'Fstn2', F(:, 2, 1, 1))
            call write_values(u, 'Fstn3', F(:, 3, 1, 1))
            call write_values(u, 'vn_history', nw_vn_history(1, 1, :))
            call write_values(u, 'vt_history', nw_vt_history(1, 1, :))
            if (c%pitch_step_at > 0) then
                call write_values(u, 'P_series', P_series)
                call write_values(u, 'T_series', T_series)
            end if
        end if
    end subroutine run_case

end program driver_original
