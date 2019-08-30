c-----------------------------------------------------------------------
c            ####  #######  #    #     #     #  #  #    #  ####
c            #   #    #     #    #     #  #  #  #  # #  #  #   #
c            #   #    #     #    #     #  #  #  #  #  # #  #   #
c            ####     #      ####       #   #   #  #    #  ####
c-----------------------------------------------------------------------
c
c     A fast SMEARING CORRECTION for the actuator line 
c
c     Author:        Alexander R. Meyer Forsting, alrf@dtu.dk
c     Contributors:  Georg Pirrung,
c                    Néstor Ramos-García,
c                    Ang Li 
c     Revisions : 
c         date   |      author     |    task
c     - 24-07-19 | A MeyerForsting | Version 1.0
c
c     Purpose of file: Correct actuator line forces for the missing
c                      induction, originating from the force smearing at
c                      the control points. This code is an 
c                      implementation of the correction presented in 
c                      Meyer Forsting et al. (2019). 
c
c     Overview of functions contained in this file:
c     - module nw_smearing_correction_module
c     - subroutine nw_smearing_correction()
c     - subroutine smearing_factor()
c     - subroutine nw_hr_components()
c     - subroutine nw_phi_components()
c     - subroutine nw_XY_components()
c     - subroutine nw_circulation()
c     - subroutine nw_calc_relax()
c     - subroutine init_nw_smearing_correction()
c     - subroutine nw_interp1d()
c
c-----------------------------------------------------------------------
c     This file contains all elements necessary for running the smearing
c     correction for the actuator line (AL). The AL code only needs to 
c     call nw_smearing_correction() with the given in- and output 
c     variables. Note that it needs to be called for each blade 
c     individually. Call init_* before calling the routine to initialize 
c     all variables. Activate the correction by setting kaya=.true. and 
c     set a convgerence limit (nw_max_resid, 1d-6 for example). When 
c     using deformation some constants might have to be recomputed which 
c     can be achieved through setting kaya_refresh of the specific rotor 
c     and blade to .true.. Variables from previous time steps are 
c     stored for all rotors and blades in single arrays (variables 
c     ending on *_history). 
c
c     This implementation is based on several scientific publications
c     related to the smearing correction and the near-wake model, which
c     are referenced throughout the code:
c
c     Meyer Forsting et al. (2019)
c        Meyer Forsting, A. R., Pirrung, G. R., and Ramos-García, N.: 
c        A vortex-based tip/smearing correction for the actuator line, 
c        Wind Energ. Sci., 4, 369-383, 
c        https://doi.org/10.5194/wes-4-369-2019, 2019. 
c
c     Pirrung et al. (2016)
c        Pirrung, G., Madsen, H. A., Kim, T., and Heinz, J.: 
c        A coupled near and far wake model for wind turbine aerodynamics
c        , Wind Energy, 19, 2053–2069, 
c        https://doi.org/10.1002/we.1969, 2016. 
c
c     Pirrung et al. (2017a)
c        Pirrung, G., Riziotis, V., Madsen, H., Hansen, M., and Kim, T.: 
c        Comparison of a coupled near- and far-wake model with a free-
c        wake vortex code, Wind Energ. Sci., 2, 15–33, 
c        https://doi.org/10.5194/wes-2-15-2017, 2017. 
c
c     Pirrung et al. (2017b)
c     Pirrung, G. R., Madsen, H. A., and Schreck, S.: 
c        Trailed vorticity modeling for aeroelastic wind turbine 
c        simulations in standstill, Wind Energ. Sci., 2, 521–532, 
c        https://doi.org/10.5194/wes-2-521-2017, 2017. 
c
c-----------------------------------------------------------------------

c=======================================================================
      module nw_smearing_correction_module 
c-----------------------------------------------------------------------
c     Variables used in the new near wake tip correction
c----------------------------------------------------------------------- 
      ! Global precision
      use params,only:idp 
      ! Switch activating the correction
      logical::kaya=.false.  
      ! Switch invocking the recomputation of constants (large changes 
      ! in the blade geometry from deformation could make it necessary)
      logical,dimension(:,:),allocatable::kaya_refresh 
      ! Definition of setup: number of rotors, nr. of blades, nr. of 
      ! vortex trailing points (in-between blade sections), nr. of 
      ! sections along blade
      integer::n_rotors=0,n_blades=0,n_v=0,n_s=0
      ! Index range of control points (CP): The missing induction is 
      ! only felt in the proximity of the shed vortices, allowing to 
      ! limit the number of CPs to loop over.
      integer,dimension(:,:,:,:),allocatable::cp_loop_ind 
      ! Pi
      real(kind=idp)::pi=4.d0*atan(1.d0)
      ! Convergence limit of induced velocity residual (default value)
      real(kind=idp)::nw_max_resid=1d-6
      ! Constants of the Near-wake model
      real(kind=idp),dimension(4,5)::nw_N,nw_P
      ! Variables from previous time step: For smooth transition between
      ! steps history needs to be stored
      real(kind=idp),dimension(:,:,:),allocatable::
     &         nw_vn_history,nw_vt_history,nw_Gamma_history,nw_W_history          
      real(kind=idp),dimension(:,:,:,:),allocatable::
     &                     nw_dyn_comp_history,nw_X_history,nw_Y_history
      ! Constants that need to be recomputed only after strong changes 
      ! in geometry
      real(kind=idp),dimension(:,:,:,:),allocatable::     
     &                             nw_phi,nw_phi_s,nw_dx,nw_dy,smear_fac  
      real(kind=idp),dimension(:,:,:,:,:),allocatable:: nw_a   
      end module 
c=======================================================================
c=======================================================================
      subroutine nw_smearing_correction(i,j,
     &                                  r_cp,r_vtp,vn_cfd,vt_cfd,
     &                                  eps_cp,omega,chord,twist,pitch,
     &                                  cl_2d,cd_2d,alpha_2d,len_2d,
     &                                  Fstn)
c------------------------------------------------------
c     The smearing function leads to a loss of vorticity,
c     which in turn leads to lower induced velocities. The
c     missing velocity is calculated by a near wake model
c     and used to correct the angle of attack.
c-----------------------------------------------------------------------
      use params,only:idp          ! Global precision 
      use ns_contr,only:timestep0, ! Global time step
     &                  grlev,     ! Grid level
     &                  density,   ! Air density
     &                  time       ! Global time 
      use nw_smearing_correction_module 
c-----------------------------------------------------------------------      
      implicit none
c-----[in]   
      integer,intent(in)::i,     ! ith rotor
     &                    j,     ! jth blade
     &                    len_2d ! Length of 2D airfoil data
      real(kind=idp),intent(in)::
     & r_cp(n_s),   ! Radial coordinate of control points (CP)
     & r_vtp(n_v),  ! Radial coordinate of vortex trailing points   
     & vn_cfd(n_s), ! Normal velocities at CPs from CFD domain
     & vt_cfd(n_s), ! Tangential velocities at CPs from CFD domain
     & eps_cp,      ! Epsilon, the smearing length scale
     & omega,       ! Rotational speed in radiens/sec
     & chord(n_s),  ! Chord distribution along blade at CPs 
     & twist(n_s),  ! Twist distribution along blade at CPs
     & pitch,       ! Blade pitch
     & cl_2d(len_2d,n_s),    ! Airfoil lift coeff. curve at CPs
     & cd_2d(len_2d,n_s),    ! Airfoil drag coeff. curve at CPs
     & alpha_2d(len_2d,n_s)  ! Alpha values of lift&drag curves
      ! Note that the radial coordiante is defined along the blade and
      ! starts at the rotor centre. The radial coordinate should be 
      ! ascending and r_vtp(1)>0.
c-----[out] Forces at CPs, 1:spanwise, 2:tangential, 3:normal 
      real(kind=idp),intent(out)::Fstn(n_s,3)
c-----[tmp]
      logical :: nw_recompute_relax,nw_converged
      integer :: nw_sub_iter, k, l
      real(kind=idp)::dt_global, h, h_r, nw_resmax, nw_relax,
     &                nw_relax_safety,core_cutoff,
     &                rtc,cpu_t1,cpu_t2      
      real(kind=idp),dimension(n_s):: 
     &               cl_inter, cd_inter, Gamma_cp, 
     &               Gamma_cp_steady, phi_cp_cfd, alpha_cp_cfd,
     &               lift_cp, drag_cp, nw_section_res, nw_vn,nw_vt,
     &               nw_W, phi_cp, alpha_cp, vrel_cp, nw_W_iter
      real(kind=idp),dimension(n_v)::
     &               phi_vtp,vrel_vtp, dGamma_vtp,nw_dbeta
      real(kind=idp)::nw_dyn_comp(n_s,3),
     &                nw_X(n_v,n_s),nw_Y(n_v,n_s),nw_phi_star(n_v,n_s)
c=======================================================================
      if(.not.kaya)return ! Switch model on/off
      ! Test whether the first coordinate is close to zero. This is a
      ! singularity.
      if(r_vtp(1)<0.1d0)then
       write(6,*) 
       write(6,*)
       write(6,"('!!!! WARNING smearing correction switched off !!!!')")
       write(6,"('First vortex trailing point is zero, please modify')")
       write(6,*)
       write(6,*)
       kaya =.false.
       return
      endif
      !#################### INITILIZE ##################################   

      !========= Starting Preamble =====================================
      nw_relax_safety = 0.2d0
      nw_recompute_relax = .true.
      nw_converged = .false.
      ! Get the global timestep of the current grid level
      dt_global = timestep0(grlev)
      ! Get velocities from the NW model from previous time step
      nw_W = nw_W_history(i,j,:)
      nw_vn =  nw_vn_history(i,j,:)
      nw_vt =  nw_vt_history(i,j,:)

      !+++++++++++++++++++++ ONLY IF REFRESH NEEDED ++++++++++++++++++++
      if (kaya_refresh(i,j))then 
        !============= Compute constant nw components ==================
        ! These components only need to be computed once at the very 
        ! first iteration or when the geometry of the blade changes 
        ! drastically.  

        ! Viscous core cut-off at 99th percentile of Gaussian. Outside 
        ! this region the missing induction from the vortex smearing 
        ! is neglegible.
        core_cutoff = 1.82d0*eps_cp   

          do k = 1, n_v
            do l = 1,n_s
              !---------------------------------------------------------
              ! Distance between vortex and control point
              h = r_vtp(k)-r_cp(l)
              h_r = h/r_vtp(k)
              ! Compute constants of the near-wake model 
              call nw_hr_components(h,h_r,nw_phi(i,j,k,l),
     &                                nw_phi_s(i,j,k,l),nw_dx(i,j,k,l),
     &                                nw_dy(i,j,k,l),nw_a(i,j,k,l,:))
              ! Compute the smearing factor, which is fixed in this 
              ! implementation i.e. beta=0
              call smearing_factor(r_vtp(k),h_r,eps_cp,
     &                             smear_fac(i,j,k,l))
              ! Finally determine which blade sections are influenced
              ! by the vtpex core. For a Gaussian the 99th percentile
              ! is reached at 1.82x the smearing factor. Further from
              ! the influence of the core is assumed negligible.
              if (h.gt.0d0)then
                if (abs(h).gt.core_cutoff)then
                  cp_loop_ind(i,j,k,1) = l+1
                endif
              else
                if (abs(h).lt.core_cutoff)then
                  cp_loop_ind(i,j,k,2) = l
                endif
              endif
              !--------------------------------------------------------- 
            enddo
          enddo  
        ! Now it's fresh
        kaya_refresh(i,j) =.false.
      endif ! refresh?    
      !######################### MAIN ##################################
      cpu_t1=rtc() ! Start timing 
          !===================  Fixed components =======================
          ! The physical flow angle does not change during a CFD
          ! iteration. Therefore there are no sub-iterations needed. 
          !----Flow/helix angle----
          ! Extracted from CFD at control points
          phi_cp_cfd = atan2(vn_cfd,omega*r_cp-vt_cfd) 
          ! Ensure that the angle is below |90| degrees
          do l = 1,n_s
            if (phi_cp_cfd(l).gt.pi/2d0) then 
              phi_cp_cfd(l)=pi-phi_cp_cfd(l) 
            endif
          enddo
          ! Trailing points
          ! Interpolate from CPs
          call nw_interp1d(r_cp,phi_cp_cfd,n_s,
     &                    r_vtp(2:n_v-1),phi_vtp(2:n_v-1),n_s-1)            
          phi_vtp(1) = phi_cp_cfd(1)
          phi_vtp(n_v) = phi_cp_cfd(n_s)

          !----Near Wake phi components---- 
          do k = 1,n_v
            do l =  cp_loop_ind(i,j,k,1),cp_loop_ind(i,j,k,2)
              ! Distance between trailing point and CP
              h = r_vtp(k)-r_cp(l)
              h_r = h/r_vtp(k)
              call nw_phi_components(h_r,nw_a(i,j,k,l,:),phi_vtp(k),
     &                               nw_phi(i,j,k,l),nw_phi_s(i,j,k,l),
     &                               nw_phi_star(k,l))          
            enddo
          enddo  

          !************* Iterative procedure ***************************
          ! Iterate until correctinon velocities converge 
          nw_sub_iter = 0 ! Sub-iteration counter
          do while (.not.nw_converged)

            !============= Velocities & Flow angles  ===================
            !----Control points----
            ! Corrected flow angle
            phi_cp = atan2(vn_cfd+nw_vn,
     &                     omega*r_cp-vt_cfd-nw_vt)
            ! Angle-of-attack at CP
            alpha_cp = phi_cp*180.d0/pi-twist-pitch

            ! Relative wind speed at the control points
            vrel_cp = sqrt((vn_cfd+nw_vn)**2 + 
     &                     (omega*r_cp-vt_cfd-nw_vt)**2)
            !----Trailing points----
            ! Interpolate the velocity at the trailing points
            call nw_interp1d(r_cp,vrel_cp,n_s,
     &                    r_vtp(2:n_v-1),vrel_vtp(2:n_v-1),n_s-1)
            ! Take velcoities from outer CPs (typical in LL codes)
            vrel_vtp(1) = vrel_cp(1)
            vrel_vtp(n_v) = vrel_cp(n_s)

            !================ Circulation ==============================
            ! Compute circulation at sections 
            call nw_circulation(dt_global,vrel_cp,alpha_cp,
     &             nw_Gamma_history(i,j,:),nw_dyn_comp_history(i,j,:,:),
     &             alpha_2d,cl_2d,len_2d,chord,
     &             Gamma_cp_steady,Gamma_cp,dGamma_vtp,nw_dyn_comp)
            !++++++++++++++++++ Velocity history +++++++++++++++++++++++
            ! Local advance in beta (rotation)
            nw_dbeta = (vrel_vtp*dt_global)/r_vtp          
            !============ Calculate new velocity =======================
            ! Determine the velocity induced by vorticity shed during  
            ! the current time step. 
            nw_W_iter = 0d0 ! Initialize
            nw_X = 0d0; nw_Y = 0d0
            do k = 1,n_v
               do l =  cp_loop_ind(i,j,k,1),cp_loop_ind(i,j,k,2)
               ! Distance between trailing point and CP
               h = r_vtp(k)-r_cp(l)
               h_r = h/r_vtp(k)

               !----Near Wake velocity components---- 
               call nw_XY_components(nw_phi_star(k,l),nw_dx(i,j,k,l),
     &                 nw_dy(i,j,k,l),nw_dbeta(k),dGamma_vtp(k),
     &                 nw_X_history(i,j,k,l),nw_Y_history(i,j,k,l),
     &                 nw_X(k,l),nw_Y(k,l))

               !----Velocity correction----
               nw_W_iter(l) = nw_W_iter(l)+ smear_fac(i,j,k,l)*
     &                       (nw_X(k,l)+nw_Y(k,l))
               enddo
            enddo
            !============== Relaxation factor ==========================
            if (nw_recompute_relax)then
              call nw_calc_relax(vrel_cp,alpha_cp,
     &                     (omega*r_cp-vt_cfd-nw_vt),(vn_cfd+nw_vn),
     &                      chord,dt_global,nw_dx,nw_dy,nw_dbeta,nw_phi,
     &                      nw_relax_safety,nw_relax)
              ! After computing it once 
              nw_recompute_relax = .false.
            endif

            !============== Advance velocities =========================
            ! Advance velocity correction with relaxation
            nw_W_iter = nw_relax*nw_W+(1d0-nw_relax)*nw_W_iter
            !=================== Residual ==============================
            ! Tip vortex less than 45 degrees away from rotor plane: 
            ! convergence determined by axial induction otherwise 
            ! convergence determined by tangential induction  
            nw_section_res=sqrt((nw_W_iter-nw_W)**2)

            if(phi_vtp(n_v).lt.pi/4d0)then 
              nw_resmax = maxval(nw_section_res*cos(phi_cp))
            else
              nw_resmax = maxval(nw_section_res*sin(phi_cp))
            endif

            !=============== Convergence ===============================
            nw_sub_iter = nw_sub_iter + 1
            !----Check----
            if ((nw_resmax.lt.nw_max_resid).or.
     &                                        (nw_sub_iter.eq.500)) then
              nw_converged = .true.
            endif

            !----100 Iterations---- 
            if (nw_sub_iter.eq.100)then
              nw_recompute_relax = .true.
            endif
            
            !----500 Iterations---- 
            ! If no convergence was reached after 500 iterations change 
            ! the relaxation and skip this time step
            if (nw_sub_iter.eq.500)then
              if(nw_relax_safety.lt.0.8d0) then
                nw_relax_safety=nw_relax_safety+0.1d0
                nw_relax=nw_relax+0.1d0*(1.0d0-nw_relax)
                nw_relax=min(nw_relax, 0.99d0)
                nw_recompute_relax=.true.
              endif
              !No convergence, use near wake induction of last time step
              nw_W = nw_W_history(i,j,:)
              nw_vn =  nw_vn_history(i,j,:)
              nw_vt =  nw_vt_history(i,j,:)
              nw_dyn_comp = nw_dyn_comp_history(i,j,:,:)
              Gamma_cp_steady = nw_Gamma_history(i,j,:)
              nw_X = nw_X_history(i,j,:,:) 
              nw_Y = nw_Y_history(i,j,:,:)

            else
            !============== While loop update===========================
              nw_W = nw_W_iter 
              nw_vn = nw_W*cos(phi_cp)
              nw_vt = nw_W*sin(phi_cp)
            endif
          enddo 
          !********************* END ITERATING  ************************

          !===================== Blade forces ==========================

          !----Flow/helix angle----
          ! Corrected 
          phi_cp = atan2(vn_cfd+nw_vn,omega*r_cp-vt_cfd-nw_vt)
          !----Angle-of-Attack----
          ! Corrected 
          alpha_cp = phi_cp*180.d0/pi-twist-pitch
          ! CFD
          alpha_cp_cfd = phi_cp_cfd*180.d0/pi-twist-pitch

          !----Lift & drag coefficients----
          do l= 1,n_s
           call nw_interp1d(alpha_2d(:,l),cl_2d(:,l),len_2d,alpha_cp(l),
     &                                        cl_inter(l),1)
           call nw_interp1d(alpha_2d(:,l),cd_2d(:,l),len_2d,alpha_cp(l),
     &                                  cd_inter(l),1)
          enddo
          !----Relative velocity----
          vrel_cp = sqrt((vn_cfd+nw_vn)**2 + 
     &                   (omega*r_cp-vt_cfd-nw_vt)**2)

          !----Aerodynamic forces----
          ! Lift 
          lift_cp = 0.5d0*density*vrel_cp**2*chord*cl_inter
          ! Drag
          drag_cp = 0.5d0*density*vrel_cp**2*chord*cd_inter

          !----Forces----
          ! These forces lead to the source terms and are defined as
          ! spanwise, tangential and normal to the blade.          
          Fstn(:,1) = 0.0d0
          Fstn(:,2) = lift_cp*sin(phi_cp)-drag_cp*cos(phi_cp)
          Fstn(:,3) = lift_cp*cos(phi_cp)+drag_cp*sin(phi_cp)

          cpu_t2 = rtc() ! Stop timing
          !=================== Output ==================================
          write(150+(j-1)*n_rotors+i,'(es15.7,i5,*(es15.7))')
     &      time,nw_sub_iter,nw_resmax,alpha_cp_cfd,alpha_cp,phi_vtp,
     &      nw_vt,nw_vn,dGamma_vtp,(cpu_t2-cpu_t1)

          !=================== Saving iteration ========================
          nw_W_history(i,j,:) = nw_W
          nw_dyn_comp_history(i,j,:,:) = nw_dyn_comp 
          nw_Gamma_history(i,j,:) = Gamma_cp_steady
          nw_vn_history(i,j,:) = nw_vn
          nw_vt_history(i,j,:) = nw_vt
          nw_X_history(i,j,:,:) = nw_X
          nw_Y_history(i,j,:,:) = nw_Y
          !=============================================================
          nw_converged = .false.
       
      end subroutine nw_smearing_correction
c=======================================================================
c=======================================================================
      subroutine smearing_factor(r,h_r,epsilon,smear_fac)   
c-----------------------------------------------------------------------
c     The Gaussian smearing results in a Lamb-Oseen like vortex. The 
c     smearing in the induced velocity can directly be determined from 
c     their analytical solution. 
c     Definition and details in Meyer Forsting et al. (2019).
c-----------------------------------------------------------------------
      use params,only:idp ! Global precision
c-----------------------------------------------------------------------      
      implicit none    
c-----[in]        
      real(kind=idp),intent(in)::r,h_r,epsilon
c-----[out]      
      real(kind=idp),intent(out)::smear_fac 
c-----[tmp]      
      real(kind=idp)::dist_orth
c=======================================================================      
      ! The definition of the smearing factor is given in Eq.(22) with 
      ! the perpendicular distance between vortex element and control 
      ! point defined in Eq.(21). Setting beta to zero, the expression 
      ! is greatly simplified. 
      ! Distance between CP and vortex trailing point
      dist_orth = (r*h_r)
      ! The smearing factor
      smear_fac = exp(-dist_orth**2/epsilon**2)

      end subroutine smearing_factor
c=======================================================================
c=======================================================================
      subroutine nw_hr_components(h,h_r,
     & 	                          nw_phi,nw_phi_s,nw_dx,nw_dy,nw_a)
c-----------------------------------------------------------------------      
c     Compute the components of the near-wake model that only rely on
c     geometric definiton of the blade. Unless the blade deforms heavily
c     these components remain unchaged. Definitions of all equations and
c     variables follow from Pirrung et al. (2016,2017b).
c-----------------------------------------------------------------------      
      use params,only:idp ! Global precision
      use nw_smearing_correction_module,only:pi,nw_P,nw_N ! Constants
c----------------------------------------------------------------------- 
      implicit none
c-----[in] Distance between CP and vortex trailing point 
      real(kind=idp),intent(in)::h_r,h
c-----[out]      
      real(kind=idp),intent(out)::nw_phi,nw_phi_s,nw_a(4),nw_dx,nw_dy
c-----[tmp]      
      real(kind=idp)::dw,root_corr
      integer::k
c=======================================================================
      ! Pirrung et al. (2016)
      !------ Wang Coton -----------------
      ! Determine the model phi following Wang and Coton 
      ! Eq.(4), note in the paper the log is erroneously missing
      if (h_r.gt.0d0) then
         nw_phi = pi/4.0d0*abs((1d0+0.5d0*h_r)*log(1d0-h_r))
      else
         nw_phi = log(1d0-h_r)/(1.5d0+log(1-0.5d0*h_r))
      endif
      !------- Amplitude term -----------------
      ! The start value for the indicial function. Split into slow and 
      ! fast decaying parts. Gamma is added later. 
      ! Eq.(25) 
      dw = 1d0/(4d0*pi*abs(h)*h_r)
      nw_dx = 1.359d0*dw
      nw_dy = 0.359d0/4d0*dw
      !-------- Root Correction ----------
      ! Implement the root correction that limits the induction to that 
      ! 90 degrees
      ! Eq.(27)
      root_corr=(nw_dx*(1d0-exp(-pi/(2d0*nw_phi)))-nw_dy* 
     &          (1d0-exp(-2d0*pi/nw_phi)))/(nw_dx-nw_dy)      
      ! Apply root correction 
      nw_dx = root_corr*nw_dx
      nw_dy = root_corr*nw_dy
      nw_phi = root_corr*nw_phi

      ! Pirrung et al. (2017b)
      !----------- Convection correction -----------
      ! The wake can be in-plane solely convecting downstream and 
      ! anything in-between. For this purpose the phi needs to be 
      ! corrected for the convection of the vortex helix.
      ! Phi for straight vortices
      ! Eq.(7)
      if (h_r.lt.0d0) then
         nw_phi_s = -0.788d0*h_r
      else
         nw_phi_s = 0.788d0*h_r
      endif 
      ! Determine the mixing ratio of straight and helical vorticity
      ! First compute the constants 
      ! Eq.(10) and (12)
      if (h_r.lt.0d0) then
         do k=1,4
            nw_a(k)=nw_N(k,1)+nw_N(k,2)*exp(nw_N(k,3)*h_r)
     &                +nw_N(k,4)*exp(nw_N(k,5)*h_r)
     &                -nw_N(k,4)-nw_N(k,2)
         enddo
      else
         do k=1,4
            nw_a(k)=nw_P(k,1)+nw_P(k,2)*h_r+nw_P(k,3)*h_r**2
     &               +nw_P(k,4)*h_r**3
         enddo
      endif

      end subroutine nw_hr_components
c=======================================================================
c=======================================================================
      subroutine nw_phi_components(h_r,nw_a,phi,nw_phi,nw_phi_s,
     &                             nw_phi_star)
c-----------------------------------------------------------------------      
c     Components of the near-wake model related to the inflow angle
c     phi. 
c-----------------------------------------------------------------------  
      use params,only:idp
      use nw_smearing_correction_module,only:pi
c-----------------------------------------------------------------------        
      implicit none
c-----[in]       
      real(kind=idp),intent(in)::h_r,nw_a(4),phi,nw_phi,nw_phi_s
c-----[out]     
      real(kind=idp),intent(out)::nw_phi_star
c-----[tmp]      
      real(kind=idp)::nw_k_phi
c=======================================================================  
      ! Pirrung et al. (2017b)     
      !----Mixing factor----
      ! Eq. (9) and (11)
      if (h_r.lt.0d0) then
         nw_k_phi = nw_a(1)+nw_a(2)*exp(nw_a(3)*(pi/2d0-phi))+
     &     nw_a(4)*exp(-8d0*(pi/2d0-phi))-nw_a(2)-nw_a(4)
      else
         nw_k_phi = nw_a(1)+nw_a(2)*phi+nw_a(3)*phi**2+nw_a(4)*phi**3
      endif
      ! The factor is only valid for angles below 89 degrees
      nw_k_phi = min(nw_k_phi,1d0)
      nw_k_phi = max(nw_k_phi,0d0)
      ! Determine the mix of straight and curved helix
      ! Eq. (8)
      nw_phi_star = nw_k_phi*nw_phi_s + (1d0-nw_k_phi)*nw_phi

      end subroutine nw_phi_components
c=======================================================================
c=======================================================================
      subroutine nw_XY_components(nw_phi_star,nw_dx,nw_dy,dbeta,dGamma,
     &                            nw_X_history,nw_Y_history,
     &                            nw_X,nw_Y)
c-----------------------------------------------------------------------      
c     Induced velocity components of the near-wake model, with iterative
c     dependancy through dbeta and dGamma. 
c-----------------------------------------------------------------------  
      use params,only:idp
c-----------------------------------------------------------------------        
      implicit none
c-----[in]       
      real(kind=idp),intent(in)::nw_phi_star,nw_dx,nw_dy,dbeta,dGamma,
     &                           nw_X_history,nw_Y_history  
c-----[out]     
      real(kind=idp),intent(out)::nw_X,nw_Y
c-----[tmp]      
      real(kind=idp)::nw_X_exp,nw_Y_exp,nw_X_new,nw_Y_new,
     &                nw_X_old,nw_Y_old
c=======================================================================  
      !----Induced velocity components----
      ! Determine the slow and fast decaying parts 
      ! Pirrung et al. (2016) Eq.(25) with modified phi as described in 
      ! Pirrung et al. (2017b)
      
      ! Exponential terms
      nw_X_exp=exp(-dbeta/nw_phi_star)
      nw_Y_exp=exp(-4d0*dbeta/nw_phi_star)
      ! New induction from shed vorticity
      nw_X_new = dGamma*nw_dx*nw_phi_star*(1d0-nw_X_exp)
      nw_Y_new = dGamma*nw_dy*nw_phi_star*(1d0-nw_Y_exp)
      ! Contribution from previously shed vortex elements,
      ! which need to be advanced in time
      nw_X_old = nw_X_history*nw_X_exp
      nw_Y_old = nw_Y_history*nw_Y_exp
      ! Total induction from vortex line, fast and slow 
      nw_X = nw_X_new + nw_X_old
      nw_Y = nw_Y_new + nw_Y_old
      end subroutine nw_XY_components
c=======================================================================
c=======================================================================
      subroutine nw_circulation(dt_global,vrel,alpha,Gamma_history,
     &                          dyn_comp_history,
     &                          alpha_2d,cl_2d,len_2d,chord,     
     &                          Gamma_steady,Gamma,dGamma,dyn_comp)

c-----------------------------------------------------------------------      
c     Determine the circulation along the lifting line, including the 
c     dynamic effects and also returns the delta in circulation. 
c----------------------------------------------------------------------- 
      use params,only:idp ! Global precision
      use nw_smearing_correction_module,only:n_s,n_v ! Constants
c-----------------------------------------------------------------------
      implicit none
c-----[in]  
      integer,intent(in)::len_2d
      real(kind=idp),intent(in)::dt_global,vrel(n_s),alpha(n_s),
     &               Gamma_history(n_s),dyn_comp_history(n_s,3),
     &               alpha_2d(len_2d,n_s),cl_2d(len_2d,n_s),chord(n_s)      
c-----[out]  
      real(kind=idp),intent(out)::Gamma_steady(n_s),Gamma(n_s),
     &                            dGamma(n_v),dyn_comp(n_s,3)
c-----[tmp]        
      integer::l
      real(kind=idp),dimension(n_s)::cl_inter,dt,T0
c=======================================================================

      !----Circulation at sections----
      ! Interpolate the lift coefficients, output is cl_inter
      do l = 1,n_s
        call nw_interp1d(alpha_2d(:,l),cl_2d(:,l),len_2d,alpha(l),
     &                                        cl_inter(l),1)
      enddo
      ! Calculate the steady circulation at each blade element
      Gamma_steady = 0.5d0*vrel*chord*cl_inter
      
      !----Time filter bound circulation----
      ! Madsen and Gaunaa (2004)
      T0=chord/(2d0*vrel)
      dt=dt_global/T0
      ! Filter functions by Mac Gaunaa
      dyn_comp(:,1)=dyn_comp_history(:,1)
     &         *exp(-0.3064d0*dt)+0.27735d0*(Gamma_steady
     &         +Gamma_history)*(1d0-exp(-0.3064d0*dt))

      dyn_comp(:,2)=dyn_comp_history(:,2)
     &         *exp(-0.0439d0*dt)+0.0914d0*(Gamma_steady
     &         +Gamma_history)*(1d0-exp(-0.0439d0*dt))
      dyn_comp(:,3)=dyn_comp_history(:,3)
     &         *exp(-3.227d0*dt)+0.13125d0*(Gamma_steady
     &         +Gamma_history)*(1d0-exp(-3.227d0*dt))
      ! Dynamic filtered circulation
      Gamma = dyn_comp(:,1)+dyn_comp(:,2)+dyn_comp(:,3)
      ! Change in Gamma between sections.At blade ends take all Gamma
      ! from neighboring point. 
      dGamma(2:n_v-1) = Gamma(2:n_s)-Gamma(1:n_s-1)
      dGamma(1) = Gamma(1)
      dGamma(n_v)= -Gamma(n_s);  

      end subroutine nw_circulation
c=======================================================================
c=======================================================================
      subroutine nw_calc_relax(vrel,alpha,vn,vt,chord,dt_global,
     &                         nw_dx,nw_dy,dbeta,nw_phi,
     &                         nw_relax_safety,nw_relax)
c-----------------------------------------------------------------------      
c     Determine the relaxation factor needed in the near wake model. 
c     The defintions of Pirrung et al. (2017a) are used and the 
c     corresponding equations referenced.
c-----------------------------------------------------------------------       
      use params,only:idp ! Global precision
      use nw_smearing_correction_module,only:pi,n_s,n_v ! Constants
c-----------------------------------------------------------------------        
      implicit none 
c-----[in]  
      real(kind=idp),intent(in)::vrel(n_s),alpha(n_s),
     &                           vt(n_s),
     &                  dt_global,nw_dx(n_v,n_s),
     &                   nw_dy(n_v,n_s),
     &                  dbeta(n_v),nw_phi(n_v,n_s),
     %                  nw_relax_safety,chord(n_s),
     &                  vn(n_s)
c-----[out]  
      real(kind=idp),intent(out)::nw_relax
c-----[tmp]     
      real(kind=idp)::A1,A2,gradw,tau(n_s),d(n_s),
     &                 nw_X(n_v,n_s),nw_Y(n_v,n_s)
      integer::l,k
c=======================================================================
      !----Initialize----   
      nw_relax = 0d0

      !----Dynamic to steady circulation ratio----
      ! Eq. (15)
      tau=(2d0*vrel)/chord
      d=1-0.5547d0/(tau*dt_global*0.3064d0)*
     &             (1-exp(-tau*dt_global*0.3064d0))-0.1828d0/
     &             (tau*dt_global*0.0439d0)*(1-exp(-tau*dt_global
     &              *0.0439d0))-0.2625d0/(tau*dt_global*3.2277d0)
     &                          *(1-exp(-tau*dt_global*3.2277d0))
      

      !----Slow and fast decaying components---- 
      ! in-plane only (no phi_star)
      do l = 1,n_s
        do k = 1,n_v
          nw_X(k,l) = nw_dx(k,l)*nw_phi(k,l)
     &                             *(1d0-exp(-dbeta(k)/nw_phi(k,l)))
          nw_Y(k,l) = nw_dy(k,l)*nw_phi(k,l)
     &                          *(1d0-exp(-4d0*dbeta(k)/nw_phi(k,l)))
        enddo
      enddo

      !----Time derivative of induced velocities----
      ! Eq. (27)
      do l = 1,n_v-1
        A1 = nw_X(l,l)+nw_Y(l,l)
        A2 = nw_X(l+1,l)+nw_Y(l+1,l)
        gradw = d(l)*pi*chord(l)*(A1-A2)*
     &                   ( (alpha(l)*vn(l))/vrel(l)+
     &            vrel(l)/(abs(vt(l))*(vn(l)/abs(vt(l)))**2+1d0) ) 
        ! Update the relaxation factor
        nw_relax = max(nw_relax,-(1d0+gradw)/(1d0-gradw))
      enddo

      !----Mix in the safety factor---- 
      nw_relax=nw_relax+nw_relax_safety*(1.0d0-nw_relax)
      nw_relax=min(nw_relax, 0.99d0)
      
      end subroutine nw_calc_relax
c=======================================================================
c=======================================================================
      subroutine init_nw_smearing_correction(n_rotors_in, n_blades_in,
     &                                       n_sections_in)
c-----------------------------------------------------------------------      
c     Initialize the arrays needed in the near-wake smearing corection. 
c-----------------------------------------------------------------------  
      use nw_smearing_correction_module ! Assign all 
      use TurbineData,only:ACLproject,aclnamelen ! Project name 
      use params,only:MyProcNum ! Get the processor nr. 
c-----------------------------------------------------------------------
      implicit none   
c-----[in]  
      integer,intent(in)::n_rotors_in,n_blades_in,n_sections_in 
c-----[tmp]              
      character(len=256)::filename
      character(len=30)::nw_crotor,cnblade     
      integer::i,j   
c=======================================================================      
      if(.not.kaya)return ! Switch 
      !----Copy for module----
      n_rotors = n_rotors_in
      n_blades = n_blades_in
      !----Nr of control and vortex trailing points----
      n_s = n_sections_in
      n_v = n_s+1

      !----Near-wake constants Pirrung et al. (2017b)----
      ! Positve values of h_r given in Eq.(14)
      nw_P(1,:)=(/-1.64636754184988d0,8.14821474772595d0, 
     &   -12.1784861715161d0,5.02652654794305d0,21.7713111885100d0/)
      nw_P(2,:)=(/-0.499006672527007d0,6.08465268675599d0,
     &   -15.1712005446736d0,14.8254135311781d0,-2.42318681495122d0/)
      nw_P(3,:)=(/3.90835606394205d0,-18.7662330478906d0,
     &    39.1243282362892d0,-29.4870086106427d0,-60.2947309276912d0/)
      nw_P(4,:)=(/-1.60622843516794d0,7.42952711674470d0,
     &    -15.8594777074799d0,11.6870163960092d0,-195.060865844227d0/)
      ! Negative values given in Eq.(13)
      nw_N(1,:)=(/1.01933206359188d0,-0.135668151345576d0,
     &    0.395524705910666d0,0.0801822997046969d0,44.8347503471217d0/)
      nw_N(2,:)=(/12.9874512870379d0,49.9999999999916d0, 
     &    0.00235345218088083d0,11.3116074964299d0,3935.34323353354d0/)
      nw_N(3,:)=(/-0.690159136848837d0,101.238775997879d0,
     &  -0.00154467589177718d0,3.99520376143011d0,0.394541495293350d0/)
      nw_N(4,:)=(/-0.269253409815015d0,49.9999999999843d0,
     &  -0.00247975579295973d0,0.403642518874381d0,1.16610462658971d0/)

      !================= Allocation ====================================

      !----Near-wake variables---- 
      allocate(nw_X_history(n_rotors,n_blades,n_v,n_s))
      nw_X_history =0d0
      allocate(nw_Y_history(n_rotors,n_blades,n_v,n_s))
      nw_Y_history =0d0
      allocate(nw_vn_history(n_rotors,n_blades,n_s))
      nw_vn_history =0d0
      allocate(nw_vt_history(n_rotors,n_blades,n_s))
      nw_vt_history =0d0
      allocate(nw_Gamma_history(n_rotors,n_blades,n_s))
      nw_Gamma_history =0d0
      allocate(nw_dyn_comp_history(n_rotors,n_blades,n_s,3))
      nw_dyn_comp_history =0d0
      allocate(nw_W_history(n_rotors,n_blades,n_s))
      nw_W_history =0d0
      allocate(nw_phi(n_rotors,n_blades,n_v,n_s))
      nw_phi =0d0
      allocate(nw_phi_s(n_rotors,n_blades,n_v,n_s))
      nw_phi_s =0d0
      allocate(nw_dx(n_rotors,n_blades,n_v,n_s))
      nw_dx =0d0
      allocate(nw_dy(n_rotors,n_blades,n_v,n_s))      
      nw_dy =0d0
      allocate(nw_a(n_rotors,n_blades,n_v,n_s,4)) 
      nw_a =0d0

      !----Smearing correction variables----
      allocate(cp_loop_ind(n_rotors,n_blades,n_v,2))
      ! Initialize with first and last index of control point
      cp_loop_ind(:,:,:,1) = 1
      cp_loop_ind(:,:,:,2) = n_s
      allocate(smear_fac(n_rotors,n_blades,n_v,n_s))
      smear_fac=0d0
      allocate(kaya_refresh(n_rotors,n_blades))
      kaya_refresh=.true.        

      !----Initialize the output----
      if(MyProcNum.eq.1)then
        filename(1:aclnamelen)=ACLproject(1:aclnamelen)
        do i= 1,n_rotors
          ! Create NW part of filename
          if (i<10)then
            nw_crotor(1:4)='.nw0'; write(nw_crotor(5:5),'(i1)')i
          else
            nw_crotor(1:3)='.nw'; write(nw_crotor(4:5),'(i2)')i
          endif
          ! Create file for each blade
          do j = 1,n_blades
            ! Blade name
            cnblade(1:4)='.bl0'; write(cnblade(5:5),'(i1)')j
            ! Compose name
            filename(aclnamelen+1:aclnamelen+5)=nw_crotor(1:5)
            filename(aclnamelen+6:aclnamelen+10)=cnblade(1:5)
            ! Open file 
            open(unit=150+(j-1)*n_rotors+i,
     &            file=filename(1:aclnamelen+10),form='formatted')
            rewind(150+(j-1)*n_rotors+i)
          enddo
        enddo
      endif

      end subroutine init_nw_smearing_correction
c=======================================================================
c=======================================================================
      subroutine nw_interp1d(x,y,n,xin,yin,nin)
c-----------------------------------------------------------------------      
c     One-dimensional linear interpolation
c-----------------------------------------------------------------------      
      use params,only:idp
c-----------------------------------------------------------------------      
      implicit none
c-----[in]  
      integer,intent(in)::n,nin
      real(kind=idp),intent(in)::x(n),y(n),xin(nin)
c-----[out]  
      real(kind=idp),intent(out)::yin(nin)
c-----[in]        
      integer::i,j
      real(kind=idp)::ratio
c=======================================================================
      do j=1,nin
c------  extrapolate
         if(xin(j).gt.x(n))then
           ratio=(xin(j)-x(n))/(x(n)-x(n-1))
           yin(j)=y(n)+ratio*(y(n)-y(n-1))
         endif
c------  extrapolate
         if(xin(j).lt.x(1))then
           ratio=(xin(j)-x(1))/(x(2)-x(1))
           yin(j)=y(1)+ratio*(y(2)-y(1))
         endif
c------  interpolate
         do i=2,n
            if(x(i).ge.xin(j).and.x(i-1).le.xin(j))then
               ratio=(xin(j)-x(i-1))/(x(i)-x(i-1))
               yin(j)=y(i-1)+ratio*(y(i)-y(i-1))
               exit
            endif
         enddo
      enddo

      end subroutine nw_interp1d
c=======================================================================