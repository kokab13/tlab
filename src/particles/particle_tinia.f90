 
module PARTICLE_TINIA
   use TLAB_CONSTANTS, only: MAX_VARS, wp, wi, sp 
   use AVG_SPATIAL
   use PARTICLE_TYPES, only: particle_dt
   use PARTICLE_VARS
   use PARTICLE_ARRAYS
   use TLAB_VARS, only: visc
   use TIME

   implicit none
   save
   private

   public :: PARTICLE_TINIA_READBLOCK
   public :: PARTICLE_TINIA_INITIALIZE
   public :: PARTICLE_TINIA_RHS_1

   ! ###################################################################
   ! Carrier Phase
   ! ###################################################################
   real(wp) :: K_f                        ! thermal conductivity of flow 
   real(wp) :: mu_f                       ! dynamic viscosity
   real(wp) :: rho_f                      ! density
   real(wp) :: dif_v                      ! diffusion coefficient of water vapor
   real(wp) :: ui_i                       ! velocity of vapor molecules
   real(wp) :: lambda_g_i                 ! mean free path
   real(wp) :: Kn_mt_i                    ! Knudsen number
   real(wp) :: rhoi_i                     ! vapor mass density
   real(wp) :: T_0                        ! ambient temperature
   real(wp) :: p_sat_i_liq                ! saturation vapor pressure
   real(wp) :: p_sat_i_ice                ! saturation vapor pressure   
   real(wp) :: di_mic                     ! diffusion coefficient of water vapor
   real(wp) :: cp_air                     ! heat coefficient of air
   real(wp) :: LvRv_l                     ! the ratio of latent heat to the gas constant
   real(wp) :: LvRv_ice                   ! the ratio of latent heat to the gas constant
   real(wp) :: R_G                        ! general gas constant
   real(wp) :: Rv                         ! indivudual gass constant for water vapor
   real(wp) :: rho_air                    ! density of air
   real(wp) :: Rd                         ! Rd is the dry-air gas constant 
   real(wp) :: pi
   ! ###################################################################
   ! Descrete Phase
   ! ###################################################################
   real(wp) :: dummyT_i                   ! 1/tau_T
   real(wp) :: dummy_mp                   ! 1/m_p
   real(wp) :: d_p                        ! particle diameter
   real(wp) :: T_p                        ! particle temperature
   real(wp) :: rho_l                      ! luiquid droplet density
   real(wp) :: rho_ice                    ! ice density
   real(wp) :: Cp_p                       ! particle heat coefficient
   real(wp) :: velMag_i                   ! magnitude of velocity
   real(wp) :: Re_p_i                     ! Re number
   real(wp) :: CD_i                       ! drag coefficient
   real(wp) :: mp_i                       ! particle mass
   real(wp) :: FD_x_i                     ! drag in x-direction
   real(wp) :: FD_y_i                     ! drag in y-direction
   real(wp) :: FD_z_i                     ! drag in z-direction
   real(wp) :: tau_T_i                    ! thermal response time
   real(wp) :: A_factor                   ! factor A, the following factors come from D. Niedermeier et al., 2020.
   real(wp) :: B1_factor                  ! factor B1
   real(wp) :: B2_factor                  ! factor B2
   real(wp) :: C_factor                   ! factor C
   real(wp) :: f_mt_i                     ! correction factor
   real(wp) :: p_0                        ! saturation vapor pressure
   real(wp) :: w_mi                       ! molar weight of water droplet
   real(wp) :: s_ten                      ! surface tension of water
   real(wp) :: nuPhi                      ! number of ions in NaCl with zero molality
   real(wp) :: rho_s                      ! mass density of NaCl
   real(wp) :: w_ms                       ! molar mass of NaCl
   real(wp) :: d_mo                       ! diameter of NaCl
   real(wp) :: Satu_star_i                ! köhler
   real(wp) :: Satu_p                     ! saturation ratio
   real(wp) :: vortMag_i                  ! magnitute of vorticity
   real(wp) :: Re_s_i                     ! shear Re number
   real(wp) :: bet_i                      ! ß
   real(wp) :: f_sl_i                     ! correction factor
   real(wp) :: CL_i                       ! lift coefficient
   real(wp) :: FL_x_i                     ! lift in x-direction
   real(wp) :: FL_y_i                     ! lift in y-direction
   real(wp) :: FL_z_i                     ! lift in z-direction
   real(wp) :: J_het_i                    ! Het. nucleation rate
   real(wp) :: np_liq                     ! number of liquid particles
   real(wp) :: np_ice                     ! number of ice particles
   real(wp) :: P_freeze                   ! probability of freezing
   real(wp) :: r                          ! random number
   real(wp) :: rad                        ! radius of each particle 
   real(wp) :: S_ice                      ! supersaturation over ice
   real(wp) :: alpha_c                    ! coming from Grabowski et al., 2011
   real(wp) :: alpha_T                    ! coming from Grabowski et al., 2011       
   real(wp) :: delta_v                    ! coming from Grabowski et al., 2011
   real(wp) :: delta_T                    ! coming from Grabowski et al., 2011
   real(wp) :: lambda_a                   ! Pruppacher and Klett, 1996, table 13-1
   real(wp) :: Ls                         ! latent heat required to deposit per unit mass of water vapor to ice  
   real(wp) :: Lv                         ! for liquid water the latent heat of vaporization                     
   real(wp) :: dp_i                       ! particle diameter 
   real(wp) :: k_prime                    ! coming from S. Chen et al., 2023
   real(wp) :: D_prime                    ! coming from S. Chen et al., 2023
   real(wp) :: dummyT0                    ! for mathematical simplicity
   real(wp) :: A_                         ! for mathematical simplicity
   real(wp) :: B_                         ! for mathematical simplicity
   real(wp) :: C_                         ! for mathematical simplicity
   real(wp) :: D_                         ! for mathematical simplicity
   real(wp) :: E_                         ! for mathematical simplicity
   real(wp) :: F_                         ! for mathematical simplicity
   real(wp) :: G_                         ! for mathematical simplicity
   real(wp) :: K_                         ! for mathematical simplicity
   real(wp) :: AA                         ! for mathematical simplicity
   real(wp) :: BB                         ! for mathematical simplicity
   real(wp) :: CC                         ! for mathematical simplicity
   real(wp) :: DD                         ! for mathematical simplicity
   real(wp) :: EE                         ! for mathematical simplicity
   real(wp) :: ffff                       ! for mathematical simplicity
   real(wp) :: e_v_d                      ! partial presure of the dry air 
   real(wp) :: mix                        ! water vapor mixing ratio
   real(wp) :: dt_lag
   real(wp) :: relaxation_time

   public :: np_ice
   contains
   
   !########################################################################
   !########################################################################

     subroutine PARTICLE_TINIA_READBLOCK(bakfile, inifile, block)
       use PROFILES
       use TLab_WorkFlow, only: TLab_Write_ASCII        
       character(len=*), intent(in) :: bakfile, inifile, block
       

       call TLab_Write_ASCII(bakfile, '#Diameter=<value>')
       call TLab_Write_ASCII(bakfile, '#ParDensity=<value>')
       call TLab_Write_ASCII(bakfile, '#ParHeatCoeff=<value>')
       call TLab_Write_ASCII(bakfile, '#mixRatio=<value>')
       call TLab_Write_ASCII(bakfile, '#Conductivity=<value>')

       call ScanFile_Real(bakfile, inifile, block, 'Diameter', '0.0', d_p)
       call ScanFile_Real(bakfile, inifile, block, 'TemperatureP', '1.0', T_p)
       call ScanFile_Real(bakfile, inifile, block, 'DropletDensity', '0.0', rho_l)
       call ScanFile_Real(bakfile, inifile, block, 'IceDensity', '0.0', rho_ice)
       call ScanFile_Real(bakfile, inifile, block, 'ParHeatCoeff', '0.0', Cp_p)
       call ScanFile_Real(bakfile, inifile, block, 'mixRatio', '0.0', mix) 
       call ScanFile_Real(bakfile, inifile, block, 'Conductivity', '0.0', K_f)   

       return
   end subroutine PARTICLE_TINIA_READBLOCK

   !########################################################################
   !########################################################################
   subroutine PARTICLE_TINIA_INITIALIZE()  
   implicit none

   integer(wi) :: j


     do j = 1, l_g%np

             l_q(j, 4) = 0.0_wp
             l_q(j, 5) = 0.0_wp
             l_q(j, 6) = 0.0_wp
             l_q(j,7) = d_p**2_wp
             l_q(j,8) = T_p
             l_q(j,9) = 0.0_wp 
             l_q(j,10) = 0.0_wp 
             l_q(j,12) = mix
             l_q(j,11) = 0.0_wp 

     
     end do
      

   return

   end subroutine PARTICLE_TINIA_INITIALIZE
   

   !########################################################################
   !########################################################################   
   
   subroutine PARTICLE_TINIA_RHS_1(l_hq)
        
      implicit none
       real(wp), intent(inout) :: l_hq(:,:)
       integer(wi) i
       !########################################################################
                   
       mu_f = 2.0e-5              ! dynamic viscosity, [kg/m.s], 1.81e-5
       p_0 = 611.3                ! vapor pressure, [Pa]
       LvRv_l = 5423.0            ! the ratio of latent heat to the gas constant, [K] → Lv = 2.5e06 [J/kg], Rv = 461 [J/K·kg]
       LvRv_ice = 6139.0          ! the ratio of latent heat to the gas constant, [K]
       T_0 = 273.15               ! [K]
       dummyT0 = 1.0_wp/273.15    ! [1/K]  
       w_mi = 0.018               ! molar weight (here, it is water), [kg/mol]
       R_G = 8.3145               ! general gas constant, [J/K.mol]
       dif_v = 0.242e-4           ! diffusion coefficient of water vapor, [m²/s]
       A_factor = 1.0             ! the following factors come from D. Niedermeier et al., 2020.
       B1_factor = 0.377
       B2_factor = 4.0/3.0
       C_factor = 4.0/3.0                    
       s_ten = 0.07288            ! surface tension of water, [N/m]
       nuPhi = 2.0                ! number of ions in NaCl with zero molality
       rho_s = 2163.0             ! mass density of NaCl, [kg/m³]
       w_ms = 0.05844             ! molar mass of NaCl, [kg/mol]
       d_mo = 185.0e-9            ! diameter of NaCl, [m]
       di_mic = dif_v * 1.0e12_wp ! diffusion coefficient of water vapor, [µm²/s]
       pi = 3.14159!2653589793238
       Ls = 2.83e6                ! [J/kg], latent heat required to deposit per unit mass of water vapor to ice
       Lv = 2.5e6                 ! [J/kg], latent heat required to deposit per unit mass of water vapor to ice
       Rv = 467.0                 ! [J/K.kg], individual gas constant for water vapor
       lambda_a = 8e-8            ! [m]
       delta_v = 1.3 * lambda_a   ! [m]
       delta_T = 2.16e-7          ! [m]
       alpha_c = 0.036
       alpha_T = 0.7
       rho_air = 1.293            ! [kg/m³], density of air
       Rd = 287.052874            ! [J/K.kg], individual gas constant for dry air
       cp_air = 1000.0            ! [J/kg.K]
       !########################################################################

      do i = 1, l_g%np
         
         
         dp_i = sqrt(l_q(i,7))*1e-6 ![m]
                  
          l_hq(i,9) = l_hq(i,9) + 0.0                                                                                 !keep the droplet always liquid
          
          velMag_i = sqrt ((l_txc(i,1)-l_q(i,4))**2 + (l_txc(i,2)-l_q(i,5))**2 + (l_txc(i,3)-l_q(i,6))**2 + 1e-34)    !magnitude of velocity, [m/s]
          
          Re_p_i = dp_i * velMag_i / visc                                                                             !particle Re
          
          IF ( Re_p_i < 0.5) THEN
             CD_i = 24.0 / Re_p_i                                                                                     !particle drag coefficient                                              
          ELSE
             CD_i = (24.0 / Re_p_i) * (1 + (0.15 * Re_p_i**0.687))
          ENDIF   
          
          mp_i = (1.0/6.0) * pi * (dp_i)**3 * rho_l                                                                   !particle mass [kg]

          FD_x_i = (3.0/4.0) * (rho_air/rho_l) * (mp_i/dp_i) * CD_i * (l_txc(i,1)-l_q(i,4)) * velMag_i                !particle drag force [N]
          FD_y_i = (3.0/4.0) * (rho_air/rho_l) * (mp_i/dp_i) * CD_i * (l_txc(i,2)-l_q(i,5)) * velMag_i
          FD_z_i = (3.0/4.0) * (rho_air/rho_l) * (mp_i/dp_i) * CD_i * (l_txc(i,3)-l_q(i,6)) * velMag_i
          
          dummy_mp = 1.0_wp / mp_i 

          ! equation dv_p/dt = 1/responseTime (u-v) in x, y, z 
          l_hq(i,4) =  l_hq(i,4) + dummy_mp * FD_x_i 
          
          l_hq(i,5) =  l_hq(i,5) + dummy_mp * FD_y_i - 9.8                                                             !gravity included
          
          l_hq(i,6) =  l_hq(i,6) + dummy_mp * FD_z_i 
          
          ! equation dx_p/dt = v_p in x, y, z 
          l_hq(i,1) = l_hq(i,1) + l_q(i,4)
          l_hq(i,2) = l_hq(i,2) + l_q(i,5)
          l_hq(i,3) = l_hq(i,3) + l_q(i,6) 
                    
          l_q(i,8) = l_txc(i,4)                                                                                       !temperature [K]
           
          l_q(i,12) = l_txc(i,5)                                                                                      !water vapor [g/kg]
          
          e_v_d = (l_q(i,12)*101325.0/1000.0) / (0.622 + (l_q(i,12)/1000.0))                                          !partial vapor pressure [Pa]
                          
          p_sat_i_liq = p_0 * exp(LvRv_l * (dummyT0 - (1/l_q(i,8))))                                                  !saturation vapor pressure → Clausius-Clapeyron, wrt liquid [Pa]
          
          Satu_p = e_v_d / p_sat_i_liq                                                                                !saturation ratio wrt liquid

 
             k_prime = (1.5e-11 * l_q(i,8)**3) - (4.8e-8 * l_q(i,8)**2) + (1.0e-4 * l_q(i,8)) - (3.9e-4)                 !thermal conductivity [W/mK]        
             
             D_prime = ((0.015 * l_q(i,8)) - 1.9) * 1.0e-5                                                               !water vapor diffusivity [m²/s]  
                  
             ! modification based on Pruppacher and Klett 1997, ch. 13
             rad = 0.5 * dp_i                                                                                            !radius of particle [m]
             A_ = ( rad ) / ( rad + delta_v )
             B_ = D_prime / ( rad * alpha_c )
             C_ =  (2*pi)/(Rv * l_q(i,8)) 
             D_ = D_prime / ( A_ + ( B_ * sqrt(C_) ) )                                                                   !water vapor diffusivity including kinetic effects [m²/s]        
             E_ = ( rad ) / ( rad + delta_T)
             F_ = k_prime / ( rad * alpha_T * rho_air * cp_air)
             G_ = sqrt( (2*pi)/(Rd * l_q(i,8)) )
             K_ = k_prime / ( E_ + (F_*G_) )                                                                             !thermal conductivity including kinetic effects [W/mK] 
                      
             ! equation of growth «The output of this equation is d², i.e., the unit is [µm²]» 
             EE = Lv/( Rv * l_q(i,8) )
             BB = EE-1.0
             CC = Lv/( K_ * l_q(i,8) )
             DD = ( Rv * l_q(i,8) )/( p_sat_i_liq * D_ )

                                 
          IF ( sqrt(l_q(i,7)) > d_mo ) THEN                                                                           !equilibrium saturation
             Satu_star_i = exp(((4.0 * w_mi * s_ten)/(R_G * rho_l * l_q(i,8) * (dp_i)))-&
             (((nuPhi * rho_s * w_mi)/(w_ms * rho_l))*((d_mo**3)/((dp_i)**3-d_mo**3 ))))
           ELSE
             Satu_star_i = 0.0
           ENDIF 
                            
             
          l_hq(i,7) = l_hq(i,7) +  ( ((8.0_wp/rho_l) * (1.0_wp/( (BB * CC) + DD )) * (Satu_p - Satu_star_i) ) ) * 1.0e12
    

         
      end  do 
    
   return
   end subroutine PARTICLE_TINIA_RHS_1   

   !########################################################################
   !######################################################################## 


end module PARTICLE_TINIA


