 
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

   
end module PARTICLE_TINIA


