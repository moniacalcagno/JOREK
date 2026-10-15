!> Complete the variables that were not in the restart file when importing it with var_order_import.
!!
!! Variables with an input profile (rho, T, T_i, T_e, rho_n, Vpar) get the n=0 profile, as in
!! initial_conditions; all other missing variables stay zero. If the total temperature T of the file
!! was stored in T_e (split_T), T_i is taken from its input profile and T_e = T - T_i for n=0, while
!! the n/=0 harmonics of T are split equally between T_i and T_e; corr_neg_temp is applied to the
!! n=0 part of T_e.
subroutine init_missing_variables(node_list, var_imported, split_T)

use constants
use data_structure
use phys_module
use equil_info
use mod_sources
use corr_neg, only: corr_neg_temp_node

implicit none

type (type_node_list), intent(inout) :: node_list
logical,               intent(in)    :: var_imported(n_var) !< .true. if the variable was read from the restart file
logical,               intent(in)    :: split_T             !< T_e holds the total temperature of the file

integer    :: i, ivar, n_neg, n_corr
real*8     :: psi, R, Z, T_e_min
real*8     :: zf, df_dpsi, df_dz, df_dpsi2, df_dz2, df_dpsi_dz, df_dpsi3, df_dpsi_dz2, df_dpsi2_dz
real*8     :: Omega, dOmega_dpsi, dOmega_dz, dOmega_dpsi2, dOmega_dz2, dOmega_dpsi_dz
logical    :: use_profile(n_var), use_velocity, corrected

write(*,*) '***************************************'
write(*,*) '*  init missing variables  (600)      *'
write(*,*) '***************************************'

! --- Variables that are initialised from their input profile
use_velocity = (abs(V_0) .ge. 1.d-19) .or. (num_rot)
use_profile(:) = .false.
do ivar = 1, n_var
  if ( var_imported(ivar) ) cycle
  use_profile(ivar) = (ivar == var_rho) .or. (ivar == var_T)    .or. (ivar == var_Ti) .or. &
                      (ivar == var_Te)  .or. (ivar == var_rhon) .or. ((ivar == var_Vpar) .and. use_velocity)
  if ( use_profile(ivar) .and. (n_order <= 3) ) then
    write(*,*) '  ', trim(variable_names(ivar)), ' : n=0 from the input profile'
  else
    write(*,*) '  ', trim(variable_names(ivar)), ' : zero'
  end if
end do

if ( n_order > 3 ) then
  write(*,*) 'WARNING: initialisation from input profiles is only implemented for n_order <= 3;'
  write(*,*) '         the missing variables stay zero.'
  use_profile(:) = .false.
end if

do i = 1, node_list%n_nodes

  psi = node_list%node(i)%values(1,1,var_psi)
  R   = node_list%node(i)%x(1,1,1)
  Z   = node_list%node(i)%x(1,1,2)

  if ( var_rho > 0 ) then
    if ( use_profile(var_rho) ) then
      call density(xpoint, xcase, Z, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, zf, df_dpsi, df_dz, df_dpsi2, df_dz2, &
                   df_dpsi_dz, df_dpsi3, df_dpsi_dz2, df_dpsi2_dz)
      call set_n0(var_rho)
    end if
  end if

  if ( var_T > 0 ) then
    if ( use_profile(var_T) ) then
      call temperature(xpoint, xcase, Z, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, zf, df_dpsi, df_dz, df_dpsi2, df_dz2, &
                       df_dpsi_dz, df_dpsi3, df_dpsi_dz2, df_dpsi2_dz)
      call set_n0(var_T)
    end if
  end if

  if ( var_Ti > 0 ) then
    if ( use_profile(var_Ti) ) then
      call temperature_i(xpoint, xcase, Z, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, zf, df_dpsi, df_dz, df_dpsi2, df_dz2, &
                         df_dpsi_dz, df_dpsi3, df_dpsi_dz2, df_dpsi2_dz)
      call set_n0(var_Ti)
    end if
  end if

  if ( var_Te > 0 ) then
    if ( use_profile(var_Te) ) then
      call temperature_e(xpoint, xcase, Z, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, zf, df_dpsi, df_dz, df_dpsi2, df_dz2, &
                         df_dpsi_dz, df_dpsi3, df_dpsi_dz2, df_dpsi2_dz)
      call set_n0(var_Te)
    end if
  end if

  if ( var_rhon > 0 ) then
    if ( use_profile(var_rhon) ) then
      call neutral_density(xpoint, xcase, Z, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, zf, df_dpsi, df_dz, df_dpsi2, df_dz2, &
                           df_dpsi_dz, df_dpsi3, df_dpsi_dz2, df_dpsi2_dz)
      call set_n0(var_rhon)
    end if
  end if

  if ( var_Vpar > 0 ) then
    if ( use_profile(var_Vpar) ) then
      if ( normalized_velocity_profile ) then
        call velocity(xpoint, xcase, Z, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, zf, df_dpsi, df_dz, df_dpsi2, df_dz2, &
                      df_dpsi_dz, df_dpsi3, df_dpsi_dz2, df_dpsi2_dz)
        call set_n0(var_Vpar)
      else
        ! --- Vpar,0 = 2piR^2/F0 * Omega_tor,0 (as in initial_conditions)
        call velocity(xpoint, xcase, Z, ES%Z_xpoint, psi, ES%psi_axis, ES%psi_bnd, Omega, dOmega_dpsi, dOmega_dz, &
                      dOmega_dpsi2, dOmega_dz2, dOmega_dpsi_dz, df_dpsi3, df_dpsi_dz2, df_dpsi2_dz)
        call set_n0_vpar_omega()
      end if
    end if
  end if

end do

! --- Split the total temperature of the file into T_i (input profile) and T_e = T - T_i
if ( split_T ) then
  write(*,*) '  T_i, T_e : n=0 T_i from the input profile, T_e = T - T_i (with corr_neg_temp); n/=0 T split equally'
  n_neg   = 0
  n_corr  = 0
  T_e_min = huge(1.d0)
  do i = 1, node_list%n_nodes
    if ( n_order <= 3 ) then
      node_list%node(i)%values(1,:,var_Te) = node_list%node(i)%values(1,:,var_Te) - node_list%node(i)%values(1,:,var_Ti)
    else
      node_list%node(i)%values(1,:,var_Ti) = 0.5d0 * node_list%node(i)%values(1,:,var_Te)
      node_list%node(i)%values(1,:,var_Te) = 0.5d0 * node_list%node(i)%values(1,:,var_Te)
      node_list%node(i)%deltas(1,:,var_Ti) = 0.5d0 * node_list%node(i)%deltas(1,:,var_Te)
      node_list%node(i)%deltas(1,:,var_Te) = 0.5d0 * node_list%node(i)%deltas(1,:,var_Te)
    end if
    if ( n_tor > 1 ) then
      node_list%node(i)%values(2:n_tor,:,var_Ti) = 0.5d0 * node_list%node(i)%values(2:n_tor,:,var_Te)
      node_list%node(i)%values(2:n_tor,:,var_Te) = 0.5d0 * node_list%node(i)%values(2:n_tor,:,var_Te)
      node_list%node(i)%deltas(2:n_tor,:,var_Ti) = 0.5d0 * node_list%node(i)%deltas(2:n_tor,:,var_Te)
      node_list%node(i)%deltas(2:n_tor,:,var_Te) = 0.5d0 * node_list%node(i)%deltas(2:n_tor,:,var_Te)
    end if
    T_e_min = min(T_e_min, node_list%node(i)%values(1,1,var_Te))
    if ( node_list%node(i)%values(1,1,var_Te) <= 0.d0 ) n_neg = n_neg + 1
    if ( n_order <= 3 ) then
      call corr_neg_temp_node(node_list%node(i)%values(1,:,var_Te), 2.d0, corrected)
      if ( corrected ) n_corr = n_corr + 1
    end if
  end do
  write(*,*) '  T_e = T - T_i: minimum before corr_neg_temp = ', T_e_min, ', corrected at ', n_corr, ' nodes'
  if ( n_neg > 0 ) then
    write(*,*) 'WARNING: T_e = T - T_i was not positive at ', n_neg, ' nodes (raised to the corr_neg_temp floor).'
    write(*,*) '         Check that the input T_i profile is below the total temperature T.'
  end if
end if

contains

  !> Set the n=0 values of ivar at node i from the profile value zf and its derivatives
  subroutine set_n0(ivar)
    integer, intent(in) :: ivar
    node_list%node(i)%values(1,1,ivar) = zf
    node_list%node(i)%values(1,2,ivar) = df_dpsi    * node_list%node(i)%values(1,2,var_psi) + df_dz * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,ivar) = df_dpsi    * node_list%node(i)%values(1,3,var_psi) + df_dz * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,ivar) = df_dpsi    * node_list%node(i)%values(1,4,var_psi) + df_dz * node_list%node(i)%x(1,4,2) &
                                       + df_dpsi2   * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%values(1,3,var_psi) &
                                       + df_dz2     * node_list%node(i)%x(1,2,2)            * node_list%node(i)%x(1,3,2)         &
                                       + df_dpsi_dz * node_list%node(i)%values(1,3,var_psi) * node_list%node(i)%x(1,2,2)         &
                                       + df_dpsi_dz * node_list%node(i)%values(1,2,var_psi) * node_list%node(i)%x(1,3,2)
    node_list%node(i)%deltas(1,:,ivar) = 0.d0
  end subroutine set_n0

  !> Set the n=0 values of Vpar at node i from the angular rotation profile Omega
  subroutine set_n0_vpar_omega()
    node_list%node(i)%values(1,1,var_Vpar) = R**2 * Omega
    node_list%node(i)%values(1,2,var_Vpar) = 2.d0 * R * Omega       * node_list%node(i)%x(1,2,1)        &
                                             + R**2   * dOmega_dpsi * node_list%node(i)%values(1,2,1)   &
                                             + R**2   * dOmega_dz   * node_list%node(i)%x(1,2,2)
    node_list%node(i)%values(1,3,var_Vpar) = 2.d0 * R * Omega       * node_list%node(i)%x(1,3,1)        &
                                             + R**2   * dOmega_dpsi * node_list%node(i)%values(1,3,1)   &
                                             + R**2   * dOmega_dz   * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,var_Vpar) =   2.d0 *     node_list%node(i)%x(1,2,1)**2  * Omega &
                                             + 2.d0 * R * node_list%node(i)%x(1,4,1)     * Omega &
                                             + 2.d0 * R * node_list%node(i)%x(1,2,1)     * dOmega_dpsi  * node_list%node(i)%values(1,3,1) &
                                             + 2.d0 * R * node_list%node(i)%x(1,2,1)     * dOmega_dz    * node_list%node(i)%x(1,3,2)
    node_list%node(i)%values(1,4,var_Vpar) = node_list%node(i)%values(1,4,var_Vpar) &
                                             + 2.d0 * R * dOmega_dpsi    * node_list%node(i)%x(1,3,1)      * node_list%node(i)%values(1,2,1) &
                                             + R**2     * dOmega_dpsi2   * node_list%node(i)%values(1,3,1) * node_list%node(i)%values(1,2,1) &
                                             + R**2     * dOmega_dpsi_dz * node_list%node(i)%x(1,3,2)      * node_list%node(i)%values(1,2,1) &
                                             + R**2     * dOmega_dpsi    * node_list%node(i)%values(1,4,1)
    node_list%node(i)%values(1,4,var_Vpar) = node_list%node(i)%values(1,4,var_Vpar) &
                                             + 2.d0 * R * dOmega_dz      * node_list%node(i)%x(1,3,1)      * node_list%node(i)%x(1,2,2) &
                                             + R**2     * dOmega_dpsi_dz * node_list%node(i)%values(1,3,1) * node_list%node(i)%x(1,3,2) &
                                             + R**2     * dOmega_dz2     * node_list%node(i)%x(1,3,2)      * node_list%node(i)%x(1,3,2) &
                                             + R**2     * dOmega_dz      * node_list%node(i)%x(1,4,2)
    node_list%node(i)%values(1,:,var_Vpar) = 2.d0 * PI / F0 * node_list%node(i)%values(1,:,var_Vpar)
    node_list%node(i)%deltas(1,:,var_Vpar) = 0.d0
  end subroutine set_n0_vpar_omega

end subroutine init_missing_variables
