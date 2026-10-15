module mod_parameters

use mod_settings
use mod_model_settings

implicit none



contains



! --- Variable names
character(len=64) pure elemental function variable_names(i_var)

implicit none

integer,           intent(in) :: i_var !< Variable index
!logical, optional, intent(in) :: long  !< Long output with description?

character(len=20) :: name
character(len=64) :: descr
logical :: long2

long2 = .false.

if ( (i_var < 1) .or. (i_var > n_var) ) then
  name  = 'NONE'
  descr = 'Illegal variable index!'
  return
else if ( i_var == var_psi ) then
  name  = 'Psi'
  descr = 'poloidal magnetic flux'
else if ( i_var == var_u ) then
  name  = 'u'
  descr = 'velocity stream function'
else if ( i_var == var_zj ) then
  name  = 'zj'
  descr = 'toroidal plasma current'
else if ( i_var == var_w ) then
  name  = 'omega'
  descr = 'toroidal vorticity'
else if ( i_var == var_rho ) then
  name  = 'rho'
  descr = 'particle density'
else if ( i_var == var_T ) then
  name  = 'T'
  descr = 'temperature; sum of electrons and ions'
else if ( i_var == var_Ti ) then
  name  = 'T_i'
  descr = 'ion temperature'
else if ( i_var == var_Te ) then
  name  = 'T_e'
  descr = 'electron temperature'
else if ( i_var == var_Vpar ) then
  name  = 'v_par'
  descr = 'parallel velocity'
else if ( i_var == var_rhon ) then
  name  = 'rho_n'
  descr = 'neutral particle density'
else if ( i_var == var_rhoimp ) then
  name  = 'rho_imp'
  descr = 'impurity density'
else if ( i_var == var_nre ) then
  name  = 'n_RE'
  descr = 'runaway electron particle density'
else if ( i_var == var_AR ) then
  name  = "AR"
  descr = "R component of magnetic vector potential"
else if ( i_var == var_AZ ) then
  name  = "AZ"
  descr = "Z component of magnetic vector potential"
else if ( i_var == var_uR ) then
  name  = "UR"
  descr = "Velocity in R direction"
else if ( i_var == var_uZ ) then
  name  = "UZ"
  descr = "Velocity in Z direction"
else if ( i_var == var_up ) then
  name  = "Uphi"
  descr = "Velocity in phi direction"
else
  write(name, "(I4)") i_var
  descr = "Please check implementation"
end if

if ( .not. long2 ) then
  variable_names = trim(name)
else
  variable_names = trim(name)//' ('//trim(descr)//')'
end if

end function variable_names



!> Index of a variable in the compiled model from its name ('var_psi', 'var_u', ...).
!! Returns 0 if the variable is unknown or does not exist in the compiled model.
integer pure function variable_index_from_name(var_name)

implicit none

character(len=*), intent(in) :: var_name

character(len=len(var_name)) :: lname
integer :: i, ic

lname = adjustl(var_name)
do i = 1, len_trim(lname)
  ic = iachar(lname(i:i))
  if ( (ic >= iachar('A')) .and. (ic <= iachar('Z')) ) lname(i:i) = achar(ic + 32)
end do

select case ( trim(lname) )
  case ('var_psi');    variable_index_from_name = var_psi
  case ('var_u');      variable_index_from_name = var_u
  case ('var_zj');     variable_index_from_name = var_zj
  case ('var_w');      variable_index_from_name = var_w
  case ('var_rho');    variable_index_from_name = var_rho
  case ('var_t');      variable_index_from_name = var_T
  case ('var_ti');     variable_index_from_name = var_Ti
  case ('var_te');     variable_index_from_name = var_Te
  case ('var_vpar');   variable_index_from_name = var_Vpar
  case ('var_rhon');   variable_index_from_name = var_rhon
  case ('var_rhoimp'); variable_index_from_name = var_rhoimp
  case ('var_nre');    variable_index_from_name = var_nre
  case ('var_ar');     variable_index_from_name = var_AR
  case ('var_az');     variable_index_from_name = var_AZ
  case ('var_a3');     variable_index_from_name = var_A3
  case ('var_ur');     variable_index_from_name = var_uR
  case ('var_uz');     variable_index_from_name = var_uZ
  case ('var_up');     variable_index_from_name = var_up
  case default;        variable_index_from_name = 0
end select

end function variable_index_from_name



#ifdef MODEL_FAMILY

!> Are all selected extensions available and compatible?
logical pure function ext_correct()

implicit none

! --- Local variables
integer :: iext1, iext2

! --- Preset to .true.
ext_correct = .true.

! --- Check binary compatibility of all extensions
do iext1 = 1, n_mod_ext
  if ( .not. with_ext(iext1) ) cycle
  ext_correct = ext_correct .and. ext_available(iext1)
  do iext2 = 1, n_mod_ext
    if ( .not. with_ext(iext2) ) cycle
      ext_correct = ext_correct .and. ext_compatible(iext1, iext2)
  end do
end do

end function ext_correct

#endif



end module mod_parameters
