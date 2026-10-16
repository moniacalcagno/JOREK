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

if ( (i_var < 1) .or. (i_var > n_eq_var) ) then
  name  = 'NONE'
  descr = 'Illegal variable index!'
  return
else if ( i_var == var_rho_frozen ) then
  name  = 'rho'
  descr = 'particle density (frozen, not evolved)'
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



!> Canonical name of a variable ('psi', 'u', 'zj', 'w', 'rho', 'T', 'Ti', 'Te', 'Vpar', ...).
!! Accepts the names of the var_* constants with or without the prefix 'var_' and the names
!! returned by variable_names (e.g. 'var_w', 'w' and 'omega'; 'var_Vpar', 'vpar' and 'v_par'),
!! case insensitive. Returns '' for an unknown name.
character(len=8) pure function variable_canonical_name(var_name)

implicit none

character(len=*), intent(in) :: var_name

character(len=len(var_name)) :: lname, sname
integer :: i, ic, n

! --- Lower case
lname = adjustl(var_name)
do i = 1, len_trim(lname)
  ic = iachar(lname(i:i))
  if ( (ic >= iachar('A')) .and. (ic <= iachar('Z')) ) lname(i:i) = achar(ic + 32)
end do

! --- Remove the prefix 'var_' and all underscores
if ( lname(1:min(4,len(lname))) == 'var_' ) lname = lname(5:)
sname = ''
n = 0
do i = 1, len_trim(lname)
  if ( lname(i:i) /= '_' ) then
    n = n + 1
    sname(n:n) = lname(i:i)
  end if
end do

select case ( trim(sname) )
  case ('psi');            variable_canonical_name = 'psi'
  case ('u');              variable_canonical_name = 'u'
  case ('zj');             variable_canonical_name = 'zj'
  case ('w', 'omega');     variable_canonical_name = 'w'
  case ('rho');            variable_canonical_name = 'rho'
  case ('t');              variable_canonical_name = 'T'
  case ('ti');             variable_canonical_name = 'Ti'
  case ('te');             variable_canonical_name = 'Te'
  case ('vpar');           variable_canonical_name = 'Vpar'
  case ('rhon');           variable_canonical_name = 'rhon'
  case ('rhoimp');         variable_canonical_name = 'rhoimp'
  case ('nre');            variable_canonical_name = 'nre'
  case ('ar');             variable_canonical_name = 'AR'
  case ('az');             variable_canonical_name = 'AZ'
  case ('a3');             variable_canonical_name = 'A3'
  case ('ur');             variable_canonical_name = 'uR'
  case ('uz');             variable_canonical_name = 'uZ'
  case ('up', 'uphi');     variable_canonical_name = 'up'
  case default;            variable_canonical_name = ''
end select

end function variable_canonical_name



!> Index of a variable in the compiled model from its name (see variable_canonical_name).
!! Returns -1 for an unknown name and 0 for a variable that does not exist in the compiled model.
!! The density refers to the frozen density (var_rho_frozen) in a model that stores but does not evolve it.
integer pure function variable_index_from_name(var_name)

implicit none

character(len=*), intent(in) :: var_name

select case ( trim(variable_canonical_name(var_name)) )
  case ('psi');    variable_index_from_name = var_psi
  case ('u');      variable_index_from_name = var_u
  case ('zj');     variable_index_from_name = var_zj
  case ('w');      variable_index_from_name = var_w
  case ('rho');    variable_index_from_name = merge(var_rho, var_rho_frozen, var_rho > 0)
  case ('T');      variable_index_from_name = var_T
  case ('Ti');     variable_index_from_name = var_Ti
  case ('Te');     variable_index_from_name = var_Te
  case ('Vpar');   variable_index_from_name = var_Vpar
  case ('rhon');   variable_index_from_name = var_rhon
  case ('rhoimp'); variable_index_from_name = var_rhoimp
  case ('nre');    variable_index_from_name = var_nre
  case ('AR');     variable_index_from_name = var_AR
  case ('AZ');     variable_index_from_name = var_AZ
  case ('A3');     variable_index_from_name = var_A3
  case ('uR');     variable_index_from_name = var_uR
  case ('uZ');     variable_index_from_name = var_uZ
  case ('up');     variable_index_from_name = var_up
  case default;    variable_index_from_name = -1
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
