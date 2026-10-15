module corr_neg

implicit none

  interface corr_neg_dens
    module procedure corr_neg_dens1, corr_neg_dens2, corr_neg_dens3
  end interface

  interface dcorr_neg_dens_drho
    module procedure dcorr_neg_dens_drho1, dcorr_neg_dens_drho2, dcorr_neg_dens_drho3
  end interface
      
  interface corr_neg_temp
    module procedure corr_neg_temp1, corr_neg_temp2, corr_neg_temp3
  end interface

  interface dcorr_neg_temp_dT
    module procedure dcorr_neg_temp_dT1, dcorr_neg_temp_dT2, dcorr_neg_temp_dT3
  end interface

  interface d2corr_neg_temp_dT2
    module procedure d2corr_neg_temp_dT21, d2corr_neg_temp_dT22, d2corr_neg_temp_dT23
  end interface


contains

! Source code moved to an include file to be able to inline function calls into the loops
! in element_matrix. This module file is kept for those models which do not use inlining.
!
! For optimized element_matrix construction subroutines, inlining the corr_neg functions improves
! the performance. But the -ipo switch, which would enable inlining across compilation units, makes
! the compilation very slow and certain other subroutines become slower with it. As an alternative
! to -ipo, we separete the function definitions into the corr_neg_include.f90 file, which can be
! included into the source files that define element matrix. This way the compiler will be able to
! inline the corr_neg functions.
!
! The include directive is used here to avoid code duplication.
#include "corr_neg_include.f90"



!> Apply corr_neg_temp to the n=0 temperature at one node, given as (value, d/ds, d/dt, d2/dsdt)
!! (n_order = 3). Use scale = 1 for the total temperature T and scale = 2 for T_i and T_e, as in
!! the element matrix of model600: corr_neg_temp(2*T_i)/2. The derivatives follow the chain rule.
subroutine corr_neg_temp_node(v, scale, corrected)

  real*8,  intent(inout) :: v(4)       !< n=0 temperature at the node
  real*8,  intent(in)    :: scale      !< 1 for T, 2 for T_i and T_e
  logical, intent(out)   :: corrected  !< .true. if the value was changed

  real*8 :: f0, f1, f2, v4

  f0 = corr_neg_temp(scale * v(1)) / scale
  f1 = dcorr_neg_temp_dT(scale * v(1))
  f2 = scale * d2corr_neg_temp_dT2(scale * v(1))

  corrected = (f0 /= v(1))
  v4   = f1 * v(4) + f2 * v(2) * v(3)
  v(1) = f0
  v(2) = f1 * v(2)
  v(3) = f1 * v(3)
  v(4) = v4

end subroutine corr_neg_temp_node

end module corr_neg
