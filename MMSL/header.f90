  !-----------------------------------------------------------------------------
  ! header in MMSL
  !
  ! header
  !
  !-----------------------------------------------------------------------------
  
  !------------------------------------------------------------
  ! Solvers
  !------------------------------------------------------------
  ! MINRES
  include '../MMSL/minr1d'
  include '../MMSL/minr1d2'
  ! Right preconditioned MINRES
  include '../MMSL/rpmr1d'
  include '../MMSL/icmr1d'
  include '../MMSL/rpmr1dtest'
  include '../MMSL/icmr1dtest'
  ! E-SSOR preconditioned MINRES
  include '../MMSL/esmr1d'
  include '../MMSL/esmr1d2'
  include '../MMSL/icmr1d2'
  include '../MMSL/rpmr1d2'
  ! ICCG
  include '../MMSL/m_iccg1d'
  include '../MMSL/m_cgic1d'
  ! SOR
  include '../MMSL/icsor1d'
  include '../MMSL/sor1d'
  ! MR
  include '../MMSL/mr1d'
  ! Right preconditioned MR
  include '../MMSL/pmr1d'
  include '../MMSL/wpmr1d'
  ! CR
  include '../MMSL/cr1d'
  ! Right preconditioned CR
  include '../MMSL/pcr1d'
  include '../MMSL/wpcr1d'
  ! MrR
  include '../MMSL/mrr1d'
  ! Right preconditioned MrR
  include '../MMSL/pmrr1d'
  include '../MMSL/wpmrr1d'
  ! GD
  include '../MMSL/gd1d'
  ! Right preconditioned GD
  include '../MMSL/pgd1d'
  include '../MMSL/wpgd1d'
  ! Variable preconditioned solvers
  include '../MMSL/vpcgsor1d'
  include '../MMSL/vpcgminr1d'
  include '../MMSL/vpcgiccg1d'
  !------------------------------------------------------------
  ! Others
  !------------------------------------------------------------
  ! inverse matrix
  ! include '../MMSL/inverse'
  ! depress
  include '../MMSL/m_conv_2_iccg'
  include '../MMSL/m_upper_from_iccg'
  ! incomplete cholesky decomposition
  include '../MMSL/m_sqrt_icdc'
  include '../MMSL/m_sqrt_icdc1d'
  ! A-1x
  include '../MMSL/m_icsl1d'
  ! Ax
  include '../MMSL/m_prod1d'
  
  include '../MMSL/m_icdec'
  include '../MMSL/m_icdc1d'
  
