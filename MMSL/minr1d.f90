!-------------------------------------------------------------------------------
! minres
!
! MINRES Method (ICCG format)
!
! by Moeto Suzuki
!
! arguments
! [input]
! kmax                           : number of nonzero elements
!                                  in strict upper triangular matrix
! n                              : size of symmetric matrix
! nzu,nzl,iu,il,au,al,ad,ksu,ksl : sparse symmetric matrix data
! b (n)                          : right-hand side vector
! x (n)                          : initial guess for solution
! eps                            : convergence determinant
! itmax                          : maximum number of iteration
! iopt                           : output residual history or not
!
! [output]
! x  (n)                         : solution of simultaneous equations
! 
!-------------------------------------------------------------------------------
subroutine minr1d(kmax,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::kmax,n
  real(8),intent(in)::au(kmax),ad(n),al(kmax),ksu(n),ksl(n)
  integer,intent(in)::iu(kmax),il(kmax),nzu(n),nzl(n)
  real(8),intent(in)::eps,b(n)
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)

  real(8)::v,r

  real(8),dimension(n)::v0vec,v1vec,v2vec,w0vec,w1vec,w2vec
  real(8)::gamma1,gamma2,eta,s0,s1,c0,c1,error, &
       alpha0,alpha1,alpha2,alpha3
  real(8),dimension(n)::ap,ax
  !
  ! Init v0,w0,w1
  !
  v0vec = 0.d0
  w0vec = 0.d0
  w1vec = 0.d0
  !
  ! Set v1vec
  !
  call m_prod1d(n,kmax,nzu,nzl,iu,il,au,al,ad,x,ap,int(ksu),int(ksl))
  v1vec = b - ap
  v = sqrt(sum(v1vec**2))
  r = v
  !
  ! Set initial values
  !
  gamma1 = sqrt(sum(v1vec**2))
  eta = gamma1
  s0 = 0.d0
  s1 = 0.d0
  c0 = 1.d0
  c1 = 1.d0
  !
  ! Roop until convergence
  !
  do j=1,itmax
     ! define v1(j)
     v1vec = v1vec / gamma1
     call m_prod1d(n,kmax,nzu,nzl,iu,il,au,al,ad,v1vec,ap,int(ksu),int(ksl))
     delta = dot_product(ap,v1vec)
     v2vec = ap - delta * v1vec - gamma1 * v0vec
     gamma2 = sqrt(sum(v2vec**2))
     !
     ! calculate values
     !
     alpha0 = c1*delta - c0 * s1 * gamma1
     alpha1 = sqrt(alpha0**2 + gamma2**2)
     alpha2 = s1*delta + c0*c1*gamma1
     alpha3 = s0 * gamma1
     c0 = c1
     c1 = alpha0 / alpha1
     s0 = s1
     s1 = gamma2 / alpha1
     !
     ! calculate w2
     !
     w2vec = (v1vec - alpha3 * w0vec - alpha2 * w1vec) / alpha1
     !
     ! update x
     !
     x = x + c1 * eta * w2vec
     !
     ! Calculate relative error
     !
     eta = -s1 * eta
     r = abs(s1) * r
     error = r/v
     
     !
     ! output residual history
     !
     if(iopt/=0)write(*,*) j,',',error
     if (isnan(error)) stop '"x" is a NaN'
     if(error<eps) exit
     if(j==itmax) write(*,*) '**** does not converge. ****'
     !
     ! Update values
     !
     gamma1 = gamma2
     w0vec = w1vec
     w1vec = w2vec
     v0vec = v1vec
     v1vec = v2vec
  end do
  return
end subroutine minr1d
