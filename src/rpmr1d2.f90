!-------------------------------------------------------------------------------
! rpmr1d2
!
! Right preconditioned MINRES Method with restart (ICCG format)
!
! arguments
! [input]
! kmax                           : number of nonzero elements (adim)
!                                  in strict upper triangular matrix
! n                              :  size of symmetric matrix
! nzu,nzl,iu,il,au,al,ad,ksu,ksl : sparse symetric matrix data(ICCG format)
! b (n)                          : right-hand side vector
! x (n)                          : initial guess for solution
! eps                            : convergence determinant
! itmax                          : maximum number of iteration
! iopt                           : output residual history or not
! mmax                           : number of nonzero elements (m1dim)
!                                  in strict upper triangular matrix
! mu,ml,md                       : preconditioner matrix data(ICCG format)
! k                              : restart parameter
!
! [output]
! x  (n)    : solution of simultaneous equations
! 
!-------------------------------------------------------------------------------
subroutine rpmr1d2(kmax,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt,&
     mu,ml,md,k)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::kmax,n,k
  real(8),intent(in)::au(kmax),ad(n),al(kmax),ksu(n),ksl(n)
  integer,intent(in)::iu(kmax),il(kmax),nzu(n),nzl(n)
  real(8),intent(in)::mu(kmax),ml(kmax),md(n)
  real(8),intent(in)::eps,b(n)
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)
  real(8)::r,v,rk,km
  real(8),dimension(n)::v0vec,v1vec,v2vec,w0vec,w1vec,w2vec,u1vec,u2vec
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
  ! Set v1vec,u1vec
  !
  call m_prod1d(n,kmax,nzu,nzl,iu,il,au,al,ad,x,ap,int(ksu),int(ksl))
  v1vec = b - ap
  !call m_icsl1d(n,kmax,nzu,nzl,iu,il,mu,ml,md,v1vec,ap,int(ksu),int(ksl))
  !call m_prod1d(n,kmax,nzu,nzl,iu,il,au,al,ad,ap,ax,int(ksu),int(ksl))
  !v     = ax
  v  = sqrt(sum(v1vec**2))
  r  = v
  rk = v
  call m_icsl1d(n,kmax,nzu,nzl,iu,il,mu,ml,md,v1vec,ap,int(ksu),int(ksl))
  u1vec = ap
  !
  ! Set initial values
  !
  gamma1 = sqrt(abs(dot_product(v1vec,u1vec)))
  eta = gamma1
  s0 = 0.d0
  s1 = 0.d0
  c0 = 1.d0
  c1 = 1.d0
  iii=0
  !
  ! Loop until convergence
  !
  do j=1,itmax
     ! Define v1,u1
     v1vec = v1vec / gamma1
     u1vec = u1vec / gamma1
     !
     ! Calculate v2,u2
     !
     call m_prod1d(n,kmax,nzu,nzl,iu,il,au,al,ad,u1vec,ap,int(ksu),int(ksl))
     delta = dot_product(u1vec,ap)
     v2vec = ap - delta * v1vec - gamma1 * v0vec
     ap =0.d0
     call m_icsl1d(n,kmax,nzu,nzl,iu,il,mu,ml,md,v2vec,ap,int(ksu),int(ksl))
     u2vec = ap
     !
     ! Calculate values
     !
     gamma2 = sqrt(abs(dot_product(v2vec,u2vec)))
     alpha0 = c1*delta - c0 * s1 * gamma1
     alpha1 = sqrt(alpha0**2 + gamma2**2)
     alpha2 = s1*delta + c0*c1*gamma1
     alpha3 = s0 * gamma1
     c0 = c1
     c1 = alpha0 / alpha1
     s0 = s1
     s1 = gamma2 / alpha1
     !
     ! Calculate w2
     !
     w2vec = (u1vec - alpha3 * w0vec - alpha2 * w1vec) / alpha1
     !
     ! Update x
     !
     x = x + c1 * eta * w2vec
     !
     ! Calculate relative error
     !
     eta = -s1 * eta
     r = abs(s1)*r
     error = r / v
     !
     ! output residual history
     !
     if(iopt/=0) write(*,*) j,',',error
     if(isnan(error)) stop '"x" is a NaN'
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
     u1vec = u2vec
     if(mod(j,k)==0)then
        !rm = r
        !if(abs(rk-rm)/abs(v) <= eps/real(k) )then
           !iii=iii+1
           !write(*,*)'Restart',iii
           v0vec = 0.d0
           w0vec = 0.d0
           w1vec = 0.d0
           call m_prod1d(n,kmax,nzu,nzl,iu,il,au,al,ad,x,ap,int(ksu),int(ksl))
           v1vec = b - ap
           call m_icsl1d(n,kmax,nzu,nzl,iu,il,mu,ml,md,v1vec,ap,int(ksu),int(ksl))
           u1vec = ap
           gamma1 = sqrt(abs(dot_product(v1vec,u1vec)))
           eta = gamma1
           s0 = 0.d0
           s1 = 0.d0
           c0 = 1.d0
           c1 = 1.d0
        !end if
        !rk = rm
     end if
  end do
  return
end subroutine rpmr1d2
