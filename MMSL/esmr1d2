!-------------------------------------------------------------------------------
! esmr1d
!
! Eisenstat SSOR for Right preconditioned MINRES
!
! by Moeto Suzuki
!
!-------------------------------------------------------------------------------
subroutine esmr1d2(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt, &
     omg)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n_nzr,n
  real(8),intent(in)::au(n_nzr),ad(n),al(n_nzr),ksu(n),ksl(n)
  integer,intent(in)::iu(n_nzr),il(n_nzr),nzu(n),nzl(n)
  real(8),intent(in)::eps,b(n),omg
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)
  real(8)::r(n),v(n),y(n),vbar(n)
  real(8),dimension(n)::v0vec,v1vec,v2vec,w0vec,w1vec,w2vec,u1vec,u2vec
  real(8)::gamma1,gamma2,eta,s0,s1,c0,c1,error, &
       alpha0,alpha1,alpha2,alpha3
  real(8),dimension(n)::ap,ax
  real(8)::zu(n_nzr),zl(n_nzr)
  real(8)::uu(n_nzr),ul(n_nzr),ud(n)
  zu = 0.d0
  zl = 0.d0
  call m_icdc1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,ksu,ksl,uu,ul,ud)
  !
  ! Initialization
  !
  v0vec = 0.d0
  w0vec = 0.d0
  w1vec = 0.d0
  
  call m_prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,x,ap,int(ksu),int(ksl))
  v = b - ap
  
  call m_icsl1d(n,n_nzr,nzu,nzl,iu,il,zu,ul,abs(ud)/omg,v,ap,int(ksu),int(ksl))
  call m_prod1d(n,n_nzr,nzu,nzl,iu,il,zu,zl,sqrt(abs(ud)),ap,ax,int(ksu),int(ksl))
  v1vec = ax

  
  omg1 = (2.d0-omg)/omg
  gamma1 = sqrt(omg1 * dot_product(v1vec,v1vec))
  eta = gamma1
  s0 = 0.d0
  s1 = 0.d0
  c0 = 1.d0
  c1 = 1.d0
  !
  ! Iterative loop
  !
  do j=1,itmax
     !
     ! Define v1,u1
     ! 
     v1vec = v1vec / gamma1
     !
     ! Eisenstat's trick
     !
     call m_prod1d(n,n_nzr,nzu,nzl,iu,il,zu,zl,sqrt(abs(ud)),v1vec,vbar,int(ksu),int(ksl))
     call m_icsl1d(n,n_nzr,nzu,nzl,iu,il,uu,zl,abs(ud)/omg,vbar,y,int(ksu),int(ksl))
     call m_prod1d(n,n_nzr,nzu,nzl,iu,il,zu,zl,((2.d0/omg*abs(ud))-ud),y,ap,int(ksu),int(ksl))
     ax = vbar - ap
     call m_icsl1d(n,n_nzr,nzu,nzl,iu,il,zu,ul,abs(ud)/omg,ax,ap,int(ksu),int(ksl))
     call m_prod1d(n,n_nzr,nzu,nzl,iu,il,zu,zl,sqrt(abs(ud)),y+ap,ax,int(ksu),int(ksl))
     u1vec = ap
     
     delta = omg1**2 * dot_product(u1vec,v1vec)

     v2vec = omg1 * u1vec - delta * v1vec - gamma1 * v0vec

     gamma2 = sqrt(omg1 * dot_product(v2vec,v2vec))

     alpha0 = c1*delta - c0*s1*gamma1
     alpha1 = sqrt(alpha0**2 + gamma2**2)
     alpha2 = s1*delta + c0*c1*gamma1
     alpha3 = s0*gamma1
     c0 = c1
     c1 = alpha0/alpha1
     s0 = s1
     s1 = gamma2 / alpha1

     w2vec = (omg1*y - alpha3*w0vec - alpha2*w1vec)/alpha1

     x = x + c1*eta * w2vec
     eta = - s1 * eta

     call m_prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,x,ap,int(ksu),int(ksl))
     r = b - ap
     !call icsl1d(n,nzu,nzl,iu,il,au,al,ad,r,ap,int(ksu),int(ksl))
     !call m_prod1d(n,kmax,nzu,nzl,iu,il,au,al,ad,ap,ax,int(ksu),int(ksl))
     error = sqrt(sum(r**2))/sqrt(sum(v**2))
     !error = sqrt(sum(ax**2))/ sqrt(sum(v**2))
     
     if(iopt/=0)write(*,*)j,',',error
     if(error<eps) return
     if(isnan(error))stop'Nan'
     if(itmax==j)write(*,*)'** does not converge **'
     
     gamma1 = gamma2
     w0vec = w1vec
     w1vec = w2vec
     v0vec = v1vec
     v1vec = v2vec
     u1vec = u2vec
  end do
  
end subroutine esmr1d2
