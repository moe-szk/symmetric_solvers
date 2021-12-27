subroutine  wpmrr1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt&
   , uu,ul, ud)
  integer,intent(in)::n_nzr,n
  integer,intent(in)::nzu(n),nzl(n),iu(n_nzr),il(n_nzr)
  real(8),intent(in)::au(n_nzr),al(n_nzr),ad(n),ksu(n_nzr),ksl(n_nzr)
  real(8),intent(in)::b(n),eps,uu(n_nzr),ul(n_nzr),ud(n)
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)
  real(8),dimension(n)::r,p,z,ar,rd,s,kr
  real(8)::mu,nu,omg,zeta,alpha,beta,err,r0
  
  call prod1d(n,nzu,nzl,iu,il,au,al,ad,x,ar,int(ksu),int(ksl))
  
  r = b - ar
  r0 = sqrt(sum(r**2))
  call icsl1d(n,nzu,nzl,iu,il,uu,ul,ud,r,kr,int(ksu),int(ksl))
  p =  -kr
  
  z = 0.d0

  do j=1,itmax

     !call prod1d(n,nzu,nzl,iu,il,au,al,ad,r,ar,int(ksu),int(ksl))
     call icsl1d(n,nzu,nzl,iu,il,uu,ul,ud,r,kr,int(ksu),int(ksl))
     call prod1d(n,nzu,nzl,iu,il,au,al,ad,kr,ar,int(ksu),int(ksl))
     
     mu = dot_product(p,p)
     nu = dot_product(p,ar)
     !nu = dot_product(y,kr)
     omg = dot_product(p,r)
     
     alpha = omg / mu
     beta =  nu / mu
     if(j==1)gamma1=0.d0;gamma2=0.d0
     
     rd = r - alpha * p
     
     s = ar - beta * p

     zeta = dot_product(rd,s)/dot_product(s,s)
     eta = alpha - zeta * beta

     !call icsl1d(n,nzu,nzl,iu,il,uu,ul,ud,ar,kr,int(ksu),int(ksl))
     
     p =  eta * p + zeta * kr
     z = eta * z - zeta * r

     x = x - z
     r = r - p
     err = sqrt(sum(r**2))/r0

     if(iopt/=0)write(*,*)j,err
     if(err<eps)return
     if(isnan(err))return

  end do
  write(*,*)'** does not converge. **'
  return
end subroutine wpmrr1d
