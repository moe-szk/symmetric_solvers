subroutine  wpmrr1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt&
   ,uu,ul,ud)
  integer,intent(in)::n_nzr,n
  integer,intent(in)::nzu(n),nzl(n),iu(n_nzr),il(n_nzr)
  real(8),intent(in)::au(n_nzr),al(n_nzr),ad(n),ksu(n_nzr),ksl(n_nzr)
  real(8),intent(in)::b(n),eps,uu(n_nzr),ul(n_nzr),ud(n)
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)
  real(8),dimension(n)::r,y,z,ar,rd,s,kr
  real(8)::mu,nu,omg,zeta,gamma1,gamma2,err,r0
  
  call prod1d(n,nzu,nzl,iu,il,au,al,ad,x,ar,int(ksu),int(ksl))
  
  r = b - ar
  r0 = sqrt(sum(r**2))
  
  call icsl1d(n,nzu,nzl,iu,il,uu,ul,ud,r,kr,int(ksu),int(ksl))
  y = - kr
  
  z = 0.d0

  do j=1,itmax

     !call icsl1d(n,nzu,nzl,iu,il,uu,ul,ud,r1,kr,ksu,ksl)
     call prod1d(n,nzu,nzl,iu,il,au,al,ad,r,ar,int(ksu),int(ksl))
     
     mu = dot_product(y,y)
     nu = dot_product(y,ar)
     omg = dot_product(y,r)
     
     gamma1 = omg / mu
     gamma2 = nu / mu
     if(j==1)gamma1=0.d0;gamma2=0.d0
     
     rd = r - gamma1 * y
     
     s = ar - ganma2 * y

     zeta = dot_product(rd,s)/dot_product(s,s)
     eta = gamma1 - zeta * gamma2

     call icsl1d(n,nzu,nzl,iu,il,uu,ul,ud,r,kr,int(ksu),int(ksl))
     !y =  eta * y + zeta * ar
     y = eta * y + zeta * kr
     z = eta * z - zeta * r

     x = x - z
     r = r - y
     err = sqrt(sum(r**2))/r0

     if(iopt/=0)write(*,*)j,err
     if(err<eps)return
     if(isnan(err))return

  end do
  write(*,*)'** does not converge. **'
  return
end subroutine wpmrr1d
