subroutine wpcr1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt,&
     uu,ul,ud)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n,n_nzr
  real(8),intent(in)::au(n_nzr),ad(n),al(n_nzr),ksu(n),ksl(n),uu(n_nzr),ul(n_nzr),ud(n)
  integer,intent(in)::iu(n_nzr),il(n_nzr),nzu(n),nzl(n)
  real(8),intent(in)::eps,b(n)
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)
  real(8)::r0
  real(8),dimension(n)::r,p,&
       kr,akr,ap,kap,&
       kr1,akr1,ap1!,kap1
  
  call prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,x,ap,int(ksu),int(ksl))
  r = b - ap
  r0 = sqrt(sum(r**2))
  
  call icsl1d(n,n_nzr,nzu,nzl,iu,il,uu,ul,ud,r,p,int(ksu),int(ksl))

  call icsl1d(n,n_nzr,nzu,nzl,iu,il,uu,ul,ud,r,kr,int(ksu),int(ksl))
  call prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,kr,akr,int(ksu),int(ksl))

  call prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,p,ap,int(ksu),int(ksl))
  call icsl1d(n,n_nzr,nzu,nzl,iu,il,uu,ul,ud,ap,kap,int(ksu),int(ksl))
  
  do j=1,itmax
     
     alpha = dot_product(akr,kr)/dot_product(kap,ap)
     
     x = x + alpha * p
     
     r = r - alpha * ap

     kr1 = kr - alpha * kap
     call prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,kr1,akr1,int(ksu),int(ksl))
     
     beta = dot_product(akr1,kr1)/dot_product(akr,kr)

     p = kr1 + beta * p

     ap1 = akr1 + beta * ap

     call icsl1d(n,n_nzr,nzu,nzl,iu,il,uu,ul,ud,ap1,kap,int(ksu),int(ksl))
     
     error = sqrt(sum(r**2)) / r0
     
     if(iopt/=0)write(*,*)j,',',error
     if(error<eps) return
     
     kr=kr1
     akr=akr1
     ap=ap1
     ! kap=kap1
  end do
  
  write(*,*)'** does not converge. **'
  return
end subroutine wpcr1d
