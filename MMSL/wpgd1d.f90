subroutine wpgd1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt,&
     uu,ul,ud)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n_nzr,n
  real(8),intent(in)::au(n_nzr),ad(n),al(n_nzr),ksu(n),ksl(n),uu(n_nzr),ul(n_nzr),ud(n)
  integer,intent(in)::iu(n_nzr),il(n_nzr),nzu(n),nzl(n)
  real(8),intent(in)::eps,b(n)
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)
  
  real(8)::r(n),kr(n),akr(n),alpha,r0


  call prod1d(n,nzu,nzl,iu,il,au,al,ad,x,r,int(ksu),int(ksl))
  r = b - r
  r0 = sum(r**2)

  do j=1,itmax
     
     call icsl1d(n,nzu,nzl,iu,il,uu,ul,ud,r,kr,int(ksu),int(ksl))
     call prod1d(n,nzu,nzl,iu,il,au,al,ad,kr,akr,int(ksu),int(ksl))

     alpha = dot_product(kr,r) / dot_product(kr,akr)
     
     x = x + alpha * kr

     r = r - alpha * akr

     error = sqrt(sum(r**2) / r0)
     
     if(iopt/=0)write(*,*)j,',',error
     if(error<eps .or. isnan(error)) return
     
  end do
  
  write(*,*) '** does not converge. **'
  return
end subroutine wpgd1d
