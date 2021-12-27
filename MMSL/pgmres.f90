subroutine pgmres(n,n_nzr,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,k,itmax,eps,x,iopt,&
     uu,ul,ud)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n,n_nzr,itmax,iopt,k,nzu(n),nzl(n),&
       iu(n_nzr),il(n_nzr),ksu(n),ksl(n)
  real(8),intent(in)::b(n),eps,au(n_nzr),al(n_nzr),ad(n),uu(n_nzr),ul(n_nzr),ud(n)
  real(8),intent(inout)::x(n)
  real(8),dimension(n)::v0,v1,v2,r,r0,u0,u1,u2,kv,akv
  real(8),dimension(k)::e
  
  r0 = b - dot_product(adim,x)
  r = r0

  v1 = r0 / sqrt(sum(r0**2))

  e(1) = sqrt(sum(r0**2))
  
  do j=1,itmax
     call m_icsl1d(n,n_nzr,nzu,nzl,iu,il,uu,ul,ud,v1,kv,ksu,ksl)
     call m_prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,kv,u2,ksu,ksl)

     
     
  end do
  write(*,*)'** does note converge. **'
  return
end subroutine pgmres
