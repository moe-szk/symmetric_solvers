subroutine cr1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n,n_nzr
  real(8),intent(in)::au(n_nzr),ad(n),al(n_nzr),ksu(n),ksl(n)
  integer,intent(in)::iu(n_nzr),il(n_nzr),nzu(n),nzl(n)
  real(8),intent(in)::eps,b(n)
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)
  real(8)::r1(n),r0(n),ap(n),ap1(n),ar(n),alpha,p(n)
  
  call m_prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,x,ap,int(ksu),int(ksl))
  r0 = b - ap
  p = r0
  
  do j=1,itmax
     call m_prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,r0,ar,int(ksu),int(ksl))
     call m_prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,p,ap,int(ksu),int(ksl))
     
     alpha = dot_product(ar,r0)/dot_product(ap,ap)
     
     x = x + alpha * p
     
     r1 = r0 - alpha * ap

     call m_prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,r1,ap1,int(ksu),int(ksl))
     
     beta = dot_product(ap1,r1)/dot_product(ar,r0)

     p = r1 + beta * p
     
     error = sqrt(sum(r1**2))
     
     if(iopt/=0)write(*,*)j,',',error
     if(error<eps) exit
     if(j==itmax)write(*,*)'** does not converge. **'
     r0 = r1
  end do
end subroutine cr1d
