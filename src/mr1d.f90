subroutine mr1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n,n_nzr
  real(8),intent(in)::au(n_nzr),ad(n),al(n_nzr),ksu(n),ksl(n)
  integer,intent(in)::iu(n_nzr),il(n_nzr),nzu(n),nzl(n)
  real(8),intent(in)::eps,b(n)
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)
  real(8)::r(n),ap(n),alpha
  
  call prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,x,ap,int(ksu),int(ksl))
  r = b - ap

  do j=1,itmax
     call prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,r,ap,int(ksu),int(ksl))
     
     alpha = dot_product(r,ap)/dot_product(ap,ap)

     x = x + alpha * r

     r = r - alpha * ap
     
     error = sqrt(sum(r**2))
     
     if(iopt/=0)write(*,*)j,',',error
     if(error<eps) exit
     if(j==itmax)write(*,*)'** does not converge. **'
  end do
end subroutine mr1d
