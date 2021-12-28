subroutine gd1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n,n_nzr
  real(8),intent(in)::au(n_nzr),ad(n),al(n_nzr),ksu(n),ksl(n)
  integer,intent(in)::iu(n_nzr),il(n_nzr),nzu(n),nzl(n)
  real(8),intent(in)::eps,b(n)
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)
  real(8)::r(n),ar(n),alpha,r0
  
  call prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,x,ar,int(ksu),int(ksl))
  r = b - ar
  r0 = sum(r**2)
  
  do j=1,itmax
     
     call prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,r,ar,int(ksu),int(ksl))
     
     alpha = dot_product(r,r) / dot_product(r,ar)
     !alpha = sqrt(sum(r**2)) / dot_product(r,ar)
     !if(iopt/=0)write(*,*)alpha
     x = x + alpha * r

     r = r - alpha * ar
     
     error = sqrt(sum(r**2) / r0)
     
     if(iopt/=0)write(*,*)j,',',error
     if(error<eps .or. isnan(error)) exit
     if(j==itmax)write(*,*)'** does not converge. **'
  end do
end subroutine gd1d
