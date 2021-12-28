subroutine prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,x,ap,ksu,ksl)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n,nzu(n),nzl(n),il(n_nzr),iu(n_nzr),ksu(n),ksl(n)
  real(8),intent(in)::au(n_nzr),al(n_nzr),ad(n),x(n)
  real(8),intent(out)::ap(n)
  real(8),parameter::zero=0.d0
  integer::k,mu,muak,nu,nuak
  real(8)::t
  do k=1,n
     t = ad(k)* x(k)
     !
     ! upper matrix
     !
     do i=1,nzu(k)
        j = ksu(k)+i
        t = t + au(j)*x(iu(j))
     end do
     !
     ! lower matrix
     !
     do i=1,nzl(k)
        j = ksl(k)+i
        t = t + al(j)*x(il(j))
     end do
     ap(k) = t
  end do
  end subroutine prod1d
