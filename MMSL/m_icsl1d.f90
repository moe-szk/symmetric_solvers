subroutine m_icsl1d(n,n_nzr,nzu,nzl,iu,il,uu,ul,ud,b,x,kus,kls)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n,n_nzr
  integer,intent(in)::nzu(n),nzl(n),il(n_nzr),iu(n_nzr),kus(n),kls(n)
  real(8),intent(in)::uu(n_nzr),ul(n_nzr),ud(n),b(n)
  real(8),intent(inout)::x(n)
  real(8),parameter::zero=0.d0
  integer::k,mu,muak,nu,nuak
  real(8)::t
  !------------------------------------------------------------
  !
  ! forward substitution
  !
  !------------------------------------------------------------
  do k=1,n
     t = b(k)
     do mu=1,nzl(k)
        muak = kls(k)+mu
        t = t - ul(muak)*x(il(muak))
     end do
     x(k)=ud(k)*t
  end do
  !------------------------------------------------------------
  !
  ! backward substitution
  !
  !------------------------------------------------------------
  do k=n,1,-1
     t = zero
     do nu=1,nzu(k)
        nuak = kus(k)+nu
        t=t+uu(nuak)*x(iu(nuak))
     end do
     x(k)=x(k)-ud(k)*t
  end do
  return
end subroutine m_icsl1d
