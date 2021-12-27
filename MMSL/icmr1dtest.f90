!-------------------------------------------------------------------------------
! icmr
!
! Incomplete Cholesky decompositioned MINRES method
!
! by Moeto Suzuki
!
! arguments
!
!-------------------------------------------------------------------------------
subroutine icmr1dtest(kmax,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt,u)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::kmax,n
  real(8),intent(in)::au(kmax),ad(kmax),al(kmax),ksu(n),ksl(n)
  integer,intent(in)::iu(kmax),il(kmax),nzu(n),nzl(n)
  real(8),intent(in)::eps,b(n),u(n)
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)
  real(8)::uu(kmax),ul(kmax),ud(n)
  call m_icdc1d(n,kmax,nzu,nzl,iu,il,au,al,ad,ksu,ksl,uu,ul,ud)
  call rpmr1dtest(kmax,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt,&
       uu,ul,abs(ud),u)
  return
end subroutine icmr1dtest
