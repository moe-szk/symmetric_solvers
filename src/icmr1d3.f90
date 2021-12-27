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
subroutine icmr1d3(kmax,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::kmax,n
  real(8),intent(in)::au(kmax),ad(n),al(kmax),ksu(n),ksl(n)
  integer,intent(in)::iu(kmax),il(kmax),nzu(n),nzl(n)
  real(8),intent(in)::eps,b(n)
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)

  real(8),allocatable::uu(:),ul(:),ud(:),kus(:),kls(:)
  integer,allocatable::ivu(:),ivl(:),mzu(:),mzl(:)

  !
  real(8),dimension(n,n)::adim,mdim
  call k_convm_from_iccg(n,kmax,nzu,iu,au,ad,adim)
  call m_icdec(n,adim,mdim)
  allocate(mzu(n),mzl(n),ivu(kmax),ivl(kmax),&
       uu(kmax),ul(kmax),ud(n),kus(n),kls(n))
  call m_conv_2_iccg(n,kmax,mdim,mzu,mzl,ivu,ivl,uu,ul,ud,kus,kls)

  call rpmr1d(kmax,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt,&
       uu,ul,ud)

  deallocate(uu,ul,ud,ivu,ivl,kus,kls)
end subroutine icmr1d3
