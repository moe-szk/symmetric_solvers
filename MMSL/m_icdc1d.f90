subroutine m_icdc1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,ksu,ksl,uu,ul,ud)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n,n_nzr
  integer,intent(in)::nzu(n),nzl(n),iu(n_nzr),il(n_nzr)
  real(8),intent(in)::au(n_nzr),al(n_nzr),ad(n),ksu(n),ksl(n)
  real(8),intent(out)::uu(n_nzr),ul(n_nzr),ud(n)
  real(8),allocatable::kus(:),kls(:)
  integer,allocatable::ivu(:),ivl(:)
  
  allocate(kus(n),kls(n),ivu(n_nzr),ivl(n_nzr))
  
  call cnonzr(n_nzr,n,iu,nzu, kus,nzl,kls,ktot)
  call crco1d(n,nzu,nzl,iu,il,au,al,ivu,kus,kls,ktot)
  call icdc1d(n,nzl,il,al,ad,ul,ud,kls)
  call crco1d(n,nzl,nzu,il,iu,ul,uu,ivl,kls,kus,ktot)
  
  deallocate(kus,kls,ivu,ivl)
  return
end subroutine m_icdc1d
