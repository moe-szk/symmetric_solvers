subroutine m_sqrt_icdc1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,ksu,ksl,uu,ul,ud)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n,n_nzr
  integer,intent(in)::nzu(n),nzl(n),iu(n_nzr),il(n_nzr)
  real(8),intent(in)::au(n_nzr),al(n_nzr),ad(n),ksu(n),ksl(n)
  real(8),intent(out)::uu(n_nzr),ul(n_nzr),ud(n)
  integer::x,y
  integer,allocatable::nzuu(:),iuu(:),ksuu(:)
  allocate(nzuu(n),iuu(n_nzr),ksuu(n))
  
  ul = al
  
  do k=1,n

     if(ad(k)<0.d0)then
        write(*,*)" **** Error in m_sqrt)icdc1d **** "
        stop "There is a negative number in diagonal element. "
     end if
     ud(k)=sqrt(ad(k))
     
     do i=1,n_nzr
        if(il(i)==k)then
           ul(i)=al(i)/ud(k)
        end if
     end do

     do j=k+1,n
        do i=j,n


           do y=1,n_nzr
              if(il(y)==k .and. (y<=ksl(i)))then
                 if((y<=ksl(j)))then
                    if(il(y)==j .and. y<=ksl(i))then
                       ikloc = ksl(i)+k
                       jkloc = ksl(j)+k
                       ul(y) = ul(y) - ul(ikloc)*ul(jkloc)
                       write(*,*)ul(y)
                    end if
                 end if
              end if
           end do
           
        end do
     end do
     
  end do
  
  call m_upper_from_iccg(n,n_nzr,nzl,il,ul,int(ksl),nzuu,iuu,uu,ksuu)
  deallocate(nzuu,iuu,ksuu)
  return
end subroutine m_sqrt_icdc1d
