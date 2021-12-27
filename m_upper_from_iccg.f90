!-------------------------------------------------------------------------------
! m_upper_from_iccg
!
! by Moeto Suzuki
!
!-------------------------------------------------------------------------------
subroutine m_upper_from_iccg(n,n_nzr,nzl,il,al,ksl,nzu,iu,au,ksu)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n,n_nzr
  integer,intent(in)::nzl(n),il(n_nzr),ksl(n)
  real(8),intent(in)::al(n_nzr)
  integer,intent(out)::nzu(n),iu(n_nzr),ksu(n)
  real(8),intent(out)::au(n_nzr)
  integer::ku,kl
  ku = 0
  ksu(1) = 0
  
  do i=1,n

     nzu(i) = 0
     ksu(i+1) =0
     do kl=1,n_nzr
        if(il(kl)==i)then
           ku = ku + 1
           nzu(i) = nzu(i) + 1
           if(i/=n) ksu(i+1) = ksu(i)+ksu(i+1)+1
           au(ku) = al(kl)
           do j=1,n
              if(ksl(j)<ku) iu(ku) = j
           end do
        end if
     end do

  end do
  
  if(ku/=n_nzr)then
     write(*,*)"**** Error in m_upper_from_iccg ****"
     write(*,*)"n_nzr : ",n_nzr,"sum of nzu",ku
     stop
  end if

  return
end subroutine m_upper_from_iccg
