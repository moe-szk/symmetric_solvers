!-------------------------------------------------------------------------------
! m_conv_2_iccg
!
! Compression Symmetric matrix by ICCG format
! 
! by Moeto Suzuki
!
! arguments
! [input]
! n
! n_nzr
! adim(n)
!
! [output]
! nzu(n)
! nzl(n)
! iu(n_nzr)
! il(n_nzr)
! au(n_nzr)
! al(n_nzr)
! ad(n_nzr)
! ksu(n)
! ksl(n)
!-------------------------------------------------------------------------------
subroutine m_conv_2_iccg(n,n_nzr,adim,nzu,nzl,iu,il,au,al,ad,ksu,ksl)
  implicit real(8)(a-h,o-z)
  integer,intent(in )::n,n_nzr
  real(8),intent(in )::adim(n,n)
  integer,intent(out)::nzu(n),nzl(n),iu(n_nzr),il(n_nzr)
  real(8),intent(out)::au(n_nzr),al(n_nzr),ad(n),ksu(n),ksl(n)

  n_nonzero = 0
  do i=1,n 
   adim(i,i) = ad(i) 
   do j=i+1,n 
      if(adim(i,j)/=0.d0)then
         n_nonzero     = n_nonzero + 1
         au(n_nonzero) = adim(i,j)
         nzu(i) = nzu(i) + 1
      end if
   end do
  end do 

  if(n_nzr/=n_nonzero) then
     write(*,*) '*** Error in s_conv_2_iccg ***'
     write(*,*) 'argument:n_nzr is wrong.'
     write(*,*) 'n_nzr :',n_nzr,'n_nonzero',n_nonzero
     stop
  end if

  n_nonzero = 0
  do i=1,n 
   do j=1,i
      if(adim(i,j)/=0.d0)then
         n_nonzero     = n_nonzero + 1
         al(n_nonzero) = adim(i,j)
         nzl(i) = nzl(i) + 1
      end if
   end do
  end do 

  ksu(1) = 0
  ksl(1) = 0
  do i=2,n
     ksu(i) = ksu(i-1) + nzu(i-1)
     ksl(i) = ksl(i-1) + nzl(i-1)
  end do
  return
end subroutine m_conv_2_iccg
