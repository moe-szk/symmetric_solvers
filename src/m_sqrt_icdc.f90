subroutine m_sqrt_icdc(n,adim,udim)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n
  real(8),intent(in)::adim(n,n)
  real(8),intent(out)::udim(n,n)
  udim = adim
  do k=1,n
     udim(k,k) = sqrt(adim(k,k))
     w = 1.d0 / udim(k,k)

     do i=k+1,n
        if(adim(i,k)/=0.d0)then
           udim(i,k) = udim(i,k)/w
        end if
     end do

     do j=k+1,n
        if(adim(j,k)/=0.d0)then
           do i=j,n
              if(adim(i,j)/=0.d0 .and. adim(i,j)/=0.d0)then
                 udim(i,j) = udim(i,j) - udim(i,k)*udim(j,k)
              end if
           end do
        end if
     end do
  end do

  do i=1,n
     do j=i+1,n
        udim(i,j) = udim(j,i)
     end do
  end do
end subroutine m_sqrt_icdc
