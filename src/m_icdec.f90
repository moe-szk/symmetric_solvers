!-------------------------------------------------------------------------------
! m_icdec
!
! Incomplete Cholesky Decomposition
!
! by Moeto Suzuki
!
! arguments
! [input]
! n    : size of a matrix 
! amat : symmetric matrix
!
! [output]
! mmat : incomplete Cholesky decomposed matrix
!        (M = LLt)
!-------------------------------------------------------------------------------
subroutine m_icdec(n,amat,mmat)
  integer,intent(in )::n
  real(8),intent(in )::amat(n,n)
  real(8),intent(out)::mmat(n,n)
  Mmat = Amat
  do i=1,n
     do j=i+1,n
        Mmat(i,j) = 0.d0
     end do
  end do

  do k=1,n
     Mmat(k,k) = sqrt(abs(Mmat(k,k)))
     do i=k+1,n
        if(Mmat(i,k)/=0.d0) then
           Mmat(i,k) = Mmat(i,k) / Mmat(k,k)
        end if
     end do
     do j=k+1,n
        do i=j,n
           if((Mmat(i,j)/=0.d0) .and. (Mmat(i,k)/=0.d0))then
              Mmat(i,j) = Mmat(i,j) - Mmat(i,k) * Mmat(j,k)
           end if
        end do
     end do
  end do
  
  do i=1,n
     do j=i+1,n
        mmat(i,j) = mmat(j,i)
     end do
  end do
  return
end subroutine m_icdec
