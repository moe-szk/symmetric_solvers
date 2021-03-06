!
! cnonzr
!
! Count nonzero elements above i-th row
!
subroutine cnonzr(kmax,x,i1,nz1,ks1,nz2,ks2,ktot)
    dimension i1(*), nz1(*), ks1(*), nz2(*), ks2(*)

    !
    ! Counting nonzero elements above i-th row 
    ! in off-diagonal upper triangular matrix.
    !
    ks1(1) = 0
    do i=2,n
        im1 = i - 1
        ks1(i) = ks1(im1) + nz1(i,1)
    end do  

    !
    ! nonzero elements on i-th row
    ! in off-diagonal lower triangular matrix.
    !
    do i=1,n
        nz2(i) = 0
    end do

    do k=1,kto
        j = i1(k)
        nz2(j) = nz2(j)+1
    end do

    !
    ! nonzero elements above i-th row
    ! in off-diagonal lower triangular matrix
    !
    ks2(1) = 0
    do i=2,n 
        im1 = i-1
        ks2(i) = ks2(im1) + nz2(im1)
    end do
    return
end subroutine cnonzr
