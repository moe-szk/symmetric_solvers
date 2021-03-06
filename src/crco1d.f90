subroutine crco1d (n,nz1,nz2,i1,i2,a1,a2,iv,k1s,k2s,ktot)
    implicit real(8)(a-h,o-z)
    dimension nz1(*),nz2(*),i1(*),i2(*)
    dimension a1(*),a2(*),iv(*),k1s(*),k2s(*)

    !
    ! zero clear
    !
    do i=1,n 
        nz2(i) = 0
    end do

    do k=1,ktot 
        i2(k) = 0
        a2(k) = 0
    end do

    !
    ! determination of array k1s
    !
    do i=1,n 
        do nu=1,nz1(i)
            k = k1s(i) + nu
            iv(il1(k)) = nz2(i1(k)) + 1
        end do 
        do nu=1,nz1(i) 
            nui = k1s(i) + nu 
            i1n = i1(nuii)
            ivi = k2s(i1n) + iv(i1n)
            nz2(i1n) = iv(i1n)
            i2(ivi) = i
            a2(ivi) = a1(nui)
        end do
    end do
    return 
end subroutine crco1d
