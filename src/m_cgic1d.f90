!-------------------------------------------------------------------------------
! m_cgic
!
! conjugate gradient method for iccg
!
! by Moeto Suzuki
!-------------------------------------------------------------------------------
subroutine m_cgic1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt, &
     uu,ul,ud)
  implicit real(8)(a-h,o-z)
  integer,intent(in)::n_nzr,n
  real(8),intent(in)::au(n_nzr),ad(n),al(n_nzr),ksu(n),ksl(n)
  integer,intent(in)::iu(n_nzr),il(n_nzr),nzu(n),nzl(n)
  real(8),intent(in)::uu(n_nzr),ul(n_nzr),ud(n)
  real(8),intent(in)::eps,b(n)
  integer,intent(in)::itmax,iopt
  real(8),intent(inout)::x(n)
  real(8),dimension(n)::ap,r,ru,p
  real(8)::error0,error1
  !
  ! Initialization
  !
  error0 = 0.d0
  call m_prod1d(n,nzu,nzl,iu,il,au,al,ad,x,ap,int(ksu),int(ksl))

  res2 = 0.d0
  do i=1,n
     r(i) = b(i) - ap(i)
     res2 = res2 + r(i)**2
  end do
  k = 0

  call m_icsl1d(n,nzu,nzl,iu,il,uu,ul,ud,r,ru,int(ksu),int(ksl))

  rur0 = 0.d0

  do i=1,n
     p(i) = ru(i)
     rur0 = rur0 + r(i)*ru(i)
  end do

  b2 = 0.d0
  do i=1,n
     b2 = b2 + b(i)**2
  end do

  !
  ! Iterative loop
  !

  do k=1,itmax
     
     call m_prod1d(n,nzu,nzl,iu,il,au,al,ad,p,ap,int(ksu),int(ksl))
     
     pap = 0.d0
     do i=1,n
        pap = pap + p(i) * ap(i)
     end do

     alpha = rur0 /pap

     res2 = 0.d0

     do i=1,n
        x(i) = x(i) + alpha*p(i)
        r(i) = r(i) - alpha*ap(i)
        res2 = res2 + r(i)**2
     end do

     error1 = sqrt(res2/b2)
     
     if(iopt/=0) then
        if(iopt==-1)then
           write(*,*)k,',',error1 - error0
        else if(mod(k,iopt)==0)then
           write(*,*)k,',',error1
        end if
     end if
     
     !
     ! Convergence check
     !
     if(error1 < eps) then
        if(iopt/=0) write(*,*)'Number of Iteration:',k
        return
     end if
     if(k==itmax) write(*,*)'** does not convergence **'
     
     call m_icsl1d(n,nzu,nzl,iu,il,uu,ul,ud,r,ru,int(ksu),int(ksl))
     
     rur1 = 0.d0

     do i=1,n
        rur1 = rur1 + r(i)* ru(i)
     end do

     beta = rur1 / rur0
     rur0 = rur1

     do i=1,n
        p(i) = ru(i) + beta*p(i)
     end do
     error0 = error1
  end do
  return
end subroutine m_cgic1d
