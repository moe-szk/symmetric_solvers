!------------------------------------------------------------
! vpcgsor1d
!
! Variable preconditioning CG method
! (innder preconditioner:SOR)
!
!programmed by Moeto Suzuki
!
!------------------------------------------------------------
subroutine vpcgsor1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt)
  implicit real(8) (a-h,o-z)
  integer,intent(in) ::n,n_nzr
  integer,intent(in)::nzu(n),nzl(n),iu(n_nzr),il(n_nzr)
  real(8),intent(in)::au(n_nzr),al(n_nzr),ad(n),ksu(n),ksl(n)
  real(8),intent(in) ::b(n)
  real(8),intent(out)::x(n)
  real(8),parameter  ::zero=0.d0,one=1.d0,pi=4.d0*atan(1.d0)
  real(8),dimension(n)::r,p,ru,ap
  !
  eps1 = 1.0e-3
  !
  call prod1d(n,nzu,nzl,iu,il,au,al,ad,x,ap,int(ksu),int(ksl))
  !
  res2 = zero
  do i=1,n
     r(i) = b(i)-ap(i)
     res2 = res2+r(i)**2
  end do
  k = 0
  
  !
  !SOR
  !
  nmax  = 100
  omega = 1.8d0
  call icsor1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,r,ru,eps1,nmax,0,omega)
  !
  !
  !
  rur0  = zero
  
  do i=1,n
     p(i)=ru(i)
     rur0=rur0+r(i)*ru(i)
  end do
  b2 = zero
  do i=1,n
     b2 = b2 + b(i)**2
  end do
  !
  !iterative loop
  !
  itmax = 1000
  do k=1,itmax
     call prod1d(n,nzu,nzl,iu,il,au,al,ad,p,ap,int(ksu),int(ksl))
     pap = zero
     do i=1,n
        pap = pap+p(i)*ap(i)
     end do
     alpha = rur0/pap
     res2  = zero

     do i=1,n
        x(i) = x(i) + alpha * p(i)
        r(i)     = r(i) - alpha * ap(i)
        res2     = res2 + r(i)**2
     end do
     !
     if(sqrt(res2/b2) <= eps2) then
        write(*,*) k,sqrt(res2/b2)
        return
     else if(k == itmax) then
        stop 'Does not converge!!!'
     end if

     write(*,*)k,sqrt(res2/b2)
     !
     !SOR
     !
     call icsor1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,r,ru,eps1,nmax,0,omega)
     !
     !
     !
     
     rur1 = zero
     do i=1,n
        rur1 = rur1 + r(i)*ru(i)
     end do
     beta = rur1/rur0
     rur0 = rur1

     do i=1,n
        p(i) = ru(i)+beta * p(i)
     end do
  end do
  
  return
  
end subroutine vpcgsor1d
