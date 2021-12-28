!------------------------------------------------------------
! vpcgminr1d
!
! Variable preconditioning CG method
! (innder preconditioner:MINRES)
!
!------------------------------------------------------------
subroutine vpcgminr1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt)
  implicit real(8) (a-h,o-z)
  integer,intent(in) ::n,n_nzr
  integer,intent(in)::nzu(n),nzl(n),iu(n_nzr),il(n_nzr)
  real(8),intent(in)::au(n_nzr),al(n_nzr),ad(n),ksu(n),ksl(n)
  real(8),intent(in) ::b(n)
  real(8),intent(out)::x(n)
  real(8),parameter  ::zero=0.d0,one=1.d0,pi=4.d0*atan(1.d0)
  real(8),dimension(n)::r,p,ru,ap
  !
  eps1 = 1.0e-5
  !
  call prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,x,ap,int(ksu),int(ksl))
  !
  res2 = zero
  do i=1,n
     r(i) = b(i)-ap(i)
     res2 = res2+r(i)**2
  end do
  k = 0
  
  !
  !MINRES
  !
  nmax  = 10000
  !omega = 1.8d0
  !call sor_fast(r,ru,n,n_nzr,nzu,nzl,iu,il,au,al,ad,int(ksu),int(ksl),nmax,omega,eps1)
  call icmr1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,r,ru,eps1,nmax,0)
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
     call prod1d(n,n_nzr,nzu,nzl,iu,il,au,al,ad,p,ap,int(ksu),int(ksl))
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
     error1 = sqrt(res2/b2)

     if(iopt/=0)then
        if(iopt==-1)then
           write(*,*)k,',',error1-error0
        else if(mod(k,iopt)==0)then
           write(*,*)k,',',error1
        end if
     end if
     
     if(error1 <= eps) then
        if(iopt/=0) write(*,*) 'Number of Iteration:',k
        return
     else if(k == itmax) then
        stop '** does not converge **'
     end if

     !
     ! MINRES
     !
     !call sor_fast(r,ru,n,n_nzr,nzu,nzl,iu,il,au,al,ad,int(ksu),int(ksl),nmax,omega,eps1)
     call icmr1d(n_nzr,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,r,ru,eps1,nmax,0)
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
     
     error0 = error1
  end do
  return
end subroutine vpcgminr1d
