!------------------------------------------------------------
! sor1d
! SOR
!
!------------------------------------------------------------
subroutine sor1d(kmax,nval,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt,omg)
  implicit real(8) (a-h,o-z)
  integer,intent(in)   ::nval,kmax
  integer,intent(in)   ::nzu(nval),nzl(nval),iu(kmax),il(kmax),itmax
  integer,intent(in)   ::ksl(nval),ksu(nval)
  real(8),intent(in)   ::b(nval),au(kmax),al(kmax),eps,ad(nval),omg
  real(8),intent(inout)::x(nval)
  real(8),parameter    ::zero=0.d0
  real(8)::xold(nval),error,xber
  xold = x
  
  do n=1,itmax
     
     do i=1,nval
        
        sum1 = zero
        do j=1,nzu(i)
           mu   = ksu(i) + j
           sum1 = sum1 + (au(mu) * x(iu(mu)))
        end do
        
        sum2 = zero
        do j=1,nzl(i)
           ml   = ksl(i) + j
           sum2 = sum2 + (al(ml) * xold(il(ml)))
        end do

        xber = (b(i) - (sum1 + sum2)) / ad(i)

        x(i) = xold(i) + omg * (xber - xold(i))
        
     end do

     error = maxval(abs(x-xold))/maxval(abs(x))
     
     if(iopt/=0)then
        if(mod(n,iopt)==0) write(*,*)n,',',error
     end if
     if(isnan(error))stop
     if (error < eps)then
        return
     end if

     
     xold = x
  end do
  return
end subroutine sor1d
