  !----------------------------------------------------------------------
  !     icdc1d
  !----------------------------------------------------------------------
  subroutine icdc1d(n,nzl,il,al,ad,ul,ud,kls)
    implicit real(8)(a-h,o-z)
    dimension nzl(*),il(*)
    dimension al(*),ad(*),ul(*),ud(*)
    dimension kls(*)
    
    data one/1.0/
    
    eps = 1.0e-6
    
    do k=1,n
       
       do mu=1,nzl(k)
          muak     = kls(k)+mu
          i        = il(muak)
          mui      = 1
          muk      = 1
          ul(muak) = al(muak)
          
10        continue
          
          muiai = kls(i)+mui
          mukak = kls(k)+muk
          
          if(mui > nzl(i)) then
             iliai = 0
          else
             iliai = il(muiai)
          end if
          
          if(muk > nzl(k)) then
             ilkak = 0
          else
             ilkak = il(mukak)
          end if
          
          if(iliai > ilkak) then
             muk = muk+1
          else if(iliai < ilkak) then
             mui = mui+1
          else
             ul(muak) = ul(muak)-ud(ilkak)*ul(muiai)*ul(mukak)
             mui      = mui+1
             muk      = muk+1
          end if
          
          if((mui <= nzl(i)) .and. (muk <= nzl(k))) go to 10
          
       end do
       
       t = ad(k)
       
       do mu=1,nzl(k)
          muak = kls(k)+mu
          t = t-ud(il(muak))*ul(muak)**2
       end do
       
       if(abs(t) <= eps) t = eps
       
       ud(k) = one/t
    end do
    
    return
  end subroutine icdc1d
  
