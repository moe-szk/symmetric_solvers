subroutine rpmr1d(kmax,n,nzu,nzl,iu,il,au,al,ad,ksu,ksl,b,x,eps,itmax,iopt,&
     uu,ul,ud)
    implicit real(8)(a-h,o-z)
    integer,intent(in)::kmax,n
    real(8),intent(in)::au(kmax),ad(n),al(kmax),ksu(n),ksl(n)
    integer,intent(in)::iu(kmax),il(kmax),nzu(n),nzl(n)
    real(8),intent(in)::uu(kmax),ul(kmax),ud(n)
    real(8),intent(in)::eps,b(n)
    integer,intent(in)::itmax,iopt
    real(8),intent(inout)::x(n)
    real(8)::r,v
    real(8),dimension(n)::v0vec,v1vec,v2vec,w0vec,w1vec,w2vec,z1vec,z2vec
    real(8)::gamma1,gamma2,eta,s0,s1,c0,c1,error0,error1, &
              alpha0,alpha1,alpha2,alpha3
    real(8),dimension(n)::kx,ax

    !
    ! Set,v0,w0,w1
    !
    v0vec = 0.d0
    w0vec = 0.d0
    w1vec = 0.d0
    !
    ! Set v1,z1
    !
    call prod1d(n,nzu,nzl,iu,il,au,al,ad,x,ax,int(ksu),int(ksl))
    v1vec = b - ax
    call icsl1d(n,nzu,nzl,iu,il,uu,ul,ud,ax,kx,int(ksu),int(ksl))
    z1vec = kx
    !
    ! gamma,v1,z1
    !
    gamma1 = sqrt(dot_product(z1vec,v1vec))
    v1vec = v1vec / gamma1
    z1vec = z1vec / gamma1

    psi 
    
    return 
  end subroutine rpmr1d
