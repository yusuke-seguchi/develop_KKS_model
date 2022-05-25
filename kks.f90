! ===================================================
! PHASE-FIELD CALCULATION PROGRAM FOR BINARY ALLOY
! FOR FREE DNDRITE ANALYSIS
! 大出真知子、鈴木俊夫、”フェーズフィールドシミュレーション入門”、鋳造工学　第73巻（2001）第5号
! ===================================================

module COMMON_V
    integer :: l=0, i=1, j=1, m, n, mc, nc
    integer :: triangle=15, modsave=1, lsave=0
    real(8) :: xm, ep, ep2, W
    real(8) :: dx, dy, dt, dt1, ds, dl, dl2, ds2
    real(8) :: cle, cse, tmpmelt, tmp, vm, c0
    real(8) :: v, yk, sigma, xme, ke, beta, comnoise
    double precision :: test;    ! variable for debag

    double precision, allocatable :: phi(:,:,:), com(:,:,:)
    double precision, allocatable :: fcl(:,:), fcc(:,:), cl(:,:), cs(:,:)

    real,parameter :: r = 8.31451

end module COMMON_V

! ===================================================
! MAIN PROGRAM
! ===================================================
program main
    use COMMON_V
    implicit none

    double precision :: pd, phix, phiy, phixx, phiyy, phixy
    real(8) :: e1, e2, e3, e4, e5, th, eta, deta
    real(8) :: p, pp, pg, gd, gg, fp, dc, c, dphi
    real(8) :: fccw, fcce, fccn, fccs, fccl
    real(8) :: xj1, xj2, xj3, xj4
    real(8) :: d1, d2, d3, d4, d5
    real(8) :: a, aa, bb, cc


! READ CALCULATION CONDITION
    call cal_cond

! SET THE INITIAL PHASE-CONDITION
    call init_cond
    call outsave
    call outsave_

    allocate( cs(0:m+1, 0:n+1) ); allocate( cl(0:m+1, 0:n+1) )
    allocate( fcl(0:m+1, 0:n+1) ); allocate( fcc(0:m+1, 0:n+1) )

! SET THE PHASE-FIELD PARAMETERS (&TIMESTEP)
    ep = sqrt(18./2.2*dx*sigma); ep2 = ep*ep
    W = 2.2*sigma/dx
    call mobilitiy
    a = cle/cse * (1.-cse)/(1.-cle)

    write(6,*)'SET ALL THE CALCULATION CONDITIONS'
    write(6,*)'NOW CALCULATING.....'

! ===================================================
! CALCULATE THE GOVERNING EQUATIONS
! ===================================================
500 l = l + 1

! CALCULATION CS & CL
    do i = 0, mc + 1
        do j = 0, nc + 1
            p = phi(0,i,j)
            pp = p**3 * (10. - 15.*p + 6.*p*p) 

            if (phi(0,i,j).lt.0.001) then
                cl(i,j) = com(0,i,j)
                cs(i,j) = cl(i,j)/(a + (1.-a)*cl(i,j))

            else if (phi(0,i,j).gt.0.999) then
                cs(i,j) = com(0,i,j)
                cl(i,j) = a*cs(i,j)/(1. + (a-1.)*cs(i,j))
            
            else
                aa = pp * (1.-a)
                bb = pp + (1.-pp)*a + com(0,i,j)*(1.-a)
                cc = com(0,i,j)

                cs(i,j) = (bb - sqrt(bb*bb - 4.*aa*cc))/(2.*aa)
                cl(i,j) = (cc - pp*cs(i,j))/(1.-pp)
            end if

            fcl(i,j) = r*tmp/vm*log( cl(i,j)/(1.-cl(i,j)) )
            fccl = r*tmp/vm/(cl(i,j)*(1.-cl(i,j)))
            fccs = r*tmp/vm/(cs(i,j)*(1.-cs(i,j)))
            fcc(i,j) = fccl*fccs/((1.-pp)*fccs + pp*fccl)

        end do
    end do

    print '(i10.5, e12.5)', lsave, com(0, 10, 10)
    print '(i10.5, e12.5)', lsave, com(1, 10, 10)

! GOVERNING EQUATIONS
    do i = 1, mc
        do j = 1, nc

            p = phi(0,i,j); c = com(0,i,j)

! TIME SAVING
            pd = (phi(0,i+1,j) + phi(0,i-1,j) + phi(0,i,j+1) + phi(0,i,j-1))/4.
            if (pd.le.1.e-5) then
                dphi = 0.
                dc = dl*( (com(0,i,j+1) + com(0,i,j-1) - 2.*c)/(dy*dy) + (com(0,i+1,j) + com(0,i-1,j) - 2.*c)/(dx*dx) )
            
            else if (pd.ge.(1.-1.e-5)) then
                dphi = 0.
                dc = ds*( (com(0,i,j+1) + com(0,i,j-1) - 2.*c)/(dy*dy) + (com(0,i+1,j) + com(0,i-1,j) - 2.*c)/(dx*dx) )

! NON-TIME SAVING (140行目)
            else
                pg = 30.*p*p*(1-p)*(1-p)
                gd = 2.*p*(1.-p)*(1.-2.*p)

                gg = pg*log( (1.-cse)/(1.-cle) * (1.-cl(i,j))/(1.-cs(i,j)) )

                fp = r*tmp/vm*gg - W*gd

            phix = (phi(0,i-1,j) - phi(0,i+1,j))/(2.*dx)
            phiy = (phi(0,i,j-1) - phi(0,i,j+1))/(2.*dy)

            ! test = phi(0,100,j) - phi(0,101,j)
            ! print '(i5,i5,e12.5)', i, j, test
        

            phixx = (phi(0,i-1,j) + phi(0,i+1,j) - 2.*p)/(dx*dx)
            phiyy = (phi(0,i,j-1) + phi(0,i,j+1) - 2.*p)/(dy*dy)
            phixy = (phi(0,i+1,j+1) + phi(0,i-1,j-1) - phi(0,i-1,j+1) - phi(0,i-1,j+1))/(2.*dx*2.*dy)
            th = atan(phiy/(phix + 1.e-20))

            eta = 1. + v*cos(yk*th)
            deta = v*yk*sin(yk*th)

            e1 = ep2*eta*eta*(phixx + phiyy)
            e2 = ep2*eta*(-deta)*(sin(2.*th))*(phiyy - phixx) + 2.*cos(2.*th)*phixy
            e3 = 0.5*ep2
            e4 = deta*deta+eta*(-v*yk*yk*cos(yk*th))
            e5 = 2.*sin(2.*th)*phixy - phixx - phiyy - cos(2.*th)*(phiyy-phixx)

            dphi = xm*(e1+e2-e3*e4*e5+fp)

            d1 = ds; d2 = ds; d3 = ds; d4 = ds; d5 = ds

            if (p.le.0.9) d1 = dl
            if (phi(0,i-1,j).le.0.9) d2 = dl; if (phi(0,i+1,j).le.0.9) d3 = dl
            if (phi(0,i,j+1).le.0.9) d4 = dl; if (phi(0,i,j-1).le.0.9) d5 = dl

            fccw = 2.*d1/fcc(i,j)* d2/fcc(i-1,j) / (d1/fcc(i,j) + d2/fcc(i-1,j))
            fcce = 2.*d1/fcc(i,j)* d3/fcc(i+1,j) / (d1/fcc(i,j) + d3/fcc(i+1,j))
            fccs = 2.*d1/fcc(i,j)* d4/fcc(i,j+1) / (d1/fcc(i,j) + d4/fcc(i,j+1))
            fccn = 2.*d1/fcc(i,j)* d5/fcc(i,j-1) / (d1/fcc(i,j) + d5/fcc(i,j-1))

            xj1 = (-fcl(i,j) + fcl(i-1,j))/dx * fccw
            xj2 = (-fcl(i,j) + fcl(i+1,j))/dx * fcce
            xj3 = (-fcl(i,j) + fcl(i,j+1))/dy * fccs
            xj4 = (-fcl(i,j) + fcl(i,j-1))/dy * fccn

            dc = (xj1 + xj2)/dx + (xj3 + xj4)/dy

            end if

        phi(1,i,j) = p + dphi*dt; com(1,i,j) = c + dc*dt
        end do
    end do

! END GOVERNING EQUAITION CALCULATIONS

! BOUNDARY CONDITION
    do i=0, m+1
        phi(1,i,0) = phi(1,i,1); phi(1,i,n+1) = phi(1,i,n)
        com(1,i,0) = com(1,i,1); com(1,i,n+1) = com(1,i,n)
    end do

    do j=0, n+1
        phi(1,0,j) = phi(1,1,j); phi(1,m+1,j) = phi(1,m,j)
        com(1,0,j) = com(1,1,j); com(1,m+1,j) = com(1,m,j)
    end do

! RENEWAL OF PHASE & CONCENTRATION FIELDS
    do i=0,m+1
        do j=0,n+1
            phi(0,i,j) = phi(1,i,j); com(0,i,j) = com(1,i,j)
        end do
    end do
    
! NOISE
    call noise

! AREA SET FOR TIME SAVING
    if (l.ge.100) call areaset

! OUTPUT
    lsave = lsave + 1
    if (mod(l,modsave).eq.0) call outsave
    if (mod(l,modsave).eq.0) call outsave_

! END CONDITION
    if (phi(0,1,n-10).le.0.5) goto 500
    write(6,*)'CALCULATION HAS FINISHED'
    
end program main

! ===================================================
! SUBRUTINE
! ===================================================
! READ CALCULATION CONDITION
subroutine cal_cond
    use COMMON_V
    implicit none

    m = 250                 ! x-direction mesh number
    n = 250                 ! y-direction mesh number
    dx = 1.e-8              ! mesh size
    dy = 1.e-8              ! mesh size
    dl = 3.e-9              ! dl
    ds = 3.e-13             ! ds
    c0 = 0.0196             ! initial solute content
    tmp = 900.0             ! initial temperature
    tmpmelt = 933.3        ! melting point
    vm = 10.547e-6          ! moler volume
    xme = 640.0             ! liquidus slope
    ke = 0.14               ! partition coefficient
    beta = 0.0              ! kinetic coefficient
    v = 0.03                ! anisotropy ep=ep(1+v*cos(yk*th))
    yk = 4.0                ! anisotropy (yk-fold)
    sigma = 0.093           ! interface energy
    comnoise = 0.1          ! noise

    nc = n; mc = m; dl2 = 2.0*dl/dx; ds2 = 2.0*ds/dx
    cle = (tmpmelt - tmp)/xme; cse = cle*ke
    
end subroutine cal_cond

! INIT_COND
subroutine init_cond
    use COMMON_V
    implicit none
    
    allocate( phi(0:1, 0:m+1, 0:n+1) )
    allocate( com(0:1, 0:m+1, 0:n+1) )

    do i = 1, m
        do j = 1, n
            phi(0,i,j) = 0.; phi(1,i,j) = 0.
            com(0,i,j) = c0; com(1,i,j) = c0
        end do
    end do

    do i = 1, triangle+5
        do j = 1, triangle+5
            if (j.lt.(-i+triangle)) then
                phi(0,i,j) = 1.; phi(1,i,j) = 1.
                com(0,i,j) = cse; com(1,i,j) = cse
            end if

            if (j.eq.(-i+triangle)) then
                phi(0,i,j) = 0.5; phi(1,i,j) = 0.5
            end if

        end do
    end do

    do i = 0, m+1
        phi(0,i,0) = phi(0,i,1); phi(0,i,n+1) = phi(0,i,n)
        com(0,i,0) = com(0,i,1); com(0,i,n+1) = com(0,i,n)
    end do

    do j = 0, n+1
        phi(0,0,j) = phi(0,1,j); phi(0,m+1,j) = phi(0,m,j)
        com(0,0,j) = com(0,1,j); com(0,m+1,j) = com(0,m,j)
    end do

    return
end subroutine init_cond

! PF MOBILITY
subroutine mobilitiy
    use COMMON_V
    implicit none

    integer :: a
    real(8) :: p1, p2, pp1, pp2, fun1, fun2, zeta, fccle, fccse, alpha

    ep2 = ep*ep; zeta = 0.

    fccle = r*tmp/vm/(cle*(1.-cle))
    fccse = r*tmp/vm/(cse*(1.-cse))

    do a = 1, 998
        p1 = a/1000.
        p2 = p1 + 0.001
        pp1 = p1**3 * (10. - 15.*p1 + 6.*p1*p1)
        pp2 = p2**3 * (10. - 15.*p2 + 6.*p2*p2)
        fun1 = pp1*(1.-pp1) / ((1.-pp1)*fccse + pp1*fccle) / (p1*(1.-p1))
        fun2 = pp2*(1.-pp2) / ((1.-pp2)*fccse + pp1*fccle) / (p2*(1.-p2)) !　pp1*fccleがあやしい。
        zeta = zeta + (fun1 + fun2)*0.001/2.
    end do

    alpha = beta*r*tmp*(1.-ke) / (vm*xme)
    xm = 1./(ep*ep/sigma * (alpha+ep/(dl*sqrt(2.*W))*zeta*fccse*fccle*(cle-cse)**2))
    print *, xm

    dt = dx**2 / (5.*xm*ep**2)
    dt1 = dx**2 / (5.*dl)
    dt = dmin1(dt, dt1)

    return
end subroutine mobilitiy


! COMNOISE
subroutine noise
    use COMMON_V
    implicit none

    real(8) :: countnoise, comtot, comdam, cnoise
    integer :: lam = 12869, c = 6925, mmu = 32768, x = 19724

    countnoise = 0.; comtot = 0.

    do i = 1, mc
        do j = 1, nc
            if (phi(0,i,j).gt.0.01.and.phi(0,i,j).le.0.5) then
                comdam = com(0,i,j)
                
                x = mod((x*lam+c), mmu)
                cnoise = (real(x)/mmu - 0.5) * comnoise

                com(0,i,j) = comdam*(1. + cnoise)
                comtot = comtot + comdam*cnoise
                countnoise = countnoise + 1
            end if
        end do
    end do

    do i = 1, mc
        do j = 1, nc
            if (phi(0,i,j).gt.0.01.and.phi(0,i,j).le.0.5) then
                com(0,i,j) = com(0,i,j) - (comtot/countnoise)
            end if
        end do
    end do

    do i = 0, m+1
        com(0,i,0) = com(0,i,1); com(0,i,n+1) = com(0,i,n)
    end do

    do j = 0, n+1
        com(0,0,j) = com(0,1,j); com(0,m+1,j) = com(0,m,j)
    end do
    
    return
end subroutine noise

! AREA SET 
subroutine areaset
    use COMMON_V
    implicit none

    do j = 1, n
        if (abs(com(0,1,j)/c0 - 1.).gt.1.e-5) nc = j
    end do

    do i = 1, m
        if (abs(com(0,i,1)/c0 - 1.).gt.1.e-5) mc = i
    end do

    nc = nc + 10; mc = mc + 10
    if (nc.gt.n) nc = n; if (mc.gt.m) mc = m

    return
end subroutine areaset

! OUTSAVE
subroutine outsave
    use COMMON_V
    implicit none

    character*3 :: out_num
    character*15 :: fpout
    integer :: one, ten, hand

    
    one = mod(lsave, 10)
    ten = mod(int(real(lsave)/10.), 10)
    hand = mod(int(real(lsave)/100.), 10)

    one = 48 + one; ten = 48 + ten; hand = 48 + hand;
    out_num = char(hand)//char(ten)//char(one)
    
    fpout = "output"//out_num//".vtk"
    open(101,file=fpout, err=1000)
    write(101,'(a)') '# vtk DataFile Version 3.0'
    write(101,'(a)') 'output.vtk'
    write(101,'(a)') 'ASCII'
    write(101,'(a)') 'DATASET STRUCTURED_POINTS'
    write(101,'(a,3i5)') 'DIMENSIONS',m,n,1
    write(101,'(a,3f4.1)')'ORIGIN' ,0.0,0.0,0.0
    write(101,'(a,3i2)')'ASPECT_RATIO',1,1,1
    write(101,'(a,1i11)')'POINT_DATA',m*n*1
    write(101,'(a)')'SCALARS concentration double'
    write(101,'(a)')'LOOKUP_TABLE default'

    do j=1,n
        do i=1,m
            write(101,300) com(0,i,j)
        end do
    end do

    write(101,'(a)')'SCALARS phase_field double'
    write(101,'(a)')'LOOKUP_TABLE default'

    do j=1,n
        do i=1,m
            write(101,300) phi(0,i,j)
        end do
    end do

    close(101)
    return

300     format(e12.5)
1000    write(6,*)'ERROR IN FILE OPEN'

end subroutine outsave

subroutine outsave_
    use COMMON_V
    implicit none

    character*3 :: out_num
    character*15 :: fpout
    integer :: one, ten, hand

    one = mod(lsave, 10)
    ten = mod(int(real(lsave)/10.), 10)
    hand = mod(int(real(lsave)/100.), 10)

    one = 48 + one; ten = 48 + ten; hand = 48 + hand;
    out_num = char(hand)//char(ten)//char(one)
    
    fpout = "output_"//out_num//".vtk"
    open(101,file=fpout, err=1000)
    write(101,'(a)') '# vtk DataFile Version 3.0'
    write(101,'(a)') 'output.vtk'
    write(101,'(a)') 'ASCII'
    write(101,'(a)') 'DATASET STRUCTURED_POINTS'
    write(101,'(a,3i5)') 'DIMENSIONS',mc,nc,1
    write(101,'(a,3f4.1)')'ORIGIN' ,0.0,0.0,0.0
    write(101,'(a,3i2)')'ASPECT_RATIO',1,1,1
    write(101,'(a,1i11)')'POINT_DATA',mc*nc*1
    write(101,'(a)')'SCALARS concentration double'
    write(101,'(a)')'LOOKUP_TABLE default'

    do j=1,n
        do i=1,m
            write(101,300) com(1,i,j)
        end do
    end do

    write(101,'(a)')'SCALARS phase_field double'
    write(101,'(a)')'LOOKUP_TABLE default'

    do j=1,n
        do i=1,m
            write(101,300) phi(1,i,j)
        end do
    end do

    close(101)
    return

! 410     format(i4, i4, i10)
300     format(e12.5)
1000    write(6,*)'ERROR IN FILE OPEN'

end subroutine outsave_


 





