
module DBL
    implicit none
    integer, parameter :: dp=kind(0.0d0)
    real(dp) :: pi = 3.141592
end module DBL

subroutine compute_intensity(q,x,y,z,ff0,Natoms,Iq,Nq,Nr,beamprofile,Ny)
    use DBL
    implicit none
    ! if beamprofile is supplied
    !f2py integer optional, intent(in) :: Ny
    integer :: Ny
    !f2py real(dp) optional dimension(Ny,2),intent(in) :: beamprofile
    real(dp), dimension(Ny,2) :: beamprofile
    ! end optional arguments
    integer :: Natoms
    integer :: i,j,k,r_i,Nq
    real(dp), dimension(Natoms) :: x,y,z,ff0
    real(dp), dimension(Nq) :: q, Eqsq
    real(dp), dimension(Nq), intent(out) :: Iq
    real(dp) :: dist, b, dr, dq, fi, fj
    real(dp) :: rem_dr, Dmax, qr, aij, qy, dy
    real(dp), dimension(:,:), allocatable :: A  
    real(dp), dimension(1000) :: r,Pr
    integer, intent(out) :: Nr
    logical :: smear = .false.

    b  = 0.23 ! gaussian approximation of form factor
    Iq = 0
    dr = 0.5
    dq = q(2)-q(1)
    Pr = 0.0 ! initialize P(r) to zeros
    Eqsq = exp(-b*q**2) ! dummy atom, represented by a gaussian
    Dmax = 0

    do i=1,Natoms
        fi = ff0(i) 
        do j=1,i
            if (i/=j) then
                fj = ff0(j)
                dist = sqrt((x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2)
                rem_dr = mod(dist,dr)               
                r_i    = int((dist-rem_dr) / dr) + 1 ! index at distance, dist
                Pr(r_i)  = Pr(r_i) + fi*fj 
                if (dist.ge.Dmax) then
                    Dmax = dist
                end if
            end if
        end do
    end do

!     write (*,101) Dmax
! 101 format("Dmax = ",f6.2)

    ! perform fourier transform and add auto-correlation part
    Nr = int( (Dmax-mod(Dmax,dr))/dr ) + 1
    do i=1,Nr
        r(i) = (i-1)*dr
    end do

    allocate(A(Nq,Nr))

    if (Ny>2) then !Ny is flag for beam profile presence
        smear = .true.
        dy = beamprofile(2,1) - beamprofile(1,1)
        ! print *,"Calculating smearing matrix ..."
        ! do i=1,Ny
        !     write(*,'(2f10.4)') beamprofile(i,1),beamprofile(i,2)
        ! end do
        do i=1,Nq
            do j=1,Nr
                aij = 0.0
                do k=1,Ny
                    qy = sqrt(q(i)**2 + beamprofile(k,1)**2) * r(j)
                    if (qy.lt.1e-15) then
                        aij = aij + beamprofile(k,2) * 1.0 * dy
                    else
                        aij = aij + beamprofile(k,2) * (sin(qy)/qy) * dy
                    end if
                end do
                A(i,j) = 4.0*pi*2.0*aij*dr ! multiply by 2 due to beam symmetry
            end do
        end do 
    else
        print *,"No smearing, continuing ... "
        do i=1,Nq
            do j=1,Nr
                qr = q(i)*r(j)
                if (qr.lt.1e-15) then
                    A(i,j) = A(i,j) + 4.0*pi*dr
                else
                    A(i,j) = A(i,j) + 4.0*pi*(sin(qr)/qr)*dr
                end if
            end do
        end do
    end if

    ! use intrinsic function 
    Iq = matmul(A,Pr(1:Nr)) * dr

    Iq = 2*Iq ! inter atomic scaling

    if (smear) then
        ! also smear the dummy form factor
        Eqsq=0.0
        do i=1,Nq
            do j=1,Ny
                qy = sqrt(q(i)**2 + beamprofile(j,1)**2)
                Eqsq(i) = Eqsq(i) + 2*beamprofile(j,2)*exp(-b*qy*qy) * dy
            end do
        end do
    end if

    do i=1,Natoms
        Iq  = Iq + Eqsq*ff0(i) 
    end do

end subroutine compute_intensity

subroutine find_neighbor(x,y,z,radlist,proberad,neighborlist,Natoms)
    ! Algorithm to find neighboring atoms within the radius of an 
    ! atom + the probe diameter
    ! ####### USING SHRAKE RUPLEY ALGORITHM, see below #########
    ! Basic algorithm is to generate a mesh of point around
    ! an atom at its VDW radius. Then simply count the number 
    ! of points that don't overlap with other atoms (including probes)
    ! the SASA is then approximated by SASA = 4*pi*r^2 * Nexposed
    use DBL
    integer :: i,j,Natoms,found_neighbor
    real(dp) :: proberad, r, range
    real(dp) :: dist ! function
    real(dp), dimension(Natoms) :: x,y,z,radlist
    integer, dimension(Natoms,Natoms), intent(out) :: neighborlist

    neighborlist = 0
    ! find neighbor atoms: within a diameter of water/probe
    ! in this case is 2*proberad (usually proberad=1.4)
    ! pre-compute this and save index for each atom i
    do i=1,Natoms
        ! for every Atom i, the neighbors should be within range
        ! defined by the diameter of the probe
        range = radlist(i) + 2*proberad
        found_neighbor = 0
        do j=1,Natoms        
            if (i/=j) then
                r = dist(x(i),y(i),z(i),x(j),y(j),z(j))
                if (r < range + radlist(j)) then ! if Atom j is within range
                    found_neighbor = found_neighbor + 1
                    neighborlist(i,found_neighbor) = j
                end if
            end if
        end do
    end do

end subroutine find_neighbor

subroutine shrakerupley(x,y,z,radlist,proberad,Nfib,sasalist,Natoms)
    use DBL
    implicit none
    integer :: i,j,n,nid,Nfib,Natoms,naccess
    real(dp), dimension(Natoms) :: x,y,z,radlist
    real(dp), dimension(Natoms), intent(out) :: sasalist
    real(dp) :: proberad, r, testdistsq, fourpi, distsq ! probe radius
    integer, dimension(Natoms,Natoms) :: neighborlist
    logical :: exposed    
    real(dp), dimension(Nfib) :: x_fib, y_fib, z_fib
    real(dp) :: xs, ys, zs
    character*72 :: dumpfile

    dumpfile = "surface_points.xyz"
    fourpi = 4.0 * pi / real(Nfib)
    sasalist = 0.0
    ! compute a fibonacci-grid on a spherical surface
    call fibonacci(Nfib,x_fib,y_fib,z_fib)
    ! find neighboring atoms
    call find_neighbor(x,y,z,radlist,proberad,neighborlist,Natoms)

    ! open diagnostic file
    open(unit=500, file=dumpfile)
    print *,"File ", trim(dumpfile), " is open for writing..."
    ! compute distance between neighboring atoms
    do i=1,Natoms ! for every atom
        ! generate fibonacci grid with a given radius of atom i
        ! start fibonacci grid exploration
        naccess = 0

        do j=1,Nfib ! for every point
            xs = x_fib(j)*(radlist(i)+proberad) + x(i)
            ys = y_fib(j)*(radlist(i)+proberad) + y(i)
            zs = z_fib(j)*(radlist(i)+proberad) + z(i)
            exposed = .true.
            ! check with the neighbors
            do n=1,Natoms 
                if(neighborlist(i,n).gt.0) then
                    nid = neighborlist(i,n)
                    r   = radlist(nid) + proberad ! neighbor radius + probe
                    testdistsq = distsq(x(nid),y(nid),z(nid),xs,ys,zs)
                    ! if distance of atomic surface to neighboring atom
                    ! is within probe radius + neighboring atom radius
                    ! then this atom is buried
                    if (testdistsq .lt. r*r) then 
                        exposed = .false.
                        exit 
                    end if
                else
                    exit 
                end if
            end do 

            if (exposed) then
                naccess = naccess + 1
                write(500,105) xs,ys,zs ! write do dump
            105 format(3f8.3)
            end if   
        end do 
    end do 
        ! end fibonacci grid exploration

        ! the total amount of of points exposed is proportional to surface area
    sasalist(i) = fourpi * real(naccess) * radlist(i)*radlist(i)

    close(500)
        ! print out for diagnostics
    !     write (*,102) i, naccess, sasalist(i)
    ! 102 format("Atom ",i5," has ",i5, " mesh points exposed. SA=", f8.3)

end subroutine shrakerupley

function dist(xi,yi,zi,xf,yf,zf) result(r)
    use DBL
    ! simple function to compute distance
    real(dp), intent(in) :: xi,yi,zi,xf,yf,zf
    reaL(dp) :: r
    r = sqrt((xf-xi)**2 + (yf-yi)**2 + (zf-zi)**2)
end function dist

function distsq(xi,yi,zi,xf,yf,zf) result(r)
    use DBL
    ! simple function to compute distance^2
    real(dp), intent(in) :: xi,yi,zi,xf,yf,zf
    reaL(dp) :: r
    r = (xf-xi)**2 + (yf-yi)**2 + (zf-zi)**2
end function distsq

subroutine fibonacci(Npts,x_fib,y_fib,z_fib)
    ! given an origin xc,yc,zc and a radius of a sphere
    ! generate Npts points on the sphere following the golden-section spiral
    ! whose coordinates are x_fib,y_fib,z_fib 
    use DBL
    integer :: i,Npts
    real(dp) :: inc,offset,r,k,phi
    real(dp), dimension(Npts), intent(out) :: x_fib, y_fib, z_fib

    inc    = pi * (3.0-sqrt(5.0))
    offset = 2.0 / real(Npts)
    ! initialize points at center 
    x_fib = 0
    y_fib = 0
    z_fib = 0

    do i=1,Npts
        k = real(i)-1.0
        phi = k * inc
        y_fib(i) = (k*offset - 1.0) + 0.5*offset
        r   = sqrt(1.0 - y_fib(i)**2)
        x_fib(i) = cos(phi) * r
        z_fib(i) = sin(phi) * r
    end do

end subroutine fibonacci