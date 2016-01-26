
module DBL
    implicit none
    integer, parameter :: dp=kind(0.0d0)
end module DBL

subroutine compute_intensity(q,x,y,z,ff0,Natoms,Iq,Nq,Nr)
    use DBL
    implicit none
    integer :: Natoms
    integer :: i,j,r_i,Nq
    real(dp), dimension(Natoms) :: x,y,z,ff0
    real(dp), dimension(Nq) :: q, Eqsq
    real(dp), dimension(Nq), intent(out) :: Iq
    real(dp) :: dist, b, dr, dq, fi, fj, rem_dr, Dmax, qr
    real(dp), dimension(1000) :: r,Pr
    integer, intent(out) :: Nr

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
    print *,"P(r) computed, computing auto-correlation part."
    write (*,101) Dmax
101 format("Dmax = ",f6.2)

    ! perform fourier transform and add auto-correlation part
    Nr = int( (Dmax-mod(Dmax,dr))/dr ) + 1
    do i=1,Nr
        r(i) = (i-1)*dr
    end do

    do i=1,Nq
        do j=1,Nr
            qr = q(i)*r(j)
            if (qr.lt.1e-8) then
                Iq(i)  = Iq(i) + Pr(j) * dr
            else
                Iq(i)  = Iq(i) + Pr(j) * dr * sin(qr)/qr
            end if
        end do
    end do
    ! compute Eq^2 (N + 2 * P(r)sinc(qr))
    ! and subtract in vacuo scattering with dummy excluded solvent
    Iq  = Eqsq*(Natoms + 2*Iq)

end subroutine compute_intensity