subroutine init(data_root)

    use geomod10b, only: init1
   
    character(*) :: data_root

    print *, 'Calling geomod10 init with data root: ', trim(data_root)
    call init1(data_root)

end subroutine init


subroutine wind_emiss(freq,tht,n,sst,wind,phir,emiss_tot)
   
    use geomod10b, only: fd_wind_emiss

    real(4)  :: freq
    real(4)  :: tht
    integer(4) :: n
    real(4),dimension(n)  :: sst,wind,phir
    real(4),dimension(n,2) :: emiss_tot

    !f2py intent(in) :: freq,tht,sst,wind,phir
    !f2py intent(hidden), depend(sst) :: n = shape(sst)
    !f2py intent(out) :: emiss_tot

    integer(4) :: i
    real(4),dimension(4)   :: emiss_phi
    real(4),dimension(2)   :: emiss_totx
    
    do i=1,n
        call fd_wind_emiss(freq,tht,sst(i),wind(i),phir(i),emiss_totx,emiss_phi)
        emiss_tot(i,:) = emiss_totx
    enddo

end subroutine wind_emiss

    