!     Forked from Frank's .f version on the O-drive 10/12/2022
! converted to module

!     1/3/2018 changed 9/13/2022.  error message in fd_scatterm_all change from
!      if(opacty.lt.0)   stop 'wind oob in fd_scatterm, pgm stopped'  to
!      if(opacty.lt.0)   stop 'opacty oob in fd_scatterm, pgm stopped'


!     6/17/2015 changed 1/3/2018.  subroutine findstokes_toa added


!     this is the same as 'geomod9b.f' dated july 25 2009 except that:
!     1. xexp for eia dependence in wind model has changed for hpol
!     2. to compensate for change 1, the linear coef for 37H was slightly increased
!     3. 19H wind emiss decreased by 1%
!     see 'memo14.txt'

!     Many changes made during the april-july 2009 time frame.  Final version data is July 24 2009
module geomod10b

      use geomod10c, only: fdem0_meissner_wentz,dielectric_meissner_wentz,cosd,sind

      integer(4), parameter :: nstoke=4
      real(4)      :: acoef(5,2,6)
      real(4)      :: bcoef(5,2,nstoke,6)
      real(4)      :: rcoef(3,0:2,6)
      real(4)      :: scatterm(91,50,26,13,2)

      logical      :: initialized = .false.

contains

    subroutine init1(data_root)
      implicit none

      character(*) :: data_root

      print *, 'initializing geomod10b with data root: ', trim(data_root)
      print *, 'loading acoef from ' // trim(data_root) // 'finetune_emiss_wind.dat'
      open(3,file=trim(data_root) // '/finetune_emiss_wind.dat',&
                   access='stream',form='unformatted',action='read')
      read(3) acoef
      close(3)
      acoef(1,2,1)=1.0015*acoef(1,2,1)  !new geomod10
      acoef(1,2,2)=1.0170*acoef(1,2,2)  !new geomod10
      acoef(1,2,3)=0.9910*acoef(1,2,3)  !new geomod10
      acoef(1,2,5)=1.0271*acoef(1,2,5)  !new geomod10


      print *, 'loading bcoef from ' // trim(data_root) // 'fit_emiss_phir_wind.dat'
      open(3,file=trim(data_root) // '/fit_emiss_phir_wind.dat',&
                   access='stream',form='unformatted',action='read')
      read(3) bcoef
      close(3)
      
      print *, 'loading rcoef from ' // trim(data_root) // 'em0_ref_freq_sst.dat'
      open(3,file=trim(data_root) // '/em0_ref_freq_sst.dat',&
                   access='stream',form='unformatted',action='read')
      read(3) rcoef
      close(3)

      print *, 'loading scatterm from ' // trim(data_root) // 'mk_scatterm_table_all.dat'
      open(3,file=trim(data_root) // '/mk_scatterm_table_all.dat',&
                   access='stream',form='unformatted',action='read')
      read(3) scatterm
      close(3)

      initialized = .true.

      end subroutine init1

    subroutine findtb_toa(freq,tht,sst,wind,phir,tran,tbdw_in,tbup_in, tbv,tbh)
            implicit none
      
            real(4) freq,tht,sst,wind,phir,tran,tbdw_in,tbup_in, tbv,tbh
            real(4) emiss(2),emiss_phi(4),xscat(2)
            real(4) tbdown(2),surtep,tcos,tbdw,tbup
            real(4) costht,path,opacty,absp
      
            call fd_tcos_eff(freq, tcos)

      !     the absp**3 adjustment comes from emiss_model folder analysis.  it makes better agreement between tasim and tamea
      !     for very high water vapor at 85 ghz channels.  it appears that the teff expression derived from the raob is a bit
      !     too high at high absp, and this code corrects this.  
            absp=1-tran
            tbup=tbup_in -2.*absp**3
            tbdw=tbdw_in -2.*absp**3

            surtep=sst+273.16
            call fd_wind_emiss(freq,tht,sst,wind,phir,   emiss,emiss_phi)

      !    see folders o:\emiss_model\emiss_model and o:\tbmodel1\geoptics
            costht=cosd(tht)
            path=1.00035/sqrt(costht*costht+7.001225e-4)   !(1+hratio)/sqrt(costht**2+hratio*(2+hratio)), hratio=.00035
            opacty=-alog(tran)/path
            call fd_scatterm_all(freq,tht,wind,opacty, xscat) 

            tbdown(1)=tran*tcos + tbdw + xscat(1)
            tbdown(2)=tran*tcos + tbdw + xscat(2)

            tbv=tbup + tran*(emiss(1)*surtep + (1-emiss(1))*tbdown(1))
            tbh=tbup + tran*(emiss(2)*surtep + (1-emiss(2))*tbdown(2))
      
            end subroutine findtb_toa
      

    subroutine findstokes_toa(freq,tht,sst,wind,phir,tran,tbdw_in,tbup_in, tb_stokes)
            implicit none
      
            real(4) freq,tht,sst,wind,phir,tran,tbdw_in,tbup_in, tb_stokes(4)
            real(4) emiss(2),emiss_phi(4),xscat(2)
            real(4) tbdown(2),surtep,tcos,tbdw,tbup,avg_tbdown
            real(4) costht,path,opacty,absp
      
            call fd_tcos_eff(freq, tcos)

      !     the absp**3 adjustment comes from emiss_model folder analysis.  it makes better agreement between tasim and tamea
      !     for very high water vapor at 85 ghz channels.  it appears that the teff expression derived from the raob is a bit
      !     too high at high absp, and this code corrects this.  
            absp=1-tran
            tbup=tbup_in -2.*absp**3
            tbdw=tbdw_in -2.*absp**3

            surtep=sst+273.16
            call fd_wind_emiss(freq,tht,sst,wind,phir,   emiss,emiss_phi)

      !     see folders o:\emiss_model\emiss_model and o:\tbmodel1\geoptics
            costht=cosd(tht)
            path=1.00035/sqrt(costht*costht+7.001225e-4)   !(1+hratio)/sqrt(costht**2+hratio*(2+hratio)), hratio=.00035
            opacty=-alog(tran)/path
            call fd_scatterm_all(freq,tht,wind,opacty, xscat) 

            tbdown(1)=tran*tcos + tbdw + xscat(1)
            tbdown(2)=tran*tcos + tbdw + xscat(2)

            avg_tbdown=0.5*(tbdown(1) + tbdown(2))

            tb_stokes(1:2) = tbup + tran*(emiss(1:2)*surtep + (1-emiss(1:2))*tbdown(1:2))
            tb_stokes(  3) = tran*(surtep - avg_tbdown)*emiss_phi(3)
            tb_stokes(  4) = tran*(surtep - avg_tbdown)*emiss_phi(4)

    end subroutine findstokes_toa

!     ============================================================================================================================== 
!     ============================================================================================================================== 
!     ============================================================================================================================== 
 
    subroutine fd_wind_emiss(freq,tht,sst,wind,phir, emiss_tot,emiss_phi)
            implicit none

            !integer(4), parameter :: nstoke=4
      
            integer(4) ifreq1,ifreq2,ipol,istoke,iharm
            real(4) freq,tht,sst,wind,phir, emiss_tot(2),emiss_phi(4)
            real(4) em0(2)
            real(4) xexp(2),xexp_phir(2,nstoke),freq0(6),u
            real(4) thtref,qtht,wt,emiss1(2),emiss2(2),emiss(2),enad,h1,h2

            real(4) phirsv
            real(4) cos1phi,cos2phi,sin1phi,sin2phi
            real(4) aharm1(2,nstoke),aharm2(2,nstoke),aharm(2,nstoke),amp1,amp2,amp,anad(2,nstoke)

            data phirsv/1.e30/
            data thtref/55.2/
      !      data xexp/4.,1.5/ !see 'O:\emiss_model\eia_dependence.doc'
            data xexp_phir/2,2, 1,4, 1,4, 2,2/  !'O:\windsat\geomod\memo2.txt'
            data freq0/6.8000,  10.7000,  18.7000,  23.8000,  37.0000, 85.5/  !now reference to windsat and ssmi rather than amsr

            if(freq.lt.6.5) stop 'pgm stopped, freq too small in oob in      fd_wind_emiss'
            if (.not. initialized) stop 'pgm stopped, geomod10b not initialized in fd_wind_emiss'

            qtht=tht
            if(qtht.gt.65) qtht=65.  !qtht is just used for extrpolation for tht>thtref and i limit it to 60 deg      

            u=0
            if(freq.gt.7) then
                  u=(freq-7.)/12.
                  if(u.gt.1) u=1
            endif
            xexp(1)=4.
            xexp(2)=2.0 +1.*u*u*(3-2*u)  !new geomod10  
      

            call fdem0_meissner_wentz(freq,tht,sst, em0)
 
            ifreq1=1
            if(freq.gt.freq0(2)) ifreq1=2
            if(freq.gt.freq0(3)) ifreq1=3  !between 18.7 and 37 ghz
            if(freq.gt.freq0(5)) ifreq1=5
 
            if(ifreq1.ne.3) then
                  ifreq2=ifreq1+1
            else
                  ifreq2=ifreq1+2
            endif
 
            wt=(freq-freq0(ifreq1))/(freq0(ifreq2)-freq0(ifreq1))
            if(freq.gt.freq0(ifreq2)) wt=1  !only occurs for freq>85.5

!            =========================================================================================
!           =============================isotropic wind-induced emissivity ==========================
!           =========================================================================================

            call get_emiss_wind(ifreq1,sst,wind, emiss1)
            call get_emiss_wind(ifreq2,sst,wind, emiss2)
            emiss=(1-wt)*emiss1 + wt*emiss2 !emiss is thtref value interpolated to input freq

            enad=0.5*(emiss(1)+emiss(2))

            do ipol=1,2
                  if(tht.le.thtref) then
                        emiss(ipol)=enad        + (emiss(ipol)-enad)*(tht/thtref)**xexp(ipol)
                  else
                        emiss(ipol)=emiss(ipol) + (emiss(ipol)-enad)*(qtht-thtref)*xexp(ipol)/thtref  
                  endif
            enddo  !ipol

!     =========================================================================================
!     ================================= find emiss_phi ========================================
!     =========================================================================================

            if(phir.lt.-998. .or. wind.le.3) then !-999. default for doing no correction
                  emiss_phi=0
            else  !find emiss_phi
                  call get_aharm_phir(ifreq1,sst,wind, aharm1)  !aharm in terms of true stokes
                  call get_aharm_phir(ifreq2,sst,wind, aharm2)  !aharm in terms of true stokes
                  aharm=(1-wt)*aharm1 + wt*aharm2  !aharm is thtref value interpolated to input freq

!                 get nadir harmonic
                  call get_aharm_phir_nad(ifreq1,freq0(ifreq1),sst,wind, amp1) 
                  call get_aharm_phir_nad(ifreq2,freq0(ifreq2),sst,wind, amp2)  
                  amp=(1-wt)*amp1 + wt*amp2

                  anad=0  !most elements are zero
                  anad(2,2)=  amp
                  anad(2,3)= -amp

                  do istoke=1,nstoke
                        do iharm=1,2
                              if(tht.le.thtref) then
                                    aharm(iharm,istoke)= anad(iharm,istoke) + &
                                    (aharm(iharm,istoke)-anad(iharm,istoke))* &
                                    (tht/thtref)**xexp_phir(iharm,istoke)
                              else
                                    aharm(iharm,istoke)=aharm(iharm,istoke) + &
                                    (aharm(iharm,istoke)-anad(iharm,istoke))* &
                                    (qtht-thtref)*xexp_phir(iharm,istoke)/thtref  
                              endif
                        enddo  !iharm
                  enddo  !istoke

!                 convert back from true stokes to v and h
                  do iharm=1,2
                        h1=aharm(iharm,1) + 0.5*aharm(iharm,2) !(v+h)/2 + (v-h)/2=v
                        h2=aharm(iharm,1) - 0.5*aharm(iharm,2) !(v+h)/2 - (v-h)/2=h
                        aharm(iharm,1)=h1
                        aharm(iharm,2)=h2
                  enddo

                  if(abs(phir-phirsv).gt.0.01) then
                        phirsv=phir
                        cos1phi=cosd(  phir)
                        cos2phi=cosd(2*phir)
                        sin1phi=sind(  phir)
                        sin2phi=sind(2*phir)
                  endif

                  emiss_phi(1:2)=aharm(1,1:2)*cos1phi + aharm(2,1:2)*cos2phi
                  emiss_phi(3:4)=aharm(1,3:4)*sin1phi + aharm(2,3:4)*sin2phi
            endif

        emiss_tot=em0 + emiss  + emiss_phi(1:2)

      return
      end subroutine fd_wind_emiss



!     ============================================================================================================================== 
!     ============================================================================================================================== 
!     ============================================================================================================================== 

      subroutine get_emiss_wind(ifreq,sst,wind, emiss)
            implicit none

            integer(4) ifreq
            !integer(4)      :: istart
            real(4) sst,wind,emiss(2)
            !real(4) :: acoef(5,2,6)
            real(4) sst_fac(2)
            real(8) xmea(5)

            ! data istart/1/

            ! if(istart.eq.1) then
            ! istart=0
            ! !call openbig(3,'o:\emiss_model\finetune_emiss_wind.dat','old')
            ! open(3,file='/mnt/oserver/o/emiss_model/finetune_emiss_wind.dat',&
            !        access='stream',form='unformatted',action='read')
            ! read(3) acoef
            ! close(3)
            ! acoef(1,2,1)=1.0015*acoef(1,2,1)  !new geomod10
            ! acoef(1,2,2)=1.0170*acoef(1,2,2)  !new geomod10
            ! acoef(1,2,3)=0.9910*acoef(1,2,3)  !new geomod10
            ! acoef(1,2,5)=1.0271*acoef(1,2,5)  !new geomod10
            ! endif

            if(ifreq.eq.4) stop 'ifreq oob in get_emiss_wind, pgm stopped'

            call  fd_xmea_win(wind, xmea)
            emiss(1)=dot_product(acoef(:,1,ifreq),xmea) 
            emiss(2)=dot_product(acoef(:,2,ifreq),xmea) 

            call get_sst_fac(ifreq,1,sst, sst_fac(1))
            call get_sst_fac(ifreq,2,sst, sst_fac(2))

                  emiss=emiss*sst_fac
            
            return
      end subroutine get_emiss_wind


!   ============================================================================================================================== 
!   ============================================================================================================================== 
!   ============================================================================================================================== 

 
    subroutine get_aharm_phir(ifreq,sst,wind, aharm)
            implicit none

            !integer(4), parameter :: nstoke=4
      
            integer(4) ifreq,istoke,iharm
            ! integer(4)      :: istart
            real(4)      sst,wind,aharm(2,nstoke)
            real(4) h1,h2
            !real(4) :: bcoef(5,2,nstoke,6)
            real(4) sst_fac(nstoke)
            real(8) xmea(5)

            ! data istart/1/
 
            ! if(istart.eq.1) then
            !       istart=0
            !       !call openbig(3,'o:\emiss_model\fit_emiss_phir_wind.dat','old')
            !       open(3,file='/mnt/oserver/o/emiss_model/fit_emiss_phir_wind.dat', &
            !              access='stream',form='unformatted',action='read')
            
            !       read(3) bcoef
            !       close(3) 
            ! endif

            if(ifreq.eq.4) stop 'ifreq oob in get_aharm_phir, pgm stopped'
  
            call  fd_xmea_win(wind, xmea)

            call get_sst_fac(ifreq,1,sst, sst_fac(1))
            call get_sst_fac(ifreq,2,sst, sst_fac(2))
            sst_fac(3:4)=0.5*(sst_fac(1)+sst_fac(2))

            do istoke=1,nstoke
                  aharm(1,istoke)=sst_fac(istoke)*dot_product(xmea,bcoef(:,1,istoke,ifreq))
                  aharm(2,istoke)=sst_fac(istoke)*dot_product(xmea,bcoef(:,2,istoke,ifreq))
            enddo  !istoke

!           convert to true stokes paramters,ie (v+h)/2 and v-h, rather than v and h in order to do tht adjustment
            do iharm=1,2
                  h1=0.5*(aharm(iharm,1)+aharm(iharm,2))
                  h2=     aharm(iharm,1)-aharm(iharm,2)
                  aharm(iharm,1)=h1
                  aharm(iharm,2)=h2
            enddo
      end subroutine get_aharm_phir





! this routine is based on 'O:\windsat\geomod\memo2.txt' 
    subroutine get_aharm_phir_nad(ifreq,freq,sst,wind, amp)
            implicit none

            integer(4) ifreq
            real(4)      freq,sst,wind
            real(4) amp_10_nad,amp,ywind,qfreq
            real(4) sst_fac

            qfreq=freq
            if(qfreq.gt.37) qfreq=37

            if(freq.lt.3) then
                  amp_10_nad=.2/290.
            else
                  amp_10_nad=2*(1. - 0.9*alog10(30./qfreq))/290.
            endif
      
            ywind=wind
            if(wind.lt. 0) ywind= 0
            if(wind.gt.15) ywind=15

            amp=amp_10_nad*ywind*(ywind - ywind**2/22.5)/55.5556
            call get_sst_fac(ifreq,0,sst, sst_fac)
            amp=amp*sst_fac
            return
      end subroutine get_aharm_phir_nad




!   rcoef values for the ratio of nadir em0(sst)/em0(sst=20) come from 'O:\emiss_model\em0_ref_freq_sst.f'
      subroutine get_sst_fac(ifreq,ipol,sst, sst_fac)  !ipol=0 denotes nadir value, sst_fac=em0(sst)/em0(sst=20)
            implicit none

            integer(4) ifreq,ipol
            ! integer(4)      :: istart
            real(4) sst,sst_fac
            real(4) xmea(3)
            ! real(4) :: rcoef(3,0:2,6)
            ! data istart/1/

            ! if(istart.eq.1) then
            !       istart=0
            !       !call openbig(3,'o:\emiss_model\em0_ref_freq_sst.dat','old')
            !       open(3,file='/mnt/oserver/o/emiss_model/em0_ref_freq_sst.dat',&
            !              access='stream',form='unformatted',action='read')
            !       read(3) rcoef
            !       close(3)
            ! endif

            xmea(1)= sst-20
            xmea(2)=xmea(1)*xmea(1)
            xmea(3)=xmea(1)*xmea(2)
            sst_fac=1 + dot_product(rcoef(:,ipol,ifreq),xmea)  

            return
      end subroutine get_sst_fac

!     ============================================================================================================================== 
!     ============================================================================================================================== 
!     ============================================================================================================================== 

      subroutine fd_xmea_win(wind, xmea)
            implicit none
      
            real(4) wind,x,dif,wcut
            real(8) xmea(5)
            data wcut/20./
            
            x=wind
            if(x.lt.0) x=0
      
            xmea(1)=x
            if(x.le.wcut) then
            xmea(2)=xmea(1)*x
            xmea(3)=xmea(2)*x
            xmea(4)=xmea(3)*x
            xmea(5)=xmea(4)*x

            else
            dif=x-wcut
            xmea(2)=2*dif*wcut       + wcut**2
            xmea(3)=3*dif*wcut**2    + wcut**3
            xmea(4)=4*dif*wcut**3    + wcut**4
            xmea(5)=5*dif*wcut**4    + wcut**5
            endif
      
            return
      end subroutine fd_xmea_win

!     ============================================================================================================================== 
!     ============================================================================================================================== 
!     ============================================================================================================================== 


    subroutine fd_scatterm_all(freq,tht,wind,opacty, xscat) 
            implicit none

            real(4) freq,tht,wind,opacty, xscat(2)
            real(4) xlog_freq,xscat1(2),xscat2(2)
            real(4) a1,a2,b1,b2,c1,c2,brief,d1,d2
            !(4)      ::  scatterm(91,50,26,13,2)

            !integer(4)      :: istart
            integer(4) i1,i2,j1,j2,k1,k2,l1,l2

            ! data istart/1/

            ! if(istart.eq.1) then
            !       istart=0
            !       !call openbig(3,'o:\tbmodel1\geoptics\mk_scatterm_table_all.dat','old')
            !       open(3,file='/mnt/oserver/o/tbmodel1/geoptics/mk_scatterm_table_all.dat',&
            !              access='stream',form='unformatted',action='read')
            !       read(3) scatterm
            !       close(3)
            ! endif

!           check inputs

            if(freq.lt.1 .or. freq.gt.200) stop 'freq oob in fd_scatterm, pgm stopped'
            if(tht .lt.0 .or.  tht.gt. 90) stop 'tht  oob in fd_scatterm, pgm stopped'
            if(wind.lt.0 .or. wind.gt.100) stop 'wind oob in fd_scatterm, pgm stopped'
      !bug  if(opacty.lt.0)                stop 'wind oob in fd_scatterm, pgm stopped'
            if(opacty.lt.0)                stop 'opacty oob in fd_scatterm, pgm stopped'

            xlog_freq=alog10(freq)
 
!           do time,lat,lon interpolation
 
            brief=tht
            if(brief.gt.89.99) brief=89.99
            i1=1+brief
            i2=i1+1
            a1=i1-brief
            a2=1.-a1
 
            brief=wind
            if(brief.gt.24.99) brief=24.99
            j1=1+brief
            j2=j1+1
            b1=j1-brief
            b2=1-b1

            brief=xlog_freq/0.2
            if(brief.gt.11.99) brief=11.99
            k1=1+brief
            k2=k1+1
            c1=k1-brief
            c2=1-c1

            brief=opacty/0.025
            if(brief.gt.48.99) brief=48.99
            l1=1+brief
            l2=l1+1
            d1=l1-brief
            d2=1-d1

        xscat1= &
                  a1*b1*(c1*scatterm(i1,l1,j1,k1,:)+c2*scatterm(i1,l1,j1,k2,:))+ &
                  a1*b2*(c1*scatterm(i1,l1,j2,k1,:)+c2*scatterm(i1,l1,j2,k2,:))+ &
                  a2*b1*(c1*scatterm(i2,l1,j1,k1,:)+c2*scatterm(i2,l1,j1,k2,:))+ &
                  a2*b2*(c1*scatterm(i2,l1,j2,k1,:)+c2*scatterm(i2,l1,j2,k2,:))

            xscat2= &
                  a1*b1*(c1*scatterm(i1,l2,j1,k1,:)+c2*scatterm(i1,l2,j1,k2,:))+ &
                  a1*b2*(c1*scatterm(i1,l2,j2,k1,:)+c2*scatterm(i1,l2,j2,k2,:))+ &
                  a2*b1*(c1*scatterm(i2,l2,j1,k1,:)+c2*scatterm(i2,l2,j1,k2,:))+ &
                  a2*b2*(c1*scatterm(i2,l2,j2,k1,:)+c2*scatterm(i2,l2,j2,k2,:))

            xscat=d1*xscat1 + d2*xscat2
 
      return
      end subroutine fd_scatterm_all

!     for these routine the term b is the flux (2*h*f**3/(c**2*(dexp(h*f/(k*t))-1))) divided by  2*k*f**2/c**2
!     see 'o:\emiss_model\planck.docx'
    subroutine fd_tcos_eff(freq, tcos_eff)
            implicit none
            real(8), parameter :: tcos=2.73
            real(8), parameter :: teff=63.  !selected to provide optimum fit over 60-300 k range
            real(8), parameter :: h=6.6260755d-34
            real(8), parameter :: k= 1.380658d-23
            real(8), parameter :: a=h/k
            real(8)  x,b1,b2
            real(4) freq,tcos_eff

            x=a*freq*1.d9
            b1=x/(dexp(x/tcos)-1)
            b2=x/(dexp(x/teff)-1)
            tcos_eff=b1-b2+teff
            return
      end subroutine fd_tcos_eff

end    module geomod10b
