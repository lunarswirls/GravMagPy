      program gravmag_sphere_brtp
c*************************************************************************
c
c Sphere is for computing z, x and y  components on sphere (w/o mainfield)
c Tiku Ravat 18 August 2000  modified for y component sign on 160525
c
c ------------------------------------------------------------------------
c Changes by Dany Waller:
c   - Student of Tiku Ravat (2017-2019); modified this program under his
c         supervision in Fall 2018 EES-395 Gravity and Magnetics course.
c   - Generalized all sphere/sphere/sphere naming to a generic 'sphere'.
c   - Removed hard-coded body radius, radius is now passed as an argument
c         in kilometers (km).
c   - Option to compute all magnetic field components at once and output
c         to combined file with header.
c ------------------------------------------------------------------------
c
c  The program does both gravity and magnetic calculations.
c
c  Compile the program at your site and then use gravmag_sphere.bash as a guide
c
c..For facilitating testing of this program at various sites three test inputs
c  and outputs have been provided. number of gaussian nodes (nr, ntheta,
c  nphi) will always determine the accuracy of calculations, but the more the
c  nodes more the time spent in calculations. these can only be 2, 4, 8, 16.
c  there is no specific rule of thumb to select nr, ntheta, nphi - these
c  depend on thickness, and n-s and e-w widths of the body and obs. elevation,
c  etc. try: 2, 8, 8 then 2, 16, 16. If results change drastically between
c  8 and 16 then you certainly need 16. If that is not enough, you may need
c  to split the body in smaller pieces. for most thicknesses of crustal
c  dimensions (sphere's continental) nr = 2 is sufficient.
c
c  (No tests yet for sphere are being provided but they have been made.)
c  the three input test files are: gtestfigV.2.in, test.gravin and test.magin
c  the three output testfiles are gtestfigV.2.outgrid test.gravoutgrid, etc.
c
c  REFERENCES: (in publications, please give appropriate credit to following
c               individuals for their
c               efforts in developing these programs)
c
c  von Frese, R.R.B., W.J. Hinze, L.W. Braile, and A.J. Luca, 1981,
c       Spherical-sphere Gravity and Magnetic Anomaly Modeling by Gauss-
c       Legendre Quadrature Integration, Journal of Geophysics.,
c       vol.49,pp.234-242.
c
c  Ravat, D., 1989, Magsat Investigations over the Greater African Region,
c       Ph.D. Dissertation, Purdue University, West Lafayette, IN, 234p.
c
c***********************************************************************
c  tape1 = gauss.data
c          (gauss.data = gauss-legendre quadrature coefficients)
c  tape2 = program output
c  tape4 = program parameters
c  tape5 = input
c  tape6 = output
c  tape10 = unformatted read/write (scratch)
c***********************************************************************
c..!!!This is a majorly modified version of program sphere!!!
c..   modified by tiku ravat  20 february 1989 (after about 30 continuou
c..   s hours of attempts to debug sphere.
c..   this program does not do any cubic spline interpolation. the
c..   program estimates the anomaly of a polygonal body with flat
c..   (meaning, constant elevation) upper
c..   and lower surfaces. the data input is slightly different than
c..   the original program (it's easier!).
c..the major advantage of using this program over the original sphere
c..is that this program works!!
c***********************************************************************
c     Program Sphere calculates the gravity or magnetic anomaly over a
c     spherical grid due to an arbitrary 3-dimensional body within the
c     sphere.  The shape of the body in the three spherical dimensions is
c     approximated by the 'mean slope' cubic spline interpolation method
c     to yield limits for numerical integration in spherical coordinates
c     by gauss-legendre quadrature.  for a uniform body (e.g., prism,
c     etc.) the user may by-pass the interpolation procedure and supply
c     the limits of integration directly.  output options are given in
c     data card #6.
c
c     To obtain the total anomaly due to multiple bodies repeat the com-
c     plete data card input sequence outlined below for each body.  when
c     this option is used data card #2 must be the same for each body
c     considered.
c
c******************* list of subroutines used  in sphereii *********************
c
c(modified to read2)** read1             ** order2 (eliminated)
c( " " to srch1)    ** search            ** rlim   (   ""     )
c( " " to trans1)   ** trans             ** change (modified to xchang)
c( " " to order1)   ** order             ** sum1   **
c (eliminated)      ** bound1            ** gravs  **
c (   ""     )      ** cubic             ** rho1   **
c (   ""     )      ** bound2            ** contr  (eliminated 3/25/2003)
c (   ""     )      ** bpick             ** zero   **
c (   ""     )      ** spline            ** mags   **
c (   ""     )      ** setup             ** sus1   **
c (   ""     )      ** store             ** geomag **
c (   ""     )      ** interp            ** fidd   (nasa) **
c (   ""     )      ** fix               ** fid    (nasa) **
c                                        ** magf   (nasa) **
c (   ""     )      ** put               ** thetlm (new subroutine)
c
c***********************************************************************
c
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      real lat,long
      character*80 title
      common /tot/ xx(544),totalf(1801,1801)
      common /body/ nb, iprint, thetas, thetan, bound(10000,2)
      common /nodes/ phi(32),theta(64),r(32)
      common /coeffs/ cr1(32),cr2(32),ct1(32),ct2(32),cp1(32),cp2(32)
      common tsum(1801,1801),rsum(1801,1801),gfield(1801,1801),grav(1801
     &,1801)
      common /fxlim/ nblim,thmax,thmin,phmax,phmin,rrmax,rrmin
      common /rem/ nrem,rf1,ri1,rd1
c ------------------------------------------------------------------------
c Added by Dany Waller:
c   When magnetic field is requested (ifield > 0 in input), the program
c   now computes and outputs *all* magnetic field components and total
c   field in one run: Br, Btheta, Bphi, Btot.
c   Combined (multi-body) output columns (lon lat Br Btheta Bphi Btot)
c   are written to unit 2 (fort.2) when nfile != 0.
c ------------------------------------------------------------------------
      common /outfields/ brtot(1801,1801),bthetatot(1801,1801),
     &                  bphitot(1801,1801),btottot(1801,1801)
      integer*4 multi_out, ifield_in, icomp
c
c ------------------------------------------------------------------------
c Sphere reference radius (km) is provided at run time as argument #1.
c Usage: gravmag_sphere_brtp <R_sphere_km>
c ------------------------------------------------------------------------
      character*64 arg1
      integer*4 ios
      call get_command_argument(1,arg1)
      if (arg1 .eq. ' ') then
         write(*,*) 'ERROR: missing required argument R_sphere_km (km)'
         write(*,*) 'Usage: gravmag_sphere_gauss2_32_brtp <R_sphere_km>'
         stop
      endif
      read(arg1,*,iostat=ios) rsphere
      if (ios .ne. 0) then
         write(*,*) 'ERROR: could not parse R_sphere_km from: ', arg1
         stop
      endif

      knk=0
      nkn=0
c  tape1 = gauss.data
c  tape2 = program output
c  tape4 = program parameters
c  tape5 = input
c  tape6 = output
c  tape10 = unformatted read/write recfm u lrecl 80
c     open(1,file='gauss.data',status='old')
c     open(2,file='sm.eps.l0.r0.e395.testg77',status='unknown')
c     open(4,file='sm.eps.l0.r0.in',status='old')
      open(1,status='old')
      open(2,status='unknown')
      open(4,status='old')
      open(10,file='sphereiiscratch',
     . form='unformatted',status='unknown')
      open(6,file='sphout',status='unknown')
c
c        read in program parameters
c
 110  read (4,830,end=640) title
  120 knk=knk+1
      read (4,840) lat,long,dlat,dlong,elvo,nlat,nlong
      read (4,860) nr,ntheta,nphi,nblim
      read (4,850) htheta,hphi,phi1,phi2,rho
      read (4,710) ifield,nrem,rf1,ri1,rd1
      read (4,1130) iprint,nfile,nopt,first,conint
c ------------------------------------------------------------------------
c Added by Dany Waller: store requested field selector and decide whether
c we will compute all magnetic components in a single run.
c   ifield_in = value read from input file (card 5).
c   multi_out = 1 => compute Br, Btheta, Bphi, Btot (magnetic case).
c   multi_out = 0 => original behavior (gravity or single field).
c ------------------------------------------------------------------------
      ifield_in = ifield
      multi_out = 0
      if (ifield_in .gt. 0) multi_out = 1
c
      if (ifield.eq.0) write (6,1010)
      if (ifield.eq.2) write (6,1030)
      if (ifield.eq.3) write (6,1040)
      if (ifield.eq.4) write (6,1050)
      write (6,830) title
      write (6,980) lat,long,elvo,dlat,dlong,nlat,nlong
      write (6,990) htheta,hphi,phi1,phi2,rho
      write (6,1000) nr,ntheta,nphi
      write (6,1140) nfile,nopt
      if (nrem.eq.1) write (6,720) nrem,rf1,ri1,rd1
c
c        convert to spherical coordinates
c
      ylati = lat
      xlongi = long
      dylat  = dlat
      dxlong = dlong
      pi=3.1415926
c..for gravity of other spherical-planetary objects, replace 6371.0 by the
c  appropriate number.
c..sphere mean equatorial radius
c     rsphere=<set via command line>
      fact=pi/180.
      ro=elvo+rsphere
      lat=(90.-lat)*fact
      long=long*fact
      dlat=-dlat*fact
      dlong=dlong*fact
      sus=rho

      ri1=ri1*fact
      rd1=-rd1*fact
c
c
c        read gauss-legendre quadrature coefficients
c
  140 write (6,900)
      read (1,*) (xx(j),j=1,544)
      rewind 1
c
      n2=nr/2
      nt=n2*(n2-1)
      do 150 i=1,nr
         cr2(i)=xx(nt+i)
c this is for 16 nodes  72 coefficients, total in xx (144): x and coeff
c        cr1(i)=xx(72+nt+i)
c this is for 32 nodes  272 coefficients, total in xx(544): x and its coeff
         cr1(i)=xx(272+nt+i)
  150 continue
      write (6,930)
      do 160 i=1,nr
         write (6,910) cr1(i),cr2(i)
  160 continue
c
      n2=ntheta/2
      nt=n2*(n2-1)
      do 170 i=1,ntheta
         ct2(i)=xx(nt+i)
c this is for 16 nodes  72 coefficients, total in xx (144): x and coeff
c        ct1(i)=xx(72+nt+i)
c this is for 32 nodes  272 coefficients, total in xx(544): x and its coeff
         ct1(i)=xx(272+nt+i)
  170 continue
      write (6,940)
      do 180 i=1,ntheta
         write (6,910) ct1(i),ct2(i)
  180 continue
c
      n2=nphi/2
      nt=n2*(n2-1)
      do 190 i=1,nphi
         cp2(i)=xx(nt+i)
c this is for 16 nodes  72 coefficients, total in xx (144): x and coeff
c        cp1(i)=xx(72+nt+i)
c this is for 32 nodes  272 coefficients, total in xx(544): x and its coeff
         cp1(i)=xx(272+nt+i)
  190 continue
      write (6,920)
      do 200 i=1,nphi
         write (6,910) cp1(i),cp2(i)
  200 continue
c
c        read in body coordinates
c
      if (nblim.eq.0) go to 210
      read (4,850) thmax,thmin,phmax,phmin,rrmax,rrmin
      write (6,730)
      write (6,740) thmax,thmin,phmax,phmin,rrmax,rrmin
      go to 270
  210 call read2 (tmax,pmax,pmin,rrmax,rrmin)
c
c        group body coordinates by increasing 'phi' variable
c
      zphi1=pmin
      j1=1
      j0=1
      nbni1=nb
      do 220 i=1,nbni1
         if (iprint.eq.2) write (6,750) i,zphi1,j1
         call srch1 (zphi1,j1)
         if (zphi1.gt.phi2) go to 230
         zphi1=zphi1+hphi
  220 continue
c
c        arrange 'phi' groups by increasing 'theta'
c
  230 n=i
      j1=1
      j0=1
      zphi1=pmin
      do 240 i=1,n
         if (iprint.eq.2) write (6,770) i,zphi1,j1
         call order1 (zphi1,j1)
         zphi1=zphi1+hphi
  240 continue
c
c        write body points
c
      write (6,1090) nb
      do 250 i=1,nb
         write (6,1100) i,(bound(i,j),j=1,2)
  250 continue
c
c
c        find gaussian nodes for the 'phi' variable
c
c
  270 if (nblim.eq.1) phi1=phmin
      if (nblim.eq.1) phi2=phmax
  280 call xchang (phi1,phi2,3,nphi,lz1,lz2)
      if (nblim.eq.1) go to 320
      if (phi(1).gt.pmin) go to 290
      phi1=phi1+hphi*0.1
      go to 280
  290 if (phi(nphi).lt.pmax) go to 300
      phi2=phi2-hphi*0.1
      go to 280
c
c
  300 write (6,1080) phi1,phi2
      if (iprint.ne.2) go to 310
      write (6,790)
      write (6,800) (phi(i),i=1,nphi)
  310 pphi=pmin
      if (iprint.eq.2) write (6,810) pphi
      ib=1
      ii=1
      istart=1
  320 continue
c ------------------------------------------------------------------------
c Added by Dany Waller: when magnetic field is requested (multi_out=1),
c compute Br, Btheta, Bphi, and Btot in one run. Implemented by looping
c over an internal component index (icomp) and reusing the original
c Gaussian integration workflow for each component.
c ------------------------------------------------------------------------
      if (multi_out .eq. 1) then
c initialize combined output arrays once
         if (knk .eq. 1) call zero_outfields(nlat,nlong)
         icomp = 1
 3410    continue
c choose component to compute
c   icomp = 1 -> Br     (ifield=2)
c   icomp = 2 -> Btheta (ifield=4)
c   icomp = 3 -> Bphi   (ifield=3)
c   icomp = 4 -> Btot   (ifield=1)
         if (icomp .eq. 1) ifield = 2
         if (icomp .eq. 2) ifield = 4
         if (icomp .eq. 3) ifield = 3
         if (icomp .eq. 4) ifield = 1
c zero working arrays for this body+component
         call zero (nlat,nlong,-1)
         do 331 n1=1,nlat
         do 331 n2=1,nlong
  331    totalf(n1,n2)=0.0
         go to 340
      endif
c ------------------------------------------------------------------------
c Original method (single magnetic field component or gravity):
c ------------------------------------------------------------------------
      call zero (nlat,nlong,-1)
      if (knk.gt.1) go to 340
      do 330 n1=1,nlat
      do 330 n2=1,nlong
  330 totalf(n1,n2)=0.0
c
c        do loop functions
c
c            1. determine limits of integration for the 'theta'
c               variable at each 'phi' node
c            2. find gaussian nodes for the 'theta' variable
c            3. find limits of integration for the 'r' variable
c               at each 'theta' node
c
  340 do 530 i=1,nphi
         phin=phi(i)
         phix=phin*fact
         if (nblim.eq.1) go to 400

c..finds theta limits at the phin
         call thetlm(phin)

  400    nlim=2
c..calculates theta nodes
         call xchang (phi1,phi2,2,ntheta,nlim,0)

c..radial limits are already fixed to rrmin and rrmax and tranferred to
c  xchang via common fxlim
         write (6,970)
         j2=1
         do 520 j1=1,nlim,2
            call zero (nlat,nlong,1)
            do 490 j=1,ntheta
               thetaj=theta(j)
               thetax=(90.-thetaj)*fact
               call zero (nlat,nlong,0)
c
c        find gaussian nodes for the 'r' variable
c
               call xchang (phi1,phi2,1,nr,lz1,j2)
c
c        perform the gaussian summation
c
               do 460 l=1,nr
                  rl=r(l)
                  write (6,890) rl,thetaj,phin
                  rx=rl+rsphere
                  if (ifield.eq.0) go to 440
               call magssphere(ro,lat,long,rl,thetaj,phin,dlat,dlong,nlat,
     1            nlong,sus,ifield,rsphere)
                  go to 450
  440             call gravs (ro,lat,long,rx,thetax,phix,dlat,dlong,nlat
     1            ,nlong,rho)
  450             call sum1 (nlat,nlong,i,j,l)
  460          continue
cc             if (nblim.eq.1) dr=0.5*(rrmax-rrmin)
cc             if (nblim.eq.0) dr=0.5*(r2(j2)-r1(j2))
c..next statement is added instead of the previous two
               dr=0.5*(rrmax-rrmin)

               do 480 l1=1,nlat
                  do 470 l2=1,nlong
                     tsum(l1,l2)=tsum(l1,l2)+dr*rsum(l1,l2)
  470             continue
  480          continue
c..do not update j2, because the radial limits are fixed
cc             j2=j2+1
  490       continue
            if (nblim.eq.1) dtheta=0.5*(thmax-thmin)*fact
c           if (nblim.eq.0) dtheta=0.5*(xlim(j1+1,1)-xlim(j1,1))*fact
c..next statement is added instead of the previous one
            if (nblim.eq.0) dtheta=0.5*(thetan-thetas)*fact

            do 510 n1=1,nlat
               do 500 n2=1,nlong
                  gfield(n1,n2)=gfield(n1,n2)+dtheta*tsum(n1,n2)
  500          continue
  510       continue
  520    continue
c
c
  530 continue
c..these are o.k.
      if (nblim.eq.1) dphi=0.5*(phmax-phmin)*fact
      if (nblim.eq.0) dphi=0.5*(phi2-phi1)*fact
      do 550 m1=1,nlat
         do 540 m2=1,nlong
c remember that all these gfield or totalf are really based on the value specified for ifield
c ifield = 0 grav, ifield = 1 total field, ifield = 2 = z-comp,  ifield = 3 = now y-component (see below change of sign)
c ifield = 4 = x-comp
            gfield(m1,m2)=dphi*gfield(m1,m2)
            totalf(m1,m2)=totalf(m1,m2)+gfield(m1,m2)
  540    continue
  550 continue
c ------------------------------------------------------------------------
c Added by Dany Waller: in magnetic multi-component mode, store the
c combined (multi-body) field for this component in the appropriate
c accumulator array and then compute the next component.
c ------------------------------------------------------------------------
      if (multi_out .eq. 1) go to 5560
      write(6,444) nfile
 444  format(1x,8hnfile = ,i5)
      if (nfile.ne.1) go to 570
c
c        write single body anomaly on #tape2#
c
      do 560 i=1,nlat
         write (2,1150) (gfield(i,j),j=1,nlong)
  560 continue
      write (2,680)
  570 continue
c
c     zero arrays for subsequent body point interpolations (i.e., for
c     knk > 1)
c
      do 602 i=1,10000
         do 603 j=1,2
            bound(i,j)=0.0
  603    continue
  602 continue
c
      go to 110
 5560 continue
c ------------------------------------------------------------------------
c Added by Dany Waller: add current component result (totalf) to the
c global accumulators. Note sign convention: the original code wrote
c -totalf for component outputs (ifield=2,3,4). I preserve that here by
c storing the physical components with the same sign that would have been
c written to fort.2 in the single-component mode.
c ------------------------------------------------------------------------
      do 5565 m1=1,nlat
         do 5564 m2=1,nlong
            if (icomp .eq. 1) brtot(m1,m2) = brtot(m1,m2) - totalf(m1,m2)
            if (icomp .eq. 2) bthetatot(m1,m2) = bthetatot(m1,m2) - totalf(m1,m2)
            if (icomp .eq. 3) bphitot(m1,m2) = bphitot(m1,m2) - totalf(m1,m2)
            if (icomp .eq. 4) btottot(m1,m2) = btottot(m1,m2) + totalf(m1,m2)
 5564    continue
 5565 continue
c next component (or finish this body)
      icomp = icomp + 1
      if (icomp .le. 4) go to 3410
c restore requested selector for logging consistency
      ifield = ifield_in
      do 600 i=1,10000
       do 600 j = 1, 2
  600 bound(i,j)=0.0
      go to 110
  640 if (knk.lt.1) go to 670
      nkn=nkn+1
c
c        write multiple body anomaly on #tape2#
c

c
      if (nfile.eq.0) go to 670
c     do 641 j = 1, nlat
c        do 642 i = 1, nlong
c           totalf(i,j) = totalf(i,j) - amean
c642     continue
c641  continue
c     do 650 i=1,nlat
c        write (2,1150) (totalf(i,j),j=1,nlong)
c 650 continue

c ------------------------------------------------------------------------
c Added by Dany Waller: final write-out. In magnetic multi-component mode
c (multi_out=1), write a 6-column XYZ table: lon lat Br Btheta Bphi Btot.
c Otherwise, preserve original behavior (single field in totalf).
c ------------------------------------------------------------------------
      if (multi_out .eq. 1) then
c lightweight header for downstream plotting scripts
         write(2,*) '# lon_deg lat_deg Br_nT Btheta_nT Bphi_nT Btot_nT'
         do 7655 i = 1, nlat
            ylat = ylati + real(i-1)*dylat
            do 7660 j = 1, nlong
               xlong = xlongi + real(j-1)*dxlong
               write(2,*) xlong, ylat, brtot(i,j), bthetatot(i,j),
     &                     bphitot(i,j), btottot(i,j)
 7660       continue
 7655    continue
      else
      do 655 i = 1, nlat
      ylat = ylati + real(i-1)*dylat
      do 660 j = 1, nlong
      xlong = xlongi + real(j-1)*dxlong
c..the original phi component was westward (ifield.eq.3) so change sign to make it eastward or Bphi
c now writing Br, Btheta, and Bphi components instead of x and z components
      if (ifield.eq.3.or.ifield.eq.2.or.ifield.eq.4) then
      write(2,*)  xlong, ylat, -totalf(i,j)
      else
      write(2,*)  xlong, ylat, totalf(i,j)
      endif

 660     continue
 655  continue
      endif


  670 continue
c
  680 format (///)
  700 format (8a10)
  710 format (2i5,5f10.0)
  720 format (/,1x,5hnrem=,i3,2x,4hrf1=,f10.3,2x,4hri1=,f10.5,2x,4hrd1=,
     1f10.5,2x,3hri=,f10.5,2x,3hrd=,f10.5)
  730 format (/,1x,26hbody limits of integration)
  740 format (/,1x,6hthmax=,f8.4,2x,6hthmin=,f8.4,2x,6hphmax=,f8.4,2x,6h
     1phmin=,f8.4,2x,6hrrmax=,f8.4,2x,6hrrmin=,f8.4)
  750 format (/,1x,2hi=,i5,2x,11hcall srch1(,f8.4,i5,1h))
  770 format (/,1x,2hi=,i5,2x,12hcall order1(,f8.4,i5,1h))
  790 format (//,1x,9hphi-nodes,/)
  800 format (8f10.4)
  810 format (//,1x,5hpphi=,f15.11)
  820 format (/,1x,3hib=,i5,5x,3hii=,i5,5x,5hpphi=,f10.6)
  830 format (a80)
  840 format (5f10.0,2i5)
  850 format (6f10.0)
  860 format (4i5)
  870 format (8f10.8)
  880 format (///,10x,53hsteps in the calculation of the limits of integ
     1ration,///,10x,6hnode =,i5,15x,19hthe value of phi is,2x,f10.4,5x,
     25hnlim=,i5,//,13x,5htheta,32x,2hr2,18x,2hr1,//)
  890 format (4(10x,f10.5))
  900 format (///,10x,26hgaussian quadrature coeff.,10x,26hcoordinates o
     1f subdivision,/)
  910 format (18x,f10.8,26x,f10.8)
  920 format (//,39x,3hphi,//)
  930 format (//,39x,1hr,//)
  940 format (//,39x,5htheta,//)
  950 format (///,10x,35hpoints used in cross interpolations,///,13x,5ht
     1heta,15x,3hphi,17x,2hr2,17x,2hr1,//)
  960 format (///,10x,26hfixed interpolation points,///,13x,5htheta,11x,
     12hr2,17x,2hr1,//)
  970 format (///,10x,35hgaussian nodes - equivalent sources,///,11x,9he
     1levation,10x,8hlatitude,11x,9hlongitude,//)
  980 format (10x,31hobservation grid specifications,///,20x,6horigin,//
     1,15x,8hlatitude,3x,f8.4,/,14x,9hlongitude,3x,f8.4,/,14x,9helevatio
     2n,3x,f8.4,//,20x,12hgrid spacing,//,19x,4hdlat,3x,f8.4,/,18x,5hdlo
     3ng,3x,f8.4,//,20x,15hgrid dimensions,//,19x,4hnlat,3x,i4,/,18x,5hn
     4long,3x,i4,///)
  990 format (////,8x,6hhtheta,3x,f8.4,/,10x,4hhphi,3x,f8.4,//,25x,28hli
     1mits of the causative body,//,10x,4hphi1,3x,f8.4,/,10x,4hphi2,3x,f
     28.4,//10x,16hdensity contrast,3x,f8.6,///)
 1000 format (10x,24hnumber of gaussian nodes,//,15x,1hr,3x,i4,/,11x,5ht
     1heta,3x,i4,/,13x,3hphi,3x,i4,//)
 1010 format (1h1,25x,7hgravity,////)
 1020 format (1h1,25x,20htotal magnetic field,////)
 1030 format (1h1,25x,29hmagnetic field - z- component,////)
 1040 format (1h1,25x,39hmagnetic field - east or  phi component,
     1////)
 1050 format (1h1,25x,29hmagnetic field - x- component,////)
 1060 format (1h1,25x,43hmagnetic field - total horizontal component,///
     1/)
 1070 format (1h1,5x,70hcoefficients for spherical harmonic expansion of
     1 the geomagnetic field,///)
 1080 format (///,25x,37hadjusted limits of the causative body,///,15x,4
     1hphi1,3x,f8.4,//,15x,4hphi2,3x,f8.4,//)
 1090 format (///,10x,23hlist of boundary points,///,15x,5htheta,11x,3hp
     1hi,12x,2hr2,13x,2hr1,8x,12hboundary no.,9x,9htotal no.,2x,i3,///)
 1100 format (2x,i3,2x,2(5x,f10.4))
 1120 format (//,10x,36hlimits of integration - *r* variable,//,13x,5hth
     1eta,32x,2hr2,18x,2hr1,//)
 1130 format (3i5,2f10.4)
 1140 format (//,11x,5hnfile,3x,i4,/,11x,4hnopt,4x,i4,//)
 1150 format (8f10.3)
c
      stop
      end

      subroutine read2 (tmax,pmax,pmin,rrmax,rrmin)
c
c
c        *********************************************************
c        ******************  read2 (modified read1) **************
c        *********************************************************
c
c
c           tmax  **  maximum value of theta (body coordinate)
c           pmax  **  maximun value of phi (body coordinate)
c           pmin  **  minimum value of phi (body coordinate)
c
c
c              n  **  total number of points
c             nb  **  number of boundary points
c
c
c
      common /body/ nb, iprint, thetas, thetan, bound(10000,2)
c
      read (4,130) n, rrmax, rrmin
      tmax=0.0
      pmax=0.0
      pmin=10000.
      nb=0
      do 120 k=1,n
         read (4,140) theta,phi
         if (theta.gt.tmax) tmax=theta
         if (phi.gt.pmax) pmax=phi
         if (phi.lt.pmin) pmin=phi
c
         nb=nb+1
         bound(nb,1)=theta
         bound(nb,2)=phi
  120 continue
      return
c
  130 format (i5,2f10.0)
  140 format (2f10.0)
c
      end
      subroutine srch1 (phi,j)
c
c
c        ********************************************************
c        *********************  srch1 (modified search)  ********
c        ********************************************************
c
c
c        this subroutine searches for boundary
c        points with the same  'phi' coordinate and groups
c        them together.
c
c             phi  **  value of 'phi' coordinate
c               n  **  total number of points
c
c               j  **  index for placing points in their proper
c                      locations
c                      starting point for search of the array
c
c        subroutines used
c             ** trans1 **
c
c
      common /body/ nb, iprint, thetas, thetan, bound(10000,2)
c
      j1=j
c
c        perform search on boundary points
c
      do 110 i=j1,nb
         abound=abs(bound(i,2)-phi)
         if (iprint.eq.2) write (6,140) i,abound,phi
         if (abound.gt.1.e-7) go to 110
         if (iprint.eq.2) write (6,150) i,j
         call trans1 (i,j)
         j=j+1
  110 continue
      return
c
c
  140 format (/,1x,2hi=,i5,2x,7habound=,f12.8,2x,4hphi=,f12.8)
  150 format (/,1x,2hi=,i5,2x,2hj=,i5)
c
      end
      subroutine trans1 (i,j)
c..modified from trans
      common /body/ nb, iprint, thetas, thetan, bound(10000,2)

c
      do 110 k=1,2
         temp=bound(i,k)
         bound(i,k)=bound(j,k)
         bound(j,k)=temp
         if (iprint.eq.2) write (6,180) j,k,bound(j,k)
  110 continue
      return
c
  180 format (/,1x,2hj=,i5,2x,2hk=,i5,2x,11hbound(j,k)=,f8.4)
c
      end
      subroutine order1 (phi,j)
c
c
c        **********************************************************
c        ********************* order1 (modified order) ************
c        **********************************************************
c
c
c        this subroutine arranges boundary             points
c        with the same 'phi' coordinate in increasing order of
c        the 'theta' coordinate.
c
c
c             phi  **  value of the 'phi' coordinate
c               j  **  location of the lowest value of 'theta'
c
c
      common /body/ nb, iprint, thetas, thetan, bound(10000,2)
c
c        determine the range of points with equal 'phi'
c        coordinates
c
      do 120 i=j,nb
         abound=abs(bound(i,2)-phi)
         if (iprint.eq.2) write (6,230) i,abound,phi
         if (abound.gt.1.e-7) go to 110
         go to 120
  110    i2=i-1
         go to 320
  120 continue
      i2=nb
      go to 320
c
c        perform the ordering
c
  320 k1=j
      i1=j
c
      do 180 i=i1,i2
         xmin=bound(i,1)
         k1=i
         do 170 k=i,i2
            if (xmin.lt.bound(k,1)) go to 170
            xmin=bound(k,1)
            k1=k
  170    continue
         call trans1 (i,k1)
  180 continue
      go to 220
c
c
  220 j=i2+1
      return
c
  230 format (/,1x,2hi=,i5,2x,7habound=,f12.8,2x,4hphi=,f12.8)
c
      end

      subroutine thetlm (phin)
c..this subroutine determines south and north theta limits for the
c..given phin
      common /body/ nb, iprint, thetas, thetan, bound(10000,2)

      do 10 i = 1, nb
         if (bound(i+1,2).ge.phin .and. bound(i,2).le.phin) go to 11
 10   continue
      write(6,*) ' something has gone wrong in thetlm, stopping exec '
      stop
 11   phi1 = bound(i,2)
      phi2 = bound(i+1,2)
      if (iprint.eq.2) write(6,*) 'phi1, phin, phi2 ', phi1, phin, phi2
      thets1 = 111.0
      thetn1 = -111.0
      thets2 = 111.0
      thetn2 = -111.0
      do 20 i = 1, nb
      if (bound(i,2).eq.phi1.and.bound(i,1).le.thets1) thets1=bound(i,1)
      if (bound(i,2).eq.phi1.and.bound(i,1).ge.thetn1) thetn1=bound(i,1)
      if (bound(i,2).eq.phi2.and.bound(i,1).le.thets2) thets2=bound(i,1)
      if (bound(i,2).eq.phi2.and.bound(i,1).ge.thetn2) thetn2=bound(i,1)
 20   continue
      thetas = thets1 + (phin-phi1) * (thets2-thets1) / (phi2-phi1)
      thetan = thetn1 + (phin-phi1) * (thetn2-thetn1) / (phi2-phi1)
      if (iprint.eq.2) write(6,*) ' thetas, thetan = ', thetas, thetan
      return
      end
      subroutine xchang (phi1,phi2,iflag,nd,nlim,iset)
c
c        *********************************************************
c        ******************   xchang (modified from change)*******
c        *********************************************************
c
c
c        this subroutine calculates the location of the gaussian
c        nodes between the limits of integration
c
c
c              iflag  **  variable option
c                             iflag=1, 'r' variable
c                             iflag=2, 'theta' variable
c                             iflag=3, 'phi' variable
c
c                 nd  **  number of gaussian nodes
c
c               nlim  **  number of limits in the 'theta'
c                         variable
c               iset  **  index of the particular set of limits
c                         for the 'r' variable
c
c                 dx  **  integration scale factor
c          phi1,phi2  **  limits of integraton in the 'phi'
c                         variable
c
c
      common /body/ nb, iprint, thetas, thetan, bound(10000,2)
      common /nodes/ phi(32),theta(64),r(32)
      common /coeffs/ cr1(32),cr2(32),ct1(32),ct2(32),cp1(32),cp2(32)
      common /fxlim/ nblim,thmax,thmin,phmax,phmin,rrmax,rrmin
      dimension dx(1), r1(1), r2(1)
c
c
      if (iflag.ne.1) go to 130
cc    if (nblim.eq.0) go to 110
      r2(iset)=rrmax
      r1(iset)=rrmin
  110 dx(1)=r2(iset)-r1(iset)
      u=0.5*(r2(iset)+r1(iset))
      v=0.5*dx(1)
      do 120 i=1,nd
         r(i)=v*cr2(i)+u
  120 continue
      return
c
c
  130 if (iflag.eq.3) go to 180
      i=1
      do 170 j=1,nlim,2
         if (nblim.eq.1) go to 140
         theta1=thetas
         theta2=thetan
         go to 150
  140    theta1=thmin
         theta2=thmax
  150    write (6,200)
         write (6,210) theta1,theta2
         dx(j)=theta2-theta1
         v=0.5*dx(j)
         u=0.5*(theta2+theta1)
         write (6,220)
         do 320 i1=1,nd
            theta(i)=v*ct2(i1)+u
            write (6,210) theta(i)
            i=i+1
  320    continue
  170 continue
      return
c
c
  180 dx(1)=phi2-phi1
      v=0.5*dx(1)
      u=0.5*(phi2+phi1)
      do 190 i=1,nd
         phi(i)=v*cp2(i)+u
  190 continue
      return
c
  200 format (/,1x,12htheta limits)
  210 format (2(5x,f10.5))
  220 format (/,1x,11htheta nodes)
c
      end
      subroutine sum1 (ntheta,nphi,i,j,k)
c
c        ***********************************************************
c        **********************  sum1  *****************************
c        ***********************************************************
c
c        this subroutine performs the gaussian summation
c
c
c           ntheta,nphi  **  dimensions of the observation grid
c                 i,j,k  **  index of the gaussian node
c
c
      common /nodes/ phi(32),theta(64),r(32)
c
      common /coeffs/ cr1(32),cr2(32),ct1(32),ct2(32),cp1(32),cp2(32)
      common tsum(1801,1801),rsum(1801,1801),gfield(1801,1801),grav(1801
     &,1801)
c
c
      do 120 m=1,ntheta
         do 110 n=1,nphi
            rsum(m,n)=rsum(m,n)+cp1(i)*ct1(j)*cr1(k)*grav(m,n)
  110    continue
  120 continue
      return
c
      end
      subroutine gravs (r,theta,phi,r1,theta1,phi1,htheta,hphi,ntheta,np
     1hi,rho)
c
c
c        ************************************************************
c        ********************        gravs       ********************
c        ************************************************************
c
c
c        *************************************************************
c        this subroutine calculates the radial component (w.r.t. sphere
c        center)  of the gravitational field on the gridded region
c        of a sphere due to a single point source
c        *************************************************************
c
c
c        *************************************************************
c        calculated from the formula for the grav field in spherical
c        coordinates
c        *************************************************************
c
c
c        *************************************************************
c           r1, theta1, phi1 **  spherical coordinates of the point
c                                source
c              r, theta, phi **  spherical coordinates of the point of
c                                observation ( kilo, rad, rad)
c
c                     htheta **  grid spacing of the theta variable
c                       hphi **  grid spacing of the phi variable
c
c                     ntheta **  number of sample points in the theta
c                                variable
c                       nphi **  number of sample points in the phi
c                                variable
c
c                       grav **  array which stores the values of the
c                                calculated gravity field
c
c                        rho **  mass density
c        **************************************************************
c
c
c        units of gravity are in milligals
c
      common tsum(1801,1801),rsum(1801,1801),gfield(1801,1801),grav(1801
     &,1801)
      g=6.67e-3
      g=g*1000.
      yyyyyy=theta
      zzzzzz=phi
      xxxxxx=phi
      cost1=cos(theta1)
      sint1=sin(theta1)
      if (abs(rho).lt.1.e-20) rho=rho1(r1,theta1,phi1)
c
c
      do 120 m=1,ntheta
         phi=xxxxxx
         do 110 n=1,nphi
c
c        calculate the vertical component of gravity
c
c
            sint=sin(theta)
            cost=cos(theta)
            cosp=cos(phi-phi1)
            cosd=cost1*cost+sint1*sint*cosp
            scale=r1*r1*sint1
            r2=r**2+r1**2-2.0*r*r1*cosd
            r15=r2**1.5
            grav(m,n)=g*rho*scale*(r-r1*cosd)/r15
c
c
            phi=phi+hphi
  110    continue
         theta=theta+htheta
  120 continue
      theta=yyyyyy
      phi=zzzzzz
      return
c
      end
      function rho1(r,theta,phi)
      rho1=r+theta+phi
      return
c
      end
      subroutine zero (m,n,iflag)
      common tsum(1801,1801),rsum(1801,1801),gfield(1801,1801),grav(1801
     &,1801)
c
      if (iflag) 110,140,170
c
  110 do 130 i=1,m
         do 120 j=1,n
            gfield(i,j)=0.0
  120    continue
  130 continue
      return
c
  140 do 320 i=1,m
         do 150 j=1,n
            rsum(i,j)=0.0
  150    continue
  320 continue
      return
c
  170 do 190 i=1,m
         do 180 j=1,n
            tsum(i,j)=0.0
  180    continue
  190 continue
      return
c
      end
      subroutine magssphere(r,theta,phi,rx,thetax,phix,htheta,hphi,
     1ntheta,nphi,k,nopt,rsphere)
c
c        ************************************************************
c        ********************       magssphere       ****************
c        ************************************************************
c
c
c        *************************************************************
c        this subroutine calculates the z, theta, or phi component on the
c        gridded region of a sphere due to a single dipole
c        *************************************************************
c
c
c        *************************************************************
c        calculated from the gradient of the scalar magnetic
c        potential in rectangular coordinates by transforming
c                  spherical # reactangular
c        *************************************************************
c
c
c        *************************************************************
c        rx, thetax, phix  **  sphere coordinates of the dipole
c                              (kilo.,deg.,deg.)
c           r, theta, phi  **  spherical coordinates of the point of
c                              observation ( kilo, rad, rad)
c
c
c                  htheta  **  grid spacing of the theta variable
c                    hphi  **  grid spacing of the phi variable
c
c             ntheta,nphi  **  number of sample points in the theta
c                              and phi variables
c
c                    tmag  **  array which stores the values of the
c                              selected magnetic field component
c
c                  k * f1 = magnetization in nT equivalent (i.e., not in emu/cc)
c
c
c                    nopt  **  magnetic field component option
c                                  nopt = 1, total field
c                                       = 2, vertical comp. ('r')
c                                       = 3, y-component ('phi')
c                                       = 4, x-component ('theta')
c                                       = 5, NOT ALLOWED
c
c          !!note!!  **  units of magnetic field are in gammas (nT)
c
c
c        **************************************************************
c
      common tsum(1801,1801),rsum(1801,1801),gfield(1801,1801),tmag(1801
     &,1801)
      common /rem/ nrem,rf1,ri1,rd1

      real i1,k
      real jr,jtheta,jphi
      real jx,jy,jz
c
c
c..sphere radius passed as argument (km)
c
c     rsphere passed as argument
c
c
  110 pi=3.1415926
      fact=pi/180.
      r1=rx+rsphere
      theta1=(90.-thetax)*fact
      phi1=phix*fact
      yyyyyy=theta
      zzzzzz=phi
      xxxxxx=phi
      cost1=cos(theta1)
      sint1=sin(theta1)
      cosp1=cos(phi1)
      sinp1=sin(phi1)
      scale=r1*r1*sint1

      i1=ri1
      d1=rd1
      f1=rf1
c
c        expression of the geomagnetic field vector in the
c        orthogonal basis er, etheta, ephi
c
      jr=f1*k*sin(i1)
      jtheta=f1*k*cos(i1)*cos(d1)
      jphi=f1*k*cos(i1)*sin(d1)
c
c
c        coordinate change of dipole position
c                  spherical # rectangular
c
      x1=r1*sint1*cosp1
      y1=r1*sint1*sinp1
      z1=r1*cost1
c
c        coordinate change of the magnetization vector by a
c        change of basis
c                  er, etheta, ephi # i, j, k
c        at the point (r1,theta1, phi1)
c
      jx=sint1*cosp1*jr+cost1*cosp1*jtheta-sinp1*jphi
      jy=sint1*sinp1*jr+cost1*sinp1*jtheta+cosp1*jphi
      jz=cost1*jr-sint1*jtheta
c
c
c     write (6,*) 'ntheta, nphi = ', ntheta, nphi

      do 190 m=1,ntheta
         phi=xxxxxx
         do 180 n=1,nphi
            sint=sin(theta)
            cost=cos(theta)
            sinp=sin(phi)
            cosp=cos(phi)
c
c        coordinate change of observation point
c                  spherical # rectangular
c
            x=r*sint*cosp
            y=r*sint*sinp
            z=r*cost
c
            xdif=x1-x
            ydif=y1-y
            zdif=z1-z
            r2=(xdif)**2+(ydif)**2+(zdif)**2
            r25=r2**2.5
c
            xmag=jx*(3.0*(xdif)**2-r2)
            xmag=xmag+jy*(xdif)*(ydif)*3.0
            xmag=xmag+jz*(xdif)*(zdif)*3.0
c
c
            ymag=jy*(3.0*(ydif)**2-r2)
            ymag=ymag+jx*(ydif)*(xdif)*3.0
            ymag=ymag+jz*(ydif)*(zdif)*3.0
c
c
            zmag=jz*(3.0*(zdif)**2-r2)
            zmag=zmag+jx*(zdif)*(xdif)*3.0
            zmag=zmag+jy*(zdif)*(ydif)*3.0
c
c
            ur=0.
            utheta=0.
            uphi=0.
            if (nopt.eq.1) go to 178
            if (nopt.eq.2) ur=1.0
            if (nopt.eq.3) uphi=1.0
            if (nopt.eq.4) utheta=1.0

            ux=sint*cosp*ur+cost*cosp*utheta-sinp*uphi
            uy=sint*sinp*ur+cost*sinp*utheta+cosp*uphi
            uz=cost*ur-sint*utheta
c
c Because the anomaly and component computations are done by converting
c spherical coordinates to rectangular.  Even for components, it is important
c to do the entire dot product as below.

            tmag(m,n)=scale*(xmag*ux+ymag*uy+zmag*uz)/r25
            go to 179

 178        ur=1.0
            ux=sint*cosp*ur+cost*cosp*utheta-sinp*uphi
            uy=sint*sinp*ur+cost*sinp*utheta+cosp*uphi
            uz=cost*ur-sint*utheta
            tz=scale*(xmag*ux+ymag*uy+zmag*uz)/r25
            ur = 0.0
            uphi=1.0
            ux=sint*cosp*ur+cost*cosp*utheta-sinp*uphi
            uy=sint*sinp*ur+cost*sinp*utheta+cosp*uphi
            uz=cost*ur-sint*utheta
            ty=scale*(xmag*ux+ymag*uy+zmag*uz)/r25
            uphi = 0.0
            utheta=1.0
            ux=sint*cosp*ur+cost*cosp*utheta-sinp*uphi
            uy=sint*sinp*ur+cost*sinp*utheta+cosp*uphi
            uz=cost*ur-sint*utheta
            tx=scale*(xmag*ux+ymag*uy+zmag*uz)/r25
            tmag(m,n) = sqrt(tx**2+ty**2+tz**2)

 179        continue
c
            phi =phi+hphi
  180    continue
         theta=theta+htheta
  190 continue
c
c
      theta=yyyyyy
      phi=zzzzzz
      return
c
      end
      function sus1(r,theta,phi)
      sus1=r+theta+phi
      return
c
      end

c ------------------------------------------------------------------------
c helper to zero multi-component output accumulators
c ------------------------------------------------------------------------
      subroutine zero_outfields(nlat,nlong)
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /outfields/ brtot(1801,1801),bthetatot(1801,1801),
     &                  bphitot(1801,1801),btottot(1801,1801)
      do 9001 i=1,nlat
         do 9000 j=1,nlong
            brtot(i,j) = 0.0
            bthetatot(i,j) = 0.0
            bphitot(i,j) = 0.0
            btottot(i,j) = 0.0
 9000    continue
 9001 continue
      return
      end