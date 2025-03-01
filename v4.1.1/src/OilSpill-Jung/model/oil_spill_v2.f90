!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!										!
!	Oil Spill Modelling Code for data format 5 (SZ). 			!
!       Author: Dr. Yun C. Jung (ycjung@hhu.ac.kr)
!
!	Routines bought from PTRACK & SELFE:					!
!	cpp, quicksearch, intersect2, signa, area_coord, levels			!
!       Warning: indices in 2D arrays are not swapped (i.e., (np,nv)).		! 
!       Also interpolation is along S-coord in tracking.	        	!
!										!
!	Input: 
!          a) hgrid.ll (if ics=2 in spill.in) or hgrid.gr3 (if ics=1 in spill.in)
!          b) vgrid.in
!          c) spill.in (see format below; also sample 'spill.in.sample')
!          d) *hvel.64, *vert.63,*elev.61, *tdff.63, *wind.62 
!          e) (for advanced users) some constants specified near the beginning of the main routine 
!	Input spill.in (example script to generate this input can be found in gen_data.f90):
!	  (1) description - info only
!	  (2) nscreen - =0: no screen outputs; =1: with screen outputs 
!         (3) istiff - particles at fixed distance from f.s. (1) or not (0). 
!                      For oil spill, istiff=0 in general;
!	  (4) ics,slam0,sfea0 - ics=1 (map projection); =2 (lat/lon); (slam0,sfea0)
!                               are lon/lat of the projection center 
!	  (5) h0,rnday,dtm,nspool,ihfskip,ndeltp - min. depth, # of days, 	!
!	      time step used in the original run, nspool and ihfskip used 	!
!             in the run, # of sub-division to be used in the tracking; 
!         (6) ihdf,hdc,horcon - if ihdf=0, the constant diffusivity is given by hdc;
!                              if ihdf=1, the diffusivity is calculated from Smagorinsky
!                              scheme with the coefficient given by horcon
!         (7) ibuoy,iwind - on/off (0/1) flags for buoyancy effects and wind effects
!         (8) pbeach - percentage of refloating for stranded particles; e.g. '20' 
!                     means 20% of refloating rate
!	  (9) nparticle - # of particles;			
!         (10) (For each particle) idp(i),st_p(i),xpar(i),ypar(i),zpar0(i) -
!                                 particle id, start time (sec),starting x,y, and 
!                                 z relative to the instant f.s. (<=0). 
!							
!	Outputs: spill.out (drogue format; see plot/plot_oil.m for details); fort.16 (run info), fort.11 (fatal errors)
!										!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ...  code description
!      oil_spill_v2.f90 modified on 01/10/2011 from oil_spill_v1.f90
!      added a routine for getting wind & elevation at spill location 
!      added a routine for computing the rising velocity based on Proctor et al
!      added the smagorinsky algorithm for hvis_e & hvis_p (@element and nodes)
!      added minor changes to get the code simplified

!Compile:
!ifort -O3 -assume byterecl -Bstatic -o oil_spill_v2 oil_spill_v2.f90

!...  Data type consts
      module kind_par
        implicit none
        integer, parameter :: sng_kind1=4
        integer, parameter :: dbl_kind1=8
        real(kind=dbl_kind1), parameter :: small1=1.e-6 
!       small is non-negative number; must be identical to that in global
      end module kind_par

!...  definition of variables
!************************************************************************
!     			mnp < mne < mns					*
!************************************************************************
      module global
        implicit none
        integer, parameter :: sng_kind=4
        integer, parameter :: dbl_kind=8
        
!...  	Dimensioning parameters
        integer, parameter :: mnp=35000
        integer, parameter :: mne=66000
        integer, parameter :: mns=100000
        integer, parameter :: mnv=85
        integer, parameter :: mnei=40 !neighbor
        integer, parameter :: nbyte=4
        integer, parameter :: mnout=100 !max. # of output files
      	integer, parameter :: mnope=6 !# of open bnd segements
        integer, parameter :: mnond=1000 !max. # of open-bnd nodes on each segment
      	integer, parameter :: mnland=50 !# of land bnd segements
      	integer, parameter :: mnlnd=10000 !max. # of land nodes on each segment
!      	integer, parameter :: mnoe=20000 !max. # of open-bnd elements on each segment
!      	integer, parameter :: mnosd=20000 !max # of open-bnd sides on each segment
!      	integer, parameter :: mnbfr=9 !# of forcing freqs.
!      	integer, parameter :: itmax=5000 !# of iteration for itpack solvers used for dimensioning
!      	integer, parameter :: nwksp=6*mne+4*itmax !available work space for itpack solvers
!      	integer, parameter :: mirec=1109000000 !997000000) !max. record # to prevent output ~> 4GB
      	real(kind=dbl_kind), parameter :: zero=1.e-5 !small postive number in lieu of 0; usually used to check areas 
        real(kind=dbl_kind), parameter :: small1=1.e-6 !small non-negative number; must be identical to that in kind_par
        real(kind=dbl_kind), parameter :: pi=3.1415926d0 

!...  	Important variables
      	integer :: np,ne,ns,nvrt,istiff,ivcor,kz,nsig
      	real(kind=dbl_kind) :: h0,rho0,dt,h_c,s_con1,theta_b,theta_f,h_s,rnday

!...    Output handles
        character(len=48) :: start_time,version
        character(len=12) :: ifile_char
        character(len=48), dimension(mnout) :: outfile,variable_nm,variable_dim
        integer :: ihot,nrec,nspool,igmp,noutgm,ifile,noutput,nparticle,ifort12(100)
        integer, dimension(mnout) :: ichan,irec,iof
!        real(kind=dbl_kind), dimension(mnout) :: vpos      

!...  1D arrays
        integer :: nne(mnp),nnp(mnp),idry(mnp),idry_e(mne),idry_e0(mne),iback(mnp)
        integer :: kbp(mnp),kbs(mns),kbe(mne),kbp00(mnp),isbnd(mnp)
        real(kind=dbl_kind), dimension(mnp) :: x,y,dp,hmod,eta1,eta2,eta3,wnx1,wnx2, &
            & wny1,wny2
        real(kind=dbl_kind), dimension(mne) :: area,radiel,xctr,yctr
        real(kind=dbl_kind), dimension(mns) :: snx,sny,distj,xcj,ycj,dps
        real(kind=dbl_kind) :: sigma(mnv),cs(mnv),dcs(mnv),ztot(mnv)
        real(kind=dbl_kind) :: zpar0(10*mnp)

!...  2D arrays
        integer :: elnode(3,mne),nx(3,2),ic3(3,mne),indel(mnei,mnp),elside(3,mne),isdel(2,mns), &
           &isidenode(2,mns),iself(mnei,mnp)
        real(kind=dbl_kind) :: ssign(3,mne),z(mnp,mnv),icum1(mnp,mnv),icum2(mnp,mnv,2)
        real(kind=dbl_kind), dimension(mnp,mnv) :: hf1,vf1,uu1,vv1,ww1,hf2,vf2,uu2,vv2,ww2
        real(kind=dbl_kind), dimension(mne,mnv) :: dudx_e,dudy_e,dvdx_e,dvdy_e,hvis_e
        real(kind=dbl_kind) :: hvis_p(mnp,mnv)
      end module global

!...  Main program
      program ptrack
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      real(kind=sng_kind) :: floatout,floatout2
      real(kind=dbl_kind), dimension(10*mnp) :: xpar,ypar,zpar,st_p,tpar,upar,vpar,wpar, &
          &dhfx,dhfy,dhfz,grdx,grdy,grdz,amas,wndx,wndy
      integer, dimension(10*mnp) :: ielpar,levpar,iabnorm,ist,inp
      character(len=25) :: idp(10*mnp)
      real(kind=dbl_kind) :: staint(4),dx(10),dy(10),dz(10)
      real xcp(3),ycp(3),lng(10),psa(10)
      dimension arco(3),ztmp(mnv),ztmp2(mnv)
      dimension nond(mnope),iond(mnope,mnond),nlnd(mnland),ilnd(mnland,mnlnd)
      data iseed/5/
! .......................................................................................... 
      do i=1,3
        do j=1,2
          nx(i,j)=i+j
          if(nx(i,j)>3) nx(i,j)=nx(i,j)-3
          if(nx(i,j)<1.or.nx(i,j)>3) then
            write(*,*)'nx wrong',i,j,nx(i,j)
            stop
          endif
        enddo !j
      enddo !i

!...  Initialize arrays and variables
      isbnd=0
!...  Initialize levpar as flag
      levpar=-99
!...  ist=0(not active), 1(active), -1(stranded), 2&3(out of boundary)
      ist=0      ! particle status flag
      inp=0      ! neighboring element # for refloating  
      amas=1.0   ! initial mass of particles
      iabnorm=0  ! abnormal tracking exit flag
      ifort12=0

!...  Open input file
      open(90,file='spill.in',status='old')
      open(95,file='spill.out',status='unknown')

!...  Read in spill data
      read(90,*)         ! read data_format
      read(90,*) nscreen
      read(90,*) istiff  !1: fixed distance from F.S.
      if(istiff/=0.and.istiff/=1) then
        write(*,*)'Wrong istiff',istiff
        stop
      endif
      read(90,*) ics,slam0,sfea0
      slam0=slam0*pi/180
      sfea0=sfea0*pi/180
      read(90,*) h0,rnday,dtm,nspool,ihfskip,ndeltp !# of sub-divisions
      if(mod(ihfskip,nspool).ne.0) then
        write(*,*)'ihfskip must be a multiple of nspool'
        stop
      endif

! ........................................................
! ... Description of parameters
!     ihdf  : turn Smagorinsky algorithm - off(0), on(1)
!     ibuoy : turn buoyancy of particle  - off(0), on(1)
!     iwind : turn wind effect - off(0), on(1)
!     pbch  : set percentage of stranding on shore
! ........................................................
      read(90,*) ihdf,hdc,horcon
      if(ihdf.eq.0) then
         hf2=hdc     ! constant diffusivity 
      endif

      read(90,*) ibuoy,iwind
      read(90,*) pbeach

! .......................................................................
! ... treatment for buoyant oil particles 
      if(ibuoy.eq.1) then
!...  compute the rising velocity(m/s) based on Proctor et al., 1994 
        gr=9.8               ! m/s^2
        rho_o=900.0d3        ! kg/m^3
        rho_w=1025.0d3       ! kg/m^3
        di=500.0d-6          ! m
        smu=1.05d-6          ! m^2/s  
! ... critical diameter, dc
        dc=9.52*smu**(2./3.)/(gr**(1./3.)*(1.-rho_o/rho_w)**(1./3.)) !m     
!... compute the rising velocity, m/s
         if(di.ge.dc) then 
           rsl=sqrt(8./3.*gr*di*(1.-rho_o/rho_w))  ! large droplet
         else
           rsl=gr*di**2*(1.-rho_o/rho_w)/(18.*smu) ! small droplet
         endif
      else if(ibuoy.eq.0) then
        rsl=0.0
      else
        write(11,*)'ibuoy is incorrect'
        stop
      endif
! .......................................................................

      read(90,*) nparticle
      if(nparticle.gt.10*mnp) then
        write(11,*)'nparticle > 10*mnp'
        stop
      endif
      dt=dtm*nspool !output time step
      st_m=rnday*86400 ! max. running time
      do i=1,nparticle
        if(ics.eq.1) then
!	  zpar0: relative to f.s.
          read(90,*)idp(i),st_p(i),xpar(i),ypar(i),zpar0(i)
        else
          read(90,*)idp(i),st_p(i),xparl,yparl,zpar0(i)
          xparl=xparl/180*pi
          yparl=yparl/180*pi
          call cpp(xpar(i),ypar(i),xparl,yparl,slam0,sfea0)
        endif
        if(st_p(i)<0.or.st_p(i)>rnday*86400) then
          write(11,*)'Starting time for particle ',i,' out of range:',st_p(i)
          stop
        endif
        if(zpar0(i)>0) then
          write(11,*)'Starting z-coord above f.s.',i
          stop
        endif
        if(st_p(i)<st_m) st_m=st_p(i)
      enddo !i
      close(90)

!...  set spill location
      x_spill=xpar(1); y_spill=ypar(1)

!...  Read in h- and v-grid and compute geometry
      open(19,file='vgrid.in',status='old')
      ivcor=2 !S only
      read(19,*) nvrt,kz,h_s !kz>=1
      if(nvrt>mnv.or.nvrt<3) then
        write(11,*)'nvrt > mnv or nvrt<3'
        stop
      endif
      if(kz<1.or.kz>nvrt-2) then
        write(11,*)'Wrong kz:',kz
        stop
      endif
      if(h_s<10) then
        write(11,*)'h_s needs to be larger:',h_s
        stop
      endif

!     # of z-levels excluding "bottom" at h_s
      read(19,*) !for adding comment "Z levels"
      do k=1,kz-1
       read(19,*)j,ztot(k)
       if(k>1.and.ztot(k)<=ztot(k-1).or.ztot(k)>=-h_s) then
         write(11,*)'z-level inverted:',k
         stop
       endif
      enddo !k

      read(19,*) !level kz     
!     In case kz=1, there is only 1 ztot(1)=-h_s
      ztot(kz)=-h_s
      nsig=nvrt-kz+1 !# of S levels (including "bottom" & f.s.)
      read(19,*) !for adding comment "S levels"
      read(19,*)h_c,theta_b,theta_f
      if(h_c<5) then !large h_c to avoid 2nd type abnormaty
        write(11,*)'h_c needs to be larger:',h_c
        stop
      endif
      if(theta_b<0.or.theta_b>1) then
        write(11,*)'Wrong theta_b:',theta_b
        stop
      endif
      if(theta_f<=0) then 
        write(11,*)'Wrong theta_f:',theta_f 
        stop
      endif

!     Pre-compute constants
      s_con1=dsinh(theta_f)
      sigma(1)=-1 !bottom
      sigma(nsig)=0 !surface
      read(19,*) !level kz
      do k=kz+1,nvrt-1
        kin=k-kz+1
        read(19,*) j,sigma(kin)
        if(sigma(kin)<=sigma(kin-1).or.sigma(kin)>=0) then
          write(11,*)'Check sigma levels at:',k,sigma(kin),sigma(kin-1)
          stop
        endif
      enddo
      read(19,*) !level nvrt
      close(19)

!     Compute C(s) and C'(s)
      do k=1,nsig
        cs(k)=(1-theta_b)*dsinh(theta_f*sigma(k))/dsinh(theta_f)+ &
     &theta_b*(dtanh(theta_f*(sigma(k)+0.5))-dtanh(theta_f*0.5))/2/dtanh(theta_f*0.5)
        dcs(k)=(1-theta_b)*theta_f*dcosh(theta_f*sigma(k))/dsinh(theta_f)+ &
     &theta_b*theta_f/2/dtanh(theta_f*0.5)/dcosh(theta_f*(sigma(k)+0.5))**2
      enddo !k=1,nvrt

      if(ics.eq.1) then
         open(14,file='hgrid.gr3',status='old')
        else
         open(14,file='hgrid.ll', status='old')
      endif    

      read(14,*) 
      read(14,*) ne,np
      if(ne.gt.mne.or.np.gt.mnp) then
        write(11,*)'Increase mne/mnp',mne,mnp,ne,np
        stop
      endif

      do i=1,np
        if(ics.eq.1) then
          read(14,*) j,x(i),y(i),dp(i)
         else
          read(14,*) j,xlon,ylat,dp(i)
          ylat=ylat/180*pi
          xlon=xlon/180*pi
          call cpp(x(i),y(i),xlon,ylat,slam0,sfea0)
        endif
        hmod(i)=dmin1(dp(i),h_s)
      enddo !i=1,np

      do i=1,ne
        read(14,*) j,l,(elnode(l,i),l=1,3)
        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        area(i)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
        if(area(i)<=0) then
          write(11,*)'Negative area at',i
          stop
        endif
      enddo !i=1,ne

!     Open bnds
      read(14,*) nope
      if(nope>mnope) then
        write(11,*) 'nope > mnope' 
        stop
      endif

      read(14,*) neta
      ntot=0
      do k=1,nope
        read(14,*) nond(k)
        if(nond(k)>mnond) then
          write(11,*) 'nond(k) > mnond'
          stop
        endif
        do i=1,nond(k)
          read(14,*) iond(k,i)
          isbnd(iond(k,i))=k
        enddo
        if(iond(k,1)==iond(k,nond(k))) then
          write(11,*)'Looped open bnd:',k
          stop
        endif
        ntot=ntot+nond(k)
      enddo

      if(neta/=ntot) then
        write(11,*)'neta /= total # of open bnd nodes',neta,ntot
        stop
      endif

!     Land bnds
      read(14,*) nland
      if(nland>mnland) then
        write(11,*) 'nland > mnland'
        stop
      endif

      read(14,*) nvel
      do k=1,nland
        read(14,*) nlnd(k)
        if(nlnd(k)>mnlnd) then
          write(11,*)'nlnd(k) > mnlnd',k,nlnd(k),mnlnd
          stop
        endif
        do i=1,nlnd(k)
          read(14,*) ilnd(k,i)
          if(isbnd(ilnd(k,i))==0) isbnd(ilnd(k,i))=-1 !overlap of open bnd
        enddo
      enddo !k=1,nland
      close(14)

!... isbnd(node#) means land(-1),inside(0),ocean bnd(2),river bnd(1)    
!...  End fort.14
!                                                                             *
!                                                                             *
!******************************************************************************
!                                                                             *
!     			Compute geometry 				      *
!                                                                             *
!******************************************************************************
!                                                                             *
!                                                                             *

!...  compute the elements connected to each node 
      do i=1,np
        nne(i)=0
      enddo

      do i=1,ne
        do j=1,3
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd).gt.mnei) then
            write(11,*)'Too many neighbors',nd
            stop
          endif
          indel(nne(nd),nd)=i
          iself(nne(nd),nd)=j
        enddo
      enddo

!...  Check hanging nodes
!...
      ihang=0
      do i=1,np
        if(nne(i).eq.0) then
          ihang=1
          write(11,*)'Hanging node',i
        endif
      enddo
      if(ihang.eq.1) then
        write(11,*)'Hanging nodes found'
        stop
      endif
      
!     Compute ball info
      do i=1,ne
        do j=1,3
          ic3(j,i)=0 !index for bnd sides
          nd1=elnode(nx(j,1),i)
          nd2=elnode(nx(j,2),i)
          do k=1,nne(nd1)
            iel=indel(k,nd1)
            if(iel/=i.and.(elnode(1,iel)==nd2.or.elnode(2,iel)==nd2 &
                     .or.elnode(3,iel)==nd2)) then
              ic3(j,i)=iel
            end if
          enddo !k
        enddo !j
      enddo !i

!... ic3(j,i) shows the element # of neighbouring element.
!... ic3(j,i) = 0 shows no neighbouring element, it is boundary.

!...  compute the sides information
!...
      ns=0 !# of sides
      do i=1,ne
        do j=1,3
          nd1=elnode(nx(j,1),i)
          nd2=elnode(nx(j,2),i)
          if(ic3(j,i)==0.or.i<ic3(j,i)) then !new sides
            ns=ns+1
            if(ns.gt.mns) then
              write(11,*)'Too many sides'
              stop
            endif
            elside(j,i)=ns
            isdel(1,ns)=i
            isidenode(1,ns)=nd1
            isidenode(2,ns)=nd2
!            x1j(ns)=x(nd1)
!            y1j(ns)=y(nd1)
!            x2j(ns)=x(nd2)
!            y2j(ns)=y(nd2)
            xcj(ns)=(x(nd1)+x(nd2))/2
            ycj(ns)=(y(nd1)+y(nd2))/2
            dps(ns)=(dp(nd1)+dp(nd2))/2
            distj(ns)=dsqrt((x(nd2)-x(nd1))**2+(y(nd2)-y(nd1))**2)
            if(distj(ns)==0) then
              write(11,*)'Zero side',ns
              stop
            endif
            thetan=datan2(x(nd1)-x(nd2),y(nd2)-y(nd1))
            snx(ns)=dcos(thetan)
            sny(ns)=dsin(thetan)
            isdel(2,ns)=ic3(j,i) !bnd element => bnd side
!       Corresponding side in element ic3(j,i) 
            if(ic3(j,i)/=0) then !old internal side
              iel=ic3(j,i)
              index=0
              do k=1,3
                if(ic3(k,iel)==i) then
                  index=k
                  exit
                endif
              enddo !k
              if(index==0) then
                write(*,*)'Wrong ball info',i,j
                stop
              endif
              elside(index,iel)=ns
            endif !ic3(j,i).ne.0
          endif !ic3(j,i)==0.or.i<ic3(j,i)
        enddo !j=1,3
      enddo !i=1,ne

      if(ns.lt.ne.or.ns.lt.np) then 
        write(11,*)'Weird grid with ns < ne or ns < np',np,ne,ns
        stop
      endif

      do i=1,ne
        do j=1,3
          jsj=elside(j,i)
          ssign(j,i)=(isdel(2,jsj)-2*i+isdel(1,jsj))/(isdel(2,jsj)-isdel(1,jsj))
        enddo
      enddo

      if(nscreen.eq.1) write(*,*) 'There are',ns,' sides in the grid...'
      write(16,*) 'There are',ns,' sides in the grid...'

!...  compute centers of each triangle 
!...
      do i=1,ne
        xctr(i)=0
        yctr(i)=0
        do j=1,3
          xctr(i)=xctr(i)+x(elnode(j,i))/3
          yctr(i)=yctr(i)+y(elnode(j,i))/3
        enddo !j
      enddo !i=1,ne

      if(nscreen.eq.1) write(*,*)'done computing geometry...'
      write(16,*)'done computing geometry...'

!... isd() represents the side# that a particle is crossed.
!... isdenode(isd,1&2) represents both node# that a particle is crossed.

!...  Read in header 
      nt=rnday*86400/dt !total # of iterations
      ifile=1 !for st_m=0
      do i=1,nt/(ihfskip/nspool)+1
        if(st_m/dt>(i-1)*(ihfskip/nspool).and.st_m/dt<=i*(ihfskip/nspool)) then
          ifile=i
          exit
        endif
      enddo
      iths=(ifile-1)*ihfskip/nspool+1 !iteration # for the 1st time step output in ifile
      write(ifile_char,'(i12)') ifile
      open(60,file=ifile_char//'_elev.61',access='direct',recl=nbyte)
      open(61,file=ifile_char//'_wind.62',access='direct',recl=nbyte)
      open(62,file=ifile_char//'_vert.63',access='direct',recl=nbyte)
      open(63,file=ifile_char//'_tdff.63',access='direct',recl=nbyte)
      open(64,file=ifile_char//'_hvel.64',access='direct',recl=nbyte)
      irec00=5*48/nbyte+5+7+nvrt+2 !elev

!     Read initial bottom index
      do m=1,np
        read(60,rec=irec00+4)kbp00(m)
        irec00=irec00+4
      enddo !m=1,np
      irec00=irec00+4*ne

!     Initialize kbp for levels()
      kbp=kbp00

!...  Compute record # offset for a node and level for 3D outputs
!...
      icount1=0
      icount2=0
      do i=1,np
        do k=max0(1,kbp00(i)),nvrt
          do m=1,2
            icount2=icount2+1
            icum2(i,k,m)=icount2
          enddo !m
          icount1=icount1+1
          icum1(i,k)=icount1
        enddo !k
      enddo !i=1,np

      irec01=irec00 ! end of the header for elev.61
      irec02=irec00 ! wind.62
      irec03=irec00 ! vert.63
      irec04=irec00 ! tdff.63
      irec05=irec00 ! hvel.64
      irec1=irec01
      irec2=irec02
      irec3=irec03
      irec4=irec04
      irec5=irec05

!...  Compute initial elements for particle tracking
!...
      lp1: do i=1,nparticle
        do k=1,ne
          aa=0
          do j=1,3
            n1=elnode(j,k)
            n2=elnode(nx(j,1),k)
            aa=aa+dabs(signa(x(n1),x(n2),xpar(i),y(n1),y(n2),ypar(i)))
          enddo !j
          ae=dabs(aa-area(k))/area(k)

          if(ae.lt.1.e-5) then
!          if(ae<small) then
            ielpar(i)=k
            cycle lp1
          endif
        enddo !k=1,ne
        write(11,*)'Cannot find init. element for particle',i
        stop
      end do lp1 !i=1,nparticle

! ... set element# of spill location
      ne_spill=ielpar(1)

      write(95,*)'Drogues'
      write(95,*) nt-iths+1
      if(nscreen.eq.1) write(*,*)'done initialization...'
      write(16,*)'done initialization...'

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!									!
!	       Time iteration						!
!									!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      it2=nt
!     Initialize for output before particle moving
      upar=0; vpar=0; wpar=0
      zpar=zpar0

      do it=iths,it2
!--------------------------------------------------------------------------
      time=it*dt
      
!...  Read in elevation and vel. info
      if(it>ihfskip/nspool*ifile) then
        ifile=ifile+1
        write(ifile_char,'(i12)') ifile
        close(60)
        close(61)
        close(62)
        close(63)
        close(64)
!
        open(60,file=ifile_char//'_elev.61',access='direct',recl=nbyte)
        open(61,file=ifile_char//'_wind.62',access='direct',recl=nbyte)
        open(62,file=ifile_char//'_vert.63',access='direct',recl=nbyte)
        open(63,file=ifile_char//'_tdff.63',access='direct',recl=nbyte)
        open(64,file=ifile_char//'_hvel.64',access='direct',recl=nbyte)
        irec1=irec01
        irec2=irec02
        irec3=irec03
        irec4=irec04
        irec5=irec05
      endif 
        irec1=irec1+np+2
        irec2=irec2+np+2
        irec3=irec3+np+2
        irec4=irec4+np+2      
        irec5=irec5+np+2

        do i=1,np
          read(60,rec=irec1+1) floatout
          irec1=irec1+1
          eta2(i)=floatout
          read(61,rec=irec2+1) floatout
          read(61,rec=irec2+2) floatout2
          irec2=irec2+2
          wnx2(i)=floatout
          wny2(i)=floatout2
          do k=max0(1,kbp00(i)),nvrt
            read(62,rec=irec3+1) floatout
            irec3=irec3+1
            ww2(i,k)=floatout
            read(63,rec=irec4+1) floatout
            irec4=irec4+1
            vf2(i,k)=floatout
            read(64,rec=irec5+1) floatout
            read(64,rec=irec5+2) floatout2
            irec5=irec5+2
            uu2(i,k)=floatout
            vv2(i,k)=floatout2
          enddo !k
        enddo !i 

! ..............................................................
! ... get wind & elevation at spill location       
      call area_coord(ne_spill,x_spill,y_spill,arco)                
      windx=0;windy=0;selev=0
!      dep=0 
      do j=1,3
         nd=elnode(j,ne_spill)
         windx=windx+wnx2(nd)*arco(j)
         windy=windy+wny2(nd)*arco(j)
         selev=selev+eta2(nd)*arco(j)   
!         dep=dep+dp(nd)*arco(j)   
      enddo !j
!         print*,ne_spill,dep
!         pause
! ..............................................................

!...  Store info for first step
      if(it==iths) then
        uu1=uu2; vv1=vv2; ww1=ww2; eta1=eta2; eta3=eta2;&
        & vf1=vf2; wnx1=wnx2; wny1=wny2
      endif

!...  Compute z-cor
      call levels

! ..............................................................
!...  compute hvis_e & hvis_e based on Smagorinsky Algorithm
      if(ihdf.eq.1) then
      do k=1,nvrt 
      do i=1,ne
          if(idry_e(i).eq.1) then
            hvis_e(i,k)=0
            goto 150
          endif 
          do j=1,3 
            xcp(j)=x(elnode(j,i))
            ycp(j)=y(elnode(j,i))        
          enddo
          dux=0;duy=0;dvx=0;dvy=0
          do j=1,3
            m1=j+1;m2=j+2
            if(m1>3) m1=m1-3
            if(m2>3) m2=m2-3            
            dux=dux+uu2(elnode(j,i),k)*(ycp(m1)-ycp(m2))
            duy=duy+uu2(elnode(j,i),k)*(xcp(m2)-xcp(m1))
            dvx=dvx+vv2(elnode(j,i),k)*(ycp(m1)-ycp(m2))
            dvy=dvy+vv2(elnode(j,i),k)*(xcp(m2)-xcp(m1))
          enddo                  
            dudx_e(i,k)=1./(2.*area(i))*dux
            dudy_e(i,k)=1./(2.*area(i))*duy
            dvdx_e(i,k)=1./(2.*area(i))*dvx
            dvdy_e(i,k)=1./(2.*area(i))*dvy
            hvis_e(i,k)=horcon*area(i)*dsqrt(dudx_e(i,k)**2+dvdy_e(i,k)**2 &
                       &+.5*(dudy_e(i,k)+dvdx_e(i,k))**2)
150       continue
      enddo !i=1,ne
      enddo !k=1,nvrt

!...  compute hvis_p(i,k) from hvis_e(i,k)
      do k=1,nvrt
      do i=1,np
         if(idry(i).eq.1) then
           hvis_p(i,k)=0
           goto 160
         endif
         nec=nne(i) ! number of elements connected to node i
         spm=0;slm=0
         do j=1,nec
           nel=indel(j,i)
           lng(j)=dsqrt((x(i)-xctr(nel))**2+(y(i)-yctr(nel))**2)
           psa(j)=hvis_e(nel,k)
           slm=slm+1./lng(j)
           spm=spm+psa(j)*1./lng(j)
         enddo  
         hvis_p(i,k)=spm/slm
160      hf2(i,k)=hvis_p(i,k)
      enddo !i=1,np
      enddo !k=1,nvrt
      endif !ihdf.eq.1

!       tmax=0.0;tmin=10.0      
!       do k=1,nvrt
!       do i=1,ne
!          hh=hvis_p(i,k)
!          if(hh>tmax) tmax=hh
!          if(hh<tmin) tmin=hh
!          if(hh>46) then
!             print*,i,k,hh
!          endif
!        enddo
!        enddo
!        print*, 'tmin,tmax : ',tmin,tmax 
!        pause

! ..............................................................

!...  Store info for diffusion term in first step
      if(it==iths) then
        hf1=hf2; grdx=0.0; grdy=0.0; grdz(i)=0.0;&
        & dhfx=hdc; dhfy=hdc; dhfz=3.0d-4
      endif

      if(nscreen.eq.1) write(*,*)'begin ptrack...'
      write(16,*)'begin ptrack...'

!...  Particle tracking
      write(95,*) time,nparticle
      
      do i=1,nparticle
        eta_p=0; dp_p=0 !for output before moving
        if(time<=st_p(i)) go to 450
        if(ist(i).eq.0) ist(i)=1
!...    no refloating when the particle is on boundary  
        if(ist(i).eq.-1 .and. inp(i).eq.0) go to 115 
 
!...    Refloat the stranded particle if re-wet
        if(ist(i).eq.-1 .and. idry_e(inp(i)).eq.0) ist(i)=1
115     continue
        if(ist(i).ne.1) go to 450

        pt=dt !tracking time step
!...    Initialize starting level
        if(levpar(i)==-99) then !just started
          pt=(time-st_p(i))
          if(pt<=0) then
            write(*,*)'Tracking step negative:',pt
            stop
          endif
          iel=ielpar(i)
          if(idry_e(iel)==1) then !dry
            levpar(i)=-1
          else !wet
            call area_coord(iel,xpar(i),ypar(i),arco)
            do k=kbe(iel),nvrt
              ztmp2(k)=0
              do j=1,3
                nd=elnode(j,iel)
                ztmp2(k)=ztmp2(k)+z(nd,k)*arco(j)
              enddo !j
            enddo !k
            zpar(i)=dmax1(zpar0(i)+ztmp2(nvrt),ztmp2(kbe(iel))) !zpar0<=0
            jlev=0
            do k=kbe(iel),nvrt-1
              if(zpar(i)>=ztmp2(k).and.zpar(i)<=ztmp2(k+1)) then
                jlev=k+1
                exit
              endif
            enddo !k
            if(jlev==0) then
              write(11,*)'Cannot find an init. level:',i,zpar(i),(ztmp2(k),k=kbe(iel),nvrt)
              stop
            endif

            levpar(i)=jlev
            upar(i)=0; vpar(i)=0; wpar(i)=0; wndx(i)=0; wndy(i)=0 
            do j=1,3
              nd=elnode(j,iel)
              upar(i)=upar(i)+uu2(nd,jlev)*arco(j)
              vpar(i)=vpar(i)+vv2(nd,jlev)*arco(j)
              wpar(i)=wpar(i)+ww2(nd,jlev)*arco(j)
              wndx(i)=wndx(i)+wnx2(nd)*arco(j)
              wndy(i)=wndy(i)+wny2(nd)*arco(j)
            enddo !j  
          endif !wet
        endif !levpar=-99

!	Wetting/drying
        if(idry_e(ielpar(i))==1) then
          levpar(i)=-1
          go to 450
        endif

        nnel=ielpar(i)
        jlev=levpar(i)
!       Rewetted elements
        if(jlev==-1) then !element nnel wet
          jlev=nvrt
          zpar(i)=(eta3(elnode(1,nnel))+eta3(elnode(2,nnel))+eta3(elnode(3,nnel)))/3
        endif
  
!	Tracking
        x0=xpar(i)
        y0=ypar(i)
        z0=zpar(i)
        nnel0=nnel
        jlev0=jlev
        dtb=pt/ndeltp

        do idt=1,ndeltp
          if(ist(i).ne.1) go to 404
          trat=real(idt)/ndeltp
 
!        upar, vpar, wpar : advection velocity at particle location
!        grdx, grdy, grdz : gradient of diffusivity
!        dhfx, dhfy, dhfz : diffusivity at half time step
! ...    wind rotation & apply drag_c
          rotate_angle=3.0*(pi/180.0)   ! coriolis effect
!          rotate_angle=0.0
          drag_c=0.03
          dir = atan2(wndx(i),wndy(i))
          speed = sqrt(wndx(i)**2+wndy(i)**2)
          dir = dir+rotate_angle
          wind_x = speed*sin(dir)*drag_c
          wind_y = speed*cos(dir)*drag_c

!...      turn off wind
          if(iwind.eq.0) then
            wind_x=0;wind_y=0; windx=0; windy=0
          endif  
          if(z0<-0.1) then  ! sub_surface particles are not influenced by wind
            wind_x=0;wind_y=0
          endif 

! ...    generating random number
          do k=1,3
            dx(k)=ran1(iseed)
            dy(k)=ran2(iseed)
            dz(k)=ran3(iseed)
          enddo

          rndx=dx(1)
          rndy=dx(2)
          rndz=dx(3)
          xadv=(upar(i)+wind_x+grdx(i))*dtb
          yadv=(vpar(i)+wind_y+grdy(i))*dtb
          zadv=(wpar(i)+rsl+grdz(i))*dtb        
          xdif=(2*rndx-1)*sqrt(6*dhfx(i)*dtb)
          ydif=(2*rndy-1)*sqrt(6*dhfy(i)*dtb)
          zdif=(2*rndz-1)*sqrt(6*dhfz(i)*dtb)   
!          xadv=0;yadv=0;zadv=0;zdif=0
!          if(i.eq.1) then
!            iel=nnel
!            k=jlev
!            nd1=elnode(1,iel)
!            nd2=elnode(2,iel)
!            nd3=elnode(3,iel)
!            print*,iel,idry_e(iel),hvis_e(iel,k)
!            print*,hvis_p(nd1,k),hvis_p(nd2,k),hvis_p(nd3,k)
!            print*,grdx(i),grdy(i),grdz(i)
!            print*,dhfx(i),dhfy(i),dhfz(i)
!            pause  
!          endif

          xt=x0+xadv+xdif
          yt=y0+yadv+ydif
          zt=z0+zadv+zdif    

! ...    rnds - random number used for stranding
         rnds=dy(1)
         call quicksearch(1,idt,i,nnel0,jlev0,dtb,x0,y0,z0,xt,yt,zt,nnel,jlev,arco,zrat,&
               & ss,nfl,eta_p,dp_p,ztmp,kbpl,ist,inp,rnds,pbeach)    

!        Interpolation data for next step
         call interp(i,nnel,jlev,ztmp,trat,zrat,arco,upar,vpar,wpar,dhfx,dhfy,dhfz,grdx,&
               & grdy,grdz,wndx,wndy)

          !if(iflqs1==1) then
          if(nfl==1) then
            iabnorm(i)=1
            go to 404
          endif

          x0=xt
          y0=yt
          z0=zt
          nnel0=nnel
          jlev0=jlev
        enddo !idt=1,ndeltp

404     xpar(i)=xt
        ypar(i)=yt
        zpar(i)=zt
        tpar(i)=real(zpar(i)-eta_p)
!        if(tpar(i).gt.0.0) tpar(i)=0.0    
        ielpar(i)=nnel
        levpar(i)=jlev

450     continue
! ...   mass decay - evaporation
        ddt=dt/ndeltp
        if(ist(i).ne.0) call evapor(ddt,time,st_p(i),amas(i))      

        if (ics.eq.2) then
           call cppinverse(xout, yout, xpar(i), ypar(i), slam0, sfea0)
           xout = xout * 180.0 / pi
           yout = yout * 180.0 / pi
        else
           xout = xpar(i)
           yout = ypar(i)
        endif

! ...   print out results
        if(ist(i).ne.0) then
          write(95,401) i,ist(i),amas(i),xout,yout,tpar(i)
         else
          write(95,401) i,ist(i),amas(i),xout,yout,zpar(i)
        endif  

401     format(i6,i3,2x,f7.5,1x,2f14.4,f10.4)
      enddo !i=1,nparticle

! ... save wind & elevation at spill location
      write(95,'(3f10.4)') windx,windy,selev
      write(15,'(4f10.4)') time,windx,windy,selev

!...  Store info for next step
      uu1=uu2; vv1=vv2; ww1=ww2; eta1=eta2; hf1=hf2;&
      & vf1=vf2; wnx1=wnx2; wny1=wny2

      if(nscreen.eq.1) write(*,101)'Time(day)=',time/(24.*60.*60.)
!      if(nscreen.eq.1) write(*,*)'Time=',time
      write(16,101)'Time(day)=',time/(24.*60.*60.)
101   format(A11,f5.2)
!--------------------------------------------------------------------------
      enddo !it
       
      close(60)
      close(61)
      close(62)
      close(63)
      close(64)  
      if(nscreen.eq.1) write(*,*)'Completed'
      write(16,*)'Completed'
      stop
      end


!******************************************************************************
!                                                                             *
!    Transform from lon,lat (lamda,phi) coordinates into CPP coordinates.     *
!    Lon,Lat must be in radians.                                              *
!                                                                             *
!******************************************************************************

      subroutine cpp(x,y,rlambda,phi,rlambda0,phi0)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)

      r=6378206.4
      x=r*(rlambda-rlambda0)*cos(phi0)
      y=phi*r

      return
      end
! .........................................................
      subroutine cppinverse(rlambda,phi,x,y,rlambda0,phi0)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)

      r=6378206.4
      rlambda=x / (r * cos(phi0)) + rlambda0
      phi=y/r

      return
      end
!
!********************************************************************************
!										*
!     Program to detect if two segments (1,2) and (3,4) have common pts   	*
!     Assumption: the 4 pts are distinctive.					*
!     The eqs. for the 2 lines are: X=X1+(X2-X1)*tt1 and X=X3+(X4-X3)*tt2.	*
!     Output: iflag: 0: no intersection or colinear; 1: exactly 1 intersection.	*
!     If iflag=1, (xin,yin) is the intersection.				*
!										*
!********************************************************************************
!
      subroutine intersect2(x1,x2,x3,x4,y1,y2,y3,y4,iflag,xin,yin,tt1,tt2)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)
      real(kind=dbl_kind1), parameter :: zero1=0.0 !small positive number or 0
      real(kind=dbl_kind1), intent(in) :: x1,x2,x3,x4,y1,y2,y3,y4
      integer, intent(out) :: iflag
      real(kind=dbl_kind1), intent(out) :: xin,yin,tt1,tt2

      tt1=-1000
      tt2=-1000
      iflag=0
      delta=(x2-x1)*(y3-y4)-(y2-y1)*(x3-x4)
      delta1=(x3-x1)*(y3-y4)-(y3-y1)*(x3-x4)
      delta2=(x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

      if(delta.ne.0.0d0) then
        tt1=delta1/delta
        tt2=delta2/delta
        if(tt1.ge.-zero1.and.tt1.le.1+zero1.and.tt2.ge.-zero1.and.tt2.le.1+zero1) then
          iflag=1
          xin=x1+(x2-x1)*tt1
          yin=y1+(y2-y1)*tt1
        endif
      endif

      return
      end

      function signa(x1,x2,x3,y1,y2,y3)
      use kind_par
      implicit real(kind=dbl_kind1)(a-h,o-z),integer(i-n)
      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2      
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                       !
!       Compute area coordinates of pt (xt,yt), which must be inside element nnel.      !
!       Impose bounds for area coordinates.                                             !
!                                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine area_coord(nnel,xt,yt,arco)                                             
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      integer, intent(in) :: nnel
      real(kind=dbl_kind), intent(in) :: xt,yt
      real(kind=dbl_kind), intent(out) :: arco(3)

      n1=elnode(1,nnel)
      n2=elnode(2,nnel)
      n3=elnode(3,nnel)
      arco(1)=signa(xt,x(n2),x(n3),yt,y(n2),y(n3))/area(nnel)
      arco(2)=signa(x(n1),xt,x(n3),y(n1),yt,y(n3))/area(nnel)
      arco(1)=dmax1(0.0d0,dmin1(1.0d0,arco(1)))
      arco(2)=dmax1(0.0d0,dmin1(1.0d0,arco(2)))
      if(arco(1)+arco(2)>1) then
        arco(3)=0
        arco(1)=dmin1(1.d0,dmax1(0.d0,arco(1)))
        arco(2)=1-arco(1)
      else
        arco(3)=1-arco(1)-arco(2)
      endif

      return
      end

!
!********************************************************************
!								    *
!	Routine to update z-coordinates and wetting and drying      *
!								    *
!********************************************************************
!
      subroutine levels
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      dimension idry_new(mnp) !,out2(12+mnv)

!...  z-coor. for nodes
!...  
      do i=1,np        
        if(dp(i)+eta3(i)<=h0) then !dry
          idry_new(i)=1
          if(dp(i)>=h_s) then
            write(11,*)'Deep depth dry:',i
            stop
          endif
          kbp(i)=0
        else !wet
          idry_new(i)=0

!         S-levels
          do k=kz,nvrt
            kin=k-kz+1
            if(hmod(i)<=h_c) then
              if(ifort12(12)==0) then
                ifort12(12)=1
                write(12,*)'Initial depth too shallow for S:',i,hmod(i),h_c
              endif
              iback(i)=1
              z(i,k)=sigma(kin)*(hmod(i)+eta3(i))+eta3(i)
            else if(eta3(i)<=-h_c-(hmod(i)-h_c)*theta_f/s_con1) then !hmod(i)>h_c>=0
              write(11,*)'Pls choose a larger h_c (1):',eta3(i),h_c
              stop
            else
              z(i,k)=eta3(i)*(1+sigma(kin))+h_c*sigma(kin)+(hmod(i)-h_c)*cs(kin)
            endif
          enddo !k=kz,nvrt

!         z-levels
          if(dp(i)<=h_s) then
            kbp(i)=kz
          else !bottom index will never change
            if(kbp(i)>=kz.or.kbp(i)<1) then
              write(11,*)'Impossible 92:',kbp(i),kz,i
              stop
            endif
            z(i,kbp(i))=-dp(i)       !!! added to ptrack.f90 !!!
            do k=kbp(i)+1,kz-1
              z(i,k)=ztot(k)
            enddo !k
          endif      

          do k=kbp(i)+1,nvrt
            if(z(i,k)-z(i,k-1)<=0) then
              write(11,*)'Inverted z-levels at:',i,k,z(i,k)-z(i,k-1),eta3(i),hmod(i)
              stop
            endif
          enddo !k
        endif !wet ot dry
      enddo !i=1,np

!...  Set wet/dry flags for element; element is "dry" if one of nodes is dry; conversely, 
!...  an element is wet if all nodes are wet (and all sides are wet as well)
!...  Weed out fake we nodes; a node is wet if and only if at least one surrounding element is wet
!...
!      idry_e0=idry_e !save
      idry=1 !dry unless wet
      kbe=0
      do i=1,ne
        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        idry_e(i)=max0(idry_new(n1),idry_new(n2),idry_new(n3))
        if(idry_e(i)==0) then
          idry(n1)=0; idry(n2)=0; idry(n3)=0
          kbe(i)=max0(kbp(n1),kbp(n2),kbp(n3))
        endif
      enddo !i

      return
      end

!
!********************************************************************************
!										*
!     Straightline search algorithm. Initially nnel0 is an element that 	*
!     encompasses (x0,y0). iloc=0: do not nudge initial pt; iloc=1: nudge.	* 
!     Input: iloc,nnel0,x0,y0,z0,xt,yt,zt,jlev0, time, and uu2,vv2,ww2 for 	*
!	abnormal cases;								*
!     Output the updated end pt (xt,yt,zt) (if so), nnel1, jlev1, area          *
!       coordinates, vertical ratio, a flag nfl, and local elevation and depth. *
!     nfl=1 if a bnd or dry element is hit and vel. there is small,		* 
!	or death trap is reached.						*
!										*
!********************************************************************************
!
      subroutine quicksearch(iloc,idt,ipar,nnel0,jlev0,time,x0,y0,z0,xt,yt,zt,nnel1,jlev1, &
     &arco,zrat,ss,nfl,etal,dp_p,ztmp,kbpl,ist,inp,str,pbeach)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      integer, intent(in) :: iloc,idt,ipar,nnel0,jlev0
      real(kind=dbl_kind), intent(in) :: time,x0,y0,z0,str
      integer, intent(out) :: nnel1,jlev1,nfl,kbpl
      integer, intent(inout) :: ist(10*mnp),inp(10*mnp)
      real(kind=dbl_kind), intent(inout) :: xt,yt,zt
      real(kind=dbl_kind), intent(out) :: arco(3),zrat,ss,etal,dp_p,ztmp(mnv)
!      dimension out2(mnv)
! ...............................................................................
      if(iloc>1) then
        write(11,*)'iloc > 1'
        stop
      endif
      if(idry_e(nnel0)==1) then
        write(11,*)'Starting element is dry'
        stop
      endif

      nfl=0
      trm=time !time remaining

!     Starting element nel
!     Try area_coord
      nel=nnel0
      aa=0
      aa1=0
      do i=1,3
        n1=elnode(i,nel)
        n2=elnode(nx(i,1),nel)
        aa=aa+dabs(signa(x(n1),x(n2),x0,y(n1),y(n2),y0))
        aa1=aa1+dabs(signa(x(n1),x(n2),xt,y(n1),y(n2),yt))
      enddo !i
      ae=dabs(aa-area(nel))/area(nel)
      if(ae>small1) then
        write(11,*)'(x0,y0) not in nnel0 initially',ae,nnel0
        stop
      endif

      ae=dabs(aa1-area(nel))/area(nel)
      if(ae<small1) then
        nnel1=nel
        go to 400
      endif

!     (xt,yt) not in nel, and thus (x0,y0) and (xt,yt) are distinctive
!     An interior pt close to (x0,y0) to prevent underflow for iloc >=1.
      if(iloc==0) then
        xcg=x0
        ycg=y0
      else if(iloc==1) then
!        weit=1./3; al=0; bet=0
!        xint=x(elnode(1,nel))*(weit-al)+x(elnode(2,nel))*(weit-bet)+x(elnode(3,nel))*(weit+al+bet)
!        yint=y(elnode(1,nel))*(weit-al)+y(elnode(2,nel))*(weit-bet)+y(elnode(3,nel))*(weit+al+bet)
!        xcg=(1-5.e-4)*x0+5.e-4*xint
!        ycg=(1-5.e-4)*y0+5.e-4*yint
        xcg=(1-1.0d-4)*x0+1.0d-4*xctr(nel)
        ycg=(1-1.0d-4)*y0+1.0d-4*yctr(nel)
      endif

      pathl=dsqrt((xt-xcg)**2+(yt-ycg)**2)
      if(pathl==0) then
        write(11,*)'Zero path',x0,y0,xt,yt,xcg,ycg
        stop
      endif

!     Starting edge nel_j
      nel_j=0
      do j=1,3
        jd1=elnode(nx(j,1),nel)
        jd2=elnode(nx(j,2),nel)
        call intersect2(xcg,xt,x(jd1),x(jd2),ycg,yt,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
        if(iflag==1) then
          nel_j=j
          exit
        endif
      enddo !j=1,3
      if(nel_j==0) then
        write(11,*)'Found no intersecting edges I:',nel,xcg,ycg,xt,yt,ae
        stop
      endif

      zin=z0 !intialize
      it=0
      loop4: do
!----------------------------------------------------------------------------------------
      it=it+1
      if(it>1000) then
        if(ifort12(3)==0) then
          ifort12(3)=1
          write(12,*)'Death trap reached',idt,id0
        endif
        nfl=1
        xt=xin
        yt=yin
        zt=zin
        nnel1=nel
        exit loop4
      endif
      md1=elnode(nx(nel_j,1),nel)
      md2=elnode(nx(nel_j,2),nel)
      
!     Compute z position 
      dist=dsqrt((xin-xt)**2+(yin-yt)**2)
      if(dist/pathl.gt.1+1.0d-4) then
        write(11,*)'Path overshot'
        stop
      endif
      zin=zt-dist/pathl*(zt-zin)
      trm=trm*dist/pathl !time remaining
      pathl=dsqrt((xin-xt)**2+(yin-yt)**2)
      if(pathl==0.or.trm==0) then
        write(11,*)'Target reached'
        stop
      endif

      lit=0 !flag
!     For horizontal exit and dry elements, compute tangential vel.,
!     update target (xt,yt,zt) and continue.
      if(ic3(nel_j,nel)==0.or.idry_e(ic3(nel_j,nel))==1) then
        lit=1
        nfl=1
        isd=elside(nel_j,nel)
        nd1=isidenode(1,isd)
        nd2=isidenode(2,isd)       
        if(nd1+nd2/=md1+md2) then
          write(11,*)'Wrong side'
          stop
        endif

!... start boundary treatment
!... particle hits the land
        if(idry_e(ic3(nel_j,nel))==1) then  ! dry element
!    let the pbeach percent of hitted particles still moved
          if(str.ge.0.0 .and. str.lt.pbeach) then
            go to 105
           else
            ist(ipar)=-1  ! stranding to shore
            inp(ipar)=ic3(nel_j,nel)
            nfl=1 
!            xt=(1-1.0d-4)*xin+1.0d-4*xctr(nel)
!            yt=(1-1.0d-4)*yin+1.0d-4*yctr(nel)        
            xt=xin
            yt=yin
            zt=zin
            nnel1=nel
            exit loop4
          endif
        endif

!... particle hits the boundary
!    for the land boundary
        if(isbnd(nd1).eq.-1 .and. isbnd(nd2).eq.-1) then
          if(str.ge.0.0 .and. str.lt.0.2) then
            go to 105
           else
            ist(ipar)=-1  ! stranding to shore permanently
            nfl=1
!            xt=(1-1.0d-4)*xin+1.0d-4*xctr(nel)
!            yt=(1-1.0d-4)*yin+1.0d-4*yctr(nel)        
            xt=xin
            yt=yin
            zt=zin
            nnel1=nel
            exit loop4
          endif
        endif

!    for the river boundary 
        if(isbnd(nd1).eq.1 .or. isbnd(nd2).eq.1) then
          ist(ipar)=2  
          nfl=1
          xt=xin
          yt=yin
          zt=zin
          nnel1=nel
          exit loop4
        endif
!    for the ocean boundary 
        if(isbnd(nd1).eq.2 .or. isbnd(nd2).eq.2) then
          ist(ipar)=3  
          nfl=1
          xt=xin
          yt=yin
          zt=zin
          nnel1=nel
          exit loop4
        endif
!... finish boundary treatment

105     continue       
!       Nudge intersect (xin,yin), and update starting pt
        xin=(1-1.0d-4)*xin+1.0d-4*xctr(nel)
        yin=(1-1.0d-4)*yin+1.0d-4*yctr(nel)
        xcg=xin
        ycg=yin
        vtan=-(uu2(md1,jlev0)+uu2(md2,jlev0))/2*sny(isd)+(vv2(md1,jlev0)+vv2(md2,jlev0))/2*snx(isd)
        xvel=-vtan*sny(isd)
        yvel=vtan*snx(isd)
        zvel=(ww2(md1,jlev0)+ww2(md2,jlev0))/2
        xt=xin-xvel*trm
        yt=yin-yvel*trm
        zt=zin-zvel*trm
        hvel=dsqrt(xvel**2+yvel**2)
        if(hvel<1.e-4) then
          nfl=1
          xt=xin
          yt=yin
          zt=zin
          nnel1=nel
          exit loop4
        endif
        pathl=hvel*trm
      endif !abnormal cases

!     Search for nel's neighbor with edge nel_j, or in abnormal cases, the same element
      if(lit==0) nel=ic3(nel_j,nel) !next front element
      aa=0
      do i=1,3
        k1=elnode(i,nel)
        k2=elnode(nx(i,1),nel)
        aa=aa+dabs(signa(x(k1),x(k2),xt,y(k1),y(k2),yt))
      enddo !i
      ae=dabs(aa-area(nel))/area(nel)
      if(ae<small1) then
        nnel1=nel
        exit loop4
      endif

!     Next intersecting edge
      do j=1,3
         jd1=elnode(nx(j,1),nel)
         jd2=elnode(nx(j,2),nel)
!        For abnormal case, same side (border side) cannot be hit again
         if(jd1==md1.and.jd2==md2.or.jd2==md1.and.jd1==md2) cycle
         call intersect2(xcg,xt,x(jd1),x(jd2),ycg,yt,y(jd1),y(jd2),iflag,xin,yin,tt1,tt2)
         if(iflag==1) then
           nel_j=j !next front edge
           cycle loop4
         endif
      enddo !j
      write(11,*)'Failed to find next edge I:',lit,xin,yin,xt,yt,nel,md1,md2,idt,id0,ae
      stop
!----------------------------------------------------------------------------------------
      end do loop4  

400   continue
!     No vertical exit from domain
      if(idry_e(nnel1)==1) then
        write(11,*)'Ending element is dry'
        stop
      endif

!     Compute area & sigma coord.
      call area_coord(nnel1,xt,yt,arco)
      n1=elnode(1,nnel1)
      n2=elnode(2,nnel1)
      n3=elnode(3,nnel1)
      etal=eta3(n1)*arco(1)+eta3(n2)*arco(2)+eta3(n3)*arco(3)
      dep=dp(n1)*arco(1)+dp(n2)*arco(2)+dp(n3)*arco(3)
      dp_p=dep
      if(etal+dep<h0) then
        write(11,*)'Weird wet element in quicksearch:',nnel1,eta3(n1),eta3(n2),eta3(n3)
        stop
      endif

      if(istiff==1) zt=etal+zpar0(ipar)

!     Compute z-levels
      do k=kz,nvrt
        kin=k-kz+1
        hmod2=dmin1(dep,h_s)
        if(hmod2<=h_c) then
          ztmp(k)=sigma(kin)*(hmod2+etal)+etal
        else if(etal<=-h_c-(dep-h_c)*theta_f/s_con1) then
          write(11,*)'Pls choose a larger h_c (2):',etal,h_c
          stop
        else
          ztmp(k)=etal*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
        endif

!       Following to prevent underflow
        if(k==kz) ztmp(k)=-hmod2
        if(k==nvrt) ztmp(k)=etal
      enddo !k

      if(dep<=h_s) then
        kbpl=kz
      else !z levels
!       Find bottom index
        kbpl=0
        do k=1,kz-1
          if(-dep>=ztot(k).and.-dep<ztot(k+1)) then
            kbpl=k
            exit
          endif
        enddo !k
        if(kbpl==0) then
          write(11,*)'Cannot find a bottom level at foot:',dep
          stop
        endif
        ztmp(kbpl)=-dep
        do k=kbpl+1,kz-1
          ztmp(k)=ztot(k)
        enddo !k
      endif

!        print*, etal,zt,zpar0(ipar)
!       do k=1,nvrt
!        print*, i,k,ztmp(k)
!       enddo
!       pause

      do k=kbpl+1,nvrt
        if(ztmp(k)-ztmp(k-1)<=0) then
          write(11,*)'Inverted z-level in quicksearch:',nnel1,etal,dep,ztmp(k)-ztmp(k-1)
          stop
        endif
      enddo !k
 
      if(zt<=ztmp(kbpl)) then
        zt=ztmp(kbpl)
        zrat=1
        jlev1=kbpl+1
      else if(zt>=ztmp(nvrt)) then
        zt=ztmp(nvrt)
        zrat=0
        jlev1=nvrt
      else
        jlev1=0
        do k=kbpl,nvrt-1
          if(zt>=ztmp(k).and.zt<=ztmp(k+1)) then 
            jlev1=k+1
            exit
          endif
        enddo !k
        if(jlev1==0) then
          write(11,*)'Cannot find a vert. level:',zt,etal,dep
          write(11,*)(ztmp(k),k=kbpl,nvrt)
          stop
        endif
        zrat=(ztmp(jlev1)-zt)/(ztmp(jlev1)-ztmp(jlev1-1))
      endif

      if(zrat<0.or.zrat>1) then
        write(11,*)'Sigma coord. wrong (4):',jlev1,zrat
        stop
      endif

      if(kbpl==kz) then !in pure S region
        ss=(1-zrat)*sigma(jlev1-kz+1)+zrat*sigma(jlev1-kz)
      else
        ss=-99
      endif

500   continue
      return
      end

      subroutine interp(i,nnel,jlev,ztmp,trat,zrat,arco,upar,vpar,wpar,dhfx,dhfy,dhfz,grdx,&
           & grdy,grdz,wndx,wndy)
!     data interpolation for tracking
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      integer, intent(in) :: i,nnel,jlev
      real(kind=dbl_kind), intent(in) :: trat,zrat,arco(3)
      real(kind=dbl_kind), intent(out) :: upar(10*mnp),vpar(10*mnp),wpar(10*mnp)
      real(kind=dbl_kind), intent(out) :: wndx(10*mnp),wndy(10*mnp)
      real(kind=dbl_kind), intent(out) :: dhfx(10*mnp),dhfy(10*mnp),dhfz(10*mnp)
      real(kind=dbl_kind), intent(out) :: grdx(10*mnp),grdy(10*mnp),grdz(10*mnp)
      real(kind=dbl_kind) :: as(4),bs(4),xcp(4),ycp(4),ztmp(mnv)
      real(kind=dbl_kind) :: vwx(4),vwy(4),vxl(4,2),vyl(4,2),vzl(4,2),val(4,2), &
           & vbl(4,2),vcl(4,2),vdl(4,2),vxn(4),vyn(4),vzn(4),van(4),vcn(4),vdn(4)
! .................................................................................
!     Interpolate in time
      do j=1,3
        nd=elnode(j,nnel)
        vwx(j)=wnx1(nd)*(1-trat)+wnx2(nd)*trat  ! windx 
        vwy(j)=wny1(nd)*(1-trat)+wny2(nd)*trat  ! windy
        do l=1,2
          lev=jlev+l-2
          vxl(j,l)=uu1(nd,lev)*(1-trat)+uu2(nd,lev)*trat
          vyl(j,l)=vv1(nd,lev)*(1-trat)+vv2(nd,lev)*trat
          vzl(j,l)=ww1(nd,lev)*(1-trat)+ww2(nd,lev)*trat
          val(j,l)=hf1(nd,lev)*(1-trat)+hf2(nd,lev)*trat      ! grdx
          vbl(j,l)=vf1(nd,lev)*(1-trat)+vf2(nd,lev)*trat      ! grdz
          vcl(j,l)=hf1(nd,lev)*(1-trat/2)+hf2(nd,lev)*trat/2  ! Dh for half-step
          vdl(j,l)=vf1(nd,lev)*(1-trat/2)+vf2(nd,lev)*trat/2  ! Dz for half-step
        enddo !l
      enddo !j
     
!     Interpolate in vertical 
      do j=1,3
        vxn(j)=vxl(j,2)*(1-zrat)+vxl(j,1)*zrat
        vyn(j)=vyl(j,2)*(1-zrat)+vyl(j,1)*zrat
        vzn(j)=vzl(j,2)*(1-zrat)+vzl(j,1)*zrat
        van(j)=val(j,2)*(1-zrat)+val(j,1)*zrat  ! grdx
        vcn(j)=vcl(j,2)*(1-zrat)+vcl(j,1)*zrat  ! Dh for half-step
        vdn(j)=vdl(j,2)*(1-zrat)+vdl(j,1)*zrat  ! Dz for half-step
      enddo !j

!     Interpolate in horizontal
        wndx(i)=0;wndy(i)=0;upar(i)=0;vpar(i)=0
        wpar(i)=0;dhfx(i)=0;dhfz(i)=0   
      do j=1,3
        wndx(i)=wndx(i)+vwx(j)*arco(j)  ! windx
        wndy(i)=wndy(i)+vwy(j)*arco(j)  ! windy
        upar(i)=upar(i)+vxn(j)*arco(j)
        vpar(i)=vpar(i)+vyn(j)*arco(j)
        wpar(i)=wpar(i)+vzn(j)*arco(j)
        dhfx(i)=dhfx(i)+vcn(j)*arco(j)  ! Dh for half-step
        dhfz(i)=dhfz(i)+vdn(j)*arco(j)  ! Dz for half-step
      enddo !j
        dhfy(i)=dhfx(i)

! ..............................................................
!     Compute the vertical gradient of Dz
      az=0 
      do j=1,3
        az=az+(vbl(j,2)-vbl(j,1))*arco(j)
      enddo !j
      lev=jlev
      grdz(i)=az/(ztmp(lev)-ztmp(lev-1)) 
!     print*, az,lev,ztmp(lev),ztmp(lev-1),grdz(i)
!     pause
! ..............................................................
!     Compute the horizontal gradient of Dh
      do j=1,3
         xcp(j)=x(elnode(j,i))
         ycp(j)=y(elnode(j,i))               
      enddo
      ax=0;ay=0   
      do j=1,3
         m1=j+1;m2=j+2
         if(m1>3) m1=m1-3
         if(m2>3) m2=m2-3            
           ax=ax+van(j)*(ycp(m1)-ycp(m2))
           ay=ay+van(j)*(xcp(m2)-xcp(m1))
      enddo                  
      grdx(i)=1./(2.*area(i))*ax
      grdy(i)=1./(2.*area(i))*ay 
! ..............................................................
      return
      end
!
! 
      subroutine evapor(ddt,time,stime,pmass)
      use global
      implicit real(kind=dbl_kind)(a-h,o-z),integer(i-n)
      integer ts,elt
      real(kind=dbl_kind) :: ddt,tt,t0,a,y0,yc,rat,stime,pmass
      real, dimension(10*mnp) :: ym
      logical first
      data first/.true./
      save first,ym
! ..............................................................
! ... compute the decay rate
      if(first) then    
        tt=rnday*24*3600      ! total time
        t0=1.0*24*3600        ! half-life
        a=-log(2.0)/t0        ! decay rate
        ts=int(tt/ddt)        ! total time step
!
        y0=100.0              ! initial mass
        rat=0.6               ! remain ratio
        yc=y0*rat             ! last mass 
        ym(1)=y0

        do k=1,ts+1
           ym(k+1)=(y0-yc)*exp(a*k*ddt)+yc       
        enddo
        first=.false.
!        do i=1,ts,50 
!          write(*,*) i,ym(i)
!        enddo
!        pause
      endif   

! ... remained particle mass after evaporation
      elt=nint((time-stime)/ddt)
      pmass=ym(elt+1)/100.
      return
      end

! ****************************************************      
! RANDOM NUMBER GENERATORS from numerical recipes       
! ****************************************************
  FUNCTION ran1(idum)
   Integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,&
      & IR2,NTAB,NDIV
   Real ran1,AM,EPS,RNMX
   Parameter(IM1=2147483563,IM2=2147483399,AM=1./&
      & IM1,IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,&
      & IQ2=52774,IR1=12211,IR2=3791,NTAB=32,NDIV=1+&
      & IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   Integer idum2,j,k,iv(NTAB),iy
   SAVE iv,iy,idum2
   DATA idum2/123456789/,iv/NTAB*0/,iy/0/
! ....................................................
   if(idum.le.0) then
     idum=max(-idum,1)
     idum2=idum
     do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if(idum.lt.0) idum=idum+IM1
        if(j.le.NTAB) iv(j)=idum
     enddo
     iy=iv(1) 
   endif
   k=idum/IQ1
   idum=IA1*(idum-k*IQ1)-k*IR1
   if(idum.lt.0) idum=idum+IM1
     k=idum2/IQ2
     idum2=IA2*(idum2-k*IQ2)-k*IR2
   if(idum2.lt.0) idum2=idum2+IM2
     j=1+iy/NDIV
     iy=iv(j)-idum2
     iv(j)=idum
   if(iy.lt.1) iy=iy+IMM1
     ran1=min(AM*iy,RNMX)  
   return
  end
      
  FUNCTION ran2(idum)
   Integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,&
      & IR2,NTAB,NDIV
   Real ran2,AM,EPS,RNMX
   Parameter(IM1=2147483563,IM2=2147483399,AM=1./&
      & IM1,IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,&
      & IQ2=52774,IR1=12211,IR2=3791,NTAB=32,NDIV=1+&
      & IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   Integer idum2,j,k,iv(NTAB),iy
   SAVE iv,iy,idum2
   DATA idum2/123456789/,iv/NTAB*0/,iy/0/
! ....................................................
   if(idum.le.0) then
     idum=max(-idum,1)
     idum2=idum
     do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if(idum.lt.0) idum=idum+IM1
        if(j.le.NTAB) iv(j)=idum
     enddo
     iy=iv(1) 
   endif
   k=idum/IQ1
   idum=IA1*(idum-k*IQ1)-k*IR1
   if(idum.lt.0) idum=idum+IM1
     k=idum2/IQ2
     idum2=IA2*(idum2-k*IQ2)-k*IR2
   if(idum2.lt.0) idum2=idum2+IM2
     j=1+iy/NDIV
     iy=iv(j)-idum2
     iv(j)=idum
   if(iy.lt.1) iy=iy+IMM1
     ran2=min(AM*iy,RNMX)  
   return
  end

  FUNCTION ran3(idum)
   Integer idum
   Integer MBIG,MSEED,MZ
   Real ran3,FAC 
   Parameter(MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
   Integer i,iff,ii,inext,inextp,k
   Integer mj,mk,ma(55) 
   SAVE iff,inext,inextp,ma
   DATA iff /0/
! ........................................................... 
   if(idum.lt.0.or.iff.eq.0) then
     iff=1
     mj=abs(MSEED-abs(idum))
     mj=mod(mj,MBIG)
     ma(55)=mj
     mk=1
     do i=1,54
        ii=mod(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        if(mk.lt.MZ) mk=mk+MBIG
          mj=ma(ii)
     enddo   
     do k=1,4
       do i=1,55
          ma(i)=ma(i)-ma(1+mod(i+30,55))
          if(ma(i).lt.MZ) ma(i)=ma(i)+MBIG
       enddo
     enddo
     inext=0
     inextp=31
     idum=1
   endif  
     inext=inext+1
     if(inext.eq.56) inext=1
       inextp=inextp+1
     if(inextp.eq.56) inextp=1
       mj=ma(inext)-ma(inextp)
     if(mj.lt.MZ) mj=mj+MBIG
       ma(inext)=mj
       ran3=mj*FAC
   return
  end
