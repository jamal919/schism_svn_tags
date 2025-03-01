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

!
!****************************************************************************************
!											*
!	read_output7b with (x,y) read in from station.bp (build pts), and output
!       results along a transect (defined below) for 3D variables. 
!       Need to manually modify ntran etc. If filtering is used, it'd be done on hourly
!       series.

!       Inputs: screen; vgrid (in this dir or ../); station.bp (build pts); binary
!       Outputs: residual.out (*.6[12] type) if ifilter=1; otherwise 
!       transect.out & transect_grd.[zr]0 (ascii; use plot_transect.m)

!       ifort -Bstatic -assume byterecl -O3 -o read_output7b_transect read_output7b_transect.f90 moving_average_filter.f90 ../../UtilLib/compute_zcor.f90 
!       pgf90 -O2 -mcmodel=medium -Mbounds -Bstatic -o read_output7b_transect read_output7b_transect.f90 moving_average_filter.f90 ../../UtilLib/compute_zcor.f90
!											*
!****************************************************************************************
!
      program read_output_transect
      character*30  :: file    !< base file name (e.g. elev.61) of files to be read 
      character*80  :: stationfile  !< input file, e.g. station.bp
      character*80  :: vgrid  !< vertical grid file path, e.g. vgrid.in
      character*80  :: outfile
      integer       :: iday1
      integer       :: iday2
      real          :: dt_output
      integer       :: ifilter
      call read_output7b_transect_input(file,vgrid,stationfile,outfile,iday1,iday2,dt_output,ifilter)
      call read_output7b_transect(file,vgrid,stationfile,outfile,iday1,iday2,dt_output,ifilter)
      end program


      subroutine read_output7b_transect_input(file,vgrid,stationfile,outfile,iday1,iday2,dt_output,ifilter)
      use argparse
      implicit none
      character(len=*)  :: file    !< base file name (e.g. elev.61) of files to be read 
      character(len=*)  :: outfile    !< base file name (e.g. elev.61) of files to be read 
      character(len=80)  :: infile  !< input file containing options for this routine
      character(len=*)  :: stationfile  !< input file, e.g. station.xyzt 
      character(len=*)  :: vgrid        !< vgrid file
      integer       :: iday1
      integer       :: iday2
      real          :: dt_output
      integer       :: ifilter
      integer       :: comcount
      logical       :: lfilter
      cmd_name = "read_output7b_transect"
      call cla_init(cmd_name)
      call cla_register('-f','--file', 'base input file (e.g. salt.63) without block number prefex', cla_char,'')
      call cla_register('-v','--vgrid', 'path to vertical grid file (typically vgrid.in or ../vgrid.in', cla_char,'vgrid.in')
      call cla_register('-i','--in', 'input file containing options (possibly overriden at command line)',cla_char,'')
      call cla_register('-s','--stations','input file listing stations (first line is #, then x y t)',cla_char,'station.bp')
      call cla_register('-b','--bgn', 'first block to read', cla_int,'1')
      call cla_register('-e','--end', 'last block to read', cla_float  ,'-1')
      call cla_register('-d','--dt','output dt)', cla_float  ,'-1')
      call cla_register('-l','--filter','filter time series', cla_flag, 'f')

      call cla_validate

      comcount = command_argument_count()
      if (comcount == 0) then
          infile = "read_output7b_transect.in"
      else
          call cla_get("--in",infile)
      end if

      if (len_trim(infile)>1) then
      open(10,file=trim(infile),status='old')
!      print*, 'Input file to read from (without *_):'
      read(10,'(a30)')file
!     Fill number used for invalid 3D variables: below bottom; dry spot; no parents
!      print*, 'Input values to be used for invalid place:'
      read(10,*) iday1, iday2
!      print*, 'Input starting corie day:'
      read(10,*)dt_output
      read(10,*)ifilter
      close(10)
      ! override with (explicitly provided) command line values)
      if (cla_key_present("--file")) call cla_get("--file",file)
      if (cla_key_present("--bgn"))  call cla_get("--bgn",iday1)
      if (cla_key_present("--end"))  call cla_get("--end",iday2)
      if (cla_key_present("--dt"))   call cla_get("--dt", dt_output)
      if (cla_key_present("--filter")) call cla_get_flag("--filter", lfilter)
      
      else
      call cla_get("--file",file)
      call cla_get("--bgn",iday1)
      call cla_get("--end",iday2)
      call cla_get("--dt", dt_output)
      call cla_get_flag("--filter", lfilter)

      end if
      if (lfilter)then
          ifilter = 1
      else 
          ifilter = 0
      end if
      call cla_get("--stations",stationfile)
      call cla_get("--vgrid",vgrid)
      end subroutine
      


      subroutine read_output7b_transect(file63,vgrid,stationfile,outfile,iday1,iday2,dt_output,ifilter)
      parameter(nbyte=4)
      character(len=*) :: file63
      character(len=*) :: stationfile
      character(len=*) :: outfile
      character(len=*) :: vgrid
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      integer,allocatable :: elnode(:,:)
      allocatable :: sigma(:),cs(:),ztot(:),x(:),y(:),dp(:),kbp00(:),kfp(:)
      allocatable:: out(:,:,:,:),out2(:,:,:),icum(:,:,:),eta2(:),node3(:,:),arco(:,:)
      allocatable :: ztmp(:),x00(:),y00(:),iep(:),out3(:,:,:),out4(:,:)
      allocatable :: ser(:),ser2(:),out5(:,:,:),nmxz(:,:),z0(:),r0(:),r00(:),z00(:)
      allocatable :: sigma_lcl(:,:),kbp(:),ztmp2(:,:)
      dimension swild(3)
      
 !     Input transect depths
      ntran=375
      allocate(z0(ntran))
      z0(1:ntran)=(/(-374+i, i=0,ntran-1) /)

      open(10,file=trim(stationfile),status='old')
      read(10,*) 
      read(10,*) nxy
      nxz=nxy*ntran
      nxz_e=(nxy-1)*(ntran-1)*2
      allocate(x00(nxy),y00(nxy),r0(nxy),r00(nxz),z00(nxz),nmxz(nxz_e,4),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do i=1,nxy
        read(10,*)j,x00(i),y00(i)
      enddo !i
      !Along transect distance
      r0(1)=0
      do i=2,nxy
        r0(i)=r0(i-1)+sqrt((x00(i)-x00(i-1))**2+(y00(i)-y00(i-1))**2)
      enddo !i
      close(10)
 
!     Output transect grid
      if(ifilter==0) then
        open(22,file='transect_grd.z0',status='replace')
        open(23,file='transect_grd.r0',status='replace')
        write(22,'(f12.3)')z0
        write(23,'(e15.6)')r0
        close(22) 
        close(23)
      endif !ifilter

!     Compute connectivity
      do i=1,nxy
        do j=1,ntran
          nd=ntran*(i-1)+j
          r00(nd)=r0(i)
          z00(nd)=z0(j)
        enddo !j
      enddo !i

      do i=1,nxy-1
        do j=1,ntran-1
          ie=(ntran-1)*(i-1)+j
          n1=ntran*(i-1)+j
          n2=ntran*i+j
          nmxz(2*ie-1,1)=n1
          nmxz(2*ie-1,2)=n2
          nmxz(2*ie-1,3)=n2+1
          nmxz(2*ie,1)=n1
          nmxz(2*ie,2)=n2+1
          nmxz(2*ie,3)=n1+1
        enddo !j
      enddo !i

      open(21,file='transect.out',status='replace')
      open(65,file='residual.out',access='direct',recl=nbyte)
      
!...  Header
!...
      open(63,file='1_'//file63,status='old',access='direct',recl=nbyte)
      irec=0
      do m=1,48/nbyte
        read(63,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m)
      enddo
      if(data_format.ne.'DataFormat v5.0') then
        print*, 'This code reads only v5.0:  ',data_format
        stop
      endif
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) version(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) start_time(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte

      write(*,'(a48)')data_format
      write(*,'(a48)')version
      write(*,'(a48)')start_time
      write(*,'(a48)')variable_nm
      write(*,'(a48)')variable_dim

      read(63,rec=irec+1) nrec
      read(63,rec=irec+2) dtout
      read(63,rec=irec+3) nspool
      read(63,rec=irec+4) ivs
      read(63,rec=irec+5) i23d

      iskip=dt_output/dtout
      if(iskip<1) then
        print*, 'Output interval too small; resetting to dtout=',dtout
        iskip=1
      endif
      nstep=(iday2-iday1+1)*nrec/iskip
!      if(nxz*ivs>6000) stop 'Increase output statement below!'
 
      write(65,rec=irec+1) nstep
      write(65,rec=irec+2) 3600.
      write(65,rec=irec+3) 1
      write(65,rec=irec+4) ivs
      write(65,rec=irec+5) 2
      irec=irec+5

      print*, 'i23d=',i23d,' nrec= ',nrec

!     Vertical grid (obsolete)
      read(63,rec=irec+1) nvrt
      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h0
      read(63,rec=irec+4) h_s
      read(63,rec=irec+5) h_c
      read(63,rec=irec+6) theta_b
      read(63,rec=irec+7) theta_f
      write(65,rec=irec+1) nvrt
      write(65,rec=irec+2) kz
      write(65,rec=irec+3) h0
      write(65,rec=irec+4) h_s
      write(65,rec=irec+5) h_c
      write(65,rec=irec+6) theta_b
      write(65,rec=irec+7) theta_f
      irec=irec+7
!      allocate(ztot(nvrt),sigma(nvrt),cs(nvrt),ztmp(nvrt),stat=istat)
!      if(istat/=0) stop 'Falied to allocate (2)'
!
      do k=1,kz-1
        read(63,rec=irec+k) tmp !ztot(k)
        write(65,rec=irec+k) tmp !ztot(k)
      enddo
      do k=kz,nvrt
        kin=k-kz+1
        read(63,rec=irec+k) tmp !sigma(kin)
        write(65,rec=irec+k) tmp !sigma(kin)
!        cs(kin)=(1-theta_b)*sinh(theta_f*sigma(kin))/sinh(theta_f)+ &
!     &theta_b*(tanh(theta_f*(sigma(kin)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
      enddo
      irec=irec+nvrt
      irec3=irec !for 65

!     Horizontal grid
      read(63,rec=irec+1) np
      read(63,rec=irec+2) ne
      irec=irec+2
      write(65,rec=irec3+1) nxz
      write(65,rec=irec3+2) nxz_e
      irec3=irec3+2
      allocate(x(np),y(np),dp(np),kbp00(np),kfp(np),elnode(4,ne),out(np,3,nvrt,2), &
     &out2(np,nvrt,2),icum(np,nvrt,2),eta2(np),node3(np,3),arco(np,3),iep(np), &
     &out3(nxz,2,nstep),out4(ntran,2),ser(nstep),ser2(nstep),out5(np,2,nstep),stat=istat)
      if(istat/=0) stop 'Falied to allocate (5)'

      do m=1,np
        read(63,rec=irec+1)x(m)
        read(63,rec=irec+2)y(m)
        read(63,rec=irec+3)dp(m)
        read(63,rec=irec+4)kbp00(m)
        irec=irec+4
      enddo !m=1,np
 
      do m=1,nxz
        write(65,rec=irec3+1)r00(m)
        write(65,rec=irec3+2)z00(m)
        write(65,rec=irec3+3)1.
        write(65,rec=irec3+4)1
        irec3=irec3+4
      enddo !m
      do m=1,ne
        read(63,rec=irec+1)i34
        irec=irec+1
        do mm=1,i34
          read(63,rec=irec+1)elnode(mm,m)
          irec=irec+1
        enddo !mm
      enddo !m
      do m=1,nxz_e
        write(65,rec=irec3+1) 3
        irec3=irec3+1
        do mm=1,3
          write(65,rec=irec3+1)nmxz(m,mm)
          irec3=irec3+1
        enddo !mm
      enddo !m
      irec0=irec

      print*, 'last element',(elnode(j,ne),j=1,3)

!     Read in vgrid.in
      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np))
      call get_vgrid(trim(vgrid),np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

      allocate(ztmp(nvrt),ztmp2(nvrt,3))

!...  Find parent element for (x00,y00)
      iep=0
      do i=1,ne
        do l=1,nxy
          aa=0
          ar=0 !area
          do j=1,3
            j1=j+1
            j2=j+2
            if(j1>3) j1=j1-3
            if(j2>3) j2=j2-3
            n0=elnode(j,i)
            n1=elnode(j1,i)
            n2=elnode(j2,i)
            swild(j)=signa(x(n1),x(n2),x00(l),y(n1),y(n2),y00(l)) !temporary storage
            aa=aa+abs(swild(j))
            if(j==1) ar=signa(x(n1),x(n2),x(n0),y(n1),y(n2),y(n0))
          enddo !j
          if(ar<=0) then
            print*, 'Negative area:',ar
            stop
          endif
          ae=abs(aa-ar)/ar
          if(ae<=1.e-5) then
            iep(l)=i
            node3(l,1:3)=elnode(1:3,i)
            arco(l,1:3)=swild(1:3)/ar
            arco(l,1)=max(0.,min(1.,arco(l,1)))
            arco(l,2)=max(0.,min(1.,arco(l,2)))
            if(arco(l,1)+arco(l,2)>1) then 
              arco(l,3)=0
              arco(l,2)=1-arco(l,1)
            else
              arco(l,3)=1-arco(l,1)-arco(l,2)
            endif
            cycle
          endif
        enddo !l; build pts

        ifl=0 !flag
        do l=1,nxy
          if(iep(l)==0) then
            ifl=1
            exit
          endif
        enddo !l
        if(ifl==0) exit
      enddo !i=1,ne

      do j=1,nxy
        if(iep(j)==0) then
          print*, 'Cannot find a parent for pt:',j,x00(j),y00(j)
          stop
        endif
      enddo !j

!...  Compute relative record # for a node and level for 3D outputs
!...
      icount=0
      do i=1,np
        do k=max0(1,kbp00(i)),nvrt
          do m=1,ivs
            icount=icount+1
            icum(i,k,m)=icount
          enddo !m
        enddo !k
      enddo !i=1,np

      if(file63(1:7).eq.'hvel.64') then
        rjunk=0 !invalid values
      else
        rjunk=-9999
      endif

!...  Time iteration
!...
      its=0 !time step counter
      do iday=iday1,iday2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      open(63,file=it_char//'_'//file63,status='old',access='direct',recl=nbyte)

      irec=irec0
      do it1=1,nrec
        its=its+1
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2

!        print*, 'time=',time/86400

        do i=1,np
          read(63,rec=irec+i) eta2(i)
        enddo !i
        irec=irec+np

        out2=0
        if(i23d.eq.2) then
          do i=1,nxy
            do j=1,3 !nodes
              nd=node3(i,j)
              do m=1,ivs
                read(63,rec=irec+(nd-1)*ivs+m) tmp
                out2(i,1,m)=out2(i,1,m)+arco(i,j)*tmp
              enddo !m
            enddo !j
          enddo !i
          irec=irec+np*ivs
!          write(65,'(e16.8,12(1x,f12.3))')time,(out2(i,1,1),i=1,nxy)
        else !i23d=3 
          do i=1,nxy
            do j=1,3 !nodes
              nd=node3(i,j)
              do k=max0(1,kbp00(nd)),nvrt
                do m=1,ivs
                  read(63,rec=irec+icum(nd,k,m)) out(i,j,k,m)
                enddo !m
              enddo !k
            enddo !j
          enddo !i
          irec=irec+icum(np,nvrt,ivs)

!         Do interpolation
          do i=1,nxy
            etal=0; dep=0; idry=0
            do j=1,3
              nd=node3(i,j)
              if(eta2(nd)+dp(nd)<h0) idry=1
              etal=etal+arco(i,j)*eta2(nd)
              dep=dep+arco(i,j)*dp(nd)
      
!             Debug
!              write(11,*)i,j,nd,dp(nd),arco(i,j)

            enddo !j
            if(idry==1) then
              out4(:,:)=rjunk
            else !element wet
              !Compute z-coordinates
              if(ivcor==1) then !localized
                do j=1,3
                  nd=node3(i,j)
                  do k=kbp(nd)+1,nvrt-1
                    ztmp2(k,j)=(eta2(nd)+dp(nd))*sigma_lcl(k,nd)+eta2(nd)
                  enddo !k
                  ztmp2(kbp(nd),j)=-dp(nd) !to avoid underflow
                  ztmp2(nvrt,j)=eta2(nd) !to avoid underflow
                enddo !j

                ztmp=0
                kbpl=minval(kbp(node3(i,1:3)))
                do k=kbpl,nvrt
                  do j=1,3
                    nd=node3(i,j)
                    ztmp(k)=ztmp(k)+arco(i,j)*ztmp2(max(k,kbp(nd)),j)
                  enddo !j
                enddo !k

              else if(ivcor==2) then !SZ
                call zcor_SZ(dep,etal,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
     &sigma,ztmp(:),idry2,kbpl)
              endif

              if(1==2) then
!--------------------------------------
              do k=kz,nvrt
                kin=k-kz+1
                hmod2=min(dep,h_s)
                if(hmod2<=h_c) then
                  ztmp(k)=sigma(kin)*(hmod2+etal)+etal
                else if(etal<=-h_c-(hmod2-h_c)*theta_f/sinh(theta_f)) then
                  write(*,*)'Pls choose a larger h_c (2):',etal,h_c
                  stop
                else
                  ztmp(k)=etal*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
                endif

!               Following to prevent underflow
                if(k==kz) ztmp(k)=-hmod2
                if(k==nvrt) ztmp(k)=etal
              enddo !k

              if(dep<=h_s) then
                kbpl=kz
              else !z levels
!               Find bottom index
                kbpl=0
                do k=1,kz-1
                  if(-dep>=ztot(k).and.-dep<ztot(k+1)) then
                    kbpl=k
                    exit
                  endif
                enddo !k
                if(kbpl==0) then
                  write(*,*)'Cannot find a bottom level:',dep,i
                  stop
                endif
                ztmp(kbpl)=-dep
                do k=kbpl+1,kz-1
                  ztmp(k)=ztot(k)
                enddo !k
              endif

              do k=kbpl+1,nvrt
                if(ztmp(k)-ztmp(k-1)<=0) then
                  write(*,*)'Inverted z-level:',etal,dep,ztmp(k)-ztmp(k-1)
                  stop
                endif
              enddo !k
!------------------------------------
              endif !1==2
       
              do k=kbpl,nvrt
                do m=1,ivs
                  do j=1,3
                    nd=node3(i,j)
                    kin=max(k,kbp00(nd))
                    out2(i,k,m)=out2(i,k,m)+arco(i,j)*out(i,j,kin,m)
                  enddo !j
                enddo !m
!                write(65,*)i,k,ztmp(k),(out2(i,k,m),m=1,ivs)     
              enddo !k

!             Interplate in vertical
              do kk=1,ntran
                k0=0
                do k=kbpl,nvrt-1
                  if(z0(kk)>=ztmp(k).and.z0(kk)<=ztmp(k+1)) then
                    k0=k
                    rat=(z0(kk)-ztmp(k))/(ztmp(k+1)-ztmp(k))
                    exit
                  endif
                enddo !k
                if(k0==0) then
!                  write(12,*)'Warning: failed to find a vertical level:',it,i
                  out4(kk,1:ivs)=rjunk
                else
                  do m=1,ivs
                    out4(kk,m)=out2(i,k0,m)*(1-rat)+out2(i,k0+1,m)*rat
        
                    !Debug
                    !if(mod(its,iskip)==0) write(99,*)time,i,kk,out4(kk,m)
                  enddo !m
                endif
              enddo !kk
            endif !dry/wet

            if(mod(its,iskip)==0) then
              iout3=its/iskip
              do kk=1,ntran
                itmp=(i-1)*ntran+kk
!               out3(nxz,2,nsteps)
                out3(itmp,1:ivs,iout3)=out4(kk,1:ivs)
              enddo !kk
            endif !mod
          enddo !i=1,nxy

          if(mod(its,iskip)==0) then
            !Along all vertical levels of each build pt
            iout3=its/iskip
            if(ifilter==0) then
              do j=1,nxz
                write(21,'(e22.12,2(1x,f10.3))')time,out3(j,1:ivs,iout3)
              enddo !j
            endif
          endif !mod
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday
 
      if(iout3/=nstep) then
        write(*,*)'Mismatch in time steps'
        stop
      endif

      if(ifilter/=0) then
!       Filter time series
        do i=1,nxz
          do j=1,ivs
            ser(1:iout3)=out3(i,j,1:iout3)
            call avstrp(ser,iout3,24,ser2,nxf)
            call avstrp(ser2,nxf,24,ser,nxf2)
            call avstrp(ser,nxf2,25,ser2,nxf3)
            out5(i,j,1:nxf3)=ser2(1:nxf3)
          enddo !j
        enddo !i

        print*, 'done filtering...'

!       Output
        do it=1,nxf3
          write(65,rec=irec3+1) it*3600.
          write(65,rec=irec3+2) it
          irec3=irec3+2
          do i=1,nxz
            write(65,rec=irec3+i) 0.
          enddo !i
          irec3=irec3+nxz
 
          do i=1,nxz
            do j=1,ivs
              write(65,rec=irec3+1) out5(i,j,it)
              irec3=irec3+1
            enddo !j
          enddo !i
        enddo !it
      endif !ifilter

      stop
      end subroutine

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

