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
!	Read in (x,y,time) from station.xyt (time in sec); e.g., casts; 
!       for 3D variables (surface values for 2D variables) DEFINED AT NODES. Interpolation in time, and
!       add extra times before and after to examine phase errors.
!       Works for mixed tri/quad outputs.
!       Inputs: (1) nc files;
!               (2) station.xyt (bp format): make sure all times (in sec) are after 1st record (to ensure interpolation in time); 
!                                pad extra days before and after if necessary.
!               (3) read_output8_xyt.in: 
!                   variable name (e.g. elev)\n 
!                   invalid value (for out of domain, dry etc)\n 
!                   window (hours) b4 and after the cast, stride (hrs) - used to 
!                   examine the phase error. If window=0, each cast is repeated twice 
!                   there are no other extra casts.
!                (4) vgrid.in (in this dir or ../)
!       Outputs: fort.1[89]; fort.11 (fatal errors); fort.12: nonfatal errors.
!                The total # of 'virtual' casts for each actual cast is 2*window/stride+2
!									
! ifort -mcmodel=medium -assume byterecl -O2 -o read_output8_xyt.WW ../UtilLib/extract_mod.f90 read_output8_xyt.f90 ../UtilLib/compute_zcor.f90 ../UtilLib/pt_in_poly.f90 -I$NETCDF/include -I$NETCDF_FORTRAN/include -L$NETCDF_FORTRAN/lib -L$NETCDF/lib -lnetcdf -lnetcdff
!****************************************************************************************
!
      program read_out
      use netcdf
      use extract_mod

      character(len=30) :: file63,varname
      character(len=12) :: it_char
      integer,allocatable :: kbp(:),node3(:,:),iep(:),iday(:,:),irecord(:,:)
      dimension swild(3),out3(2,2)
      real,allocatable :: sigma(:),cs(:),ztot(:), &
    &out(:,:,:),out2(:,:,:),eta2(:),arco(:,:),ztmp(:),x00(:),y00(:), &
    &out4(:,:),t00(:),times(:,:),sigma_lcl(:,:),ztmp2(:,:),outvar(:,:,:)
      real*8,allocatable :: timeout(:)

      integer :: nodel(3),dimids(100),idims(100), &
     &start_2d(2),start_3d(3),start_4d(4), &
     &count_2d(2),count_3d(3),count_4d(4)
      !long int for large files
      !integer(kind=8) :: irec
      
      open(10,file='read_output8_xyt.in',status='old')
      read(10,'(a30)')varname
!     Junk used for 3D variables: below bottom; dry spot; no parents
      read(10,*)rjunk
      read(10,*)window,wstride !in hours
      close(10)
      varname=adjustl(varname); len_var=len_trim(varname)

      if(wstride==0) stop 'wstride=0'
      nextra=2*window/wstride+1 !extra casts in addition to each cast (for phase error)
      
!...  Header
      call readheader('schout_1.nc')
      !Returned vars: ne,np,ns,nrec,start_time,[x y dp ](np),
      !elnode,i34,nvrt,itime_id,ielev_id,h0,dtout

      print*, 'After header:',ne,np,nrec,i34(ne),itime_id,ielev_id, &
     &elnode(1:i34(ne),ne),nvrt,h0,x(np),y(np),dp(np),start_time

!     Read in station.xyt
      open(10,file='station.xyt',status='old')
      read(10,*) 
      read(10,*) nxy
      nxy2=nxy*(1+nextra)
      allocate(x00(nxy2),y00(nxy2),t00(nxy2),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do i=1,nxy
        read(10,*)k,xtmp,ytmp,ttmp
!       Check if time is before first record
        if(ttmp<dtout) then
          write(11,*)'Time before first record; try to pad extra day:',i,ttmp
          stop
        endif

        indx=(i-1)*(1+nextra)+1
        x00(indx)=xtmp
        y00(indx)=ytmp
        t00(indx)=ttmp !time in sec
        !Add extra casts
        do j=1,nextra
          indx=indx+1
          x00(indx)=xtmp
          y00(indx)=ytmp
          t00(indx)=max(dtout,ttmp-window*3600+(j-1)*wstride*3600) !also need to ensure it's <last time
        enddo !j
      enddo !i
      close(10)
      nxy=nxy2

!      print*, 'i23d=',i23d,' nrec= ',nrec

!     Read in vgrid.in
      allocate(timeout(nrec),ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np),outvar(2,nvrt,np), &
    &node3(nxy,3),arco(nxy,3),iep(nxy))
      call get_vgrid('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

      allocate(ztmp(nvrt),ztmp2(nvrt,3),out(3,nvrt,2),out2(2,nvrt,2),out4(nvrt,2),iday(2,nxy), &
    &irecord(2,nxy),times(2,nxy),stat=istat)
      if(istat/=0) stop 'Falied to allocate (2)'

!     Calculate kbp00
      if(ivcor==1) then
        kbp00=kbp
      else
        do i=1,np
          !Use large eta to get true bottom
          call zcor_SZ(dp(i),1.e8,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma, &
     &ztmp(:),idry2,kbp00(i))
        enddo !i
      endif !ivcor

!...  Find parent element for (x00,y00)
      iep=0
      arco=1./3 !initialize for pts without parents
      do l=1,nxy
        node3(l,1:3)=elnode(1:3,1) !initialize for pts without parents
      enddo !l

      do i=1,ne
        do l=1,nxy
          if(iep(l)/=0) cycle

          call pt_in_poly(i34(i),x(elnode(1:i34(i),i)),y(elnode(1:i34(i),i)),x00(l),y00(l),inside,arco(l,1:3),nodel)
          if(inside==1) then
            iep(l)=i
            !print*, 'Found:',l,arco(l,1:3),nodel
            node3(l,1:3)=elnode(nodel(1:3),i)
          endif !inside
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
          write(12,*)'Cannot find a parent for pt:',j,x00(j),y00(j)
!          stop
        endif
      enddo !j

!...  Compute stack and record # for each pt
      do i=1,nxy
!       Check if time is before first record
!        if(t00(i)<dtout) then
!          write(11,*)'Time before first record; try to padd extra day (0):',i,t00(i)
!!'
!          stop
!        endif

!       Lower and upper bound stacks and record #s for t00(i)
        iday(1,i)=(t00(i)-dtout)/nrec/dtout+1
        if(iday(1,i)<1) then
          write(11,*)'Make sure cast time is after 1st record:',i,t00(i) 
          stop
        else
          irecord(1,i)=(t00(i)-(iday(1,i)-1)*nrec*dtout)/dtout
          !Bounding record time just b4 the cast time, corresponding to record
          !irecord(1,i) in stack iday(1,i)
          times(1,i)=((iday(1,i)-1)*nrec+irecord(1,i))*dtout 
          iday(2,i)=t00(i)/nrec/dtout+1
          irecord(2,i)=(t00(i)-(iday(2,i)-1)*nrec*dtout)/dtout+1
          !Bounding record time just after the cast time, corresponding to record 
          !irecord(2,i) in stack iday(2,i). Note that irecord(2,i)
          !may<irecord(1,i) (e.g. t00 before 1st record of iday(2,i))
          times(2,i)=((iday(2,i)-1)*nrec+irecord(2,i))*dtout 
        endif

        if(irecord(1,i)>nrec.or.irecord(2,i)>nrec) then
          write(11,*)'Record # overflow: ',i,irecord(:,i)
          stop
        endif
        if(t00(i)<times(1,i).or.t00(i)>times(2,i)) then
          write(11,*)'Wrong time bounds:',i,t00(i),times(:,i),iday(:,i),irecord(:,i)
          stop
        endif
      enddo !i=1,nxy

!...  Time iteration
!...
      do i=1,nxy
        loop1: do l=1,2 !2 times
          write(it_char,'(i12)')iday(l,i)
          it_char=adjustl(it_char); leng=len_trim(it_char)
          file63='schout_'//it_char(1:leng)//'.nc'
          iret=nf90_open(trim(adjustl(file63)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid)
          !time is double
          iret=nf90_inq_varid(ncid,'time',itime_id)
          iret=nf90_get_var(ncid,itime_id,timeout,(/1/),(/nrec/))
!      print*, 'time=',timeout,trim(adjustl(file63))

          irec=irecord(l,i)
          !Get elev
          start_2d(1)=1; start_2d(2)=irec
          count_2d(1)=np; count_2d(2)=1
          iret=nf90_get_var(ncid,ielev_id,eta2,start_2d,count_2d)

          iret=nf90_inq_varid(ncid,varname(1:len_var),ivarid1)
          if(iret/=nf90_NoErr) stop 'Var not found'
          iret=nf90_Inquire_Variable(ncid,ivarid1,ndims=ndims,dimids=dimids)
          if(ndims>100) stop 'increase dimension of dimids & idims'
          do ii=1,ndims
            iret=nf90_Inquire_Dimension(ncid,dimids(ii),len=idims(ii))
          enddo !i
          if(idims(ndims-1)/=np) stop 'can only handle node-based'
          if(idims(ndims)/=nrec) stop 'last dim is not time'

          iret=nf90_get_att(ncid,ivarid1,'i23d',i23d)
          iret=nf90_get_att(ncid,ivarid1,'ivs',ivs)
          !print*, 'i23d:',i23d,ivs,idims(1:ndims)
          if(ivs==1) then !scalar
            if(mod(i23d-1,3)==0) then !2D
              start_2d(1)=1; start_2d(2)=irec
              count_2d(1)=np; count_2d(2)=1
              iret=nf90_get_var(ncid,ivarid1,outvar(1,1,:),start_2d,count_2d)
            else !3D
              start_3d(1:2)=1; start_3d(3)=irec
              count_3d(2)=np; count_3d(1)=nvrt; count_3d(3)=1
              iret=nf90_get_var(ncid,ivarid1,outvar(1,:,:),start_3d,count_3d)
            endif
          else !vector
            if(mod(i23d-1,3)==0) then !2D
              start_3d(1:2)=1; start_3d(3)=irec
              count_3d(2)=np; count_3d(1)=2; count_3d(3)=1
              iret=nf90_get_var(ncid,ivarid1,outvar(1:2,1,:),start_3d,count_3d)
            else if(ndims-1==3) then !3D vector
              start_4d(1:3)=1; start_4d(4)=irec
              count_4d(3)=np; count_4d(2)=nvrt; count_4d(1)=2; count_4d(4)=1
              iret=nf90_get_var(ncid,ivarid1,outvar(:,:,:),start_4d,count_4d)
            else
              stop 'Unknown type(2)'
            endif
          endif !ivs

          !Available now: outvar(2,nvrt,np)


!          open(63,file=it_char//'_'//file63,status='old',access='direct',recl=nbyte)
!
!          if(i23d==2) then
!            irec=irec0+(irecord(l,i)-1)*(2+np+np*ivs)
!          else
!            irec=irec0+(irecord(l,i)-1)*(2+np+icum(np,nvrt,ivs))
!          endif
!
!!          read(63,rec=irec+1) time
!!          read(63,rec=irec+2) it
!          irec=irec+2
!
!!         Only read in 3 nodes for efficiency
!          do j=1,3
!            nd=node3(i,j)
!            read(63,rec=irec+nd) eta2(nd)
!          enddo !j
!!          do ii=1,np
!!            read(63,rec=irec+ii) eta2(ii)
!!          enddo !ii
!          irec=irec+np
!
          out2(l,:,:)=0
          if(mod(i23d-1,3)==0) then !2D
            do j=1,3 !nodes
              nd=node3(i,j)
              do m=1,ivs
!                read(63,rec=irec+(nd-1)*ivs+m) tmp
                out2(l,1,m)=out2(l,1,m)+arco(i,j)*outvar(m,1,nd)
              enddo !m
            enddo !j
          else !3D
!            do j=1,3 !nodes
!              nd=node3(i,j)
!              do k=max0(1,kbp00(nd)),nvrt
!                do m=1,ivs
!                  read(63,rec=irec+icum(nd,k,m)) out(j,k,m)
!                enddo !m
!              enddo !k
!            enddo !j

!           Do interpolation
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
              out2(:,:,:)=rjunk
              exit loop1
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
       
              do k=kbpl,nvrt
                do m=1,ivs
                  do j=1,3
                    nd=node3(i,j)
                    kin=max(k,kbp00(nd))
                    out2(l,k,m)=out2(l,k,m)+arco(i,j)*outvar(m,kin,nd) !out(j,kin,m)
                  enddo !j
                enddo !m
              enddo !k
            endif !dry/wet
          endif !i23d
        enddo loop1 !l=1,2; 2 times

!       Interpolate in time
        trat=(t00(i)-times(1,i))/(times(2,i)-times(1,i)) !must be [0,1]
        if(mod(i23d-1,3)==0) then
          if(iep(i)==0) then !no parents
            out4(1,1:ivs)=rjunk
          else
            out4(1,1:ivs)=out2(1,1,1:ivs)*(1-trat)+out2(2,1,1:ivs)*trat
          endif
          write(18,'(e16.8,2(1x,f12.3))')t00(i)/86400,out4(1,1:ivs)
        else !3D
          if(iep(i)==0) then !no parents
            out4(:,1:ivs)=rjunk
          else
            out4(kbpl:nvrt,1:ivs)=out2(1,kbpl:nvrt,1:ivs)*(1-trat)+out2(2,kbpl:nvrt,1:ivs)*trat
            !Extend
            do k=1,kbpl-1
              out4(k,1:ivs)=out4(kbpl,1:ivs)
              ztmp(k)=ztmp(kbpl)
            enddo !k
          endif

          do k=nvrt,1,-1
            !First of each cast suite is at the actual cast time (followed by b4 and after)
            write(18,'(i6,4(1x,f12.3))')i,out4(k,1),ztmp(k)-ztmp(nvrt),ztmp(k),t00(i)/86400
            if(ivs==2) write(19,'(i6,4(1x,f12.3))')i,out4(k,2),ztmp(k)-ztmp(nvrt),ztmp(k),t00(i)/86400
          enddo !k
        endif
      enddo !i=1,nxy

      print*, 'Finished!'

      stop
      end

