!====================================================================
!  Routines to calculate various types of z-coord. in _3D_ SELFE
!  zcor_SZ
!  zcor_PWS
!  zcor_VQS
!  get_vgrid

!  ifort -Bstatic -O3 -c compute_zcor.f90
!  pgf90 -O2 -mcmodel=medium  -Bstatic -c compute_zcor.f90
!====================================================================
!====================================================================

      subroutine zcor_SZ(dp,eta,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma,zcor,idry,kbp)
!     Routine to compute z coordinates for SELFE's SZ vertical system
!     Inputs:
!             dp: depth;
!             eta: elevation;
!             h0: min. depth;
!             h_s: transition depth between S and Z layers;
!             h_c: transition depth between S and sigma
!             theta_b, theta_f: spacing const. in S coordinate system;
!             nvrt: total # of vertical levels (S+Z);
!             kz: # of Z levels (1 if pure S);
!             ztot(1:kz):  z coordinates for Z levels; note that ztot(kz)=-h_s;
!             sigma(1:nvrt-kz+1): sigma coordinates for S (or sigma) levels; note that sigma(1)=-1, sigma(nvrt-kz+1)=0;
!     Outputs:
!             idry: wet (0) or dry (1) flag;
!             kbp: bottom index (0 if dry);
!             zcor(kbp:nvrt): z coordinates (junk if dry);    
!      implicit real*8(a-h,o-z)
      integer, intent(in) :: kz,nvrt
      real, intent(in) :: dp,eta,h0,h_s,h_c,theta_b,theta_f,ztot(nvrt),sigma(nvrt)
      integer, intent(out) :: idry,kbp
      real, intent(out) :: zcor(nvrt)

      !Local
      real :: s_con1,cs(nvrt)

!     Sanity check done before
!     Pre-compute constants
      s_con1=sinh(theta_f)
      nsig=nvrt-kz+1 !# of S levels 
      do k=1,nsig
        cs(k)=(1-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
     &theta_b*(tanh(theta_f*(sigma(k)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
      enddo !k

      zcor=-99
      if(eta<=h0-h_s) then
        stop 'Deep depth dry'
      else if(eta+dp<=h0) then
        idry=1; kbp=0
      else !wet
!       S-levels
        idry=0
        hmod=min(h_s,dp)
        do k=kz,nvrt
          kin=k-kz+1
          if(hmod<=h_c) then
            zcor(k)=sigma(kin)*(hmod+eta)+eta
          else if(eta<=-h_c-(hmod-h_c)*theta_f/s_con1) then !hmod(i)>h_c>=0
            print*, 'Pls choose a larger h_c (2):',eta,h_c
            stop
          else
            zcor(k)=eta*(1+sigma(kin))+h_c*sigma(kin)+(hmod-h_c)*cs(kin)
          endif
        enddo !k=kz,nvrt

!       z-levels
        if(dp<=h_s) then
          kbp=kz
        else !bottom index 
          kbp=0 !flag
          do k=1,kz-1
            if(-dp>=ztot(k).and.-dp<ztot(k+1)) then
              kbp=k
              exit
            endif
          enddo !k
          if(kbp==0) then
            print*, 'Cannot find a bottom level for node (3):',i
            stop
          endif

          if(kbp>=kz.or.kbp<1) then
            print*, 'Impossible 92:',kbp,kz
            stop
          endif
          zcor(kbp)=-dp
          do k=kbp+1,kz-1
            zcor(k)=ztot(k)
          enddo !k
        endif

        do k=kbp+1,nvrt
          if(zcor(k)-zcor(k-1)<=0) then
            write(*,*)'Inverted z-levels at:',i,k,zcor(k)-zcor(k-1),eta,hmod
            stop
          endif
        enddo !k
      endif !wet ot dry

      end subroutine zcor_SZ

!====================================================================
      subroutine zcor_PWS(etal,dep,h0,nvrt,m_pws,h_tran1,h_tran2,a_pws,a_pws2,kkm,hsm,sigma,ztmp,idry,kbpl)
!     Piece-wise sigma - not used at the moment
!      implicit real*8(a-h,o-z)
      integer, intent(in) :: nvrt,m_pws,kkm(0:m_pws+1)
      real, intent(in) :: etal,dep,h0,h_tran1,h_tran2,a_pws,a_pws2,hsm(m_pws),sigma(nvrt)
      integer, intent(out) :: idry,kbpl
      real, intent(out) :: ztmp(nvrt)

      !Local
      real :: z_sigma(nvrt),z_pws(nvrt),sigma2(nvrt)

      if(etal<=-hsm(m_pws)) then
        write(*,*)'ZCOOR: elev<hsm:',etal
        stop
      endif

      !Pre-compute sigma2
      do k=1,nvrt
        sigma2(k)=a_pws*sigma(k)*sigma(k)+(1+a_pws)*sigma(k) !transformed sigma
      enddo !k

      ztmp=-99 !init
      if(dep+etal<=h0) then
        idry=1; kbpl=0
      else !wet
        idry=0; kbpl=1

!       z_sigma (used in shallow)
        do k=1,nvrt
          tmp=a_pws2*sigma(k)*sigma(k)+(1+a_pws2)*sigma(k) !transformed sigma
          z_sigma(k)=(etal+dep)*tmp+etal
        enddo !k

        z_pws(kbpl)=-dep !to avoid underflow
        z_pws(nvrt)=etal !to avoid underflow
        do m=m_pws,0,-1
          z0=z_pws(kkm(m+1)+1) !upper bound for the stencil
          if(m==0) then
            z_1=-dep
          else !m>0
            tmp=-hsm(m)*tanh(dble(dep/hsm(m))) !proposed lower bound for the stencil
            if(dep<=0.1d0.or.tmp>=z0-0.1d0) then
              z_1=-min(hsm(m),dep) !lower bound for the stencil
            else
              z_1=tmp
            endif
          endif !m
          !assert?
          if(-z_1>dep.or.z_1>=z0) then
            write(*,*)'ZCOOR: Below bottom:',m,z_1,dep,z0
            stop
          endif
          do k=kkm(m)+1,kkm(m+1)
            sp=(sigma2(k)-sigma2(kkm(m+1)+1))/(sigma2(kkm(m+1)+1)-sigma2(kkm(m)))
            if(sp<=-1.or.sp>=0) then
              write(*,*)'ZCOOR: Out of bound (1):',m,k,sp
              stop
            endif
            z_pws(k)=(z0-z_1)*sp+z0 
          enddo !k
        enddo !m

        ztmp(kbpl)=-dep !to avoid underflow
        ztmp(nvrt)=etal !to avoid underflow
        do k=2,nvrt-1
          tmp=(1+tanh(dble((dep+etal-h_tran1)/h_tran2)))/2 !blend ratio
          ztmp(k)=z_sigma(k)+(z_pws(k)-z_sigma(k))*tmp
        enddo !k

        do k=2,nvrt
          if(ztmp(k)-ztmp(k-1)<=0) then
            write(*,*)'Inverted z-levels at (3):',k
            stop
          endif
        enddo !k
      endif !dep+etal<=h0

      end subroutine zcor_PWS

!====================================================================
      subroutine zcor_VQS(inode,etal,dep,h0,sigma_lcl,nvrt,kbp,ztmp,idry,kbpl)
!     Localized sigma - not used at the moment
!      implicit real*8(a-h,o-z)
      integer, intent(in) :: inode,nvrt,kbp
      real, intent(in) :: etal,dep,h0,sigma_lcl
      integer, intent(out) :: idry,kbpl
      real, intent(out) :: ztmp(nvrt)

      if(etal+dep<=h0) then
        kbpl=0; idry=1
      else !wet
        kbpl=kbp; idry=0
        do k=kbpl+1,nvrt-1
          ztmp(k)=(etal+dep)*sigma_lcl+etal
        enddo !k
        ztmp(kbpl)=-dep !to avoid underflow
        ztmp(nvrt)=etal !to avoid underflow
      endif !wet/dry

      end subroutine zcor_VQS

!====================================================================
      subroutine get_vgrid(vgrid,np,nvrt1,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)
!     Routine to read in vgrid.in (either in the current dir or ../) 
!     Not all outputs have meaningful values, depending on ivcor.
!      implicit real*8(a-h,o-z)
      character(len=*), intent(in) :: vgrid
      integer, intent(in) :: np,nvrt1
      integer, intent(out) :: ivcor,kz,kbp(np)
      real, intent(out) :: h_s,h_c,theta_b,theta_f,ztot(nvrt1),sigma(nvrt1),sigma_lcl(nvrt1,np)

      !local
      logical :: lexist
      
      inquire(file=vgrid,exist=lexist)
      if(lexist) then
        open(19,file=vgrid,status='old')
      else
        inquire(file='../'//vgrid,exist=lexist)
        if(lexist) then
          open(19,file='../'//vgrid,status='old')
        else
          write(*,*)'Unable to find vgrid.in'
          stop
        endif
      endif

      read(19,*)ivcor
      select case(ivcor)
        case(2) !SZ
          read(19,*) nvrt,kz,h_s !kz>=1
          if(nvrt<3) stop 'nvrt<3'
          if(kz<1.or.kz>nvrt-2) then
            write(*,*)'Wrong kz:',kz
            stop
          endif
          if(h_s<6) then
            write(*,*)'h_s needs to be larger:',h_s
            stop
          endif

          ! Allocate vertical layers arrays
          !allocate(ztot(nvrt),sigma(nvrt),stat=istat)
          !if(istat/=0) stop 'get_vgrid: ztot allocation failure'

          ! # of z-levels excluding "bottom" at h_s
          read(19,*) !for adding comment "Z levels"
          do k=1,kz-1
            read(19,*)j,ztot(k)
            if(ztot(k)>=-h_s) then
              print*, 'Illegal Z level:',k
              stop
            endif
            if(k>1) then; if(ztot(k)<=ztot(k-1)) then
              print*, 'z-level inverted:',k
              stop
            endif; endif
          enddo !k
          read(19,*) !level kz       
          ! In case kz=1, there is only 1 ztot(1)=-h_s
          ztot(kz)=-h_s

          nsig=nvrt-kz+1 !# of S levels (including "bottom" & f.s.)
          read(19,*) !for adding comment "S levels"
          read(19,*)h_c,theta_b,theta_f
          if(h_c<5.or.h_c>=h_s) then !large h_c to avoid 2nd type abnormality
            print*, 'h_c needs to be larger avoid 2nd type abnormality; &
     &do u want to continue? Enter 1 to continue, or ctrl-C to abort:'
            read*, itmp
          endif
          if(theta_b<0.or.theta_b>1) then
            write(*,*)'Wrong theta_b:',theta_b
            stop
          endif
          if(theta_f<=0) then
            write(*,*)'Wrong theta_f:',theta_f
            stop
          endif

          sigma(1)=-1 !bottom
          sigma(nsig)=0 !surface
          read(19,*) !level kz
          do k=kz+1,nvrt-1
            kin=k-kz+1
            read(19,*) j,sigma(kin)
            if(sigma(kin)<=sigma(kin-1).or.sigma(kin)>=0) then
              write(*,*)'Check sigma levels at:',k,sigma(kin),sigma(kin-1)
              stop
            endif
          enddo !k
          read(19,*) !level nvrt

        case(1) !localized sigma
          read(19,*)nvrt
          do i=1,np
            read(19,*)j,kbp(i),sigma_lcl(kbp(i):nvrt,i)
          enddo !i
        case default
          write(*,*)'Unknown ivcor:',ivcor
          stop
      end select
      close(19)      

      end subroutine get_vgrid

!====================================================================
