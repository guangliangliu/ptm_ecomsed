      program efloats
      use hydro_mod
      use particle_mod
      use param_mod
      implicit none
      integer  ::   idum
          
      call scan
      
      do layer=1,kb
            do iint=irelst,irelend
      
                  if (iint.ge.irelst) then
                        if (iint.eq.irelst) idum=-1
                        call input
                        call partrak(idum)
                  endif
      
                  if (mod(iint-irelst+1,floor(3600.e0/dti))==0.or.iint==irelend) then
                        call archive
                  endif
      
            enddo
      enddo
      
      end program
      
      subroutine scan
      use particle_mod
      use hydro_mod
      implicit none
      integer     ::    ixm, jym, kbm
      integer     ::    iugrd, iufloats, iucor
      character (len=80)      :: head
      integer     ::    fileend
      iugrd =     55
      iufloats =  56
      iucor =     57
      
      
      !
      open(iugrd,file='../../Share/model_grid')
      read(iugrd,'(a80)') head
      read(iugrd,'(a80)') head
      read(iugrd, *) kbm

      do k=1,kbm
            read(iugrd,*) z(k)
      enddo

      do k=1,kbm-1
            dz(k)=z(k)-z(k+1)
            zz(k)=.5*(z(k)+z(k+1))
      enddo
      dz(kb)=dz(kb-1)

      do k=1,kbm-2
            dzz(k)=zz(k)-zz(k+1)
      enddo
      dzz(kb-1)=dzz(kb-2)
	  dzz(kb)=dzz(kb-2)

      read(iugrd,'(a80)') head
      read(iugrd,*)  ixm, jym
      h1   = 0.0
      h2   = 0.0
      h     = 0.0
      ang   = 0.0
      cor  = 0.0
      datum= 0.0
      ygrid = 0.0
      fileend='true'
      do
        read(iugrd,*,iostat=fileend) &
         i,j,h1(i,j),h2(i,j),h(i,j),ang(i,j),ygrid(i,j),xgrid(i,j), & 
         datum(i,j)
            if (fileend<0.e0) then
                  exit
            endif
      enddo
      fileend='true'
      close(iugrd)
      !
      
      open (unit=iufloats,file='partrack.inp',form='formatted')
      read (iufloats,'(a80)') head
      read (iufloats,*) nfreq,npart,nconv,irelst,irelend,npclass,nsource
      do mm=1,nsource
            read (iufloats,*) isource(mm),jsource(mm),ksource(mm)
      enddo
      close(iufloats)

      open (iucor,file='../../Share/corner_loc',form='formatted',status='old')
      do 
        read(iucor,*,iostat=fileend)    i,j,xcor(i,j),ycor(i,j)
        if (fileend<0.e0) then
           exit
        endif
      enddo
      
      end subroutine scan
 
      subroutine partrak(idum)
!
!*************************************************************************
!     ecomsed model
!     version 1.3
!     february 2002
!*************************************************************************
!               copyright (c) 2002, hydroqual, inc.                      *
!                                                                        *
!  these coded instructions, statements, and computer programs  contain  *
!  unpublished  proprietary  information of hydroqual, inc., and         *
!  are protected by federal copyright law.  they  may  not be disclosed  *
!  to  third  parties  or copied or duplicated in any form, in whole or  *
!  in part, without the prior written consent of hydroqual, inc.         *
!                                                                        *
! point of contact: ecomsed-support@hydroquap.com                        *
!*************************************************************************
!
      use particle_mod
      use hydro_mod
      use param_mod
      implicit none
      real  ::    xstart, ystart, zstart
      real  ::    xend, yend, zend
      integer     ::    intip
      integer     ::    iout
      integer        ::   idum
!
      if (iint.lt.irelst) return
      intip=iint-irelst
!
!---> calculate how many groups of particles in the system
!
      ngradeloop=intip/nfreq + 1
      if (ngradeloop.ge.nconv) ngradeloop=nconv
!
!---> advance all groups to a higher grade before new group is introduced
!

      if (mod(intip,nfreq).ne.0.or.iint.gt.irelend) then
      else
            if(ngradeloop.eq.1) then

            elseif (ngradeloop.gt.1) then
                  do ll=1,nsource
                        do nn=ngradeloop,2,-1
                              do mm=1,npart
                                    xp(ll,mm,nn)=   xp(ll,mm,nn-1)
                                    yp(ll,mm,nn)=   yp(ll,mm,nn-1)
                                    zp(ll,mm,nn)=   zp(ll,mm,nn-1)
                                    inout(ll,mm,nn)=inout(ll,mm,nn-1)
                              enddo
                        enddo
                  enddo
            endif
                  !
      
                  !---> introduce a new group every nfreq timestep
                  !
            do ll=1,nsource
                  do mm=1,npart
                        xp(ll,mm,1)= isource(ll)+0.5
                        yp(ll,mm,1)= jsource(ll)+0.5
	                  !liang 20120617
                        if (layer==1) then
                              zp(ll,mm,1)= z(1)-1.e-3
                        elseif (layer>1) then
                              zp(ll,mm,1)= zz(layer-1)
                        endif
				!liang 20120617
                        inout(ll,mm,1)  = 0
                  enddo
            enddo

      endif
!

          if (iint==irelst) then
            call archive
          endif

		  continue
!---> advance all particles in all possible groups at the current time
!
        ! to calculate the momentum balance
        ! liang 20120322
        ! call momentum
        
        do ll=1,nsource
            do nn=1,ngradeloop
                  do mm=1,npart
                        !
                        if(inout(ll,mm,nn).eq.1) then
                        else
                              !
                              xstart=xp(ll,mm,nn)
                              ystart=yp(ll,mm,nn)
                              zstart=zp(ll,mm,nn)
                              iout=inout(ll,mm,nn)
                              
                              !
                              call onepart(xstart,ystart,zstart,idum,&
                                  dti,xend,yend,zend,iout)
                              ! 

                              ! 
                              xp(ll,mm,nn)=xend
                              yp(ll,mm,nn)=yend
                              zp(ll,mm,nn)=zend
                              inout(ll,mm,nn)=iout
                        endif
                  enddo
            enddo
        enddo
! 
      return
      end
!
!**************************************************************************
      subroutine onepart(xstart,ystart,zstart,&
           idum,deltat,xend,yend,zend,iout)
!*************************************************************************
!     ecomsed model
!     version 1.3
!     february 2002
!*************************************************************************
!               copyright (c) 2002, hydroqual, inc.                      *
!                                                                        *
!  these coded instructions, statements, and computer programs  contain  *
!  unpublished  proprietary  information of hydroqual, inc., and         *
!  are protected by federal copyright law.  they  may  not be disclosed  *
!  to  third  parties  or copied or duplicated in any form, in whole or  *
!  in part, without the prior written consent of hydroqual, inc.         *
!                                                                        *
! point of contact: ecomsed-support@hydroquap.com                        *
!*************************************************************************
      use hydro_mod
      use particle_mod
      use param_mod
      implicit none

      real, intent(in) :: xstart, ystart, zstart, deltat
      real, intent(out):: xend, yend, zend
      integer, intent(out)    ::    iout
      real  ::  hscnu, vscnu
      integer     ::    idum
      integer     ::    ic, jc, kc
      
      integer     ::    i0u, j0u, k0u, i1u, j1u, k1u
      integer     ::    i0v, j0v, k0v, i1v, j1v, k1v
      integer     ::    i0w, j0w, k0w, i1w, j1w, k1w
      integer     ::    i0h, j0h, k0h, i1h, j1h, k1h
      
      real  ::    xc, yc, zc
      real  ::    xlocal, ylocal, zlocal
      real  ::    xdist, ydist, zdist
      real  ::    deltax, deltay, deltaz
      
      
      real  ::    ui0k0, ui1k0, uk0, ui0k1, ui1k1, uk1, up
      real  ::    vj0k0, vj1k0, vk0, vj0k1, vj1k1, vk1, vp
      real  ::    wj0k0, wj1k0, wk0, wj0k1, wj1k1, wk1, wp
      real  ::    h1j0,  h1j1, h1pxyz, h2j0, h2j1, h2pxyz
      real  ::    dj0, dj1, dpxyz
      real  ::    aamj0k0, aamj1k0, aamk0, aamj0k1, aamj1k1, aamk1, aamp
      real  ::    khj0k0, khj1k0, khk0, khj0k1, khj1k1, khk1, khp
      real  ::    udj0k0, udj1k0, udk0, udj0k1, udj1k1, udk1, udp
      real  ::    vdi0k0, vdi1k0, vdk0, vdi0k1, vdi1k1, vdk1, vdp
      real  ::    wdp
      real  ::    upmax,vpmax,wpmax
      real  ::    deltaxmax,deltaymax,deltazmax
      real  ::    udpmax,vdpmax,wdpmax 
      
      integer     ::    loop
      integer     ::    icnew, jcnew, kcnew
      real  ::    xcnew, ycnew, zcnew
      real  ::    icwall, xwall, xtowall, xbounce  
      real  ::    jcwall, ywall, ytowall, ybounce
      real  ::    zwall, ztowall, zbounce
      integer     ::    icnewest, jcnewest

      real  ::    gasdev

!
      hscnu=hprnu
      vscnu=vprnu

!     absolute coordinates in tranformed grid
!     for cell i:  x=[i,i+1) left close '[', right open ')' 
!     for cell j:  y=[j,j+1) left close '[', right open ')'
!     for cell k:  z=[z(k),z(k+1)] 

!     find cell indieces
      ic= xstart
      jc= ystart

      do k=1,kbm1
            if(zstart.le.z(k).and.zstart.ge.z(k+1)) then
                  kc=k
                  exit
            else
                  kc=0
            endif
      enddo
      if (kc==0) then
        write(*,*) 'trouble !', zstart, iint, kc    
      endif

!     coordinates of cell center
      xc= ic+0.5
      yc= jc+0.5
      zc= zz(kc)


!---> interpolate u  
!     find i-range for u 
      i0u   = ic
      i1u   = ic+1
      xlocal= xstart-xc
      xdist = abs(xlocal)

!     find j-range for u
      j0u   = jc
      ylocal= ystart-yc
      ydist = abs(ylocal)
      if(ylocal.gt.0.) j1u=j0u+1
      if(ylocal.le.0.) j1u=j0u-1

!     free-slip on land boundary
      if(fsm(ic,j1u).eq.0.) j1u=j0u

!     find k-range for u
      k0u   = kc
      zlocal= zstart-zc
      zdist = abs(zlocal)
      if(zlocal.gt.0.) k1u=k0u-1
      if(zlocal.le.0.) k1u=k0u+1

!     free-slip on land boundary
      if(k1u.lt.1) k1u=k0u
      if(k1u.gt.kbm1) k1u=k0u

!     interpolate u in k0 plane
      ui0k0=u(i0u,j0u,k0u)+ydist*(u(i0u,j1u,k0u)-u(i0u,j0u,k0u))
      ui1k0=u(i1u,j0u,k0u)+ydist*(u(i1u,j1u,k0u)-u(i1u,j0u,k0u))
      uk0  =ui0k0*(0.5-xlocal)+ui1k0*(0.5+xlocal)

!     interpolate u in k1 plane
      ui0k1=u(i0u,j0u,k1u)+ydist*(u(i0u,j1u,k1u)-u(i0u,j0u,k1u))
      ui1k1=u(i1u,j0u,k1u)+ydist*(u(i1u,j1u,k1u)-u(i1u,j0u,k1u))
      uk1  =ui0k1*(0.5-xlocal)+ui1k1*(0.5+xlocal)

!     interpolate u in z direction
      up=uk0 + (zdist/abs(0.5*dz(k1u)+0.5*dz(k0u))) * (uk1-uk0)

!---> interpolate v

!     find i-range for v
      i0v   = ic
      xlocal= xstart-xc
      xdist = abs(xlocal)
      if(xlocal.gt.0.) i1v=i0v+1
      if(xlocal.le.0.) i1v=i0v-1

!     free-slip on land
      if(fsm(i1v,jc).eq.0.) i1v=i0v

!     find j-range for v
      j0v   = jc
      j1v   = j0v+1
      ylocal= ystart-yc
      ydist = abs(ylocal)

!     find k-range for v
      k0v   = kc
      zlocal= zstart-zc
      zdist = abs(zlocal)
      if(zlocal.gt.0.) k1v=k0v-1
      if(zlocal.le.0.) k1v=k0v+1

!     free-slip on land boundary
      if(k1v.lt.1) k1v=k0v
      if(k1v.gt.kbm1)k1v=k0v

!     interpolate v in k0 plane
      vj0k0=v(i0v,j0v,k0v)+xdist*(v(i1v,j0v,k0v)-v(i0v,j0v,k0v))
      vj1k0=v(i0v,j1v,k0v)+xdist*(v(i1v,j1v,k0v)-v(i0v,j1v,k0v))
      vk0  =vj0k0*(0.5-ylocal)+vj1k0*(0.5+ylocal)

!     interpolate v in k1 plane
      vj0k1=v(i0v,j0v,k1v)+xdist*(v(i1v,j0v,k1v)-v(i0v,j0v,k1v))
      vj1k1=v(i0v,j1v,k1v)+xdist*(v(i1v,j1v,k1v)-v(i0v,j1v,k1v))
      vk1  =vj0k1*(0.5-ylocal)+vj1k1*(0.5+ylocal)

!     interpolate v in z direction
      vp=vk0 + (zdist/abs(0.5*dz(k1v)+0.5*dz(k0v))) * (vk1-vk0)


!---> interpolate w

!     find i-range for w
      i0w   = ic
      xlocal= xstart-xc
      xdist = abs(xlocal)
      if(xlocal.gt.0.) i1w=i0w+1
      if(xlocal.le.0.) i1w=i0w-1

!     free-slip on land
      if(fsm(i1w,jc).eq.0.)i1w=i0w

!     find j-range for w
      j0w   = jc
      ylocal= ystart-yc
      ydist = abs(ylocal)
      if(ylocal.gt.0.) j1w=j0w+1
      if(ylocal.le.0.) j1w=j0w-1

!     free-slip on land
      if(fsm(ic,j1w).eq.0.)j1w=j0w

      if(fsm(i1w,j1w).eq.0.)then
      i1w=i0w
      j1w=j0w
      endif

!     find k-range for w
      k0w   = kc
      k1w   = kc+1
      zlocal= zstart-zc
      zdist = abs(zlocal)

!     interpolate w in k0 plane
      wj0k0= w(i0w,j0w,k0w) + xdist * (w(i1w,j0w,k0w)-w(i0w,j0w,k0w))
      wj1k0= w(i0w,j1w,k0w) + xdist * (w(i1w,j1w,k0w)-w(i0w,j1w,k0w))
      wk0  = wj0k0          + ydist * (wj1k0-wj0k0)

!     interpolate w in k1 plane
      wj0k1= w(i0w,j0w,k1w) + xdist * (w(i1w,j0w,k1w)-w(i0w,j0w,k1w))
      wj1k1= w(i0w,j1w,k1w) + xdist * (w(i1w,j1w,k1w)-w(i0w,j1w,k1w))
      wk1  = wj0k1          + ydist * (wj1k1-wj0k1)

!     interpolate w in z direction
      wp=wk0*(0.5+zlocal/dz(kc)) + wk1*(0.5-zlocal/dz(kc))

!---> interpolate h1 and h2

!     find i-range for h1 and h2
      i0h   = ic
      xlocal= xstart-xc
      xdist = abs(xlocal)
      if(xlocal.gt.0.) i1h=i0h+1
      if(xlocal.le.0.) i1h=i0h-1

!     use value at one element
      if(fsm(i1h,jc).eq.0.)i1h=i0h

!     find j-range for h1 and h2
      j0h   = jc
      ylocal= ystart-yc
      ydist = abs(ylocal)
      if(ylocal.gt.0.) j1h=j0h+1
      if(ylocal.le.0.) j1h=j0h-1

!     use value at one element
      if(fsm(ic,j1h).eq.0.) j1h=j0h

      if(fsm(i1h,j1h).eq.0.)then
        i1h=i0h
        j1h=j0h
      endif

!     interpolate h1
      h1j0   = h1(i0h,j0h) + xdist * (h1(i1h,j0h)-h1(i0h,j0h))
      h1j1   = h1(i0h,j1h) + xdist * (h1(i1h,j1h)-h1(i0h,j1h))
      h1pxyz = h1j0        + ydist * (h1j1-h1j0)

!     interpolate h2
      h2j0   = h2(i0h,j0h) + xdist * (h2(i1h,j0h)-h2(i0h,j0h))
      h2j1   = h2(i0h,j1h) + xdist * (h2(i1h,j1h)-h2(i0h,j1h))
      h2pxyz = h2j0        + ydist * (h2j1-h2j0)

!     interpolate d
      dj0    = d(i0h,j0h) + xdist * (d(i1h,j0h)-d(i0h,j0h))
      dj1    = d(i0h,j1h) + xdist * (d(i1h,j1h)-d(i0h,j1h))
      dpxyz  = dj0        + ydist * (dj1-dj0)

      if(dpxyz.le.0) then
        write(*,*)dpxyz,xdist,ydist
        write(*,*)i0h,i1h,j0h,j1h
        write(*,*)d(i0h,j0h),d(i1h,j0h),d(i0h,j1h),d(i1h,j1h)
      endif

!     find k-range for aam, kh
      k0h   = kc
      zlocal= zstart-zc
      zdist = abs(zlocal)
      if(zlocal.gt.0.) k1h=k0h-1
      if(zlocal.le.0.) k1h=k0h+1

!     free-slip on land boundary
      if(k1h.lt.1) k1h=k0h
      if(k1h.gt.kbm1)k1h=k0h

!     interpolate aam in k0 plane
      aamj0k0=aam(i0h,j0h,k0h)+xdist*(aam(i1h,j0h,k0h)-aam(i0h,j0h,k0h))
      aamj1k0=aam(i0h,j1h,k0h)+xdist*(aam(i1h,j1h,k0h)-aam(i0h,j1h,k0h))
      aamk0  =aamj0k0         +ydist*(aamj1k0-aamj0k0)

!     interpolate aam in k1 plane
      aamj0k1=aam(i0h,j0h,k1h)+xdist*(aam(i1h,j0h,k1h)-aam(i0h,j0h,k1h))
      aamj1k1=aam(i0h,j1h,k1h)+xdist*(aam(i1h,j1h,k1h)-aam(i0h,j1h,k1h))
      aamk1  =aamj0k1         +ydist*(aamj1k1-aamj0k1)

!     interpolate aam in z direction
      aamp=aamk0 + (zdist/abs(0.5*dz(k1h)+0.5*dz(k0h))) * (aamk1-aamk0)


!     interpolate kh in k0 plane
      khj0k0=kh(i0h,j0h,k0h)+xdist*(kh(i1h,j0h,k0h)-kh(i0h,j0h,k0h))
      khj1k0=kh(i0h,j1h,k0h)+xdist*(kh(i1h,j1h,k0h)-kh(i0h,j1h,k0h))
      khk0  =khj0k0         +ydist*(khj1k0-khj0k0)

!     interpolate kh in k1 plane
      khj0k1=kh(i0h,j0h,k1h)+xdist*(kh(i1h,j0h,k1h)-kh(i0h,j0h,k1h))
      khj1k1=kh(i0h,j1h,k1h)+xdist*(kh(i1h,j1h,k1h)-kh(i0h,j1h,k1h))
      khk1  =khj0k1         +ydist*(khj1k1-khj0k1)

!     interpolate kh in z dirction
      khp=khk0 + (zdist/abs(0.5*dz(k1h)+0.5*dz(k0h))) * (khk1-khk0)



!---> interpolate ud and vd

!     interpolate ud at k0 plane
      udj0k0= h2(i1h,j0h)/h1(i1h,j0h)*aam(i1h,j0h,k0h)*d(i1h,j0h)- &
              h2(i0h,j0h)/h1(i0h,j0h)*aam(i0h,j0h,k0h)*d(i0h,j0h)
      udj1k0= h2(i1h,j1h)/h1(i1h,j1h)*aam(i1h,j1h,k0h)*d(i1h,j1h)- &
              h2(i0h,j1h)/h1(i0h,j1h)*aam(i0h,j1h,k0h)*d(i0h,j1h)
      udk0  = udj0k0        + ydist * (udj1k0-udj0k0)

!     interpolate ud at k1 plane
      udj0k1= h2(i1h,j0h)/h1(i1h,j0h)*aam(i1h,j0h,k1h)*d(i1h,j0h)- &
              h2(i0h,j0h)/h1(i0h,j0h)*aam(i0h,j0h,k1h)*d(i0h,j0h)
      udj1k1= h2(i1h,j1h)/h1(i1h,j1h)*aam(i1h,j1h,k1h)*d(i1h,j1h)- &
              h2(i0h,j1h)/h1(i0h,j1h)*aam(i0h,j1h,k1h)*d(i0h,j1h)
      udk1  = udj0k1        + ydist * (udj1k1-udj0k1)

!     interpolate ud in z-direction
      udp   =udk0 +(zdist/abs(0.5*dz(k1h)+0.5*dz(k0h)))*(udk1-udk0)
      udp   =udp/hscnu / (h1pxyz*h2pxyz*dpxyz)
      if(xlocal.le.0.0) udp= -udp


!     interpolate vd at k0 plane
      vdi0k0= h1(i0h,j1h)/h2(i0h,j1h)*aam(i0h,j1h,k0h)*d(i0h,j1h)- &
              h1(i0h,j0h)/h2(i0h,j0h)*aam(i0h,j0h,k0h)*d(i0h,j0h)
      vdi1k0= h1(i1h,j1h)/h2(i1h,j1h)*aam(i1h,j1h,k0h)*d(i1h,j1h)- &
              h1(i1h,j0h)/h2(i1h,j0h)*aam(i1h,j0h,k0h)*d(i1h,j0h)
      vdk0  = vdi0k0        + xdist * (vdi1k0-vdi0k0)

!     interpolate vd at k1 plane
      vdi0k1= h1(i0h,j1h)/h2(i0h,j1h)*aam(i0h,j1h,k1h)*d(i0h,j1h)- &
              h1(i0h,j0h)/h2(i0h,j0h)*aam(i0h,j0h,k1h)*d(i0h,j0h)
      vdi1k1= h1(i1h,j1h)/h2(i1h,j1h)*aam(i1h,j1h,k1h)*d(i1h,j1h)- &
              h1(i1h,j0h)/h2(i1h,j0h)*aam(i1h,j0h,k1h)*d(i1h,j0h)
      vdk1  = vdi0k1        + xdist * (vdi1k1-vdi0k1)

!     interpolate vd in z-direction
      vdp   =vdk0 + (zdist/abs(0.5*dz(k1h)+0.5*dz(k0h)))*(vdk1-vdk0)
      vdp   =vdp/hscnu / (h1pxyz*h2pxyz*dpxyz)
      if(ylocal.le.0.0) vdp= -vdp

!     interpolate wd at k0 plane
      wdp   =(khk1-khk0)/abs(0.5*dz(k1h)+0.5*dz(k0h))
      wdp   =wdp/vscnu / (dpxyz*dpxyz) 
      if(zlocal.le.0.0) wdp= -wdp

      if(abs(up).gt.abs(upmax))upmax=up
      if(abs(vp).gt.abs(vpmax))vpmax=vp
      if(abs(wp).gt.abs(wpmax))wpmax=wp
      if(abs(udp).gt.abs(udpmax))udpmax=udp
      if(abs(vdp).gt.abs(vdpmax))vdpmax=vdp
      if(abs(wdp).gt.abs(wdpmax))wdpmax=wdp

         
!---> advection + random walk + pseudo velocity
!       deltax=deltat * (up/h1pxyz)
!      .   +gasdev(idum)*sqrt(2.*deltat*(aamp/hscnu)/(h1pxyz*h1pxyz))
!      .   +deltat*udp 
!       deltay=deltat*vp/h2pxyz
!      .   +gasdev(idum)*sqrt(2.*deltat*(aamp/hscnu)/(h2pxyz*h2pxyz))
!      .   +deltat*vdp 
!       deltaz=deltat*wp/ dpxyz
!      .   +gasdev(idum)*sqrt(2.*deltat*( khp/vscnu)/( dpxyz* dpxyz))
!      .   +deltat*wdp
      deltax=deltat * up/h1pxyz
      deltay=deltat * vp/h2pxyz
      deltaz=deltat * wp/dpxyz

      if(abs(deltax).gt.abs(deltaxmax))deltaxmax=deltax
      if(abs(deltay).gt.abs(deltaymax))deltaymax=deltay
      if(abs(deltaz).gt.abs(deltazmax))deltazmax=deltaz

!---> update particle location
      xend=xstart + deltax
      yend=ystart + deltay
      zend=zstart + deltaz

!     indices for new cell
      icnew = xend
      jcnew = yend
      do k=1,kbm1
            if(zend.le.z(k).and.zend.ge.z(k+1)) then
                  kcnew=k
                  exit
            endif
      enddo

      if(abs(ic-icnew).gt.1.or.abs(jc-jcnew).gt.1) then
        write(*,*)      iint,ic,icnew,jc,jcnew
        write(*,*)      deltax,deltay,deltaz
      endif

!     check if it is on open boundary
      if(icnew<im-1 .and. icnew>2 .and. jcnew>2 .and. jcnew<jm-1) then
        iout=0
      else
        iout=1
        write(*,*) 'i am out from open boundary'
      endif

!---> reflection from landboundary, surface and bottom
      if(ireflect.eq.0) then
        return
      else
            
!     center coordinates for new cell 
      xcnew = icnew+0.5
      ycnew = jcnew+0.5
      zcnew = zz(kcnew)

!---> new position is on land
      if(fsm(icnew,jcnew).eq.0.0) then

!---> reflection from +x direction
      if(deltax.gt.0.) then
            do loop=ic+1,icnew
                  if(fsm(loop,jc).eq.0.0) then
                        icwall =loop
                        xwall  =icwall
                        xtowall=xwall-xstart
                        xbounce=deltax-xtowall
                        if(xbounce.ge.0.) xend=xwall-xbounce
                  else
                  endif
            enddo
      endif

!---> reflection from -x direction
      if(deltax.lt.0.) then
            do loop=ic-1,icnew,-1
                  if(fsm(loop,jc).eq.0.0) then
                        icwall =loop
                        xwall  =icwall+1
                        xtowall=xstart-xwall
                        xbounce=abs(deltax)-xtowall
                        if(xbounce.ge.0.) xend=xwall+xbounce
                  else
                  endif
            enddo
      endif

!---> reflection from +y direction
      if(deltay.gt.0.) then
            do loop=jc+1,jcnew
                  if(fsm(ic,loop).eq.0.0) then
                        jcwall =loop
                        ywall  =jcwall
                        ytowall=ywall-ystart
                        ybounce=deltay-ytowall
                        if(ybounce.gt.0.) yend=ywall-ybounce
                  else
                  endif
            enddo
      endif

!---> reflection from -y direction
      if(deltay.lt.0.) then
            do loop=jc-1,jcnew,-1
                  if(fsm(ic,loop).eq.0.0) then
                        jcwall =loop
                        ywall  =jcwall+1
                        ytowall=ystart-ywall
                        ybounce=abs(deltay)-ytowall
                        if(ybounce.ge.0.) yend=ywall+ybounce
                  else
                  endif
            enddo
      endif

      endif
!     end of horizontal reflection 


!---> reflection from surface
      if(zend.gt.0)then
      zwall  =0.
      ztowall=zwall-zstart
      zbounce=deltaz-ztowall
      zend   =zwall-zbounce
      endif

!---> reflection from bottom
      if(zend.lt.-1.0)then
      zwall  =-1.0
      ztowall=zstart-zwall
      zbounce=abs(deltaz)-ztowall
      zend   =zwall+zbounce
      endif
   
!---> 'big incremental distance' due to advection and/or
!     diffusion might bounce the particle out of domain. 
!     force particle to stay if still on land.

!     also force it to stay if it on edge of the outmost cell
      if(istay.eq.1)then
            icnewest = xend
            jcnewest = yend
            if(fsm(icnewest,jcnewest).eq.0.0) then
                  xend=xstart
                  yend=ystart
                  zend=zstart
            endif
            if(zend.le.-1.0.or.zend.gt.0.0) then
                  xend=xstart
                  yend=ystart
                  zend=zstart
            endif
      endif
      endif
      return
      end
!***************************************************************************
!
      function gasdev(idum)
      data iset/0/
      if( iset.eq.0 )then
            v1=2.0*ran2(idum)-1.
            v2=2.0*ran2(idum)-1.
            r=v1*v1+v2*v2
            do
                  if (r.ge.1.0) then
                        v1=2.0*ran2(idum)-1.
                        v2=2.0*ran2(idum)-1.
                        r=v1*v1+v2*v2
                        cycle
                  else
                        exit
                  endif
            enddo
            fac=sqrt(-2.0*log(r)/r)
            gset=v1*fac
            gasdev=v2*fac
            iset=1
      else
            gasdev=gset
            iset=0
      endif
      !
      return
      end
!
!**************************************************************************
!
      function ran2(idum)
      parameter(m=714025, ia=1366, ic=150889, rm=1.0/m)
      dimension ir(97)
      data iff/0/
      if(idum.lt.0.or.iff.eq.0)then
            iff=1
            idum=mod(ia-idum,m)
            do j=1,97
                  idum=mod(ia*idum+ic,m)
                  ir(j)=idum
            enddo
            idum=mod(ia*idum+ic,m)
            iy=idum
      endif
      j=1+(97*iy)/m
      if(j.gt.97.or.j.lt.1)then
            print *, 'stop in ran2 ',j
            stop
      endif
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      return
      end
!*************************************************************************
      subroutine archive
      !
      !*************************************************************************
      !     ecomsed model
      !     version 1.3
      !     february 2002
      !*************************************************************************
      !               copyright (c) 2002, hydroqual, inc.                      *
      !                                                                        *
      !  these coded instructions, statements, and computer programs  contain  *
      !  unpublished  proprietary  information of hydroqual, inc., and         *
      !  are protected by federal copyright law.  they  may  not be disclosed  *
      !  to  third  parties  or copied or duplicated in any form, in whole or  *
      !  in part, without the prior written consent of hydroqual, inc.         *
      !                                                                        *
      ! point of contact: ecomsed-support@hydroqual.com                        *
      !*************************************************************************
      use hydro_mod
      use particle_mod
      implicit none

      
      integer     ::    iout
      integer     ::    ic, jc, kc
      integer     ::    i0h, j0h, k0h, i1h, j1h, k1h

      real  ::    xstart, ystart, zstart
      real  ::    xc, yc, zc
      real  ::    xlocal, ylocal, zlocal
      real  ::    xdist, ydist, zdist
      real  ::    dj0, dj1, dp
      real  ::    elj0, elj1, elp
      real  ::    cpi  
      real  ::    xoutp(npartm),youtp(npartm),zoutp(npartm)

      if(iint==irelst .and. layer==1) then
            open(32,file='part_location',status='new', action='write')
            ngradeloop=1
      endif
      
      cpi=acos(-1.)/180.
      do ll=1,nsource
        do nn=1,ngradeloop
          do mm=1,npart

            xoutp(mm)=0.0
            youtp(mm)=0.0
            zoutp(mm)=0.0

            iout=inout(ll,mm,nn)

            xstart=xp(ll,mm,nn)
            ystart=yp(ll,mm,nn)
            zstart=zp(ll,mm,nn)

            ic=xstart
            jc=ystart

            do k=1,kbm1
                  if(zstart.le.z(k).and.zstart.ge.z(k+1)) then
                        kc=k
                        exit
                  else
                        kc=0
                  endif
            enddo
            if (kc==0) then
                  write(*,*) 'trouble !', zstart, iint, kc    
            endif

            xc= ic+0.5
            yc= jc+0.5
            zc= zz(kc)

            !---> interpolate d
            !     find i-range for d
            !
            i0h   = ic
            xlocal= xstart-xc
            xdist = abs(xlocal)
            if(xlocal.gt.0.) i1h=i0h+1
            if(xlocal.le.0.) i1h=i0h-1
            !
            !     use value at one element
            !
            if(fsm(i1h,jc).eq.0.)i1h=i0h
            !
            !     find j-range for d
            !
            j0h   = jc
            ylocal= ystart-yc
            ydist = abs(ylocal)
            if(ylocal.gt.0.) j1h=j0h+1
            if(ylocal.le.0.) j1h=j0h-1
            !
            !     use value at one element
            !
            if(fsm(ic,j1h).eq.0.)j1h=j0h
            !
            if(fsm(i1h,j1h).eq.0) then
                  i1h=i0h
                  j1h=j0h
            endif
            !
            !     interpolate d
            !
            dj0= d(i0h,j0h) + xdist * (d(i1h,j0h)-d(i0h,j0h))
            dj1= d(i0h,j1h) + xdist * (d(i1h,j1h)-d(i0h,j1h))
            dp = dj0        + ydist * (dj1-dj0)
            !
            !     interpolate el
            !
            elj0= el(i0h,j0h) + xdist * (el(i1h,j0h)-el(i0h,j0h))
            elj1= el(i0h,j1h) + xdist * (el(i1h,j1h)-el(i0h,j1h))
            elp = elj0        + ydist * (elj1-elj0)
            !
            xoutp(mm)=xcor(ic,jc)      * (0.5-xlocal)*(0.5-ylocal) &
                      +xcor(ic+1,jc)   * (0.5+xlocal)*(0.5-ylocal) &
                      +xcor(ic+1,jc+1) * (0.5+xlocal)*(0.5+ylocal) &
                      +xcor(ic,jc+1)   * (0.5-xlocal)*(0.5+ylocal)
            !
            youtp(mm)=ycor(ic,jc)      * (0.5-xlocal)*(0.5-ylocal) &
                      +ycor(ic+1,jc)   * (0.5+xlocal)*(0.5-ylocal) &
                      +ycor(ic+1,jc+1) * (0.5+xlocal)*(0.5+ylocal) &
                      +ycor(ic,jc+1)   * (0.5-xlocal)*(0.5+ylocal)
            !
            zoutp(mm)=elp+zstart*dp
            !
          enddo

          write(32,'(4i5,100f15.5)') ll,nn,iout,iint-irelst,zstart,elp,dp, &
              (xoutp(np),np=1,npart),(youtp(np),np=1,npart) &
              ,(zoutp(np),np=1,npart)

        enddo
      enddo
      !write (32,930)ll,nn,npart,time
      !write (32,940) (xoutp(np),np=1,npart)
      !write (32,940) (youtp(np),np=1,npart)
      !write (32,940) (zoutp(np),np=1,npart)
      if(iint==irelend .and. layer==kb) then
		      close(32)
      endif
      return
      end subroutine archive

      subroutine input
      !read hydrodynamic files      !liang 20090727
      use hydro_mod
      implicit none
      character*30 filenow,filenext
      integer     fileend
      !!

        write(*,*) iint

        if(mod(iint,1118)==0) then
            write(filenow,'(I4)') 1000+1118+445
        else
            write(filenow,'(I4)') 1000+mod(iint,1118)+445
        endif
        if(mod(iint,1118)==0) then
            write(filenext,'(I4)') 1000+1118+446
        else
            write(filenext,'(I4)') 1000+mod(iint,1118)+446
        endif
      !!
	  !!
	  open(3570,file='/home/liugl/ouc/JZB_T_SW_TS/H_T&
/hydro/hydro'//trim(filenext)//'.dat',form='unformatted',status='old')
        do
            read(3570,iostat=fileend) i,j,elf(i,j),dt(i,j),&
           (uf(i,j,k),k=1,kb),     (vf(i,j,k),k=1,kb),&
           (w(i,j,k),k=1,kb),      (aam(i,j,k),k=1,kb),&
           (kh(i,j,k),k=1,kb),     (km(i,j,k),k=1,kb)
           ! (s(i,j,k),k=1,kb)
            if(fileend<0) then
                  fileend='true'
                  exit
            endif
        enddo
        close(3570)

      !!  
	  open(3570,file='/home/liugl/ouc/JZB_T_SW_TS/H_T&
/hydro/hydro'//trim(filenow)//'.dat',form='unformatted',status='old')
      !!

        fsm(i,j)=0
		d(i,j)=0.05
        do
            read(3570,iostat=fileend) i,j,el(i,j),d(i,j),&
           (u(i,j,k),k=1,kb),      (v(i,j,k),k=1,kb),&
           (w(i,j,k),k=1,kb),      (aam(i,j,k),k=1,kb),&
           (kh(i,j,k),k=1,kb),     (km(i,j,k),k=1,kb)
           ! (s(i,j,k),k=1,kb)
            if(fileend<0) then
                  fileend='true'
                  exit
            endif      
        enddo
        close(3570)

      !!
        do i=1,im;      do j=1,jm
            el(i,j)=el(i,j)*fsm(i,j)
			if(d(i,j)<=0.05) then
              fsm(i,j)=0
			  d(i,j)=0.05
            else
              fsm(i,j)=1
            endif
        enddo;          enddo
      !!
        do i=1,im
            do j=1,jm
              do k=1,kb
                  u(i,j,k)=u(i,j,k)*fsm(i,j)
                  v(i,j,k)=v(i,j,k)*fsm(i,j)
                  w(i,j,k)=w(i,j,k)*fsm(i,j)
                  aam(i,j,k)=aam(i,j,k)*fsm(i,j)
                  km(i,j,k)=km(i,j,k)*fsm(i,j)
                  kh(i,j,k)=kh(i,j,k)*fsm(i,j)
              enddo
            enddo
        enddo
      !!
      end subroutine input
