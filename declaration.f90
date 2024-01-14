      module param_mod
      implicit none
      
      !reflection from landboundary, surface and bottom
      integer, parameter :: ireflect = 1
      
      ! 'big incremental distance' due to advection and/or
      ! diffusion might bounce the particle out of domain.
      ! force particle to stay if still on land.
      ! also force it to stay if it on edge of the outmost cell
      integer, parameter :: istay  = 1
      
      end module param_mod
      
      module particle_mod
      !  for particle tracking
      !
      !  nsrcm  = maximum no. of particle sources
      !  npartm = maximum no. of particles/release
      !  nconvm = maximum no. of particle releases before conversion to
      !           concentration
      !
      !  note:  max. no. of particles = npartm * nconvm
      !
      integer, parameter :: nsrcm=40000,npartm=1,nconvm=10
      
      integer :: ll, mm, nn, np
      integer :: nfreq,npart,nconv,irelst,irelend,npclass,nsource,ngradeloop      
      
      real :: isource(nsrcm),jsource(nsrcm),ksource(nsrcm)
      
      real :: xp(nsrcm,npartm,nconvm), yp(nsrcm,npartm,nconvm), zp(nsrcm,npartm,nconvm)
      
      real :: inout(nsrcm,npartm,nconvm)
      
      
      end module particle_mod
      
      
      module hydro_mod
      
      integer,  parameter :: im=312, jm=213, kb=17
      integer, parameter :: kbm1=kb-1, kbm2=kb-2
      integer, parameter :: numebc=459

      real, parameter :: dti=40

      real, parameter :: hprnu=0.100e+01
      !horizontal prandtl number - ratio of horizontal
      !viscosity to horizontal diffusivity (momentum
      !mixing/dispersive mixing)
      
      real, parameter :: vprnu=0.100e+01
      !vertical prandtl number - ratio of vertical viscosity
      !to vertical diffusivity (momentum
      !mixing/diffusive mixing)
      
      
      
      integer :: iint
      integer :: i,j,k,layer

      real, allocatable, dimension(:) :: ieta, jeta
      !---------------- 1-d arrays -------------------------------------------
      real, dimension(kb) :: z, zz, dz, dzz
      
      !---------------- 2-d arrays -------------------------------------------
      
      real :: h(im,jm),h1(im,jm),h2(im,jm),d(im,jm),dt(im,jm),datum(im,jm),ang(im,jm),cor(im,jm)
      integer :: dum(im,jm),dvm(im,jm),fsm(im,jm)
      real :: xgrid(im,jm),ygrid(im,jm)
      real :: xcor(im+1,jm+1),ycor(im+1,jm+1)
      real :: art(im,jm),aru(im,jm),arv(im,jm),cbc(im,jm)
      real :: wusurf(im,jm),wvsurf(im,jm),wubot(im,jm),wvbot(im,jm)
      real :: wtsurf(im,jm),wssurf(im,jm),tps(im,jm)
      
      real :: elf(im,jm),el(im,jm),elb(im,jm)
      real :: uaf(im,jm),ua(im,jm),uab(im,jm),vaf(im,jm),va(im,jm),vab(im,jm)
      
      !
      !---------------- 3-d arrays -------------------------------------------
      !
      real :: uf(im,jm,kb), u(im,jm,kb), ub(im,jm,kb)
      real :: vf(im,jm,kb), v(im,jm,kb), vb(im,jm,kb)
      real :: w(im,jm,kb)
      
      real :: km(im,jm,kb), kh(im,jm,kb), aam(im,jm,kb)
      
      real :: t(im,jm,kb), tb(im,jm,kb)
      real :: s(im,jm,kb), sb(im,jm,kb),  rho(im,jm,kb), rmean(im,jm,kb)
      
      end module hydro_mod
       
