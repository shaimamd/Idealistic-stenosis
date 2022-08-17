#define RAD  1.0         !pipe radius
#define ZLENPIPE 54.0   !pipe length
!-----------------------------------------------------------------------
!
!     user subroutines required by nek5000
!
!     Parameters used by this set of subroutines:
!
!-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'NEKUSE'          ! UDIFF, UTRANS

      UDIFF =0.0
      UTRANS=0.0

      return
      end
!-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      
      include 'SIZE'
      include 'NEKUSE'          ! FF[XYZ]
      include 'PARALLEL'

      integer ix,iy,iz,ieg

      ffx = 0.0
      ffy = 0.0
      ffz = 0.0

      return
      end
!-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'NEKUSE'          ! QVOL

      QVOL   = 0.0

      return
      end


!-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'           ! PI
      include 'NEKUSE'          ! UX, UY, UZ, TEMP, Z

!     argument list
      integer ix,iy,iz,ieg

!     local variables
      real C, k, kth, kz,alpha,beta,cons,eps,Ri,rR,xL1,zL1
      real rp,ReTau, th, u_axial, u_th, u_rho, kx, dr
      integer meanProfType

c     Define the initial mean profile
      meanProfType=1;   !1: parabolic velocity profile
                        !2: Reichardt law at a given ReTau

 
c     Parabolic profile ub=1.0
      if (meanProfType.eq.1) then
         r = sqrt(y*y+x*x)
         Ri = 1.0      !pipe radius
         rR = (r*r)/(Ri*Ri)
         cons = 2.0
         uz = cons*(1-rR)
      endif

c     perturb
      eps  = 0
      kz   = 23
      kth  = 13

      th = atan2(y,x)

      xL1=(2.0*PI)
      zL1=ZLENPIPE
      alpha = kz * 2*PI/zL1
      beta  = kth * 2*PI/xL1
      dr = 0.2

      u_axial = eps*beta*sin(alpha*z)*cos(beta*th)
      u_th = -eps*alpha*r*cos(alpha*z)*sin(beta*th)
      if (r>dr) then
         u_rho = eps*alpha*sin(alpha*z)*sin(beta*th)/r
      else
         u_rho = eps*alpha*sin(alpha*z)*sin(beta*th)/dr
      endif

      ux = cos(th)*u_rho - r*sin(th)*u_th
      uy = sin(th)*u_rho + r*cos(th)*u_th
      uz = uz + u_axial

      temp = 0.0

      return
      end
!=======================================================================

      subroutine usrdat
      implicit none
      include 'SIZE'

      return
      end
      
!-----------------------------------------------------------------------
      subroutine userchk
      implicit none
      include 'SIZE'            !
      include 'TSTEP'           ! ISTEP, lastep, time
      include 'INPUT'           ! IF3D, PARAM
      
!     start framework
      if (ISTEP.eq.0) call frame_start

!     monitor simulation
      call frame_monitor

!     post processing step
      call pstat3d_main()
      
!     finish simulation
      LASTEP=1

!     finalise framework
      if (ISTEP.eq.NSTEPS.or.LASTEP.eq.1) then
         call frame_end
      endif
      
      call exitt0
     
      return
      end
!--------------------------------------------------
      subroutine userbc(ix,iy,iz,iside,ieg)
      implicit none
      integer ix,iy,iz,iside,eg,e,ieg,i
      real U0,delta,sn(3),S0,vn,cons,rR,Ri
      
      include 'SIZE'
      include 'SOLN'
      include 'GEOM'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'NEKNEK'

      e=gllel(ieg)
      pa=0.0
      
      if(cbu.eq.'o  ') then    !Dong BC
         U0=1.0
         delta=1.0
         ux=vx(ix,iy,iz,e)
         uy=vy(ix,iy,iz,e)
         uz=vz(ix,iy,iz,e)
         call getSnormal(sn,ix,iy,iz,iside,e)
         vn=ux*sn(1)+uy*sn(2)+uz*sn(3)
         S0=0.5*(1.0-tanh(uz/U0/delta))
         pa=-0.5*(ux*ux+uy*uy+uz*uz)*S0
      elseif(cbu.eq.'v  ') then
         ux=0.0
         uy=0.0
          r = sqrt(y*y+x*x)
         Ri = 1.0      !pipe radius
         rR = (r*r)/(Ri*Ri)
         cons = 2.0
         uz = cons*(1-rR)
      elseif(cbu.eq.'W  ') then
         ux=0.0
         uy=0.0
         uz=0.0        
      endif
      temp=0.0

      return
      end
!--------------------------------------------------
      subroutine usrdat2
      !implicit none
      include 'SIZE'
      include 'TOTAL'

      do iel=1,nelv
         do ifc=1,2*ndim
            id_face=bc(5,ifc,iel,1)
            if(id_face.eq.1)then
              cbc(ifc,iel,1)='v  '
            elseif(id_face.eq.2)then
              cbc(ifc,iel,1)='o  '   !'O' sponge and 'o' dong
            elseif(id_face.eq.3) then
              cbc(ifc,iel,1)='W  '
            endif
          enddo
      enddo

      return
      end

!--------------------------------------------------
      subroutine usrdat3
      implicit none
      include 'SIZE'
      include 'INPUT'

      return
      end
!======================================================================


!======================================================================
!> @brief Register user specified modules
      subroutine frame_usr_register
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------
!     register modules
      call io_register()
      call pstat3d_register()

      return
      end subroutine
!======================================================================
!> @brief Initialise user specified modules
      subroutine frame_usr_init
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------
!     initialise modules
      call pstat3d_init()

      return
      end subroutine
!======================================================================
!> @brief Finalise user specified modules
      subroutine frame_usr_end
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------
!     finalise modules
      
      return
      end subroutine
!======================================================================

c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
