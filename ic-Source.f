      subroutine ic(t,Nx,Nc,beta,gamma,ro,cells)
      
      implicit none

      double precision t
      integer Nx,Nc

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      double precision beta(Nc),gamma(Nx),ro(Nc)
      double precision gamma0(10),beta0(10),ro0(10)
      double precision dke(Nx),dsigma(Nx)
      double precision cells(Nc)
      integer nfix, i, j
      integer grid(Nx)
      integer factor

!     %%%%% Initial State for ke=9.5 sigma=0.55
!         nfix=1
!         gamma0(1)=0.006
!         beta0(1)=0.3167
!         ro0(1)=0.8953

!     %%%%% Initial State for ke=7.0 sigma=0.55
         nfix=1
         gamma0(1)=0.359
         beta0(1)=13.91
         ro0(1)=0.285

!     %%%%% Initial High state
!         nfix=1
!         gamma0(1)=0.5
!         beta0(1)=0.3167
!         ro0(1)=0.8953

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do i=1,nfix
      write(6,*)  'gamma=',gamma0(i),'  beta=',beta0(i),'  ro=',ro0(i)
      enddo

      call FromCellToGrid(Nx,Nc,cells,grid)

      do i=1,Nx
           if(grid(i) .gt. 0.5)then
            gamma(i)=gamma0(1)
            else
!            gamma(i)=0.0
            gamma(i)=gamma0(1)
            endif
      enddo

      do i=1,Nc
            ro(i)=ro0(1)
            beta(i)=beta0(1)
      enddo



         gamma01=gamma0(1)
         beta01=beta0(1)
         ro01=ro0(1)

      return

      end

