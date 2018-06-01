      subroutine out(t,Nx,Nc,gamma,ro,beta,cells)

      implicit none

      double precision t,dls
      integer Nx, Nc,i

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx),ro(Nc),beta(Nc)
      doubleprecision betagrid(Nx), rhogrid(Nx)
      double precision cells(Nc)
      integer grid (Nx)

      call GridifyingBetaRho(Nx,Nc,cells,grid,beta,ro,betagrid
     . ,rhogrid)

      dls=dk1/(dke0*Diffgamma)**0.5

      do i=1,Nx
            write(10,*) t/dk1,i*dx/dls,gamma(i),rhogrid(i),
     .       betagrid(i),grid(i)
         enddo

      write(10,*)
      return
      end





