      subroutine functionLap(Nx,gamma,gLaplace,xgradeC)

      implicit none
      integer Nx,i

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision gamma(Nx)
      double precision gLaplace(Nx),xgradeC(Nx)
      double precision gammaim2,gammaim1,gammaip1
      double precision gLapX,gLapY,thetai,thetaim1,psii,psiim1





       do i=1,Nx
!       No-Flux boundary condition
       if(i .eq. 1) then
        gammaim2=gamma(i+1)
        gammaim1=gamma(i+1)
        gammaip1=gamma(i+1)

!        gammaim2=0.070
!        gammaim1=0.070

       elseif(i .eq. 2) then
        gammaim2=-gamma(i)+2*gamma(i-1)
        gammaim1=gamma(i-1)
        gammaip1=gamma(i+1)

!        gammaim2=0.070

       elseif(i .eq. Nx) then
        gammaim2=gamma(i-2)
        gammaim1=gamma(i-1)
        gammaip1=gamma(i-1)
       else
        gammaim2=gamma(i-2)
        gammaim1=gamma(i-1)
        gammaip1=gamma(i+1)

       endif

!      Periodic Boundary
!       if(i .eq. 1) then
!        gammaim2=gamma(Nx-1)
!        gammaim1=gamma(Nx)
!        gammaip1=gamma(2)
!
!       elseif(i .eq. 2) then
!        gammaim2=gamma(Nx)
!        gammaim1=gamma(i-1)
!        gammaip1=gamma(i+1)
!
!
!       elseif(i .eq. Nx) then
!        gammaim2=gamma(i-2)
!        gammaim1=gamma(i-1)
!        gammaip1=gamma(1)
!       else
!        gammaim2=gamma(i-2)
!        gammaim1=gamma(i-1)
!        gammaip1=gamma(i+1)
!
!       endif


!     Non Flux Boundary Cond
!       if(Ny .eq. 1) then
!        gammajm1=gamma(i,j)
!        gammajp1=gamma(i,j)
!       elseif(j .eq. 1) then
!        gammajm1=gamma(i,j+1)
!        gammajp1=gamma(i,j+1)
!       elseif(j .eq. Ny) then
!        gammajp1=gamma(i,j-1)
!        gammajm1=gamma(i,j-1)
!       else
!        gammajp1=gamma(i,j+1)
!        gammajm1=gamma(i,j-1)
!       endif





        gLaplace(i)=(gammaip1+gammaim1-2*gamma(i))/(dx**2)


!
        if(gammaip1 .eq. gamma(i)) then
        thetai=1.d-10
        else
        thetai=(gamma(i)-gammaim1)/(gammaip1-gamma(i))+1.d-10
        endif

        if(gamma(i) .eq. gammaim1)then
        thetaim1=1.d-10
        else
        thetaim1=(gammaim1-gammaim2)/(gamma(i)-gammaim1)+1.d-10
        endif

        psii=max(0.0,min(1.0,1.0/3.0+thetai/6.0,thetai))
        psiim1=max(0.0,min(1.0,1.0/3.0+thetaim1/6.0,thetaim1))

      xgradeC(i)=(1.0-psiim1+psii/thetai)*(-gammaim1+gamma(i))/(dx)


!        xgradeC(i)=(gammaip1-gammaim1)/(2*dx)
!        ygradeC(i,j)=(gammajp1-gamma(i,j))/(dy)

      enddo

      return
      end
