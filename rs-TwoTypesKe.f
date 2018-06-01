      subroutine rs(t,Nx,Nc,beta,gamma,ro,betaprime,gammaprime,roprime
     . ,cells,vdx)

      implicit none
      double precision t, factor, aux
      integer Nx, i, j,ID
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision beta(Nc),gamma(Nx),ro(Nc)
      double precision gammalist(Nc)
      double precision betaprime(Nc),gammaprime(Nx),roprime(Nc)
      double precision betagrid(Nx), rhogrid(Nx)
      double precision f1,f2,Phi,Y,r1,r2,noise
!      double precision keB, keU
      double precision gLaplace(Nx)
      double precision xgradeC(Nx),ygradeC(Nx)
      double precision vdx(Nx),vdy, percentage
      double precision dke(Nx),dsigma(Nx), cells(Nx)
      integer grid(Nx), Nc

      call functionLap(Nx,gamma,gLaplace,xgradeC)
      call GridifyingBetaRho(Nx,Nc,cells,grid,beta,ro,betagrid
     . ,rhogrid)
      call GammaToList(Nx,Nc,cells,gamma,gammalist)

       percentage=real(Nc)/(real(Nx))

!      keU=2.0
!      keB=6.0
       do i=1,Nc
!       Extra variables calculation
         if(cells(i) .lt. 0) then
           betaprime(i)=0
           roprime(i)=0
         endif
            vdy=0.0
            aux=gammalist(i)
            f1=(1.d0+dk*aux)/(1.d0+aux)
            f2=(dL1+dk*dL2*dc*aux)/(1.d0+dc*aux)
            Y=ro(i)*aux/(1.d0+aux)
            Phi=(dlambda1+Y*Y)/(dlambda2+Y*Y)

!           Right hand side
            betaprime(i)=(s1*Phi-beta(i))
     .                        /depsilonp
            roprime(i)=(-f1*ro(i)+f2*(1.d0-ro(i)))
       enddo

       do i=1,Nx
!       Right hand side
        if(grid(i) .gt. 0.5)then
            gammaprime(i)=1.0/depsilon*s2*betagrid(i)
     .                 -(keB*grid(i)+keU)/dke0/depsilon*gamma(i)
     .                      +depsilon*gLaplace(i)
        else
            gammaprime(i)=-(keU/dke0)/depsilon*gamma(i)
     .                      +depsilon*gLaplace(i)
        endif
      enddo

      return

      end
!      ***********************************************************
!     *************************************************
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
!        if(gammaip1 .eq. gamma(i,j)) then
!        thetai=1.d-10
!        else
!        thetai=(gamma(i,j)-gammaim1)/(gammaip1-gamma(i,j))+1.d-10
!        endif
!
!        if(gamma(i,j) .eq. gammaim1)then
!        thetaim1=1.d-10
!        else
!        thetaim1=(gammaim1-gammaim2)/(gamma(i,j)-gammaim1)+1.d-10
!        endif
!
!        psii=max(0.0,min(1.0,1.0/3.0+thetai/6.0,thetai))
!        psiim1=max(0.0,min(1.0,1.0/3.0+thetaim1/6.0,thetaim1))
!
!      xgradeC(i,j)=(1.0-psiim1+psii/thetai)*(-gammaim1+gamma(i,j))/(dx)


        xgradeC(i)=(gammaip1-gammaim1)/(2*dx)
!        ygradeC(i,j)=(gammajp1-gamma(i,j))/(dy)

      enddo

      return
      end



!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine UpdateDiscreteCells(t,h,Nx,Nc, gamma, ro, beta,cells)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      integer Nx, Ny,Nc, i, j, k, oldi, oldj
      double precision cells(Nc), ro(Nc), beta(Nc)
      double precision gamma(Nx), betagrid(Nx), rhogrid(Nx)
      double precision gLaplace(Nx),xgradeC(Nx)
      double precision grad, angle, h,speed, newx, newy,t
      integer grid(Nx)

      call GridifyingBetaRho(Nx,Nc,cells,grid,beta,ro,betagrid
     . ,rhogrid)

      call functionLap(Nx,gamma,gLaplace,xgradeC)

      call FromCellToGrid(Nx,Nc,cells,grid)


!        speed=0.02/sqrt(dke0*0.024)
        speed=0.0;

!      do k= 1,Nc
!        i=ceiling(cells(k)/dx)
!        grad=abs(xgradeC(i))
!        if(ro(k) .ge. 0.6 .and. grad .ge. 1.0
!     .   .and. t/dk1 .gt. 50)then
!!          write(6,*) 'Atempt move'
!            angle=xgradeC(i)/grad
!            newx=cells(k)+h*speed*angle
!            oldi=ceiling(cells(k)/dx)
!            i=ceiling(newx/dx)
!            if (oldi .eq. i )then
!                cells(k)=newx
!            else
!
!                    if(i .gt.Nx)then
!                        newx=newx-Nx*dx
!                        i=ceiling(newx/dx)
!                    endif
!
!                    if(i .lt.1)then
!                        newx=newx+Nx*dx
!                        i=ceiling(newx/dx)
!                    endif
!
!!                    if(grid(i) .eq. 0) then
!                    if(grid(i) .lt. 5) then
!!                  write(6,*) 'Actual move'
!                        cells(k)=newx
!                        grid(i)=grid(i)+1
!                        grid(oldi)=grid(oldi)-1
!                    endif
!            endif
!
!
!        endif
!      enddo

      return
      end
