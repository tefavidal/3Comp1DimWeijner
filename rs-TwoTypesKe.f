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
     .                     -vdx(i)*xgradeC(i)
        else
            gammaprime(i)=-(keU/dke0)/depsilon*gamma(i)
     .                      +depsilon*gLaplace(i)
     .                     -vdx(i)*xgradeC(i)
        endif
      enddo

      return

      end
!      ***********************************************************
!     *************************************************



