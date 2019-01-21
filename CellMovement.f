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

      do k= 1,Nc
        i=ceiling(cells(k)/dx)
        grad=abs(xgradeC(i))
        if(rhogrid(i) .ge. 0.6 .and. grad .ge. 1.0
     .   .and. t/dk1 .gt. 50)then
!          write(6,*) 'Atempt move'
            angle=xgradeC(i)/grad
            newx=cells(k)+h*speed*angle
            oldi=ceiling(cells(k)/dx)
            i=ceiling(newx/dx)
            if (oldi .eq. i )then
                cells(k)=newx
            else

                    if(i .gt.Nx)then
                        newx=newx-Nx*dx
                        i=ceiling(newx/dx)
                    endif

                    if(i .lt.1)then
                        newx=newx+Nx*dx
                        i=ceiling(newx/dx)
                    endif

!                    if(grid(i) .eq. 0) then
                    if(grid(i) .lt. 5) then
!                  write(6,*) 'Actual move'
                        cells(k)=newx
                        grid(i)=grid(i)+1
                        grid(oldi)=grid(oldi)-1
                    endif
            endif


        endif
      enddo

      return
      end

!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine ScalatedMovement(t,h,Nx,Nc, gamma, ro, beta,cells)

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
      integer MaxPerGrid

      call GridifyingBetaRho(Nx,Nc,cells,grid,beta,ro,betagrid
     . ,rhogrid)

      call functionLap(Nx,gamma,gLaplace,xgradeC)

      call FromCellToGrid(Nx,Nc,cells,grid)


!      speed=0.02/sqrt(dke0*0.024)
      speed=0.0;
      MaxPerGrid=1
      do k= 1,Nc
         i=ceiling(cells(k)/dx)
         grad=abs(xgradeC(i))
         if(rhogrid(i) .ge. 0.6 .and. grad .ge. 1.0
     .   .and. t/dk1 .gt. 50)then
            if(t/dk1 .gt.100)then
               MaxPerGrid=2
            endif
!            write(6,*) 'Atempt move'
            angle=xgradeC(i)/grad
            newx=cells(k)+h*speed*angle
            oldi=ceiling(cells(k)/dx)
            i=ceiling(newx/dx)

            if (oldi .eq. i )then
                cells(k)=newx

            else
               if(i .gt.Nx)then
                  newx=newx-Nx*dx
                  i=ceiling(newx/dx)
               endif

               if(i .lt.1)then
                  newx=newx+Nx*dx
                  i=ceiling(newx/dx)
               endif

!               if(grid(i) .eq. 0) then
               if(grid(i) .lt. MaxPerGrid .and. grid(i) .ge. 0) then
                  cells(k)=newx
                  grid(i)=grid(i)+1
                  grid(oldi)=grid(oldi)-1
               endif
            endif


        endif
      enddo

      return
      end
