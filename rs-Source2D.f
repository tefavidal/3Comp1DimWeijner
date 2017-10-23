      subroutine rs(t,Nx,beta,gamma,ro,betaprime,gammaprime,roprime
     . ,grid,vdx)

      implicit none
      double precision t, factor, aux
      integer Nx, i, j,ID
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision beta(Nx),gamma(Nx),ro(Nx)
      double precision betaprime(Nx),gammaprime(Nx),roprime(Nx)
      double precision f1,f2,Phi,Y,r1,r2,noise
      double precision gLaplace(Nx)
      double precision xgradeC(Nx),ygradeC(Nx)
      double precision vdx(Nx),vdy, percentage, factor2
      double precision dke(Nx),dsigma(Nx)
      integer grid(Nx), Nc
!  ! ----- variables for portable seed setting -----
!      INTEGER :: i_seed
!      INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
!      INTEGER, DIMENSION(1:8) :: dt_seed
!  ! ----- end of variables for seed setting -----
!
! ! ----- Set up random seed portably -----
!      CALL RANDOM_SEED(size=i_seed)
!      ALLOCATE(a_seed(1:i_seed))
!      CALL RANDOM_SEED(get=a_seed)
!      CALL DATE_AND_TIME(values=dt_seed)
!      a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)
!     .*dt_seed(6)
!!       write(6,*) 'seed=',a_seed(i_seed)
!      CALL RANDOM_SEED(put=a_seed)
!      DEALLOCATE(a_seed)
!  ! ----- Done setting up random seed ----

        Nc=0
      do i=1,Nx
         if(grid(i) .gt. 0.5)then
            Nc=Nc+1
         endif
      enddo

      call functionLap(Nx,gamma,gLaplace,xgradeC)

       percentage=real(Nc)/(real(Nx))

       do i=1,Nx
!       Extra variables calculation
        if(grid(i) .gt. 0.5)then
            vdy=0.0
            aux=gamma(i)
            dke(i)=1.0
            dsigma(i)=1.0
            f1=(1.d0+dk*aux)/(1.d0+aux)
            f2=(dL1+dk*dL2*dc*aux)/(1.d0+dc*aux)
            Y=ro(i)*aux/(1.d0+aux)
            Phi=(dlambda1+Y*Y)/(dlambda2+Y*Y)
!            call random_number(r1)
!            call random_number(r2)
!            noise=sqrt(-2*log(r1))*cos(2*Pi*r2);
!
!           Right hand side
            betaprime(i)=(s1*Phi*dsigma(i)-beta(i))
     .                        /depsilonp
            roprime(i)=(-f1*ro(i)+f2*(1.d0-ro(i)))
            gammaprime(i)=1.0/depsilon*s2*beta(i)
     .                 -1.0/depsilon*gamma(i)
     .                      +depsilon*gLaplace(i)
     .              -  (vdx(i)*xgradeC(i))
!     .              + 10*noise
        else
!       Right hand side
            betaprime(i)=0.0
            roprime(i)=0.0
            gammaprime(i)=-1.0/depsilon*gamma(i)
     .                      +depsilon*gLaplace(i)
     .              -  (vdx(i)*xgradeC(i))
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
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

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

!        gammaim2=0.0
!        gammaim1=0.0

       elseif(i .eq. 2) then
        gammaim2=-gamma(i)+2*gamma(i-1)
        gammaim1=gamma(i-1)
        gammaip1=gamma(i+1)

!        gammaim2=0.0

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
      subroutine UpdateDiscreteCells(h,Nx,Nc, gamma, ro, cells)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      integer Nx, Ny,Nc, i, j, k, oldi, oldj
      double precision cells(Nc)
      double precision gamma(Nx), ro(Nx)
      double precision gLaplace(Nx),xgradeC(Nx)
      double precision grad, angle, h,speed, newx, newy
      integer grid(Nx)


      call functionLap(Nx,gamma,gLaplace,xgradeC)

      call FromCellToGrid(Nx,Nc,cells,grid)


!        speed=0.05/sqrt(dke0*0.024)
        speed=0.0;

      do k= 1,Nc
        i=ceiling(cells(k)/dx)
        grad=abs(xgradeC(i))
        if(ro(i) .ge. 0.6 .and. grad .ge. 1.0)then
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

                    if(grid(i) .eq. 0) then
!                  write(6,*) 'Actual move'
                        cells(k)=newx
                        grid(i)=1
                        grid(oldi)=0
                    endif
            endif


        endif
      enddo

      return
      end
