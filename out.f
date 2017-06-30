      subroutine out(t,Nx,Nc,gamma,ro,beta,cells)

      implicit none

      double precision t,dls
      integer Nx, Nc,i

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx),ro(Nx),beta(Nx)
      double precision cells(Nc)
      integer grid (Nx)

      call FromCellToGrid(Nx,Nc,cells,grid)

      dls=dk1/(dke0*Diffgamma)**0.5

      do i=1,Nx
            write(10,*) t/dk1,i*dx/dls,gamma(i),ro(i),beta(i),grid(i)
         enddo

      write(10,*)
      return
      end


!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine outFinal(t,Nx,Nc,gamma,ro,beta,cells)

      implicit none
      integer Nx, Nc, i
      double precision t

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision gamma(Nx),ro(Nx),beta(Nx)
      double precision cells(Nc)



      do i=1,Nx
            write(42,*) t/dk1,gamma(i),ro(i),beta(i)
      enddo
      write(42,*)


      do i=1,Nc
        write(43,*) cells(i)
      enddo
      write(43,*)


      return
      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine loadState(t,Nx,Nc,gamma,ro,beta,cells)

      implicit none

      integer Nx, Nc, i
      double precision t, aux

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      double precision gamma(Nx),ro(Nx), beta(Nx)
      double precision cells(Nc)

      do i=1,Nx
            read(7,*) aux,gamma(i),ro(i),beta(i)
      enddo
      close(7)
      t=aux*dk1

      do i=1,Nc
        read(8,*) cells(i)
      enddo
      close(8)


      return
      end

