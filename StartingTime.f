!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FromCellToGrid(Nx,Nc,cells,grid)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      integer Nx, Nc,i, j,k
      double precision cells(Nc)
      double precision d
      integer grid(Nx)


        d=1.0

      do i=1,Nx
         grid(i)=0
      enddo

      do k=1,Nc
        i=ceiling(cells(k)/dx)
        if(i .gt.Nx)then
            cells(k)=cells(k)-Nx*dx
            i=ceiling(cells(k)/dx)
        endif

        if(i .lt.1)then
            cells(k)=cells(k)+Nx*dx
            i=ceiling(cells(k)/dx)
        endif

        grid(i)=grid(i)+1
      enddo

      return
      end

!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine initialDiscreteDistribution(Nx,Nc,cells)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      integer Nx, Nc,i,k
      double precision cells(Nc), percentage, d
      real aux
      integer grid(Nx), auxInt

  ! ----- variables for portable seed setting -----
      INTEGER :: i_seed
      INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
      INTEGER, DIMENSION(1:8) :: dt_seed
  ! ----- end of variables for seed setting -----

 ! ----- Set up random seed portably -----
      CALL RANDOM_SEED(size=i_seed)
      ALLOCATE(a_seed(1:i_seed))
      CALL RANDOM_SEED(get=a_seed)
      CALL DATE_AND_TIME(values=dt_seed)
      a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)
     .*dt_seed(6)
       write(6,*) 'seed=',a_seed(i_seed)
      CALL RANDOM_SEED(put=a_seed)
      DEALLOCATE(a_seed)
  ! ----- Done setting up random seed ----

        d=1.0;
      do i=1,Nx
         grid(i)=0
      enddo
            k=1;
            percentage=real(Nc)/(real(Nx))
      do i=1,Nx
         call random_number(aux)
         if(aux .lt. percentage .and. grid(i) .ge. 0.0)then
                if (k .gt. Nc)then
                    exit
                endif
                cells(k)=(i-0.5)*dx
                k=k+1
         endif
      enddo

      auxInt=k-1

      do i=k,Nc
        cells(i)=cells(1)
      enddo

      call FromCellToGrid(Nx,Nc,cells,grid)


      do while (k .le. Nc)
        call random_number(aux)
        i=ceiling(aux*Nx)

        if(grid(i) .eq. 0)then
            cells(k)=(i-0.2)*dx
            k=k+1
            grid(i)=1
        endif
      enddo




      return
      end
