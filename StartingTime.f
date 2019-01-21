!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GammaToList(Nx,Nc,cells,gamma,gammalist)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      integer Nx, Nc,i, j,k
      double precision cells(Nc)
      double precision d
      double precision gamma(Nx), gammalist(Nc)


      do k=1,Nc
        gammalist(k)=0.0
        if(cells(k) .lt. 0)then
           cycle
        endif
        i=ceiling(cells(k)/dx)
        if(i .gt.Nx)then
            cells(k)=cells(k)-Nx*dx
            i=ceiling(cells(k)/dx)
        endif

        if(i .lt.1)then
            cells(k)=cells(k)+Nx*dx
            i=ceiling(cells(k)/dx)
        endif
        gammalist(k)=gamma(i)
      enddo

      return
      end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GridifyingBetaRho(Nx,Nc,cells,grid,beta,ro,betagrid
     . ,rhogrid)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      integer Nx, Nc,i, j,k
      double precision cells(Nc)
      double precision beta(Nc), ro(Nc)
      double precision betagrid(Nx), rhogrid(Nx)
      double precision d
      integer grid(Nx)


        d=1.0

      do i=1,Nx
         grid(i)=0
         betagrid(i)=0
         rhogrid(i)=0
      enddo

      do k=1,Nc
        if(cells(k) .lt. 0)then
           cycle
        endif
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
        betagrid(i)=beta(k)+betagrid(i)
        rhogrid(i)=ro(k)+rhogrid(i)
      enddo
!     Setting beta and rho as the average of the cells in a spot

      do i=1,Nx
         if (grid(i).gt.0)then
            rhogrid(i)=rhogrid(i)/grid(i)
         endif
      enddo



      return
      end
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FromCellToGrid(Nx,Nc,cells,grid)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      integer Nx, Nc,i, j,k
      double precision cells(Nc)
      double precision d
      integer grid(Nx)


        d=1.0

      do i=1,Nx
         grid(i)=0
      enddo

      do k=1,Nc
        if(cells(k) .lt. 0)then
           cycle
        endif
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
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

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
!      a_seed(i_seed)=dt_seed(8)
      a_seed(i_seed)=dt_seed(8)+dt_seed(7)*1000
      a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
       write(6,*) 'seed size=',i_seed
       write(6,*) 'date seed=',dt_seed
       write(6,*) 'seed=',a_seed(1:i_seed)
!       write(6,*) 'seed=',a_seed(i_seed)
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


!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine NonRandomDistribution(Nx,Nc,cells,ClusterSize)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      integer Nx, Nc,i,k
      double precision cells(Nc), percentage, d, ClusterSize
      real aux
      integer grid(Nx), auxInt


        d=1.0;
      do i=1,Nx
         grid(i)=0
      enddo
            k=1;
      do i=1,Nx
         if(i .gt. (Nx-ClusterSize)/2)then
                if (k .gt. ClusterSize)then
                    exit
                endif
                cells(k)=(i-0.5)*dx
                k=k+1
         endif
      enddo


      do i=k,Nc
        cells(i)=-1.0
      enddo

      return
      end


!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine DispersionRelationTest(Nx,Nc,cells)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      integer Nx, Nc,i,k
      double precision cells(Nc), percentage, d
      real aux
      integer grid(Nx), auxInt


        d=1.0;
      do i=1,Nx
         grid(i)=0
      enddo
            k=1;
      do i=1,Nx
!%%%%%%  10%
!         if(mod(i,10) .eq. 2)then
!
!%%%%%%  20%
!         if(mod(i,5) .eq. 2)then
!
!%%%%%%% 30 %
!         if(mod(i,10) .eq. 3 .or. mod(i,10) .eq. 6
!     .    .or. mod(i,10) .eq. 9)then

!%%%%%%  40%
         if(mod(i,5) .eq. 2 .or. mod(i,5) .eq. 4)then

!%%%%%%%% 50%
!         if(mod(i,2) .eq. 1)then

!%%%%%%  60%
!      if(mod(i,5) .eq. 1 .or. mod(i,5) .eq. 3 .or. mod(i,5)
!     .  .eq. 0)then

!%%%%%%% 70 %
!         if(mod(i,10) .ne. 3 .and. mod(i,10) .ne. 6
!     .    .and. mod(i,10) .ne. 9)then
                if (k .gt. Nc)then
                    exit
                endif
                cells(k)=(i-0.5)*dx
                k=k+1
         endif
      enddo

!      do i=201,210
!                cells(k)=(i-0.5)*dx
!                k=k+1
!      enddo
!
!      do i=211,Nx
!!         if(mod(i,5) .eq. 2 .or. mod(i,5) .eq. 4)then
!         if(mod(i,2) .eq. 1)then
!                if (k .gt. Nc)then
!                    exit
!                endif
!                cells(k)=(i-0.5)*dx
!                k=k+1
!         endif
!      enddo



      return
      end

!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine DispersionRelationRandom(Nx,Nc,cells)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

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
      do i=1,200
         call random_number(aux)
         if(aux .lt. 0.4 .and. grid(i) .ge. 0.0)then
                if (k .gt. Nc)then
                    exit
                endif
                cells(k)=(i-0.5)*dx
                k=k+1
         endif
      enddo

      do i=201,300
                cells(k)=(i-0.5)*dx
                k=k+1
      enddo

      do i=301,Nx
         call random_number(aux)
         if(aux .lt. 0.4 .and. grid(i) .ge. 0.0)then
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


