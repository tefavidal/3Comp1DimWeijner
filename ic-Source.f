      subroutine ic(t,Nx,Nc,beta,gamma,ro,cells)
      
      implicit none

      double precision t
      integer Nx,Nc

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      double precision beta(Nc),gamma(Nx),ro(Nc)
      double precision gamma0(10),beta0(10),ro0(10)
      double precision dke(Nx),dsigma(Nx)
      double precision cells(Nc),noise,r1,r2
      integer nfix, i, j
      integer grid(Nx)
      integer factor

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
!     %%%%% Initial State for ke=9.5 sigma=0.55
!         nfix=1
!         gamma0(1)=0.006
!         beta0(1)=0.3167
!         ro0(1)=0.8953

!     %%%%% Initial State for ke=7.0 sigma=0.55
!         nfix=1
         gamma0(1)=0.359
         beta0(1)=13.91
         ro0(1)=0.285

!     %%%%% Initial State for ke=5.0 sigma=0.55
!         nfix=1
!         gamma0(1)=0.4455
!         beta0(1)=112.5272
!         ro0(1)=0.2299

!     %%%%% Initial High state
!         nfix=1
!         gamma0(1)=0.5
!         beta0(1)=0.3167
!         ro0(1)=0.8953

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      call FromCellToGrid(Nx,Nc,cells,grid)

      do i=1,Nx
           if(grid(i) .gt. 0.5)then
             call random_number(r1)
             call random_number(r2)
             noise=sqrt(-2*log(r1))*cos(2*Pi*r2);
            gamma(i)=gamma0(1)+0.05*noise
            else
            gamma(i)=0.0
!            gamma(i)=gamma0(1)
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

