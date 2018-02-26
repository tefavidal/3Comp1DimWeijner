      program main

      implicit none


      INTEGER, PARAMETER :: Nx=2000
      INTEGER, PARAMETER :: Nc=20
!     Nc decided as Nx*percentageofcells

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision beta(Nx),gamma(Nx),ro(Nx)
      double precision t
      double precision cells(Nc)
      integer pgbeg, i, counter, factor
      integer grid(Nx)
      character(len=4) ct1
      character(len=52) ct2
!      character(len=22) ct2
      character(len=56) ct3

      double precision gamma0(10),ro0(10),beta0(10)

      t=0.d0
      counter=0;
      call anfang(t,Nx,Nc,beta,gamma,ro,cells)
      ct2='/data.lfpn/evidal/3Comp1DimWeijner/OutputData1D/data'
!      ct2='/scratch01/evidal/data'
      open(10,file =ct2 ,status = 'unknown',form = 'formatted')
      call out(t,Nx,Nc,gamma,ro,beta,cells)


!      ct2='OutputData2D/data'

 5    continue


      call ODE(t,Nx,Nc,beta,gamma,ro,cells)
      write(6,*) 'real t= '
      write(6,'(F6.2)') t/dk1
      counter=counter+1
      write(ct1,'(I4)') counter
      call out(t,Nx,Nc,gamma,ro,beta,cells)

!      if(mod(counter,85) .eq. 10)then
!         write(6,*) 'Perturbing'
!         do i=201,250
!               gamma(i)=gamma(i)+3.0
!         enddo
!      endif


      if (t+dt .lt. tend) then
         goto 5

      endif

      close(10)

!

!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     WRITES FINAL STATE
!      open(42,file ='/data.lfpn/evidal/3Comp1DimWeijner/OutputData1D/
!     ./Final-State',status = 'unknown',form = 'formatted')
!
!      open(43,file ='/data.lfpn/evidal/3Comp1DimWeijner/OutputData1D/
!     ./Final-Positions',status = 'unknown',form = 'formatted')
!
!      call outFinal(t,Nx,Nc,gamma,ro,beta,cells)
!      close(42)
!      close(43)

      end



