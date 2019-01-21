      program main

      implicit none


      INTEGER, PARAMETER :: Nx=2000
      INTEGER, PARAMETER :: Nc=1000
!     Nc decided as Nx*percentageofcells

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision beta(Nc),gamma(Nx),ro(Nc)
      double precision t
      double precision cells(Nc)
      integer pgbeg, i, counter, factor
      integer grid(Nx)
      character(len=64) ct1
      character(len=48) ct2
!      character(len=22) ct2
      character(len=16) ct3

      double precision gamma0(10),ro0(10),beta0(10)

      t=0.d0
      counter=0;
      call anfang(t,Nx,Nc,beta,gamma,ro,cells)
      write(ct3,'(A,F3.1,A,F3.0,A,F4.2)')
     . 'ke', keU,'_',isf,'_ep',
     .  0.2575/dlambda2
      ct2='/data.lfpn/evidal/3Comp1DimWeijner/OutputData1D/'
!      ct2='/scratch01/evidal/data'
       ct1=ct2 // ct3
!       i=index(ct1,'.')
!       ct1(i:i)='_'
      open(10,file =ct1 ,status = 'unknown',form = 'formatted')
      call out(t,Nx,Nc,gamma,ro,beta,cells)



 5    continue


      call ODE(t,Nx,Nc,beta,gamma,ro,cells)
      write(6,*) 'real t= '
      write(6,'(F6.2)') t/dk1
      counter=counter+1
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

      end



