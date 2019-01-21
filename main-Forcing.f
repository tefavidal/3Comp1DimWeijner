      program main

      implicit none


      INTEGER, PARAMETER :: Nx=2000
      INTEGER, PARAMETER :: Nc=800
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
      integer pgbeg, i, counter, factor, ForcingPeriod
      integer grid(Nx)
      character(len=57) ct1
      character(len=48) ct2
!      character(len=22) ct2
      character(len=9) ct3
      character(len=100) ct
      double precision gamma0(10),ro0(10),beta0(10)

      t=0.d0
      counter=0;
      call anfang(t,Nx,Nc,beta,gamma,ro,cells)
!%%%%%%%%% Reading imput parameter
      if (command_argument_count() .eq. 0)then
        write(6,*) 'No forcing period given'
            write(6,*) 'Emergency Exit'
            call EXIT(0)
      else

        CALL get_command_argument(1, ct)
        read(ct,*) ForcingPeriod
        write(6,*) 'Using Forcing Period, T=',ForcingPeriod
      endif
!%%%%%%%%%% end reading imput

!%%%%% Output file creation and first line
      write(ct3,'(A,F3.1,A,I2)')
     . 'ke', dke0,'_F',ForcingPeriod
      ct2='/data.lfpn/evidal/3Comp1DimWeijner/OutputData1D/'
       ct1=ct2 // ct3
      open(10,file =ct1 ,status = 'unknown',form = 'formatted')
      call out(t,Nx,Nc,gamma,ro,beta,cells)



 5    continue


      call ODE(t,Nx,Nc,beta,gamma,ro,cells)
      write(6,*) 'real t= '
      write(6,'(F6.2)') t/dk1
      counter=counter+1
      call out(t,Nx,Nc,gamma,ro,beta,cells)

      if(mod(counter,ForcingPeriod) .eq. 10)then
         write(6,*) 'Perturbing'
         do i=201,250
               gamma(i)=gamma(i)+3.0
         enddo
      endif


      if (t+dt .lt. tend) then
         goto 5

      endif

      close(10)

!

      end



