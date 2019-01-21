      subroutine ODE(t,Nx,Nc,beta,gamma,ro,cells)
      
      implicit none
      double precision t
      integer Nx,Nc,i

      double precision beta(Nc),gamma(Nx),ro(Nc)
      double precision beta0(Nc),gamma0(Nx),ro0(Nc)
      double precision betak1(Nc),gammak1(Nx),rok1(Nc)
      double precision betak2(Nc),gammak2(Nx),rok2(Nc)
      double precision betak3(Nc),gammak3(Nx),rok3(Nc)
      double precision betak4(Nc),gammak4(Nx),rok4(Nc)
      double precision betak5(Nc),gammak5(Nx),rok5(Nc)
      double precision b1(Nc),g1(Nx),r1(Nc)
      double precision vdx(Nx), vdy(Nx)
      double precision cells(Nc)
      integer grid(Nx)

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob,keB,keU

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision tau,h,err
      integer iteration, index

      tau=0.d0
      h=dt
!      call flow(t,Nx,vdx,vdy)
             do i=1,Nx
                vdx(i)=vd
                vdy(i)=0
                enddo

 13   call FromCellToGrid(Nx,Nc,cells,grid)

       do i=1,Nc
        beta0(i)=beta(i)
        ro0(i)=ro(i)
       enddo
       do i=1,Nx
        gamma0(i)=gamma(i)
       enddo
      iteration=0

      call rs(t,Nx,Nc,beta0,gamma0,ro0,betak1,gammak1,rok1,cells,vdx)
!     Runge-Kutta-Merson Method

 16      do i=1,Nx
        gamma(i)=gamma0(i)+h*gammak1(i)/3
       enddo
       do i=1,Nc
        beta(i)=beta0(i)+h*betak1(i)/3
        ro(i)=ro0(i)+h*rok1(i)/3
       enddo
      call rs(t+h/3,Nx,Nc,beta,gamma,ro,betak2,gammak2,rok2,cells,vdx)


       do i=1,Nc
        beta(i)=beta0(i)+h*(betak1(i)+betak2(i))/6
        ro(i)=ro0(i)+h*(rok1(i)+rok2(i))/6
       enddo
       do i=1,Nx
        gamma(i)=gamma0(i)+h*(gammak1(i)+gammak2(i))/6
       enddo

      call rs(t+h/3,Nx,Nc,beta,gamma,ro,betak3,gammak3,rok3,cells,vdx)


       do i=1,Nc
        beta(i)=beta0(i)+h*(betak1(i)+3*betak3(i))/8
        ro(i)=ro0(i)+h*(rok1(i)+3*rok3(i))/8
       enddo
       do i=1,Nx
        gamma(i)=gamma0(i)+h*(gammak1(i)+3*gammak3(i))/8
       enddo

      call rs(t+h/2,Nx,Nc,beta,gamma,ro,betak4,gammak4,rok4,cells,vdx)


       do i=1,Nc
        beta(i)=beta0(i)+h*(betak1(i)-3*betak3(i)
     .   +4*betak4(i))/2
        ro(i)=ro0(i)+h*(rok1(i)-3*rok3(i)
     .   +4*rok4(i))/2
       enddo
       do i=1,Nx
        gamma(i)=gamma0(i)+h*(gammak1(i)-3*gammak3(i)
     .   +4*gammak4(i))/2
       enddo

      call rs(t+h,Nx,Nc,beta,gamma,ro,betak5,gammak5,rok5,cells,vdx)


       do i=1,Nc
        beta(i)=beta0(i)+h*(betak1(i)+4*betak4(i)
     .   +betak5(i))/6
        ro(i)=ro0(i)+h*(rok1(i)+4*rok4(i)
     .   +rok5(i))/6
       enddo
       do i=1,Nx
        gamma(i)=gamma0(i)+h*(gammak1(i)+4*gammak4(i)
     .   +gammak5(i))/6
       enddo


       do i=1,Nc
        b1(i)=beta(i)-h*(betak1(i)+betak2(i)+betak3(i)
     .   +betak4(i)+betak5(i))/5-beta0(i)
        r1(i)=ro(i)-h*(rok1(i)+rok2(i)+rok2(i)+rok3(i)
     .   +rok4(i)+rok5(i))/5-ro0(i)
       enddo
       do i=1,Nx
        g1(i)=gamma(i)-h*(gammak1(i)+gammak2(i)+gammak3(i)
     .   +gammak4(i)+gammak5(i))/5-gamma0(i)
       enddo

      err=0.d0
      index=0

       do i=1,Nx
        err=max(abs(g1(i)),err)
        if (gamma(i) .lt. 0)then
         index=1
        endif
       enddo

       do i=1,Nc
        err=max(abs(b1(i)),abs(r1(i)),err)
        if (beta(i) .lt. 0 .or. ro(i) .lt. 0)then
         index=1
        endif
       enddo

      if (err .gt. tol .or. index .eq. 1) then
      h=h/2
      iteration=iteration+1

        if (iteration .gt. 2) then
            write(6,*) 't =',t,' index =',index, 'iteration=',iteration
        endif
        if (iteration .gt. 40) then
            write(6,*) 'Emergency Exit'
            call EXIT(0)
        endif

      go to 16
      endif


      t=t+h
      tau=tau+h

      call UpdateDiscreteCells(t,h,Nx,Nc, gamma, ro, beta, cells)
!      call ScalatedMovement(t,h,Nx,Nc, gamma, ro, beta, cells)
      h=dt

      if (tau + h .le. tout+tol*dt) then


       go to 13
      elseif(tau .lt. tout-tol*dt)then
         h = tout - tau

         go to 13
      endif


      return
      end




