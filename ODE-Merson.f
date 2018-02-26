      subroutine ODE(t,Nx,Nc,beta,gamma,ro,cells)
      
      implicit none
      double precision t
      integer Nx,Nc,i

      double precision beta(Nx),gamma(Nx),ro(Nx)
      double precision beta0(Nx),gamma0(Nx),ro0(Nx)
      double precision betak1(Nx),gammak1(Nx),rok1(Nx)
      double precision betak2(Nx),gammak2(Nx),rok2(Nx)
      double precision betak3(Nx),gammak3(Nx),rok3(Nx)
      double precision betak4(Nx),gammak4(Nx),rok4(Nx)
      double precision betak5(Nx),gammak5(Nx),rok5(Nx)
      double precision b1(Nx),g1(Nx),r1(Nx)
      double precision vdx(Nx), vdy(Nx)
      double precision cells(Nc)
      integer grid(Nx)

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

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

       do i=1,Nx
        beta0(i)=beta(i)
        gamma0(i)=gamma(i)
        ro0(i)=ro(i)
       enddo

      iteration=0

      call rs(t,Nx,beta0,gamma0,ro0,betak1,gammak1,rok1,grid,vdx)
!     Runge-Kutta-Merson Method

 16      do i=1,Nx
        beta(i)=beta0(i)+h*betak1(i)/3
        gamma(i)=gamma0(i)+h*gammak1(i)/3
        ro(i)=ro0(i)+h*rok1(i)/3
       enddo


      call rs(t+h/3,Nx,beta,gamma,ro,betak2,gammak2,rok2,grid,vdx)


       do i=1,Nx
        beta(i)=beta0(i)+h*(betak1(i)+betak2(i))/6
        gamma(i)=gamma0(i)+h*(gammak1(i)+gammak2(i))/6
        ro(i)=ro0(i)+h*(rok1(i)+rok2(i))/6
       enddo


      call rs(t+h/3,Nx,beta,gamma,ro,betak3,gammak3,rok3,grid,vdx)


       do i=1,Nx
        beta(i)=beta0(i)+h*(betak1(i)+3*betak3(i))/8
        gamma(i)=gamma0(i)+h*(gammak1(i)+3*gammak3(i))/8
        ro(i)=ro0(i)+h*(rok1(i)+3*rok3(i))/8
       enddo


      call rs(t+h/2,Nx,beta,gamma,ro,betak4,gammak4,rok4,grid,vdx)


       do i=1,Nx
        beta(i)=beta0(i)+h*(betak1(i)-3*betak3(i)
     .   +4*betak4(i))/2
        gamma(i)=gamma0(i)+h*(gammak1(i)-3*gammak3(i)
     .   +4*gammak4(i))/2
        ro(i)=ro0(i)+h*(rok1(i)-3*rok3(i)
     .   +4*rok4(i))/2
       enddo


      call rs(t+h,Nx,beta,gamma,ro,betak5,gammak5,rok5,grid,vdx)


       do i=1,Nx
        beta(i)=beta0(i)+h*(betak1(i)+4*betak4(i)
     .   +betak5(i))/6
        gamma(i)=gamma0(i)+h*(gammak1(i)+4*gammak4(i)
     .   +gammak5(i))/6
        ro(i)=ro0(i)+h*(rok1(i)+4*rok4(i)
     .   +rok5(i))/6
       enddo



       do i=1,Nx
        b1(i)=beta(i)-h*(betak1(i)+betak2(i)+betak3(i)
     .   +betak4(i)+betak5(i))/5-beta0(i)
        g1(i)=gamma(i)-h*(gammak1(i)+gammak2(i)+gammak3(i)
     .   +gammak4(i)+gammak5(i))/5-gamma0(i)
        r1(i)=ro(i)-h*(rok1(i)+rok2(i)+rok2(i)+rok3(i)
     .   +rok4(i)+rok5(i))/5-ro0(i)
       enddo


      err=0.d0
      index=0

       do i=1,Nx
       err=max(abs(b1(i)),abs(g1(i)),abs(r1(i)),err)
      if (beta(i) .lt. 0 .or. gamma(i) .lt. 0 .or. ro(i) .lt. 0)
     . then
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

      call UpdateDiscreteCells(h,Nx,Nc, gamma, ro, cells)

      h=dt

      if (tau + h .le. tout+tol*dt) then


       go to 13
      elseif(tau .lt. tout-tol*dt)then
         h = tout - tau

         go to 13
      endif


      return
      end




