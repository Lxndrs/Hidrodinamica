program burguers
  implicit none
  
  !Objetivo: Resolver las ecuaciones de Euler en una dimension

  ! Declaracion de variables -----------------------------------------------------------

  real, parameter :: ncou=0.5, tlim=1.5, xmin=0., xmax=1., gamma=5./3, Nx=100, rl=0.8, ru=0.2
  real :: c, dt, h, t
  real :: U1(0:101), Up1(0:101), U2(0:101), Up2(0:101), U3(0:101), Up3(0:101)
  real :: rho(0:101), p(0:101), eint(0:101), E(0:101), u(0:101), a(0:101)
  real :: A21(0:101), A22(0:101), A23(0:101), A31(0:101), A32(0:101), A33(0:101)
  integer :: i, 

  ! Significado de las variables:
  ! Nx: Resolucion de la malla
  ! ncou: numero de Courant
  ! c: variable adimensional
  ! dt: salto en el tiempo
  ! h: salto espacial
  ! tlim: tiempo limite de la simulacion
  ! ul: valor de la funcion al lado izquierdo del escalon
  ! ur: valor de la funcion al lado derecho del escalon
  ! xmin: extremo izquierdo de la malla espacial
  ! xmax: extremo derecho de la malla espacial
  ! t: tiempo actual
  ! i: contador espacial
  ! nstep: numero de pasos que le toma a la simulacion terminarse
  ! nmax: numero de pasos maximo para parar el programa en caso de fallo
  ! U123: Funcion solucion al tiempo t de la variable conservada 1, 2 o 3
  ! Up123: Funcion solucion al tiempo t+dt de la variable conservada 1, 2 o 3
  ! Call problem: condiciones iniciales y de frontera -----------------------------------

  open(1, file="CI-euler.dat", status="unknown")
  write(1, *)"step ", "x ", "rho ", "u ", "p ", "e ", "E "
  h = (xmax - xmin)/Nx ! definir malla espacial
 
  do i=1,Nx ! condiciones iniciales: u0=0, p y rho son iguales a la funcion escalon
     u(i)= 0
     if (i.lt.30) then
        rho(i)=rl
     else
        rho(i)=ru
     end if
     p(i)=rho(i)
     eint(i)=p(i)/(gamma-1)/rho(i)
     E(i)=rho(i)*(0.5*u(i)**2+eint(i))
  end do
  u(0)=u(1) ! Primera actualizacion de condiciones de frontera
  u(Nx+1)=u(Nx)
  rho(0) = rho(1)
  rho(Nx+1) = rho(Nx)
  p(1)=p(0)
  p(Nx+1)=p(Nx)
  eint(0)= eint(1)
  eint(Nx+1) = eint(Nx)
  E(0) = E(1)
  E(Nx+1)=E(Nx)

  !Guardar en archivo condiciones iniciales
  do i=0, Nx+1
     write(1, *)i, i*h, rho(i), u(i), p(i), eint(i), E(i)
     U1(i)=rho(i)   ! Cambio a cantidades conservadas
     U2(i)=rho(i)*u(i)
     U3(i)=E(i)
     a(i)=sqrt(gamma*p(i)/rho(i))
  end do
  close(1)


  
  
  
  ! Call Timestep -----------------------------------------------------------------------
  t = 0
  dt = ncou*h ! Medida temporal
  !Solve ---------------------------------------------------------------------------------
  open(2, file="euler-final.dat", status="unknown")
  open(3, file="euler-middle.dat", status="unknown")
  write(2, *)"step ", "x ", "rho ", "u ", "p ", "e ", "E "
  write(3, *)"step ", "x ", "rho ", "u ", "p ", "e ", "E "

  do while(t.le.tlim)
     t = t + dt ! Salto en el tiempo
     c = dt/h
     
     do i=1, Nx
        A21(i) = 0.5*(gamma-3)*u(i)**2
        A22(i) = (3-gamma)*u(i)
        A23(i) = (gamma-1)
        A31(i) = 0.5*(gamma-2)*u(i)**3 - (a(i)**2*u(i))/(gamma-1)
        A32(i) = u(i)**2/(gamma-1) + 0.5*(3-2*gamma)*u(i)**2
        A33(i) = gamma-1
        Up1(i) = U1(i) -c*(U2(i)-U2(i-1))
        Up2(i) = U2(i) -c*(A21(i)*(U1(i)-U1(i-1))+A22(i)*(U2(i)-U2(i-1))+ A23(i)*(U3(i)-U3(i-1)))
        Up3(i) = U3(i) -c*(A31(i)*(U1(i)-U1(i-1))+ A32(i)*(U2(i)-U2(i-1))+ A33(i)*(U3(i)-U3(i-1)))
        !Actualizacion de variables primitivas
        rho(i)=Up1
        u(i)=Up2(i)/rho(i)
        E(i) = Up3(i)
        eint(i) = E(i)/rho(i) - 0.5*u(i)**2
        p(i) = (gamma-1)*rho(i)*eint(i)
        a(i) = sqrt(gamma*p(i)/rho(i))
        if  ((t.lt.tlim/2. + 0.5*dt).and.(t.gt.tlim/2. - 0.5*dt)) then 
           write(3, *)i, i*h, rho(i), u(i), p(i), eint(i), E(i) 
        end if
     end do
     ! Condiciones de frontera
     Up1(0) = Up1(1)
     Up1(Nx+1) = Up1(Nx)
     Up2(0) = Up1(1)
     Up2(Nx+1) = Up1(Nx)
     Up3(0) = Up1(1)
     Up3(Nx+1) = Up1(Nx)

     !Actualizacion de variables conservadas
     U1=Up1
     U2=Up2
     U3=Up3

     
  end do
  do i=0, Nx+1
     write(2, *)i, i*h, rho(i), u(i), p(i), eint(i), E(i)
  end do
  print*, "Done, hope everything is all right."
end program burguers

!subroutine frontier_update()
