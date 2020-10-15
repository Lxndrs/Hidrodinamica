program advection
  implicit none
  
  !Objetivo: Resolver la ecuacion u_t + au_x =0

  ! Declaracion de variables -----------------------------------------------------------

  real :: ncou=0.5, c, dt, h, a=0.2, tlim=1.5, ul=0.75, ur=0.2, xmin=0., xmax=1., t
  real :: U(0:101), Up(0:101)
  integer :: i, nstep, Nx=100, nmax=1000

  ! Significado de las variables:
  ! Nx: Resolucion de la malla
  ! ncou: numero de Courant
  ! c: variable adimensional
  ! dt: salto en el tiempo
  ! h: salto espacial
  ! a: velocidad
  ! tlim: tiempo limite de la simulacion
  ! ul: valor de la funcion al lado izquierdo del escalon
  ! ur: valor de la funcion al lado derecho del escalon
  ! xmin: extremo izquierdo de la malla espacial
  ! xmax: extremo derecho de la malla espacial
  ! t: tiempo actual
  ! i: contador espacial
  ! nstep: numero de pasos que le toma a la simulacion terminarse
  ! nmax: numero de pasos maximo para parar el programa en caso de fallo
  ! U: Funcion solucion al tiempo t
  ! Up: Funcion solucion al tiempo t+dt
  ! Call problem: condiciones iniciales y de frontera -----------------------------------

  open(1, file="CI-advection.dat", status="unknown")
  write(1, *)"step ", "x ", "U(x,0)"
  h = (xmax - xmin)/Nx ! definir malla espacial
  U(0) = ul ! Evitar que en el archivo de condiciones de frontera contengan basura
  U(Nx+1) = ur
  do i=1,Nx ! condiciones iniciales
     if(i.le.30) then
        U(i)=ul
     else
        U(i)=ur
     end if
  end do

  !Guardar en archivo condiciones iniciales
  do i=0, Nx+1
     write(1, *)i, i*h, U(i)
  end do
  close(1)
  
  ! Call Timestep -----------------------------------------------------------------------
  t = 0
  dt = ncou*h/a
  !Solve ---------------------------------------------------------------------------------
  open(2, file="advection-final.dat", status="unknown")
  open(3, file="advection-middle.dat", status="unknown")
  write(2, *)"step " ,"x ","U(x,tlim)"
  write(3, *)"step " ,"x ","U(x,t)"
 ! nstep=1 ! inicializar contador de pasos
  do while(t.le.tlim)!.or.(nstep.ge.nmax))
     t = t + dt ! Salto en el tiempo
     c = a*dt/h
     do i=1, Nx        
        Up(i) = U(i) - c*(U(i)-U(i-1))
        if  ((t.lt.tlim/2. + 0.5*dt).and.(t.gt.tlim/2. - 0.5*dt)) then 
           write(3, *)i, i*h, Up(i)
        end if
     end do
     !condiciones de frontera
     Up(0) = Up(1)
     Up(Nx+1) = Up(Nx)
     U=Up
!     nstep = nstep+1
  end do
  do i=0, Nx+1
     write(2, *)i, i*h, U(i)
  end do
  print*, "Done, hope everything is all right."
end program advection

