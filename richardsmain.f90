program richardsmain
  ! Description: Main program solving theta-based richards equation
  ! Reference: R.L. Zarba, A general Mass-Conservative Numerical Solution
  !            for the Unsaturated Flow Equation, Water Resources Research, 
  !            Vol.26, No.7, Pages:1483-1496, 1990
  ! USES:
  use soitexmod, only : prep, evap
  use saxtonscheme
  use tridiagonalmod
  use iomod
  implicit none
  !
  ! Variables:
  ! model setup options
  integer, parameter :: nstep = 96      ! number of total model steps
  real, parameter :: dt = 0.5          ! time per step (unit = h)
  integer, parameter :: nodes = 40     ! number of total layers
  real, parameter :: dz = 5.0          ! depth per layer (unit = cm)
  ! soil state variables
  real, dimension(nodes,1) :: theta       ! volumetric water content current step
  real, dimension(nodes,1) :: thetalast   ! volumetric water content previous step
  real, dimension(nodes,1) :: phi         ! soil matric water potential
  ! soil water characteristic parameters
  real, dimension(nodes,1) :: D           ! soil water diffusivity
  real, dimension(nodes,1) :: K           ! soil water conductivity
  real, dimension(nodes,1) :: vwcwp       ! volumetric soil water content at
                                          ! wilting point
  real, dimension(nodes,1) :: vwcfc       ! volumetric soil water content at
                                          ! field capacity
  real, dimension(nodes,1) :: vwcsat      ! saturated volumetric soil water 
                                          ! content
  ! other variables
  logical :: flag, notfirst               ! flag for picard iteration
  integer :: i, j, t                   ! indices
  integer :: iteration                 ! counter of iteration
  real, dimension(nodes,1) :: Kfactor      ! Vector of (K(lev+1)-K(lev))
  real, dimension(nodes,nodes) :: IM       ! Identity matrix
  real, dimension(nodes,nodes) :: Dplus    ! Diagonal matrix of D(lev+0.5)
  real, dimension(nodes,nodes) :: Dminus   ! Diagonal matrix of D(lev-0.5)
  real, dimension(nodes,nodes) :: deltaplus   ! Coefficients matrix of delta
  real, dimension(nodes,nodes) :: deltaminus  ! Coefficients matrix of delta
  real, dimension(nodes,nodes) :: A       ! tri-diagonal matrix for solving
                                          ! theta
  real, dimension(nodes,1) :: RHS     ! Residuals when using picard iteration
  real, dimension(nodes,1) :: delta   ! Difference between current theta and
                                      ! previous theta when iterating
!
!--------------------------------End of Definition-----------------------------
!
  ! Determine initial condition
  call moiregress(vwcwp(:,1), vwcfc(:,1), vwcsat(:,1), nodes)
  theta = vwcfc
  thetalast = vwcfc
!  write(6,*) 'theta: ', theta

  ! Determine matrix of delta and I
  do i = 1, nodes
    do j = 1, nodes
      if(i == j) then
        deltaplus(i, j) = -1.0
        deltaminus(i, j) = 1.0
        IM(i, j) = 1.0
      else if((i + 1) == j) then
        deltaplus(i, j) = 1.0
        IM(i ,j) = 0.0
      else if(i == (j + 1)) then
        deltaminus(i, j)= -1.0
        IM(i, j) = 0.0
      else
        deltaplus(i, j) = 0.0
        deltaminus(i, j) = 0.0
        IM(i, j) = 0.0
      end if
    end do
  end do
  deltaminus(1, 1) = 0.0
  deltaplus(nodes, nodes) = 0.0
  
  ! Start loop of water movement simulation
  do t = 1, nstep
    iteration = 0

    ! Get soil water characteristics through theta
    call watconduct(theta(:,1), vwcwp(:,1), vwcfc(:,1), vwcsat(:,1), &
                    K(:,1), nodes)
!    write(6,*) 'K: ', K
    call watdiffusiv(theta(:,1), vwcwp(:,1), vwcfc(:,1), vwcsat(:,1), &
                    K(:,1), D(:,1), nodes)
!    write(6,*) 'D: ', D
    
    ! Add Upper boundary condition
    ! Upper boundary condition is qflux = (prep - evap)/10
    ! derive new theta(1) by the eqation below
    ! theta(1) = dz*(qflux/D -1) + theta(1)
    theta(1,1) = dz*(0.1*(prep(t) - evap(t) - K(1,1)) / D(1,1)) + theta(1,1)
 
    ! Add Lower boundary condition
    ! Lower boundary condition is qflux = 0
    ! So D=0
    D(nodes, 1) = 0.0
      
    flag = .true.
    notfirst = .false.
    
    ! Picard iteration
    do while(flag)

      iteration = iteration + 1
      
      if(notfirst) then
        ! Get soil water characteristics through theta
        call watconduct(theta(:,1), vwcwp(:,1), vwcfc(:,1), vwcsat(:,1), &
                        K(:,1), nodes)
        call watdiffusiv(theta(:,1), vwcwp(:,1), vwcfc(:,1), vwcsat(:,1), &
                        K(:,1), D(:,1), nodes)
    
        D(nodes, 1) = 0.0
      end if

      ! Set up Dplus and Dminus
      ! Dplus and Dminus save the mean diffusivity between two layers
      ! by evaluate simple average
      do i = 1, nodes
        do j = 1, nodes
          if(i==j) then
            if(i == 1) then
              Dplus(i, j) = 0.5*D(i, 1) + 0.5*D(i + 1, 1)
              Dminus(i, j) = D(i, 1)
            else if(i == nodes) then
              Dplus(i, j) = D(i, 1)
              Dminus(i, j) = 0.5*D(i - 1, 1) + 0.5*D(i, 1)
            else
              Dplus(i, j) = 0.5*D(i, 1) + 0.5*D(i + 1, 1)
              Dminus(i, j) = 0.5*D(i - 1, 1) + 0.5*D(i, 1)
            end if
          else
            Dplus(i, j) = 0.0
            Dminus(i, j) = 0.0
          end if
        end do
      end do

      ! Set up Kfactor
      Kfactor(1, 1) = 0.5*K(2, 1) - 0.5*K(1, 1)
      Kfactor(nodes, 1) = 0.5*K(nodes, 1) - 0.5*K(nodes - 1, 1)
      j = nodes - 1
      do i = 2, j
        Kfactor(i, 1) = 0.5*K(i + 1, 1) - 0.5*K(i - 1, 1)
      end do

      ! Evaluate tri-diagonal matrix A
      A = 1.0/dt*IM - 1.0/(dz*dz)*(matmul(Dplus, deltaplus) - &
          matmul(Dminus, deltaminus))

      ! Evaluate Residual vector RHS
      RHS = 1.0/(dz*dz)*(matmul(Dplus, matmul(deltaplus, theta)) - &
            matmul(Dminus, matmul(deltaminus, theta))) + 1.0/dz*Kfactor - &
            1.0/dt*(theta - thetalast)
!      write(6,*) 'RHS: ', RHS
      ! Slove delta
      call tridiagonal(nodes, A, RHS, delta)
!      write(6,*) 'deltaplus: ', deltaplus
!      pause
!      write(6,*) 'deltaminus: ', deltaminus
!      pause
!      write(6,*) 'Kfactor: ', Kfactor
!      pause
!      write(6,*) 'A: ', A
!      pause
!      write(6,*) 'delta: ', delta
!      pause
!      write(6,*) 'maxabs(delta): ', maxval(abs(delta(:,1)))
!      pause
      ! If delta < 0.001, stop iteration and write out the result
      if(maxval(abs(delta(:, 1))) < 0.0002) then
        flag = .false.
        notfirst = .true.
        thetalast = theta
        theta = thetalast + delta
      else
        notfirst = .true.
        theta = theta + delta
      end if
    
    end do
    call soitension(theta(:,1), vwcwp(:,1), vwcfc(:,1), vwcsat(:,1), &
                    phi(:,1), nodes)
    call wrtout(theta, phi, nodes, t, iteration)
    
  end do

end program richardsmain
