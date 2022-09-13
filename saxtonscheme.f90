module saxtonscheme
! Description: Determine parameters by using the scheme from Saxton 
! (Saxton, 2006)
  implicit none
  public :: moiregress
  public :: soitension
  public :: watconduct
  public :: watdiffusiv
!
contains
!
  subroutine moiregress(vwcwp, vwcfc, vwcsat, siz)
  ! Description: Determine vwc thresholds being used in water conductivity 
  ! solution
  ! USES:
    use soitexmod, only : S, C, OM
    implicit none
  !
    ! Arguments
    ! Volumetric soil water content at wilting point (1500kPa) for each layer
    real, dimension(:), intent(out) :: vwcwp
    ! Volumetric soil water content at field capacity (33kPa) for each layer
    real, dimension(:), intent(out) :: vwcfc
    ! Saturated volumetric soil water content
    real, dimension(:), intent(out) :: vwcsat
    ! total number of layers
    integer, intent(in) :: siz
  !
    ! Local variables
    integer :: i      ! indice
    real, dimension(1:siz) :: temp1, temp2    ! temporary variables
  !
    ! Calculate vwcwp
    do i = 1,siz
      temp1(i) = -0.024*S(i) + 0.487*C(i) + 0.006*OM + 0.005*S(i)*OM - &
                 0.013*C(i)*OM + 0.068*S(i)*C(i) + 0.031
      vwcwp(i) = 1.14*temp1(i) - 0.02
    end do
    ! Calculate vwcfc
    do i = 1,siz
      temp1(i) = -0.251*S(i) + 0.195*C(i) + 0.011*OM + 0.006*S(i)*OM - &
                 0.027*C(i)*OM + 0.452*S(i)*C(i) + 0.299
      vwcfc(i) = 1.283*temp1(i)*temp1(i) + 0.626*temp1(i) - 0.015
    end do
    ! Calculate vwcsat
    do i = 1,siz
      temp1(i) = 0.278*S(i) + 0.034*C(i) + 0.022*OM - 0.018*S(i)*OM - &
                 0.027*C(i)*OM - 0.584*S(i)*C(i) + 0.078
      temp2(i) = 1.636*temp1(i) - 0.107
      vwcsat(i) = vwcfc(i) + temp2(i) - 0.097*S(i) + 0.043
    end do
    
  end subroutine moiregress
!
!------------------------------------------------------------------------------
!
  subroutine soitension(vwc, vwcwp, vwcfc, vwcsat, soiten, siz)
  ! Description: solve soil water tension for each layer
  ! USES:
    use soitexmod, only : S, C, OM
    implicit none
    save
  !
    ! Arguments
    ! Volumetric soil water content at current step
    real, dimension(:), intent(in) :: vwc
    ! Volumetric soil water content at wilting point (1500kPa) for each layer
    real, dimension(:), intent(in) :: vwcwp
    ! Volumetric soil water content at field capacity (33kPa) for each layer
    real, dimension(:), intent(in) :: vwcfc
    ! Saturated volumetric soil water content
    real, dimension(:), intent(in) :: vwcsat
    ! Soil tension (in cm)
    real, dimension(:), intent(out) :: soiten
    ! total number of layers
    integer, intent(in) :: siz
  !
    ! Local variables
    integer :: i      ! indice
    real, dimension(1:siz) :: temp1, temp2    ! temporary variables
  !
    ! Calculate soil moisture tension
    do i = 1,siz
      if(vwc(i) < vwcfc(i)) then       ! soiten > 33kPa
        temp2 = (log(1500.0) - log(33.0)) / (log(vwcfc(i)) - log(vwcwp(i)))
        temp1 = exp(log(33.0)+ log(vwcfc(i))*temp2(i))
        soiten(i) = temp1(i)*(vwc(i)**(-temp2(i)))
      else if(vwc(i) < vwcsat(i)) then    ! soiten < 33kPa and vwc < vwcsat
        temp1(i) = 0.278*S(i) + 0.034*C(i) + 0.022*OM - 0.018*S(i)*OM - &
                   0.027*C(i)*OM - 0.584*S(i)*C(i) + 0.078
        temp1(i) = 1.636*temp1(i) - 0.107
        temp2(i) = -21.67*S(i) - 27.93*C(i) - 81.97*temp1(i) + &
                   71.12*S(i)*temp1(i) + 8.29*C(i)*temp1(i) + &
                   14.05*S(i)*C(i) + 27.16
        soiten(i) = 0.02*temp2(i)*temp2(i) + 0.887*temp2(i) - 0.07
        soiten(i) = 33.0 - (vwc(i) - vwcfc(i))*(33.0 - soiten(i)) / &
                    (vwcsat(i) - vwcfc(i))
      else       ! vwc = vwcsat
        temp1(i) = 0.278*S(i) + 0.034*C(i) + 0.022*OM - 0.018*S(i)*OM - &
                   0.027*C(i)*OM - 0.584*S(i)*C(i) + 0.078
        temp1(i) = 1.636*temp1(i) - 0.107
        temp2(i) = -21.67*S(i) - 27.93*C(i) - 81.97*temp1(i) + &
                   71.12*S(i)*temp1(i) + 8.29*C(i)*temp1(i) + &
                   14.05*S(i)*C(i) + 27.16
        soiten(i) = 0.02*temp2(i)*temp2(i) + 0.887*temp2(i) - 0.07
      end if
      ! Change unit from 'kPa' to 'cm'
      soiten(i) = 0.01*1030.3*soiten(i)
    end do

  end subroutine soitension
!
!------------------------------------------------------------------------------
!
  subroutine watconduct(vwc, vwcwp, vwcfc, vwcsat, conduct, siz)
  ! Description: solve soil water conductivity for each layer
  ! USES:
    implicit none
    save
  !
    ! Arguments
    ! Volumetric soil water content at current step
    real, dimension(:), intent(in) :: vwc
    ! Volumetric soil water content at wilting point (1500kPa) for each layer
    real, dimension(:), intent(in) :: vwcwp
    ! Volumetric soil water content at field capacity (33kPa) for each layer
    real, dimension(:), intent(in) :: vwcfc
    ! Saturated volumetric soil water content
    real, dimension(:), intent(in) :: vwcsat
    ! Water conductivity (in cm/h)
    real, dimension(:), intent(out) :: conduct
    ! total number of layers
    integer, intent(in) :: siz
  !
    ! Local variables
    integer :: i      ! indice
    real, dimension(1:siz) :: temp1, temp2    ! temporary variables
  !
    ! Calculate water conductivity
    do i = 1,siz
      if(vwc(i) < vwcsat(i)) then    ! Unsaturated water
        temp2 = (log(1500.0) - log(33.0)) / (log(vwcfc(i)) - log(vwcwp(i)))
        temp1(i) = 1930.0*(vwcsat(i) - vwcfc(i))**(3.0 - 1.0 / temp2(i))
        conduct(i) = temp1(i)*(vwc(i) / vwcsat(i))**(3.0 + 2.0*temp2(i))
      else       ! Saturated water
        temp2 = (log(1500.0) - log(33.0)) / (log(vwcfc(i)) - log(vwcwp(i)))
        conduct(i) = 1930.0*(vwcsat(i) - vwcfc(i))**(3.0 - 1.0 / temp2(i))
      end if
      ! Change unit from 'mm/h' to 'cm/h'
      conduct(i) = 0.1*conduct(i)
    end do

  end subroutine watconduct
!
!------------------------------------------------------------------------------
!
  subroutine watdiffusiv(vwc, vwcwp, vwcfc, vwcsat, conduct, diffusiv, siz)
  ! Description: solve soil water diffusivity for each layer
  ! D(theta) = K(theta) * dh / d(theta)
  ! D ------- soil water diffusivity
  ! K ------- soil water conductivity
  ! dh / d(theta) ------- soil water-retention relation, h is soil matric 
  !                       potential (cm)
  ! USES:
    use soitexmod, only : S, C, OM
    implicit none
    save
  !
    ! Arguments
    ! Volumetric soil water content at current step
    real, dimension(:), intent(in) :: vwc
    ! Volumetric soil water content at wilting point (1500kPa) for each layer
    real, dimension(:), intent(in) :: vwcwp
    ! Volumetric soil water content at field capacity (33kPa) for each layer
    real, dimension(:), intent(in) :: vwcfc
    ! Saturated volumetric soil water content
    real, dimension(:), intent(in) :: vwcsat
    ! Soil water conductivity (in cm/h)
    real, dimension(:), intent(in) :: conduct
    ! Soil water difusivity (in cm^2/h)
    real, dimension(:), intent(out) :: diffusiv
    ! total number of layers
    integer, intent(in) :: siz
  !
    ! Local variables
    integer :: i      ! indice
    real, dimension(1:siz) :: temp1, temp2    ! temporary variables
    real, dimension(1:siz) :: phie            ! soitension at air entry
  !
    do i = 1,siz
      ! Calculate soitension at air entry
      temp1(i) = 0.278*S(i) + 0.034*C(i) + 0.022*OM - 0.018*S(i)*OM - &
                   0.027*C(i)*OM - 0.584*S(i)*C(i) + 0.078
      temp1(i) = 1.636*temp1(i) - 0.107
      temp2(i) = -21.67*S(i) - 27.93*C(i) - 81.97*temp1(i) + &
                   71.12*S(i)*temp1(i) + 8.29*C(i)*temp1(i) + &
                   14.05*S(i)*C(i) + 27.16
      phie(i) = 0.02*temp2(i)*temp2(i) + 0.887*temp2(i) - 0.07
      ! Calculate parameter A and B
      temp2 = (log(1500.0) - log(33.0)) / (log(vwcfc(i)) - log(vwcwp(i)))
      temp1 = exp(log(33.0)+ log(vwcfc(i))*temp2(i))
      ! Calculate diffusivity 
      if(vwc(i) < vwcfc(i)) then
        diffusiv(i) = 0.01*1030.3*temp1(i)*temp2(i)*(vwc(i)**(-temp2(i)-1))*conduct(i)
      else       ! Saturated water
        diffusiv(i) = -0.01*1030.3*(phie(i)-33.0)/(vwcsat(i)-vwcfc(i))*conduct(i)
      end if
    end do

  end subroutine watdiffusiv
  
end module saxtonscheme
