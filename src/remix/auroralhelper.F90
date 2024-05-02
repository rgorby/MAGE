!Module to hold routines to help calculate auroral conductance

module auroralhelper
    use kdefs
    use mixdefs
    use mixtypes
    
    implicit none

    contains

!Functions for auroral conductance from precipitation values
    !Returns Robinson's SigP from eavg [kEv] and eflux [ergs/cm^2]
    elemental function SigmaP_Robinson(eavg,eflux) result(SigP)
        real(rp), intent(in) :: eavg,eflux
        real(rp) :: SigP
        SigP = 40.0*eavg*sqrt(eflux)/(16.0 + eavg**2.0)
    end function SigmaP_Robinson

    !Returns Robinson's SigH from eavg [kEv] and eflux [ergs/cm^2]
    !NOTE: Extra correction from Fedder
    elemental function SigmaH_Robinson(eavg,eflux) result(SigH)
        real(rp), intent(in) :: eavg,eflux
        real(rp) :: SigH
        real(rp) :: SigP

        SigP = SigmaP_Robinson(eavg,eflux)
        !Use to use Fedder correction: SigH = 0.45*SigP*(eavg**0.85)/(1.0 + 0.0025*eavg**2.0)
        !Use Kaeppler+ 15 correction to SigH/SigP from Robinson
        SigH = 0.57*SigP*(eavg**0.53)
    end function SigmaH_Robinson

end module auroralhelper
