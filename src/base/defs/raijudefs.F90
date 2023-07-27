! kaimag definitions/constants

module raijudefs
    use kdefs

    implicit none

    !------
    ! Enumerators
    !------

    ! Grid settings
    enum, bind(C)
        enumerator :: G_UNISPH, G_SHGRID
    endenum

    ! Species
    enum, bind(c)
        enumerator :: F_PSPH=0,F_HOTE,F_HOTP
            !! These flavors have reserved numbers
    endenum

    ! Topology
    enum, bind(C)
        enumerator :: RAIJUOPEN=0, RAIJUCLOSED
            !! Whether the field line corresponding to grid point is open or closed
    endenum

    ! Active/buffer/inactive cells
    enum, bind(C)
        enumerator :: RAIJUINACTIVE=-1, RAIJUBUFFER, RAIJUACTIVE
            !! Helps determine how RAIJU is going to treat the grid point
    endenum

    ! Species type
    enum, bind(C)
        enumerator :: RAIJUNSPC=0,RAIJUELE,RAIJUHPLUS,RAIJUOPLUS
                    ! Null species, electron, h+, o+
    endenum

    !------
    ! Defaults
    !------

    ! Units
    real(rp) :: sclEta = 1.0e9  ! [1/nT -> 1/T on DkT2eta conversion]
    real(rp) :: sclIntens = 1.e-4*sqrt(ev2J/(8.0*dalton))/PI ! code eta to intensity [1/(s*sr*keV*cm^2)]

    ! Settings
    integer :: nSpacesDef = 4
        !! Number of i spaces between last good value and active i for species
    real(rp) :: fracWorthyDef = 0.001
        !! Fraction that a lambda channel must contribute to total pressure or density in order to be worthy of being evolved

end module raijudefs