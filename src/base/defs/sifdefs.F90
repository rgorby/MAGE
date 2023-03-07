! kaimag definitions/constants

module sifdefs
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
        enumerator :: SIFOPEN, SIFCLOSED
    endenum

    ! Active/buffer/inactive cells
    enum, bind(C)
        enumerator :: SIFACTIVE, SIFBUFFER, SIFINACTIVE
    endenum

    !------
    ! Defaults
    !------


end module sifdefs