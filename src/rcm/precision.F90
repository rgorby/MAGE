module rcm_precision
    ! preventing rcm namespace pollution
    use kdefs, ONLY: kip => ip, krp => rp
    private kip,krp
  ! use gamera precision
    INTEGER, PARAMETER :: iprec = kip
    INTEGER, PARAMETER :: rprec = krp
end module rcm_precision
