! MIX definitions/constants

module mixdefs
  use kdefs
  implicit none

  integer, parameter :: nVars = 16 ! change together wiht the enumerator below
  enum, bind(C)
     enumerator :: POT=1,FAC,SIGMAP,SIGMAH,SOUND_SPEED,DENSITY,AVG_ENG,NUM_FLUX,NEUTRAL_WIND,EFIELD,IM_EAVG,IM_EFLUX,IM_IAVG,IM_IFLUX,Z_EAVG,Z_NFLUX
  end enum

  ! enumerator for MHD->MIX variables
  enum, bind(C)
     enumerator :: MHDJ=1, MHDD, MHDC
  end enum

  ! enumerator for MIX->MHD variable(s)
  enum, bind(C)
     enumerator :: MHDPSI=1
  end enum

  enum, bind(C)
     enumerator :: NORTH=1,SOUTH
  end enum

  enum, bind(C)
     enumerator :: AMIE=1,MOEN_BREKKE
  end enum

  enum, bind(C)
     enumerator :: FEDDER=1,ZHANG,RCMONO,RCMFED
  end enum

  enum, bind(C)
     enumerator :: NONE=1,ADHOC,CMIT
  end enum

end module mixdefs
