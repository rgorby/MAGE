module testHelper

    ! limiting namespace pollution from pFUnit
    use pFUnit, only: SourceLocation,anyExceptions,&
                      assertTrue,assertFalse,assertEqual,assertLessThanOrEqual, &
                      assertLessThan,assertGreaterThan,assertGreaterThanOrEqual, &
                      assertAny,assertAll,assertNotAll,assertIsFinite, &
                      assertNone,assertExceptionRaised,assertSameShape,assertIsNaN

    implicit none

    contains

end module testHelper

