module raijuELossWM
    use ioh5
    use files
    
    use raijudefs
    use raijutypes

    implicit none

    integer, parameter, private :: MAXIOVAR = 50

    contains

    subroutine initEWM(eWM, configFname)
        type(eLossWM_T), intent(inout) :: eWM
        character(len=strLen) :: configFname

        type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
        type(IOVar_T) :: tauVAR

        call ClearIO(IOVars)

        if(.not. ioExist(configfname, "Tau", "waveModel")) then
            write(*,*) "This config file not structured for RAIJU, use genRAIJU.py. Good day."
            stop
        endif

        !Chorus wave
        call AddInVar(IOVars,"Kp")  
        call AddInVar(IOVars,"MLT") 
        call AddInVar(IOVars,"L") 
        call AddInVar(IOVars,"Ek") 
        call AddInVar(IOVars,"Tau")
        call ReadVars(IOVars,.false.,configFname, "waveModel")

        eWM%Nkp  = IOVars(FindIO(IOVars, "Kp"))%N
        eWM%Nmlt = IOVars(FindIO(IOVars, "MLT"))%N
        eWM%Nl   = IOVars(FindIO(IOVars, "L"))%N
        eWM%Ne   = IOVars(FindIO(IOVars, "Ek"))%N

        tauVAR = IOVars(FindIO(IOVars, "Tau"))

        ! Do some dimension checks
        if (tauVAR%Nr /= 4) then
            write(*,*) "tauDim:",tauVAR%Nr
            write(*,*) 'Currently only support tau model files in the form tau(Kp,MLT,L,Ek)'
            stop
        endif

        if (eWM%Nkp /=  tauVAR%dims(1) .or. eWM%Nmlt /=  tauVAR%dims(2) .or. eWM%Nl /=  tauVAR%dims(3) .or. eWM%Ne /=  tauVAR%dims(4)) then
            write(*,*) "tauDims:",tauVAR%dims,"Nk:",eWM%Nkp,"Nm:",eWM%Nmlt,"Nl:",eWM%Nl,"Ne:",eWM%Ne
            write(*,*) 'Dimensions of tau arrays are not compatible'
            stop
        endif

        ! If still here we're good to allocate and store data

        allocate(eWM%Kp1D    (eWM%Nkp ))
        allocate(eWM%MLT1D   (eWM%Nmlt))
        allocate(eWM%L1D     (eWM%Nl  ))
        allocate(eWM%Energy1D(eWM%Ne  ))
        allocate(eWM%Tau4D(eWM%Nkp,eWM%Nmlt,eWM%Nl,eWM%Ne))

        call IOArray1DFill(IOVars,"Kp" , eWM%Kp1D    )
        call IOArray1DFill(IOVars,"MLT", eWM%MLT1D   )
        call IOArray1DFill(IOVars,"L"  , eWM%L1D     )
        call IOArray1DFill(IOVars,"Ek" , eWM%Energy1D)
        call IOArray4DFill(IOVars,"Tau", eWM%Tau4D   )

        call ClearIO(IOVars)

        !Array order check: array is in acsending order
           !Chorus
        if(eWM%Kp1D(1) > eWM%Kp1D(eWM%Nkp)) then
            write(*,*) "Kp: ",eWM%Kp1D
            write(*,*) "reorder wave model so Kp is in ascending order"
            stop
         end if

         if(eWM%MLT1D(1) > eWM%MLT1D(eWM%Nmlt)) then
            write(*,*) "MLT: ",eWM%MLT1D
            write(*,*) "reorder wave model so MLT is in ascending order"
            stop
         end if

         if(eWM%L1D(1) > eWM%L1D(eWM%Nl)) then
            write(*,*) "L: ",eWM%L1D
            write(*,*) "reorder wave model so L shell is in ascending order"
            stop
         end if

         if(eWM%Energy1D(1) > eWM%Energy1D(eWM%Ne)) then
            write(*,*) "Ek: ",eWM%Energy1D
            write(*,*) "reorder wave model so Ek is in ascending order"
            stop
         end if 

    end subroutine initEWM


end module raijuELossWM