program Kaimag
   use clocks
   use xml_input
   use kaimagio
   use files
   use strings
   use kaimagmain

   implicit none

   type(kaimagApp_T) :: kaimagApp, kaimagIC
   character(len=strLen) :: optFilename
   procedure(kaimagupdateV_T), pointer :: updateV => NULL()
   write(*,*) 'Beware the Big Bad Wolf!'

   call initClocks()
   call Tic("Omega")

   call Tic("Init")

   call getIDeckStr(optFilename)
   call init_kaimag(kaimagApp,optFilename,updateV)
   call Toc("Init")

   call Tic("Loop")

   kaimagApp%oState = kaimagApp%State
   kaimagIC%State = kaimagApp%State
   ! Write the initial conditions
   if (kaimagApp%Model%IO%doOutput(kaimagApp%Model%t)) then 
   call fOutputKaimag(kaimagApp%Model,kaimagApp%Grid,kaimagApp%State)
   endif

   if (kaimagApp%Model%IO%doRestart(kaimagApp%Model%t)) then
   call resOutputKaimag(kaimagApp%Model,kaimagApp%Grid,kaimagApp%State)
   endif 

   !Now do the main computational loop
   do while (kaimagApp%Model%t < kaimagApp%Model%tFin)

   call stepKaimag(kaimagApp,updateV)

   if (kaimagApp%Model%IO%doConsole(kaimagApp%Model%ts)) then 
     call consoleOutputKaimag(kaimagApp%Model,kaimagApp%Grid,kaimagApp%State)
   endif

   if (kaimagApp%Model%IO%doOutput(kaimagApp%Model%t)) then 
     call fOutputKaimag(kaimagApp%Model,kaimagApp%Grid,kaimagApp%State)
   endif

   if (kaimagApp%Model%IO%doRestart(kaimagApp%Model%t)) then
     call resOutputKaimag(kaimagApp%Model,kaimagApp%Grid,kaimagApp%State)
   endif

   enddo ! while loop

   call calcError(kaimagApp%Model,kaimagApp%Grid,kaimagApp%State,kaimagIC%State)

   call Toc("Loop")
   call Toc("Omega")
   call printClocks()


end program Kaimag
