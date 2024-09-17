submodule (volttypes) raijuCplTypesSub

    use raijjuCplHelpers


    module subroutine raiCplInitModel(App, xml)
        class(raijuCoupler_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml
    end subroutine raiCplInitModel


    module subroutine raiCplInitIO(App, xml)
        class(raijuCoupler_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml
    end subroutine raiCplInitIO
    
end submodule raijuCplTypesSub