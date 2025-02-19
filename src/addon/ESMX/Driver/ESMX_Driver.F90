!==============================================================================
! ESMX_Driver (Earth System Model eXecutable Driver)
!   This file contains the ESMX top level driver.
!==============================================================================
#define FILENAME "src/addon/ESMX/Driver/ESMX_Driver.F90"
!==============================================================================

module ESMX_Driver

  use ESMF
  use NUOPC
  use NUOPC_Driver,             driverSS    => SetServices

  include "compUse.inc"

  implicit none

  private

  public SetServices, SetVM, HConfigCreateFoundNode

#if defined (__NVCOMPILER)
!TODO: remove once NVHPC and PGI compilers work correctly w/o work-around

  abstract interface
    recursive subroutine SetServicesRoutine(gridcomp, rc)
      use ESMF
      implicit none
      type(ESMF_GridComp)        :: gridcomp ! must not be optional
      integer, intent(out)       :: rc       ! must not be optional
    end subroutine
    recursive subroutine SetVMRoutine(gridcomp, rc)
      use ESMF
      implicit none
      type(ESMF_GridComp)        :: gridcomp ! must not be optional
      integer, intent(out)       :: rc       ! must not be optional
    end subroutine
  end interface

  type type_CompDef
    procedure(SetServicesRoutine), pointer, nopass :: ssPtr => null()
    procedure(SetVMRoutine),       pointer, nopass :: svPtr => null()
    character(ESMF_MAXSTR)                  :: name = "__uninitialized__"
  end type

#else
  type type_CompDef
    procedure(SetServices), pointer, nopass :: ssPtr => null()
    procedure(SetVM),       pointer, nopass :: svPtr => null()
    character(ESMF_MAXSTR)                  :: name = "__uninitialized__"
  end type
#endif

  include "compCnt.inc"

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    ! Derive from NUOPC_Driver
    call NUOPC_CompDerive(driver, driverSS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

    ! Specialize SetModelServices
    call NUOPC_CompSpecialize(driver, specLabel=label_SetModelServices, &
      specRoutine=SetModelServices, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

    ! Specialize SetRunSequence
    call NUOPC_CompSpecialize(driver, specLabel=label_SetRunSequence, &
      specRoutine=SetRunSequence, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine SetModelServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    integer                         :: i, j, componentCount, ompNumThreads
    integer, allocatable            :: petList(:), devList(:)
    type(ESMF_GridComp)             :: comp
    type(ESMF_HConfig)              :: hconfig, hconfigNode, hconfigNode2
    character(:), allocatable       :: configKey(:)
    character(:), allocatable       :: componentList(:)
    character(:), allocatable       :: compLabel
    character(:), allocatable       :: model
    character(:), allocatable       :: string1, string2
    type(ESMF_Info)                 :: info
    type(type_CompDef), allocatable :: CompDef(:)
    logical                         :: inCompDef, isFlag

    rc = ESMF_SUCCESS

    ! Look for config in the component
    call ESMF_GridCompGet(driver, hconfigIsPresent=isFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

    if (isFlag) then
      ! Get hconfig from component
      call ESMF_GridCompGet(driver, hconfig=hconfig, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
    else
      ! Attempt to open hconfig from default file "esmxRun.yaml", set on Comp
      call ESMF_GridCompSet(driver, configFile="esmxRun.yaml", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      ! Get hconfig from component
      call ESMF_GridCompGet(driver, hconfig=hconfig, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      ! Find hconfig node that holds driver level attributes, conditionally ingest
      configKey = ["ESMX      ", "Driver    ", "attributes"]
      hconfigNode2 = HConfigCreateFoundNode(hconfig, configKey=configKey, &
        foundFlag=isFlag, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (isFlag) then
        call NUOPC_CompAttributeIngest(driver, hconfigNode2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      endif
      call ESMF_HConfigDestroy(hconfigNode2, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
    endif

    ! Find hconfigNode that holds driver level settings according to configKey
    configKey = ["ESMX  ", "Driver"]
    hconfigNode = HConfigCreateFoundNode(hconfig, configKey=configKey, &
      foundFlag=isFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out
    if (isFlag) then
      ! Validate hconfigNode against ESMX/Driver controlled key vocabulary
      isFlag = ESMF_HConfigValidateMapKeys(hconfigNode, &
        vocabulary=["attributes   ", &  ! ESMX_Driver option
                    "componentList", &  ! ESMX_Driver option
                    "runSequence  "  &  ! ESMX_Driver option
                    ], badKey=string1, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      if (.not.isFlag) then
        call ESMF_LogSetError(ESMF_RC_ARG_WRONG, &
          msg="An invalid key was found in config under ESMX/Driver (maybe a typo?): "//string1, &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif
      ! Ingest the generic component label list
      isFlag = ESMF_HConfigIsDefined(hconfigNode, keyString="componentList", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      isFlag = isFlag .and. &
        .not.ESMF_HConfigIsNull(hconfigNode, keyString="componentList", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (isFlag) then
        componentList = ESMF_HConfigAsStringSeq(hconfigNode, stringLen=32, &
          keyString="componentList", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      endif
      call ESMF_HConfigDestroy(hconfigNode, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
    endif

    ! Determine componentCount
    componentCount = 0  ! no components for unallocated componentList
    if (allocated(componentList)) componentCount = size(componentList)

    ! Setup CompDef structure
    allocate(CompDef(componentDefCount))
    include "compDef.inc"

    ! Determine information for each component and add to the driver
    do i=1, componentCount
      ! compLabel
      compLabel=trim(componentList(i))

      ! Find hconfigNode that holds component level settings
      hconfigNode = HConfigCreateFoundNode(hconfig, configKey=[compLabel], &
        foundFlag=isFlag, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (.not.isFlag) then
        call ESMF_LogSetError(ESMF_RC_ARG_INCOMP, &
          msg="Must provide settings for component: "//compLabel, &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return  ! bail out
      endif

      ! Cannot validate hconfigNode here, because valid keys depend on model.
      ! Must rely on model implementation to do the key validation!

      ! Set model
      isFlag = ESMF_HConfigIsDefined(hconfigNode, keyString="model", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (.not.isFlag) then
        call ESMF_LogSetError(ESMF_RC_ARG_INCOMP, &
          msg="Must provide `model` setting for component: "//compLabel, &
          line=__LINE__, file=FILENAME, rcToReturn=rc)
        return  ! bail out
      endif
      model = ESMF_HConfigAsString(hconfigNode, keyString="model", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out

      ! set up petList
      isFlag = ESMF_HConfigIsDefined(hconfigNode, keyString="petList", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (isFlag) then
        hconfigNode2 = ESMF_HConfigCreateAt(hconfigNode, keyString="petList", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call NUOPC_IngestPetList(petList, hconfigNode2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_HConfigDestroy(hconfigNode2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      else
        allocate(petList(0))
      endif

      ! set up devList
      isFlag = ESMF_HConfigIsDefined(hconfigNode, keyString="devList", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (isFlag) then
        hconfigNode2 = ESMF_HConfigCreateAt(hconfigNode, keyString="devList", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call NUOPC_IngestPetList(devList, hconfigNode2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_HConfigDestroy(hconfigNode2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      else
        allocate(devList(0))
      endif

      ! Set NUOPC hint for OpenMP
      info = ESMF_InfoCreate(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      isFlag = ESMF_HConfigIsDefined(hconfigNode, &
        keyString="ompNumThreads", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (isFlag) then
        ompNumThreads = ESMF_HConfigAsI4(hconfigNode, &
          keyString="ompNumThreads", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_InfoSet(info, key="/NUOPC/Hint/PePerPet/MaxCount", &
          value=ompNumThreads, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      endif

      ! Search for an entry for this component model inside CompDef
      inCompDef = .false.
      do j=1, componentDefCount
        if (trim(CompDef(j)%name)=="__uninitialized__") exit
        ! case insensitive string comparison
        string1 = ESMF_UtilStringLowerCase(trim(CompDef(j)%name), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        string2 = ESMF_UtilStringLowerCase(trim(model), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        if (string1==string2) then
          inCompDef = .true.
          exit
        endif
      enddo

      if (inCompDef) then
        ! add child component with SetVM and SetServices in CompDef
#if defined (__INTEL_LLVM_COMPILER) || defined (__NVCOMPILER) || defined (NAGFOR)
!TODO: remove once IFX, NVHPC, and NAG compilers work correctly w/o work-around
        call NUOPC_DriverAddGridCompPtr(driver, trim(compLabel), hconfig=hconfig, &
          compSetServicesRoutine=CompDef(j)%ssPtr, compSetVMRoutine=CompDef(j)%svPtr, &
          info=info, petList=petList, devList=devList, comp=comp, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, &
          msg="Unable to add component '"//trim(compLabel)// &
            "' to driver via Fortran module.", &
          line=__LINE__, file=FILENAME)) return  ! bail out
#else
        call NUOPC_DriverAddComp(driver, trim(compLabel), hconfig=hconfig, &
          compSetServicesRoutine=CompDef(j)%ssPtr, compSetVMRoutine=CompDef(j)%svPtr, &
          info=info, petList=petList, devList=devList, comp=comp, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, &
          msg="Unable to add component '"//trim(compLabel)// &
            "' to driver via Fortran module.", &
          line=__LINE__, file=FILENAME)) return  ! bail out
#endif
      else
        ! add child component with SetVM and SetServices in shared object
        call NUOPC_DriverAddComp(driver, trim(compLabel), hconfig=hconfig, &
          sharedObj=trim(model), info=info, petList=petList, devList=devList, &
          comp=comp, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, &
          msg="Unable to add component '"//trim(compLabel)// &
            "' to driver via shared object: "//trim(model), &
          line=__LINE__, file=FILENAME)) return  ! bail out
      endif

      ! Find hconfig node that holds model level attributes, conditionally ingest
      isFlag = ESMF_HConfigIsDefined(hconfigNode, keyString="attributes", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (isFlag) then
        hconfigNode2 = ESMF_HConfigCreateAt(hconfigNode, &
          keyString="attributes", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call NUOPC_CompAttributeIngest(comp, hconfigNode2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_HConfigDestroy(hconfigNode2, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      endif

      ! clean-up
      deallocate(petList, devList)
      call ESMF_InfoDestroy(info, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      call ESMF_HConfigDestroy(hconfigNode, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
    enddo

    ! Ensure clock is set (by parent)
    call ESMF_GridCompGet(driver, clockIsPresent=isFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out
    if (.not.isFlag) then
      call ESMF_LogSetError(ESMF_RC_ARG_INCOMP, &
        msg="Parent must set the driver clock!", &
        line=__LINE__, file=FILENAME, rcToReturn=rc)
      return  ! bail out
    endif

  end subroutine SetModelServices

  !-----------------------------------------------------------------------------

  subroutine SetRunSequence(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_HConfig)              :: hconfig, hconfigNode
    character(:), allocatable       :: configKey(:)
    logical                         :: isFlag

    rc = ESMF_SUCCESS

    ! Query the driver for hconfig
    call ESMF_GridCompGet(driver, hconfig=hconfig, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out

    ! Find hconfigNode that holds driver level settings according to configKey
    configKey = ["ESMX  ", "Driver"]
    hconfigNode = HConfigCreateFoundNode(hconfig, configKey=configKey, &
      foundFlag=isFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) return  ! bail out
    if (isFlag) then
      isFlag = ESMF_HConfigIsDefined(hconfigNode, keyString="runSequence", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
      if (isFlag) then
        hconfig = ESMF_HConfigCreateAt(hconfigNode, keyString="runSequence", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call NUOPC_DriverIngestRunSequence(driver, hconfig, &
          autoAddConnectors=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, &
          msg="Unable to ingest RunSequence!", &
          line=__LINE__, file=FILENAME)) return  ! bail out
        call ESMF_HConfigDestroy(hconfig, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, file=FILENAME)) return  ! bail out
      endif
      ! clean-up
      call ESMF_HConfigDestroy(hconfigNode, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) return  ! bail out
    endif

  end subroutine SetRunSequence

  !-----------------------------------------------------------------------------

  function HConfigCreateFoundNode(hconfig, configKey, foundFlag, rc)
    type(ESMF_HConfig)    :: HConfigCreateFoundNode
    type(ESMF_HConfig)    :: hconfig
    character(*)          :: configKey(:)
    logical, intent(out)  :: foundFlag
    integer, intent(out)  :: rc

    ! local variables
    type(ESMF_HConfig)    :: hconfigNodePrev
    integer               :: i
    logical               :: isFlag

    rc = ESMF_SUCCESS

    HConfigCreateFoundNode = ESMF_HConfigCreate(hconfig, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=FILENAME)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    foundFlag = .true.
    do i=1, size(configKey)
      isFlag = ESMF_HConfigIsMap(HConfigCreateFoundNode, &
        keyString=configKey(i), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      if (.not.isFlag) then
        ! configKey must be a map
        foundFlag = .false.
        exit  ! break out of loop
      endif
      hconfigNodePrev = HConfigCreateFoundNode
      HConfigCreateFoundNode = ESMF_HConfigCreateAt(hconfigNodePrev, &
        keyString=configKey(i),rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_HConfigDestroy(hconfigNodePrev, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=FILENAME)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
    enddo

  end function

  !-----------------------------------------------------------------------------

end module ESMX_Driver
