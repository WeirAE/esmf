! $Id$
!
! Example/test code which shows User Component calls.

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!BOP
!
! !DESCRIPTION:
!  User-supplied Coupler
!
!
!\begin{verbatim}

    module user_coupler

    ! ESMF Framework module
    use ESMF
    
    implicit none
    
    public usercpl_register
        
    ! global data
    type(ESMF_RouteHandle), save :: routehandle

    contains

!-------------------------------------------------------------------------
!   !  The Register routine sets the subroutines to be called
!   !   as the init, run, and finalize routines.  Note that these are
!   !   private to the module.
 
    subroutine usercpl_register(comp, rc)
      type(ESMF_CplComp) :: comp
      integer, intent(out) :: rc

      rc = ESMF_SUCCESS
      print *, "in user setservices routine"

      ! Register the callback routines.
      call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, user_init, rc=rc)
      if(rc/=ESMF_SUCCESS) return
      call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_RUN, user_run, rc=rc)
      if(rc/=ESMF_SUCCESS) return
      call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, user_final, rc=rc)
      if(rc/=ESMF_SUCCESS) return

      print *, "Registered Initialize, Run, and Finalize routines"

    end subroutine

!-------------------------------------------------------------------------
!   !User Comp Component created by higher level calls, here is the
!   ! Initialization routine.
 
    
    subroutine user_init(comp, importState, exportState, clock, rc)
      type(ESMF_CplComp) :: comp
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      ! Local variables
      integer :: itemcount
      type(ESMF_Field) :: humidity1, humidity2
      type(ESMF_VM) :: vm
      integer(ESMF_KIND_I4), allocatable          :: factorList(:)
      integer, allocatable                        :: factorIndexList(:,:)
      integer       :: lpe

      rc = ESMF_SUCCESS
      print *, "User Coupler Init starting"

      call ESMF_CplCompGet(comp, vm=vm, rc=rc)
      call ESMF_VMGet(vm, localPET=lpe, rc=rc)

      ! uncomment the following when locstream supports reconcile

      ! Since the components we are coupling between are running concurrently,
      ! they have each separately created ESMF objects.   We are planning to
      ! use a communications call (SMM) here, so first we must make a new
      ! call to reconcile the object lists in all the import and export states.

      call ESMF_StateReconcile(importState, vm=vm, rc=rc)
      if(rc/=ESMF_SUCCESS) return

      call ESMF_StateReconcile(exportState, vm=vm, rc=rc)
      if(rc/=ESMF_SUCCESS) return

      call ESMF_StateGet(importState, itemcount=itemcount, rc=rc)
      if(rc/=ESMF_SUCCESS) return
      print *, "Import State contains ", itemcount, " items."
      call ESMF_StateGet(exportState, itemcount=itemcount, rc=rc)
      if(rc/=ESMF_SUCCESS) return
      print *, "Export State contains ", itemcount, " items."
       
      ! Get input data
      call ESMF_StateGet(importState, "humidity", humidity1, rc=rc)
      if(rc/=ESMF_SUCCESS) return
      ! call ESMF_FieldPrint(humidity1, rc=rc)

      ! Get location of output data
      call ESMF_StateGet(exportState, "humidity", humidity2, rc=rc)
      if(rc/=ESMF_SUCCESS) return
      ! call ESMF_FieldPrint(humidity2, rc=rc)

      ! Get VM from coupler component to send down to SMM Store routine
      call ESMF_CplCompGet(comp, vm=vm, rc=rc)
      if(rc/=ESMF_SUCCESS) return

      ! Call ESMF_FieldRegridStore() to calculate routehandle between different locstreams
      call ESMF_FieldRegridStore(srcField=humidity1, srcMaskValues=(/1/), &
           dstField=humidity2, dstMaskValues=(/2/), &
           regridmethod=ESMF_REGRIDMETHOD_NEAREST_STOD, &
           routeHandle=routehandle, rc=rc)
      if(rc/=ESMF_SUCCESS) return

      print *, "User Coupler Init returning"
   
    end subroutine user_init


!-------------------------------------------------------------------------
!   !  The Run routine where data is coupled.
!   !
 
    subroutine user_run(comp, importState, exportState, clock, rc)
      type(ESMF_CplComp) :: comp
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      ! Local variables
      type(ESMF_Field) :: humidity1, humidity2

      rc = ESMF_SUCCESS
      print *, "User Coupler Run starting"

      ! Get input data
      call ESMF_StateGet(importState, "humidity", humidity1, rc=rc)
      if(rc/=ESMF_SUCCESS) return
      ! call ESMF_FieldPrint(humidity1, rc=rc)

      ! Get location of output data
      call ESMF_StateGet(exportState, "humidity", humidity2, rc=rc)
      if(rc/=ESMF_SUCCESS) return
      ! call ESMF_FieldPrint(humidity2, rc=rc)

      ! Call ESMF_FieldRegrid() to move  data from src locstream to destination locstream
      call ESMF_FieldRegrid(humidity1, humidity2, routehandle, rc=rc)
      if(rc/=ESMF_SUCCESS) return

      ! Data is moved directly to the field in the output state, so no
      ! "put" is needed here.
 
      print *, "User Coupler Run returning"

    end subroutine user_run


!-------------------------------------------------------------------------
!   !  The Finalization routine where things are deleted and cleaned up.
!   !
 
    subroutine user_final(comp, importState, exportState, clock, rc)
      type(ESMF_CplComp) :: comp
      type(ESMF_State) :: importState, exportState
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      ! Local variables

      rc = ESMF_SUCCESS
      print *, "User Coupler Final starting"
   
      ! Release resources stored for the SMM
      call ESMF_FieldRegridRelease(routehandle, rc=rc)
      if(rc/=ESMF_SUCCESS) return

      print *, "User Coupler Final returning"

    end subroutine user_final


    end module user_coupler
    
!\end{verbatim}
    
