! $Id$
!
! Earth System Modeling Framework
! Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!
!==============================================================================
!
program ESMF_FileRegridUTest

!------------------------------------------------------------------------------

#include "ESMF_Macros.inc"
#include "ESMF.h"

!==============================================================================
!BOP
! !PROGRAM: ESMF_FileRegridUTest - Test two ESMF_FileRegrid() APIs
!
! !DESCRIPTION:
!
! The code in this file tests the two ESMF_FileRegrid() APIs -- one generates a weight file and the other
! generate a route handle.
!
!-----------------------------------------------------------------------------
! !USES:
  use ESMF_TestMod     ! test methods
  use ESMF_FileRegridMod
  use ESMF_FileRegridCheckMod

  implicit none

!------------------------------------------------------------------------------
! The following line turns the CVS identifier string into a printable variable.
  character(*), parameter :: version = &
    '$Id$'
!------------------------------------------------------------------------------
    
  ! cumulative result: count failures; no failures equals "all pass"
  integer :: result = 0

  ! individual test result code
  integer :: localrc, rc, petCount,localPet
  ! individual test failure message
  character(ESMF_MAXSTR) :: name, failMsg

  logical :: correct
  type(ESMF_VM) :: vm
  character(len=256) :: srcfile, dstfile
  type(ESMF_PoleMethod_Flag) :: pole
  type(ESMF_FileFormat_Flag) :: srcFileType, dstFileType
  type(ESMF_RegridMethod_Flag) :: methodflag
  type(ESMF_UnmappedAction_Flag) :: unmappedaction
  type(ESMF_RouteHandle) :: routehandle

  rc = ESMF_SUCCESS

  !-----------------------------------------------------------------------------
  call ESMF_TestStart(ESMF_SRCLINE, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  !-----------------------------------------------------------------------------

  ! get global VM
  call ESMF_VMGetGlobal(vm, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
  !----------------------------------------------------------------------------
  !NEX_UTest
  write(name, *) "FileRegrid bilinear Test using UGRID"
  write(failmsg, *) "Did not return ESMF_SUCCESS"

  srcfile = 'data/BT42_ugrid.nc'
  dstfile = 'data/BT42_ugrid_dual.nc'

  pole = ESMF_POLEMETHOD_ALLAVG
  methodflag = ESMF_REGRIDMETHOD_BILINEAR
  unmappedaction = ESMF_UNMAPPEDACTION_ERROR

  call ESMF_FileRegrid(srcfile, dstfile, "nodedata", "nodedata", &
       dstLoc='node', regridmethod=methodflag, &
       polemethod = pole, unmappedaction = unmappedaction, &
       verboseFlag = .true., rc=rc)

! The above call needs PIO, so check output based on that
#ifdef ESMF_PIO
  call ESMF_Test(((rc.eq.ESMF_SUCCESS)), name, failMsg, result, ESMF_SRCLINE)
#else
  write(failMsg, *) "Did not return ESMF_RC_LIB_NOT_PRESENT"
  call ESMF_Test(((rc==ESMF_RC_LIB_NOT_PRESENT) .or. (rc==ESMC_RC_LIB_NOT_PRESENT)), &
       name, failMsg, result, ESMF_SRCLINE) 
#endif

  !NEX_UTest
  write(name, *) "FileRegridCheck bilinear Test using UGRID"
  write(failmsg, *) "Did not return ESMF_SUCCESS"


  ! Only call to check if PIO is present and a file has been written above
#ifdef ESMF_PIO
  call ESMF_FileRegridCheck(dstfile, "nodedata", regridmethod=methodflag, &
       rc=rc)  

  call ESMF_Test(((rc.eq.ESMF_SUCCESS)), name, failMsg, result, ESMF_SRCLINE)
#else
  ! Otherwise just pass
  call ESMF_Test((.true.),  name, failMsg, result, ESMF_SRCLINE) 
#endif

  !----------------------------------------------------------------------------
  !NEX_UTest
  write(name, *) "FileRegrid neareststod Test using UGRID files with data on node"
  write(failmsg, *) "Did not return ESMF_SUCCESS"

  methodflag = ESMF_REGRIDMETHOD_NEAREST_STOD

  call ESMF_FileRegrid(srcfile, dstfile, "nodedata", "nodedata", &
       dstLoc = "node", regridmethod=methodflag, &
       polemethod = pole, unmappedaction = unmappedaction, &
       verboseFlag = .true., rc=rc)

! The above call needs NETCDF, so check output based on that
#ifdef ESMF_NETCDF
  call ESMF_Test(((rc.eq.ESMF_SUCCESS)), name, failMsg, result, ESMF_SRCLINE)
#else
  write(failMsg, *) "Did not return ESMF_RC_LIB_NOT_PRESENT"
  call ESMF_Test(((rc==ESMF_RC_LIB_NOT_PRESENT) .or. (rc==ESMC_RC_LIB_NOT_PRESENT)), &
       name, failMsg, result, ESMF_SRCLINE) 
#endif

  !NEX_UTest
  write(name, *) "FileRegridCheck neareststod Test using UGRID with data no node"
  write(failmsg, *) "Did not return ESMF_SUCCESS"

  ! Only check file if NETCDF is present and test above ran
#ifdef ESMF_NETCDF
  call ESMF_FileRegridCheck(dstfile, "nodedata", regridmethod=methodflag, &
       rc=rc)  

  call ESMF_Test(((rc.eq.ESMF_SUCCESS)), name, failMsg, result, ESMF_SRCLINE)
#else
  ! Otherwise just pass
  call ESMF_Test((.true.),  name, failMsg, result, ESMF_SRCLINE) 
#endif

  !----------------------------------------------------------------------------
  !NEX_UTest
  write(name, *) "FileRegrid conserve using UGRID files with data on element"
  write(failmsg, *) "Did not return ESMF_SUCCESS"

  pole = ESMF_POLEMETHOD_NONE
  methodflag = ESMF_REGRIDMETHOD_CONSERVE
  unmappedaction = ESMF_UNMAPPEDACTION_IGNORE

  call ESMF_FileRegrid(srcfile, dstfile, "elmtdata", "elmtdata", &
       regridmethod=methodflag, &
       unmappedaction = unmappedaction, &
       verboseFlag = .true., rc=rc)

! The above call needs PIO, so check output based on that
#ifdef ESMF_PIO
  call ESMF_Test(((rc.eq.ESMF_SUCCESS)), name, failMsg, result, ESMF_SRCLINE)
#else
  write(failMsg, *) "Did not return ESMF_RC_LIB_NOT_PRESENT"
  call ESMF_Test(((rc==ESMF_RC_LIB_NOT_PRESENT) .or. (rc==ESMC_RC_LIB_NOT_PRESENT)), &
       name, failMsg, result, ESMF_SRCLINE) 
#endif

  !NEX_UTest
  write(name, *) "FileRegridCheck conserve Test using UGRID with data on element"
  write(failmsg, *) "Did not return ESMF_SUCCESS"

  ! Only call to check if PIO is present and a file has been written above
#ifdef ESMF_PIO
  call ESMF_FileRegridCheck(dstfile, "elmtdata", regridmethod=methodflag, &
       rc=rc)  

  call ESMF_Test(((rc.eq.ESMF_SUCCESS)), name, failMsg, result, ESMF_SRCLINE)
#else
  ! Otherwise just pass
  call ESMF_Test((.true.),  name, failMsg, result, ESMF_SRCLINE) 
#endif

  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  call ESMF_TestEnd(ESMF_SRCLINE)
  !-----------------------------------------------------------------------------

end program ESMF_FileRegridUTest
