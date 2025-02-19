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

!==============================================================================
!ESMF_MULTI_PROC_EXAMPLE        String used by test script to count examples.
!==============================================================================

program ESMF_ArrayScatterGatherEx
#include "ESMF.h"

  use ESMF
  use ESMF_TestMod
  
  implicit none
  
  ! local variables
  integer:: rc, petCount, localPet
  type(ESMF_VM):: vm
  integer:: i, j, k
  type(ESMF_DistGrid):: distgrid
  type(ESMF_Array):: array
  type(ESMF_ArraySpec):: arrayspec
  integer(ESMF_KIND_I4), allocatable:: farray(:,:,:)
  real(ESMF_KIND_R8), pointer:: myFarray2D(:,:), myFarray2D2(:,:)
  integer :: finalrc, result
  character(ESMF_MAXSTR) :: testname
  character(ESMF_MAXSTR) :: failMsg

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

  write(failMsg, *) "Example failure"
  write(testname, *) "Example ESMF_ArrayScatterGatherEx"
  
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
  finalrc = ESMF_SUCCESS
  call ESMF_Initialize(vm=vm,defaultlogfilename="ArrayScatterGatherEx.Log", &
                    logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
  if (petCount /= 4) then
    finalrc = ESMF_FAILURE
    goto 10
  endif
  
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

!BOE
!
! \subsubsection{Communication -- Scatter and Gather}
! \label{Array:ScatterGather}
! 
! It is a common situation, particularly in legacy code, that an ESMF Array
! object must be filled with data originating from a large Fortran array stored
! on a single PET.
!EOE
!BOC
  if (localPet == 0) then
    allocate(farray(10,20,30))
    do k=1, 30
      do j=1, 20
        do i=1, 10
          farray(i, j, k) = k*1000 + j*100 +  i
        enddo
      enddo
    enddo
  else
    allocate(farray(0,0,0))
  endif
!EOC
!BOC
  distgrid = ESMF_DistGridCreate(minIndex=(/1,1,1/), maxIndex=(/10,20,30/), &
    rc=rc)
!EOC  
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  call ESMF_ArraySpecSet(arrayspec, typekind=ESMF_TYPEKIND_I4, rank=3, rc=rc)
!EOC  
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  array = ESMF_ArrayCreate(arrayspec=arrayspec, distgrid=distgrid, rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOE
! The {\tt ESMF\_ArrayScatter()} method provides a convenient way of scattering
! array data from a single root PET across the DEs of an ESMF Array object.
!EOE
!BOC
  call ESMF_ArrayScatter(array, farray=farray, rootPet=0, rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  deallocate(farray)
!EOC
!BOE
! The destination of the ArrayScatter() operation are all the DEs of a single
! tile. For multi-tile Arrays the destination tile can be specified. The 
! shape of the scattered Fortran array must match the shape of the destination
! tile in the ESMF Array.
!
! Gathering data decomposed and distributed across the DEs of an ESMF Array
! object into a single Fortran array on root PET is accomplished by calling
! {\tt ESMF\_ArrayGather()}.
!EOE
!BOC
  if (localPet == 3) then
    allocate(farray(10,20,30))
  else
    allocate(farray(0,0,0))
  endif
  
  call ESMF_ArrayGather(array, farray=farray, rootPet=3, rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  deallocate(farray)
!EOC
!BOE
! The source of the ArrayGather() operation are all the DEs of a single
! tile. For multi-tile Arrays the source tile can be specified. The 
! shape of the gathered Fortran array must match the shape of the source
! tile in the ESMF Array.
!EOE
  call ESMF_ArrayDestroy(array, rc=rc) ! destroy the Array object
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_DistGridDestroy(distgrid, rc=rc) ! destroy the DistGrid object
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
!BOE
! The {\tt ESMF\_ArrayScatter()} operation allows to fill entire replicated
! Array objects with data coming from a single root PET.
!EOE
!BOC
  distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/5,5/), &
    regDecomp=(/2,3/), rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  call ESMF_ArraySpecSet(arrayspec, typekind=ESMF_TYPEKIND_R8, rank=2, rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  array = ESMF_ArrayCreate(arrayspec=arrayspec, distgrid=distgrid, &
    distgridToArrayMap=(/0,0/), undistLBound=(/11,21/), &
    undistUBound=(/14,22/), rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOE
! The shape of the Fortran source array used in the Scatter() call must be
! that of the contracted Array, i.e. contracted DistGrid dimensions do not
! count. For the {\tt array} just created this means that the source array
! on {\tt rootPet} must be of shape 4 x 2.
!EOE
!BOC
  if (localPet == 0) then
    allocate(myFarray2D(4,2))
    do j=1,2
      do i=1,4
        myFarray2D(i,j) = i * 100.d0 + j * 1.2345d0 ! initialize
      enddo
    enddo
  else
    allocate(myFarray2D(0,0))
  endif
  
  call ESMF_ArrayScatter(array, farray=myFarray2D, rootPet=0, rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  deallocate(myFarray2D)
!EOC
!BOE
! This will have filled each local 4 x 2 Array piece with the replicated
! data of {\tt myFarray2D}.
!EOE
!  call ESMF_ArrayPrint(array, rc=rc)
!  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  call ESMF_ArrayDestroy(array, rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  call ESMF_DistGridDestroy(distgrid, rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
  
!!!!!!!!!!  
  
!BOE
! As a second example for the use of Scatter() and Gather() consider the
! following replicated Array created from existing local Fortran arrays.
!EOE
!BOC
  allocate(myFarray2D(3,10))
  distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/40,10/), rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  array = ESMF_ArrayCreate(farray=myFarray2D, distgrid=distgrid, &
    indexflag=ESMF_INDEX_DELOCAL, distgridToArrayMap=(/0,2/), rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOE
! The {\tt array} object associates the 2nd DistGrid dimension with the 2nd
! Array dimension. The first DistGrid dimension is not associated with any
! Array dimension and will lead to replication of the Array along the DEs of
! this direction. Still, the local arrays that comprise the {\tt array} 
! object refer to independent pieces of memory and can be initialized 
! independently.
!EOE
!BOC
  myFarray2D = localPet ! initialize
!EOC
!  call ESMF_ArrayPrint(array, rc=rc)
!  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOE
! However, the notion of replication becomes visible when an array of shape
! 3 x 10 on root PET 0 is scattered across the Array object.
!EOE
!BOC
  if (localPet == 0) then
    allocate(myFarray2D2(5:7,11:20))
  
    do j=11,20
      do i=5,7
        myFarray2D2(i,j) = i * 100.d0 + j * 1.2345d0 ! initialize
      enddo
    enddo
  else
    allocate(myFarray2D2(0,0))
  endif
  
  call ESMF_ArrayScatter(array, farray=myFarray2D2, rootPet=0, rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  deallocate(myFarray2D2)
!EOC
!BOE
! The Array pieces on every DE will receive the same source data, resulting
! in a replication of data along DistGrid dimension 1.
!EOE
!  call ESMF_ArrayPrint(array, rc=rc)
!  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOE
! When the inverse operation, i.e. {\tt ESMF\_ArrayGather()}, is applied to
! a replicated Array an intrinsic ambiguity needs to be considered. ESMF 
! defines the gathering of data of a replicated Array as the collection of data
! originating from the numerically higher DEs. This means that data in
! replicated elements associated with numerically lower DEs will be ignored
! during {\tt ESMF\_ArrayGather()}. For the current example this means that
! changing the Array contents on PET 1, which here corresponds to DE 1,
!EOE
!BOC
  if (localPet == 1) then
    myFarray2D = real(1.2345, ESMF_KIND_R8)
  endif
!EOC
!BOE
! will {\em not} affect the result of
!EOE
  call ESMF_ArrayPrint(array, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  allocate(myFarray2D2(3,10))
  myFarray2D2 = 0.d0    ! initialize to a known value
  call ESMF_ArrayGather(array, farray=myFarray2D2, rootPet=0, rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  print *, "localPet = ", localPet, "myFarray2D2: ", myFarray2D2
  
!BOE
! The result remains completely defined by the unmodified values of Array in 
! DE 3, the numerically highest DE. However, overriding the DE-local Array
! piece on DE 3
!EOE
!BOC
  if (localPet==3) then
    myFarray2D = real(5.4321, ESMF_KIND_R8)
  endif
!EOC
  call ESMF_ArrayPrint(array, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOE
! will change the outcome of
!EOE
!BOC
  call ESMF_ArrayGather(array, farray=myFarray2D2, rootPet=0, rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOE
! as expected.
!EOE
  print *, "localPet = ", localPet, "myFarray2D2: ", myFarray2D2
!BOC
  deallocate(myFarray2D2)

  call ESMF_ArrayDestroy(array, rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!BOC
  call ESMF_DistGridDestroy(distgrid, rc=rc)
!EOC
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  deallocate(myFarray2D)
  
  
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
10 continue

  ! IMPORTANT: ESMF_STest() prints the PASS string and the # of processors in the log
  ! file that the scripts grep for.
  call ESMF_STest((finalrc.eq.ESMF_SUCCESS), testname, failMsg, result, ESMF_SRCLINE)

  call ESMF_Finalize(rc=rc)
  
  if (rc/=ESMF_SUCCESS) finalrc = ESMF_FAILURE
  if (finalrc==ESMF_SUCCESS) then
    print *, "PASS: ESMF_ArrayScatterGatherEx.F90"
  else
    print *, "FAIL: ESMF_ArrayScatterGatherEx.F90"
  endif
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------
  
end program
