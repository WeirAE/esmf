! $Id$
!
! Earth System Modeling Framework
! Copyright 2002-2021, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!
!==============================================================================
!
program ESMF_AltGridCreateUTest

!------------------------------------------------------------------------------

#include "ESMF.h"
#include "ESMF_Macros.inc"

!==============================================================================
!BOP
! !PROGRAM: ESMF_AltGridCreateTest - Check Grid Create Routines
!
! !DESCRIPTION:
!
! The code in this file drives F90 Grid Create unit tests.
!
!-----------------------------------------------------------------------------
! !USES:
  use ESMF_TestMod     ! test methods
  use ESMF

  implicit none

!------------------------------------------------------------------------------
! The following line turns the CVS identifier string into a printable variable.
  character(*), parameter :: version = &
    '$Id$'
!------------------------------------------------------------------------------
    
  ! cumulative result: count failures; no failures equals "all pass"
  integer :: result = 0

  ! individual test result code
  integer :: localrc, rc, localPet, petCount
  logical :: correct, xercesNotPresent
  type(ESMF_TypeKind_Flag) :: typekind

  ! individual test failure message
  character(ESMF_MAXSTR) :: failMsg
  character(ESMF_MAXSTR) :: name, grid_name

  type(ESMF_Grid) :: grid, grid2, gridAlias,grid_multi
  type(ESMF_VM) :: vm
  type(ESMF_DistGrid) :: distgrid, distgrid2, distgrid_multi
  type(ESMF_Array) :: array
  integer :: coordDimMap(2,2), dimCount, undistLBound(3), undistUBound(3)
  type(ESMF_Index_Flag) :: indexflag
  integer :: distgridToGridMap(2), coordDimCount(2)
  integer :: distgridToArrayMap(3)
  integer :: coordDimCount2(3),coordDimMap2(3,3)
  integer :: gridEdgeLWidth(3),gridEdgeUWidth(3),gridAlign(3)
  integer :: exlbnd(3),exubnd(3)
  integer :: clbnd(3),cubnd(3)
  character, pointer :: buf(:)
  real(ESMF_KIND_R8), pointer :: fptr2D(:,:)
  integer :: bufCount, offset, localDECount, rank, i1,i2,lDE, i, j
  type(ESMF_StaggerLoc)          :: staggerloc8
  integer :: minIndex(3), maxIndex(3) 
  integer :: celw(3),ceuw(3)
  logical :: isLBound(2),isUBound(2)
  integer :: petMap2D(2,2,1)
  integer :: petMap2x3(2,3,1)
  real(ESMF_KIND_R8), pointer :: fptr(:,:), fptr1(:,:), fptr2(:,:)
  real(ESMF_KIND_R4), pointer :: fptr1R4(:,:), fptr2R4(:,:)
  logical:: gridBool
  logical:: isCreated
  type(ESMF_GridStatus_Flag) :: status
  ! test the AttributeGet for Grid info
  type(ESMF_TypeKind_Flag) :: attrValue
  type(ESMF_CoordSys_Flag) :: coordSys
  integer :: minIndexPTile(2,4), maxIndexPTile(2,4)
  integer :: regDecompPTile4(2,4)
  type(ESMF_Decomp_Flag) :: decompPTile(2,4)
  integer :: tile
  integer(ESMF_KIND_I4) :: regDecompPTile(2,6), deLabelList(6)
  type(ESMF_Decomp_Flag) :: decompFlagPTile(2,6)
  real(ESMF_KIND_R8), pointer :: lonPtrR8(:,:), latPtrR8(:,:)
  real(ESMF_KIND_R4), pointer :: lonPtrR4(:,:), latPtrR4(:,:)
  real(ESMF_KIND_R8), allocatable :: lonDiff(:,:), latDiff(:,:)
  real(ESMF_KIND_R8), allocatable :: mean(:)
  real(ESMF_KIND_R8) :: lonmin, latmin, lonmax, latmax, lonmean, latmean
  real(ESMF_KIND_R8) :: threshhold
  type(ESMF_DELayout) :: delayout
  integer :: total, s,  decount, localDe
  type(ESMF_Staggerloc) :: staggerLocList(2)
  type(ESMF_CubedSphereTransform_Args) :: transformArgs
  character(100) :: mosaicIn   
 
  !-----------------------------------------------------------------------------
  call ESMF_TestStart(ESMF_SRCLINE, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  !-----------------------------------------------------------------------------

  ! get global VM
  call ESMF_VMGetGlobal(vm, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! prepare DistGrid
  distgrid=ESMF_DistGridCreate(minIndex=(/1,1/),maxIndex=(/10,10/), rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  !-----------------------------------------------------------------------------
  
  !NEX_UTest
  write(name, *) "Testing ATM GridCreateMosaic with different CoordTypeKind "
  write(failMsg, *) "Returns incorrect results"

  rc = ESMF_SUCCESS
  write(mosaicIn, *) "./data/C48_mosaic.nc"

  staggerLocList(1) = ESMF_STAGGERLOC_CENTER
  staggerLocList(2) = ESMF_STAGGERLOC_CORNER
  threshhold = 1.0E-5
  ! Create cubed sphere grid with coordTypeKind == ESMF_TYPEKIND_R4
  print *, 'Reading R4 file ', mosaicIn
  grid = ESMF_GridCreateMosaic(filename=mosaicIn, &
                staggerLocList= staggerLocList, &
                coordTypeKind = ESMF_TYPEKIND_R4, &
                tileFilePath='./data/', rc=localrc)

#ifndef ESMF_NETCDF
  write(failMsg, *) "Did not return ESMF_RC_LIB_NOT_PRESENT"
  call ESMF_Test((localrc==ESMF_RC_LIB_NOT_PRESENT), name, failMsg, result, ESMF_SRCLINE) 
#else
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  call ESMF_GridGet(grid, distgrid = distgrid, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  call ESMF_DistGridGet(distgrid, delayout = delayout, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
 
  call ESMF_DELayoutGet(delayout, localDeCount = decount, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  ! Create cubed sphere grid with coordTypeKind == ESMF_TYPEKIND_R8
  print *, 'Reading R8 file ', mosaicIn
  !call ESMF_MeshWrite(mosaicIn, 'ATMmosaicIn.vtk', rc=localrc)
  grid2 = ESMF_GridCreateMosaic(filename=mosaicIn, &
                staggerLocList= staggerLocList, &
                coordTypeKind = ESMF_TYPEKIND_R8, &
                tileFilePath='./data/', rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
  
  !call ESMF_MeshWrite(grid2, 'ATMmosaic2grid.vtk', rc=localrc)
  
  do s = 1, 2
    do localDe = 0, decount-1  
      call ESMF_GridGetCoord(grid2, coordDim=1, localDe=localDe, &
         staggerloc=staggerLocList(s), farrayPtr=lonPtrR8, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
      call ESMF_GridGetCoord(grid2, coordDim=2, localDe=localDe, &
         staggerloc=staggerLocList(s), farrayPtr=latPtrR8, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
      call ESMF_GridGetCoord(grid, coordDim=1, localDe=localDe, &
       staggerloc=staggerLocList(s), farrayPtr=lonPtrR4, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
      call ESMF_GridGetCoord(grid, coordDim=2, localDe=localDe, &
         staggerloc=staggerLocList(s), farrayPtr=latPtrR4, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

      allocate(lonDiff(size(lonPtrR8,1), size(lonPtrR8,2)))
      allocate(latDiff(size(lonPtrR8,1), size(lonPtrR8,2)))
      total = size(lonPtrR8,1)*size(lonPtrR8,2)
      lonDiff = abs(lonPtrR8-lonPtrR4)
      latDiff = abs(latPtrR8-latPtrR4)

      ! Find the max/min/mean errors
      lonmin = minval(lonDiff)
      latmin = minval(latDiff)
      lonmax = maxval(lonDiff)
      latmax = maxval(latDiff)
      allocate(mean(size(lonDiff,2)))
      do  j=1, size(lonDiff,2)
         mean(j) = sum(lonDiff(:,j))
      enddo
      lonmean = sum(mean)/total
      do j=1, size(latDiff,2)
         mean(j) = sum(latDiff(:,j))
      enddo
      latmean = sum(mean)/total
      
      deallocate(lonDiff, latDiff, mean)

      print *, localPet, localDe, 'min/max/mean:', lonmin, latmin, lonmax, latmax, lonmean, latmean
      if (lonmean > threshhold .or. latmean > threshhold) rc = ESMF_FAILURE
    enddo
  enddo

  ! destroy grid
  call ESMF_GridDestroy(grid,rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
  call ESMF_GridDestroy(grid2,rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  call ESMF_Test((rc==ESMF_SUCCESS), name, failMsg, result, ESMF_SRCLINE)
#endif

 !-----------------------------------------------------------------------------
  
  !NEX_UTest
  write(name, *) "Testing ocean GridCreateMosaic with different CoordTypeKind "
  write(failMsg, *) "Returns incorrect results"

  rc = ESMF_SUCCESS
  write(mosaicIn, *) "./data/ocean_mosaic.nc"

  staggerLocList(1) = ESMF_STAGGERLOC_CENTER
  staggerLocList(2) = ESMF_STAGGERLOC_CORNER
  threshhold = 1.0E-5
  ! Create cubed sphere grid with coordTypeKind == ESMF_TYPEKIND_R4
  print *, 'Reading file ', mosaicIn
  grid = ESMF_GridCreateMosaic(filename=mosaicIn, &
                staggerLocList= staggerLocList, &
                coordTypeKind = ESMF_TYPEKIND_R4, &
                tileFilePath='./data/', rc=localrc)

#ifndef ESMF_NETCDF
  write(failMsg, *) "Did not return ESMF_RC_LIB_NOT_PRESENT"
  call ESMF_Test((localrc==ESMF_RC_LIB_NOT_PRESENT), name, failMsg, result, ESMF_SRCLINE) 
#else
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  call ESMF_GridGet(grid, distgrid = distgrid, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  call ESMF_DistGridGet(distgrid, delayout = delayout, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
 
  call ESMF_DELayoutGet(delayout, localDeCount = decount, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  ! Create cubed sphere grid with coordTypeKind == ESMF_TYPEKIND_R8
  print *, 'Reading file ', mosaicIn
  !call ESMF_MeshWrite(mosaicIn, 'OCNmosaicIn.vtk', rc=localrc)
  grid2 = ESMF_GridCreateMosaic(filename=mosaicIn, &
                staggerLocList= staggerLocList, &
                coordTypeKind = ESMF_TYPEKIND_R8, &
                tileFilePath='./data/', rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
  
  !call ESMF_MeshWrite(grid2, 'OCNmosaic2grid.vtk', rc=localrc)
  
  do s = 1, 2
    do localDe = 0, decount-1  
      call ESMF_GridGetCoord(grid2, coordDim=1, localDe=localDe, &
         staggerloc=staggerLocList(s), farrayPtr=lonPtrR8, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
      call ESMF_GridGetCoord(grid2, coordDim=2, localDe=localDe, &
         staggerloc=staggerLocList(s), farrayPtr=latPtrR8, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
      call ESMF_GridGetCoord(grid, coordDim=1, localDe=localDe, &
       staggerloc=staggerLocList(s), farrayPtr=lonPtrR4, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
      call ESMF_GridGetCoord(grid, coordDim=2, localDe=localDe, &
         staggerloc=staggerLocList(s), farrayPtr=latPtrR4, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

      allocate(lonDiff(size(lonPtrR8,1), size(lonPtrR8,2)))
      allocate(latDiff(size(lonPtrR8,1), size(lonPtrR8,2)))
      total = size(lonPtrR8,1)*size(lonPtrR8,2)
      lonDiff = abs(lonPtrR8-lonPtrR4)
      latDiff = abs(latPtrR8-latPtrR4)

      ! Find the max/min/mean errors
      lonmin = minval(lonDiff)
      latmin = minval(latDiff)
      lonmax = maxval(lonDiff)
      latmax = maxval(latDiff)
      allocate(mean(size(lonDiff,2)))
      do  j=1, size(lonDiff,2)
         mean(j) = sum(lonDiff(:,j))
      enddo
      lonmean = sum(mean)/total
      do j=1, size(latDiff,2)
         mean(j) = sum(latDiff(:,j))
      enddo
      latmean = sum(mean)/total
      
      deallocate(lonDiff, latDiff, mean)

      print *, localPet, localDe, 'min/max/mean:', lonmin, latmin, lonmax, latmax, lonmean, latmean
      if (lonmean > threshhold .or. latmean > threshhold) rc = ESMF_FAILURE
    enddo
  enddo

  ! destroy grid
  call ESMF_GridDestroy(grid,rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
  call ESMF_GridDestroy(grid2,rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  call ESMF_Test((rc==ESMF_SUCCESS), name, failMsg, result, ESMF_SRCLINE)
#endif

  !-----------------------------------------------------------------------------
  !NEX_UTest
  write(name, *) "Test ATM ESMF_GridCreate with different coordTypeKind with GRIDSPEC supergrid tile file"
  write(failMsg, *) "Did not return ESMF_SUCCESS"

  grid=ESMF_GridCreate('data/C48_grid.tile6.nc', &
    ESMF_FILEFORMAT_GRIDSPEC, coordTypeKind=ESMF_TYPEKIND_R4, rc=rc)

#ifdef ESMF_NETCDF

  call ESMF_GridGet(grid, distgrid = distgrid, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  call ESMF_DistGridGet(distgrid, delayout = delayout, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
 
  call ESMF_DELayoutGet(delayout, localDeCount = decount, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  grid2=ESMF_GridCreate('data/C48_grid.tile6.nc', &
    ESMF_FILEFORMAT_GRIDSPEC, coordTypeKind=ESMF_TYPEKIND_R8, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
  
  threshhold = 1.0E-5

  do localDe = 0, decount-1  
      call ESMF_GridGetCoord(grid2, coordDim=1, localDe=localDe, &
         staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=lonPtrR8, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
      call ESMF_GridGetCoord(grid2, coordDim=2, localDe=localDe, &
         staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=latPtrR8, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
      call ESMF_GridGetCoord(grid, coordDim=1, localDe=localDe, &
       staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=lonPtrR4, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
      call ESMF_GridGetCoord(grid, coordDim=2, localDe=localDe, &
         staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=latPtrR4, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

      allocate(lonDiff(size(lonPtrR8,1), size(lonPtrR8,2)))
      allocate(latDiff(size(lonPtrR8,1), size(lonPtrR8,2)))
      total = size(lonPtrR8,1)*size(lonPtrR8,2)
      lonDiff = abs(lonPtrR8-lonPtrR4)
      latDiff = abs(latPtrR8-latPtrR4)

      ! Find the mean errors
      allocate(mean(size(lonDiff,2)))
      do  j=1, size(lonDiff,2)
         mean(j) = sum(lonDiff(:,j))
      enddo
      lonmean = sum(mean)/total
      do j=1, size(latDiff,2)
         mean(j) = sum(latDiff(:,j))
      enddo
      latmean = sum(mean)/total
      
      deallocate(lonDiff, latDiff, mean)

      if (lonmean > threshhold .or. latmean > threshhold) rc = ESMF_FAILURE
  enddo

  call ESMF_GridDestroy(grid,rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
  call ESMF_GridDestroy(grid2,rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  call ESMF_Test((rc==ESMF_SUCCESS), name, failMsg, result, ESMF_SRCLINE)
#else
  write(failMsg, *) "Did not return ESMF_RC_LIB_NOT_PRESENT"
  call ESMF_Test((rc==ESMF_RC_LIB_NOT_PRESENT), name, failMsg, result, ESMF_SRCLINE) 
#endif

  !-----------------------------------------------------------------------------
    !NEX_UTest
  write(name, *) "Test OCN ESMF_GridCreate with different coordTypeKind with GRIDSPEC supergrid tile file"
  write(failMsg, *) "Did not return ESMF_SUCCESS"

  grid=ESMF_GridCreate('data/ocean_hgrid.nc', &
    ESMF_FILEFORMAT_GRIDSPEC, coordTypeKind=ESMF_TYPEKIND_R4, rc=rc)

#ifdef ESMF_NETCDF

  call ESMF_GridGet(grid, distgrid = distgrid, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  call ESMF_DistGridGet(distgrid, delayout = delayout, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
 
  call ESMF_DELayoutGet(delayout, localDeCount = decount, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  grid2=ESMF_GridCreate('data/ocean_hgrid.nc', &
    ESMF_FILEFORMAT_GRIDSPEC, coordTypeKind=ESMF_TYPEKIND_R8, rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
  
  threshhold = 1.0E-5

  do localDe = 0, decount-1  
      call ESMF_GridGetCoord(grid2, coordDim=1, localDe=localDe, &
         staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=lonPtrR8, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
      call ESMF_GridGetCoord(grid2, coordDim=2, localDe=localDe, &
         staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=latPtrR8, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
      call ESMF_GridGetCoord(grid, coordDim=1, localDe=localDe, &
       staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=lonPtrR4, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
      call ESMF_GridGetCoord(grid, coordDim=2, localDe=localDe, &
         staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=latPtrR4, rc=localrc)
      if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

      allocate(lonDiff(size(lonPtrR8,1), size(lonPtrR8,2)))
      allocate(latDiff(size(lonPtrR8,1), size(lonPtrR8,2)))
      total = size(lonPtrR8,1)*size(lonPtrR8,2)
      lonDiff = abs(lonPtrR8-lonPtrR4)
      latDiff = abs(latPtrR8-latPtrR4)

      ! Find the mean errors
      allocate(mean(size(lonDiff,2)))
      do  j=1, size(lonDiff,2)
         mean(j) = sum(lonDiff(:,j))
      enddo
      lonmean = sum(mean)/total
      do j=1, size(latDiff,2)
         mean(j) = sum(latDiff(:,j))
      enddo
      latmean = sum(mean)/total
      
      deallocate(lonDiff, latDiff, mean)

      if (lonmean > threshhold .or. latmean > threshhold) rc = ESMF_FAILURE
  enddo

  call ESMF_GridDestroy(grid,rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE
  call ESMF_GridDestroy(grid2,rc=localrc)
  if (localrc .ne. ESMF_SUCCESS) rc=ESMF_FAILURE

  call ESMF_Test((rc==ESMF_SUCCESS), name, failMsg, result, ESMF_SRCLINE)
#else
  write(failMsg, *) "Did not return ESMF_RC_LIB_NOT_PRESENT"
  call ESMF_Test((rc==ESMF_RC_LIB_NOT_PRESENT), name, failMsg, result, ESMF_SRCLINE) 
#endif

  call ESMF_TestEnd(ESMF_SRCLINE)
  !-----------------------------------------------------------------------------

end program ESMF_AltGridCreateUTest
