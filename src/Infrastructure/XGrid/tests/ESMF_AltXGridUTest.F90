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
      program ESMF_AltXGridUTest

!------------------------------------------------------------------------------

#include "ESMF.h"

!==============================================================================
!BOPI
! !PROGRAM: ESMF_XGridUTest - Unit tests for Field Create and Get methods
!
! !DESCRIPTION:
!
! The code in this file drives F90 Field Create and Get unit tests.
! The companion folder Fieldsrc contains the definitions for the
! Field methods.
!EOPI
!-----------------------------------------------------------------------------
! !USES:
    use ESMF_TestMod     ! test methods
    use ESMF
    implicit none

!------------------------------------------------------------------------------
! The following line turns the CVS identifier string into a printable variable.
    character(*), parameter :: version = &
      '$Id$'

    type(ESMF_XGrid)                    :: xgrid

    ! cumulative result: count failures; no failures equals "all pass"
    integer :: result = 0

    ! individual test result code
    integer :: rc = 1

    ! individual test failure message
    character(ESMF_MAXSTR) :: failMsg
    character(512) :: name

    logical :: isCreated

    call ESMF_TestStart(ESMF_SRCLINE, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
 
    !------------------------------------------------------------------------

    !NEX_UTest
    ! Create an XGrid in 2D from Meshes with user supplied area
    print *, 'Starting test7'
    call test7(rc)
    write(failMsg, *) ""
    write(name, *) "Creating an XGrid in 2D with user supplied area"
    call ESMF_Test((rc.eq.ESMF_SUCCESS), name, failMsg, result, ESMF_SRCLINE)

    !------------------------------------------------------------------------


    call ESMF_TestEnd(ESMF_SRCLINE)
  
contains 
#define ESMF_METHOD "ESMF_TESTS"

  subroutine test7(rc)
    integer, intent(out)                :: rc
    integer                             :: localrc, i, npet
    integer                             :: sideACount, sideBCount, XgridCount  
    character(100)                      :: mosaicA, mosaicB
    type(ESMF_XGrid)                    :: xgrid
    type(ESMF_Grid)                     :: sideA, sideB
    type(ESMF_VM)                       :: vm
    real(ESMF_KIND_R8)                  :: xgrid_area(12), B_area(2,2)
    type(ESMF_Staggerloc)               :: staggerLocList(2)

    rc = ESMF_SUCCESS
    localrc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm=vm, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
      ESMF_ERR_PASSTHRU, &
      ESMF_CONTEXT, rcToReturn=rc)) return

    call ESMF_VMGet(vm, petcount=npet, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
      ESMF_ERR_PASSTHRU, &
      ESMF_CONTEXT, rcToReturn=rc)) return

    ! set up input files
    write(mosaicA, *) "./data/C48_mosaic.nc"
    write(mosaicB, *) "./data/ocean_mosaic.nc"   
    print *, "Read input mosaics and build grids"
    print *, "SideA:",mosaicA
    print *, "SideB:",mosaicB
    
  staggerLocList(1) = ESMF_STAGGERLOC_CENTER
  staggerLocList(2) = ESMF_STAGGERLOC_CORNER
    sideA = ESMF_GridCreateMosaic(filename=mosaicA, &
       staggerLocList= staggerLocList, &
       coordTypeKind = ESMF_TYPEKIND_R8, &
       tileFilePath='./data/', rc=localrc)
    !call ESMF_GridWriteVTK(sideA,staggerloc=ESMF_STAGGERLOC_CORNER, &
       !filename="ATM_grid_in.vtk", &
       !rc=localrc) 
    sideB = ESMF_GridCreateMosaic(filename=mosaicB, &
       staggerLocList= staggerLocList, &
       coordTypeKind = ESMF_TYPEKIND_R8, &
       tileFilePath='./data/', rc=localrc)
    !call ESMF_GridWriteVTK(sideB,staggerloc=ESMF_STAGGERLOC_CORNER, &
       !filename="OCN_grid_in.vtk", &
       !rc=localrc)
    
    
    if (ESMF_LogFoundError(localrc, &
      ESMF_ERR_PASSTHRU, &
      ESMF_CONTEXT, rcToReturn=rc)) return
      
   ! up, down
    print *, "Creating XGrid from Input Mosaics"
    grid = ESMF_XGridCreate(sideAGrid=(/sideA/), &
      sideBGrid=(/sideB/), &
      storeOverlay = .true., &
      rc=localrc)
    if (ESMF_LogFoundError(localrc, &
      ESMF_ERR_PASSTHRU, &
      ESMF_CONTEXT, rcToReturn=rc)) return
      
    ! record the number of cells
      call ESMF_GridGet(grid=sideA, nodeCount=sideACount, rc=localrc)
      call ESMF_GridGet(grid=sideB, nodeCount=sideBCount, rc=localrc)
      call ESMF_XGridGet(xgrid=xgrid, elementCount=XgridCount, rc=localrc)
      print *, "Num cells A/B/X: ",sideACount,"/",sideBCount,"/",XgridCount

    call flux_exchange_sph(xgrid, area_adj=0.95, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
      ESMF_ERR_PASSTHRU, &
      ESMF_CONTEXT, rcToReturn=rc)) return

  end subroutine test7


  !------------------------------------------------------------------------
  ! Utility methods
  !------------------------------------------------------------------------

  function make_grid(atm_nx, atm_ny, atm_dx, atm_dy, atm_sx, atm_sy, field, rc)

    ! return value
    type(ESMF_Grid)                           :: make_grid
    ! arguments
    integer, intent(in)                       :: atm_nx, atm_ny
    real(ESMF_KIND_R4), intent(in)            :: atm_dx, atm_dy
    real(ESMF_KIND_R4), intent(in)            :: atm_sx, atm_sy
    type(ESMF_Field), intent(out), optional   :: field
    integer, intent(out), optional            :: rc

    ! local variables
    integer                                   :: localrc, i, j
    real(ESMF_KIND_R8), pointer               :: coordX(:), coordY(:)
    real(ESMF_KIND_R8)                        :: startx, starty

    make_grid = ESMF_GridCreateNoPeriDim(maxIndex=(/atm_nx, atm_ny/), &
      coordSys=ESMF_COORDSYS_CART, &
      indexflag=ESMF_INDEX_GLOBAL, &
      coordDep1=(/1/), &
      coordDep2=(/2/), &
      rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    call ESMF_GridAddCoord(make_grid, staggerloc=ESMF_STAGGERLOC_CENTER, &
        rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    call ESMF_GridAddCoord(make_grid, staggerloc=ESMF_STAGGERLOC_CORNER, &
        rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    ! global indexing
    ! atm grid is not decomposed in the y direction
    !startx = lpet*atm_nx/npet*atm_dx
    startx = atm_sx
    starty = atm_sy
    ! compute coord
    ! X center
    !call ESMF_GridGetCoord(make_grid, localDE=0, staggerLoc=ESMF_STAGGERLOC_CENTER, &
    !    coordDim=1, farrayPtr=coordX, rc=localrc)
    !if (ESMF_LogFoundError(localrc, &
    !    ESMF_ERR_PASSTHRU, &
    !    ESMF_CONTEXT, rcToReturn=rc)) return
    !do i = lbound(coordX,1), ubound(coordX,1)
    !  coordX(i) = startx + atm_dx/2. + (i-1)*atm_dx
    !enddo
    !print *, 'coordX: ', coordX
    ! X corner
    call ESMF_GridGetCoord(make_grid, localDE=0, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        coordDim=1, farrayPtr=coordX, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    do i = lbound(coordX,1), ubound(coordX,1)
      coordX(i) = startx + (i-1)*atm_dx
    enddo
    !print *, 'startx: ', startx, lbound(coordX, 1), 'coordX: ', coordX
    ! Y center
    !call ESMF_GridGetCoord(make_grid, localDE=0, staggerLoc=ESMF_STAGGERLOC_CENTER, &
    !    coordDim=2, farrayPtr=coordY, rc=localrc)
    !if (ESMF_LogFoundError(localrc, &
    !    ESMF_ERR_PASSTHRU, &
    !    ESMF_CONTEXT, rcToReturn=rc)) return
    !do i = lbound(coordY,1), ubound(coordY,1)
    !  coordY(i) = starty + atm_dy/2. + (i-1)*atm_dy
    !enddo
    ! Y corner
    call ESMF_GridGetCoord(make_grid, localDE=0, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        coordDim=2, farrayPtr=coordY, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    do i = lbound(coordY,1), ubound(coordY,1)
      coordY(i) = starty + (i-1)*atm_dy
    enddo

    if(present(field)) then
      field = ESMF_FieldCreate(make_grid, typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    endif
  
    if(present(rc)) rc = ESMF_SUCCESS

  end function make_grid

  function make_grid_sph(atm_nx, atm_ny, atm_dx, atm_dy, atm_sx, atm_sy, area_adj, tag, scheme, rc)

    ! return value
    type(ESMF_Grid)                           :: make_grid_sph
    ! arguments
    integer, intent(in)                       :: atm_nx, atm_ny
    real(ESMF_KIND_R4), intent(in)            :: atm_dx, atm_dy
    real(ESMF_KIND_R4), intent(in)            :: atm_sx, atm_sy
    real(ESMF_KIND_R4), intent(in), optional  :: area_adj
    character(len=*), intent(in), optional    :: tag
    integer, intent(in) , optional            :: scheme
    integer, intent(out), optional            :: rc

    ! local variables
    integer                                   :: localrc, i, j
    real(ESMF_KIND_R8), pointer               :: coordX(:,:), coordY(:,:)
    real(ESMF_KIND_R8), pointer               :: f_area(:,:), f_area_m(:), o_area(:,:)
    real(ESMF_KIND_R8)                        :: startx, starty
    integer                                   :: l_scheme
    type(ESMF_Mesh)                           :: mesh
    type(ESMF_Field)                          :: field

    l_scheme = ESMF_REGRID_SCHEME_REGION3D
    if(present(scheme)) l_scheme = scheme

    if(l_scheme == ESMF_REGRID_SCHEME_FULL3D) then
      make_grid_sph = ESMF_GridCreate1PeriDim(maxIndex=(/atm_nx, atm_ny/), &
        indexflag=ESMF_INDEX_GLOBAL, &
        gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), &
        !regDecomp=(/npet, 1/), &
        rc=localrc)
    else
      make_grid_sph = ESMF_GridCreateNoPeriDim(maxIndex=(/atm_nx, atm_ny/), &
        indexflag=ESMF_INDEX_GLOBAL, &
        gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/1,1/), &
        !regDecomp=(/npet, 1/), &
        rc=localrc)
    endif 
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    call ESMF_GridAddCoord(make_grid_sph, staggerloc=ESMF_STAGGERLOC_CENTER, &
        rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    call ESMF_GridAddCoord(make_grid_sph, staggerloc=ESMF_STAGGERLOC_CORNER, &
        rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    ! global indexing
    ! atm grid is not decomposed in the y direction
    !startx = lpet*atm_nx/npet*atm_dx
    startx = atm_sx
    starty = atm_sy
    ! compute coord
    ! X center
    call ESMF_GridGetCoord(make_grid_sph, localDE=0, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        coordDim=1, farrayPtr=coordX, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    ! Y center
    call ESMF_GridGetCoord(make_grid_sph, localDE=0, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        coordDim=2, farrayPtr=coordY, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    do i = lbound(coordX,1), ubound(coordX,1)
      do j = lbound(coordX, 2), ubound(coordX, 2)
        coordX(i,j) = startx + atm_dx/2. + (i-1)*atm_dx
        coordY(i,j) = starty + atm_dy/2. + (j-1)*atm_dy
      enddo
    enddo
    !print *, 'startx: ', startx, lbound(coordX, 1), ubound(coordX, 1), 'coordX: ', coordX(:,1)
    ! X corner
    call ESMF_GridGetCoord(make_grid_sph, localDE=0, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        coordDim=1, farrayPtr=coordX, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    ! Y corner
    call ESMF_GridGetCoord(make_grid_sph, localDE=0, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        coordDim=2, farrayPtr=coordY, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    do i = lbound(coordX,1), ubound(coordX,1)
      do j = lbound(coordX, 2), ubound(coordX, 2)
        coordX(i,j) = startx + (i-1)*atm_dx
        coordY(i,j) = starty + (j-1)*atm_dy
      enddo
    enddo

    if(present(area_adj)) then
      ! retrieve area

      !mesh = ESMF_GridToMesh(make_grid_sph, &
      !  ESMF_STAGGERLOC_CORNER, 0, &
      !  regridConserve=ESMF_REGRID_CONSERVE_ON, rc=localrc)
      !if (ESMF_LogFoundError(localrc, &
      !    ESMF_ERR_PASSTHRU, &
      !    ESMF_CONTEXT, rcToReturn=rc)) return

      !allocate(f_area_m(mesh%NumOwnedElements))
      !call ESMF_MeshGetElemArea(mesh,  arealist=f_area_m, rc=localrc)
      !if (ESMF_LogFoundError(localrc, &
      !    ESMF_ERR_PASSTHRU, &
      !    ESMF_CONTEXT, rcToReturn=rc)) return
      !deallocate(f_area_m)

      ! find out original Grid cell area
      field = ESMF_FieldCreate(make_grid_sph, typekind=ESMF_TYPEKIND_R8, &
        staggerloc=ESMF_STAGGERLOC_CENTER, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridGetArea(field, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(field, farrayPtr=o_area, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

      ! add area to Grid
      call ESMF_GridAddItem(make_grid_sph, ESMF_GRIDITEM_AREA, &
        staggerloc=ESMF_STAGGERLOC_CENTER,  rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

      call ESMF_GridGetItem(make_grid_sph, ESMF_GRIDITEM_AREA, &
        staggerloc=ESMF_STAGGERLOC_CENTER, farrayptr=f_area, &
        rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

      ! adjust Grid area
      f_area = area_adj*o_area

    endif

    if(present(rc)) rc = ESMF_SUCCESS

  end function make_grid_sph

  function make_mesh_sph(atm_nx, atm_ny, atm_dx, atm_dy, atm_sx, atm_sy, tag, scheme, rc)

    ! return value
    type(ESMF_Mesh)                           :: make_mesh_sph
    ! arguments
    integer, intent(in)                       :: atm_nx, atm_ny
    real(ESMF_KIND_R4), intent(in)            :: atm_dx, atm_dy
    real(ESMF_KIND_R4), intent(in)            :: atm_sx, atm_sy
    character(len=*), intent(in), optional    :: tag
    integer, intent(in) , optional            :: scheme
    integer, intent(out), optional            :: rc

    ! local variables
    integer                                   :: localrc, i, j
    real(ESMF_KIND_R8), pointer               :: coordX(:,:), coordY(:,:)
    real(ESMF_KIND_R8)                        :: startx, starty
    integer                                   :: l_scheme
    type(ESMF_Grid)                           :: make_grid_sph

    l_scheme = ESMF_REGRID_SCHEME_REGION3D
    if(present(scheme)) l_scheme = scheme

    if(l_scheme == ESMF_REGRID_SCHEME_FULL3D) then
      make_grid_sph = ESMF_GridCreate1PeriDim(maxIndex=(/atm_nx, atm_ny/), &
        indexflag=ESMF_INDEX_GLOBAL, &
        gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,1/), &
        !regDecomp=(/npet, 1/), &
        rc=localrc)
    else
      make_grid_sph = ESMF_GridCreateNoPeriDim(maxIndex=(/atm_nx, atm_ny/), &
        indexflag=ESMF_INDEX_GLOBAL, &
        gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/1,1/), &
        !regDecomp=(/npet, 1/), &
        rc=localrc)
    endif 
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    call ESMF_GridAddCoord(make_grid_sph, staggerloc=ESMF_STAGGERLOC_CENTER, &
        rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    call ESMF_GridAddCoord(make_grid_sph, staggerloc=ESMF_STAGGERLOC_CORNER, &
        rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    ! global indexing
    ! atm grid is not decomposed in the y direction
    !startx = lpet*atm_nx/npet*atm_dx
    startx = atm_sx
    starty = atm_sy
    ! compute coord
    ! X center
    call ESMF_GridGetCoord(make_grid_sph, localDE=0, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        coordDim=1, farrayPtr=coordX, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    ! Y center
    call ESMF_GridGetCoord(make_grid_sph, localDE=0, staggerLoc=ESMF_STAGGERLOC_CENTER, &
        coordDim=2, farrayPtr=coordY, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    do i = lbound(coordX,1), ubound(coordX,1)
      do j = lbound(coordX, 2), ubound(coordX, 2)
        coordX(i,j) = startx + atm_dx/2. + (i-1)*atm_dx
        coordY(i,j) = starty + atm_dy/2. + (j-1)*atm_dy
      enddo
    enddo
    !print *, 'startx: ', startx, lbound(coordX, 1), ubound(coordX, 1), 'coordX: ', coordX(:,1)
    ! X corner
    call ESMF_GridGetCoord(make_grid_sph, localDE=0, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        coordDim=1, farrayPtr=coordX, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    ! Y corner
    call ESMF_GridGetCoord(make_grid_sph, localDE=0, staggerLoc=ESMF_STAGGERLOC_CORNER, &
        coordDim=2, farrayPtr=coordY, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    do i = lbound(coordX,1), ubound(coordX,1)
      do j = lbound(coordX, 2), ubound(coordX, 2)
        coordX(i,j) = startx + (i-1)*atm_dx
        coordY(i,j) = starty + (j-1)*atm_dy
      enddo
    enddo

    make_mesh_sph = ESMF_GridToMesh(make_grid_sph, &
      ESMF_STAGGERLOC_CORNER, 0, &
      regridConserve=ESMF_REGRID_CONSERVE_ON, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    if(present(rc)) rc = ESMF_SUCCESS

  end function make_mesh_sph

  subroutine flux_exchange(xgrid, srcField, dstField, rc)

    type(ESMF_XGrid), intent(inout)           :: xgrid
    type(ESMF_Field), intent(inout)           :: srcField(:)
    type(ESMF_Field), intent(inout)           :: dstField(:)
    integer, intent(out), optional            :: rc


    integer                                   :: localrc, i, j, nsrc, ndst, lpet, npet
    type(ESMF_Field)                          :: f_xgrid
    type(ESMF_Grid), allocatable              :: srcGrid(:)
    type(ESMF_Field), allocatable             :: srcFrac(:), srcArea(:)
    type(ESMF_Grid), allocatable              :: dstGrid(:)
    type(ESMF_Field), allocatable             :: dstFrac(:), dstArea(:)
    type(ESMF_Field), allocatable             :: srcFrac2(:), dstFrac2(:)
    type(ESMF_RouteHandle), allocatable       :: s2d_rh(:,:)
    type(ESMF_RouteHandle), allocatable       :: d2s_rh(:,:)
    type(ESMF_RouteHandle), allocatable       :: s2x_rh(:), x2s_rh(:)
    type(ESMF_RouteHandle), allocatable       :: d2x_rh(:), x2d_rh(:)
    real(ESMF_KIND_R8), pointer               :: src(:,:), dst(:,:), exf(:)
    real(ESMF_KIND_R8), pointer               :: src_area(:,:), dst_area(:,:), exf_area(:)
    real(ESMF_KIND_R8), pointer               :: src_frac(:,:), dst_frac(:,:), exf_frac(:)
    real(ESMF_KIND_R8), pointer               :: src_frac2(:,:), dst_frac2(:,:)
    real(ESMF_KIND_R8)                        :: srcsum(3), allsrcsum(3), scale=2.0, exf_tarea, exf_tflux
    type(ESMF_VM)                             :: vm

    call ESMF_VMGetCurrent(vm=vm, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    call ESMF_VMGet(vm, localpet=lpet, petCount=npet, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    !------------------------------------
    ! build Fields on the Grids
    !------------------------------------

    ! create a Field on the xgrid
    f_xgrid = ESMF_FieldCreate(xgrid=xgrid, TYPEKIND=ESMF_TYPEKIND_R8, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    call ESMF_FieldGet(f_xgrid, farrayPtr=exf, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    nsrc = size(srcField)
    ndst = size(dstField)
    allocate(srcGrid(nsrc), srcFrac(nsrc), srcFrac2(nsrc), srcArea(nsrc))
    allocate(dstGrid(ndst), dstFrac(ndst), dstFrac2(ndst), dstArea(ndst))
    do i = 1, size(srcField)
      call ESMF_FieldGet(srcField(i), grid=srcGrid(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      srcFrac(i) = ESMF_FieldCreate(srcGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      srcFrac2(i) = ESMF_FieldCreate(srcGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      srcArea(i) = ESMF_FieldCreate(srcGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridGetArea(srcArea(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    do i = 1, size(dstField)
      call ESMF_FieldGet(dstField(i), grid=dstGrid(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      dstFrac(i) = ESMF_FieldCreate(dstGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      dstFrac2(i) = ESMF_FieldCreate(dstGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(dstFrac(i), localDe=0, farrayPtr=dst_frac, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      dstArea(i) = ESMF_FieldCreate(dstGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridGetArea(dstArea(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    allocate(s2d_rh(size(srcField), size(dstField)), d2s_rh(size(dstField), size(srcField)))
    allocate(s2x_rh(size(srcField)), x2s_rh(size(srcField)))
    allocate(d2x_rh(size(dstField)), x2d_rh(size(dstField)))

    do i = 1, size(srcField)
      do j = 1, size(dstField)
        call ESMF_FieldRegridStore(srcField=srcField(i), dstField=dstField(j), &
          regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
          routehandle=s2d_rh(i,j), &
          unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
          srcFracField=srcFrac(i), dstFracField=dstFrac(j), & 
          rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return

        call ESMF_FieldRegridStore(srcField=dstField(j), dstField=srcField(i), &
          regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
          routehandle=d2s_rh(j,i), &
          unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
          srcFracField=dstFrac(j), dstFracField=srcFrac(i), & 
          rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
      enddo
    enddo

    do i = 1, size(srcField)
      call ESMF_FieldRegridStore(xgrid, srcField=srcField(i), dstField=f_xgrid, &
        routehandle=s2x_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridStore(xgrid, srcField=f_xgrid, dstField=srcField(i), &
        routehandle=x2s_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo
    do i = 1, size(dstField)
      call ESMF_FieldRegridStore(xgrid, srcField=dstField(i), dstField=f_xgrid, &
        routehandle=d2x_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridStore(xgrid, srcField=f_xgrid, dstField=dstField(i), &
        routehandle=x2d_rh(i), dstFracField=dstFrac(i), dstMergeFracField=dstFrac2(i), &
        rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    !----------------------------------------------------
    ! Compute flux integrals
    ! Initialize src flux to constant
    !----------------------------------------------------
    exf = 0.
    do i = 1, size(srcField)
      call ESMF_FieldGet(srcField(i), farrayPtr=src, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      src = scale
    enddo

    ! Perform flux exchange
    do i = 1, size(srcField)
      call ESMF_FieldRegrid(srcField=srcField(i), dstField=f_xgrid, &
        routehandle=s2x_rh(i), zeroregion=ESMF_REGION_EMPTY, &
        rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    ! make sure flux is conserved on XGrid
    allocate(exf_area(lbound(exf,1):ubound(exf,1)))
    allocate(exf_frac(lbound(exf,1):ubound(exf,1)))
    call ESMF_XGridGet(xgrid, area=exf_area, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    exf_frac = 1.0
    call compute_flux1D(vm, exf, exf_area, exf_frac, allsrcsum, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    if(lpet == 0) print *, ' xgrid flux and area: ', allsrcsum
    if(abs(allsrcsum(1) - allsrcsum(2)*scale) .gt. 1.e-10) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
         msg="- inconsistent flux and area found", &
         ESMF_CONTEXT, rcToReturn=rc) 
      return
    endif
    exf_tflux = allsrcsum(1)
    exf_tarea = allsrcsum(2)
    deallocate(exf_area, exf_frac)

    !make sure flux is conserved on dst Fields
    do i = 1, size(dstField)
      call ESMF_FieldRegrid(srcField=f_xgrid, dstField=dstField(i), &
        routehandle=x2d_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(dstField(i), farrayPtr=dst, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

      ! fraction
      call ESMF_FieldGet(dstFrac(i), farrayPtr=dst_frac, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(dstFrac2(i), farrayPtr=dst_frac2, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

      ! area
      call ESMF_FieldGet(dstArea(i), farrayPtr=dst_area, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

      call compute_flux2D(vm, dst, dst_area, dst_frac, dst_frac2, allsrcsum, dstflux=.true., rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      if(lpet == 0) print *, 'dst flux and area: ', allsrcsum
      if((abs(exf_tarea - allsrcsum(2)) .gt. 1.e-10) .or. &
         (abs(exf_tflux - allsrcsum(1)) .gt. 1.e-10)) then
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
           msg="- inconsistent flux and area found", &
           ESMF_CONTEXT, rcToReturn=rc) 
        return
      endif
    enddo

    do i = 1, size(dstField)
      call ESMF_FieldRegrid(srcField=dstField(i), dstField=f_xgrid, &
        routehandle=d2x_rh(i), &
        rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo
    do i = 1, size(srcField)
      call ESMF_FieldRegrid(srcField=f_xgrid, dstField=srcField(i), &
        routehandle=x2s_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(srcField(i), farrayPtr=src, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    !----------------------------------------------------
    ! clean up
    !----------------------------------------------------
    do i = 1, size(srcField)
      call ESMF_FieldDestroy(srcField(i), rc=localrc)
      call ESMF_FieldDestroy(srcArea(i), rc=localrc)
      call ESMF_FieldDestroy(srcFrac(i), rc=localrc)
      call ESMF_FieldDestroy(srcFrac2(i), rc=localrc)
      call ESMF_RoutehandleRelease(s2x_rh(i), rc=localrc)
      call ESMF_RoutehandleRelease(x2s_rh(i), rc=localrc)
    enddo
    do i = 1, size(dstField)
      call ESMF_FieldDestroy(dstField(i), rc=localrc)
      call ESMF_FieldDestroy(dstArea(i), rc=localrc)
      call ESMF_FieldDestroy(dstFrac(i), rc=localrc)
      call ESMF_FieldDestroy(dstFrac2(i), rc=localrc)
      call ESMF_RoutehandleRelease(d2x_rh(i), rc=localrc)
      call ESMF_RoutehandleRelease(x2d_rh(i), rc=localrc)
    enddo
    do i = 1, size(srcField)
      do j = 1, size(dstField)
        call ESMF_RoutehandleRelease(s2d_rh(i,j), rc=localrc)
        call ESMF_RoutehandleRelease(d2s_rh(j,i), rc=localrc)
      enddo
    enddo

    deallocate(srcArea, srcFrac, dstArea, dstFrac)
    deallocate(s2d_rh, d2s_rh)
    deallocate(s2x_rh, x2s_rh)
    deallocate(d2x_rh, x2d_rh)

    call ESMF_XGridDestroy(xgrid,rc=localrc)

    if(present(rc)) rc = ESMF_SUCCESS

  end subroutine flux_exchange

!----------------------------------------------------
  subroutine compute_flux1D(vm, flux_density, area, fraction, allsum, rc)
    type(ESMF_VM), intent(in)        :: vm
    real(ESMF_KIND_R8), pointer      :: flux_density(:) 
    real(ESMF_KIND_R8), pointer      :: area(:) 
    real(ESMF_KIND_R8), pointer      :: fraction(:) 
    real(ESMF_KIND_R8), intent(out)  :: allsum(3)
    integer, intent(out), optional   :: rc

    real(ESMF_KIND_R8)               :: sum(3)
    integer                          :: i,j, localrc

    if(present(rc)) rc = ESMF_SUCCESS

    sum = 0.
    do i = lbound(flux_density, 1), ubound(flux_density, 1)
      sum(1) = sum(1) + flux_density(i)*area(i)*fraction(i)
      sum(2) = sum(2) +                 area(i)*fraction(i)
      sum(3) = sum(3) +                 area(i)
    enddo

    call ESMF_VMAllReduce(vm, sum, allsum, 3, ESMF_REDUCE_SUM, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

  end subroutine compute_flux1D

  subroutine compute_flux2D(vm, flux_density, area, fraction, fraction2, allsum, dstflux, rc)
    type(ESMF_VM), intent(in)        :: vm
    real(ESMF_KIND_R8), pointer      :: flux_density(:,:) 
    real(ESMF_KIND_R8), pointer      :: area(:,:) 
    real(ESMF_KIND_R8), pointer      :: fraction(:,:) 
    real(ESMF_KIND_R8), pointer      :: fraction2(:,:) 
    real(ESMF_KIND_R8), intent(out)  :: allsum(3)
    logical, intent(in),  optional   :: dstflux
    integer, intent(out), optional   :: rc

    real(ESMF_KIND_R8)               :: sum(3)
    integer                          :: i,j, localrc, npet, lpet
    logical                          :: l_dstflux

    if(present(rc)) rc = ESMF_SUCCESS
    l_dstflux = .false.
    if(present(dstflux)) l_dstflux = dstflux

    call ESMF_VMGet(vm, petCount=npet, localPet=lpet, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
  
    !if(lpet == 0) write(*, '(A, 4I5)') 'compute_flux2D bounds: ', &
    !  lbound(flux_density, 1), ubound(flux_density, 1), &
    !  lbound(flux_density, 2), ubound(flux_density, 2)

    sum = 0.
    do i = lbound(flux_density, 1), ubound(flux_density, 1)
      do j = lbound(flux_density, 2), ubound(flux_density, 2)
        if(l_dstflux) then
          sum(1) = sum(1) + flux_density(i,j)*area(i,j)*fraction2(i,j)
        else
          sum(1) = sum(1) + flux_density(i,j)*area(i,j)*fraction(i,j)*fraction2(i,j)
        endif
        sum(2) = sum(2) +                 area(i,j)*fraction(i,j)
        sum(3) = sum(3) +                 area(i,j)
      enddo
    enddo

    call ESMF_VMAllReduce(vm, sum, allsum, 3, ESMF_REDUCE_SUM, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

  end subroutine compute_flux2D

  subroutine flux_exchange_sph(xgrid, scheme, area_adj, rc)

    type(ESMF_XGrid), intent(inout)           :: xgrid
    integer, intent(in),  optional            :: scheme
    real(ESMF_KIND_R8), pointer               :: coordX(:), coordY(:)
    real(ESMF_KIND_R4), intent(in), optional  :: area_adj
    integer, intent(out), optional            :: rc


    integer                                   :: localrc, i, j, nsrc, ndst, lpet, npet
    type(ESMF_Field)                          :: f_xgrid
    type(ESMF_Grid), allocatable              :: srcGrid(:)
    type(ESMF_Field), allocatable             :: srcFrac(:), srcArea(:)
    type(ESMF_Grid), allocatable              :: dstGrid(:)
    type(ESMF_Field), allocatable             :: dstFrac(:), dstArea(:)
    type(ESMF_Field), allocatable             :: srcFrac2(:)
    type(ESMF_Field), allocatable             :: dstFrac2(:)
    type(ESMF_RouteHandle), allocatable       :: s2d_rh(:,:)
    type(ESMF_RouteHandle), allocatable       :: d2s_rh(:,:)
    type(ESMF_RouteHandle), allocatable       :: s2x_rh(:), x2s_rh(:)
    type(ESMF_RouteHandle), allocatable       :: d2x_rh(:), x2d_rh(:)
    real(ESMF_KIND_R8), pointer               :: src(:,:), dst(:,:), exf(:)
    real(ESMF_KIND_R8), pointer               :: src_area(:,:), dst_area(:,:), exf_area(:)
    real(ESMF_KIND_R8), pointer               :: src_frac(:,:), dst_frac(:,:), exf_frac(:)
    real(ESMF_KIND_R8), pointer               :: src_frac2(:,:), dst_frac2(:,:)
    real(ESMF_KIND_R8)                        :: srcsum(3), allsrcsum(3), scale=2.0, exf_tarea, exf_tflux
    type(ESMF_VM)                             :: vm
    type(ESMF_Field), allocatable             :: srcField(:)
    type(ESMF_Field), allocatable             :: dstField(:)
    integer                                   :: l_scheme
    integer                                   :: sideAGC, sideBGC, sideAMC, sideBMC
    real(ESMF_KIND_R8)                        :: global_sum, l_area_adj
    character(len=1)                          :: cids(10) = (/'0','1','2','3','4','5','6','7','8','9'/)

    l_scheme = ESMF_REGRID_SCHEME_REGION3D
    if(present(scheme)) l_scheme = scheme
    l_area_adj = 1.0
    if(present(area_adj)) l_area_adj = area_adj

    call ESMF_VMGetCurrent(vm=vm, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    call ESMF_VMGet(vm, localpet=lpet, petCount=npet, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    !------------------------------------
    ! build Fields on the Grids
    !------------------------------------

    ! create a Field on the xgrid
    f_xgrid = ESMF_FieldCreate(xgrid=xgrid, TYPEKIND=ESMF_TYPEKIND_R8, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    call ESMF_FieldGet(f_xgrid, farrayPtr=exf, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    call ESMF_XGridGet(xgrid, &
        sideAGridCount=sideAGC, sideBGridCount=sideBGC, &
        sideAMeshCount=sideAMC, sideBMeshCount=sideBMC, &
        rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    nsrc = sideAGC
    ndst = sideBGC
    allocate(srcGrid(nsrc), srcField(nsrc), srcFrac(nsrc), srcFrac2(nsrc), srcArea(nsrc))
    allocate(dstGrid(ndst), dstField(ndst), dstFrac(ndst), dstFrac2(ndst), dstArea(ndst))

    call ESMF_XGridGet(xgrid, sideAGrid=srcGrid, sideBGrid=dstGrid, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    do i = 1, nsrc
      srcField(i) = ESMF_FieldCreate(srcGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      srcFrac(i) = ESMF_FieldCreate(srcGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      srcFrac2(i) = ESMF_FieldCreate(srcGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      srcArea(i) = ESMF_FieldCreate(srcGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridGetArea(srcArea(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    do i = 1, ndst
      dstField(i) = ESMF_FieldCreate(dstGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      dstFrac(i) = ESMF_FieldCreate(dstGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      dstFrac2(i) = ESMF_FieldCreate(dstGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(dstFrac(i), localDe=0, farrayPtr=dst_frac, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      dstArea(i) = ESMF_FieldCreate(dstGrid(i), typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridGetArea(dstArea(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    allocate(s2d_rh(size(srcField), size(dstField)), d2s_rh(size(dstField), size(srcField)))
    allocate(s2x_rh(size(srcField)), x2s_rh(size(srcField)))
    allocate(d2x_rh(size(dstField)), x2d_rh(size(dstField)))

    do i = 1, size(srcField)
      do j = 1, size(dstField)
        call ESMF_FieldRegridStore(srcField=srcField(i), dstField=dstField(j), &
          regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
          routehandle=s2d_rh(i,j), &
          unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
          srcFracField=srcFrac(i), dstFracField=dstFrac(j), & 
          rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return

        call ESMF_FieldRegridStore(srcField=dstField(j), dstField=srcField(i), &
          regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
          routehandle=d2s_rh(j,i), &
          unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
          srcFracField=dstFrac(j), dstFracField=srcFrac(i), & 
          rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
      enddo
    enddo

    do i = 1, size(srcField)
      call ESMF_FieldRegridStore(xgrid, srcField=srcField(i), dstField=f_xgrid, &
        routehandle=s2x_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridStore(xgrid, srcField=f_xgrid, dstField=srcField(i), &
        routehandle=x2s_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      !call ESMF_GridWrite(srcGrid(i), cids(i)//'_srcmesh.vtk', rc=localrc)
      !if (ESMF_LogFoundError(localrc, &
      !    ESMF_ERR_PASSTHRU, &
      !    ESMF_CONTEXT, rcToReturn=rc)) return
    enddo
    do i = 1, size(dstField)
      call ESMF_FieldRegridStore(xgrid, srcField=dstField(i), dstField=f_xgrid, &
        routehandle=d2x_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridStore(xgrid, srcField=f_xgrid, dstField=dstField(i), &
        routehandle=x2d_rh(i), dstFracField=dstFrac(i), dstMergeFracField=dstFrac2(i), &
        rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      !call ESMF_GridWrite(dstGrid(i), cids(i)//'_dstmesh.vtk', rc=localrc)
      !if (ESMF_LogFoundError(localrc, &
      !    ESMF_ERR_PASSTHRU, &
      !    ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    !----------------------------------------------------
    ! Compute flux integrals
    ! Initialize src flux to constant
    !----------------------------------------------------
    exf = 0.
    do i = 1, size(srcField)
      call ESMF_FieldGet(srcField(i), farrayPtr=src, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      src = scale
    enddo

    ! Perform flux exchange
    do i = 1, size(srcField)
      call ESMF_FieldRegrid(srcField=srcField(i), dstField=f_xgrid, &
        routehandle=s2x_rh(i), zeroregion=ESMF_REGION_EMPTY, &
        rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    ! make sure flux is conserved on XGrid
    !call ESMF_MeshWrite(xgrid%xgtypep%mesh, 'xgrid.vtk', rc=localrc)
    !if (ESMF_LogFoundError(localrc, &
    !    ESMF_ERR_PASSTHRU, &
    !    ESMF_CONTEXT, rcToReturn=rc)) return
    allocate(exf_area(lbound(exf,1):ubound(exf,1)))
    allocate(exf_frac(lbound(exf,1):ubound(exf,1)))
    call ESMF_XGridGet(xgrid, area=exf_area, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    exf_frac = 1.0
    call compute_flux1D(vm, exf, exf_area, exf_frac, allsrcsum, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    if(lpet == 0) print *, ' xgrid flux and area: ', allsrcsum
    if(abs(allsrcsum(1) - allsrcsum(2)*scale*l_area_adj) .gt. 1.e-10) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
         msg="- inconsistent flux and area found", &
         ESMF_CONTEXT, rcToReturn=rc) 
      return
    endif
    exf_tflux = allsrcsum(1)
    exf_tarea = allsrcsum(2)
    deallocate(exf_area, exf_frac)

    !make sure flux is conserved on dst Fields
    global_sum = 0.
    do i = 1, size(dstField)
      call ESMF_FieldRegrid(srcField=f_xgrid, dstField=dstField(i), &
        routehandle=x2d_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(dstField(i), farrayPtr=dst, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

      ! fraction
      call ESMF_FieldGet(dstFrac(i), farrayPtr=dst_frac, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(dstFrac2(i), farrayPtr=dst_frac2, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

      ! area
      call ESMF_FieldGet(dstArea(i), farrayPtr=dst_area, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

      call compute_flux2D(vm, dst, dst_area, dst_frac, dst_frac2, allsrcsum, dstflux=.true., rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      if(lpet == 0) print *, 'dst flux and area: ', allsrcsum
      if(ndst == 1) then
        if((abs(exf_tarea*l_area_adj - allsrcsum(2)) .gt. 1.e-10) .or. &
           (abs(exf_tflux - allsrcsum(1)) .gt. 1.e-10)) then
          call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
             msg="- inconsistent flux and area found", &
             ESMF_CONTEXT, rcToReturn=rc) 
          return
        endif
      else
        call compute_flux2D(vm, dst, dst_area, dst_frac, dst_frac2, allsrcsum, dstflux=.true., rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
        if(lpet == 0) print *, 'dst flux and area using frac2: ', allsrcsum
        global_sum = global_sum + allsrcsum(1)
      endif

    enddo

    ! make sure going to multiple Grids also conserve global flux
    if(ndst .gt. 1) then
        if ((abs(exf_tflux - global_sum) .gt. 1.e-10)) then
        call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
           msg="- inconsistent flux and area found", &
           ESMF_CONTEXT, rcToReturn=rc) 
        return
      endif
    endif

    do i = 1, size(dstField)
      call ESMF_FieldRegrid(srcField=dstField(i), dstField=f_xgrid, &
        routehandle=d2x_rh(i), &
        rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo
    do i = 1, size(srcField)
      call ESMF_FieldRegrid(srcField=f_xgrid, dstField=srcField(i), &
        routehandle=x2s_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(srcField(i), farrayPtr=src, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    !----------------------------------------------------
    ! clean up
    !----------------------------------------------------
    do i = 1, size(srcField)
      call ESMF_FieldDestroy(srcField(i), rc=localrc)
      call ESMF_FieldDestroy(srcArea(i), rc=localrc)
      call ESMF_FieldDestroy(srcFrac(i), rc=localrc)
      call ESMF_FieldDestroy(srcFrac2(i), rc=localrc)
      call ESMF_RoutehandleRelease(s2x_rh(i), rc=localrc)
      call ESMF_RoutehandleRelease(x2s_rh(i), rc=localrc)
    enddo
    do i = 1, size(dstField)
      call ESMF_FieldDestroy(dstField(i), rc=localrc)
      call ESMF_FieldDestroy(dstArea(i), rc=localrc)
      call ESMF_FieldDestroy(dstFrac(i), rc=localrc)
      call ESMF_FieldDestroy(dstFrac2(i), rc=localrc)
      call ESMF_RoutehandleRelease(d2x_rh(i), rc=localrc)
      call ESMF_RoutehandleRelease(x2d_rh(i), rc=localrc)
    enddo
    do i = 1, size(srcField)
      do j = 1, size(dstField)
        call ESMF_RoutehandleRelease(s2d_rh(i,j), rc=localrc)
        call ESMF_RoutehandleRelease(d2s_rh(j,i), rc=localrc)
      enddo
    enddo

    deallocate(srcArea, srcFrac, dstArea, dstFrac)
    deallocate(s2d_rh, d2s_rh)
    deallocate(s2x_rh, x2s_rh)
    deallocate(d2x_rh, x2d_rh)

    call ESMF_XGridDestroy(xgrid,rc=localrc)

    if(present(rc)) rc = ESMF_SUCCESS

  end subroutine flux_exchange_sph

  subroutine flux_exchange_sph_mesh(xgrid, scheme, rc)

    type(ESMF_XGrid), intent(inout)           :: xgrid
    integer, intent(in),  optional            :: scheme
    integer, intent(out), optional            :: rc


    integer                                   :: localrc, i, j, nsrc, ndst, lpet, npet
    type(ESMF_Field)                          :: f_xgrid
    type(ESMF_Mesh), allocatable              :: srcGrid(:)
    type(ESMF_Field), allocatable             :: srcFrac(:), srcArea(:)
    type(ESMF_Mesh), allocatable              :: dstGrid(:)
    type(ESMF_Field), allocatable             :: dstFrac(:), dstArea(:)
    type(ESMF_Field), allocatable             :: srcFrac2(:)
    type(ESMF_Field), allocatable             :: dstFrac2(:)
    type(ESMF_RouteHandle), allocatable       :: s2d_rh(:,:)
    type(ESMF_RouteHandle), allocatable       :: d2s_rh(:,:)
    type(ESMF_RouteHandle), allocatable       :: s2x_rh(:), x2s_rh(:)
    type(ESMF_RouteHandle), allocatable       :: d2x_rh(:), x2d_rh(:)
    real(ESMF_KIND_R8), pointer               :: src(:), dst(:), exf(:)
    real(ESMF_KIND_R8), pointer               :: src_area(:), dst_area(:), exf_area(:)
    real(ESMF_KIND_R8), pointer               :: src_frac(:), dst_frac(:), exf_frac(:)
    real(ESMF_KIND_R8), pointer               :: src_frac2(:), dst_frac2(:)
    real(ESMF_KIND_R8)                        :: srcsum(3), allsrcsum(3), scale=2.0, exf_tarea, exf_tflux
    type(ESMF_VM)                             :: vm
    type(ESMF_Field), allocatable             :: srcField(:)
    type(ESMF_Field), allocatable             :: dstField(:)
    integer                                   :: l_scheme
    integer                                   :: sideAGC, sideBGC, sideAMC, sideBMC
    real(ESMF_KIND_R8)                        :: global_sum
    character(len=1)                          :: cids(10) = (/'0','1','2','3','4','5','6','7','8','9'/)

    l_scheme = ESMF_REGRID_SCHEME_REGION3D
    if(present(scheme)) l_scheme = scheme

    call ESMF_VMGetCurrent(vm=vm, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    call ESMF_VMGet(vm, localpet=lpet, petCount=npet, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    !------------------------------------
    ! build Fields on the Grids
    !------------------------------------

    ! create a Field on the xgrid
    f_xgrid = ESMF_FieldCreate(xgrid=xgrid, TYPEKIND=ESMF_TYPEKIND_R8, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    call ESMF_FieldGet(f_xgrid, farrayPtr=exf, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    call ESMF_XGridGet(xgrid, &
        sideAGridCount=sideAGC, sideBGridCount=sideBGC, &
        sideAMeshCount=sideAMC, sideBMeshCount=sideBMC, &
        rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    nsrc = sideAMC
    ndst = sideBMC
    allocate(srcGrid(nsrc), srcField(nsrc), srcFrac(nsrc), srcFrac2(nsrc), srcArea(nsrc))
    allocate(dstGrid(ndst), dstField(ndst), dstFrac(ndst), dstFrac2(ndst), dstArea(ndst))

    call ESMF_XGridGet(xgrid, sideAMesh=srcGrid, sideBMesh=dstGrid, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return

    do i = 1, nsrc
      srcField(i) = ESMF_FieldCreate(srcGrid(i), meshloc=ESMF_MESHLOC_ELEMENT, &
        typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      srcFrac(i) = ESMF_FieldCreate(srcGrid(i), meshloc=ESMF_MESHLOC_ELEMENT, &
        typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      srcFrac2(i) = ESMF_FieldCreate(srcGrid(i), meshloc=ESMF_MESHLOC_ELEMENT, &
        typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      srcArea(i) = ESMF_FieldCreate(srcGrid(i), meshloc=ESMF_MESHLOC_ELEMENT, &
        typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridGetArea(srcArea(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    do i = 1, ndst
      dstField(i) = ESMF_FieldCreate(dstGrid(i), meshloc=ESMF_MESHLOC_ELEMENT, &
        typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      dstFrac(i) = ESMF_FieldCreate(dstGrid(i), meshloc=ESMF_MESHLOC_ELEMENT, &
        typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      dstFrac2(i) = ESMF_FieldCreate(dstGrid(i), meshloc=ESMF_MESHLOC_ELEMENT, &
        typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(dstFrac(i), localDe=0, farrayPtr=dst_frac, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      dstArea(i) = ESMF_FieldCreate(dstGrid(i), meshloc=ESMF_MESHLOC_ELEMENT, &
        typekind=ESMF_TYPEKIND_R8, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridGetArea(dstArea(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    allocate(s2d_rh(size(srcField), size(dstField)), d2s_rh(size(dstField), size(srcField)))
    allocate(s2x_rh(size(srcField)), x2s_rh(size(srcField)))
    allocate(d2x_rh(size(dstField)), x2d_rh(size(dstField)))

    do i = 1, size(srcField)
      do j = 1, size(dstField)
        call ESMF_FieldRegridStore(srcField=srcField(i), dstField=dstField(j), &
          regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
          routehandle=s2d_rh(i,j), &
          unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
          srcFracField=srcFrac(i), dstFracField=dstFrac(j), & 
          rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return

        call ESMF_FieldRegridStore(srcField=dstField(j), dstField=srcField(i), &
          regridmethod=ESMF_REGRIDMETHOD_CONSERVE, &
          routehandle=d2s_rh(j,i), &
          unmappedaction = ESMF_UNMAPPEDACTION_IGNORE, &
          srcFracField=dstFrac(j), dstFracField=srcFrac(i), & 
          rc=localrc)
        if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
      enddo
    enddo

    do i = 1, size(srcField)
      call ESMF_FieldRegridStore(xgrid, srcField=srcField(i), dstField=f_xgrid, &
        routehandle=s2x_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridStore(xgrid, srcField=f_xgrid, dstField=srcField(i), &
        routehandle=x2s_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      !call ESMF_MeshWrite(srcGrid(i), cids(i)//'_srcmesh.vtk', rc=localrc)
      !if (ESMF_LogFoundError(localrc, &
      !    ESMF_ERR_PASSTHRU, &
      !    ESMF_CONTEXT, rcToReturn=rc)) return
    enddo
    do i = 1, size(dstField)
      call ESMF_FieldRegridStore(xgrid, srcField=dstField(i), dstField=f_xgrid, &
        routehandle=d2x_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldRegridStore(xgrid, srcField=f_xgrid, dstField=dstField(i), &
        routehandle=x2d_rh(i), dstFracField=dstFrac(i), dstMergeFracField=dstFrac2(i), &
        rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      !call ESMF_MeshWrite(dstGrid(i), cids(i)//'_dstmesh.vtk', rc=localrc)
      !if (ESMF_LogFoundError(localrc, &
      !    ESMF_ERR_PASSTHRU, &
      !    ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    !----------------------------------------------------
    ! Compute flux integrals
    ! Initialize src flux to constant
    !----------------------------------------------------
    exf = 0.
    do i = 1, size(srcField)
      call ESMF_FieldGet(srcField(i), farrayPtr=src, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      src = scale
    enddo

    ! Perform flux exchange
    do i = 1, size(srcField)
      call ESMF_FieldRegrid(srcField=srcField(i), dstField=f_xgrid, &
        routehandle=s2x_rh(i), zeroregion=ESMF_REGION_EMPTY, &
        rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    ! make sure flux is conserved on XGrid
    !call ESMF_MeshWrite(xgrid%xgtypep%mesh, 'xgrid.vtk', rc=localrc)
    !if (ESMF_LogFoundError(localrc, &
    !    ESMF_ERR_PASSTHRU, &
    !    ESMF_CONTEXT, rcToReturn=rc)) return
    allocate(exf_area(lbound(exf,1):ubound(exf,1)))
    allocate(exf_frac(lbound(exf,1):ubound(exf,1)))
    call ESMF_XGridGet(xgrid, area=exf_area, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    exf_frac = 1.0
    call compute_flux1D(vm, exf, exf_area, exf_frac, allsrcsum, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
        ESMF_ERR_PASSTHRU, &
        ESMF_CONTEXT, rcToReturn=rc)) return
    if(lpet == 0) print *, ' xgrid flux and area: ', allsrcsum
    !if(abs(allsrcsum(1) - allsrcsum(2)*scale) .gt. 1.e-10) then
    !  call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
    !     msg="- inconsistent flux and area found", &
    !     ESMF_CONTEXT, rcToReturn=rc) 
    !  return
    !endif
    exf_tflux = allsrcsum(1)
    exf_tarea = allsrcsum(2)
    deallocate(exf_area, exf_frac)

    !make sure flux is conserved on dst Fields
    global_sum = 0.
    do i = 1, size(dstField)
      call ESMF_FieldRegrid(srcField=f_xgrid, dstField=dstField(i), &
        routehandle=x2d_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(dstField(i), farrayPtr=dst, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

      ! fraction
      call ESMF_FieldGet(dstFrac(i), farrayPtr=dst_frac, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(dstFrac2(i), farrayPtr=dst_frac2, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

      ! area
      call ESMF_FieldGet(dstArea(i), farrayPtr=dst_area, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return

      !call compute_flux2D(vm, dst, dst_area, dst_frac, dst_frac2, allsrcsum, dstflux=.true., rc=localrc)
      !if (ESMF_LogFoundError(localrc, &
      !    ESMF_ERR_PASSTHRU, &
      !    ESMF_CONTEXT, rcToReturn=rc)) return
      !if(lpet == 0) print *, 'dst flux and area: ', allsrcsum
      !if(ndst == 1) then
      !  if((abs(exf_tarea - allsrcsum(2)) .gt. 1.e-10) .or. &
      !     (abs(exf_tflux - allsrcsum(1)) .gt. 1.e-10)) then
      !    call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
      !       msg="- inconsistent flux and area found", &
      !       ESMF_CONTEXT, rcToReturn=rc) 
      !    return
      !  endif
      !else
      !  call compute_flux2D(vm, dst, dst_area, dst_frac, dst_frac2, allsrcsum, dstflux=.true., rc=localrc)
      !  if (ESMF_LogFoundError(localrc, &
      !      ESMF_ERR_PASSTHRU, &
      !      ESMF_CONTEXT, rcToReturn=rc)) return
      !  if(lpet == 0) print *, 'dst flux and area using frac2: ', allsrcsum
      !  global_sum = global_sum + allsrcsum(1)
      !endif

    enddo

    ! make sure going to multiple Grids also conserve global flux
    !if(ndst .gt. 1) then
    !    if ((abs(exf_tflux - global_sum) .gt. 1.e-10)) then
    !    call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, & 
    !       msg="- inconsistent flux and area found", &
    !       ESMF_CONTEXT, rcToReturn=rc) 
    !    return
    !  endif
    !endif

    do i = 1, size(dstField)
      call ESMF_FieldRegrid(srcField=dstField(i), dstField=f_xgrid, &
        routehandle=d2x_rh(i), &
        rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo
    do i = 1, size(srcField)
      call ESMF_FieldRegrid(srcField=f_xgrid, dstField=srcField(i), &
        routehandle=x2s_rh(i), rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
      call ESMF_FieldGet(srcField(i), farrayPtr=src, rc=localrc)
      if (ESMF_LogFoundError(localrc, &
          ESMF_ERR_PASSTHRU, &
          ESMF_CONTEXT, rcToReturn=rc)) return
    enddo

    !----------------------------------------------------
    ! clean up
    !----------------------------------------------------
    do i = 1, size(srcField)
      call ESMF_FieldDestroy(srcField(i), rc=localrc)
      call ESMF_FieldDestroy(srcArea(i), rc=localrc)
      call ESMF_FieldDestroy(srcFrac(i), rc=localrc)
      call ESMF_FieldDestroy(srcFrac2(i), rc=localrc)
      call ESMF_RoutehandleRelease(s2x_rh(i), rc=localrc)
      call ESMF_RoutehandleRelease(x2s_rh(i), rc=localrc)
    enddo
    do i = 1, size(dstField)
      call ESMF_FieldDestroy(dstField(i), rc=localrc)
      call ESMF_FieldDestroy(dstArea(i), rc=localrc)
      call ESMF_FieldDestroy(dstFrac(i), rc=localrc)
      call ESMF_FieldDestroy(dstFrac2(i), rc=localrc)
      call ESMF_RoutehandleRelease(d2x_rh(i), rc=localrc)
      call ESMF_RoutehandleRelease(x2d_rh(i), rc=localrc)
    enddo
    do i = 1, size(srcField)
      do j = 1, size(dstField)
        call ESMF_RoutehandleRelease(s2d_rh(i,j), rc=localrc)
        call ESMF_RoutehandleRelease(d2s_rh(j,i), rc=localrc)
      enddo
    enddo

    deallocate(srcArea, srcFrac, dstArea, dstFrac)
    deallocate(s2d_rh, d2s_rh)
    deallocate(s2x_rh, x2s_rh)
    deallocate(d2x_rh, x2d_rh)

    call ESMF_XGridDestroy(xgrid,rc=localrc)

    if(present(rc)) rc = ESMF_SUCCESS

  end subroutine flux_exchange_sph_mesh

end program ESMF_AltXGridUTest
