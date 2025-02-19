// $Id$
//
// Earth System Modeling Framework
// Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
// Massachusetts Institute of Technology, Geophysical Fluid Dynamics
// Laboratory, University of Michigan, National Centers for Environmental
// Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
// NASA Goddard Space Flight Center.
// Licensed under the University of Illinois-NCSA License.
//
//==============================================================================
#define ESMC_FILENAME "ESMCI_Regrid_F.C"
//==============================================================================
//
// This file contains the Fortran interface code to link F90 and C++.
//
//------------------------------------------------------------------------------
// INCLUDES
//------------------------------------------------------------------------------
#include "ESMCI_Macros.h"
#include "ESMCI_VM.h"
#include "ESMCI_LogErr.h"
#include "ESMCI_Grid.h"
#include "ESMCI_GridToMesh.h"
#include "ESMC_Util.h"
#include "ESMCI_Array.h"
#include "ESMCI_PointList.h"
#include "Mesh/include/ESMCI_Mesh.h"
#include "Mesh/include/ESMCI_MeshCap.h"
#include "Mesh/include/Regridding/ESMCI_Integrate.h"
#include "Mesh/include/Regridding/ESMCI_ExtrapolationPoleLGC.h"
#include "Mesh/include/Legacy/ESMCI_MeshRead.h"
#include "Mesh/include/Legacy/ESMCI_Exception.h"

//------------------------------------------------------------------------------
//BOP
// !DESCRIPTION:
//
//
//EOP
//-------------------------------------------------------------------------


using namespace ESMCI;


extern "C" void FTN_X(c_esmc_regrid_create)(MeshCap **meshsrcpp,
                                            ESMCI::Array **arraysrcpp,
                                            ESMCI::PointList **plsrcpp,
                                            int *src_pl_used,
                                            MeshCap **meshdstpp,
                                            ESMCI::Array **arraydstpp,
                                            ESMCI::PointList **pldstpp,
                                            int *dst_pl_used,
                                            int *regridMethod,
                                            int *map_type,
                                            int *norm_type,
                                            int *_vectorRegrid, 
                                            int *regridPoleType, int *regridPoleNPnts,
                                            int *extrapMethod,
                                            int *extrapNumSrcPnts,
                                            ESMC_R8 *extrapDistExponent,
                                            int *extrapNumLevels,
                                            int *extrapNumInputLevels,
                                            int *unmappedaction, int *_ignoreDegenerate,
                                            int *srcTermProcessing, int *pipelineDepth,
                                            ESMCI::RouteHandle **rh, int *has_rh, int *has_iw,
                                            int *nentries, ESMCI::TempWeights **tweights,
                                            int *has_udl, int *_num_udl, ESMCI::TempUDL **_tudl,
                                            int *has_statusArray, ESMCI::Array **statusArray,
                                            int *checkFlag, 
                                            int*rc) {
#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_regrid_create()"

  // Nullify Mesh or PointList based on usage
  if (*src_pl_used==1) {
    *meshsrcpp=NULL;
  } else {
    *plsrcpp=NULL;
  }

  if (*dst_pl_used==1) {
    // TODO: figure out how to include ESMC_EXTRAPMETHOD_CREEP here without MOAB
    if ((*extrapMethod != ESMC_EXTRAPMETHOD_CREEP) && (*extrapMethod != ESMC_EXTRAPMETHOD_CREEP_NRST_D)) {
      *meshdstpp=NULL;    
    }
  } else {
    *pldstpp=NULL;
  }

MeshCap::regrid_create(meshsrcpp, arraysrcpp, plsrcpp,
                       meshdstpp, arraydstpp, pldstpp,
                       regridMethod,
                       map_type,
                       norm_type,
                       _vectorRegrid,                        
                       regridPoleType, regridPoleNPnts,
                       extrapMethod,
                       extrapNumSrcPnts,
                       extrapDistExponent,
                       extrapNumLevels, 
                       extrapNumInputLevels, 
                       unmappedaction, _ignoreDegenerate,
                       srcTermProcessing, pipelineDepth,
                       rh, has_rh, has_iw,
                       nentries, tweights,
                       has_udl, _num_udl, _tudl,
                       has_statusArray, statusArray,
                       checkFlag, 
                       rc);
}

extern "C" void FTN_X(c_esmc_regrid_getiwts)(Grid **gridpp,
                   MeshCap **meshpp, ESMCI::Array **arraypp, int *staggerLoc,
                   int *rc) {
#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_regrid_getiwts()"
  MeshCap::regrid_getiwts(gridpp,
                          meshpp, arraypp, staggerLoc,
                          rc);
}


extern "C" void FTN_X(c_esmc_regrid_getarea)(Grid **gridpp,
                   MeshCap **meshpp, ESMCI::Array **arraypp, int *staggerLoc,
                   int *rc) {
#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_regrid_getarea()"
  MeshCap::regrid_getarea(gridpp,
                          meshpp, arraypp, staggerLoc,
                          rc);
}


extern "C" void FTN_X(c_esmc_regrid_getfrac)(Grid **gridpp,
                   MeshCap **meshpp, ESMCI::Array **arraypp, int *staggerLoc,
                   int *rc) {
#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_regrid_getfrac()"
  MeshCap::regrid_getfrac(gridpp,
                          meshpp, arraypp, staggerLoc,
                          rc);
}


#undef  ESMC_METHOD
