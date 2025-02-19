// $Id$
//
// Earth System Modeling Framework
// Copyright (c) 2002-2024, University Corporation for Atmospheric Research,    
// Massachusetts Institute of Technology, Geophysical Fluid Dynamics
// Laboratory, University of Michigan, National Centers for Environmental
// Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
// NASA Goddard Space Flight Center.
// Licensed under the University of Illinois-NCSA License.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// ESMCI Test method implementation (body) file
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//  
// !DESCRIPTION:
//  
// The code in this file implements the C++ Test methods declared
// in the companion file ESMCI_Test.h
//
//-----------------------------------------------------------------------------
//
// insert any higher level, 3rd party or system includes here
#include <stdio.h>
#if !defined (ESMF_OS_MinGW)
#include <sys/time.h>
#else
#include <windows.h>
#endif
#include "ESMCI.h"
#include "ESMC.h"

// associated class definition file
#include "ESMCI_Test.h"

timeval  start_time;
int PETnum;
//-----------------------------------------------------------------------------
// leave the following line as-is; it will insert the cvs ident string
// into the object file for tracking purposes.
static const char *const version = "$Id$";
//-----------------------------------------------------------------------------

namespace ESMCI {

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//
// This section includes all the Test routines
//
// TODO: these all use printf, but should use log write routines so output
//  goes to the log, not stdout.
//

//-----------------------------------------------------------------------------
//BOP
// !IROUTINE:  ESMCI::Test() - Prints whether a test passed or failed
//
// !INTERFACE:   
int Test(
//
// !RETURN VALUE:
//    ESMF_SUCCESS or ESMF_FAILURE
//
// !ARGUMENTS:
  int condition,       // in - the test pass/fail condition
  const char *name,    // in - the test name
  const char *failMsg, // in - optional message printed on test failure
  int *result,         // in/out - cumulative failure count
  const char *file,    // in - test filename
  int line,            // in - test line number in test filename
  int only) {          // in - if set to 0, print on stderr also
//
// !DESCRIPTION:
//    Prints PASS/FAIL based on passed-in condition.  If FAIL, prints
//    optional failure message and increments failure result counter.
//    If {\tt only} is zero, also print same message to stderr as well
//    as the normal output on stdout.  The default for {\tt only} is 1.
//
//EOP
//-----------------------------------------------------------------------------
  std::stringstream msgbuf;
  ESMCI::LogErr *whichLog;

  // TODO: this should be settable by the user
  whichLog = &ESMC_LogDefault;

  if (name == NULL || result == NULL || failMsg == NULL || file == NULL) {
    msgbuf << "FAIL " << __FILE__ << ", line " << __LINE__
           << ", null pointer(s) passed to ESMCI::Test()";
    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());
    return(ESMF_FAILURE);
  }

  if (condition) {
    msgbuf << "PASS " << name << ", " << file << ", line " << line;
    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());
  }else {
    msgbuf << "FAIL " << name << ", " << file << ", line " << line << ": "
      << failMsg;
    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());
    (*result)++; // count total failures; 0 = all pass
  }

  return(ESMF_SUCCESS);

}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//BOP
// !IROUTINE:  ESMCI::TestEnd() - Print a standard end test message
//
// !INTERFACE:   
int TestEnd(
//
// !RETURN VALUE:
//    ESMF_SUCCESS or ESMF_FAILURE
//
// !ARGUMENTS:
  const char *file, // in - test filename
  int line,         // in - test line number in test filename
  int only) {       // in - if set to 0, print on stderr also
// 
// !DESCRIPTION:
//    Prints summary message about total failures, and standard exit message.
//    If {\tt only} is zero, also print same message to stderr as well
//    as the normal output on stdout.  The default for {\tt only} is 1.
//
//EOP
//-----------------------------------------------------------------------------
  int rc;
  std::stringstream msgbuf;
  ESMCI::LogErr *whichLog;
  timeval  end_time;
  double  elapsed_time;

  // TODO: this should be settable by the user
  whichLog = &ESMC_LogDefault;

  if (file == NULL) {
    msgbuf << "FAIL " << __FILE__ << ", line " << __LINE__ 
           << ", null filename passed to ESMCI::TestEnd()";
    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());
    return(ESMF_FAILURE);
  }

  msgbuf << "Ending Test, file " << file << ", line " << line;
  whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
  if (!only)
    fprintf(stderr, "%s\n", msgbuf.str().c_str());

  rc = ESMCI_Finalize();
  if (rc != ESMF_SUCCESS) {
    msgbuf << "FAIL: " << file << ", line" << line << ", Finalizing ESMF";
    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());
    return(rc);
  }

  // Calculate & print test elapsed time.
#if !defined (ESMF_OS_MinGW)
  gettimeofday(&end_time, NULL);
  elapsed_time = (end_time.tv_sec - start_time.tv_sec) ;
  elapsed_time += (end_time.tv_usec - start_time.tv_usec) / 1000.0;   // us to ms
#else
  elapsed_time = 0.0;
#endif
  msgbuf << " PET " << PETnum << " Test Elapsed Time " 
         << elapsed_time << " msec.";
  fprintf(stdout, "%s\n", msgbuf.str().c_str());
 
  return(ESMF_SUCCESS);

}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//BOP
// !IROUTINE:  ESMCI::TestMaxPETs() - Verify there are not too many PETs
//
// !INTERFACE:   
bool TestMaxPETs(
//
// !RETURN VALUE:
//    true or false
//
// !ARGUMENTS:
  int petCount,   // in - the maximum acceptable number of PETs
  char *file,     // in - test filename
  int line,       // in - test line number in test filename
  int only) {     // in - if set to 0, print on stderr also
// 
// !DESCRIPTION:
//    Returns true if there are not more than the specified number of PETs.
//    If {\tt only} is zero, also print same message to stderr as well
//    as the normal output on stdout.  The default for {\tt only} is 1.
//
//EOP
//-----------------------------------------------------------------------------
  int rc;
  ESMCI::VM *globalVM;
  std::stringstream msgbuf, failMsg;
  int numPETs;
  ESMCI::LogErr *whichLog;

  // TODO: this should be settable by the user
  whichLog = &ESMC_LogDefault;

  if (file == NULL) {
    msgbuf << "FAIL " << __FILE__ << ", line " << __LINE__ 
           << ", null filename passed to ESMCI::TestMaxPETs()";
    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());
    return (false);
  }

  globalVM = ESMCI::VM::getGlobal(&rc);
  if ((globalVM == NULL) || (rc != ESMF_SUCCESS)) {
    msgbuf << "FAIL  rc=" << rc << ", " << file << ", line " << line
           << ", Unable to get GlobalVM";

    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());

    return (false);
  }

  numPETs = globalVM->getPetCount();

  if (numPETs > petCount) {
    failMsg << "These tests must not run on more than " << petCount 
            << " processors.";
    msgbuf << "SKIP  " << failMsg.str() << ", " << file << ", line " << line;
    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());

    return (false);
  }

  return (true);

}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//BOP
// !IROUTINE:  ESMCI::TestMinPETs() - Verify there are not too few PETs
//
// !INTERFACE:   
bool TestMinPETs(
//
// !RETURN VALUE:
//    true or false
//
// !ARGUMENTS:
  int petCount,   // in - the minimum acceptable number of PETs
  char *file,     // in - test filename
  int line,       // in - test line number in test filename
  int only) {     // in - if set to 0, print on stderr also
// 
// !DESCRIPTION:
//    Returns true if there are at least petCount PETs.
//    If {\tt only} is zero, also print same message to stderr as well
//    as the normal output on stdout.  The default for {\tt only} is 1.
//
//EOP
//-----------------------------------------------------------------------------
  int rc;
  ESMCI::VM *globalVM;
  std::stringstream msgbuf, failMsg;
  int numPETs;
  ESMCI::LogErr *whichLog;

  // TODO: this should be settable by the user
  whichLog = &ESMC_LogDefault;

  if (file == NULL) {
    msgbuf << "FAIL " << __FILE__ << ", line " << __LINE__ 
           << ", null filename passed to ESMCI::TestMinPETs()";
    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());
    return (false);
  }

  globalVM = ESMCI::VM::getGlobal(&rc);
  if ((globalVM == NULL) || (rc != ESMF_SUCCESS)) {
    msgbuf << "FAIL  rc=" << rc << ", " << file << ", line " << line
           << ", Unable to get GlobalVM";

    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());

    return (false);
  }

  numPETs = globalVM->getPetCount();

  if (numPETs < petCount) {
    failMsg << "These tests must not run on less than " << petCount 
            << " processors.";
    msgbuf << "SKIP  " << failMsg.str() << ", " << file << ", line " << line;
   whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());

    return (false);
  }

  return (true);

}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//BOP
// !IROUTINE:  ESMCI::TestNumPETs() - Verify there are exactly this number of PETs
//
// !INTERFACE:   
bool TestNumPETs(
//
// !RETURN VALUE:
//    true or false
//
// !ARGUMENTS:
  int petCount,   // in - the exact acceptable number of PETs
  char *file,     // in - test filename
  int line,       // in - test line number in test filename
  int only) {     // in - if set to 0, print on stderr also
// 
// !DESCRIPTION:
//    Returns true only if there are exactly the number of requested PETs.
//    If {\tt only} is zero, also print same message to stderr as well
//    as the normal output on stdout.  The default for {\tt only} is 1.
//
//EOP
//-----------------------------------------------------------------------------
  int rc;
  ESMCI::VM *globalVM;
  std::stringstream msgbuf, failMsg;
  int numPETs;
  ESMCI::LogErr *whichLog;

  // TODO: this should be settable by the user
  whichLog = &ESMC_LogDefault;

  if (file == NULL) {
    msgbuf << "FAIL " << __FILE__ << ", line " << __LINE__ 
           << ", null filename passed to ESMCI::TestNumPETs()";
    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());
    return (false);
  }

  globalVM = ESMCI::VM::getGlobal(&rc);
  if ((globalVM == NULL) || (rc != ESMF_SUCCESS)) {
    msgbuf << "FAIL  rc=" << rc << ", " << file << ", line " << line
           << ", Unable to get GlobalVM";

    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());

    return (false);
  }

  numPETs = globalVM->getPetCount();

  if (numPETs != petCount) {
    failMsg << "These tests must not run on exactly " << petCount 
            << " processors.";
    msgbuf << "SKIP  " << failMsg.str() << ", " << file << ", line " << line;
    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());

    return (false);
  }

  return (true);

}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//BOP
// !IROUTINE:  ESMCI::TestStart() - Initialize the framework, print a standard msg
//
// !INTERFACE:   
int TestStart(
//
// !RETURN VALUE:
//    ESMF_SUCCESS or ESMF_FAILURE
//
// !ARGUMENTS:
  const char *file, // in - test filename
  int line,         // in - test line number in test filename
  int only) {       // in - if set to 0, print on stderr also
// 
// !DESCRIPTION:
//    Initializes the framework, prints out the standard messages needed
//    by the testing scripts.
//    If {\tt only} is zero, also print same message to stderr as well
//    as the normal output on stdout.  The default for {\tt only} is 1.
//
//EOP
//-----------------------------------------------------------------------------
  int rc, rc2;
  ESMCI::VM *globalVM;
  std::stringstream msgbuf;
  int numPETs;
  ESMCI::LogErr *whichLog;

  // TODO: this should be settable by the user
  whichLog = &ESMC_LogDefault;


  if (file == NULL) {
    msgbuf << "FAIL " << __FILE__ << ", line " << __LINE__ 
           << ", null filename passed to ESMCI::TestStart()";
    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());
    return(ESMF_FAILURE);
  }

  const char *underScore = strchr(file, '_');
  if (underScore == NULL) underScore = file-1;
  const char *period = strrchr(file, '.');
  int numChars = period - underScore;
  std::string logFileName = std::string(underScore+1, numChars) + "Log";
  
  rc = ESMC_Initialize(&rc2,
    ESMC_InitArgDefaultConfigFilename(NULL),
    ESMC_InitArgDefaultCalKind(ESMC_CALKIND_NOCALENDAR),
    ESMC_InitArgLogFilename(logFileName.c_str()),
    ESMC_InitArgLogKindFlag(ESMC_LOGKIND_MULTI),
    ESMC_ArgLast);
  if (rc2 != ESMF_SUCCESS)
    rc = rc2;
  if (rc != ESMF_SUCCESS) {
    msgbuf << "FAIL  rc=" << rc << ", " << file << ", line " << line
           << ", Unable to initialize ESMF";
    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());
    return(rc);
  }
  
  ESMC_LogSet(1);

  // Get test start time
#if !defined (ESMF_OS_MinGW)
  gettimeofday(&start_time, NULL);
#endif

  globalVM = ESMCI::VM::getGlobal(&rc);
  if ((globalVM == NULL) || (rc != ESMF_SUCCESS)) {
    msgbuf << "FAIL  rc=" << rc << ", " << file << ", line " << line
           << ", Unable to get GlobalVM";


    whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
    if (!only)
      fprintf(stderr, "%s\n", msgbuf.str().c_str());

    return (false);
  }

  numPETs = globalVM->getPetCount();
  PETnum = globalVM->getLocalPet();

  msgbuf << "Beginning Test, file " << file << ", line " << line;
  whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
  if (!only)
    fprintf(stderr, "%s\n", msgbuf.str().c_str());

  msgbuf << "NUMBER_OF_PROCESSORS " << numPETs;
  whichLog->Write(msgbuf, ESMC_LOGMSG_INFO);
  if (!only)
    fprintf(stderr, "%s\n", msgbuf.str().c_str());
 
  return(ESMF_SUCCESS);

} // end ESMC_TestStart
//-----------------------------------------------------------------------------

} // namespace ESMCI


//-----------------------------------------------------------------------------
extern "C"{
  void FTN_X(c_esmc_printpassflush)(){
    printf("PASS: \n");
    fflush(stdout);
  }
}
//-----------------------------------------------------------------------------

