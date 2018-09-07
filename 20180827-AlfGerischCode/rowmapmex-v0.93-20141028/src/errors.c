/* description of errors, included in rowmapM.c */

/* general Errors */
case 1:
  mexPrintf("Only 0,2,3 or 4 output arguments are possible, but %i were requested.\n",i1);
  msg="Invalid number of output arguments";break;

case 2:
  mexPrintf("Only 3 or 4 input arguments are possible, but %i were passed.\n",i1);
  msg="Invalid number of input arguments";break;

case 3:
  mexPrintf("1st argument has to be a string (function name), a function handle, or an inline function.\n");
  msg="1. Arg: Function (String,handle,inline) expected";break;

case 4:
  mexPrintf("1st argument has characters, but needs to be ONE string. ");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if (i2!=1) 
      mexPrintf("Not a matrix containing strings.");
  mexPrintf("\n");
  msg="1. Arg: only ONE string expected";break;

case 5:
  mexPrintf("2nd argument has to be a double row vector.\n");
  msg="2. Arg: double row vector expected";break;

case 6:
  mexPrintf("2nd argument has to be a double row vector. ");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if (i2!=1) mexPrintf("Not a double matrix");
  mexPrintf("\n");
  msg="2. Arg: (1,k) double vector expected";break;

case 7:
  mexPrintf("2nd argument has to be a double vector with at least 2 components.\n");
  mexPrintf("The length of the vector found is %i.\n",i1);
  msg="2. Arg: (1,k) double vector (k>=2) expected";break;

case 8:
  mexPrintf("3rd argument has to be a double column vector.\n");
  msg="3. Arg: double column vector expected";break;

case 9:
  mexPrintf("3rd argument has to be a double column vector.\n");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if (i2!=1) mexPrintf("Not a double matrix");
  mexPrintf("\n");
  msg="3. Arg: (d,1) double vector expected";break;

case 10:
  mexPrintf("3rd argument has to be a double column vector with at least one component.\n");
  mexPrintf("The length of the vector found was %i.\n",paramGlobal.d);
  msg="3. Arg: (d,1) double vector (d>=1) expected";break;    

case 11:
  mexPrintf("4th argument has to be a struct.\n");
  msg="4. Arg: struct expected";break;

case 12:
  mexPrintf("3rd arg: initial time is greater or equal final time.\n");
  msg="3. Arg: t0>=tend";break;  

case 14:
  mexPrintf("2nd argument has to be strictly ordered.\n");
  if (i1>0)
    mexPrintf("Because of t0<tend, the 2nd argument has to be strictly increasing.\n");
  else
    mexPrintf("Because of t0>tend, the 2nd argument has to be strictly descending.\n");
  msg="2. Arg: has to be strictly ordered";break;

case 15:
  mexPrintf("1st argument contained inline-funcs or function handles, but I need only ONE. ");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if ((i2!=1) || (i3!=1)) 
      mexPrintf("Not a matrix containing inline-funcs or function handles.");
  mexPrintf("\n");
  msg="1. Arg: only ONE function expected";break;

case 16:
  mexPrintf("Concering Option '%s':\n",OPT_FUNCCALLMETHOD);
  mexPrintf("Requirement: '%s'==0 or '%s'==1\n",OPT_FUNCCALLMETHOD,OPT_FUNCCALLMETHOD);
  msg="Invalid function call method";break;

case 17:
  mexPrintf("Concering Option '%s':\n",OPT_FUNCCALLMETHOD);
  mexPrintf("'%s'==0 was chosen. Hence all Matlab-functions must be given as string.\n",
    OPT_FUNCCALLMETHOD);
  mexPrintf("But the RHSFun was not supplied as string.\n");
  msg="function call method 0 => all user-supplied functions must be given as strings";break;

case 18:
  mexPrintf("Concering Option '%s':\n",OPT_FUNCCALLMETHOD);
  mexPrintf("'%s'==0 was chosen. Hence all Matlab-functions must be given as string.\n",
    OPT_FUNCCALLMETHOD);
  mexPrintf("But the OutputFun was not supplied as string.\n");
  msg="function call method 0 => all user-supplied functions must be given as strings";break;

case 19:
  mexPrintf("Concering Option '%s':\n",OPT_FUNCCALLMETHOD);
  mexPrintf("'%s'==0 was chosen. Hence all Matlab-functions must be given as string.\n",
    OPT_FUNCCALLMETHOD);
  mexPrintf("But the PostStepFun was not supplied as string.\n");
  msg="function call method 0 => all user-supplied functions must be given as strings";break;

case 20:
  mexPrintf("Concering Option '%s':\n",OPT_FUNCCALLMETHOD);
  mexPrintf("'%s'==0 was chosen. Hence all Matlab-functions must be given as string.\n",
    OPT_FUNCCALLMETHOD);
  mexPrintf("But the JacVFun was not supplied as string.\n");
  msg="function call method 0 => all user-supplied functions must be given as strings";break;

case 21:
  mexPrintf("Concering Option '%s':\n",OPT_FUNCCALLMETHOD);
  mexPrintf("'%s'==0 was chosen. Hence all Matlab-functions must be given as string.\n",
    OPT_FUNCCALLMETHOD);
  mexPrintf("But the FdtFun was not supplied as string.\n");
  msg="function call method 0 => all user-supplied functions must be given as strings";break;


case 81:
  mexPrintf("Concerning Option '%s':\n",OPT_INITIALSS);
  mexPrintf("Requirement: 0<'%s'.\n",OPT_INITIALSS);
  mexPrintf("But I found '%s'=%e.\n",OPT_INITIALSS,d1);
  msg="invalid initial step size";break;

case 82:
  mexPrintf("Concerning Option '%s':\n",OPT_CONTOUTPUTINTERP);
  mexPrintf("Requirement: '%s'= 0 or 1.\n",OPT_CONTOUTPUTINTERP);
  mexPrintf("But I found '%s'=%i.\n",OPT_CONTOUTPUTINTERP,i1);
  msg="invalid interpolation method";break;

/* Errors concering WORK and IWORK-options */
case 101:
  mexPrintf("Concerning Option '%s':\n",OPT_MAXSTEPS);
  mexPrintf("Requirement: 0<='%s'.\n",OPT_MAXSTEPS);
  mexPrintf("But I found '%s'=%i.\n",OPT_MAXSTEPS,i1);
  msg="invalid maximal number of steps";break;

case 102:
  mexPrintf("Concerning Option '%s':\n",OPT_ROWMAPCOEFFS);
  mexPrintf("Requirement: '%s'=1,2,3,4,5, or 6.\n",OPT_ROWMAPCOEFFS);
  mexPrintf("But I found '%s'=%i.\n",OPT_ROWMAPCOEFFS,i1);
  msg="invalid identifier for ROWMAP coefficient set";break;

case 103:
  mexPrintf("Concerning Option '%s':\n",OPT_EPS);
  mexPrintf("Requirement: 1e-35<'%s'<1.0.\n",OPT_EPS);
  mexPrintf("But I found '%s'=%e.\n",OPT_EPS,d1);
  msg="invalid precision eps";break;

case 104:
  mexPrintf("Concerning Option '%s':\n",OPT_RHO);
  mexPrintf("Requirement: 1e-4<'%s'<1.\n",OPT_RHO);
  mexPrintf("But I found '%s'=%e.\n",OPT_RHO,d1);
  msg="invalid value for safety factor in step size prediction";break;

case 105:
  mexPrintf("Concerning Option '%s':\n",OPT_SSMINSEL);
  mexPrintf("Requirement: 0<'%s'.\n",OPT_SSMINSEL);
  mexPrintf("But I found '%s'=%e.\n",OPT_SSMINSEL,d1);
  msg="invalid lower bound for step size prediction";break;

case 106:
  mexPrintf("Concerning Option '%s':\n",OPT_SSMAXSEL);
  mexPrintf("Requirement: 0<'%s'.\n",OPT_SSMAXSEL);
  mexPrintf("But I found '%s'=%e.\n",OPT_SSMAXSEL,d1);
  msg="invalid upper bound for step size prediciton";break;

case 121:
  mexPrintf("Concerning Option '%s':\n",OPT_MAXKRYLOV);
  mexPrintf("Requirement: 0<'%s'.\n",OPT_MAXKRYLOV);
  mexPrintf("But I found '%s'=%i.\n",OPT_MAXKRYLOV,i1);
  msg="invalid maximal number of Krylov steps";break;

case 122:
  mexPrintf("Concerning Option '%s':\n",OPT_KTOL);
  mexPrintf("Requirement: 1e-12<'%s'<=1.\n",OPT_KTOL);
  mexPrintf("But I found '%s'=%e.\n",OPT_KTOL,d1);
  msg="invalid value for Krylov tolerance";break;

case 123:
  mexPrintf("Concerning Option '%s':\n",OPT_MAXSTP);
  mexPrintf("Requirement: '%s'>=0.\n",OPT_MAXSTP);
  mexPrintf("But I found '%s'=%e.\n",OPT_MAXSTP,d1);
  msg="invalid value for maximum step size";break;

/* Errors concerning RHSFun() */
case 301:
  mexPrintf("Problem calling RHSFun().\n");
  mexPrintf("I have not got return values.\n");
  msg="RHSFun(): without return values";break;

case 302:
  mexPrintf("Problem calling RHSFun().\n");
  mexPrintf("Return value was not a double vector.\n");
  msg="RHSFun(): return value was not a double vector";break;

case 303:
  mexPrintf("Problem calling RHSFun().\n");
  mexPrintf("Return value was not a (d,1) or (1,d) double Vector.\n");
  mexPrintf("A (%i,%i) double matrix was returned.\n",i1,i2);
  msg="RHSFun(): return value was not a (d,1) or (1,d) double vector";break;

/* Errors concerning JacVFun() */
case 311:
  mexPrintf("Problem calling JacVFun().\n");
  mexPrintf("I have not got return values.\n");
  msg="JacVFun(): without return values";break;

case 312:
  mexPrintf("Problem calling JacVFun().\n");
  mexPrintf("Return value was not a double vector.\n");
  msg="JacVFun(): return value was not a double vector";break;

case 313:
  mexPrintf("Problem calling JacVFun().\n");
  mexPrintf("Return value was not a (d,1) or (1,d) double Vector.\n");
  mexPrintf("A (%i,%i) double matrix was returned.\n",i1,i2);
  msg="JacVFun(): return value was not a (d,1) or (1,d) double vector";break;

/* Errors concerning FdtFun() */
case 321:
  mexPrintf("Problem calling FdtFun().\n");
  mexPrintf("I have not got return values.\n");
  msg="FdtFun(): without return values";break;

case 322:
  mexPrintf("Problem calling FdtFun().\n");
  mexPrintf("Return value was not a double vector.\n");
  msg="FdtFun(): return value was not a double vector";break;

case 323:
  mexPrintf("Problem calling FdtFun().\n");
  mexPrintf("Return value was not a (d,1) or (1,d) double Vector.\n");
  mexPrintf("A (%i,%i) double matrix was returned.\n",i1,i2);
  msg="FdtFun(): return value was not a (d,1) or (1,d) double vector";break;

/* Errors concerning PostStepFun() */
case 331:
  mexPrintf("Problem calling PostStepFun().\n");
  mexPrintf("I have not got return value 2.\n");
  msg="PostStepFun(): without return value 2";break;

case 332:
  mexPrintf("Problem calling PostStepFun().\n");
  mexPrintf("Return value 2 was not a double vector.\n");
  msg="PostStepFun(): return value 2 was not a double vector";break;

case 333:
  mexPrintf("Problem calling PostStepFun().\n");
  mexPrintf("Return value 2 was not a (d,1) or (1,d) double Vector.\n");
  mexPrintf("A (%i,%i) double matrix was returned.\n",i1,i2);
  msg="PostStepFun(): return value 2 was not a (d,1) or (1,d) double vector";break;

case 341:
  mexPrintf("Problem calling PostStepFun().\n");
  mexPrintf("I have not got return value 3.\n");
  msg="PostStepFun(): without return value 3";break;

case 342:
  mexPrintf("Problem calling PostStepFun().\n");
  mexPrintf("Return value 3 was not a double vector.\n");
  msg="PostStepFun(): return value 3 was not a double vector";break;

case 343:
  mexPrintf("Problem calling PostStepFun().\n");
  mexPrintf("Return value 3 was not a (d,1) or (1,d) double Vector.\n");
  mexPrintf("A (%i,%i) double matrix was returned.\n",i1,i2);
  msg="PostStepFun(): return value 3 was not a (d,1) or (1,d) double vector";break;

/* Errors concering optional user supplied function */
case 401:
  mexPrintf("Concering option '%s':\n",OPT_OUTPUTFUNCTION);
  mexPrintf("Option contained inline-funcs or function handles, but I need only ONE. ");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if ((i2!=1) || (i3!=1)) 
      mexPrintf("Not a matrix containing inline-funcs or function handles.");
  mexPrintf("\n");
  msg="only ONE function expected";break;
case 402:
  mexPrintf("Concering option '%s':\n",OPT_POSTSTEPFUNCTION);
  mexPrintf("Option contained inline-funcs or function handles, but I need only ONE. ");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if ((i2!=1) || (i3!=1)) 
      mexPrintf("Not a matrix containing inline-funcs or function handles.");
  mexPrintf("\n");
  msg="only ONE function expected";break;
case 403:
  mexPrintf("Concering option '%s':\n",OPT_JACVFUNCTION);
  mexPrintf("Option contained inline-funcs or function handles, but I need only ONE. ");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if ((i2!=1) || (i3!=1)) 
      mexPrintf("Not a matrix containing inline-funcs or function handles.");
  mexPrintf("\n");
  msg="only ONE function expected";break;
case 404:
  mexPrintf("Concering option '%s':\n",OPT_FDTFUNCTION);
  mexPrintf("Option contained inline-funcs or function handles, but I need only ONE. ");
  if (i1!=2) 
    mexPrintf("Not a 0-, 1- or >=3-dimensional thing.");
  else
    if ((i2!=1) || (i3!=1)) 
      mexPrintf("Not a matrix containing inline-funcs or function handles.");
  mexPrintf("\n");
  msg="only ONE function expected";break;
case 405:
  mexPrintf("Concering option '%s':\n",OPT_RETURNMODE);
  mexPrintf("Allowed values are '%s'=0,1,2,3 or 4 but '%s'=%i supplied.\n",
	    OPT_RETURNMODE,OPT_RETURNMODE,i1);
  msg="invalied value for return mode";break;
case 406:
  mexPrintf("Concering option '%s':\n",OPT_OUTPUTCALLMODE);
  mexPrintf("Allowed values are '%s'=0,2,3 or 4 but '%s'=%i supplied.\n",
	    OPT_OUTPUTCALLMODE,OPT_OUTPUTCALLMODE,i1);
  msg="invalied value for output call mode";break;
case 407:
  mexPrintf("Concering option '%s':\n",OPT_POSTSTEPCALLMODE);
  mexPrintf("Allowed values are '%s'=0 or 3 but '%s'=%i supplied.\n",
	    OPT_POSTSTEPCALLMODE,OPT_POSTSTEPCALLMODE,i1);
  msg="invalied value for post-step call mode";break;

/* Errors that should never occur */
case 1001:
  msg="internal error: that should not have happened.";
  break;

case 1002:
  msg="internal error: unknown funcCallMethod";
  break;

case 1101:
  mexPrintf("ERROR: ROWMAP terminated with IDID = %i .\n", i1);
  msg="internal error: unknown IDID value from ROWMAP.";
  break;
