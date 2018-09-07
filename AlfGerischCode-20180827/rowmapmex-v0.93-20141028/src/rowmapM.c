/*
 * ROWMAP-MEX-Interface by A. Gerisch
 *   Adapted from: DOPRI5MEX-Interface by C. Ludwig (Ludwig_C@gmx.de)
 *                 Version:  Id: dopri5Mex.c 278 2006-01-11 16:46:40Z luchr 
 */
#define ROWMAPMexVersion "08 June 2009 (v. 0.93)"
/*
 *   Questions, requests, problems, notes,
 *   remarks and bugs to
 *   alf.gerisch@mathematik.uni-halle.de
 */
#include <math.h>
#include "mex.h"
#include "string.h"
#include "inttype.h"
#include "options.h"
#include "rowmapM.h"


static SOptionSettings optSet;
static SParameterGlobal paramGlobal;
static SParameterOptions paramOptions;
static SParameterROWMAP paramROWMAP;
static SParameterRHSFun paramRHSFun;
static SParameterJacVFun paramJacVFun;
static SParameterFdtFun paramFdtFun;
static SParameterIOutputFun paramIOutputFun;


/* 
 * Helper Functions
 */
static char ismxArrayString (const mxArray *arr) {
  if (arr==NULL) return (char)0;
  if (mxIsChar(arr)) return (char)1;
  return (char)0;
}

static char ismxArrayFunction (const mxArray *arr) {
  if (arr==NULL) return (char)0;
  return (mxGetClassID(arr)==mxFUNCTION_CLASS)?(char)1:(char)0;
}


static char ismxArrayInline (const mxArray *arr) {
  if (arr==NULL) return (char)0;
  return (strcmp(mxGetClassName(arr),"inline")==0)?(char)1:(char)0;
}

static void clearList (PListElement current) {
  PListElement next;
  
  while (current!=NULL) {
    next=current->next;
    if (current->values!=NULL)
      {mxFree(current->values);current->values=NULL;}
    mxFree(current);

    current=next;
  }
}


static void initVars (void) {
  /* Option settings */
  optSet.warnMiss=0;
  optSet.warnType=1;
  optSet.warnSize=1;     
  
  /* global parameters */
  paramGlobal.tspanPointer=NULL;

  /* Parameters for Options*/
  paramOptions.opt=NULL;
  paramOptions.optCreated=0;    
  
  /* parameters for ROWMAP */
  paramROWMAP.yStart=NULL;
  paramROWMAP.RTOL=NULL;
  paramROWMAP.ATOL=NULL;
  paramROWMAP.WORK=NULL;
  paramROWMAP.IWORK=NULL;
  paramROWMAP.RPAR=NULL;
  paramROWMAP.IPAR=NULL;
  
  /* parameters for RHSFun */
  paramRHSFun.RHSFun=NULL;
  paramRHSFun.RHSFunH=NULL;
  paramRHSFun.tArg=NULL;
  paramRHSFun.yArg=NULL;

  /* parameters for JacVFun */
  paramJacVFun.JacVFun=NULL;
  paramJacVFun.JacVFunH=NULL;
  paramJacVFun.tArg=NULL;
  paramJacVFun.yArg=NULL;
  paramJacVFun.vArg=NULL;

  /* parameters for FdtFun */
  paramFdtFun.FdtFun=NULL;
  paramFdtFun.FdtFunH=NULL;
  paramFdtFun.tArg=NULL;
  paramFdtFun.yArg=NULL;

  /* parameters for internal output function */
  paramIOutputFun.tyList.values=NULL;
  paramIOutputFun.tyList.next=NULL;
  paramIOutputFun.tyLastElement=&paramIOutputFun.tyList;
  paramIOutputFun.tyNoOfElements=0;
  paramIOutputFun.OutputFun=NULL;
  paramIOutputFun.OutputFunH=NULL;
  paramIOutputFun.PostStepFun=NULL;
  paramIOutputFun.PostStepFunH=NULL;
  paramIOutputFun.tArg=NULL;
  paramIOutputFun.yArg=NULL;
  paramIOutputFun.fArg=NULL;
}


static void doneVars (void) {

  /* global parameters */
  /* It's not allowed to free paramGlobal.tPointer, i.e. the tspan argument 
     of the rowmap call because it belogns to the caller */

  /* Parameters for Options*/
  if ((paramOptions.optCreated) && (paramOptions.opt!=NULL)) {
    mxDestroyArray((mxArray*)paramOptions.opt);
    paramOptions.opt=NULL;
    paramOptions.optCreated=0;
  }
        
  /* parameters for ROWMAP */
  if (paramROWMAP.yStart!=NULL)
    {mxFree(paramROWMAP.yStart);paramROWMAP.yStart=NULL;}
  if (paramROWMAP.RTOL!=NULL)
    {mxFree(paramROWMAP.RTOL);paramROWMAP.RTOL=NULL;}
  if (paramROWMAP.ATOL!=NULL)
    {mxFree(paramROWMAP.ATOL);paramROWMAP.ATOL=NULL;}
  if (paramROWMAP.WORK!=NULL)
    {mxFree(paramROWMAP.WORK);paramROWMAP.WORK=NULL;}
  if (paramROWMAP.IWORK!=NULL)
    {mxFree(paramROWMAP.IWORK);paramROWMAP.IWORK=NULL;}
  if (paramROWMAP.RPAR!=NULL)
    {mxFree(paramROWMAP.RPAR);paramROWMAP.RPAR=NULL;}
  if (paramROWMAP.IPAR!=NULL)
    {mxFree(paramROWMAP.IPAR);paramROWMAP.IPAR=NULL;}
    
  /* parameters for RHSFun */
  if (paramRHSFun.RHSFun!=NULL)
    {mxFree(paramRHSFun.RHSFun);paramRHSFun.RHSFun=NULL;}
  paramRHSFun.RHSFunH=NULL; /* belongs to caller => don't free */
  if (paramRHSFun.tArg!=NULL)
    {mxDestroyArray(paramRHSFun.tArg);paramRHSFun.tArg=NULL;}
  if (paramRHSFun.yArg!=NULL)
    {mxDestroyArray(paramRHSFun.yArg);paramRHSFun.yArg=NULL;}
    
  /* parameters for JacVFun */
  if (paramJacVFun.JacVFun!=NULL)
    {mxFree(paramJacVFun.JacVFun);paramJacVFun.JacVFun=NULL;}
  paramJacVFun.JacVFunH=NULL; /* belongs to caller => don't free */
  if (paramJacVFun.tArg!=NULL)
    {mxDestroyArray(paramJacVFun.tArg);paramJacVFun.tArg=NULL;}
  if (paramJacVFun.yArg!=NULL)
    {mxDestroyArray(paramJacVFun.yArg);paramJacVFun.yArg=NULL;}
  if (paramJacVFun.vArg!=NULL)
    {mxDestroyArray(paramJacVFun.vArg);paramJacVFun.vArg=NULL;}
 
  /* parameters for FdtFun */
  if (paramFdtFun.FdtFun!=NULL)
    {mxFree(paramFdtFun.FdtFun);paramFdtFun.FdtFun=NULL;}
  paramFdtFun.FdtFunH=NULL; /* belongs to caller => don't free */
  if (paramFdtFun.tArg!=NULL)
    {mxDestroyArray(paramFdtFun.tArg);paramFdtFun.tArg=NULL;}
  if (paramFdtFun.yArg!=NULL)
    {mxDestroyArray(paramFdtFun.yArg);paramFdtFun.yArg=NULL;}

  /* parameters for internal output function */
  clearList(paramIOutputFun.tyList.next);
  if (paramIOutputFun.tyList.values!=NULL){
    mxFree(paramIOutputFun.tyList.values);
    paramIOutputFun.tyList.values=NULL;
  }
  paramIOutputFun.tyList.next=NULL;
  paramIOutputFun.tyLastElement=&paramIOutputFun.tyList;
  paramIOutputFun.tyNoOfElements=0;

  if (paramIOutputFun.OutputFun!=NULL)
    {mxFree(paramIOutputFun.OutputFun);paramIOutputFun.OutputFun=NULL;}
  paramIOutputFun.OutputFunH=NULL; /* belongs to caller => don't free */
  if (paramIOutputFun.PostStepFun!=NULL)
    {mxFree(paramIOutputFun.PostStepFun);paramIOutputFun.PostStepFun=NULL;}
  paramIOutputFun.PostStepFunH=NULL; /* belongs to caller => don't free */

  if (paramIOutputFun.tArg!=NULL)
    {mxDestroyArray(paramIOutputFun.tArg);paramIOutputFun.tArg=NULL;}
  if (paramIOutputFun.yArg!=NULL)
    {mxDestroyArray(paramIOutputFun.yArg);paramIOutputFun.yArg=NULL;}
  if (paramIOutputFun.fArg!=NULL)
    {mxDestroyArray(paramIOutputFun.fArg);paramIOutputFun.fArg=NULL;}
}


static void stopMexFunction (INTTYPE errNo,
  INTTYPE i1, INTTYPE i2, INTTYPE i3, INTTYPE i4, double d1) {
  char *msg;

  mexPrintf("Error (%i) [Version: %s]:\n",errNo,ROWMAPMexVersion);

  switch (errNo) {
    #include "errors.c"      
    default: msg="Unknown errornumber";break;
  }
  
  doneVars();
  mexErrMsgTxt(msg);
}

static void addTYtoList (double t, double *y) {
  /* function adds a new entry to the list and stores data there */
  double* dpointer;
  PListElement target;
  
  target=mxMalloc((mwSize)sizeof(SListElement));
  target->next=NULL;
  paramIOutputFun.tyLastElement->next=target;
  paramIOutputFun.tyLastElement=target;
  paramIOutputFun.tyNoOfElements++;

  target->values=mxMalloc((mwSize)((paramGlobal.d+1)*sizeof(double)));
  dpointer=target->values;
  *dpointer=t;
  dpointer++;
  memcpy(dpointer,y,paramGlobal.d*sizeof(double));  
}





static void checkNumberOfArgs (INTTYPE nlhs, INTTYPE nrhs) {
  if ((nlhs == 1) || (nlhs>4)) 
    stopMexFunction(1,nlhs,0,0,0,0);    
  paramIOutputFun.nlhs = nlhs;

  if ((nrhs!=3) && (nrhs!=4)) 
    stopMexFunction(2,nrhs,0,0,0,0);
}

static void processArgs (INTTYPE nrhs, const mxArray* prhs[]) {
  INTTYPE buflen;
  double *dpointer;
  
  /* 1st arg: right-hand side function */
  if (ismxArrayString(prhs[0])) {
    if ((mxGetNumberOfDimensions(prhs[0])!=2) || (mxGetM(prhs[0])!=1))
      stopMexFunction(4,mxGetNumberOfDimensions(prhs[0]),mxGetM(prhs[0]),0,0,0);
    buflen=mxGetN(prhs[0])*sizeof(mxChar)+1;
    paramRHSFun.RHSFun=mxMalloc(buflen);
    mxGetString(prhs[0],paramRHSFun.RHSFun,buflen);
    paramRHSFun.RHSFunH=prhs[0];
  } else 
    if ( (ismxArrayFunction(prhs[0])) || (ismxArrayInline(prhs[0])) ) {
      if ((mxGetNumberOfDimensions(prhs[0])!=2) || 
      (mxGetM(prhs[0])!=1) || (mxGetN(prhs[0])!=1)) 
    stopMexFunction(15,mxGetNumberOfDimensions(prhs[0]),
            mxGetM(prhs[0]),mxGetN(prhs[0]),0,0);
      if (paramRHSFun.RHSFun!=NULL) { /* clear if a string is stored here */
    mxFree(paramRHSFun.RHSFun); 
    paramRHSFun.RHSFun=NULL;
      }
      paramRHSFun.RHSFunH=prhs[0];
    } else
      stopMexFunction(3,0,0,0,0,0);

  
  /* 2nd arg: row-vector containing tspan values */
  if (!mxIsDouble(prhs[1])) 
    stopMexFunction(5,0,0,0,0,0);
  if ((mxGetNumberOfDimensions(prhs[1])!=2) || (mxGetM(prhs[1])!=1))
    stopMexFunction(6,mxGetNumberOfDimensions(prhs[1]),mxGetM(prhs[1]),0,0,0);
  paramGlobal.tspanLength=mxGetN(prhs[1]);
  if (paramGlobal.tspanLength<2) 
    stopMexFunction(7,paramGlobal.tspanLength,0,0,0,0);
  dpointer=mxGetPr(prhs[1]);
  paramGlobal.tspanPointer=dpointer;
  paramROWMAP.tStart=dpointer[0];
  /*mexPrintf(" set tStart = %e \n", paramROWMAP.tStart);*/
  paramROWMAP.tEnd=dpointer[paramGlobal.tspanLength-1];
  if (paramROWMAP.tStart>=paramROWMAP.tEnd) 
    stopMexFunction(12,0,0,0,0,0);
  paramGlobal.direction = 1.0;
  for (buflen=1; buflen<paramGlobal.tspanLength; buflen++,dpointer++)
    if (paramGlobal.direction*(dpointer[0]-dpointer[1])>=0)
      stopMexFunction(14,paramGlobal.direction,0,0,0,0);        
  
  /* 3rd arg: initial vector */
  if (!mxIsDouble(prhs[2])) 
    stopMexFunction(8,0,0,0,0,0);
  if ((mxGetNumberOfDimensions(prhs[2])!=2) || (mxGetN(prhs[2])!=1))
    stopMexFunction(9,mxGetNumberOfDimensions(prhs[2]),mxGetN(prhs[2]),0,0,0);
  paramGlobal.d=mxGetM(prhs[2]);
  if (paramGlobal.d<1) 
    stopMexFunction(10,0,0,0,0,0);
  dpointer=mxGetPr(prhs[2]);
  paramROWMAP.yStart=mxMalloc((mwSize)(paramGlobal.d*sizeof(double)));
  memcpy(paramROWMAP.yStart,dpointer,paramGlobal.d*sizeof(double));
  /* Remark: A COPY of the startvector is made, because
     it will be passed to Fortran code. To be sure, that
     it will not be changed there, a copy will be passed.
     y0 belongs to the caller and it is not allowed 
     (due to Matlab-Contract) to change it. */
  
  /* 4th arg: struct with options */
  if (nrhs==4) {
    if (!mxIsStruct(prhs[3])) 
      stopMexFunction(11,0,0,0,0,0);
    
    paramOptions.opt=prhs[3];
    paramOptions.optCreated=0;
  } else {
    paramOptions.opt=mxCreateStructMatrix(1,1,0,NULL);
    paramOptions.optCreated=1;
  }  
}


static void extractOptionSettings (void) {
  optSet.warnMiss=opt_getIntFromOpt(paramOptions.opt,&optSet,OPT_WARNMISS,0);
  optSet.warnType=opt_getIntFromOpt(paramOptions.opt,&optSet,OPT_WARNTYPE,1);
  optSet.warnSize=opt_getIntFromOpt(paramOptions.opt,&optSet,OPT_WARNSIZE,1);
}

static void prepareWorkArrays (void) {
  INTTYPE i, mx;
  
  /* extract maximum Krylov dimension */
  mx=opt_getIntFromOpt(paramOptions.opt,&optSet,OPT_MAXKRYLOV,70);
  if (mx<=0) stopMexFunction(121,mx,0,0,0,0);

  /* Allocate WORK */
  paramROWMAP.LWORK=10+paramGlobal.d*(mx+11)+mx*(mx+4);
  paramROWMAP.WORK=mxMalloc((mwSize)(paramROWMAP.LWORK*sizeof(double)));
  /* Std-Values */
  for (i=1-1; i<=10-1; i++) paramROWMAP.WORK[i]=0.0;
  
  /* Allocate IWORK */
  paramROWMAP.LIWORK=mx+20;
  paramROWMAP.IWORK=mxMalloc((mwSize)(paramROWMAP.LIWORK*sizeof(INTTYPE)));  
  /* Std-Values */
  for (i=1-1; i<=20-1; i++) paramROWMAP.IWORK[i]=0;

  /* Set maximum Krylov dimension */
  paramROWMAP.IWORK[3-1]=mx; 
  
}


static void extractTOLs (void) {
  INTTYPE m1,n1,m2,n2,res1,res2;
  INTTYPE takeScalar;
  
  res1=opt_getSizeOfOptField(paramOptions.opt,OPT_RTOL,&m1,&n1);
  res2=opt_getSizeOfOptField(paramOptions.opt,OPT_ATOL,&m2,&n2);
  
  takeScalar=0;
  if (((res1!=0) || (res2!=0)) ||
      ((res1==0) && (m1==1)) ||
      ((res2==0) && (m2==1))) {
    takeScalar=1; 
  } else {
    if ((n1!=1) || (n2!=1)) takeScalar=1;
    if (m1!=m2) takeScalar=1;
    if (m1!=paramGlobal.d) takeScalar=1;
  }

  if (takeScalar) {
    paramROWMAP.RTOL=mxMalloc((mwSize)(1*sizeof(double)));
    paramROWMAP.ATOL=mxMalloc((mwSize)(1*sizeof(double)));
    *paramROWMAP.RTOL=
      opt_getDoubleFromOpt(paramOptions.opt,&optSet,OPT_RTOL,1e-3);
    *paramROWMAP.ATOL=
      opt_getDoubleFromOpt(paramOptions.opt,&optSet,OPT_ATOL,1e-6);
    paramROWMAP.ITOL=0;
  } else {
    paramROWMAP.RTOL=mxMalloc((mwSize)(paramGlobal.d*sizeof(double)));
    paramROWMAP.ATOL=mxMalloc((mwSize)(paramGlobal.d*sizeof(double)));
    opt_getDoubleVectorFromOpt(paramOptions.opt,&optSet,OPT_RTOL,
      paramGlobal.d,1,paramROWMAP.RTOL);
    opt_getDoubleVectorFromOpt(paramOptions.opt,&optSet,OPT_ATOL,
      paramGlobal.d,1,paramROWMAP.ATOL);
    paramROWMAP.ITOL=1;
  }
}

static void extractHs (void) {
    double d;
    d = opt_getDoubleFromOpt(paramOptions.opt,&optSet,OPT_INITIALSS,0.0);
    if (!(d>=0.0)) stopMexFunction(81,0,0,0,0,d);

    paramROWMAP.hs = d;
}

static void extractIfcn (void) {
    INTTYPE i;
    i = opt_getIntFromOpt(paramOptions.opt,&optSet,OPT_RHSNONAUTONOMOUS,1);
    if (i==0)
      paramROWMAP.IFCN=0; /* ODE is autonomous */
    else
      paramROWMAP.IFCN=1; /* ODE is non-autonomous */
}

static void extractContOutputInterp (void) {
    INTTYPE i;
    i = opt_getIntFromOpt(paramOptions.opt,&optSet,OPT_CONTOUTPUTINTERP,0);
    if ((i<0) || (i>1))
      stopMexFunction(81,i,0,0,0,0);
    paramIOutputFun.ContOutputInterp=i;
}

static void extractIWORKOpt (void) {
  INTTYPE i;
  
  /* maximum number of steps */
  i=opt_getIntFromOpt(paramOptions.opt,&optSet,OPT_MAXSTEPS,10000);
  if (i<=0) 
    stopMexFunction(101,i,0,0,0,0);
  paramROWMAP.IWORK[1-1]=i;
  
  /* select coefficient set for ROWMAP */
  i=opt_getIntFromOpt(paramOptions.opt,&optSet,OPT_ROWMAPCOEFFS,2);
  if ((i<1) || (i>6))
    stopMexFunction(102,i,0,0,0,0);
  paramROWMAP.IWORK[2-1]=i; 


  /* maximum Krylov dimension*/
  /* IS ALREADY DONE IN prepareWorkArrays !! */
  /*
    i=opt_getIntFromOpt(paramOpt.opt,&optSet,OPT_MAXKRYLOV,70);
    if (i<=0) stopMexFunction(121,i,0,0,0,0);
    paramROWMAP.IWORK[3-1]=i; 
  */

  /* ROWMAP: please don't try to write */
  paramROWMAP.IWORK[4-1]=-1; 

}

static void extractWORKOpt (void) {
  double d;
  
  d=opt_getDoubleFromOpt(paramOptions.opt,&optSet,OPT_SSMINSEL,0.25);
  if (!(d>0.0)) 
    stopMexFunction(105,0,0,0,0,d);
  paramROWMAP.WORK[1-1]=d;

  d=opt_getDoubleFromOpt(paramOptions.opt,&optSet,OPT_SSMAXSEL,2.0);
  if (!(d>0.0)) 
    stopMexFunction(106,0,0,0,0,d);
  paramROWMAP.WORK[2-1]=d;

  d=opt_getDoubleFromOpt(paramOptions.opt,&optSet,OPT_RHO,0.8);
  if (!((d>1e-4) && (d<1.0))) 
    stopMexFunction(104,0,0,0,0,d);
  paramROWMAP.WORK[3-1]=d;

  d=opt_getDoubleFromOpt(paramOptions.opt,&optSet,OPT_EPS,2.3e-16);
  if (!((d>=1e-35) && (d<1.0))) 
    stopMexFunction(103,0,0,0,0,d);
  paramROWMAP.WORK[4-1]=d;
  
  d=opt_getDoubleFromOpt(paramOptions.opt,&optSet,OPT_KTOL,1.0e-1);
  if (!((d>=1e-12) && (d<=1.0))) 
    stopMexFunction(122,0,0,0,0,d);
  paramROWMAP.WORK[5-1]=d;

  d=opt_getDoubleFromOpt(paramOptions.opt,&optSet,OPT_MAXSTP,0.0e0);
  if (!(d>=0.0e0)) 
    stopMexFunction(123,0,0,0,0,d);
  paramROWMAP.WORK[6-1]=d;

}


static void extractOptFuns (void) {
  mxArray *arr;
  INTTYPE i, defaultValue;

  /* extract optional user-supplied Jacobian times vector function */
  if (paramJacVFun.JacVFun!=NULL)
    {mxFree(paramJacVFun.JacVFun);paramJacVFun.JacVFun=NULL;}
  paramJacVFun.JacVFunH=NULL; /* belongs to caller => don't free */
  arr=mxGetField(paramOptions.opt,0,OPT_JACVFUNCTION);
  if ((arr!=NULL) && (!mxIsEmpty(arr))) {
    /* something is provided */
    if ( (ismxArrayFunction(arr)) || (ismxArrayInline(arr)) ) {
      /* we have function handles or inline functions */
      if ((mxGetNumberOfDimensions(arr)!=2) || 
          (mxGetM(arr)!=1) || (mxGetN(arr)!=1)) 
        stopMexFunction(403,mxGetNumberOfDimensions(arr),mxGetM(arr),
            mxGetN(arr),0,0);
      paramJacVFun.JacVFunH=arr;
    }
    else {
      /* we have strings */
      paramJacVFun.JacVFun=opt_getStringFromOpt(paramOptions.opt,&optSet,
                             OPT_JACVFUNCTION,NULL);
      if (paramJacVFun.JacVFun!=NULL) {
        if (strlen(paramJacVFun.JacVFun)==0) {
          mxFree(paramJacVFun.JacVFun);
      paramJacVFun.JacVFun=NULL;
        } else {
          paramJacVFun.JacVFunH=arr;
        }
      }
    }
  }

  if (paramJacVFun.JacVFunH != NULL)
    paramROWMAP.IJACV=1;
  else
    paramROWMAP.IJACV=0;

  /* extract optional user-supplied Fdt function */
  if (paramFdtFun.FdtFun!=NULL)
    {mxFree(paramFdtFun.FdtFun);paramFdtFun.FdtFun=NULL;}
  paramFdtFun.FdtFunH=NULL; /* belongs to caller => don't free */
  arr=mxGetField(paramOptions.opt,0,OPT_FDTFUNCTION);
  if ((arr!=NULL) && (!mxIsEmpty(arr))) {
    /* something is provided */
    if ( (ismxArrayFunction(arr)) || (ismxArrayInline(arr)) ) {
      /* we have function handles or inline functions */
      if ((mxGetNumberOfDimensions(arr)!=2) || 
          (mxGetM(arr)!=1) || (mxGetN(arr)!=1)) 
        stopMexFunction(404,mxGetNumberOfDimensions(arr),mxGetM(arr),
            mxGetN(arr),0,0);
      paramFdtFun.FdtFunH=arr;
    }
    else {
      /* we have strings */
      paramFdtFun.FdtFun=opt_getStringFromOpt(paramOptions.opt,&optSet,
                             OPT_FDTFUNCTION,NULL);
      if (paramFdtFun.FdtFun!=NULL) {
        if (strlen(paramFdtFun.FdtFun)==0) {
          mxFree(paramFdtFun.FdtFun);
      paramFdtFun.FdtFun=NULL;
        } else {
          paramFdtFun.FdtFunH=arr;
        }
      }
    }
  }
  if (paramFdtFun.FdtFunH != NULL)
    paramROWMAP.IFDT=1;
  else
    paramROWMAP.IFDT=0;


  /* extract optional user-supplied output function */
  if (paramIOutputFun.OutputFun!=NULL)
    {mxFree(paramIOutputFun.OutputFun);paramIOutputFun.OutputFun=NULL;}
  paramIOutputFun.OutputFunH=NULL; /* belongs to caller => don't free */
  arr=mxGetField(paramOptions.opt,0,OPT_OUTPUTFUNCTION);
  if ((arr!=NULL) && (!mxIsEmpty(arr))) {
    /* something is provided */
    if ( (ismxArrayFunction(arr)) || (ismxArrayInline(arr)) ) {
      /* we have function handles or inline functions */
      if ((mxGetNumberOfDimensions(arr)!=2) || 
          (mxGetM(arr)!=1) || (mxGetN(arr)!=1)) 
        stopMexFunction(401,mxGetNumberOfDimensions(arr),mxGetM(arr),
            mxGetN(arr),0,0);
      paramIOutputFun.OutputFunH=arr;
    }
    else {
      /* we have strings */
      paramIOutputFun.OutputFun=opt_getStringFromOpt(paramOptions.opt,&optSet,
                             OPT_OUTPUTFUNCTION,NULL);
      if (paramIOutputFun.OutputFun!=NULL) {
        if (strlen(paramIOutputFun.OutputFun)==0) {
          mxFree(paramIOutputFun.OutputFun);
      paramIOutputFun.OutputFun=NULL;
        } else {
          paramIOutputFun.OutputFunH=arr;
        }
      }
    }
  }

  /* extract optional user-supplied post-step function */
  if (paramIOutputFun.PostStepFun!=NULL)
    {mxFree(paramIOutputFun.PostStepFun);paramIOutputFun.PostStepFun=NULL;}
  paramIOutputFun.PostStepFunH=NULL; /* belongs to caller => don't free */
  arr=mxGetField(paramOptions.opt,0,OPT_POSTSTEPFUNCTION);
  if ((arr!=NULL) && (!mxIsEmpty(arr))) {
    /* something is provided */
    if ( (ismxArrayFunction(arr)) || (ismxArrayInline(arr)) ) {
      /* we have function handles or inline functions */
      if ((mxGetNumberOfDimensions(arr)!=2) || 
          (mxGetM(arr)!=1) || (mxGetN(arr)!=1)) 
        stopMexFunction(402,mxGetNumberOfDimensions(arr),mxGetM(arr),
            mxGetN(arr),0,0);
      paramIOutputFun.PostStepFunH=arr;
    }
    else {
      /* we have strings */
      paramIOutputFun.PostStepFun=opt_getStringFromOpt(paramOptions.opt,&optSet,
                             OPT_POSTSTEPFUNCTION,NULL);
      if (paramIOutputFun.PostStepFun!=NULL) {
        if (strlen(paramIOutputFun.PostStepFun)==0) {
          mxFree(paramIOutputFun.PostStepFun);
      paramIOutputFun.PostStepFun=NULL;
        } else {
          paramIOutputFun.PostStepFunH=arr;
        }
      }
    }
  }

  /* extract function call method and check whether we can call 
     all user supplied functions accordingly.
  */
  paramGlobal.funcCallMethod=
    opt_getIntFromOpt(paramOptions.opt,&optSet,OPT_FUNCCALLMETHOD,1);
  switch (paramGlobal.funcCallMethod) {
    case 0: /* use mexCallMATLAB user supplied functions directly, they 
           must be provided as strings  */
      /* check RHSFun() */
      if (paramRHSFun.RHSFun==NULL) 
    stopMexFunction(17,0,0,0,0,0);

      /* check JacVFun() */
      if ((paramJacVFun.JacVFunH!=NULL) && (paramJacVFun.JacVFun==NULL))
    stopMexFunction(20,0,0,0,0,0);

      /* check FdtFun() */
      if ((paramFdtFun.FdtFunH!=NULL) && (paramFdtFun.FdtFun==NULL))
    stopMexFunction(21,0,0,0,0,0);

      /* check OutputFun() */
      if ((paramIOutputFun.OutputFunH!=NULL) && (paramIOutputFun.OutputFun==NULL))
        stopMexFunction(18,0,0,0,0,0);

      /* check PostStepFun() */
      if ((paramIOutputFun.PostStepFunH!=NULL) && (paramIOutputFun.PostStepFun==NULL))
        stopMexFunction(19,0,0,0,0,0);
      break;
    case 1: /* use mexCallMATLAB to call feval */
      break;
    default:
      stopMexFunction(16,0,0,0,0,0);
      break;
  }

  /* extract ReturnMode */
  if (paramIOutputFun.nlhs == 0)
    /* ignore option value ReturnMode */
    paramIOutputFun.ReturnMode = 0;
  else {
    if (paramGlobal.tspanLength == 2)
      defaultValue = 3;
    else
      defaultValue = 2;
    i = opt_getIntFromOpt(paramOptions.opt,&optSet,OPT_RETURNMODE,defaultValue);
    if ((i<0) || (i>4))
      stopMexFunction(405,i,0,0,0,0);
    paramIOutputFun.ReturnMode = i;
  }
  
  /* extract OutputCallMode */
  if (paramIOutputFun.OutputFunH == NULL)
    /* we have no output function to call */
    paramIOutputFun.OutputCallMode = 0;
  else {
    i = opt_getIntFromOpt(paramOptions.opt,&optSet,OPT_OUTPUTCALLMODE,2);
    if ((i<0) || (i==1) || (i>4))
      stopMexFunction(406,i,0,0,0,0);
    paramIOutputFun.OutputCallMode = i;
  }
  

  /* extract PostStepCallMode */
  if (paramIOutputFun.PostStepFunH == NULL)
    /* we have no post-step function to call */
    paramIOutputFun.PostStepCallMode = 0;
  else {
    i = opt_getIntFromOpt(paramOptions.opt,&optSet,OPT_POSTSTEPCALLMODE,3);
    if (!((i==0) || (i==3)))
      stopMexFunction(407,i,0,0,0,0);
    paramIOutputFun.PostStepCallMode = i;  
  }

}


static void extractOptions (void) {
  extractTOLs();
  extractHs();
  extractIfcn();
  extractContOutputInterp();
  extractIWORKOpt();
  extractWORKOpt();
  extractOptFuns();
 }

static void prepareHelpMxArrays (void) {
  INTTYPE d;

  d=paramGlobal.d;
  
  paramRHSFun.tArg=mxCreateDoubleMatrix(1,1,mxREAL);
  paramRHSFun.yArg=mxCreateDoubleMatrix(d,1,mxREAL);

  if (paramROWMAP.IJACV) {
    paramJacVFun.tArg=mxCreateDoubleMatrix(1,1,mxREAL);
    paramJacVFun.yArg=mxCreateDoubleMatrix(d,1,mxREAL);
    paramJacVFun.vArg=mxCreateDoubleMatrix(d,1,mxREAL);
  }

  if (paramROWMAP.IFDT) {
    paramFdtFun.tArg=mxCreateDoubleMatrix(1,1,mxREAL);
    paramFdtFun.yArg=mxCreateDoubleMatrix(d,1,mxREAL);
  }
    
  if ((paramIOutputFun.OutputFunH != NULL) || 
      (paramIOutputFun.PostStepFunH != NULL) ||
      (paramIOutputFun.ReturnMode != 0)) {
    paramIOutputFun.tArg=mxCreateDoubleMatrix(1,1,mxREAL);
    paramIOutputFun.yArg=mxCreateDoubleMatrix(paramGlobal.d,1,mxREAL);    
    paramIOutputFun.fArg=mxCreateDoubleMatrix(paramGlobal.d,1,mxREAL);    
  }
}


void ROWMAP_RHSFun (INTTYPE *n, double *t,
            double *y, double *f, 
            double *rpar, INTTYPE *ipar) {

  INTTYPE d;
  mxArray *rhs[3];
  mxArray *lhs[1];
  
  d=*n;lhs[0]=NULL;

  /* Call by value */
  /* Use always the SAME Matlab tArg and yArg */
  /* IMPORTANT (if the right-hand side is also a MEX-File)
     The right-hand side function must not change the passed 
     mxArrays: the size must not be changed and the memory 
     must not be freed. 
     Hence: the right-hand side function has to take care, 
     that the mxArrays passed have the same size and memory
     when returning. The values in the memory(-block)
     may be overwritten.
     If the right-hand side function is an m-file, MATLAB 
     obeys this restriction automatically. 
  */

  /* set t-value */
  *mxGetPr(paramRHSFun.tArg)= *t;
  /* set y-value */
  memcpy(mxGetPr(paramRHSFun.yArg),y,d*sizeof(double));

  /* call user supplied RHSFun */
  /* mexPrintf("Calling RHSFun().\n"); */
  switch (paramGlobal.funcCallMethod) {
    case 0: 
      rhs[0]=paramRHSFun.tArg; 
      rhs[1]=paramRHSFun.yArg;
      mexCallMATLAB(1,lhs,2,rhs,paramRHSFun.RHSFun);
      /* mexPrintf("Called RHSFun() direct.\n"); */
      break;
    case 1:
      rhs[0]=(mxArray*)paramRHSFun.RHSFunH;
      rhs[1]=paramRHSFun.tArg;
      rhs[2]=paramRHSFun.yArg;
      mexCallMATLAB(1,lhs,3,rhs,"feval"); 
      /* mexPrintf("Called RHSFun() feval.\n"); */
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
  }
  
  /* check return values */
  if (lhs[0]==NULL) 
    stopMexFunction(301,0,0,0,0,0);

  if ((!mxIsDouble(lhs[0])) || (mxGetNumberOfDimensions(lhs[0])!=2))
    stopMexFunction(302,0,0,0,0,0);

  if (!(((mxGetM(lhs[0])==d) && (mxGetN(lhs[0])==1)) ||
        ((mxGetM(lhs[0])==1) && (mxGetN(lhs[0])==d))))
    stopMexFunction(303,mxGetM(lhs[0]),mxGetN(lhs[0]),0,0,0);
    
  /* copy back */
  memcpy(f,mxGetPr(lhs[0]),d*sizeof(double));
  
  /* free memory */
  mxDestroyArray(lhs[0]);  
}

void ROWMAP_JacVFun(INTTYPE *n, double *t,
            double *y, double *v, double *z, 
            double *rpar, INTTYPE *ipar) {

  INTTYPE d;
  mxArray *rhs[4];
  mxArray *lhs[1];
  
  d=*n;lhs[0]=NULL;

  /* Call by value */
  /* Use always the SAME Matlab tArg, yArg, vArg */
  /* IMPORTANT (if the JacVFun() is also a MEX-File)
       --> see ROWMAP_RHSFun
  */

  /* set t-value */
  *mxGetPr(paramJacVFun.tArg)= *t;
  /* set y-value */
  memcpy(mxGetPr(paramJacVFun.yArg),y,d*sizeof(double));
  /* set v-value */
  memcpy(mxGetPr(paramJacVFun.vArg),v,d*sizeof(double));

  /* call user supplied JacVFun */
  /* mexPrintf("Calling JacVFun().\n"); */
  switch (paramGlobal.funcCallMethod) {
    case 0: 
      rhs[0]=paramJacVFun.tArg; 
      rhs[1]=paramJacVFun.yArg;
      rhs[2]=paramJacVFun.vArg;
      mexCallMATLAB(1,lhs,3,rhs,paramJacVFun.JacVFun);
      /* mexPrintf("Called JacVFun() direct.\n"); */
      break;
    case 1:
      rhs[0]=(mxArray*)paramJacVFun.JacVFunH;
      rhs[1]=paramJacVFun.tArg;
      rhs[2]=paramJacVFun.yArg;
      rhs[3]=paramJacVFun.vArg;
      mexCallMATLAB(1,lhs,4,rhs,"feval");
      /* mexPrintf("Called JacVFun() feval.\n"); */
     break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
  }

  /* check return values */
  if (lhs[0]==NULL) 
    stopMexFunction(311,0,0,0,0,0);

  if ((!mxIsDouble(lhs[0])) || (mxGetNumberOfDimensions(lhs[0])!=2))
    stopMexFunction(312,0,0,0,0,0);

  if (!(((mxGetM(lhs[0])==d) && (mxGetN(lhs[0])==1)) ||
        ((mxGetM(lhs[0])==1) && (mxGetN(lhs[0])==d))))
    stopMexFunction(313,mxGetM(lhs[0]),mxGetN(lhs[0]),0,0,0);
    
  /* copy back */
  memcpy(z,mxGetPr(lhs[0]),d*sizeof(double));
  
  /* free memory */
  mxDestroyArray(lhs[0]);  
}

void ROWMAP_FdtFun(INTTYPE *n, double *t,
           double *y, double *ft,
           double *rpar, INTTYPE *ipar) {

  INTTYPE d;
  mxArray *rhs[3];
  mxArray *lhs[1];
  
  d=*n;lhs[0]=NULL;

  /* Call by value */
  /* Use always the SAME Matlab tArg and yArg */
  /* IMPORTANT (if theFdtFun() is also a MEX-File)
       --> see ROWMAP_RHSFun
  */

  /* set t-value */
  *mxGetPr(paramFdtFun.tArg)= *t;
  /* set y-value */
  memcpy(mxGetPr(paramFdtFun.yArg),y,d*sizeof(double));

  /* call user supplied FdtFun */
  /* mexPrintf("Calling FdtFun().\n"); */
  switch (paramGlobal.funcCallMethod) {
    case 0: 
      rhs[0]=paramFdtFun.tArg; 
      rhs[1]=paramFdtFun.yArg;
      mexCallMATLAB(1,lhs,2,rhs,paramFdtFun.FdtFun);
      /* mexPrintf("Called FdtFun() direct.\n"); */
      break;
    case 1:
      rhs[0]=(mxArray*)paramFdtFun.FdtFunH;
      rhs[1]=paramFdtFun.tArg;
      rhs[2]=paramFdtFun.yArg;
      mexCallMATLAB(1,lhs,3,rhs,"feval");
      /* mexPrintf("Called FdtFun() feval.\n"); */
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
  }
  
  /* check return values */
  if (lhs[0]==NULL) 
    stopMexFunction(321,0,0,0,0,0);

  if ((!mxIsDouble(lhs[0])) || (mxGetNumberOfDimensions(lhs[0])!=2))
    stopMexFunction(322,0,0,0,0,0);

  if (!(((mxGetM(lhs[0])==d) && (mxGetN(lhs[0])==1)) ||
        ((mxGetM(lhs[0])==1) && (mxGetN(lhs[0])==d))))
    stopMexFunction(323,mxGetM(lhs[0]),mxGetN(lhs[0]),0,0,0);
    
  /* copy back */
  memcpy(ft,mxGetPr(lhs[0]),d*sizeof(double));
  
  /* free memory */
  mxDestroyArray(lhs[0]);  
}


static void callOutputAndPostStepFun_Init(void) {
  INTTYPE offset;
  double *dpointer;
  mxArray* rhs[5];
  mxArray* lhs[3];

  /* store initial (t,y) value pair */
  switch (paramIOutputFun.ReturnMode) {
  case 2:
  case 3:
  case 4:
    /* if==2,3 or 4: add initial time and initial values to tyList */
    addTYtoList(paramROWMAP.tStart,paramROWMAP.yStart);
    break;
  } /* end of switch (paramROWMAP.ReturnMode) */
    
  /* call OutputFun() */
  switch (paramIOutputFun.OutputCallMode) {
  case 2:
  case 3:
  case 4:
    /* if==2,3 or 4: init call to OutputFun */
    offset=0;
    switch (paramGlobal.funcCallMethod) {
    case 0: break;
    case 1:
      rhs[offset++]=(mxArray*)paramIOutputFun.OutputFunH;
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
    }
    /* store [t0,tend] in first OutputFun() argument */
    rhs[offset]=mxCreateDoubleMatrix(1,2,mxREAL);
    dpointer=mxGetPr(rhs[offset]);
    *dpointer=paramROWMAP.tStart;
    dpointer++;
    *dpointer=paramROWMAP.tEnd;
    offset++;
    /* store y0 in second OutputFun() argument */
    rhs[offset]=mxCreateDoubleMatrix(paramGlobal.d,1,mxREAL);
    memcpy(mxGetPr(rhs[offset]),paramROWMAP.yStart,
       paramGlobal.d*sizeof(double));        
    offset++;
    /* store string 'init' in third OutputFun() argument */
    rhs[offset++]=mxCreateString("init");
    
    switch (paramGlobal.funcCallMethod) {
    case 0: 
      mexCallMATLAB(0,lhs,offset,rhs,paramIOutputFun.OutputFun); 
      while (offset>0) {offset--;mxDestroyArray(rhs[offset]);}
      break;
    case 1:
      mexCallMATLAB(0,lhs,offset,rhs,"feval");
      while (offset>1) {offset--;mxDestroyArray(rhs[offset]);}
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
    }
    break;
  } /* end of switch (paramROWMAP.OutputCallMode) */
  
  /* call PostStepFun() */
  switch (paramIOutputFun.PostStepCallMode) {
  case 3:
    /* if==3: init call to PostStepFun */
    offset=0;
    switch (paramGlobal.funcCallMethod) {
    case 0: break;
    case 1:
      rhs[offset++]=(mxArray*)paramIOutputFun.PostStepFunH;
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
    }
    /* store [t0,tend] in first PostStepFun() argument */
    rhs[offset]=mxCreateDoubleMatrix(1,2,mxREAL);
    dpointer=mxGetPr(rhs[offset]);
    *dpointer=paramROWMAP.tStart;
    dpointer++;
    *dpointer=paramROWMAP.tEnd;
    offset++;
    /* store y0 in second PostStepFun() argument */
    rhs[offset]=mxCreateDoubleMatrix(paramGlobal.d,1,mxREAL);
    memcpy(mxGetPr(rhs[offset]),paramROWMAP.yStart,
       paramGlobal.d*sizeof(double));        
    offset++;
    /* store a [] in third PostStepFun() argument */
    rhs[offset++]=mxCreateDoubleMatrix(0,0,mxREAL);
    /* store string 'init' in fourth PostStepFun() argument */
    rhs[offset++]=mxCreateString("init");
    
    switch (paramGlobal.funcCallMethod) {
    case 0: 
      mexCallMATLAB(0,lhs,offset,rhs,paramIOutputFun.PostStepFun);
      while (--offset>=0) {mxDestroyArray(rhs[offset]);}
      break;
    case 1:
      mexCallMATLAB(0,lhs,offset,rhs,"feval");
      while (--offset>=1) {mxDestroyArray(rhs[offset]);}
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
    }
    break;
  } /* end of switch (paramROWMAP.PostStepCallMode) */
} /* end of static INTTYPE callOutputAndPostStepFun_Init(void) */



static void callOutputAndPostStepFun_Done(void) {
  INTTYPE offset;
  mxArray* rhs[5];
  mxArray* lhs[3];

  /*  nothing to do for paramIOutputFun.ReturnMode */

  /* call OutputFun() */
  switch (paramIOutputFun.OutputCallMode) {
  case 2:
  case 3:
  case 4:
    /* if==2,3 or 4: done call to OutputFun */
    offset=0;
    switch (paramGlobal.funcCallMethod) {
    case 0: break;
    case 1:
      rhs[offset++]=(mxArray*)paramIOutputFun.OutputFunH;
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
    }
    /* store [] in first and second OutputFun() argument */
    rhs[offset++]=mxCreateDoubleMatrix(0,0,mxREAL);
    rhs[offset++]=mxCreateDoubleMatrix(0,0,mxREAL);
    /* store string 'done' in third OutputFun() argument */
    rhs[offset++]=mxCreateString("done");
    
    switch (paramGlobal.funcCallMethod) {
    case 0: 
      mexCallMATLAB(0,lhs,offset,rhs,paramIOutputFun.OutputFun);
      while (offset>0) {offset--;mxDestroyArray(rhs[offset]);}
      break;
    case 1:
      mexCallMATLAB(0,lhs,offset,rhs,"feval");
      while (offset>1) {offset--;mxDestroyArray(rhs[offset]);}
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
    }
    break;
  } /* end of switch (paramROWMAP.OutputCallMode) */
  
  /* call PostStepFun() */
  switch (paramIOutputFun.PostStepCallMode) {
  case 3:
    /* if==3: done call to PostStepFun */
    offset=0;
    switch (paramGlobal.funcCallMethod) {
    case 0: break;
    case 1:
      rhs[offset++]=(mxArray*)paramIOutputFun.PostStepFunH;
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
    }
    
    /* store [] in first, second and third PostStepFun() argument */
    rhs[offset++]=mxCreateDoubleMatrix(0,0,mxREAL);
    rhs[offset++]=mxCreateDoubleMatrix(0,0,mxREAL);
    rhs[offset++]=mxCreateDoubleMatrix(0,0,mxREAL);
    /* store string 'done' in fourth PostStepFun() argument */
    rhs[offset++]=mxCreateString("done");
    
    switch (paramGlobal.funcCallMethod) {
    case 0: 
      mexCallMATLAB(0,lhs,offset,rhs,paramIOutputFun.PostStepFun);
      while (--offset>=0) {mxDestroyArray(rhs[offset]);}
      break;
    case 1:
      mexCallMATLAB(0,lhs,offset,rhs,"feval");
      while (--offset>=1) {mxDestroyArray(rhs[offset]);}
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
    }
    break;
  } /* end of switch (paramROWMAP.PostStepCallMode) */    
} /* end of static INTTYPE callOutputAndPostStepFun_Done(void) */

void ROWMAP_OutputFun(INTTYPE *n, 
              double *told, double *tnew, 
              double *yold, double *ynew,
              double *fold, double *fnew, 
              double *ycon, INTTYPE* intr,
              double *rpar, INTTYPE *ipar){ 
  double tout;
  double *dpointer;
  INTTYPE offset;
  INTTYPE tnewIsInTspan;
  INTTYPE status;
  mxArray* rhs[5];
  mxArray* lhs[3];

  /* If this is the first call to output function then told equals  
   * the initial time. The next output point is the one following
   * the initial time, i.e. paramGlobal.tPointer[1]. The 
   * initial time and the initial values are already stored, if required,
   * in the tyList (during the callOutputAndPostStepFun_Init() call).
   * So we only set paramOutput.tspanPos=1 to indicate the next output 
   * time point. We do not call a user supplied output function 
   * for the initial value because this has been done with the 'init' 
   * call already before the integration started (this is constistent 
   * with the behaviour of Matlab's ode45). 
   */ 
  if (*told==paramGlobal.tspanPointer[0]) {
      paramIOutputFun.tspanPos=1;
  }

  while ((paramIOutputFun.tspanPos < paramGlobal.tspanLength) &&
     (paramGlobal.tspanPointer[paramIOutputFun.tspanPos] < *tnew)) {
    /* we have another possible output time tout in (told,tnew) */
    tout = paramGlobal.tspanPointer[paramIOutputFun.tspanPos];
    /* check if we need to output something at tout */
    if ((paramIOutputFun.ReturnMode == 2) ||
	(paramIOutputFun.ReturnMode == 4) ||
	(paramIOutputFun.OutputCallMode == 2) ||
	(paramIOutputFun.OutputCallMode == 4)) {
      /* interpolate solution to tout */
      dpointer=mxGetPr(paramIOutputFun.yArg);
      switch (paramIOutputFun.ContOutputInterp) {
      case 0:
	ROWCON_(n, &tout, 
		told, tnew,
		yold, ynew,
		fold, fnew, 
		dpointer);
	break;
      case 1:
	ROWLIN_(n, &tout, 
		told, tnew,
		yold, ynew,
		dpointer);
	break;
      default:
	stopMexFunction(1001,0,0,0,0,0);
	break;
      }
    }
    if ((paramIOutputFun.ReturnMode == 2) ||
	(paramIOutputFun.ReturnMode == 4)) {
      addTYtoList(tout,dpointer);
    }
    if ((paramIOutputFun.OutputCallMode == 2) ||
    (paramIOutputFun.OutputCallMode == 4)) {
      offset=0;
      switch (paramGlobal.funcCallMethod) {
      case 0: 
	break;
      case 1:
	rhs[offset++]=(mxArray*)paramIOutputFun.OutputFunH;
	break;
      default: 
	stopMexFunction(1002,0,0,0,0,0);break;
      }

      /* store tout in first OutputFun() argument */
      *mxGetPr(paramIOutputFun.tArg) = tout;
      rhs[offset++]=paramIOutputFun.tArg;
      /* store corresponding y(tout) in second OutputFun() argument; 
     y(tout) is stored in paramIOutputFun.yArg already 
      */
      rhs[offset++]=paramIOutputFun.yArg;
      /* store empty string in third OutputFun() argument */
      rhs[offset++]=mxCreateString("");
      /* call output function */
      switch (paramGlobal.funcCallMethod) {
      case 0: 
    mexCallMATLAB(1,lhs,offset,rhs,paramIOutputFun.OutputFun);
    mxDestroyArray(rhs[offset-1]);
    break;
      case 1:
	mexCallMATLAB(1,lhs,offset,rhs,"feval");
	mxDestroyArray(rhs[offset-1]);
	break;
      default: 
	stopMexFunction(1002,0,0,0,0,0);break;
      }
      /* deal with return value */
      if (lhs[0]==NULL) {
	status=0; /* continue integration */
      } else {
    status=mxGetScalar(lhs[0]);
    mxDestroyArray(lhs[0]);
      }
      if (status!=1) status=0;
      if (status==1) { /* Stop, NOW */
    *intr=-1; 
    return;
      }
    } /* end of call to OutputFun() */
    /* increment to next tspan output time */
    paramIOutputFun.tspanPos++;
  } /* end of 'we have another possible output time tout in (told,tnew)' */

  /* now deal with solution at tnew */
  tnewIsInTspan = 0;
  if ((paramIOutputFun.tspanPos < paramGlobal.tspanLength) &&
      (paramGlobal.tspanPointer[paramIOutputFun.tspanPos] == *tnew)) {
    /* tnew is tspan value */
    tnewIsInTspan = 1;
    paramIOutputFun.tspanPos++;
  }
  /* store solution depending on paramIOutputFun.ReturnMode */
  switch (paramIOutputFun.ReturnMode) {
  case 0:
    break;
  case 1:
    if (*tnew == paramGlobal.tspanPointer[paramGlobal.tspanLength-1])
      addTYtoList(*tnew,ynew);
    break;
  case 2:
    if (tnewIsInTspan == 1)
      addTYtoList(*tnew,ynew);
    break;
  case 3:
  case 4:
    addTYtoList(*tnew,ynew);
    break;
  }
  /* call OutputFun() depending on paramIOutputFun.OutputCallMode */
  if (((paramIOutputFun.OutputCallMode == 2) && (tnewIsInTspan == 1)) ||
      (paramIOutputFun.OutputCallMode == 3) ||
      (paramIOutputFun.OutputCallMode == 4)) {
    offset=0;
    switch (paramGlobal.funcCallMethod) {
    case 0: break;
    case 1:
      rhs[offset++]=(mxArray*)paramIOutputFun.OutputFunH;
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
    }

    /* store tnew in first OutputFun() argument */
    *mxGetPr(paramIOutputFun.tArg) = *tnew;
    rhs[offset++]=paramIOutputFun.tArg;
    /* store ynew in second OutputFun() argument */
    memcpy(mxGetPr(paramIOutputFun.yArg),ynew,paramGlobal.d*sizeof(double));
    rhs[offset++]=paramIOutputFun.yArg;
    /* store empty string in third OutputFun() argument */
    rhs[offset++]=mxCreateString("");
    /* call output function */
    switch (paramGlobal.funcCallMethod) {
    case 0: 
      mexCallMATLAB(1,lhs,offset,rhs,paramIOutputFun.OutputFun);
      mxDestroyArray(rhs[offset-1]);
      break;
    case 1:
      mexCallMATLAB(1,lhs,offset,rhs,"feval");
      mxDestroyArray(rhs[offset-1]);
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
    }
    /* deal with return value */
    if (lhs[0]==NULL) {
      status=0; /* continue integration */
    } else {
      status=mxGetScalar(lhs[0]);
      mxDestroyArray(lhs[0]);
    }
    if (status!=1) status=0;
    if (status==1) { /* Stop, NOW */
      *intr=-1; 
      return;
    }
  } /* end of 'call OutputFun() depending on paramIOutputFun.OutputCallMode' */

  /* call PostStepFun() depending on paramIOutputFun.PostStepCallMode */
  if (paramIOutputFun.PostStepCallMode == 3) {
    offset=0;
    switch (paramGlobal.funcCallMethod) {
    case 0: break;
    case 1:
      rhs[offset++]=(mxArray*)paramIOutputFun.PostStepFunH;
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);break;
    }

    /* store tnew in first PostStepFun() argument */
    *mxGetPr(paramIOutputFun.tArg) = *tnew;
    rhs[offset++]=paramIOutputFun.tArg;
    /* store ynew in second PostStepFun() argument */
    memcpy(mxGetPr(paramIOutputFun.yArg),ynew,paramGlobal.d*sizeof(double));
    rhs[offset++]=paramIOutputFun.yArg;
    /* store fnew in third PostStepFun() argument */
    memcpy(mxGetPr(paramIOutputFun.fArg),fnew,paramGlobal.d*sizeof(double));
    rhs[offset++]=paramIOutputFun.fArg;
    /* store empty string in fourth PostStepFun() argument */
    rhs[offset++]=mxCreateString("");
    /* call post-step function */
    lhs[0] = NULL;
    lhs[1] = NULL;
    lhs[2] = NULL;
    switch (paramGlobal.funcCallMethod) {
    case 0: 
      mexCallMATLAB(3,lhs,offset,rhs,paramIOutputFun.PostStepFun);
      mxDestroyArray(rhs[offset-1]);
      break;
    case 1:
      mexCallMATLAB(3,lhs,offset,rhs,"feval");
      mxDestroyArray(rhs[offset-1]);
      break;
    default: 
      stopMexFunction(1002,0,0,0,0,0);
      break;
    }
    /* deal with return value(s) */
    if (lhs[0]==NULL) {
      stopMexFunction(1002,0,0,0,0,0);
    } else {
      status=mxGetScalar(lhs[0]);
      mxDestroyArray(lhs[0]);
    }
    if (status==1) { /* Stop, NOW */
      if (lhs[1] != NULL) mxDestroyArray(lhs[1]);
      if (lhs[2] != NULL) mxDestroyArray(lhs[2]);
      *intr=-1; 
      return;
    }
    if ((status!=2) &&(status!=3)) { /* just continue */
      if (lhs[1] != NULL) mxDestroyArray(lhs[1]);
      if (lhs[2] != NULL) mxDestroyArray(lhs[2]);
      *intr=1; 
      return;
    }
    if (status==2) { /* use only ynew return value */
      /* check return value ynew */
      if (lhs[1]==NULL) 
    stopMexFunction(331,0,0,0,0,0);
      if ((!mxIsDouble(lhs[1])) || (mxGetNumberOfDimensions(lhs[1])!=2))
    stopMexFunction(332,0,0,0,0,0);
      if (!(((mxGetM(lhs[1])==paramGlobal.d) && (mxGetN(lhs[1])==1)) ||
        ((mxGetM(lhs[1])==1) && (mxGetN(lhs[1])==paramGlobal.d))))
    stopMexFunction(333,mxGetM(lhs[1]),mxGetN(lhs[1]),0,0,0);
      /* copy back */
      memcpy(ynew,mxGetPr(lhs[1]),paramGlobal.d*sizeof(double));
      /* free memory */
      mxDestroyArray(lhs[1]);  
      if (lhs[2] != NULL) mxDestroyArray(lhs[2]);
      /* call RHSFun() to update fnew */
      ROWMAP_RHSFun (n, tnew, ynew, fnew, rpar, ipar);
      *intr=1; 
      return;
    }
    if (status==3) { /* use ynew and fnew return value */
      /* check return value ynew */
      if (lhs[1]==NULL) 
    stopMexFunction(331,0,0,0,0,0);
      if ((!mxIsDouble(lhs[1])) || (mxGetNumberOfDimensions(lhs[1])!=2))
    stopMexFunction(332,0,0,0,0,0);
      if (!(((mxGetM(lhs[1])==paramGlobal.d) && (mxGetN(lhs[1])==1)) ||
        ((mxGetM(lhs[1])==1) && (mxGetN(lhs[1])==paramGlobal.d))))
    stopMexFunction(333,mxGetM(lhs[1]),mxGetN(lhs[1]),0,0,0);
      /* copy back */
      memcpy(ynew,mxGetPr(lhs[1]),paramGlobal.d*sizeof(double));
      /* check return value fnew */
      if (lhs[2]==NULL) 
    stopMexFunction(341,0,0,0,0,0);
      if ((!mxIsDouble(lhs[2])) || (mxGetNumberOfDimensions(lhs[2])!=2))
    stopMexFunction(342,0,0,0,0,0);
      if (!(((mxGetM(lhs[2])==paramGlobal.d) && (mxGetN(lhs[2])==1)) ||
        ((mxGetM(lhs[2])==1) && (mxGetN(lhs[2])==paramGlobal.d))))
    stopMexFunction(343,mxGetM(lhs[2]),mxGetN(lhs[2]),0,0,0);
      /* copy back */
      memcpy(fnew,mxGetPr(lhs[2]),paramGlobal.d*sizeof(double));

      /* free memory */
      mxDestroyArray(lhs[1]);  
      mxDestroyArray(lhs[2]);  
      *intr=1; 
      return;
    }


    
  } /* end of 'call PostStepFun() depending on paramIOutputFun.PostStepCallMode' */

  *intr=1;
}





static void createTYArrays (mxArray* plhs[]) {
    INTTYPE i,d,count, nout;
    PListElement current;
    double *tPointer, *yPointer, *vPointer;
    
    d=paramGlobal.d;
    nout=paramIOutputFun.tyNoOfElements;
    plhs[0]=mxCreateDoubleMatrix(nout,1,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(nout,d,mxREAL);
    tPointer=mxGetPr(plhs[0]);
    yPointer=mxGetPr(plhs[1]);
    
    /* in the first element of the list nothing is stored!! */
    current=paramIOutputFun.tyList.next;
    count=0;
    while (current!=NULL) {
      vPointer=current->values;
      *tPointer= *vPointer; 
      tPointer++; 
      vPointer++;
      /* the y matrix is stored column-wise as a vector; we fill y row by row */
      for (i=0; i<d; i++,vPointer++) 
    yPointer[count+i*nout]= *vPointer;
      current=current->next;
      count++;
    }
  }
static void createStatVector (mxArray* plhs[]) {
    double* dpointer;
    plhs[2]=mxCreateDoubleMatrix(1,5,mxREAL);
    dpointer=mxGetPr(plhs[2]);
    *dpointer=paramROWMAP.IDID;        dpointer++;
    *dpointer=paramROWMAP.IWORK[5-1];  dpointer++;
    *dpointer=paramROWMAP.IWORK[6-1];  dpointer++;
    *dpointer=paramROWMAP.IWORK[7-1];  dpointer++;
    *dpointer=paramROWMAP.IWORK[8-1];
  }
static void createHPred (mxArray* plhs[]) {  
    plhs[3]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[3]) = paramROWMAP.hs;
  }
  
static void createOutput (INTTYPE nlhs, mxArray* plhs[]) {
  switch (nlhs) {
  case 0:
    break;
  case 2:
    createTYArrays(plhs);
    break;
  case 3:
    createTYArrays(plhs);
    createStatVector(plhs);
    break;
  case 4:
    createTYArrays(plhs);
    createStatVector(plhs);
    createHPred(plhs);
    break;
  default:
    stopMexFunction(1,nlhs,0,0,0,0);
  }
}


void mexFunction (int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[]) {
  initVars();
  checkNumberOfArgs(nlhs,nrhs);
  processArgs(nrhs,prhs);
  extractOptionSettings();
  prepareWorkArrays();
  extractOptions();
  prepareHelpMxArrays();

  /* make init call to output and post-step function
     and store initial data in tyList if required
   */
  callOutputAndPostStepFun_Init();

  /* always call output function from rowmap */
  paramROWMAP.IOUT = 1;
  /* call rowmap time integration */
  paramROWMAP.IDID=1;

  ROWMAP_(&paramGlobal.d,
      &ROWMAP_RHSFun, &paramROWMAP.IFCN,
      &paramROWMAP.tStart, paramROWMAP.yStart, &paramROWMAP.tEnd,
      &paramROWMAP.hs,
      paramROWMAP.RTOL,paramROWMAP.ATOL,&paramROWMAP.ITOL,
      &ROWMAP_JacVFun, &paramROWMAP.IJACV,
      &ROWMAP_FdtFun, &paramROWMAP.IFDT,
      &ROWMAP_OutputFun, &paramROWMAP.IOUT,
      paramROWMAP.WORK,&paramROWMAP.LWORK,
      paramROWMAP.IWORK,&paramROWMAP.LIWORK,
      paramROWMAP.RPAR,paramROWMAP.IPAR,&paramROWMAP.IDID);

  /* make done call to output and post-step function */
  callOutputAndPostStepFun_Done();
  
  switch (paramROWMAP.IDID) {
  case 1:
    break;
  case 2:
    mexPrintf("ROWMAP: computation interrupted by OutputFun() or PostStepFun().\n");
    break;
  case -1:
    mexPrintf("ERROR: ROWMAP terminated because step-size too small.\n");
    break;
  case -2:
    mexPrintf("ERROR: ROWMAP terminated because too many time steps.\n");
    break;
  case -3:
    mexPrintf("ERROR: ROWMAP terminated because input is not consistent.\n"); 
    break;
  case -4:
    mexPrintf("ERROR: ROWMAP terminated because internal Krylov matrix is repeatedly singular.\n");
    break;
  default: 
    stopMexFunction(1101,paramROWMAP.IDID,0,0,0,0);
  }
  
  createOutput(nlhs,plhs);
  doneVars();
}
