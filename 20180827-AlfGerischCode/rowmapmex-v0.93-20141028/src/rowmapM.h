#ifndef rowmaph
#define rowmaph

/* (A) Option settings */
#define OPT_WARNMISS "OptWarnMiss" /* extracted in: extractOptionSettings() */
#define OPT_WARNTYPE "OptWarnType" /* extracted in: extractOptionSettings() */
#define OPT_WARNSIZE "OptWarnSize" /* extracted in: extractOptionSettings() */

/* (B) Initial step size and step-size selection parameters */
#define OPT_RTOL      "RelTol"      /* extracted in: extractTOLs() */
#define OPT_ATOL      "AbsTol"      /* extracted in: extractTOLs() */
#define OPT_INITIALSS "InitialStep" /* extracted in: extractHs() */
#define OPT_SSMINSEL  "StepSizeMinSelection" /* extracted in: extractWORKOpt() */
#define OPT_SSMAXSEL  "StepSizeMaxSelection" /* extracted in: extractWORKOpt() */
#define OPT_RHO       "StepSizeSafetyFactor" /* extracted in: extractWORKOpt() */
#define OPT_MAXSTP    "MaxStepSize"          /* extracted in: extractWORKOpt() */
#define OPT_MAXSTEPS  "MaxNumberOfSteps"     /* extracted in: extractIWORKOpt() */

/* (C) Optional functions and call modes*/
#define OPT_JACVFUNCTION     "JacVFun"          /* extracted in: extractOptFuns() */
#define OPT_FDTFUNCTION      "FdtFcn"           /* extracted in: extractOptFuns() */
#define OPT_OUTPUTFUNCTION   "OutputFun"        /* extracted in: extractOptFuns() */
#define OPT_POSTSTEPFUNCTION "PostStepFun"      /* extracted in: extractOptFuns() */
#define OPT_RETURNMODE       "ReturnMode"       /* extracted in: extractOptFuns() */
#define OPT_FUNCCALLMETHOD   "FuncCallMethod"   /* extracted in: extractOptFuns() */
#define OPT_OUTPUTCALLMODE   "OutputCallMode"   /* extracted in: extractOptFuns() */
#define OPT_POSTSTEPCALLMODE "PostStepCallMode" /* extracted in: extractOptFuns() */

/* (D) Krylov process */
#define OPT_MAXKRYLOV "MaxKrylovSteps" /* extracted in: prepareWorkArrays() */
#define OPT_KTOL      "KTol"           /* extracted in: extractWORKOpt() */

/* (E) Other options */
#define OPT_ROWMAPCOEFFS "ROWMAPCoefficientSet" /* extracted in: extractIWORKOpt() */
#define OPT_RHSNONAUTONOMOUS "RHSnonautonomous" /* extracted in: extractIfcn() */
#define OPT_CONTOUTPUTINTERP "ContOutputInterp" /* extracted in: extractContOutputInterp()  */
#define OPT_EPS              "eps"              /* extracted in: extractWORKOpt() */

struct ListElement
{ /* linked list with (t,y) pairs; all vectors values have length d+1 */
  double* values;
  struct ListElement *next;
};
typedef struct ListElement SListElement;
typedef SListElement* PListElement;

struct ParameterGlobal 
{ /* global variables  */
  INTTYPE d;              /* dimension of the ODE system */
  INTTYPE tspanLength;    /* length of tspan */  
  double* tspanPointer;   /* pointer into tspan */
  double direction;   /* sign(tEnd-tStart), currently only = 1.0 supported */  
  INTTYPE funcCallMethod; /* see OPTION  FuncCallMethod (extractOptFuns) */ 
};
typedef struct ParameterGlobal SParameterGlobal;

struct ParameterOptions
{ /* parameters for user supplied options */    
  const mxArray *opt;  /* options array */
  INTTYPE optCreated;      /* flag==1 if options array is created because it was 
			  not supplied with the rowmap call; ==0 otherwise */
};
typedef struct ParameterOptions SParameterOptions;

struct ParameterROWMAP
{ /* parameter for call to FORTRAN ROWMAP */
  double tStart;    /* initial time (processArgs) */
  double tEnd;      /* final time (processArgs) */
  double hs;        /* initial step size (extractHs) */
  double *yStart;   /* initial values (processArgs) */
  double *RTOL;     /* relative tolerance (extractTOLs) */
  double *ATOL;     /* absolute tolerance (extractTOLs) */
  INTTYPE ITOL;         /* switch for RTOL und ATOL (extractTOLs) */
  INTTYPE IFCN;         /* switch for autonomous/nonautonomous (extractIfcn) */
  INTTYPE IJACV;        /* switch for Jacobian x vector function (extractOptFuns) */
  INTTYPE IFDT;         /* switch for FDT function (extractOptFuns) */
  INTTYPE IOUT;         /* switch for Output function (always set to 1) */
  double *WORK;     /* double work array (prepareWorkArrays) */
  INTTYPE LWORK;        /* length of WORK (prepareWorkArrays) */
  INTTYPE *IWORK;       /* integer work array (prepareWorkArrays) */
  INTTYPE LIWORK;       /* length of IWORK (prepareWorkArrays) */
  double *RPAR;     /* additional RPAR double array */
  INTTYPE *IPAR;        /* additional IPAR integer array */
  INTTYPE IDID;         /* return status */
};
typedef struct ParameterROWMAP SParameterROWMAP;

struct ParameterRHSFun
{ /* parameter for call to right-hand side function */
  char *RHSFun;                  /* function name of RHSFun() if provided as 
				    string and FuncCallMethod==1 */
  const mxArray *RHSFunH;        /* function handle or inline function for 
				    RHSFun() */
  mxArray *tArg;                 /* To call RHSFun(): t */
  mxArray *yArg;                 /* To call RHSFun(): y */
};
typedef struct ParameterRHSFun SParameterRHSFun;

struct ParameterJacVFun
{ /* parameter for Jacobian x vector function */
  char *JacVFun;                 /* function name of JacVFun() if provided as 
				    string and FuncCallMethod==1 */  
                                    
  const mxArray *JacVFunH;       /* function handle or inline function for 
				    JacVFun() */
  mxArray *tArg;                 /* To call JacVFun(): t */
  mxArray *yArg;                 /* To call JacVFun(): y */
  mxArray *vArg;                 /* To call JacVFun(): v */
};
typedef struct ParameterJacVFun SParameterJacVFun;


struct ParameterFdtFun
{ /* parameter for call to time derivative function */
  char *FdtFun;                  /* function name of FdtFun() if provided as 
				    string and FuncCallMethod==1 */
  const mxArray *FdtFunH;        /* function handle or inline function for 
				    FdtFun() */
  mxArray *tArg;                 /* To call FdtFun(): t */
  mxArray *yArg;                 /* To call FdtFun(): y */
};
typedef struct ParameterFdtFun SParameterFdtFun;




struct ParameterIOutputFun
{ /* parameter for the internal output function */
  SListElement tyList;        /* start of (t,y) linked list*/
  PListElement tyLastElement; /* pointer to last entry in (t,y) linked list */
  INTTYPE tyNoOfElements;         /* number of elements in (t,y) linked list */
  INTTYPE tspanPos;               /* position of the next output time point in tspan vector */

  INTTYPE nlhs;                   /* number of return arguments of rowmap call */
  INTTYPE ReturnMode;             /* see OPTION ReturnMode */
  INTTYPE OutputCallMode;         /* see OPTION OutputCallMode */
  INTTYPE PostStepCallMode;       /* see OPTION PostStepCallMode */
  INTTYPE ContOutputInterp;       /* see OPTION ContOutputInterp (extractContOutputInterp) */

  char *OutputFun;            /* function name of OutputFun() if provided as 
				 string and FuncCallMethod==1 */
  const mxArray *OutputFunH;  /* function handle or inline function for 
				 OutputFun() */
  char *PostStepFun;          /* function name of PostStepFun() if provided 
				 as string and FuncCallMethod==1 */
  const mxArray *PostStepFunH;/* function handle or inline function for 
				 PostStepFun() */
 
  mxArray *tArg;              /* To call OutputFun()/PostStepFun(): t */
  mxArray *yArg;              /* To call OutputFun()/PostStepFun(): y */  
  mxArray *fArg;              /* To call PostStepFun(): f(t,y) */  
};
typedef struct ParameterIOutputFun SParameterIOutputFun;



typedef void (*T_ROWMAP_RHSFun)(INTTYPE *n, double *t,
				double *x, double *f, 
				double *rpar, INTTYPE *ipar);

typedef void (*T_ROWMAP_JacVFun)(INTTYPE *n, double *t,
				 double *y, double *v, double *z, 
				 double *rpar, INTTYPE *ipar);

typedef void (*T_ROWMAP_FdtFun)(INTTYPE *n, double *t,
				double *y, double *ft,
				double *rpar, INTTYPE *ipar);

typedef void (*T_ROWMAP_OutputFun)(INTTYPE *n, 
				   double *told, double *tnew, 
				   double *yold, double *ynew,
				   double *fold, double *fnew, 
				   double *ycon, INTTYPE* intr,
				   double *rpar, INTTYPE *ipar);


#ifdef FORTRANNOUNDER
/* Fortran functions without underscore */
#ifdef FORTRANUPP
/* Fortran functions without underscore  & UPPERCASE letters */
#define ROWMAP_ ROWMAP
#define ROWCON_ ROWCON
#define ROWLIN_ ROWLIN
#else
/* Fortran functions without underscore  & lowercase letters */
#define ROWMAP_ rowmap
#define ROWCON_ rowcon
#define ROWLIN_ rowlin
#endif
#else
/* Fortran functions with underscore */
#ifdef FORTRANUPP
/* Fortran functions with underscore & UPPERCASE letters */
#else
/* Fortran functions with underscore & lowercase letters */
#define ROWMAP_ rowmap_
#define ROWCON_ rowcon_
#define ROWLIN_ rowlin_
#endif
#endif


extern void ROWMAP_ (INTTYPE *n,
		     T_ROWMAP_RHSFun fcn, INTTYPE *ifcn,
		     double *tStart, double *x, double *tEnd,
		     double *hs,
		     double  *rtol, double *atol, INTTYPE *itol,
		     T_ROWMAP_JacVFun jacv, INTTYPE *ijacv,
		     T_ROWMAP_FdtFun fdt, INTTYPE *ifdt,
		     T_ROWMAP_OutputFun solout, INTTYPE *iout,
		     double *work, INTTYPE *lwork, INTTYPE *iwork, INTTYPE *liwork,
		     double *rpar, INTTYPE *ipar, INTTYPE *idid);
        
extern double ROWCON_ (INTTYPE *n, double *s, 
		       double *told, double *tnew,
		       double *xold, double *xnew,
		       double *fold, double *fnew, 
		       double *xcon);

extern double ROWLIN_ (INTTYPE *n, double *s, 
		       double *told, double *tnew,
		       double *xold, double *xnew,
		       double *xcon);


#endif
