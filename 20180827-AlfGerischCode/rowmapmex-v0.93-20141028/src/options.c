#include <string.h>
#include "mex.h"
#include "inttype.h"
#include "options.h"

INTTYPE opt_getSizeOfOptField (const mxArray *opt, const char *name, INTTYPE *m, INTTYPE *n)
{
  mxArray *field;
  
  field=mxGetField(opt,0,name);
  if (field==NULL) return 1;
  
  if (mxGetNumberOfDimensions(field)!=2) return 2;  
  
  if (mxIsEmpty(field)) return 3;
  
  *m=mxGetM(field);*n=mxGetN(field);
  return 0;
}

double opt_getDoubleFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, double defaultValue)
{
  mxArray *field;
  
  field=mxGetField(opt,0,name);
  if ((field==NULL) ||  (mxIsEmpty(field)))
  {
    if (optSet->warnMiss) 
      mexPrintf("Option '%s' not provided. Using default %e\n",name,defaultValue);    
    return defaultValue;
  }    
  if (!mxIsDouble(field))
  {
    if (optSet->warnType)
      mexPrintf("Option '%s' has wrong type. Using default %e\n",name,defaultValue);
    return defaultValue;
  }
  if ((mxGetNumberOfDimensions(field)!=2) ||
      (mxGetM(field)!=1) || (mxGetN(field)!=1))
  {
    if (optSet->warnSize)
      mexPrintf("Option '%s' has wrong size. Using default %e\n",name,defaultValue);
    return defaultValue;
  }
  return mxGetScalar(field);
}

INTTYPE opt_getDoubleVectorFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, INTTYPE m, INTTYPE n, double *dpointer)
{
  INTTYPE i,l;
  mxArray *field;
  double *dataPointer;
  
  mxAssert(m==1 || n==1,"getDoubleVectorFromOpt: m!=1 && n!=1");

  field=mxGetField(opt,0,name);
  if ((field==NULL) ||  (mxIsEmpty(field)))
  {
    if (optSet->warnMiss) mexPrintf("Option '%s' not provided.\n",name);    
    return 1;
  }    
  if (!mxIsDouble(field))
  {
    if (optSet->warnType) mexPrintf("Option '%s' has wrong type.\n",name);
    return 2;
  }
  if ((mxGetM(field)!=m) || (mxGetN(field)!=n))
  {
    if (optSet->warnSize) mexPrintf("Option '%s' has wrong size.\n");
    return 3;
  }
  if (m>n) l=m; else l=n;
  
  dataPointer=mxGetPr(field);
  for (i=0; i<l; i++, dpointer++,dataPointer++)
    *dpointer= *dataPointer;
  
  return 0;
}

INTTYPE opt_getIntFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, INTTYPE defaultValue)
{
  mxArray *field;
  mxClassID classID;
  
  field=mxGetField(opt,0,name);
  if ((field==NULL) ||  (mxIsEmpty(field)))
  {
    if (optSet->warnMiss) 
      mexPrintf("Option '%s' not provided. Using default %i\n",name,defaultValue);    
    return defaultValue;
  }    
  
  classID=mxGetClassID(field);
  switch (classID)
  {
    case mxDOUBLE_CLASS: case mxINT8_CLASS:   case mxUINT8_CLASS:
    case mxINT16_CLASS:  case mxUINT16_CLASS: case mxINT32_CLASS:
    case mxUINT32_CLASS: case mxINT64_CLASS:  case mxUINT64_CLASS:
    break;
    default:
    if (optSet->warnType)
      mexPrintf("Option '%s' has wrong type. Using default %i\n",
        name,defaultValue);
    return defaultValue;    
  }
    
  if ((mxGetNumberOfDimensions(field)!=2) ||
     (mxGetM(field)!=1) || (mxGetN(field)!=1))
  {
    if (optSet->warnSize)
      mexPrintf("Option '%s' has wrong size. Using default %i\n",
        name,defaultValue);
    return defaultValue;
  } 
  
  switch (classID)
  {
    case mxDOUBLE_CLASS: return (INTTYPE)*((double*)mxGetData(field));break;
    case mxINT8_CLASS:   return (INTTYPE)*((char*)mxGetData(field));break;
    case mxUINT8_CLASS:  return (INTTYPE)*((unsigned char*)mxGetData(field));break;
    case mxINT16_CLASS:  return (INTTYPE)*((short*)mxGetData(field));break;
    case mxUINT16_CLASS: return (INTTYPE)*((unsigned short*)mxGetData(field));break;
    case mxINT32_CLASS:  return (INTTYPE)*((int*)mxGetData(field));break;
    case mxUINT32_CLASS: return (INTTYPE)*((unsigned int*)mxGetData(field));break;
    case mxINT64_CLASS:  return (INTTYPE)*((long long*)mxGetData(field));break;
    case mxUINT64_CLASS: return (INTTYPE)*((unsigned long long*)mxGetData(field));break;
  }

  mxAssert(0,"getIntFromOpt: should not happen");
  return 0;
}

char* opt_getStringFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, char* defaultValue)
{
  static char* noString = "<NULL>";
  
  mxArray *field;
  INTTYPE buflen,takeDefault;
  char* erg;

  takeDefault=0;
  field=mxGetField(opt,0,name);
  if ((field==NULL) ||  (mxIsEmpty(field)))
  {
    if (optSet->warnMiss) 
      mexPrintf("Option '%s' not provided. Using default '%s'\n",name,
        defaultValue==NULL?noString:defaultValue);    
    takeDefault=1;
  } 
  
  if ((!takeDefault) && (!mxIsChar(field)))
  {
    if (optSet->warnType)
      mexPrintf("Option '%s' has wrong type. Using default %s\n",
        name,defaultValue==NULL?noString:defaultValue);
    takeDefault=1;
  }
  
  if (!takeDefault)
  {
    if ((mxGetNumberOfDimensions(field)!=2) || (mxGetM(field)!=1))
    {
      if (optSet->warnSize)
        mexPrintf("Option '%s' has wrong size. Using default %s\n",
          name,defaultValue==NULL?noString:defaultValue);
      takeDefault=1;
    }
  } 
  
  if (takeDefault)
  {
    if (defaultValue==NULL) return NULL;
    buflen=strlen(defaultValue)+1;
    erg=mxMalloc(buflen);
    strcpy(erg,defaultValue);
  } else
  {
    buflen=mxGetN(field)*sizeof(mxChar)+1;
    erg=mxMalloc(buflen);
    mxGetString(field,erg,buflen);
  }
  return erg;
}
