#ifndef optionsh
#define optionsh

/* PREFIX: opt */

struct OptionSettings 
{ /* moegliche Einstellungen (Warnstufe, usw.) bei Optionen */
  INTTYPE warnMiss;  /* Warnung, wenn Optionseintrag fehlt */
  INTTYPE warnType;  /* Warnung, wenn Optionseintrag falschen Typ */
  INTTYPE warnSize;  /* Warnung, wenn Optionseintrag falsche Groesse */
};
typedef struct OptionSettings SOptionSettings;
typedef SOptionSettings* POptionSettings;

INTTYPE opt_getSizeOfOptField (const mxArray *opt, const char *name, INTTYPE *m, INTTYPE *n);

double opt_getDoubleFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, double defaultValue);
  
INTTYPE opt_getDoubleVectorFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, INTTYPE m, INTTYPE n, double *dpointer);
  
INTTYPE opt_getIntFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, INTTYPE defaultValue);

char* opt_getStringFromOpt (const mxArray *opt, POptionSettings optSet,
  const char *name, char* defaultValue);


#endif
