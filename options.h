
#ifndef _OPTIONS_H
#define _OPTIONS_H
typedef enum { OPTION_SCALAR, OPTION_INTEGER, OPTION_III, OPTION_SSS,OPTION_SS,OPTION_II,OPTION_ISPAIR,OPTION_IIPAIR,OPTION_SPECIAL1,OPTION_PSO,OPTION_GRIDTYPE,OPTION_MULTISELECT} option_type;

PetscErrorCode csvOptions(Options *, Materials *);
PetscErrorCode declare_option( const char *pattern, option_type opt, void *option_ptr, const char *default_value );
PetscErrorCode parse_option( const char *key, const char *value , int );
void print_options( FILE *fp );
void parse_multiselect(MultiSelect *, const char *, int);
int multiselect_compare(MultiSelect *ms, char *value2 );

#endif
