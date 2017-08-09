/*
 * MATLAB Compiler: 4.2 (R14SP2)
 * Date: Wed Jan 10 16:06:55 2007
 * Arguments: "-B" "macro_default" "-m" "-W" "main" "-T" "link:exe" "mar.m" 
 */

#include <stdio.h>
#include "mclmcr.h"
#ifdef __cplusplus
extern "C" {
#endif
extern const unsigned char __MCC_mar_public_data[];
extern const char *__MCC_mar_name_data;
extern const char *__MCC_mar_root_data;
extern const unsigned char __MCC_mar_session_data[];
extern const char *__MCC_mar_matlabpath_data[];
extern const int __MCC_mar_matlabpath_data_count;
extern const char *__MCC_mar_classpath_data[];
extern const int __MCC_mar_classpath_data_count;
extern const char *__MCC_mar_lib_path_data[];
extern const int __MCC_mar_lib_path_data_count;
extern const char *__MCC_mar_mcr_runtime_options[];
extern const int __MCC_mar_mcr_runtime_option_count;
extern const char *__MCC_mar_mcr_application_options[];
extern const int __MCC_mar_mcr_application_option_count;
#ifdef __cplusplus
}
#endif

static HMCRINSTANCE _mcr_inst = NULL;


static int mclDefaultPrintHandler(const char *s)
{
    return fwrite(s, sizeof(char), strlen(s), stdout);
}

static int mclDefaultErrorHandler(const char *s)
{
    int written = 0, len = 0;
    len = strlen(s);
    written = fwrite(s, sizeof(char), len, stderr);
    if (len > 0 && s[ len-1 ] != '\n')
        written += fwrite("\n", sizeof(char), 1, stderr);
    return written;
}


/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_mar_C_API 
#define LIB_mar_C_API /* No special import/export declaration */
#endif

LIB_mar_C_API 
bool marInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler
)
{
    if (_mcr_inst != NULL)
        return true;
    if (!mclmcrInitialize())
        return false;
    if (!mclInitializeComponentInstance(&_mcr_inst, __MCC_mar_public_data,
                                        __MCC_mar_name_data,
                                        __MCC_mar_root_data,
                                        __MCC_mar_session_data,
                                        __MCC_mar_matlabpath_data,
                                        __MCC_mar_matlabpath_data_count,
                                        __MCC_mar_classpath_data,
                                        __MCC_mar_classpath_data_count,
                                        __MCC_mar_lib_path_data,
                                        __MCC_mar_lib_path_data_count,
                                        __MCC_mar_mcr_runtime_options,
                                        __MCC_mar_mcr_runtime_option_count,
                                        true, NoObjectType, ExeTarget, NULL,
                                        error_handler, print_handler))
        return false;
    return true;
}

LIB_mar_C_API 
bool marInitialize(void)
{
    return marInitializeWithHandlers(mclDefaultErrorHandler,
                                     mclDefaultPrintHandler);
}

LIB_mar_C_API 
void marTerminate(void)
{
    if (_mcr_inst != NULL)
        mclTerminateInstance(&_mcr_inst);
}

int main(int argc, const char **argv)
{
    int _retval;
    if (!mclInitializeApplication(__MCC_mar_mcr_application_options,
                                  __MCC_mar_mcr_application_option_count))
        return 0;
    
    if (!marInitialize())
        return -1;
    _retval = mclMain(_mcr_inst, argc, argv, "mar", 1);
    if (_retval == 0 /* no error */) mclWaitForFiguresToDie(NULL);
    marTerminate();
    mclTerminateApplication();
    return _retval;
}
