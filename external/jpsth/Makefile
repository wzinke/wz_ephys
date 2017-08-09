# COMSPEC is defined only on Windows. 
ifdef COMSPEC
	PATH_SEP := ;
	EXE_EXT := .exe
	MATLAB_ARGS := -nosplash -nojvm -automation
	RM := del
else
	PATH_SEP := :
	EXE_EXT :=
	MATLAB_ARGS := -nosplash -nojvm
	RM := rm -rf
endif
THE_DATE = $(shell date +%Y%m%d)
SOURCES := equation3.c covariogramBrody.c
MATLAB = matlab$(EXE_EXT)
TAR = tar cvf
TARBALL_NAME = jpsth-$(THE_DATE).tar.gz
.SILENT: compile-mex
compile-mex: $(SOURCES)
	$(MATLAB) $(MATLAB_ARGS) -r "mex('equation3.c');mex('covariogramBrody.c');exit"

.SILENT: tarball-release-code
tarball-release-code:
	$(TAR) $(TARBALL_NAME) *

all: compile-mex tarball-release-code



