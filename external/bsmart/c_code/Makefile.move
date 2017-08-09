opssmov : opssmov.o MARfit.o EEGrealA.o MARgetpfpb.o matmisc.o matmisc2.o matmisc5.o EEG_tpool.o
	cc -lm -o opssmov opssmov.o MARfit.o EEGrealA.o MARgetpfpb.o matmisc.o matmisc2.o matmisc5.o EEG_tpool.o
opssmov.o : opssmov.c EEGdef.h EEGmat.h
	cc -c opssmov.c
MARfit.o : MARfit.c EEGdef.h EEGmat.h
	cc -c MARfit.c
EEGrealA.o : EEGrealA.c EEGdef.h EEGmat.h
	cc -c EEGrealA.c
MARgetpfpb.o : MARgetpfpb.c EEGdef.h EEGmat.h
	cc -c MARgetpfpb.c
matmisc.o : matmisc.c
	cc -c matmisc.c
matmisc2.o : matmisc2.c
	cc -c matmisc2.c
matmisc5.o : matmisc5.c
	cc -c matmisc5.c
EEG_tpool.o : EEG_tpool.c
	cc -c EEG_tpool.c

