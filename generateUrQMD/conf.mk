NAME=urqmd3p9gev
NUM=50
CLONE_LIST=rdlist

TOOLCHAIN=gnu
#TOOLCHAIN=intel

.PHONY: makeseed cls
CLONE_DEP = makeseed

PROJ_INCS='tables.dat'

ifeq ($(TOOLCHAIN), intel)
LINK_CMD=source /softwares/ifort/intel/composer_xe_2013.5.192/bin/compilervars.sh intel64 && ifort
LINK_EXT_FLAG= -lstdc++ -nofor_main
else
LINK_CMD=g++
LINK_EXT_FLAG= -lgfortran
endif

mysrc_lib_file := $(wildcard ./src/mysrc/*.f)

./src/mysrc/libcurqmd.a: $(mysrc_lib_file)
	cd src/mysrc && make

_build/urqmd: urqmd.cc src/TUrQMD.h ./src/mysrc/libcurqmd.a
	@if ! [ -d _build ]; then \
		mkdir _build; \
	fi;
	g++ -O3 -c $(CFLAGS) $(RUN).cc -o _build/$(RUN).o
	$(LINK_CMD) -O3 -o _build/$(RUN) _build/$(RUN).o $(LDFLAGS) -L./src/mysrc -Wl,-Bstatic -lcurqmd -Wl,-Bdynamic $(LIBS) $(LINK_EXT_FLAG)

OTHERCLEAN = mysrc_lib_clean

mysrc_lib_clean:
	cd src/mysrc && make clean

makeseed:
	rm -f rdlist
	cd src && python extrd.py $(NUM) > ../rdlist

cls:
	cd ./src/mysrc && make clean
