REXEC = R
export REXEC

#BUILDARGS = --force-biarch --dsym
#BUILDARGS = --dsym

TESTFILE = inst/tools/testModels.R

# subdirectories
RSOURCE = R
RDOCUMENTS = man
RDATA = data

# file types
RFILES = $(wildcard R/*.R)

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo ""	
	@echo "INSTALL"
	@echo ""	
	@echo "  install       install dynr"
	@echo ""	
	@echo "CLEANING"
	@echo ""	
	@echo "  clean      remove all files from the build directory"
	@echo "TESTING"
	@echo ""	
	@echo "  test               run the test suite"
	@echo "  torture       run the test suite with gctorture(TRUE)"


r-libs-user-dir:
	./inst/tools/mk-r-libs-user-dir

install:
	MAKEFLAGS="$(INSTALLMAKEFLAGS)" $(REXEC) CMD INSTALL $(BUILDARGS) .

#MAKEFLAGS="$(INSTALLMAKEFLAGS)" $(REXEC) CMD SHLIB $(BUILDARGS) src/wrappernegloglike.c src/estimation_nloptR.c src/PANAmodel.c src/functions/*.c

clean:
	mkdir -p build
	-rm build/dynr_*.tar.gz
	-rm src/*.o
	-rm src/*.so
	-rm src/*.dll
	-rm src/Makevars
	-rm config.log config.status

test:
	$(REXEC) --vanilla --slave -f $(TESTFILE)

torture:
	$(REXEC) --vanilla --slave -f $(TESTFILE) --args gctorture

