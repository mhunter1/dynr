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
	@echo ""
	@echo "BUILDS"
	@echo ""
	@echo "  build         create a dynr binary for the local system"
	@echo "  srcbuild      create a dynr source release"



r-libs-user-dir:
	./inst/tools/mk-r-libs-user-dir

build-prep:
	@if [ $$(git status --short --untracked-files=no 2> /dev/null | wc -l) != 0 ]; then \
	  echo '***'; echo "*** UNCOMMITTED CHANGES IGNORED ***"; \
	  echo '***'; echo "*** Use 'git diff' to see what is uncommitted"; \
          echo '***'; fi
	-[ -d build ] && rm -r ./build
	mkdir build
	git archive --format=tar HEAD | (cd build; tar -xf -)

build: build-prep
	cd build && ./util/prep && $(REXEC) CMD INSTALL $(BUILDARGS) --build .
	egrep -v '@[A-Z]+@' DESCRIPTION.in > DESCRIPTION

srcbuild: build-prep
	cd build && ./util/prep npsol && $(REXEC) CMD build .
	egrep -v '@[A-Z]+@' DESCRIPTION.in > DESCRIPTION
	@echo 'To generate a PACKAGES file, use:'
	@echo '  echo "library(tools); write_PACKAGES('"'.', type='source'"')" | R --vanilla'

install:
	./util/prep
	MAKEFLAGS="$(INSTALLMAKEFLAGS)" $(REXEC) CMD INSTALL $(BUILDARGS) .
	egrep -v '@[A-Z]+@' DESCRIPTION.in > DESCRIPTION

#MAKEFLAGS="$(INSTALLMAKEFLAGS)" $(REXEC) CMD SHLIB $(BUILDARGS) src/wrappernegloglike.c src/estimation_nloptR.c src/PANAmodel.c src/functions/*.c

clean:
	mkdir -p build
	-rm build/dynr_*.tar.gz
	-rm build/dynr_*.zip
	-rm src/*.o
	-rm src/*.so
	-rm src/*.dll
	-rm src-i386/*.o
	-rm src-i386/*.so
	-rm src-i386/*.dll
	-rm src-x64/*.o
	-rm src-x64/*.so
	-rm src-x64/*.dll
	-rm inst/models/passing/*.o
	-rm inst/models/passing/*.so
	-rm inst/models/passing/*.dll
	-rm src/Makevars
	-rm config.log config.status
	egrep -v '@[A-Z]+@' DESCRIPTION.in > DESCRIPTION

test:
	$(REXEC) --vanilla --slave -f $(TESTFILE)

torture:
	$(REXEC) --vanilla --slave -f $(TESTFILE) --args gctorture

