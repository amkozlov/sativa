# Wrapper to build RAxML versions
# Based on ARB RAxML8 Makefile, by Elmar Pruesse, http://www.arb-home.de/

SATIVA_HOME=..

# Versions of RAxML to build
VERSIONS:=SSE3.PTHREADS PTHREADS 

ifeq ("$(USE_AVX)", "yes")
  VERSIONS+=AVX.PTHREADS
endif

ifeq ("$(USE_AVX2)", "yes")
  VERSIONS+=AVX2.PTHREADS
endif

TARFILE=$(wildcard sativa-RAxML-8.*.*gz)

ifeq ("$(DARWIN)","1")
  MAKE_SUFFIX=mac
  SED=sed -i ''
else
  MAKE_SUFFIX=gcc
  SED=sed -i
endif

RAXML_BIN=$(SATIVA_HOME)/raxml/raxmlHPC8-

## targets called by "upsteam makefile":

all: install

clean_builddir:
	rm -rf builddir.*
	rm -f unpack.*.stamp

clean: clean_builddir
	rm -f $(RAXML_BIN)*

## internal targets

install: build clean_builddir

build:  $(VERSIONS:%=$(RAXML_BIN)%)

# (MAKEFLAGS need to be unset as ARB sets "-r", suppressing built-in rules,
#  which are needed by RAxML)

$(RAXML_BIN)%: unpack.%.stamp
	MAKEFLAGS= $(MAKE) -C $(<:unpack.%.stamp=builddir.%) -f Makefile.$(@:$(RAXML_BIN)%=%).$(MAKE_SUFFIX) -j 2 CC=$(CC)
	cp $(<:unpack.%.stamp=builddir.%)/raxmlHPC-* $@

unpack.%.stamp: $(TARFILE)
	(DIR=$(@:unpack.%.stamp=builddir.%); \
	rm -rf $$DIR &&\
	mkdir $$DIR &&\
	tar -C $$DIR --strip-components=1 -xzf $(TARFILE) &&\
	$(SED) 's/-march=native//' $$DIR/Makefile*)
	touch $@ -r $(TARFILE)	

%.stamp:
	echo $@

# prevent make from deleting intermediate targets:
.SECONDARY:

.PHONY: clean all install build clean_builddir
