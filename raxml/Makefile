# Wrapper to build RAxML versions

SATIVA_HOME=..

# Versions of RAxML to build
VERSIONS:=SSE3.PTHREADS PTHREADS 

# AVX needs gcc >= 4.6.0 
ifeq ("$(USE_GCC_46_OR_HIGHER)", "yes")
  VERSIONS+=AVX.PTHREADS
endif

TARFILE=$(wildcard sativa-RAxML-8.*.*gz)

ifeq (($DARWIN),1)
  MAKE_SUFFIX=mac	
else
  MAKE_SUFFIX=gcc
endif

RAXML_BIN=$(SATIVA_HOME)/raxml/raxmlHPC8-

## targets called by "upsteam makefile":

all: install

clean:
	rm -rf builddir.*
	rm -f $(RAXML_BIN)*
	rm -f unpack.*.stamp

## internal targets

install: build
	clean

build:  $(VERSIONS:%=$(RAXML_BIN)%)

# (MAKEFLAGS need to be unset as ARB sets "-r", suppressing built-in rules,
#  which are needed by RAxML)

$(RAXML_BIN)%: unpack.%.stamp
	MAKEFLAGS= $(MAKE) -C $(<:unpack.%.stamp=builddir.%) -f Makefile.$(@:$(RAXML_BIN)%=%).$(MAKE_SUFFIX)
	cp $(<:unpack.%.stamp=builddir.%)/raxmlHPC-* $@

unpack.%.stamp: $(TARFILE)
	(DIR=$(@:unpack.%.stamp=builddir.%); \
	rm -rf $$DIR &&\
	mkdir $$DIR &&\
	tar -C $$DIR --strip-components=1 -xzf $(TARFILE) &&\
	sed -i 's/-march=native//' $$DIR/Makefile*)
	touch $@ -r $(TARFILE)	

%.stamp:
	echo $@

# prevent make from deleting intermediate targets:
.SECONDARY:

.PHONY: clean all install build
