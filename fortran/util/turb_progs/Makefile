# $Id: Makefile,v 1.9 2009-12-01 19:23:49 eric Exp $

#	Copyright 2011 Johns Hopkins University
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

OSARCH := $(shell uname -sp)

ifeq ($(OSARCH),Darwin i386)
# Compile both 32- and 64-bit code under MacOS X for Intel
# ARCH_FLAGS = -arch i386 -arch x86_64
else
	ARCH_FLAGS =
endif

RM     = rm -f
CFLAGS = -Wall
LDLIBS = -lm
CP     = cp
MKDIR  = mkdir -p

CC     = icx -O0$(ARCH_FLAGS)
FC     = ifort $(ARCH_FLAGS)

OBJ =	soapC.o \
	soapClient.o \
	stdsoap2.o \
        turblib.o

all: turb_VG

turb_VG : $(OBJ) turb_VG.o
	 $(FC) -o $@ $(OBJ) turb_VG.o $(LDLIBS)

turb_VG.o : turb_VG.f90
	 $(FC) -c turb_VG.f90

stdsoap2.o: stdsoap2.c
	$(CC) $(CFLAGS) -c $<

static_lib: $(OBJ)
	ar rcs libJHTDB.a $(OBJ)

install: static_lib
	$(MKDIR) $(JHTDB_PREFIX)/include
	$(MKDIR) $(JHTDB_PREFIX)/lib
	$(CP) *.h $(JHTDB_PREFIX)/include/
	$(CP) libJHTDB.a $(JHTDB_PREFIX)/lib/

# Regenerate the gSOAP interfaces if required
TurbulenceService.h : wsdl

# Update the WSDL and gSOAP interfaces
wsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

testwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://test.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

mhdtestwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://mhdtest.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

devwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://dev.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

mhddevwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://mhddev.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

prodtestwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://prodtest.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

clean:
	$(RM) *.o *.exe turb_VG

spotless: clean
	$(RM) soapClient.c TurbulenceServiceSoap.nsmap soapH.h TurbulenceServiceSoap12.nsmap soapStub.h soapC.c TurbulenceService.h

.SUFFIXES: .o .c .x

.c.o:
	$(CC) $(CFLAGS) -c $<

.PHONY: force
compiler_flags: force
	echo '$(CFLAGS)' | cmp -s - $@ || echo '$(CFLAGS)' > $@

$(OBJ): compiler_flags

