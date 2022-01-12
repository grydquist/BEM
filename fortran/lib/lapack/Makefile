FC	   = ifort
FFLAGS = -O3 -recursive
OBJDIR = obj/
SRCDIR = src/
SRCS   = $(wildcard $(SRCDIR)*.f)
OBJS   = $(patsubst $(SRCDIR)%.f,$(OBJDIR)%.o,$(SRCS))
LN     = liblapack.a
AR 	   = ar rv
RANLIB = ar -ts

$(OBJDIR)%.o : $(SRCDIR)%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -o $@ -c $<

$(LN)        : $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

$(OBJS)      : | $(OBJDIR) $

$(OBJDIR)  :
	mkdir -p $@

.PHONY : ech clean
ech:
	@echo $(OBJS)

clean:
	rm -f obj/*
	rm -f $(LN)
