FC	   = mpiifort
CC	   = mpiicc
FFLAGS = -O3 -recursive
OBJDIR = obj/
SRCDIR = src/
INCDIR = inc/
INCDIR2= inc/blank/
CFLAGS = -O3 -DAdd_ -I$(INCDIR) -I$(INCDIR2)
SRCS   = $(wildcard $(SRCDIR)*.f)
SRCSC  = $(wildcard $(SRCDIR)*.c)
OBJS   = $(patsubst $(SRCDIR)%.f,$(OBJDIR)%.o,$(SRCS))
OBJSC  = $(patsubst $(SRCDIR)%.c,$(OBJDIR)%.o,$(SRCSC))
LN     = libscalapack.a
AR 	   = ar rv
RANLIB = ar -ts

$(OBJDIR)%.o : $(SRCDIR)%.c | $(OBJDIR)
	$(CC) $(CFLAGS) -o $@ -c $<

$(OBJDIR)%.o : $(SRCDIR)%.f | $(OBJDIR)
	$(FC) $(FFLAGS) -o $@ -c $<

$(LN)        : $(OBJS) $(OBJSC)
	$(AR) $@ $(OBJS) $(OBJSC)
	$(RANLIB) $@

$(OBJSC)     : | $(OBJDIR) $

$(OBJS)      : | $(OBJDIR) $

$(OBJDIR)  :
	mkdir -p $@

.PHONY : ech clean
ech:
	@echo $(OBJS)

clean:
	rm -f $(OBJDIR)*
	rm -f $(LN)
	rm -f $(SRCDIR)*.o
