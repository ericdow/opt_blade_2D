FC	= gfortran
CC	= gcc
LD	= $(FC)
F77	= $(FC)
SWP	= swplist
RM	= /bin/rm -f
PROC	= ip27
OLEVEL	= -O0
FOPTS   =  
FFLAGS	= $(OLEVEL) $(FOPTS) -I/usr/include
CFLAGS	= $(COPTS)
LIBS	= 
LDFLAGS	= -lcgns $(FFLAGS)

FOBJS   = mesh_var.o read_cgns_2D.o write_cgns_2D.o mod_coord_2D.o sort3.o

OBJS	= $(FOBJS)

DEFAULT:	perturb_mesh_2D

perturb_mesh_2D:	$(OBJS) perturb_mesh_2D.o
	$(LD) -o $@ $(LDFLAGS) perturb_mesh_2D.o $(OBJS) $(LIBS)

clean:
	$(RM) $(EXEC) $(OBJS) perturb_mesh_2D perturb_mesh_2D.o mesh_var.mod

.SUFFIXES:
.SUFFIXES: .o .F .c .f .swp

.F.o:
	$(FC)  -c $(FFLAGS) $(QV_OPT) $(DEFINES) $<

.f.o:
	$(FC)  -c $(FFLAGS) $(QV_OPT) $(DEFINES) $<

.c.o:
	$(CC)  -c $(CFLAGS) $(QV_OPT) $(DEFINES) $<

.F.swp:
	$(SWP) -c $(FFLAGS) $(QV_OPT) $(DEFINES) -WK,-cmp=$*.m $<

.f.swp:
	$(SWP) -c $(FFLAGS) $(QV_OPT) $(DEFINES) -WK,-cmp=$*.m $<

.c.swp:
	$(SWP) -c $(CFLAGS) $(QV_OPT) $(DEFINES) $<

