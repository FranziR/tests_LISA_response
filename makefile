CC=g++

WarningFLAGS=-Wall -g

IDIR=/usr/local/include/
LDIR=/usr/local/lib/
OBJDIR=$(PWD)/src

LIBS=-lgsl -lm -lgslcblas

OBJ=ldc_code
EXE=ldc_exe
VPATH=src

all: $(OBJDIR)/$(OBJ).o $(EXE)

.PHONY : all

$(OBJDIR)/$(OBJ).o: $(OBJ).c
	$(CC) $(WarningFLAGS) -I$(IDIR) -c $< -o $@


$(EXE): $(OBJDIR)/$(OBJ).o
	$(CC) -L$(LDIR) $< -o $@ $(LIBS)


.PHONY : clean
clean : 
	-rm $(OBJDIR)/$(OBJ).o $(EXE)

