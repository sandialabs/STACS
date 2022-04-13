CHARMC  = charmc $(OPTS)
OUT    := stacs

SRC_CI := $(wildcard *.ci)
SRC_C  := $(wildcard *.C) \
          $(wildcard models/*.C)
INC    := -I./ -I./models

DECL   := $(SRC_CI:.ci=.decl.h)
DEF    := $(SRC_CI:.ci=.def.h)
OBJ    := $(SRC_C:.C=.o)

CFLAGS     = -O3 -Wall
CHARMFLAGS = -module CkMulticast -language charm++
PROJFLAGS  = -tracemode projections -tracemode summary
YARPFLAGS  = -DSTACS_WITH_YARP

LIB        = -std=c++11
LDLIB      = -lm -lyaml-cpp -lfftw3
YARPLDLIB  = -lYARP_init -lyarpcar -lYARP_OS \
             -lYARP_sig -lYARP_math -lYARP_dev \
             -lYARP_name -lACE

.PHONY: all projections clean

all: $(OUT)

$(OUT): $(DECL) $(DEF) $(OBJ)
	$(CHARMC) $(CFLAGS) $(CHARMFLAGS) $(OBJ) $(LDLIB) -o $(OUT)

projections: $(DECL) $(DEF) $(OBJ)
	$(CHARMC) $(CFLAGS) $(PROJFLAGS) $(CHARMFLAGS) $(OBJ) $(LDLIB) -o $(OUT)

%.o: %.C
	$(CHARMC) $(CFLAGS) $(CHARMFLAGS) $(INC) $(LIB) $< -o $@

$(DECL) $(DEF): $(SRC_CI)
	$(CHARMC) $(CFLAGS) $(SRC_CI)

clean:
	rm -f  $(DECL) $(DEF) $(OBJ) $(OUT) charmrun

# YARP related
ifeq ($(WITHYARP), 1)
  CFLAGS += $(YARPFLAGS)
  LDLIB += $(YARPLDLIB)
endif

# OS related
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LIB += -stdlib=libc++ -mmacosx-version-min=10.15
endif
