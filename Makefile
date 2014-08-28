CC = gcc
DEBUG = -g
CFLAGS?= $(DEBUG) -Wall -D_FILE_OFFSET_BITS=64 -pthread  -DDEBUG
LIBS =  -lm -lz -lpthread
BIN = casindel

.SUFFIXES: .c .o
.PHONY:clean 

SAMTOOLS_DIR = samtools-1.0
HSTLIB_DIR =  $(SAMTOOLS_DIR)/htslib-1.0
INCLUDE =  -I. -I$(SAMTOOLS_DIR) -I$(HSTLIB_DIR)

%.o: %.c
	$(CC)  -c $(CFLAGS) $(INCLUDE) -o $@  $^


obj = casindel.o  parse.o


all:$(BIN)

casindel:$(obj)  $(SAMTOOLS_DIR)/libbam.a  $(HSTLIB_DIR)/libhts.a
	$(CC) -o $@ $^ $(LIBS) $(SAMTOOLS_DIR)/libbam.a $(HSTLIB_DIR)/libhts.a


clean:
	rm $(BIN) *o  *~
