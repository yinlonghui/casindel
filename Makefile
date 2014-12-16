CC = gcc
DEBUG = -g  -O2
CFLAGS?= $(DEBUG) -Wall -D_FILE_OFFSET_BITS=64 -pthread  -DDEBUG  
LIBS =  -lm -lz -lpthread
BIN = casindel filter readsPairCheck sv2vcf  distHTfilter

.SUFFIXES: .c .o
.PHONY:clean 

SAMTOOLS_DIR = samtools-1.0
HSTLIB_DIR =  $(SAMTOOLS_DIR)/htslib-1.0
INCLUDE =  -I. -I$(SAMTOOLS_DIR) -I$(HSTLIB_DIR)

%.o: %.c
	$(CC)  -c $(CFLAGS) $(INCLUDE) -o $@  $^


obj = casindel.o  parse.o

filterobj = filter.o 

reads-obj = readsPairCheck.o 

vcf-obj = sv2vcf.o

distHTfilter-obj =  distHTfilter.o


all:$(BIN)

casindel:$(obj)  $(SAMTOOLS_DIR)/libbam.a  $(HSTLIB_DIR)/libhts.a
	$(CC) -o  $@ $^ $(LIBS) $(SAMTOOLS_DIR)/libbam.a $(HSTLIB_DIR)/libhts.a  

filter:$(filterobj)  $(SAMTOOLS_DIR)/libbam.a  $(HSTLIB_DIR)/libhts.a
	$(CC) -o  $@ $^ $(LIBS) $(SAMTOOLS_DIR)/libbam.a $(HSTLIB_DIR)/libhts.a  

readsPairCheck:$(reads-obj)  $(SAMTOOLS_DIR)/libbam.a  $(HSTLIB_DIR)/libhts.a
	$(CC) -o  $@ $^ $(LIBS) $(SAMTOOLS_DIR)/libbam.a $(HSTLIB_DIR)/libhts.a  

sv2vcf:$(vcf-obj)  $(SAMTOOLS_DIR)/libbam.a  $(HSTLIB_DIR)/libhts.a
	$(CC) -o  $@ $^ $(LIBS) $(SAMTOOLS_DIR)/libbam.a $(HSTLIB_DIR)/libhts.a  

distHTfilter:$(distHTfilter-obj)  $(SAMTOOLS_DIR)/libbam.a  $(HSTLIB_DIR)/libhts.a
	$(CC) -o  $@ $^ $(LIBS) $(SAMTOOLS_DIR)/libbam.a $(HSTLIB_DIR)/libhts.a  


clean:
	rm $(BIN) *o  *~
