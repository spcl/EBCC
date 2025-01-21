# Makefile for compiling image.c

CC = g++
#CFLAGS = -Wall -Wextra -std=c11 -g -Og
CFLAGS = -Wall -Wextra -std=c++11 -g  -Og #-flto -pg
TARGET = image
SRC = image.c spiht_re.c spiht.cc $(SPIHT)/Encoder.cc \
	$(SPIHT)/Decoder.cc $(SPIHT)/Image.cc \
	$(SPIHT)/SPIHT.cc $(SPIHT)/BitIOStream.cc \
	$(SPIHT)/Exception.cc $(SPIHT)/FileIOStream.cc
SPIHT = ../imshrinker/src
INCLUDE = -I. -I../wavelib/header -I$(SPIHT)
LIBS = $(PWD)/../build/lib64/libwavelib.a -lm

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(INCLUDE)  -o $(TARGET) $(SRC) $(LIBS)

clean:
	rm -f $(TARGET)