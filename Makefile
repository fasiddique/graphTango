# See LICENSE.txt for license details.

C = g++
CFLAGS = -O2 -g -std=c++11 -fpermissive

CXX = g++
CXXFLAGS = -std=c++11 -g -pthread -fopenmp -O3
#CXXFLAGS = -std=c++11 -g -pthread -fopenmp -O0
#CXXFLAGS = -std=c++11 -g -pthread -fopenmp
#CXXFLAGS = -std=c++11 -g -pthread

DYN_PREFIX := d_

DYN_DIR := src/dynamic
UTL_DIR := src/common
OBJ_DIR := obj
BIN_DIR := bin

DYN_SRC := $(wildcard $(DYN_DIR)/*.cc)
DYN_SRC += $(wildcard $(DYN_DIR)/*.c)
DYN_SRC += $(wildcard $(UTL_DIR)/*.cc)
DYN_HDR := $(wildcard $(DYN_DIR)/*.h)
DYN_HDR += $(wildcard $(UTL_DIR)/*.h)

DYN_OBJ := $(addprefix $(OBJ_DIR)/$(DYN_PREFIX),$(notdir $(patsubst %.c,%.o,$(wildcard $(DYN_DIR)/*.c))))
DYN_OBJ += $(addprefix $(OBJ_DIR)/$(DYN_PREFIX),$(notdir $(patsubst %.cc,%.o,$(wildcard $(DYN_DIR)/*.cc))))

.PHONY : all
all : frontEnd 

frontEnd : $(DYN_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJ_DIR)/$(DYN_PREFIX)%.o : $(DYN_DIR)/%.cc $(DYN_HDR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/$(DYN_PREFIX)%.o : $(DYN_DIR)/%.c $(DYN_HDR)
	$(C) $(CFLAGS) -c $< -o $@

.PHONY : clean

clean:
	rm -f frontEnd
	rm -f $(OBJ_DIR)/*.o
	rm -f $(BIN_DIR)/*
