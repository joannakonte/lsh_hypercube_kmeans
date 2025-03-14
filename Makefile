# ================== | Paths | ================== #

SRC = ./src
MODULES = ./modules
HEADERS = ./headers
BUILD = ./build

# =========== | Compiler and Linker | =========== #

CC = g++

# Compile Options
CPPFLAGS = -std=c++11 -g3 -O2 -I$(HEADERS)

# Linker Options
LDFLAGS =

# =============== | Executables | =============== #

EXEC_LSH = lsh
EXEC_CUBE = cube
EXEC_CLUSTER = cluster


# ================ | Arguments | ================= #

LSH_ARG = 	-d ../input/train-images.idx3-ubyte 		\
			-q ../input/t10k-images.idx3-ubyte 			\
			-o output_lsh.txt	

CUBE_ARG =	-d ../input/train-images.idx3-ubyte 		\
			-q ../input/t10k-images.idx3-ubyte 			\
			-o output_cube.txt -N 5 -R 2000 -M -1 -probes 2 -k 5

CLUSTER_ARG = -i ../input/train-images.idx3-ubyte 		\
			  -c cluster.conf							\
			  -o output_cluster.txt						\
			  -m LSH

# =============== | Object Files | =============== #

OBJ_LSH = $(BUILD)/main_lsh.o $(BUILD)/lsh.o $(BUILD)/brute_force.o $(BUILD)/common.o $(BUILD)/dataset_input.o
OBJ_CUBE = $(BUILD)/main_cube.o $(BUILD)/cube_projection.o $(BUILD)/brute_force.o $(BUILD)/common.o $(BUILD)/dataset_input.o
OBJ_CLUSTER = $(BUILD)/main_cluster.o $(BUILD)/cluster.o $(BUILD)/common.o  $(BUILD)/lsh.o $(BUILD)/cube_projection.o $(BUILD)/dataset_input.o

# =============== | Compile | =============== #

# Create the build folder if it doesn't exist
$(shell mkdir -p $(BUILD))


all: $(EXEC_LSH) $(EXEC_CUBE) 

# LSH
$(EXEC_LSH): $(OBJ_LSH)
	$(CC) $(LDFLAGS) $(OBJ_LSH) -o $(EXEC_LSH) 

# CUBE
$(EXEC_CUBE): $(OBJ_CUBE)
	$(CC) $(LDFLAGS) $(OBJ_CUBE) -o $(EXEC_CUBE) 

# CLUSTER
$(EXEC_CLUSTER): $(OBJ_CLUSTER)
	$(CC) $(LDFLAGS) $(OBJ_CLUSTER) -o $(EXEC_CLUSTER) 

# Pattern rule for building object files from SRC
$(BUILD)/%.o: $(SRC)/%.cpp 
	$(CC) $(CPPFLAGS) -c $< -o $@

# Pattern rule for building object files from MODULES
$(BUILD)/%.o: $(MODULES)/%.cpp 
	$(CC) $(CPPFLAGS) -c $< -o $@

# =================== | Run | ==================== #

run_lsh: $(EXEC_LSH)
	./$(EXEC_LSH) $(LSH_ARG)

run_cube: $(EXEC_CUBE)
	./$(EXEC_CUBE) $(CUBE_ARG)

run_cluster: $(EXEC_CLUSTER)
	./$(EXEC_CLUSTER) $(CLUSTER_ARG)

# =================== | Clean | ================== #

# Remove object files and executables
clean:
	rm -f $(OBJ_LSH) $(EXEC_LSH) $(OBJ_CUBE) $(EXEC_CUBE) $(OBJ_CLUSTER) $(EXEC_CLUSTER) output*.txt
