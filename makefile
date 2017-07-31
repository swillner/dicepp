GCC_OPTIONS := -ffp-contract=off -fassociative-math -fno-signed-zeros -fno-trapping-math -fno-signaling-nans -fno-strict-overflow
GCC_WARNINGS := -Wall -Wextra -Werror \
    -Wcast-align -Wpointer-arith -Wstrict-overflow=5 -Wundef -Wcast-qual -Wunreachable-code \
    -Wlogical-op -Wredundant-decls -Waddress -Wuninitialized -Winit-self -Wdisabled-optimization \
    -Wno-sign-compare -Wno-reorder -Wno-ignored-qualifiers -Wno-unknown-pragmas -Wfloat-equal -Wno-unused
YAML_CPP_FILES := $(wildcard lib/settingsnode/lib/yaml-cpp/src/*.cpp)
YAML_OBJ_FILES := $(patsubst lib/settingsnode/lib/yaml-cpp/src/%.cpp,bin/yaml/%.o,$(YAML_CPP_FILES))
CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(patsubst src/%.cpp,bin/%.o,$(CPP_FILES)) bin/borg/borg.o bin/borg/mt19937ar.o
.SECONDARY: OBJ_FILES
LD_FLAGS := -lstdc++ -lnetcdf_c++4 -lnetcdf -lnlopt -lipopt
CC_FLAGS := -std=c++11 -pthread -I include -I include/climate -I include/damage -I lib/settingsnode/include -I lib/settingsnode/lib/yaml-cpp/include -I lib/cpp-library -I lib/borg  -I /usr/include/eigen3 -I /usr/include/coin/ -DHAVE_CSTDDEF
DICEPP_VERSION := $(shell git describe --tags --dirty --always | sed -e 's/v\([0-9]\+\.[0-9]\+\)\.\(0-\([0-9]\+\)\)\?\(.*\)/\1.\3\4/')
DICEPP_FLAGS := $(GCC_WARNINGS) -DDICEPP_VERSION='"$(DICEPP_VERSION)"'
HEADERS_REQUIRING_REBUILD := $(wildcard include/*.h) $(wildcard include/**/*.h)

.PHONY: all fast debug gcc clean cleanall format

all: fast

fast: CC_FLAGS += -O3 -flto
fast: LD_FLAGS += -flto
fast: gcc

debug: CC_FLAGS += -g -DDEBUG
debug: gcc

gcc: CXX = g++
gcc: CC_FLAGS += $(GCC_OPTIONS)
gcc: dice

clean:
	@rm -f $(OBJ_FILES) dice

cleanall: clean
cleanall:
	@rm -rf bin

format:
	@clang-format -style=file -i $(wildcard src/*.cpp) $(wildcard src/*.h) $(wildcard test/*.cpp)

dice: $(YAML_OBJ_FILES) $(OBJ_FILES)
	@echo Linking to $@
	@$(CXX) -o $@ $^ $(CC_FLAGS) $(LD_FLAGS)

bin/%.o: src/%.cpp include/%.h $(HEADERS_REQUIRING_REBUILD)
	@mkdir -p bin
	@echo Compiling $<
	@$(CXX) $(CC_FLAGS) $(DICEPP_FLAGS) -c -o $@ $<

bin/borg/%.o: lib/borg/%.c lib/borg/%.h
	@mkdir -p bin/borg
	@echo Compiling $<
	@$(CXX) $(CC_FLAGS) $(DICEPP_FLAGS) -c -o $@ $<

bin/main.o: src/main.cpp $(HEADERS_REQUIRING_REBUILD)
	@mkdir -p bin
	@echo Compiling $<
	@$(CXX) $(CC_FLAGS) $(DICEPP_FLAGS) -c -o $@ $<

bin/yaml/%.o: lib/settingsnode/lib/yaml-cpp/src/%.cpp lib/settingsnode/lib/yaml-cpp/src/*.h lib/settingsnode/lib/yaml-cpp/include/yaml-cpp/*.h
	@mkdir -p bin/yaml
	@echo Compiling $<
	@$(CXX) $(CC_FLAGS) -c -o $@ $<
