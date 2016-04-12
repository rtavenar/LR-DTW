# List to '.PHONY' all fake targets, those that are neither files nor folders.
# "all" and "clean" are good candidates.
.PHONY: all, clean

# Define the final program name
PROGNAME := test_optima

# CFLAGS is used for C compilation options.
CFLAGS := -Wall -O0

# CXXFLAGS is used for C++ compilation options.
CXXFLAGS += -Wall -O0

# LDFLAGS is used for linker (-g enables debug symbols)
LDFLAGS  += -g

# List the project' sources to compile or let the Makefile recognize
# them for you using 'wildcard' function.
#
SOURCES = $(wildcard *.cc)

# List the project' headers or let the Makefile recognize
# them for you using 'wildcard' function.
#
HEADERS = $(wildcard %.h)

# Construct the list of object files based on source files using
# simple extension substitution.
OBJECTS = $(SOURCES:%.cc=%.o)

#
# Now declare the dependencies rules and targets
#
# Starting with 'all' make it  becomes the default target when none
# is specified on 'make' command line.
all : test_optima

# Declare that the final program depends on all objects and the Makfile
test_optima : test_optima.o p_dtw.o Makefile
	$(CXX) -o $@ $(LDFLAGS) test_optima.o p_dtw.o


%.o: %.cc $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(OUTPUT_OPTION) $<

# Simple clean-up target
# notes:
# - the '@' before 'echo' informs make to hide command invocation.
# - the '-' before 'rm' command to inform make to ignore errors.
clean :
	@echo "Clean."
	-rm -f *.o $(PROGNAME)