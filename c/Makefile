# List to '.PHONY' all fake targets, those that are neither files nor folders.
# "all" and "clean" are good candidates.
.PHONY: all, clean

# CFLAGS is used for C compilation options.
CFLAGS := -Wall -O3

# CXXFLAGS is used for C++ compilation options.
CXXFLAGS += -Wall -O3

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
all : test_ucr test_optima test_lr_dtw test_1nn

test_optima : test_optima.o lr_dtw.o Makefile
	$(CXX) -o $@ $(LDFLAGS) test_optima.o lr_dtw.o
test_ucr : test_ucr.o ucr.o Makefile
	$(CXX) -o $@ $(LDFLAGS) test_ucr.o ucr.o
test_p_dtw : test_lr_dtw.o lr_dtw.o Makefile
	$(CXX) -o $@ $(LDFLAGS) test_lr_dtw.o lr_dtw.o
test_1nn : test_1nn.o lr_dtw.o ucr.o 1nn.o Makefile
	$(CXX) -o $@ $(LDFLAGS) test_1nn.o lr_dtw.o ucr.o 1nn.o


%.o: %.cc $(HEADERS) Makefile
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $(OUTPUT_OPTION) $<

# Simple clean-up target
# notes:
# - the '@' before 'echo' informs make to hide command invocation.
# - the '-' before 'rm' command to inform make to ignore errors.
clean :
	@echo "Clean."
	-rm -f *.o