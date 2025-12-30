# Makefile for raytracerX11
# Cross-platform build for Linux, macOS, and WSL

CC = clang
CFLAGS = -O3 -fopenmp -Wall -Wextra -march=native
LDFLAGS =
LIBS = -lX11 -lm
TARGET = main
SRC = main.c

# Platform detection
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    # macOS with Homebrew LLVM (for OpenMP support)
    CC = /opt/homebrew/opt/llvm/bin/clang
    CFLAGS += -I/opt/X11/include -I/opt/homebrew/include
    LDFLAGS += -L/opt/X11/lib -L/opt/homebrew/lib
    LIBS += -lomp
else ifeq ($(UNAME_S),Linux)
    # Linux (native GCC or Clang)
    CC = gcc
    CFLAGS += -fopenmp
    LIBS += -lgomp
else
    # WSL or other Unix-like
    CC = gcc
    CFLAGS += -fopenmp
    LIBS += -lgomp
endif

# Build rules
all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(LDFLAGS) $(SRC) $(LIBS) -o $(TARGET)
	@echo "Build complete: ./$(TARGET)"
	@echo "Platform: $(UNAME_S)"

clean:
	rm -f $(TARGET)
	rm -rf *.dSYM .idea

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run

# GIT HELPER

MESSAGE = .

push: clean add commit
	git push

add:
	git add .

commit:
	git commit -a -m "$(MESSAGE)"

update:
	git submodule update --remote --merge
