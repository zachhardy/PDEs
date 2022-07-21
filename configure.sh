#!/bin/sh

DO_CLEAN="NO"   # remove the build directory first
DO_CMAKE="YES"  # run the CMake generator

# Loop over argmuents
for arg in "$@"; do

  # Parse clean argument
  if [ "$arg" = "clean" ]; then
    DO_CLEAN="YES"

    # If clean is the only argument, no CMake
    if [ $# -eq 1 ]; then
      DO_CMAKE="NO"
    fi
  fi

  # Parse reset argument
  if [ "$arg" = "reset" ]; then
    DO_CLEAN="YES"
    DO_CMAKE="YES"
  fi
done

# Perform actions
if [ $DO_CLEAN = "YES" ]; then
  if [ -d "build" ]; then
    rm -r build
  fi
  if [ -d "bin" ]; then
    rm -r bin
  fi

  cd Problems || exit
  for problemset in *; do
    cd "$problemset" || exit
    if [ -d "bin" ]; then
      rm -r bin
    fi
    cd ..
  done
  cd ..
fi

if [ $DO_CMAKE = "YES" ]; then
  mkdir -p build && mkdir -p bin

  cd Problems || exit
  for problemset in *; do
    if [ ! -d "bin" ]; then
      mkdir -p "${problemset}/bin"
    fi
  done

  cd ../build || exit
  cmake .. && cd ..
fi
