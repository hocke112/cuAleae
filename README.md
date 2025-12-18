# cuAleae

By Jeremiah Hockett, Chiemeka Nwakama, Aditya Parida, and Cheo Cedillo.

A modern C++/CUDA C hybrid port (or rewrite) of Aleae, a chemical reaction network simulator written by Marc Riedel. It adds GPU parallelism to individual trial, adding significant speedup to each trial. This was done as the final project for EE 5351.

This was written with pernission from Marc Riedel.

## Motivation

Aleae is rather slow and is known to have memory-related crashes. Although cuAleae was written to address the former, reliability was also a concern. A different CRN simulator, [MARlea](https://github.com/nadaso8/MARlea/), achieves a speedup via CPU parallelism across trials while ensuring reliability.

As stated before, cuAleae achieves its speedup by GPU parallelism within trials. Crashes and other memory-related issues are eliminated in cuAleae.

## Command-line Input

The format of cuAleae's command-line input is identical to Aleae's, which is to maintain compability to Aleae. Said format is as follows:

```./cuAleae <.in file> <.r file> <num trials> <max time> <verbosity bits> <max steps>```

The first two are the input files to cuAleae, while the rest are explained below:
* max time: a positive integer time limit or -1 for no max time
* max steps: a positive integer step limit or -1 for not max steps
* verbosity: an argument made up of 4-bit setting (0 to 15 inclusive)
    * PRINT_TRIALS - 1
    * PRINT_TERMINAL - 2
    * PRINT_TRACE - 4
    * PRINT_STATES - 8

## System Requirements
### Operating System
* Windows 10/11
* Linux
### Hardware
* An nVidia GPU that supports CUDA 9.0 or higher (Makefile can be tweaked to support lower versions like CUDA 8.0)
### Toolchain
* A C++ compiler that can support C++17 or higher
* A CUDA C Compiler that can support sm_70 or higher (CUDA 9.0 or higher)

## Compilation
To compile the source code, simply type `make build` to compile. Alternatively, you can type `make debug` for compile for debugging.