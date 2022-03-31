# tests_LISA_response
We want to compare different LISA response functions. In particular we want to compare the signals that are generated by different available softwares (e.g. LDC implementation of FastGB vs implementation of SyntheticLISA). Finally, we adapt the tri-linear approximation to each implementation and test its accuracy.

# Manual
## What do we want to do?

1. Replicate the _LDC implementation of FastGB_ and regenerate LDC data
2. Compare the signals to our _reference implementation_ (based on SyntheticLISA) that is adapted to LDC implementation
3. Adapte the _tri-linear approximation_ (TLA) to the LDC data

### 1. Regenerating LDC data

- `src/ldc_code.c` and `src/ldc_code.h` replicate the large LDC implementation of FastGB
- run `make` in command line to compile C code **&rarr;** `ldc_exe` is generated
- run `./ldc_exe M` to simulate `M` signals **&rarr;** random parameters and TDIs are generated (stored in `matlab/parameters.bin`, `matlab/X.bin`, `matlab/Y.bin`, `matlab/Z.bin`)

### 2.+3. Comparing to reference implementation and TLA

- script `matlab/test_TLA.m` compares the generated LDC data to ...
  - ... our reference model (adapted to LDC implementation of FastGB)
  - ... TLA (adapted to LDC implementation of FastGB)

# Comments

- all signals are generated with carrier frequency 0.001 Hz
- before generating new data with `ldc_exe`, the files `matlab/X.bin`, `matlab/Y.bin`, `matlab/Z.bin` and `matlab/parameters.bin` have to be deleted (in `scr/ldc_code.c` the data is continuously appended to the files)
