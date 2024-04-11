# TODO

- consider passing context to population and individual
  - make population a factory for individuals, setting all necessary parameters
- consider if we need to break from state machine earlier,
  since last generation is discarded
- add more stats to print_result
  - runtime
  - generations
  - gen/s
- add printing genes to console in table / in columns
- bug: 2 newlines added after print_result with progress bar
- bug: 1 newline is added before print_result with progress bar
- add multi-threading
  - try "#pragma omp parallel for"
- add style guide
  - something standard and simple
  - best - with autoformatting
    - clang-format?
    - add to Makefile
