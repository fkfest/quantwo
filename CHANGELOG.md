# Release notes

### Changed

*added permutations to ElemCo.jl generation.

### Fixed

### Added

## Version [v1.0.2] - 2024.08.05

### Changed

* separated diagram from expression.
* improved handling of input terms.
* changed to load4idx() calls in ElemCo.jl generation.
* bract definition in params.reg.

### Fixed

* commented out equation block before actual equation block in the input file with a very long line is now handled correctly.

### Added

* DC-CCSDT closed-shell generation test file.
* penalized virtual orbitals in @tensoropt calls.

## Version [v1.0.1] - 2024.02.27

### Changed

* removed _dummy variable from Tensor
* renamed TensorsSet to TensorsList
* improved test Makefile

### Fixed

* Linebreaks in input equations with explspin = 1
* implemented pointer comparison for set of pointers _residuals to get consistent orders with different gcc versions

## Version [v1.0.0] - 2024.02.21

Initial release
