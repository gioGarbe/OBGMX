# Changelog

## [1.1] - 2025-10-24

The _Nobel prize_ release.  
This release should enable compilation of the `obgmx` executable on modern platforms.

Tested on Ubuntu 24.03 LTS and macOS 26 "Tahoe".

### Added

* updated patch
* installation script
* `examples/` directory
* `Ubuntu 16.04 LTS` directory with the old release
	
## [1.0] - 2023-12-05

Initial release.
Notice that it might not compile on modern distributions (and it usually does not).

### Added

* patches to openbabel-2.3.2 (`ob-gg-2.3.2.patch`)
* source code of obgmx (`obgmx.cpp` and `bondedparams.h`)
* statically linked executable, to be run in an Ubuntu 16.04 LTS environment
  - `share` directory to be copied to `/usr/local/ob-gg-2.3.2/`
