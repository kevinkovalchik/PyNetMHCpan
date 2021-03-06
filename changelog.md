# Changelog

This changelog is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [0.1.5] 2020-12-14
### Fixed
- Wrong path was used to make temp directory.

## [0.1.4] 2020-12-12
### New
- A unique temp directory is used for every analysis. This makes it possible to run multiple instances at once,
e.g. on a computing cluster.

## [0.1.3] 2020-12-11
### Fixed
- Config file was _still_ missing when installed from PyPI because I messed up the filename of MANIFEST.in

## [0.1.2] 2020-12-11
### Fixed
- Config file was missing when installed from PyPI

## [0.1.1] 2020-12-11
### New
- Added `annotate_file` function to `NetMHCpan.Helper` class. This reads an input file (e.g. a PSM file), makes 
binding predictions based on the contents of the indicated peptide column (default is `Peptide`), and in a new 
 compy of the file writes the binding percent rank from NetMHCpan as a new column.

## [0.1.0] 2020-12-10
### New
- Everything! First release of the project.