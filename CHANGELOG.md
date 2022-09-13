# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2.0.1] - 2022-09-14
### Added
- Ability to calculate random effects estimates separately for each permutation.
- Simulate 'A' effect in FEMA_run_on_synthetic_data.m.
- light_save implementation in FEMA_wrapper.m

### Fixed
- Fixed the call to FEMA_synthesize discussed in [issue #1](cmig-research-group/cmig_tools#1)
- FEMA_fit can now run with only two random efefcts, as discussed in [issue #2](https://github.com/cmig-research-group/cmig_tools/issues/2)

## [2.0.0] - 2022-04-12
### Added
- Initial release of public repo.

[Unreleased]: https://github.com/cmig-research-group/cmig_tools_internal/compare/beta...public
[2.0.0]: https://github.com/cmig-research-group/cmig_tools/releases/tag/2.0
