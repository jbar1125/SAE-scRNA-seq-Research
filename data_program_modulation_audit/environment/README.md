# Validation-instance environment

`requirements.txt` records the installed packages from the Python 3.12 Linux
virtualenv used for the original-suite validation, counterexamples and calibration.
`validation_environment.json` records its platform, principal versions and the
freeze-file hash. This is a resolved instance record, not a claim of portability
to arbitrary Python versions, operating systems or CUDA devices.

PyTorch was installed from its official CPU wheel index before the repository
requirements. The reusable cloud setup uses that same order. The downstream
counterexamples and identity calibration do not train SAEs; their own metadata
also records the specific NumPy/SciPy/statsmodels versions and source hashes used.
Normal CPU development was limited to four BLAS/OpenMP threads; the numerical
simulation scripts used one thread. No credentials are included.
