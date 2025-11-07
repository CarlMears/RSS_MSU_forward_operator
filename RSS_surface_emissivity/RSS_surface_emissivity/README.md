# Read Me for `geomod10`

This is a simple Python wrapper about Frank's `geomod10b` and `geomod10c` routines.

At this point, only the wind emissivity part has been attached to.
(in `geomod10py.f90`)
Other attachment points could be written.

## Building

The [build](https://pypa-build.readthedocs.io/en/stable/) tool is used to
generate a Python wheel. Note that a Fortran compiler (such as `gfortran`) and
the Python development headers are required to build the Python extension.

```bash
# Install build
pip install build

# Build the wheel, the output is in dist/
python3 -m build --wheel

# The wheel can then be installed using:
pip install dist/geomod10-*.whl
```

The GitLab CI is configured to build `manylinux_2_17` (aka `manylinux2014`)
wheels, which work with a [wide variety of common Linux
distributions](https://github.com/mayeut/pep600_compliance). The wheels are
built for the `x86_64` architecture for Python versions 3.9, 3.10, 3.11, and
3.12. Following [the numpy ecosystem
policy](https://numpy.org/neps/nep-0029-deprecation_policy.html#drop-schedule),
Python 3.9 is the minimum version supported.

## Python API

~~~
emiss = wind_emiss(freq,tht,sst,wind,phir)
~~~

`freq` and `tht` are scalars

`sst`, `wind`, and `phir` are 1-D numpy arrays

## Python test code

Some simple tests are in `tests`. After installing the `geomod10` package, the tests can be run using:

```
pytest -v
```

## Modifications I needed to make to the fortran code
* changed comments to "free" format
* wrapped some lines that were too long
* changed the `open` statements to reflect linux mounts and standard fortran (the original `openbig` statements are commented out)
* added in the `cosd` and `sind` routines to `geomod10c.f90`

## Requirements
* numpy 
* The O: drive on my machine is mounted as `/mnt/oserver/o`
