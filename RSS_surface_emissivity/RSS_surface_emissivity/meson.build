project('RSS_surf_emiss', 'c', 'fortran',
  version : '0.0.1',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2'],
)

py = import('python').find_installation(pure: false, modules : ['numpy'])
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

f2py_interface = custom_target('RSS_surf_emissmodule.c',
  input : ['src/geomod10/geomod10py.f90'],
  output : ['RSS_surf_emissmodule.c', 'RSS_surf_emiss-f2pywrappers.f'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'RSS_surf_emiss', 'only:', 'wind_emiss', 'init', ':']
)

py.extension_module('RSS_surf_emiss',
  'src/geomod10/geomod10py.f90',
  'src/geomod10/geomod10b.f90',
  'src/geomod10/geomod10c.f90',
  f2py_interface,
  incdir_f2py / 'fortranobject.c',
  fortran_args : ['-fno-automatic', '-std=f2008'],
  include_directories: inc_np,
  dependencies : py_dep,
  install : true,
  subdir: 'RSS_surf_emiss',
)

py.install_sources(
    ['src/geomod10/data/em0_ref_freq_sst.dat',
     'src/geomod10/data/fit_emiss_phir_wind.dat',
     'src/geomod10/data/finetune_emiss_wind.dat',
     'src/geomod10/data/mk_scatterm_table_all.dat'],
    subdir: 'RSS_surf_emiss/data'
)
