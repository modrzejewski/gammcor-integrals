project(
  'gammcor-integrals',
  'fortran',
  version : '1.0',
  default_options : [
    'warning_level=3',
    'buildtype=plain',
    'backend=ninja',
  ]
)

if meson.backend() != 'ninja'
  error('This project only supports the Ninja backend.')
endif

subdir('src')

# Static library for Gammcor
cholesky_lib = static_library(
  'cholesky',
  library_sources,
  name_prefix : '',
  install : false,
)

# Test executable (named test_bin internally to avoid collision with Meson's reserved 'test' target)
test_exe = executable(
  'test_bin',
  test_source,
  link_with : cholesky_lib,
  install : false,
)

python = import('python').find_installation('python3')
post_build = files('.meson/scripts/post_build.py')

post_build_target = custom_target(
  'post_build_tasks',
  input : [cholesky_lib, test_exe],
  output : 'post_build.done',
  command : [
    python, post_build,
    '@INPUT0@', '@INPUT1@',
    meson.current_build_dir(),
    meson.project_source_root(),
    '@OUTPUT@',
  ],
  build_by_default : true,
)

