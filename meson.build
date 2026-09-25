project(
  'gammcor-integrals',
  'fortran',
  version : '1.0',
  meson_version : '>=1.1',
  default_options : [
    'warning_level=3',
    'buildtype=plain',
    'backend=ninja',
  ]
)

if meson.backend() != 'ninja'
  error('This project only supports the Ninja backend.')
endif

python = import('python').find_installation('python3')

cpu = get_option('cpu')
cpu_info = cpu
if cpu == 'auto'
  cpu_lines = run_command(python, files('.meson/scripts/detect_cpu.py'), check : true).stdout().strip().split('\n')
  cpu = cpu_lines[0]
  cpu_info = cpu + ' (auto: ' + cpu_lines[1] + ')'
endif
cpu_flags = {
  'avx512'     : ['-xCORE-AVX512', '-qopt-zmm-usage=high'],
  'avx512-amd' : ['-march=x86-64-v4'],
  'avx2'       : ['-march=x86-64-v3'],
  'native'     : ['-xHost'],
  'dispatch'   : ['-march=x86-64-v3', '-axCORE-AVX512', '-qopt-zmm-usage=high'],
}[cpu]
add_project_arguments(cpu_flags, language : 'fortran')
add_project_link_arguments(cpu_flags, language : 'fortran')

subdir('src')

# Static library for Gammcor
cholesky_lib = static_library(
  'cholesky',
  library_sources,
  name_prefix : '',
  install : true,
)

# Test executable (named test_bin internally to avoid collision with Meson's reserved 'test' target)
test_exe = executable(
  'test_bin',
  test_source,
  link_with : cholesky_lib,
  install : false,
)

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

summary({'cpu' : cpu_info}, section : 'gammcor-integrals')
