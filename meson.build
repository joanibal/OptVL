project(
  'optvl',
  'c',
  version: '1.2.1',
  license: 'GPL-3.0',
  meson_version: '>= 0.64.0',
  default_options: [
    'buildtype=debugoptimized',
    'c_std=c99',
    'fortran_std=legacy',
  ],
)

cc = meson.get_compiler('c')

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)
ff = meson.get_compiler('fortran')
# if ff.has_argument('-Wno-conversion')
#   add_project_arguments('-Wno-conversion', language: 'fortran')
# endif

is_windows = host_machine.system() == 'windows'


# Platform detection
is_mingw = is_windows and cc.get_id() == 'gcc'
if is_mingw and ff.get_id() != 'gcc'
  error('If you are using GCC on Windows, you must also use GFortran! Detected ' + ff.get_id())
endif

cython_c_args = ['-DCYTHON_CCOMPLEX=0'] # see gh-18975 for why we need this
if is_mingw
  # For mingw-w64, link statically against the UCRT.
  gcc_link_args = ['-lucrt', '-static']
  add_project_link_arguments(gcc_link_args, language: ['c', 'cpp', 'fortran'])
  # Force gcc to float64 long doubles for compatibility with MSVC
  # builds, for C only.
  add_project_arguments('-mlong-double-64', language: 'c')
  # Make fprintf("%zd") work (see https://github.com/rgommers/scipy/issues/118)
  add_project_arguments('-D__USE_MINGW_ANSI_STDIO=1', language: ['c', 'cpp'])
  # Silence warnings emitted by PyOS_snprintf for (%zd), see
  # https://github.com/rgommers/scipy/issues/118.
  # Use as c_args for extensions containing Cython code
  cython_c_args += ['-Wno-format-extra-args', '-Wno-format']
  # Flag needed to work around BLAS and LAPACK Gfortran dependence on
  # undocumented C feature when passing single character string arguments. See:
  #   https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90329
  #   https://github.com/wch/r-source/blob/838f9d5a7be08f2a8c08e47bcd28756f5d0aac90/src/gnuwin32/MkRules.rules#L121
  add_project_arguments('-fno-optimize-sibling-calls', language: ['fortran'])
endif


# _linker_script = meson.project_source_root() / '_build_utils/link-version-pyinit.map'
# link_args = ['-Wl,--version-script=' + _linker_script]
# link_args = ['-Wl,--image-base,0x10000000']
# link_args = ['-mcmodel=large']
link_args = []
# -shared -fPIC
# gfortran -o myprogram mysource.f90 -Wl,-stack_size,8388608

# https://mesonbuild.com/Python-module.html
py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
py3_dep = py3.dependency()


incdir_numpy = run_command(py3,
    [
    '-c',
    'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
    ],
    check: true
).stdout().strip()
inc_np = include_directories(incdir_numpy)

incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)


# Fortran warning flags

_fflag_Wno_maybe_uninitialized = ff.get_supported_arguments('-Wno-maybe-uninitialized')
_fflag_Wno_unused_variable = ff.get_supported_arguments('-Wno-unused-variable')
_fflag_Wno_align_commons = ff.get_supported_arguments('-Wno-align-commons')
_fflag_Wno_unused_label = ff.get_supported_arguments('-Wno-unused-label')
_fflag_Wno_character_truncation = ff.get_supported_arguments('-Wno-character-truncation')
_fflag_Wno_unused_dummy_argument = ff.get_supported_arguments('-Wno-unused-dummy-argument')

# The default list of warnings to ignore from Fortran code. There is a lot of
# old, vendored code that is very bad and we want to compile it silently (at
# least with GCC and Clang)
ff_args = ff.get_supported_arguments(
 _fflag_Wno_maybe_uninitialized,
 _fflag_Wno_unused_variable,
 _fflag_Wno_align_commons, 
 _fflag_Wno_unused_label,
 _fflag_Wno_character_truncation,
 _fflag_Wno_unused_dummy_argument
)


ff_args += [
  '-ffixed-line-length-80', 
  '-std=legacy',
  '-fdefault-real-8',
  '-fdefault-double-8',
  '-fPIC',
  '-O2'
]

numpy_nodepr_api = '-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION'

cc_args = [numpy_nodepr_api]

# # Warning: Array 'aict' at (1) is larger than limit set by '-fmax-stack-var-size=', moved from stack to static storage. This makes the procedure unsafe when called recursively, or concurrently from multiple threads. Consider using '-frecursive', or increase the '-fmax-stack-var-size=' limit, or change the code to use an ALLOCATABLE 

if is_windows
  ff_args += ff.get_supported_arguments('-fomit-frame-pointer')
  cc_args += cc.get_supported_arguments('-fomit-frame-pointer')
  # ff_args += ['-fmax-stack-var-size=6553']
endif


fortranobject_c = incdir_f2py / 'fortranobject.c'
# Share this object across multiple modules.
fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  c_args: cc_args,
  dependencies: py3_dep,
  include_directories: [inc_np, inc_f2py],
)
fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [inc_np, inc_f2py],
)

c = run_command('src' / 'grab-all-fortran-files.py', check: true)
avl_source_files = c.stdout().strip().split('\n')

# Specify the directory containing your Meson build script
incs_dir = 'src/includes'
avl_inc_dir = include_directories(incs_dir)

# Define the input and output file paths
inc_in = join_paths(incs_dir, 'AVL.INC.in')
inc_out = join_paths(incs_dir, 'AVL.INC')

run_command('src' / 'create_avl_inc.py', check: true)



avl_c_wrapper = custom_target('libavlmodule.c',
input : ['src/f2py/libavl.pyf'],
output : ['libavlmodule.c', 'libavl-f2pywrappers.f'],
command: [py3, '-m', 'numpy.f2py',  '@INPUT@', ]
)

py3.extension_module('libavl',
    avl_source_files,
    avl_c_wrapper,
    fortranobject_c,
    include_directories: [inc_np, inc_f2py, avl_inc_dir],
    dependencies: [fortranobject_dep],
    subdir: 'optvl',
    link_args: link_args,
    fortran_args: ff_args,
    c_args: cc_args,
    install : true,
    link_language: 'fortran'
                   
)

py3.install_sources([
  'optvl/__init__.py',
  'optvl/optvl_class.py',
  'optvl/MExt.py',
  'optvl/om_wrapper.py',
],
  subdir: 'optvl'
)


        
