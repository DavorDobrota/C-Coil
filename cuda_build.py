import os
from distutils.errors import DistutilsExecError, CompileError
from distutils.msvccompiler import MSVCCompiler
from distutils.unixccompiler import UnixCCompiler
from os.path import join
from turtle import st
import types
from typing import List, Optional, Union

from setuptools.command.build_ext import build_ext


def find_in_path(names: List[str], path: str) -> Optional[str]:
    for dir in path.split(os.pathsep):
        print(f"Looking for NVCC in: {dir}")
        for name in names:
            binpath = join(dir, name)
            if os.path.exists(binpath):
                return os.path.abspath(binpath)
    return None


def locate_CUDA() -> dict:
    if 'CUDAHOME' in os.environ:
        home = os.environ['CUDAHOME']
        nvcc = join(home, 'bin', 'nvcc')
    else:
        nvcc = find_in_path(['nvcc', 'nvcc.exe'], os.environ['PATH'])
        if not nvcc:
            raise EnvironmentError(
                'The nvcc binary could not be located in your $PATH. Either add it to your path, or set $CUDAHOME'
            )
        home = os.path.dirname(os.path.dirname(nvcc))

    cuda_config = {
        'home': home, 'nvcc': nvcc,
        'include': join(home, 'include'),
        'lib64': join(home, 'lib64'), 'lib': join(home, 'lib')
    }

    for k, v in cuda_config.items():
        if not os.path.exists(v):
            print('Warning: The CUDA %s path could not be located in %s' % (k, v))

    return cuda_config


def customize_compiler_for_nvcc(
        compiler: Union[UnixCCompiler, MSVCCompiler], CUDA: dict, cuda_compile_args: List[str],
        cuda_include_dirs_generated: List[str], cuda_macros_generated: List[str]
) -> None:
    def setup_unix_compiler():
        default_compiler_so = compiler.compiler_so
        _super = compiler._compile

        def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
            # Set compiler to NVCC for .cu files
            if os.path.splitext(src)[1] == '.cu':
                compiler.set_executable('compiler_so', CUDA['nvcc'])
                extra_postargs = cuda_compile_args + ["--compiler-options=-fpic"]

            _super(obj, src, ext, cc_args, extra_postargs, pp_opts)

            # Reset compiler to default
            compiler.compiler_so = default_compiler_so

        compiler._compile = _compile

    def setup_msvc_compiler():
        compiler.nvcc = CUDA['nvcc']
        compiler.cuda_compile_args = cuda_compile_args
        compiler.cuda_include_dirs = cuda_include_dirs_generated
        compiler.cuda_macros = cuda_macros_generated
        compiler.compile = types.MethodType(_MSVC_CUDA_compile, compiler)

    # Determine the compiler type and set it up
    compiler.src_extensions.append('.cu')

    if compiler.compiler_type == "unix":
        print("Setting up Unix compiler for CUDA building.")
        setup_unix_compiler()
    elif compiler.compiler_type == "msvc":
        print("Setting up MSVC compiler for CUDA building.")
        setup_msvc_compiler()
    else:
        raise Exception("Trying to build CUDA code with an unsupported compiler.")


def _MSVC_CUDA_compile(
    compiler: MSVCCompiler,
    sources: List[str], output_dir: Optional[str] = None,
    macros: Optional[List[tuple]] = None, include_dirs: Optional[List[str]] = None,
    debug: bool = False,
    extra_preargs: Optional[List[str]] = None, extra_postargs: Optional[List[str]] = None,
    depends: Optional[List[str]] = None,
):
    if not compiler.initialized:
        compiler.initialize()
    compile_info = compiler._setup_compile(
        output_dir, macros, include_dirs, sources, depends, extra_postargs
    )
    macros, objects, extra_postargs, pp_opts, build = compile_info

    compile_opts = extra_preargs or []
    compile_opts.append('/c')
    if debug:
        compile_opts.extend(compiler.compile_options_debug)
    else:
        compile_opts.extend(compiler.compile_options)

    for obj in objects:
        try:
            src, ext = build[obj]
        except KeyError:
            continue
        if debug:
            # pass the full pathname to MSVC in debug mode,
            # this allows the debugger to find the source file
            # without asking the user to browse for it
            src = os.path.abspath(src)

        if ext in compiler._c_extensions:
            input_opt = "/Tc" + src
        elif ext in compiler._cpp_extensions:
            input_opt = "/Tp" + src
        elif ext in compiler._rc_extensions:
            # compile .RC to .RES file
            input_opt = src
            output_opt = "/fo" + obj
            try:
                compiler.spawn([compiler.rc] + pp_opts + [output_opt] + [input_opt])
            except DistutilsExecError as msg:
                raise CompileError(msg)
            continue
        elif ext in compiler._mc_extensions:
            # Compile .MC to .RC file to .RES file.
            h_dir = os.path.dirname(src)
            rc_dir = os.path.dirname(obj)
            try:
                # first compile .MC to .RC and .H file
                compiler.spawn([compiler.mc] + ['-h', h_dir, '-r', rc_dir] + [src])
                base, _ = os.path.splitext(os.path.basename(src))
                rc_file = os.path.join(rc_dir, base + '.rc')
                # then compile .RC to .RES file
                compiler.spawn([compiler.rc] + ["/fo" + obj] + [rc_file])

            except DistutilsExecError as msg:
                raise CompileError(msg)
            continue
        elif ext == '.cu':
            # Compile CUDA file
            input_opt = src
            output_opt = f'-o={obj}'
            try:
                compiler.spawn(
                    [compiler.nvcc] + ["-c"] + compiler.cuda_compile_args 
                    + compiler.cuda_include_dirs + compiler.cuda_macros
                    + [output_opt] + [input_opt])
            except DistutilsExecError as msg:
                raise CompileError(msg)
            continue
        else:
            raise CompileError(
                "Don't know how to compile {} to {}".format(src, obj)
            )

        output_opt = "/Fo" + obj
        try:
            compiler.spawn(
                [compiler.cc]
                + compile_opts
                + pp_opts
                + [input_opt, output_opt]
                + extra_postargs
            )
        except DistutilsExecError as msg:
            raise CompileError(msg)

    return objects


class CudaBuildExt(build_ext):
    CUDA: dict
    cuda_compile_args: List[str]
    cuda_include_dirs_generated: List[str]
    cuda_macros_generated: List[str]

    def build_extensions(self):
        customize_compiler_for_nvcc(
            self.compiler, CudaBuildExt.CUDA, CudaBuildExt.cuda_compile_args,
            CudaBuildExt.cuda_include_dirs_generated, CudaBuildExt.cuda_macros_generated
        )
        build_ext.build_extensions(self)

    @staticmethod
    def setup(CUDA: dict, cuda_compile_args: List[str], cuda_include_dirs: List[str], cuda_macros: List[tuple]):
        CudaBuildExt.CUDA = CUDA
        CudaBuildExt.cuda_compile_args = cuda_compile_args
        CudaBuildExt.cuda_include_dirs_generated = [f"-I{inc_dir}" for inc_dir in cuda_include_dirs]
        CudaBuildExt.cuda_macros_generated = [f"-D{key}={value}" for key, value in cuda_macros]
