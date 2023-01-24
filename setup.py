import os
import re
import sys
import platform
import subprocess
import shutil

from distutils.version import LooseVersion
import distutils.command.build
from setuptools import setup, find_packages, Extension, Command
import setuptools.command.develop
import setuptools.command.install
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                         out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
                cfg.upper(),
                extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        if platform.system() == "Darwin" :
            if shutil.which('brew') == '/usr/local/bin/brew':
                cmake_args += ['-DCMAKE_C_COMPILER=/usr/local/bin/gcc-12', '-DCMAKE_CXX_COMPILER=/usr/local/bin/g++-12']
            else:
                cmake_args += ['-DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-12', '-DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-12']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        exedir = 'build/executable'
        if not os.path.exists(exedir):
            os.makedirs(exedir)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=exedir, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=exedir)
        print()  # Add an empty line for cleaner output



# This function downloads and builds the shared-library
def run_privateer_install_script():
    build_privateer_cmd = ['bash', 'DependencyInstaller.sh']
    if subprocess.call(build_privateer_cmd) != 0:
        sys.exit("Failed to build Privateer's dependencies")

# I make a new command that will build the shared-library
class build_privateerdeps(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        run_privateer_install_script()

# I subclass install so that it will call my new command
class install(setuptools.command.install.install):
    def run(self):
        self.run_command('build_privateerdeps')
        setuptools.command.install.install.run(self)

# I do the same for build...
class build(distutils.command.build.build):
    sub_commands = [
        ('build_privateerdeps', lambda self: True),
        ] + distutils.command.build.build.sub_commands

# ...and the same for develop
class develop(setuptools.command.develop.develop):
    def run(self):
        self.run_command('build_privateerdeps')
        setuptools.command.develop.develop.run(self)


setup(
    name='privateer',
    version='0.4',
    author='Jon Agirre',
    author_email='jon.agirre@york.ac.uk',
    description='The Swiss Army knife for carbohydrate structure validation, refinement and analysis',
    long_description='',
    packages=find_packages('src'),
    package_dir={'':'src'},
    ext_modules=[CMakeExtension('privateer/privateer')],
    cmdclass={
        'install': install,
        'build_privateerdeps': build_privateerdeps, 'build': build,
        'develop': develop,
        'build_ext': CMakeBuild,
    },
    zip_safe=False,
)

