#  Copyright 2019 UniversalQCompiler (https://github.com/Q-Compiler/UniversalQCompiler)
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

from distutils.core import setup, Extension

# change these
sys = 'Linux-x86-64'
arch_bit = '64' # replace with 32 if you have a 32bit architecture
version = '11.0'

# this should be working for you
mathematica_folder = '/usr/local/Wolfram/Mathematica/' + version \
                     + '/SystemFiles/Links/WSTP/DeveloperKit/' + sys \
                     + '/CompilerAdditions'

module1 = Extension('MathematicaLink',
                    include_dirs = [mathematica_folder],
                    library_dirs = [mathematica_folder],
                    libraries = [
                        'WSTP' + arch_bit + 'i4',
                        'm', 'pthread', 'rt', 'stdc++', 'dl', 'uuid'
                    ],
                    sources = ['mlinkmodule.c'])

setup (name = 'MathematicaLink',
       version = '1.0',
       description = 'This is a package to call Mathematica functions',
       ext_modules = [module1])
