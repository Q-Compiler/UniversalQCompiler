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

from MathematicaLink import open, close, eval, import_package
from CircuitParser import Gates
from CircuitParser.CircuitParser import mathematica_to_python

from os.path import join, isfile
from os import getcwd
from numpy import ndarray

def fulldecomp(mat):
    # FIX HERE !!!! only take the numbers and dont parse the gates

    def mathematica_setup():
        path_to_package = join(getcwd(), 'QDecomposition/Decompositions.m')
        if not isfile(path_to_package):
            raise FileNotFoundError

        # connect to Mathematica
        open()

        # import Decompositions.m
        import_package(path_to_package)

    try:
        mathematica_setup()
        gates = mathematica_to_python(query('QSD', mat))
    finally:
        close()
        print('Connection closed')
    return gates

def query(func, *args, suffix=''):
    """
    builds a string from list-type arguments to send to Mathematica,
    which will interpret it as an expression.
    
    Returns the result from Mathematica, formatted as list
    """
    command = func + '['
    
    # `first` serves to know when to add commas in the expression
    first = True
    for x in args:
        if not first: command += ','
        if type(x) is ndarray: x = x.tolist()
        command += array_to_string(x)
        first = False
    command += ']'
    command += suffix
    cmd = eval(command)
    print(cmd)
    return mathematica_to_python(cmd)

def array_to_string(xs):
    s = str(xs)
    s = s.replace('[', '{')
    s = s.replace(']', '}')

    s = s.replace('j', '*I')
    s = s.replace('e', '*^')
    return s
