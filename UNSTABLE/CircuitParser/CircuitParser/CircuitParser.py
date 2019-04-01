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

from projectq import ops
from projectq.cengines import MainEngine

from CircuitParser.Gates import OneQubitGate, Rx, Ry, Rz, CNot, ComposedGate

def mathematica_to_python(s):

    def _parser(ls):
        gate = [Rx, Ry, Rz]
        a,b,c = round(ls[0].real), float(ls[1].real), float(ls[2].real)
        if a == 0:
            return CNot(control=int(b-1), target=int(c-1))
        elif a > 0:
            return gate[a-1](channel=int(c-1), angle=float(-b))
        else:
            raise KeyError

    def _string_to_array(s):
        s = s.replace(' ', '')
        s = s.replace('\n', '')
        s = s.replace('*I', 'j')
        s = s.replace('*^', 'e')
        assert s.find('^') < 0
        assert s.find('*') < 0
        fst, empty = _parse_array(s)
        assert empty == ''
        return fst

    assert len(s.split('=')) <= 2, 'The Mathematica output should only contain one variable'

    # ignore what precedes the = symbol
    s = s.split('=')[-1]
    return list(map(_parser, _string_to_array(s)))

def python_to_projectq(gates, eng = MainEngine()):
    gates = gates.copy()
    gates.reverse()

    qubits = gates.max_channel() + 1

    circuit = eng.allocate_qureg(qubits)

    own2projectq = {Rx: ops.Rx, Ry: ops.Ry, Rz: ops.Rz, CNot: ops.CNOT}

    for g in gates:
        for gate_own, gate_projq in own2projectq.items():
            if isinstance(g, gate_own):
                if isinstance(g, OneQubitGate):
                    gate_projq(g.angle) | circuit[g.channel]
                else:
                    gate_projq | (circuit[g.channel2], circuit[g.channel1])
                break
        else:
            raise KeyError

    for q in range(qubits):
        ops.Measure | circuit[q]

    eng.flush()

    qasm = eng.backend.qasm
    qasm = ("\ninclude \"qelib1.inc\";\nqreg q[{nq}];\ncreg c[{nq}];{qasm}\n").format(
            nq=qubits, qasm=qasm)

    return qasm, circuit

def mathematica_to_qasm(s, eng=None):
    if not eng:
        eng = MainEngine()

    gates = ComposedGate(mathematica_to_python(s))
    qubits = gates.max_channel() + 1
    qasm, qureg = python_to_projectq(gates, eng)

    return qasm, qureg

def _parse_array(s):
    if s[0] == '{':
        s = s[1:]
        ret = []
        i = 1
        while s[0] != '}':
            fst, s = _parse_array(s)
            ret.append(fst)
        s = s[1:].lstrip(',')
        return ret, s
    else:
        fst, rst = _splitAt(',}', s)
        rst = rst.lstrip(',')
        fst = fst.strip('"')

        try:
            return complex(fst), rst
        except ValueError:
            raise ValueError("Mathematica return value not understood: "+fst)

def _splitAt(sep, s):
    i = 0
    while i < len(s) and s[i] not in sep: i += 1
    return s[:i], s[i:]
