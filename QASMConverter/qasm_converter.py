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

import argparse
import sys
import numpy as np

from projectq import ops
from projectq.setups import ibm, default
from projectq.cengines import BasicEngine, MainEngine
from projectq.backends import IBMBackend
from projectq.meta import get_control_count, LogicalQubitIDTag

from parsimonious.grammar import Grammar
from parsimonious.nodes   import NodeVisitor, RegexNode, Node
from parsimonious.exceptions import VisitationError

from sys import stderr

def warning(s):
    print('Warning:', s, file=stderr)

# ==== Parser code ====

# grammar defining acceptable Mathematica input
# according to manual
# int = float, but different in parsing
grammar = Grammar(r"""
    circuit = open gates close
    gates   = (gate komma gates) / gate
    gate    = open (DIAG / CZ / CNOT / R / MST / TROUT / ANC / PS) close
    qubits  = (qubit komma qubits) / qubit
    qubit   = int
    floats  = (float komma floats) / float
    float   = ~"\s*([+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)\s*"
    int     = ~"\s*([+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)\s*"
    open    = ~"\s*{\s*"
    close   = ~"\s*}\s*"
    komma   = ~"\s*,\s*"

    DIAG    = ~"\s*(-2(\.\d*)?)" komma open floats close komma open qubits close
    CZ      = ~"\s*(-1(\.\d*)?)" komma qubit komma qubit
    CNOT    = ~"\s*(0(\.\d*)?)"  komma qubit komma qubit
    R       = ~"\s*([1-3](\.\d*)?)" komma float komma qubit
    MST     = ~"\s*(4(\.\d*)?)"  komma "1" komma qubit
    TROUT   = ~"\s*(4(\.\d*)?)"  komma "0" komma qubit
    ANC     = ~"\s*(5(\.\d*)?)"  komma ~"([0-1])" komma qubit
    PS      = ~"\s*(6(\.\d*)?)"  komma ~"([0-1])" komma qubit
""")

# some helper functions

# Having a nested list of the form [(),(),()...] flattens the tupels if they contain sub-tupels
def flatten_touples(list_in):
    lst = []
    for i in range(len(list_in)):
        if (type(list_in[i][0]) !=  tuple):
            lst.append(list_in[i])
        else:
            for j in range(len(list_in[i])):
                lst.append(list_in[i][j])
    return lst



def max_qubit(gate_list):
    # look at qubit numbers in gate_list to find largest qubit
    return max(y for x in map(
        lambda x: x[1],
        gate_list)
        for y in x
    )

def build_circuit(gate_list, eng):
    # creates circuit for tuple list of gates
    # cannot be done earlier, as the number of qubits must be
    # known a priori
    nb_qubits = max_qubit(gate_list) + 1

    qubits = eng.allocate_qureg(nb_qubits)

    for (gate, qs) in gate_list:
        qs = tuple(map(lambda q: qubits[q], qs))
        gate | qs

    return nb_qubits

def get_qasm(eng, nb_qubits):
    # returns qasm code from engine

    # flush
    eng.flush()

    qasm = eng.backend.qasm

    # adds cosmetic line at the beginning (necessary for IBM)
    qasm = ("\ninclude \"qelib1.inc\";\nqreg q[{nq}];\ncreg c[{nq}];{qasm}\n").format(
            nq=nb_qubits, qasm=qasm)
    return qasm

# helper function for all grammar rules that are sequences of type xs = (x '+' xs) / x
def visit_list(xs):
    # remove parenthesese
    xs = xs[0]
    # hack to differentiate between one-el list and multi-el
    if type(xs) != list:
        return [xs]
    else:
        return [xs[0]] + xs[2]

class InvalidRotationGateError(ValueError):
    pass

# transforms parsed tree into a gate list
class ToGateListVisitor(NodeVisitor):

    # for every type of node, run one of these parsers:
    def visit_circuit(self, node, circuit):
        return circuit[1]

    def visit_gates(self, node, gs):
        return visit_list(gs)

    def visit_gate(self, node, g):
        # remove parentheses
        g = g[1]
        return g[0]

    def visit_qubit(self, node, q):
        # doesnt seem to be used (== int)
        return q[0]

    def visit_qubits(self, node, qs):
        return visit_list(qs)

    def visit_float(self, node, f):
        return float(node.match.group(1))

    def visit_floats(self, node, fs):
        return visit_list(fs)

    def visit_int(self, node, i):
        return round(float(node.match.group(1)))

    def visit_DIAG(self, node, params):
        warning("Diagonal gates cannot be translated. Please decompose it using UniversalQCompiler with the method DecDiagGate")
        return ("DIAG", (-1,))

    def visit_CZ(self, node, params):
        ctr,targ = params[2], params[4]
        return (ops.CZ, (ctr,targ))

    def visit_CNOT(self, node, params):
        ctr,targ = params[2], params[4]
        return (ops.CNOT, (ctr,targ))

    def visit_R(self, node, params):
        # type parameter
        t = params[0]
        typ = round(float(t))
        # other params
        angle = params[2]
        q   = params[4]

        # different angle convention
        angle = -angle

        # differentiate with type
        if typ == 1:
            return (ops.Rx(angle), (q,))
        elif typ == 2:
            return (ops.Ry(angle), (q,))
        elif typ == 3:
            return (ops.Rz(angle), (q,))
        else:
            raise InvalidRotationGateError

    def visit_MST(self, node, params):
        return (ops.Measure, (params[4],))

    def visit_TROUT(self, node, params):
        warning("Trace out operations were replaced by measurement operations.")
        return (ops.Measure, (params[4],))

    def visit_ANC(self, node, params):
        warning("Marking ancilla qubits is not supported at the moment. Please remove the ancillas from your input gate sequence.")
        return ('ANC',(-1,))

    def visit_PS(self, node, params):
        warning('Post selection operations in the gate list are not supported.')
        return ('PS',(-1,))

    # dummies for debuging
    def visit_open(self, node, params): return "SHOULDNT BE USED"
    def visit_close(self, node, params): return "SHOULDNT BE USED"
    def visit_komma(self, node, params): return "SHOULDNT BE USED"

    # fall through for parenthesised and regex nodes
    def generic_visit(self, node, params):
        # differentiate parenthesised from regex
        if type(node) == RegexNode:
            # regex: return matched el
            return node.match.group(1)
        else:
            # parenthesised: return children
            return params

# ==== MAIN script code ====

# The following class is inspired by the IBMBackend of ProjectQ, which is released under the Apache License, Version 2
# (The NOTICE file of ProjectQ is provided in the folder QASMConverter)

#  Copyright 2017 ProjectQ-Framework (www.projectq.ch)
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

# Note that the following code of ProjectQ was adapted by Luca Mondada

class GetQASMBackend(BasicEngine):
    def __init__(self, *args, **kwargs):
        BasicEngine.__init__(self, *args, **kwargs)

        self._reset()

    def is_available(self, cmd):
        g = cmd.gate
        if g == ops.NOT and get_control_count(cmd) <= 1:
            return True
        if get_control_count(cmd) == 0:
            if g in (ops.T, ops.Tdag, ops.S, ops.Sdag, ops.H, ops.Y, ops.Z):
                return True
            if isinstance(g, (ops.Rx, ops.Ry, ops.Rz)):
                return True
        if g in (ops.Measure, ops.Allocate, ops.Deallocate, ops.Barrier):
            return True
        return False

    def _store(self, cmd):

        if self._clear:
            self._clear = False
            self.qasm = ""
            self._allocated_qubits = set()

        gate = cmd.gate

        if gate == ops.Allocate:
            self._allocated_qubits.add(cmd.qubits[0][0].id)
            return

        if gate == ops.Deallocate:
            return

        if gate == ops.Measure:
            assert len(cmd.qubits) == 1 and len(cmd.qubits[0]) == 1
            qb_id = cmd.qubits[0][0].id
            logical_id = None
            for t in cmd.tags:
                if isinstance(t, LogicalQubitIDTag):
                    logical_id = t.logical_qubit_id
                    break
            assert logical_id is not None
            self._measured_ids += [logical_id]
        elif gate == ops.NOT and get_control_count(cmd) == 1:
            ctrl_pos = cmd.control_qubits[0].id
            qb_pos = cmd.qubits[0][0].id
            self.qasm += "\ncx q[{}], q[{}];".format(ctrl_pos, qb_pos)
        elif gate == ops.Barrier:
            qb_pos = [qb.id for qr in cmd.qubits for qb in qr]
            self.qasm += "\nbarrier "
            qb_str = ""
            for pos in qb_pos:
                qb_str += "q[{}], ".format(pos)
            self.qasm += qb_str[:-2] + ";"
        elif isinstance(gate, (ops.Rx, ops.Ry, ops.Rz)):
            assert get_control_count(cmd) == 0
            qb_pos = cmd.qubits[0][0].id
            u_strs = {'Rx': 'u3({}, -pi/2, pi/2)', 'Ry': 'u3({}, 0, 0)',
                      'Rz': 'u1({})'}
            gate = u_strs[str(gate)[0:2]].format(gate.angle)
            self.qasm += "\n{} q[{}];".format(gate, qb_pos)
        else:
            assert get_control_count(cmd) == 0
            if str(gate) in self._gate_names:
                gate_str = self._gate_names[str(gate)]
            else:
                gate_str = str(gate).lower()

            qb_pos = cmd.qubits[0][0].id
            self.qasm += "\n{} q[{}];".format(gate_str, qb_pos)

    def _run(self):
        if self.qasm == "":
            return
        # finally: add measurements (no intermediate measurements are allowed)
        for measured_id in self._measured_ids:
            qb_loc = self.main_engine.mapper.current_mapping[measured_id]
            self.qasm += "\nmeasure q[{}] -> c[{}];".format(qb_loc,
                                                            qb_loc)

    def _reset(self):
        self._measured_ids = []
        self._clear = True

    def receive(self, command_list):
        for cmd in command_list:
            if not cmd.gate == ops.FlushGate():
                self._store(cmd)
            else:
                self._run()
                self._reset()

    """
    Mapping of gate names from our gate objects to the IBM QASM representation.
    """
    _gate_names = {str(ops.Tdag): "tdg", str(ops.Sdag): "sdg"}

def _ask(msg, default=None):
    ans = input(msg)
    while True:
        ans = ans if ans != '' else default

        if ans.lower().startswith('y'):
            return True
        elif ans.lower().startswith('n'):
            return False
        else:
            print('Please answer with yes or no: ')

def _define_args():
    args = argparse.ArgumentParser(
            description='A converter for quantum circuits from Mathematica format to OpenQASM code.'
    )
    args.add_argument(
            '--ibm',
            '-q',
            default=None,
            const='ibmq_qasm_simulator',
            metavar='dev',
            nargs='?',
            help='additionally convert to OpenQASM code that is compatible with the IBM Q Computer `dev` (default: ibmq_qasm_simulator)'
    )
    args.add_argument(
            '--use-hardware',
            '-w',
            dest='on_hardware',
            action='store_true',
            help='directly send code to IBM QE and run code on real hardware (default: false)'
    )
    args.add_argument(
            '--number-iterations',
            '-n',
            metavar='N',
            dest='n',
            default=1024,
            help='number of iterations that are run on the IBM Q Computer to obtain output measurements. Only with option `--use-hardware (-w) or --simulator (-s)` (default: 1024)',
            type=int
    )
    args.add_argument(
            '--simulator',
            '-s',
            dest='s',
            action='store_true',
            help='run code on IBM simulator (instead of real IBM hardware)`'
    )
    args.add_argument(
            'infile',
            metavar='input file',
            type=argparse.FileType('r'),
            nargs='?',
            default=sys.stdin,
            help='Name of input file in Mathematica format'
    )
    args.add_argument(
            'outfile',
            metavar='output file',
            type=argparse.FileType('w'),
            nargs='?',
            default=sys.stdout,
            help='Name of output file in OpenQASM format'
    )

    return args

def _check(args):
    if args.on_hardware and args.s:
        print('You cannot use flags --use-hardware (-w) and --simulator (-s) at the same time')
        msg = 'For faster results, you can run your code on the IBM simulator. Do you still want to run the code on real hardware? (y/N): '
        if _ask(msg, default='n'):
            args.s = False
        else:
            args.on_hardware = False

    if (args.on_hardware or args.s) and not args.ibm:
        print('You cannot send your code onto IBM hardware without converting your code into IBM compatible code.')
        msg = 'Do you want to convert your code into IBM compatible code? (Y/n): '
        if _ask(msg, default='y'):
            args.ibm = True
            return True
        else:
            sys.exit(1)

def main():

    args = _define_args()
    
    args = args.parse_args()

    _check(args)

    fin = args.infile
    fout = args.outfile

    # by default use IBMQ simulator to generate engines
    dev = 'ibmq_qasm_simulator'

    # remove the CNot flipper
    if not args.ibm:
        engines = engines[:-2] + engines[-1:]
        print('Converting...')
    else:
        print('Converting for IBM Q Experience...')
        dev = ibm.args

    engines = ibm.get_engine_list(device=dev)

    if args.on_hardware:
        print('We will run the converted code on the IBM device {}\n.'.format(args.ibm))
        backend = IBMBackend(use_hardware=True, device=args.ibm, num_runs=args.n)
    elif args.s:
        print('We will run the converted code on the IBM simulator\n.')
        backend = IBMBackend(use_hardware=False, device=args.ibm, num_runs=args.n)
    else:
        # backend similar to IBMBackend creating qasm code but without running it
        backend = GetQASMBackend()

    eng = MainEngine(backend=backend, engine_list=engines)

    s = fin.read()
    fin.close()

    # parse input
    root = grammar.parse(s)
    # create visitor that transforms tree into gate list
    visitor = ToGateListVisitor()
    # pass root of parsed tree to visitor
    gate_list = visitor.visit(root)
    #gate_list=flatten_touples(gate_list)

    # create circuit from gate_list
    nb_qubits = build_circuit(gate_list, eng)
    # create qasm code
    qasm      = get_qasm(eng, nb_qubits)

    if fout == sys.stdout:
        print(5*'=' + ' OpenQASM code ' + 5*'=')

    fout.write(qasm)

    if fout != sys.stdout:
        print('OpenQASM output saved to {}'.format(fout.name))
    print()

    if hasattr(backend, 'get_probabilities'):
        print(5*'=' + ' Results from {} run on hardware '.format(args.n) + 5*'=')
        for k,v in sorted(backend.get_probabilities(qureg).items()):
            print('Probability for state |{}>: {}'.format(k,v))

if __name__ == "__main__":
    
    main()

