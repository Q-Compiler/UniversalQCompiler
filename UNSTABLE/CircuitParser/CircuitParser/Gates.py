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

import numpy as np
from numpy import sin,cos
from functools import reduce
from copy import copy

def single_qubitify(mat, channel, n):
    """ given a unitary 2x2 matrix and the channel to apply it on
    returns the associated 2**n x 2**n matrix"""

    assert mat.shape == (2,2)
    return reduce(np.kron, [
        np.identity(2**(n-channel-1)),
        mat,
        np.identity(2**channel)
    ])

def change_channel_order(mat, frm,to):
    """returns the matrix corresponding to mat when the tensored basis is reordered: frm <-> to"""

    def swap_bits(n,i,j):
        """swaps i-th and j-th bit in n"""
        bits_are_different = int((n >> i) % 2 != (n >> j) % 2)
        return n ^ ((bits_are_different << i) | (bits_are_different << j))

    n,m = mat.shape
    return np.array(
        [[mat[swap_bits(i,frm,to), swap_bits(j,frm,to)] for j in range(m)] for i in range(n)]
    )

class Gate:
    def compute(self, vec):
        raise NotImplementedError("compute() for "+type(self)+' not defined')

    def get_matrix(self, dim):
        raise TypeError("Matrix of "+type(self)+' not defined')
    
    def inverse(self):
        raise TypeError('Inverse of '+type(self)+' not defined')

class SimpleGate(Gate):
    # uses the representation defined in the concrete gates
    def __repr__(self, fmt=''):
        cpy = copy(self)
        for attr in [k for k,v in self.parameters.items() if v is float]:
            if type(getattr(cpy, attr)) is not float:
                # parametrised template
                continue
            setattr(cpy, attr, round(getattr(cpy,attr), self.str_digit_rounding))
        return self.representation.replace('{', '{0.').format(cpy)
    
    def __init__(self, channel_offset=0, **kwargs):
        self.str_digit_rounding = 2

        for k,v in kwargs.items():
            try:
                # handle synonyms
                if hasattr(self, 'synonyms') and k in self.synonyms:
                    k = self.synonyms[k]
                # this is for one-off channel count in Mathematica
                if self.parameters[k] is int:
                    setattr(self, k, self.parameters[k](v) - channel_offset)
                else:
                    setattr(self, k, self.parameters[k](v))
            except KeyError:
                raise TypeError("'" + k + "' is an invalid keyword argument for this function")

    def compute(self, vec):
        n = int(np.log2(len(vec)))
        return np.dot(self.get_matrix(n), vec)

    @property
    def channels(self):
        try:
            return [self.channel]
        except AttributeError:
            return [self.channel1, self.channel2]

    def max_channel(self):
        return max(self.channels)

    # always keep the angle between 0 and 2pi
    @property
    def angle(self):
        return self.__angle

    @angle.setter
    def angle(self, angle):
        try:
            self.__angle = angle % (2 * np.pi)
        except TypeError:
            # parametrised angles
            self.__angle = angle
            pass

class ComposedGate(Gate, list):
    def compute(self, vec):
        for g in reversed(self):
            vec = g.compute(vec)
        return vec

    def get_matrix(self, dim):
        ret = np.identity(2**dim)
        for g in self:
            ret = np.dot(ret, g.get_matrix(dim))
        return ret
    
    def inverse(self):
        return list(map(lambda gate: gate.inverse(), reversed(self)))
        
    def copy(self):
        return ComposedGate(super(ComposedGate, self).copy())
    def of_qubit(self, qubit):
        """
        iterator only through gates concerning qubit `qubit`
        """
        for g in self:
            if isinstance(g, OneQubitGate):
                if g.channel == qubit:
                    yield g
            elif isinstance(g, TwoQubitGate):
                if g.channel1 == qubit or g.channel2 == qubit:
                    yield g
            else:
                raise NotImplementedError

    def max_channel(self):
        return max(map(SimpleGate.max_channel, self))
     
class OneQubitGate(SimpleGate):
    def get_matrix(self, dim=1):
        if dim == 1:
            return self.__matrix__()
        else:
            return single_qubitify(self.__matrix__(), self.channel, dim)

class TwoQubitGate(SimpleGate):
    def get_matrix(self, dim=2):
        """applies kronecker product to gate matrix for the n-spin situation"""
        if dim == 2:
            ret = self.__matrix__()
        else:
            ret = np.kron(np.identity(2**(dim-2)), self.__matrix__())
        # you need to be careful not to swap order twice
        if self.channel1 == 1 and self.channel2 == 0:
            ret = change_channel_order(ret, 0, 1)
        elif self.channel1 == 1:
            ret = change_channel_order(ret, 1, self.channel2)
            ret = change_channel_order(ret, 0, self.channel1)
        elif self.channel2 == 0:
            ret = change_channel_order(ret, 0, self.channel1)
            ret = change_channel_order(ret, 1, self.channel2)
        elif self.channel1 == 0 and self.channel2 == 1:
            pass
        else:
            ret = change_channel_order(ret, 0, self.channel1)
            ret = change_channel_order(ret, 1, self.channel2)
        return ret

class Rx(OneQubitGate):
    representation = 'Rx({channel})({angle})'
    parameters = {'channel': int, 'angle': float}
    
    def inverse(self):
        return Rx(channel=self.channel, angle=-self.angle)

    def __matrix__(self):
        return np.array(
            [[cos(self.angle/2), -1j*sin(self.angle/2)],
            [-1j*sin(self.angle/2), cos(self.angle/2)]])

class Ry(OneQubitGate):
    representation = 'Ry({channel})({angle})'
    parameters = {'channel':int, 'angle':float}

    def inverse(self):
        return Ry(channel=self.channel, angle=-self.angle)
        
    def __matrix__(self):
        return np.array(
            [[cos(self.angle/2), -sin(self.angle/2)],
            [sin(self.angle/2), cos(self.angle/2)]])

class Rz(OneQubitGate):
    representation = 'Rz({channel})({angle})'
    parameters = {'channel':int, 'angle':float}

    def inverse(self):
        return Rz(channel=self.channel, angle=-self.angle)
        
    def __matrix__(self):
        return np.array(
            [[np.exp(-1j*self.angle/2), 0],
            [0, np.exp(1j*self.angle/2)]])

class R(OneQubitGate):
    representation = 'R({channel})({theta},{phi})'
    parameters = {'channel':int, 'theta':float, 'phi':float}
    
    def __matrix__(self):
        return np.array(
            [[np.cos(self.theta/2), -1j * np.exp(-1j*self.phi)*np.sin(self.theta/2)],
            [-1j*np.exp(1j*self.phi)*np.sin(self.theta/2), np.cos(self.theta/2)]])

class CNot(TwoQubitGate):
    representation = 'C({channel2})({channel1})'
    parameters = {'channel1':int, 'channel2':int}
    synonyms  = {'control': 'channel2', 'target': 'channel1'}
        
    def inverse(self):
        return CNot(channel1=self.channel1, channel2=self.channel2)
        
    def __matrix__(self):
        return np.array(
            [[1,0,0,0],
             [0,1,0,0],
             [0,0,0,1],
             [0,0,1,0]])

class Swap(TwoQubitGate):
    """
    This gate swaps the two qubits.

    This is useful for some template identity. All instances of Swap are removed during CNot_simplify
    """
    representation = 'S({channel1})({channel2})'
    parameters = {'channel1': int, 'channel2':int}

    def inverse(self):
        return Swap(channel1=self.channel1, channel2=self.channel2)

    def __matrix__(self):
        return np.array(
            [[1,0,0,0],
             [0,0,1,0],
             [0,1,0,0],
             [0,0,0,1]])

class XX(TwoQubitGate):
    representation = 'XX({channel1})({channel2})({angle})'
    parameters = {'channel1': int, 'channel2': int, 'angle': float}

    def __matrix__(self):
        return np.array(
            [[cos(self.angle), 0, 0, -1j*sin(self.angle)],
             [0, cos(self.angle), -1j*sin(self.angle), 0],
             [0, -1j*sin(self.angle), cos(self.angle), 0],
             [-1j*sin(self.angle), 0, 0, cos(self.angle)]])
