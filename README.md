# UniversalQCompiler - An opensource software package for decomposing generic quantum computations

UniversalQCompiler provides a Mathematica package that allows to decompose generic quantum computations into sequences of C-NOT gates and arbitrary single-qubit rotations. In particular, the package allows for the following:

*  Decomposing  isometries from *m* to *n* qubits, and hence, in particular decomposing arbitrary unitaries on *m* qubits and to do state preparation
*  Decomposing  quantum channels from *m* to *n* qubits
*  Decomposing POVMs on *m* qubits
*  Simplifying gate sequences by merging single-qubit rotations and cancelling C-NOTs that are next to each other
*  Drawing quantum circuits within Mathematica
*  Exporting quantum circuits to LATEX
*  Running quantum circuits on the IBM Q Experience (using the OpenQASM converter)

A detailed documentation of the Mathematica package can be found on our [webpage](http://www-users.york.ac.uk/~rc973/UniversalQCompiler.html). The notebook Examples.nb helps the user to get started quickly and provides a short overview over the methods provided by UniversalQCompiler.

This project also contains a converter to translate gate sequences from the Mathematica package UniversalQCompiler to gate sequences in Python using the QASM language of ProjectQ.

Moreover, we provide bindings to directly link Python to Mathematica. This is however provided without any guarantees or support in the directory `UNSTABLE`

## Getting started

To use our Mathematica package UniversalQCompiler.m, you need to have installed Wolfram Mathematica (we tested the package for Mathematica 11.1 and 11.3). The code relies on the package QI.m which can be downloaded from https://github.com/rogercolbeck/QI. The package can then be loaded in any Mathematica notebook (see our [documentation](https://docs.google.com/forms/d/e/1FAIpQLSc_QmF_qFwp25f8fsrVWiMKGkKbtPZeSNbOWLFU357tLpKNVw/viewform) for more details).

Note that the Makefile is only used for the unstable bindings to Python. You do not need any Makefile for the Mathematica package or the QASM converter.

## Please cite
An overview over the package UniversalQCompiler and some theoretical background about the decomposition methods that it uses can be found on arXiv:1904.01072.

## References

The code is mainly based on the following papers:

* Raban Iten, Roger Colbeck, Ivan Kukuljan, Jonathan Home, Matthias Christandl, Phys. Rev. A 93, 032318 (2016).
* M. Plesch and Č. Brukner, Phys. Rev. A 83, 032302 (2011).
* O. Giraud, M. Žnidarič, and B. Georgeot, Phys. Rev. A 80, 042309 (2009).
* V. V. Shende, S. S. Bullock, and I. L. Markov, IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems 25, 1000 (2006).
* V. Bergholm, J. J. Vartiainen, M. Möttönen, and M. M.Salomaa, Phys. Rev. A 71, 052330 (2005).
* V. V. Shende, S. S. Bullock, and I. L. Markov, Phys. Rev. A 70, 012310 (2004).
* E. Knill, LANL report LAUR-95-2225, arXiv:quant-ph/9508006 (1995).

## Authors

The first release of UniversalQCompiler (v0.1) was developed by Raban Iten (ETH Zurich) and Roger Colbeck (University of York) from 2016-19 with contributions from Luca Mondada, Oliver Reardon-Smith, Ethan Redmond and Ravjot Singh Kohli:

* Luca Mondada added the Python interface for converting the circuits to QASM.
* Oliver Reardon-Smith added functionality to cope with small cases, wrote the first version of the visualisation code and implemented several bug fixes.
* Ethan Redmond added functionality to compute the Knill decomposition and do state preparation.
* Ravjot Singh Kohli wrote the first version of the code for decomposing isometries with the Column-by-Column decomposition.

Luca Mondada worked on the code as a semester student at the Institute for Theoretical Physics at ETH Zürich under supervision of Raban Iten.

Ravjot Singh Kohli and Ethan Redmond worked on the code while summer students in the Department of Mathematics at the University of York under supervision of Roger Colbeck.

## Mathematica to OpenQASM converter
A simple-to-use python script converts the Mathematica output gate sequences in list format (not containing diagonal gates, ancilla markers or post selection operations) into OpenQASM, the Quantum Assembly Language used among others by the IBM Q Experience. The converted outputs can be used directly as input for the IBM quantum computers. Note that this script requires Python version 3 and the projectq package.

The converter is composed of a single file, `qasm_converter.py` that can be found
in the `QASMConverter` folder.

### Installation
The script requires Python 3 to be installed as well as the packages  ``numpy``, ``projectq`` and ``parsimonious`` for parsing.
If you are using pip, you can install any of these packages with the command
```shell
pip3 install PKG_NAME
```

### Usage
First make sure the Mathematica output is numerical. Then save it to the file "data"
```shell
out = PrepareForQASM[gateList];
strm = OpenWrite["data"];
Write[strm, out];
Close[strm];
```

This file can then be transformed into QASM using the qasm\_converter.py script:
```
python qasm_converter.py data OUTPUT.qasm
```
The QASM code will be found in OUTPUT.qasm.

The script can also convert Mathematica code into assembly code (OpenQASM) that can be direclty uploaded into the IBM Quantum Experience.
To convert Mathematica code into OpenQASM code compatible with the IBM Quantum computers, use the option ``--ibm``, or ``-q``.
You can also send the code direclty to the IBM Quantum computers using the option ``--simulator`` (``-s``) to run on the IBM simulator or ``--use-hardware`` (``-w``) to run on real hardware. See the help for more information.
```shell
python qasm_converter.py -h
```

## Unstable Mathematica-Python bindings [not supported]
These bindings give the user a python interface to the methods that are written in Mathematica. This allows to seamleesly use the Mathematica code within any Python project. The code worked well for us, but we do not support any issues with it.

Please note that you need to have Mathematica installed on your machine, together with a valid license.

### Installation
Everything can be compiled easily using the Makefile. It is highly recommended to use virtualenv to keep the installed packages on this project only.

```shell
virtualenv SOME_ENV
source SOME_ENV/bin/activate
```

To link Mathematica to this code, you need to open the file MathematicaLink/setup.py and change the parameters at the top of the file (your system architecture and the version of Mathematica you are using. Double check the mathematica\_folder to make sure it corresponds to your actual path).

Once this is done, you can install the package easily with the following commands.
Note: this Makefile uses pip as python package management system. You can also install everything manually if you are not using pip.
```
make unstable
make unstable-install
```

## License
UniversalQCompiler is released under the Apache 2 license.
