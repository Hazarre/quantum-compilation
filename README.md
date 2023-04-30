# Compiling Quantum Programs
This is the repo for Li-Heng Henry Chang's Senior project "Compiling Qauntum Programs" at Bard College, NY. 

## Objectives
The goal of this repository is to use simple code to demonstrate the implementation of the Universal Decomposition and the Solovay-Kitaev algorithm without going into the theories that require high-level math. This will allow us to build a simple compiler as the important conceptual peices used to build a quantum compiler. 

## UDSK Quantum compiler 
This repository builds a compiler with two components: 
- Universal Decomposition decomposes a quantum program $U$ into CNOT and single qubit gates (two-level) unitaries. 
- The single qubit gates are further approximated using the Solovay-Kitaev algorithm. In the compiler we use the Qiskit built-in Solovay-Kitaev algorithm since Qiskit provides a easy way to represent quantum gates and plot nice circuit diagrams. 

## Standalone ```.ipynb``` Python notebooks 
- ```UDSKcompiler.ipynb``` is the final compiler file that is used. 
- ```SKalgo.ipynb``` is a file that implements the Solovay-Kitaev thoerem from scratch.
- ```qiskitSK.ipynb``` is a example file for using the [Solovay-Kitaev algorithm inside the Qiskit library](https://qiskit.org/documentation/stubs/qiskit.transpiler.passes.SolovayKitaevSynthesis.html).

## Packages 
- ```numpy``` (1.23.3) for matrix operations
- ```scipy.stats``` for generating random unitary matrix 
- ```qiskit``` (0.42.0) for circuit diagram and SolovayKitaevSynthesis


## Todo 
1) Sort the final compiler file. 
2) Find main references 

## References 
We follow the implementation of [this repo](https://github.com/fedimser/quantum_decomp) to decompose a general unitary matrix into two-level ones. This gives us a product of two-level matrices, which is further approximated using the Solovay-Kitaev algorithm.  



