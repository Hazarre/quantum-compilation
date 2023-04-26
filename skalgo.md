# Todo
- Combine SK  
- test random unitary 
- plot seq length, runtime 
- Use Hadamard and T gate 

$$
T = \left(\begin{array}{cc}
1 & 0 \\
0 & \exp \left(\frac{i \pi}{4}\right)
\end{array}\right)
$$
- The standard gate set: Hadamard, phase, controlled-NOT and $\pi/8$ gates.
- Hadamard, phase, controlled-NOT and the Tiffoli gate.


## Preprocessing 
A general single qubit gate (2-level quantum system) can be parameterized by 3 real numbers. 
$$
U(\theta, \phi, \lambda)=\left[\begin{array}{cc}
\cos \left(\frac{\theta}{2}\right) & -e^{i \lambda} \sin \left(\frac{\theta}{2}\right) \\
e^{i \phi} \sin \left(\frac{\theta}{2}\right) & e^{i(\phi+\lambda)} \cos \left(\frac{\theta}{2}\right)
\end{array}\right]
$$
or 
$$
\begin{pmatrix}
a & b \\ -b^*e^{i\phi} & a^* e^{i\phi}
\end{pmatrix}, \text{ where } |a|^2 + |b|^2 = 1 .
$$




# Implementations 
Input: 
- error rate $\epsilon_0$ 
- target gate $U$
- universal gate set $G$ 

Output: 
- a sequence of gates 

Test: 
- compiled length
- compiled time 

Gate representation: 
- Circuit Diagram / QASM
- Strings of sequences for single qubit operations 
- Tensors / np.matrices 


Use [Qiskit Circuit Library](https://qiskit.org/documentation/stable/0.19/apidoc/circuit_library.html) to generate common gates and `Gate.to_matrix()` to obtain their numpy representation. 

Use QuTip for common gates with `Qobj().full()` to obtain their numpy represetation. 

We use the spectral/operator norm to evaluate the error between the approximation $A$ and the target matrix $U$.

$$
\epsilon = \|D\|_2=\sqrt{\lambda_{\max }\left(D^* D\right)}=\sigma_{\max }(D), \text{ where } D = A - U
$$ 

This is the Frobenius norm in [scipy.linalg.norm](https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.linalg.norm.html)


# Concept & Pseudo Code 
## Basic-Approximation 
Find an $\epsilon_0$-approx to $U \in SU(2)$ with $l_0$ sequence from the instruction set $G$. 

> Since $\epsilon_0$ is a constant, in principle this preprocessing stage may be accomplished simply by enumerating and storing a large number of instruction sequences from $G$, say up to some sufficiently large (but fixed) length $l0$, and then providing a lookup routine which, given $U$, returns the closest sequence. Appropriate values for ǫ0 and l0 will be discussed later in this section.


$O\left(\ln ^{\ln 3 / \ln (3 / 2)}(1 / \epsilon)\right)$. 

> An important practical caveat concerns the difficulty of constructing the lookup table used to obtain the basic $\epsilon_0$-approximations. This is done by enumerating all possible gate sequences up to some sufficient length. $S U(d)$ is a manifold of dimension $d^2-1$, so if we wish to approximate every gate in $S U(d)$ to within $\epsilon_0$ then we generate $O\left(1 / \epsilon_0^{d^2-1}\right)$ sequences. For an instruction set $\mathcal{G}$ there are $O\left(|\mathcal{G}|^l\right)$ sequences of length $\leq l$, some fraction of which may be redundant. We will therefore need to enumerate all sequences up to a length $l_0$ satisfying.
> $$
l_0 \geq O\left(\frac{d^2-1}{\log |\mathcal{G}|} \log \left(1 / \epsilon_0\right)\right) $$



## Group Commutator Decomposition

Suppose $U \in SU(2)$ such that $d(I,V)< \epsilon$. The objective is to find the balanced group commutator decomposition $VWV^\dagger W^\dagger=U$ such that  $d(I,V), d(I,W)< C_{gc} \sqrt{\epsilon}$ for some constant $C_{gc} > 0 $.  

Suppose now that $U$ is a rotation by an arbitrary angle $\theta$ about an arbitrary axis $\hat{p}$ on the Bloch sphere. 

Solve for $\phi$ using the equation:
$$
\sin{(\theta /2)} = 2 \sin^2{(\phi/2)} \sqrt{1-sin^4{(\phi/2)}}
$$
Define $V$ and $W$ to be rotations by $\phi$ about the $\hat{x}$ and $\hat{y}$ axis of the Bloch Sphere. 

Verify that  $C_{gc} \approx \sqrt{\epsilon/2}$. 
$$
\phi = 2\sin^{-1}([\frac{1}{2} ( 1 \pm \cos(\theta/2)  ]
^\frac{1}{4})
$$

For the given value of $\theta$ and define $V$ and $W$ to be rotations by $\phi$ about the $\hat{x}$ and $\hat{y}$ axes of the Bloch sphere. This will give the value of the group commutator decomposition.  


# Universality 
An arbitrarry unitary matrix on a $d$-dimensional Hilbert space may be written as a product of at most $2^{n-1}(2^n-1)$ two-level unitary matrices.  
Single qubit and CNOT gates can be used to implement an arbitrary two-level unitary operation on the state space of $n$ qubits.  


## Notes 
In high-dimensional spaces, the curse of dimensionality causes the algorithm to need to visit many more branches than in lower-dimensional spaces. In particular, when the number of points is only slightly higher than the number of dimensions, the algorithm is only slightly better than a linear search of all of the points. As a general rule, if the dimensionality is $k$, the number of points in the data $n$, should be $n \gg 2^k$. Otherwise, when k-d trees are used with high-dimensional data, most of the points in the tree will be evaluated and the efficiency is no better than exhaustive search. (from wikipedia)


## References: 
- [Simple example using QuTip in python by PQCLab](https://github.com/PQCLab/SolovayKitaev) 
- [Python by cryptogoth](https://github.com/cryptogoth/skc-python/blob/master/skc/kdtree.py): Well commented code but many lines 
- [in Python by ssmi1975](https://github.com/ssmi1975/solovay-kitaev-py): functional style, 1000 lines, a bit hard to parse, has nice plots.
- [Official Implementation by Dawson in C++](https://github.com/cmdawson/sk)
- [Viz and numerical analysis by sesajad](https://github.com/sesajad/solovay-kitaev-simulation)
- [Qiskit Transpiler Plugin](https://qiskit.org/documentation/stubs/qiskit.transpiler.passes.SolovayKitaevSynthesis.html)