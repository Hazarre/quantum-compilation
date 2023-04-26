import numpy as np
from utils import PAULI_X, is_unitary, IDENTITY_2x2

class Gate2: # 2x2 single qubit gates
    """Represents gate acting on one qubit.
    Definitions:
    Ry(a) = exp(0.5*i*a*sigma_y)
    Rz(a) = exp(0.5*i*a*sigma_z)
    R1(a) = diag(1, exp(i*a))
    T: pi-8th gate 
    H: Hadamard gate
    """

    def __init__(self, name, arg=None):
        assert name in ['Ry', 'Rz', 'R1', 'X', 'T', 'H']
        self.name = name
        self.arg = arg

    def to_matrix(self):
        if self.name == 'Ry':
            return np.array([[np.cos(self.arg / 2), np.sin(self.arg / 2)],
                             [-np.sin(self.arg / 2), np.cos(self.arg / 2)]])
        elif self.name == 'Rz':
            return np.diag([np.exp(0.5j * self.arg), np.exp(-0.5j * self.arg)])
        elif self.name == 'R1':
            return np.diag([1.0, np.exp(1j * self.arg)])
        elif self.name == 'X':
            return PAULI_X
        elif self.name == 'T':
            return np.array([[1, 0],
                             [0, np.exp(-1j*np.pi/4) ]])
        elif self.name == 'H':
            return np.array([[1, 1],
                             [1, -1]])/np.sqrt(2)

    def is_identity(self):
        return np.allclose(self.to_matrix(), np.eye(2))

    def __repr__(self):
        if self.arg is not None:
            return self.name + "(" + str(self.arg) + ")"
        else:
            return self.name

class Gate:
    """Represents gate acting on register of qubits."""
    pass

class GateSingle(Gate):
    """Represents gate acting on a single qubit in a register."""
    def __init__(self, gate2, qubit_id):
        self.gate2 = gate2
        self.qubit_id = qubit_id

    def to_matrix(self, qubits_count):
        """Tensor product I x I x ... x `gate2.to_matrix()` x I x ... x I."""
        matrix = self.gate2.to_matrix()
        tile_size = 2 ** (self.qubit_id + 1)
        ts2 = tile_size // 2  # Half tile size.

        if (self.qubit_id == 0):
            tile = matrix
        else:
            tile = np.zeros((tile_size, tile_size), dtype=np.complex128)
            subtile = np.eye(tile_size // 2)
            for i in range(2):
                for j in range(2):
                    tile[i * ts2:(i + 1) * ts2,
                         j * ts2:(j + 1) * ts2] = subtile * matrix[i, j]

        matrix_size = 2 ** qubits_count
        ret = np.zeros((matrix_size, matrix_size), dtype=np.complex128)
        for i in range(2 ** (qubits_count - self.qubit_id - 1)):
            ret[i * tile_size:(i + 1) * tile_size,
                i * tile_size:(i + 1) * tile_size] = tile

        return ret

    def __repr__(self):
        return str(self.gate2) + " on bit " + str(self.qubit_id)

    def type(self):
        return self.gate2.name + "-single"

class GateFC(Gate):
    """ Represents fully contolled gate."""

    def __init__(self, gate2, qubit_id):
        self.gate2 = gate2
        self.qubit_id = qubit_id

    def to_matrix(self, qubits_count):
        matrix_size = 2 ** qubits_count
        index2 = (matrix_size - 1)
        index1 = index2 - 2 ** self.qubit_id
        matrix = TwoLevelUnitary(
            self.gate2.to_matrix(),
            matrix_size,
            index1,
            index2)
        return matrix.get_full_matrix()

    def __repr__(self):
        return "%s on bit %d, fully controlled" % (
            str(self.gate2), self.qubit_id)

    def type(self):
        return self.gate2.name + "-FC"

def gates_to_matrix(gates, qubits_count):
    """Converts gate sequence to matrix implemented by this sequence."""
    result = np.eye(2 ** qubits_count)
    for gate in gates:
        assert isinstance(gate, Gate)
        result = gate.to_matrix(qubits_count) @ result
    return result

def apply_on_qubit(gates, qubit_id):
    """Converts Gate2 gates to GateSingle gates acting on the same qubit."""
    return [GateSingle(gate, qubit_id) for gate in gates]

class TwoLevelUnitary:
    """Represents two-level unitary matrix.

    Two-level unitary matrix is a unitary matrix obtained from the identity
    matrix by changing a 2x2 principal submatrix.
    """

    def __init__(self, matrix2x2, matrix_size, index1, index2):
        assert index1 != index2
        assert index1 < matrix_size and index2 < matrix_size
        assert matrix2x2.shape == (2, 2)
        assert is_unitary(matrix2x2)

        self.matrix_size = matrix_size
        self.index1 = index1
        self.index2 = index2
        self.matrix_2x2 = matrix2x2
        self.order_indices()

    def __repr__(self):
        self.order_indices()
        return "%s on (%d, %d)" % (
            str(self.matrix_2x2), self.index1, self.index2)

    def order_indices(self):
        if self.index1 > self.index2:
            self.index1, self.index2 = self.index2, self.index1
            self.matrix_2x2 = PAULI_X @ self.matrix_2x2 @ PAULI_X

    def get_full_matrix(self):
        matrix_full = np.eye(self.matrix_size, dtype=np.complex128)
        matrix_full[self.index1, self.index1] = self.matrix_2x2[0, 0]
        matrix_full[self.index1, self.index2] = self.matrix_2x2[0, 1]
        matrix_full[self.index2, self.index1] = self.matrix_2x2[1, 0]
        matrix_full[self.index2, self.index2] = self.matrix_2x2[1, 1]
        return matrix_full

    def multiply_right(self, A):
        """M.multiply_right(A) is equivalent to A = A @ M.get_full_matrix()."""
        idx = (self.index1, self.index2)
        A[:, idx] = A[:, idx] @ self.matrix_2x2

    def inv(self):
        return TwoLevelUnitary(self.matrix_2x2.conj().T,
                               self.matrix_size,
                               self.index1,
                               self.index2)

    def apply_permutation(self, perm):
        assert(len(perm) == self.matrix_size)
        self.index1 = perm[self.index1]
        self.index2 = perm[self.index2]

    def is_identity(self):
        return np.allclose(self.matrix_2x2, IDENTITY_2x2)
