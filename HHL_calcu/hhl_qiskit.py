import numpy as np
from linear_solvers import NumPyLinearSolver, HHL
from qiskit.quantum_info import Statevector
from linear_solvers.matrices.tridiagonal_toeplitz import TridiagonalToeplitz
from matplotlib import pyplot as plt
# matrix = np.array([ [30, -10], [-10, 20]])
vector = np.array([-0.5,-0.8])
matrix = np.array([ [30, 20], [20, 20]])
# vector = np.array([0, 1])
# matrix = np.array([ [1, -1/3], [-1/3, 1]])
eigens=np.linalg.eigvals(matrix)
print('矩阵特征值为：', eigens)
print('比值为',max(eigens)/min(eigens))

#tridi_matrix = TridiagonalToeplitz(1, 20, -10)

naive_hhl_solution = HHL().solve(matrix/(np.linalg.norm(vector)), vector/(np.linalg.norm(vector)))
classical_solution = NumPyLinearSolver().solve(matrix, vector)
#tridi_solution = HHL().solve(tridi_matrix, vector)
def get_solution_vector(solution):
    """Extracts and normalizes simulated state vector
    from LinearSolverResult."""
    solution_vector = Statevector(solution.state).data[16:18].real
    norm = solution.euclidean_norm
    return norm * solution_vector / np.linalg.norm(solution_vector)
print('HHL求解为：',get_solution_vector(naive_hhl_solution))
print('经典解为：', classical_solution.state)
#print('classical Euclidean norm:', classical_solution.euclidean_norm)
#print('naive Euclidean norm:', naive_hhl_solution.euclidean_norm)
#print('full tridi solution vector:', get_solution_vector(tridi_solution)*np.linalg.norm(vector))