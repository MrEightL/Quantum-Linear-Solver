import numpy as np
from linear_solvers import NumPyLinearSolver, HHL
from qiskit.quantum_info import Statevector
from matplotlib import pyplot as plt
vol0=[ 1, 1, 1, 1, 1]
theta0=[ 0, 0, 0, 0, 0]
matrix_B1 = np.array([ [30, -10, 0, 0],[-10, 20, 0, -10],[0, 0, 30, -10],[0, -10, -10, 20] ])
matrix_B2 = np.array([[30, -10], [-10, 20]])
Y=np.array([ [-40, 20, 0, 20, 0], [ 20, -30, 10, 0, 0],[0, 10, -20, 0, 10],[20, 0, 0, -30, 10],[0, 0, 10, 10, -20] ]) * 1j
print('导纳矩阵为',Y)
#matrix_B = np.array([ [2, 0, 0, 0],[0, 4, 0, 0],[0, 0, 6, 0],[0, 0, 0, 5] ])
vector_P = np.array([-1, -1, -1, 1])
vector_Q = np.array([-0.5,-0.5])
#print(np.linalg.norm(vector_P))
classical_solution1 = NumPyLinearSolver().solve(matrix_B1, vector_P)
classical_solution2 = NumPyLinearSolver().solve(matrix_B2, vector_Q)
print('初始vol', vol0)
print('初始θ：', theta0)
#print('naive state:')
#print(naive_hhl_solution.state)
def get_solution_vector2(solution):
    """Extracts and normalizes simulated state vector
    from LinearSolverResult."""
    solution_vector = Statevector(solution.state).data[64:68].real
    #print('naive solution vector:', solution_vector)
    norm = solution.euclidean_norm
    return norm * solution_vector / np.linalg.norm(solution_vector)
def get_solution_vector1(solution):
    """Extracts and normalizes simulated state vector
    from LinearSolverResult."""
    solution_vector = Statevector(solution.state).data[16:18].real
    #print('naive solution vector:', solution_vector)
    norm = solution.euclidean_norm
    return norm * solution_vector / np.linalg.norm(solution_vector)
#设定迭代初始值
P0=vector_P#P初值
Q0=vector_Q#Q初值
delta_P=P0#ΔP初值
delta_Q=Q0#ΔQ初值
P_m1=P0#上一次迭代的初值，用来求解ΔP
Q_m1=Q0
for i in range(10):
    # deltaP、Q为去除平衡节点1后的P
    # deltaP、Q归一化
    P_iter = delta_P/vol0[1:5]
    P_num = np.linalg.norm(P_iter)
    P_iter = P_iter / P_num
    Q_iter = delta_Q/vol0[1:3]
    Q_num = np.linalg.norm(Q_iter)
    Q_iter = Q_iter / Q_num
    #B'和B"矩阵是去除平衡节点1/平衡节点1&PV节点后的Y矩阵
    B1 = matrix_B1/ P_num
    B2 = matrix_B2 / Q_num
    # print(B1)
    # print('特征值为', np.linalg.eigvals(B1))
    #开始潮流计算
    d_theta = HHL().solve(B1, P_iter)#计算V*Δθ
    d_val = HHL().solve(B2, Q_iter)#计算ΔU
    d_val0 = get_solution_vector1(d_val)
    d_theta0 = get_solution_vector2(d_theta)
    # d_theta = NumPyLinearSolver().solve(B1, P_iter)
    # d_val = NumPyLinearSolver().solve(B2, Q_iter)
    # d_val0 = np.array(d_val.state)
    # d_theta0 = np.array(d_theta.state)
    vol0[1:3] = vol0[1:3]+np.array(d_val0)#更新节点电压幅值
    theta0[1:5] = theta0[1:5]+(np.array(d_theta0)/vol0[1:5] )#更新节点相角
    V_exp = np.array(vol0)*np.exp(1j*np.array(theta0))
    #print('复数电压为',V_exp)
    I_exp = np.dot(Y, V_exp)
    #print('复数电流为', I_exp)
    P0=np.real(V_exp*np.conj(I_exp))
    #print('复数电流为', I_exp)
    Q0=np.imag(V_exp*np.conj(I_exp))
    delta_P=-P0[1:5]+vector_P
    delta_Q=-Q0[1:3]+vector_Q
    P_m1=P0[1:5]#更新上次迭代值
    Q_m1=Q0[1:3]
    theta_deg=np.array(theta0)/np.pi*180
    print('第', i+1, '次迭代电压', np.array(vol0))
    print('第', i+1,'次迭代θ：', theta_deg)

# print('full naive solution vector:',get_solution_vector1(d_val))
# print('classical state:', classical_solution2.state)