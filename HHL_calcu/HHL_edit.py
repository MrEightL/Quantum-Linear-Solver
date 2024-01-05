# Imports for Qiskit
from qiskit import QuantumCircuit, execute, Aer, IBMQ
from qiskit.compiler import transpile, assemble
from qiskit.extensions import UnitaryGate
from qiskit.tools.jupyter import *
from qiskit.circuit.library import *

from linear_solvers import NumPyLinearSolver
from qiskit.visualization import *
from qiskit import *
from qiskit import IBMQ, Aer
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.visualization import plot_histogram
from qiskit.circuit.library import QFT
import scipy as sp
from qiskit.circuit.library import QFT
import numpy as np
import math
# Various imports
import numpy as np

from copy import deepcopy
from matplotlib import pyplot as plt

#逆QFT 在Clock寄存器上，qpe的一部分
def qft_dagger(circ, q, n):
    circ.h(clock[1]);
    for j in reversed(range(n)):
        for k in reversed(range(j+1,n)):
            circ.cp(-np.pi/float(2**(k-j)), q[k], q[j]);
    circ.h(clock[0]);
    circ.swap(clock[0], clock[1]);
#QFT 在Clock寄存器上，qpe_inv的一部分
def qft(circ, q, n):
    circ.swap(clock[0], clock[1]);
    circ.h(clock[0]);
    for j in reversed(range(n)):
        for k in reversed(range(j+1,n)):
            circ.cp(np.pi/float(2**(k-j)), q[k], q[j]);
    circ.h(clock[1]);
#创建Hamilton量模拟的U门
def CGateBuild(matrix,evolution_time,label):
    Ue = QuantumCircuit(1)
    evolved = sp.linalg.expm(1j * matrix * evolution_time)
    print('Hamilton模拟矩阵为',evolved)
    Ue.unitary(evolved, [0])
    custom = Ue.to_gate(label=label).control(1)
    return custom

#QPE，在clock和b上，制备特征值
def qpe(circ, clock, target, evolution_time,na,nb,matrix):
    circuit.barrier()
    # e^{i*A*t}
    j=0
    for i in range(na):#使用cU门，哈密顿模拟
        U1=CGateBuild(matrix,(2**i)*evolution_time,'U')
        circ.append(U1,[1+j,na+nb])
        j=j+1
    circ.barrier()
    # Perform an inverse QFT on the register holding the eigenvalues
    qft_dagger(circuit, clock, na)


def inv_qpe(circ, clock, target, evolution_time, na, nb, matrix):
    # Perform a QFT on the register holding the eigenvalues
    qft(circuit, clock, na)
    circuit.barrier()
    j=na-1
    for i in reversed(range(na)):
        U1=CGateBuild(matrix, -(2**i)*evolution_time,'invU')
        circ.append(U1,[1+j, na+nb])
        j=j-1
    circuit.barrier()

def Contro_rotate(circ, na, ancilla, clock, theta,breakflag=0):
    for i in range(na):
        circuit.cry(theta[i], clock[i], ancilla)
    circ.barrier()
    circ.measure(ancilla, measurement[0])
    #with circ.if_test((measurement[0], True)):
    #    pass
    #with circ.if_test((measurement[0], False)):
     #   for i in range(na):
      #      circuit.cry(theta[i], clock[i], ancilla)


def hhl(circ, ancilla, clock, input, measurement,na,nb,matrix):
    eigenval = np.abs(np.linalg.eigvals(matrix))
    eigvec=np.linalg.eig(matrix)
    lambda_min = min(np.abs(np.linalg.eigvals(matrix)))
    lambda_max = max(np.abs(np.linalg.eigvals(matrix)))
    #esti_max=2**na-1.5
    evol_time = 1/((lambda_min * 2 ** na) / (2 * np.pi))
    #evol_time = esti_max / ((lambda_max * 2 ** na) / (2 * np.pi))
    print('t=',evol_time)
    #e_diagr=np.diag(np.array(np.real(np.e**(1j*eigenval*evol_time))))
    #e_diagi=np.diag(np.array(np.imag(np.e**(1j*eigenval*evol_time))))
    #e_diag=e_diagr+1j*e_diagi
    #print(e_diag)
    #V_trans=eigvec[1]
    #V=np.transpose(V_trans)
    #esim=V@e_diag@V_trans
    #print('哈密顿结果应为：',esim)

    qpe(circ, clock, input, evol_time,na,nb,matrix)
    circuit.barrier()
    eigenval=(eigenval * evol_time * 2 ** na) / (2 * np.pi)
    print('lambda_tude=',eigenval)
    # This section is to test and implement C = 1
    lambda_mintude = (lambda_min * evol_time * 2 ** na) / (2 * np.pi)
    C = lambda_mintude
    print('C=',C)
    theta= 2 * np.arcsin(C / eigenval)
    print('theta=', theta)
    #受控旋转
    Contro_rotate(circuit, na, ancilla, clock, theta)
    #plt.show()
    circuit.barrier()
    inv_qpe(circuit, clock, input, evol_time,na,nb,matrix)

#程序开始
#A = np.array([ [5, -1], [-1, 5] ])
A = np.array([ [1, -1/3], [-1/3, 1] ])
b = np.array([0, 1])
#A= np.array([ [19.98, -10], [-10, 19.98] ])
# A= np.array([ [20, -10], [-10, 20] ])
# b= np.array([2.8653, -0.6344])

print('A eigenvalue:', np.linalg.eigvals(A))
print('A eigenvector:', np.linalg.eig(A)[1])
b_norm=np.linalg.norm(b)
b=b/b_norm
print('b',b)
#IBMQ.save_account('Put your token')
#provider = IBMQ.load_account()
#IBMQ.get_provider(hub='ibm-q', group='open', project = 'main')
# Create the various registers needed
na = len(A)
nb = int(math.log(2 , len(b)))
clock = QuantumRegister(na, name='clock')
input = QuantumRegister(nb, name='b')
ancilla = QuantumRegister(1, name='ancilla')
measurement = ClassicalRegister(nb+1, name='c')

# Create an empty circuit with the specified registers
circuit = QuantumCircuit(ancilla, clock, input, measurement)
# State preparation.
intial_state = b
circuit.initialize(intial_state, na+1)#初始化b
circuit.barrier()

# Perform a Hadamard Transform
circuit.h(clock)#构建纠缠
hhl(circuit, ancilla, clock, input, measurement,na,nb,A)

# Perform a Hadamard Transform
circuit.h(clock)
circuit.barrier()
circuit.measure(input, measurement[1])
# circuit.draw('mpl',scale=1)
# plt.show()

# Execute the circuit using the simulator
simulator = BasicAer.get_backend('qasm_simulator')
x=np.array([0.375])
y=np.array([1.125])
for i in range(20):
    job = execute(circuit, backend=simulator, shots=5000)
#Get the result of the execution
    result = job.result()
# Get the counts, the frequency of each answer
    counts = result.get_counts(circuit)
    ans_con = counts['11']+counts['01']
    max_count=max(counts['11'],counts['01'])
    min_count=min(counts['11'],counts['01'])
    print('第',i+1,'次结果为：')
    print('概率比为：',max_count/min_count)
    classical_solution = NumPyLinearSolver().solve(A, b)
    print('经典求解为：', classical_solution.state)
    final_ans=np.sqrt(np.array([counts['01']/ans_con, counts['11']/ans_con]))
    norm=classical_solution.euclidean_norm
    ans_full=np.array(final_ans*norm/np.linalg.norm(final_ans))
    print(ans_full)
    x=np.append(x,ans_full[0])
    y=np.append(y,ans_full[1])
    print('HHL求解为：',final_ans*norm/np.linalg.norm(final_ans))
print(x)
plt.scatter(x, y, s=10)  # 绘图，s为点的大小
plt.show()
# Display the results
# plot_histogram(counts)
# plt.show()


#o_state_result = result.get_statevector(circuit, decimals=3)
#print(o_state_result)
#plot_histogram(o_state_result)

#provider.backends()

# Choose the backend on which to run the circuit
#backend = provider.get_backend('ibmq_santiago')

from qiskit.tools.monitor import job_monitor

# Execute the job
#job_exp = execute(circuit, backend=backend, shots=8192)

# Monitor the job to know where we are in the queue
#job_monitor(job_exp, interval = 2)
# Get the results from the computation
#results = job_exp.result()

# Get the statistics
#answer = results.get_counts(circuit)

# Plot the results
#plot_histogram(answer)