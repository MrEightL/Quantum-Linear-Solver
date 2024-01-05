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
import numpy as np
import math
# Various imports
import numpy as np

from copy import deepcopy
from matplotlib import pyplot as plt
def qft_rotation(circuit, n_qubits):
    if n_qubits == 0:
        return circuit
    n_qubits -= 1 #从n_qubits-1开始（0123共4个数）
    circuit.h(n_qubits)  # 添加H门
    for qubit in range(n_qubits):
        circuit.cp(np.pi / 2 ** (n_qubits - qubit), qubit, n_qubits)  # 添加旋转门
    qft_rotation(circuit, n_qubits)  # 递归调用本函数，依次添加低位量子比特所需要的量子门

def swap_register(circuit, num):  # 翻转末态
    for qubit in range(num // 2):
        circuit.swap(qubit, num - qubit - 1)  # 使用swap门执行翻转
    return circuit


def qft(circuit, q_num):
    qft_rotation(circuit, q_num)
    swap_register(circuit, q_num)#由于qiskit表示为|qn……q0>，因此需要翻转变回我们熟悉的|q0……qn>
    return circuit

def QFT_dagger(circ, register, n):
    # for j in reversed(range(n)):
    #     for k in reversed(range(j+1,n)):
    #         circ.cp(-np.pi/float(2**(k-j)), register[k], register[j]);
    #     circ.h(register[j])
    circuit = qft(QuantumCircuit(n), n)
    inv_circ = circuit.inverse()
    circ.append(inv_circ, register[:n])
        #circ.draw('mpl')
        #plt.show()

#QFT 在Clock寄存器上，qpe_inv的一部分
def QFT(circ, register, n):
    # 构建QFT电路，由于qiskit是|jn……j0>,因此反过来
    qft_circ = qft(QuantumCircuit(n), n)
    circ.append(qft_circ, register[:n])
    # for j in reversed(range(n)):
    #     circ.h(register[j])
    #     for k in reversed(range(j+1,n)):
    #        circ.cp(np.pi/float(2**(k-j)), register[k], register[j]);


#创建Hamilton量模拟的U门
def CGateBuild(matrix,evolution_time,label):
    n=int(math.log(len(matrix), 2))
    Ue = QuantumCircuit(n)
    evolved = sp.linalg.expm(1j * matrix * evolution_time)
    print('Hamilton模拟矩阵为',evolved)
    Ue.unitary(evolved, Ue.qubits)
    custom = Ue.to_gate(label=label).control(1)
    return custom

#QPE，在clock和b上，制备特征值
def qpe(circ, clock,input,evolution_time,na,nb,matrix):
    circuit.barrier()
    # e^{i*A*t}
    j=0
    b_in = []
    for x in range(nb):
        b_in.append(input[x])
    for i in range(na):#使用cU门，哈密顿模拟
        U1=CGateBuild(matrix,(2**i)*evolution_time,'U')
        final_in = []
        final_in = [clock[j]] + b_in
        circ.append(U1, final_in)
        j=j+1
    circ.barrier()
    # Perform an inverse QFT on the register holding the eigenvalues
    QFT_dagger(circuit, clock, na)


def inv_qpe(circ, clock, input, evolution_time, na, nb, matrix):
    # Perform a QFT on the register holding the eigenvalues
    QFT(circuit, clock, na)
    circuit.barrier()
    j=na-1
    b_in=[]
    print(type(b_in))
    for x in range(nb):
        b_in.append(input[x])
    for i in reversed(range(na)):
        U1=CGateBuild(matrix, -(2**i)*evolution_time,'invU')
        final_in=[]
        final_in=[clock[j]]+b_in
        circ.append(U1,final_in)
        j=j-1
    circuit.barrier()


def hhl(circ, ancilla, clock, input, measurement,na,nb,matrix):
    eigenvalue = np.abs(np.linalg.eigvals(matrix))
    flag=0
    eigenval = []
    for i in range(len(eigenvalue)):
        if eigenvalue[i] != 0:
            eigenval.append(eigenvalue[i])
    eigenval=np.array(eigenval)
    eigvec=np.linalg.eig(matrix)
    print('非0特征值为：',eigenval)
    lambda_min = min(eigenval)
    lambda_max = max(eigenval)
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
    for i in range(na):
        circuit.cry(theta[i], clock[i], ancilla)
    circuit.barrier()
    circuit.measure(ancilla, measurement[0])
    #circuit.draw('mpl')
    #plt.show()
    circuit.barrier()
    inv_qpe(circ, clock, input, evol_time,na,nb,matrix)

#程序开始
#A = np.array([ [5, -1], [-1, 5] ])
#b = np.array([0, 1])
A = np.array([ [1, -1/3, 0, 0], [-1/3, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1] ])
b = np.array([0, 1, 0, 0])
#A= np.array([ [19.98, -10], [-10, 19.98] ])
#b= np.array([2.8653, -0.6344])
#A = np.array([ [30, -10, -10, 0],[-10, 20, 0, -10],[-10, 0, 30, -10],[0, -10, -10, 20] ])
# A = np.array([ [1, 0, 0, 0],[0, 2, 0, 0],[0, 0, 3, 0],[0, 0, 0, 4] ])
# b = np.array([1, 1, 1, -1])
print('A eigenvalue:', np.linalg.eigvals(A))
print('A eigenvector:', np.linalg.eig(A)[1])
b_norm=np.linalg.norm(b)
b=b/b_norm
print('b=',b)
#IBMQ.save_account('Put your token')
#provider = IBMQ.load_account()
#IBMQ.get_provider(hub='ibm-q', group='open', project = 'main')
# Create the various registers needed
na = len(A)
nb = int(math.log(len(b),2))
clock = QuantumRegister(na, name='clock')
input = QuantumRegister(nb, name='b')
ancilla = QuantumRegister(1, name='ancilla')
measurement = ClassicalRegister(nb+1, name='c')

# Create an empty circuit with the specified registers
circuit = QuantumCircuit(ancilla, clock, input, measurement)
# State preparation.
intial_state = b
b_input = QuantumCircuit(2)
b_input.initialize(intial_state, b_input.qubits)#初始化b
print(nb)
print([na+1,na+nb])
circuit.append(b_input,circuit.qubits[na+1:na+nb+1])
circuit.barrier()

# Perform a Hadamard Transform
circuit.h(clock)#构建纠缠
hhl(circuit, ancilla, clock, input, measurement,na,nb,A)

# Perform a Hadamard Transform
circuit.h(clock)
circuit.barrier()
meas=[]
for i in range(1, nb+1):
    meas.append(measurement[i])
circuit.measure(input, meas)
circuit.draw('mpl',scale=1)
plt.show()
# Execute the circuit using the simulator
simulator = BasicAer.get_backend('qasm_simulator')
job = execute(circuit, backend=simulator, shots=1000)
#Get the result of the execution
result = job.result()

# Get the counts, the frequency of each answer
counts = result.get_counts(circuit)
plot_histogram(counts)
plt.show()
ans_con = counts['011']+counts['001']+counts['101']+counts['111']
max_count=max(counts['011'],counts['001'])
min_count=min(counts['011'],counts['001'])
print(max_count/min_count)
classical_solution = NumPyLinearSolver().solve(A, b)
print('classical state:', classical_solution.state)
final_ans=np.sqrt(np.array([counts['001']/ans_con, counts['011']/ans_con,counts['101']/ans_con, counts['111']/ans_con]))
norm=classical_solution.euclidean_norm
print(final_ans*norm/np.linalg.norm(final_ans))
# Display the results



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