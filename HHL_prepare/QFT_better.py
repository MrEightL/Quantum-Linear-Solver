import numpy as np
import matplotlib.pyplot as plt
from qiskit import Aer
from numpy import pi
from qiskit import QuantumCircuit, transpile, assemble, Aer
from qiskit import transpile
from qiskit.visualization import plot_histogram
from qiskit.circuit.library import QFT
import matplotlib.pyplot as plt
def qft_rotation(circuit, n_qubits):
    if n_qubits == 0:
        return circuit
    n_qubits -= 1 #从n_qubits-1开始（0123共4个数）
    circuit.h(n_qubits)  # 添加H门
    for qubit in range(n_qubits):
        circuit.cp(pi / 2 ** (n_qubits - qubit), qubit, n_qubits)  # 添加旋转门
    qft_rotation(circuit, n_qubits)  # 递归调用本函数，依次添加低位量子比特所需要的量子门

def swap_register(circuit, num):  # 翻转末态
    for qubit in range(num // 2):
        circuit.swap(qubit, num - qubit - 1)  # 使用swap门执行翻转
    return circuit


def qft(circuit, q_num):
    qft_rotation(circuit, q_num)
    swap_register(circuit, q_num)#由于qiskit表示为|qn……q0>，因此需要翻转变回我们熟悉的|q0……qn>
    return circuit

def inverse_qft(circuit,n):
    circ = qft(QuantumCircuit(n),n)
    inv_circ = circ.inverse()
    circuit.append(inv_circ,circuit.qubits[:n])
    return circuit.decompose()

def numberX(circ,q_num,number):#直接将原来|0>转换为|+>基底下数字
    # 添加H门，形成干涉，基底从|0>->|+>
    for x in range(q_num):
        circ.h(x)
    # 对于初始状态进行旋转，number即为最后输出结果，可以理解为改变基底
    for i in range(q_num):
        circ.p(number * np.pi / 2 ** i, q_num-1-i)
    return circ

# 测试四比特QFT电路
qc = QuantumCircuit(4)
numberX(qc,4,12)
#qc.draw("mpl")
#plt.show()
#qft(qc, 4)
qc=qft(qc, 4)
qc.measure_all()
qc.draw("mpl")
plt.show()
#qc.draw("mpl")
#plt.show()
backend = Aer.get_backend('statevector_simulator')
qc_complied = transpile(qc,backend)
job = backend.run(qc_complied,shots=10000)
result = job.result()
counts = result.get_counts(qc_complied)
print(counts)
plot_histogram(counts)
plt.show()


