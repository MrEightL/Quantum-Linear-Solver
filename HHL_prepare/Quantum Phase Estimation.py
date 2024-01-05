#initialization
import matplotlib.pyplot as plt
import numpy as np
import math

# importing Qiskit
from qiskit import IBMQ, Aer, transpile, assemble
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister

# import basic plot tools and circuits
from qiskit.visualization import plot_histogram
from qiskit.circuit.library import QFT
esti=4 #确定精度参数
qpe = QuantumCircuit(esti, esti-1)#有精度参数确定qubit个数，qubit-1个经典bit
qpe.x(esti-1) #本处cp(pi/4)门的特征向量为|1>，因此需要加X门来制备特征向量
qpe.draw()
#添加h门
for qubit in range(esti-1):
    qpe.h(qubit)
#QPE对每一位φ求取要施加2^n个U，从1开始
repetitions = 1
for counting_qubit in range(esti-1):
    for i in range(repetitions):
        qpe.cp(math.pi/3, counting_qubit, esti-1) # controlled-T
        #qpe.cx(counting_qubit,esti-1)
    repetitions *= 2
qpe.barrier()
# Apply inverse QFT
qpe = qpe.compose(QFT(esti-1,inverse=True), range(esti-1))
# Measure
qpe.barrier()
for n in range(esti-1):
    qpe.measure(n,n)

qpe.draw(output='mpl')
plt.show()
aer_sim = Aer.get_backend('aer_simulator')
shots = 10000
t_qpe = transpile(qpe, aer_sim)
results = aer_sim.run(t_qpe, shots=shots).result()
answer = results.get_counts()

plot_histogram(answer)
plt.show()
sorted_answer = dict(sorted(answer.items(), key=lambda item: item[1],reverse=True))
print(sorted_answer)