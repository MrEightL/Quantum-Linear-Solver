import numpy as np
from qiskit import *
import matplotlib.pyplot as plt
matrix = np.array([ [1, -1/3], [-1/3, 1] ])
vector = np.array([1/2, 1/2, 1/2, 1/2])
qc = QuantumCircuit(2) # Must redefine qc
qc.initialize(vector, [0,1]) # Initialize the 0th qubit in the state `initial_state`
qc.h(0)
qc.draw(output='mpl')
plt.show()
# 指定模拟器statevector_simulator
backend = BasicAer.get_backend('statevector_simulator')
# 执行线路
job = execute(qc, backend)
# 获取状态向量，get_statevector和statevector_simulator是对应的
qc_state = job.result().get_statevector(qc)
print(qc_state)