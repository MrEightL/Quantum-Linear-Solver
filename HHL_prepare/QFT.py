# 这是一个示例 Python 脚本。

# 按 Shift+F10 执行或将其替换为您的代码。
# 按 双击 Shift 在所有地方搜索类、文件、工具窗口、操作和设置。
import numpy as np
from qiskit import QuantumCircuit
import matplotlib.pyplot as plt
from qiskit import Aer
from qiskit.quantum_info import Statevector
from qiskit import *
from qiskit.visualization import plot_histogram
from qiskit.quantum_info import Statevector
#创建QFT函数
def QFT(n_qubits, number):
    # 建立n个量子比特的电路
    circ = QuantumCircuit(n_qubits)
    # 初始状态形成叠加
    for x in range(n_qubits):
        circ.h(x)
    # 对于初始状态进行旋转，number即为最后输出结果，可以理解为改变基底
    for i in range(n_qubits):
        circ.p(number*np.pi/2**i, i)
    # 构建QFT电路，由于qiskit是|jn……j0>,因此反过来
    for j in range(n_qubits):
        circ.h(j)
        for m in range(n_qubits-1, j, -1):
            circ.cp(np.pi/2**(m-j), m, j)
    return circ

qc = QFT(4, 12)
qc.measure_all()
qc.draw(output='mpl')
plt.show()
# 在statevector后端编译器上运行量子电路
backend = Aer.get_backend('statevector_simulator')
result = backend.run(transpile(qc, backend), shots=1000).result()
counts = result.get_counts(qc)
print(counts)
plot_histogram(counts)
plt.show()



