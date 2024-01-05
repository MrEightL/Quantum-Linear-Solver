import qiskit
import math
from math import pi
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, IBMQ, execute
from qiskit import Aer, transpile, assemble
import math
import random
import numpy as np
from scipy.optimize import minimize
from qiskit.visualization import circuit_drawer, plot_histogram
from math import sqrt
# 一个单元格有多个输出
from IPython.core.interactiveshell import InteractiveShell
import matplotlib.pyplot as plt
InteractiveShell.ast_node_interactivity = 'all'


# 定义函数fixed-ansatz,用于搭建门电路V，V|000>=|x(alpha)>
def apply_fixed_ansatz(circ, qubits, parameters, count_type):
    # count_type列表指定了每一层单旋转比特门的组成，列表parameters指定每个单比特旋转门的角度
    if ((len(qubits) == 3) or (len(qubits) == 6)):
        for iz in range(0, 3):
            # 随机数为1:旋转门x,随机数为2：旋转门y,随机数为3：旋转门z
            count_gate = count_type[0][iz]
            if count_gate == 1:
                circ.rx(parameters[0][iz], qubits[iz])
            if count_gate == 2:
                circ.ry(parameters[0][iz], qubits[iz])
            if count_gate == 3:
                circ.rz(parameters[0][iz], qubits[iz])

        circ.cz(qubits[0], qubits[1])
        circ.cz(qubits[2], qubits[0])

        for iz in range(0, 3):
            # 生成3个随机数,随机数为1:旋转门x,随机数为2：旋转门y,随机数为3：旋转门z
            count_gate = count_type[1][iz]
            if count_gate == 1:
                circ.rx(parameters[1][iz], qubits[iz])
            if count_gate == 2:
                circ.ry(parameters[1][iz], qubits[iz])
            if count_gate == 3:
                circ.rz(parameters[1][iz], qubits[iz])

        circ.cz(qubits[1], qubits[2])
        circ.cz(qubits[2], qubits[0])

        for iz in range(0, 3):
            # 随机数为1:旋转门x,随机数为2：旋转门y,随机数为3：旋转门z
            count_gate = count_type[2][iz]
            if count_gate == 1:
                circ.rx(parameters[2][iz], qubits[iz])
            if count_gate == 2:
                circ.ry(parameters[2][iz], qubits[iz])
            if count_gate == 3:
                circ.rz(parameters[2][iz], qubits[iz])

            # Creates the Hadamard test,for calculating <psi|psi>


def had_test(gate_type, qubits, auxiliary_index, parameters, count_type):
    # 创建局部变量，使得每执行一次循环，就会生成一个新的电路
    circ = QuantumCircuit(4, 4)
    circ.h(auxiliary_index)
    apply_fixed_ansatz(circ, qubits, parameters, count_type)

    for ie in range(0, len(gate_type[0])):
        if (gate_type[0][ie] == 1):
            circ.cz(auxiliary_index, qubits[ie])
    #         if (gate_type[0][ie] == 0):
    #             circ.cx(auxiliary_index, qubits[ie])
    #             circ.cx(auxiliary_index, qubits[ie])

    for ie in range(0, len(gate_type[1])):
        if (gate_type[1][ie] == 1):
            circ.cz(auxiliary_index, qubits[ie])
    #         if (gate_type[1][ie] == 0):
    #             circ.cx(auxiliary_index, qubits[ie])
    #             circ.cx(auxiliary_index, qubits[ie])

    circ.h(auxiliary_index)

    # 测量辅助比特
    circ.measure(0, 0)
    return circ


# #验证函数had_test()是否正确
# circ1 = had_test([[0,0,0], [0,0,0]],[1, 2, 3], 0, [[1,1,1],[1,1,1],[1,1,1]],[[1,3,3],[1,2,1],[2,3,1]])
# #电路可视化
# circ1.draw(output='mpl')
# simulator = Aer.get_backend('qasm_simulator')
# job = execute(circ1,simulator,shots=1000)
# result = job.result()
# # Returns counts
# counts = result.get_counts(circ1)
# # 统计辅助比特测得1的概率
# overall_sum_1=0
# if ('0001' in counts.keys()):
#     m_sum = float(counts["0001"])/1000
# else:
#     m_sum = 0

# overall_sum_1+=(1-2*m_sum)
# print('the final outcome of had_test: '+str(overall_sum_1))
# # 输出直方图
# plot_histogram(counts)


# Creates the swap test, for calculating |<b|psi>|^2

def swap_test(gate_type, qubits, auxiliary_index, parameters, count_type):
    circ = QuantumCircuit(7, 7)
    #     #制备|b>
    #     circ.x(6)
    #     circ.h(5)
    # swap test左侧H门
    circ.h(auxiliary_index)
    # 实现V门，用于制备|x(alpha)>
    apply_fixed_ansatz(circ, qubits, parameters, count_type)

    # 实现qubit1,2,3上的受控A
    for ie in range(0, len(gate_type[0])):
        if (gate_type[0][ie] == 1):
            circ.cz(auxiliary_index, qubits[ie])
    #         if (gate_type[0][ie] == 0):
    #             circ.cx(auxiliary_index, qubits[ie])
    #             circ.cx(auxiliary_index, qubits[ie])

    # 实现qubit4,5,6上的受控A
    for ie in range(0, len(gate_type[1])):
        if (gate_type[1][ie] == 1):
            circ.cz(auxiliary_index, qubits[ie + 3])
    #         if (gate_type[1][ie] == 0):
    #             circ.cx(auxiliary_index, qubits[ie+3])
    #             circ.cx(auxiliary_index, qubits[ie+3])

    # 多qubit,CSWAP的实现
    circ.cswap(0, 1, 4)
    circ.cswap(0, 2, 5)
    circ.cswap(0, 3, 6)

    # swap test辅助比特右侧H门
    circ.h(auxiliary_index)
    circ.measure(0, 0)
    return circ

# #验证函数swap_test()是否正确
# circ1 = swap_test([[0,1,0], [1,0,0]],[1, 2, 3,4,5,6], 0, [[1,1,1],[1,1,1],[1,1,1]],[[1,3,3],[1,2,1],[2,3,1]])
# #电路可视化
# circ1.draw(output='mpl')
# simulator = Aer.get_backend('qasm_simulator')
# job = execute(circ1,simulator,shots=1000)
# result = job.result()
# # Returns counts
# counts = result.get_counts(circ1)
# plot_histogram(counts)
# # 统计辅助比特测得1的概率
# overall_sum_2=0
# if ('0000001' in counts.keys()):
#     m_sum = float(counts["0000001"])/1000
# else:
#     m_sum = 0

# overall_sum_2+=(1-2*m_sum)
# print('the final outcome of swap_test: '+str(overall_sum_2))
# Implements the entire cost function on the quantum circuit (sampling, 1000 shots)


def calculate_had_test(parameters, count_type):
    global opt

    # 初始化参数，用于存储<psi|psi>的计算结果
    overall_sum_1 = 0
    # 对传过来的列表x0进行切片
    parameters = [parameters[0:3], parameters[3:6], parameters[6:9]]
    count_type = [count_type[0:3], count_type[3:6], count_type[6:9]]

    for i in range(0, len(gate_set)):
        for j in range(0, len(gate_set)):
            backend = Aer.get_backend('aer_simulator')

            # 统计两个受控门的系数乘积
            multiply = coefficient_set[i] * coefficient_set[j]

            # 接收had_test()函数返回的电路
            circ1 = had_test([gate_set[i], gate_set[j]], [1, 2, 3], 0, parameters, count_type)
            simulator = Aer.get_backend('qasm_simulator')
            job = execute(circ1, simulator, shots=100)
            result = job.result()
            # Returns counts
            outputstate = result.get_counts(circ1)

            # 统计辅助比特测得1的概率
            if ('0001' in outputstate.keys()):
                m_sum = float(outputstate["0001"]) / 100
            else:
                m_sum = 0

            # 内积
            overall_sum_1 += (abs(multiply)) * (1 - 2 * m_sum)
    #     print('calculate_had_test:   '+str(overall_sum_1))
    return overall_sum_1


# 计算内积的模方


def calculate_swap_test(parameters, count_type):
    global opt
    overall_sum_2 = 0
    parameters = [parameters[0:3], parameters[3:6], parameters[6:9]]
    count_type = [count_type[0:3], count_type[3:6], count_type[6:9]]

    for i in range(0, len(gate_set)):
        for j in range(0, len(gate_set)):
            multiply = coefficient_set[i] * coefficient_set[j]
            mult = 1
            backend = Aer.get_backend('aer_simulator')
            circ = swap_test([gate_set[i], gate_set[j]], [1, 2, 3, 4, 5, 6], 0, parameters, count_type)
            job = execute(circ, backend, shots=100)
            result = job.result()
            # Returns counts
            outputstate = result.get_counts(circ)

            # 计算辅助比特测得1的概率

            if ('0000001' in outputstate.keys()):
                m_sum = float(outputstate["0000001"]) / 100
            else:
                m_sum = 0

            mult = mult * (1 - 2 * m_sum)
            overall_sum_2 += (abs(multiply)) * mult
    #             print('calculate_swap_test:   '+str(overall_sum_2))
    return overall_sum_2


def calculate_cost_function2(parameters, count_type):
    global opt

    # 获取计算cost所需数值
    overall_sum_1 = calculate_had_test(parameters, count_type)
    overall_sum_2 = calculate_swap_test(parameters, count_type)

    # 计算代价函数
    outcome = 1 - float(overall_sum_2 / overall_sum_1)

    # 确保代价函数的数值大于等于0
    final_outcome = abs(outcome)

    # 返回代价函数值
    return final_outcome


# 主函数
coefficient_set = [0.3, 0.5, 0.2]
gate_set = [[0, 0, 0], [0, 1, 0], [1, 0, 0]]

for i in range(0, 1):
    # 随机生成9个参数，用于搭建ansatz,用列表count_type来存储生成的随机数
    count_type = []
    for iz in range(0, 9):
        # 生成3个随机数,随机数为1:旋转门x,随机数为2：旋转门y,随机数为3：旋转门z
        count_gate = random.randint(1, 3)
        count_type.append(count_gate)

    # 创建列表存储代价函数值,用于存储同一种ansatz下不同参数对应的代价函数值

    cost_function_value2 = []

    # 定义字典，存储随机生成的parameters,最小代价函数值和最小代价函数所对应的parameters构成了一对键值对

    her_parameters2 = {}

    for j in range(0, 100):
        x0 = [float(random.randint(0, 3000)) / 1000 for i in range(0, 9)]

        out2 = minimize(calculate_cost_function2, x0, count_type, method="COBYLA", options={'maxiter': 10})

        cost_function_value2.append(out2['fun'])
        # 记录代价函数值和参数之间的对应关系
        her_parameters2[out2['fun']] = out2['x']

    # 当前ansatz下代价函数的最小值，将代价函数值从小到大排序，输出最小值
    cost_function_value2.sort()

    print('VQLS代价函数的收敛效果，归一化：')
    print(cost_function_value2)

    # 代价函数最小值
    print(cost_function_value2[0])

    # 输出代价函数最小时所对应的参数
    print(her_parameters2[cost_function_value2[0]])



