from HW3_utils import FKHW3
import numpy as np
import math
from math import cos, sin
# file สำหรับเขียนคำตอบ
# ในกรณีที่มีการสร้าง function อื่น ๆ ให้ระบุว่า input-output คืออะไรด้วย
'''
ชื่อ_รหัส(ธนวัฒน์_6461)
1. ธีรสิทธิ บรรทัพ_6532
2. พาทัส แจ่มใจ_6544
3. -
'''
q = [0,0,0] # can change this
w = [0,0,0,0,0,0] # can change this
#=============================================<คำตอบข้อ 1>======================================================#
#code here
def endEffectorJacobianHW3(q:list[float])->list[float]:
    # set up computation
    num_round = 6
    R,P,R_e,p_e = FKHW3(q)

    # define variable from Robot
    p_0i = [P[:,0],P[:,1],P[:,2]]
    z = [np.array([0,0,1]),np.array([-sin(q[0]),cos(q[0]),0]),np.array([-sin(q[0]),cos(q[0]),0])]

    # creat blank list for receive Jacobian
    Jv = [None] * len(q)
    Jw = [None] * len(q)
    Jacobian = np.zeros((6,3))

    # computation
    for i in range(len(q)):
        # Jvi = zi x (0Pe - 0Pi)
        Jv[i] = np.cross(z[i],(p_e-p_0i[i]))

        # Jwi = zi
        Jw[i] = z[i]

        # J = (Jv Jw)T
        Jacobian[0:3,i] = Jv[i]
        Jacobian[3:6,i] = Jw[i]
    return Jacobian
#==============================================================================================================#
J_e = endEffectorJacobianHW3(q)
#=============================================<คำตอบข้อ 2>======================================================#
#code here
def checkSingularityHW3(q:list[float])->bool:
    # define epsilon
    epsilon = 0.001

    # reduce Jacobian
    Jac_reduced = endEffectorJacobianHW3(q)[0:3,:]

    # compute determinant of reduce Jacobian
    det_jac = np.linalg.det(Jac_reduced)

    # check singularity
    if abs(det_jac) < epsilon:
        return True
    elif abs(det_jac) >= epsilon:
        return False
#==============================================================================================================#
flag = checkSingularityHW3(q)
#=============================================<คำตอบข้อ 3>======================================================#
#code here
def computeEffortHW3(q:list[float], w:list[float])->list[float]:
    # Transpose matrix
    J_e = endEffectorJacobianHW3(q)
    J_eT = np.transpose(J_e)

    # w from input is [[moment:size3][force:size3]]
    # but w in formula is [[force:size3][moment:size3]]
    Wforce = w[3:6]
    Wmoment = w[0:3]
    We = Wforce + Wmoment

    # Tau = Jt*w
    return J_eT @ We
#==============================================================================================================#
tau = computeEffortHW3(q,w)