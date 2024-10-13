from FRA333_HW3_6532_6544 import endEffectorJacobianHW3 , checkSingularityHW3 , computeEffortHW3
from HW3_utils import FKHW3
import numpy as np
import math
from math import cos,sin,pi
import roboticstoolbox as rtb
from spatialmath import SE3
# file สำหรับตรวจคำตอบ
# ในกรณีที่มีการสร้าง function อื่น ๆ ให้ระบุว่า input-output คืออะไรด้วย
'''
ชื่อ_รหัส(ธนวัฒน์_6461)
1. ธีรสิทธิ บรรทัพ_6532
2. พาทัส แจ่มใจ_6544
3. -
'''
# SI unit
q = [0,0,0]
w = [5,5,9,8,4,2]

R,P,R_e,p_e = FKHW3(q)
R,P,R_e,p_e = R.round(6),P.round(6),R_e.round(6),p_e.round(6)
print("\n#===========================================<สร้าง DH โดย RBTB>====================================================#\n")
#===========================================<สร้าง DH โดย RBTB>====================================================#
T_3e = SE3(-0.47443,-0.093,0.109) * SE3.Ry(-90, unit='deg')
robot = rtb.DHRobot(
    [
        rtb.RevoluteMDH(d = 0.0892,offset=pi),
        rtb.RevoluteMDH(alpha = pi/2),
        rtb.RevoluteMDH(a = -0.425)

    ],tool=T_3e,
    name = "3DOF_Robot"
)
T_0e = robot.fkine(q)
print(f"Tranformation at frame e Robotic toolbox :\n{T_0e}\n-----------------------------------------------------\nfrom FKHW3 :\nR at frame e\n{R_e}\nP at frame e\n{p_e}")
#==============================================================================================================#
print("\n#===========================================<ตรวจคำตอบข้อ 1>====================================================#\n")
#===========================================<ตรวจคำตอบข้อ 1>====================================================#
#code here
jac = robot.jacob0(q).round(6)
print("Jacobian from roboticstoolbox :\n",jac)
print("-------------------------------------------------------\n")
print("Jacobian from HW3 :\n",endEffectorJacobianHW3(q))
#==============================================================================================================#
print("\n#===========================================<ตรวจคำตอบข้อ 2>====================================================#\n")
#===========================================<ตรวจคำตอบข้อ 2>====================================================#
#code here
epsilon = 0.001
jac_red = robot.jacob0(q).round(6)[0:3,:]
det_Jacobian = np.linalg.det(jac_red)
print(f"q1 : {q[0]}  q2 : {q[1]}  q3 : {q[2]}\ndet = {det_Jacobian}")
if abs(det_Jacobian) < epsilon:
    print("State : Singularity\n")
elif abs(det_Jacobian) >= epsilon:
    print("State : Not singularity\n")

J_e = endEffectorJacobianHW3(q)[0:3,:]
print(f"det from HW3 : {np.linalg.det(J_e)}")
print(f"Check singular from HW3 : {checkSingularityHW3(q)}")
#==============================================================================================================#
print("\n#===========================================<ตรวจคำตอบข้อ 3>====================================================#\n")
#===========================================<ตรวจคำตอบข้อ 3>====================================================#
#code here
J_e_form_Rob = robot.jacob0(q).round(6)
J_e_form_Rob_Tran = np.transpose(J_e_form_Rob)

Wforce = w[3:6]
Wmoment = w[0:3]
We = Wforce + Wmoment

Tau_func = computeEffortHW3(q,w)
Tau_check = J_e_form_Rob_Tran @ We
print(f"from func :\n{Tau_func}\n\n\nfrom Robotic toolbox :\n{Tau_check}")
#==============================================================================================================#
print("\n#==============================================================================================================#\n")