# การตรวจคำตอบ

## สมาชิก
- ธีรสิทธิ บรรทัพ 65340500032
- พาทัส แจ่มใจ 65340500044

## สารบัญ
  - [MDH Robotic toolbox](#mdh-robotic-toolbox)
  - [ข้อ 1 หา Jacobian](#ข้อ-1-หา-jacobian)
    - [Proof ข้อ 1](#proof-ข้อ-1)
  - [ข้อ 2 Check singularity](#ข้อ-2-check-singularity)
    - [Proof ข้อ 2](#proof-ข้อ-2)
  - [ข้อ 3 หา effort ของแต่ละข้อต่อ](#ข้อ-3-หา-effort-ของแต่ละข้อต่อ)
    - [Proof ข้อ 3](#proof-ข้อ-3)

## MDH Robotic toolbox
สร้าง Robot ที่มีลักษณะเดียวกับหุ่นยนต์ที่โจทย์ตั้งมา และทำการเปรียบเทียบกับ HW3_utils

```python
T_3e = SE3(-0.47443,-0.093,0.109) * SE3.Ry(-90, unit='deg')
robot = rtb.DHRobot(
    [
        rtb.RevoluteMDH(d = 0.0892,offset=pi),
        rtb.RevoluteMDH(alpha = pi/2),
        rtb.RevoluteMDH(a = -0.425)

    ],tool=T_3e,
    name = "3DOF_Robot"
    f"Tranformation at frame e Robotic toolbox :\n{T_0e}"
)
```
#### เปรียบเทียบกับ FKHW3
```python
R,P,R_e,p_e = FKHW3(q)
R,P,R_e,p_e = R.round(6),P.round(6),R_e.round(6),p_e.round(6)
print(f"R at frame e\n{R_e}\nP at frame e\n{p_e}")
```
### Output `q = [0,0,0]`
```
Tranformation at frame e Robotic toolbox :
   0         0         1         0.8994
   1         0         0         0.109
   0         1         0        -0.0038
   0         0         0         1
```
```
from FKHW3 :
R at frame e
[[ 0.  0.  1.]
 [ 1. -0. -0.]
 [ 0.  1. -0.]]
P at frame e
[ 0.89943  0.109   -0.0038 ]
```


## ข้อ 1 หา Jacobian
สร้างฟังก์ชันคำนวณเมทริกซ์ Jacobian ที่ตำแหน่ง end-effector โดยรับ input เป็นค่า `q` ซึ่งเป็นค่ามุมของแต่ละข้อต่อ
```python
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
```
### Output `q = [0,0,0]`
```
 [[-0.109   -0.093   -0.093  ]
 [ 0.89943  0.       0.     ]
 [-0.      -0.89943 -0.47443]
 [ 0.       0.       0.     ]
 [-0.       1.       1.     ]
 [ 1.       0.       0.     ]]
```

## Proof ข้อ 1
จาก MDH ของ Robot ที่ได้ทำไว้ข้างต้นจึงเรียกใช้งานและใช้ function คำนวณ Jacobian ของตัว Robot เทียบกับ Base หรือ Frame 0 `jacob0()`
```python
jac = robot.jacob0(q).round(6)
print("Jacobian from roboticstoolbox :\n",jac)
```
### Output `q = [0,0,0]`
```
Jacobian from roboticstoolbox :
 [[-0.109   -0.093   -0.093  ]
 [ 0.89943  0.       0.     ]
 [-0.      -0.89943 -0.47443]
 [ 0.       0.       0.     ]
 [-0.       1.       1.     ]
 [ 1.       0.       0.     ]]
```

## ข้อ 2 Check singularity
ฟังก์ชั่นในการตรวจสอบความเป็นสภาวะ Singularity โดยรับ input เป็นค่า `q` ซึ่งเป็นค่ามุมของแต่ละข้อต่อ

หากเป็นจะ Return `True` หากไม่เป็นจะ Return `False`

โดยใช้ค่า `ε = 0.001` หาก determinant ของ Jacobian มีขนาดน้อยกว่า `ε` จะถือว่าเป็นสภาวะ Singularity
```python
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
```
### Output `q = [0,0,0]`
```
False
```

## Proof ข้อ 2
จาก MDH ของ Robot ที่ได้ทำไว้ข้างต้นจึงเรียกใช้งานและหา Jacobian ของตัว Robot เช่นเดียวกับ `Proof ข้อ 1` และทำการเปียบเทียบค่า `determinant` ของ Robot และจาก func `HW3`
```python
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
```

### Output 
`q = [0,0,0]`
```
q1 : 0  q2 : 0  q3 : 0
det = 0.03554997074999999
State : Not singularity

det from HW3 : 0.03554997075000001
Check singular from HW3 : False
```

`q = [0,0,2.9]`
```
q1 : 0  q2 : 0  q3 : 2.9
det = -0.0005711105868000003
State : Singularity

det from HW3 : -0.0005711080707850093
Check singular from HW3 : True
```

## ข้อ 3 หา effort ของแต่ละข้อต่อ
สร้างฟังก์ชั่นในการหา effort ของแต่ละข้อต่อ โดย Output จะ represent ถึงค่า tau ในแต่ละข้อต่อ

โดยรับ input เป็นค่า `q` ซึ่งเป็นค่ามุมของแต่ละข้อต่อ และค่า `w` ซึ่งเป็นค่า wrench ที่กระทำกับจุดกึ่งกลางของ Fream e ที่ได้จาก `Force Sensor รุ่น FT300` 

โดยค่าที่จาก `Force Sensor รุ่น FT300` มีลักษณะดังนี้
![inputw](pic\pic2.png)

```python
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
```
### Output `q = [0,0,0]` `w = [5,5,9,8,4,2]`
```
[11.72572  2.45714  3.30714]
```

## Proof ข้อ 3
จาก MDH ของ Robot ที่ได้ทำไว้ข้างต้นจึงเรียกใช้งาน `jacob0()` ในการหา Jacobian และทำการ Transpose และคูณกับ Matrix input `w`
```python
J_e_form_Rob = robot.jacob0(q).round(6)
J_e_form_Rob_Tran = np.transpose(J_e_form_Rob)

Wforce = w[3:6]
Wmoment = w[0:3]
We = Wforce + Wmoment

Tau_func = computeEffortHW3(q,w)
Tau_check = J_e_form_Rob_Tran @ We
print(f"from func :\n{Tau_func}\n\n\nfrom Robotic toolbox :\n{Tau_check}")
```
### Output
`q = [0,0,0]` `w = [5,5,9,8,4,2]`
```
from func :
[11.72572  2.45714  3.30714]


from Robotic toolbox :
[11.72572  2.45714  3.30714]
```

`q = [0.84,2.4,2]` `w = [5,5,9,8,4,2]`
```
from func :
[9.31191195 1.96075079 3.72190876]


from Robotic toolbox :
[9.311916 1.96075  3.72191 ]
```