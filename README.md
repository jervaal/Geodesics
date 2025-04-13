"""distance between separable states and Bell states. 
and the distance for GHZ, W states and the base state for three qubits.""
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
from scipy.optimize import minimize

def state(t0, p0, l0, t1, p1, l1):
    """ función para crear el estado U(t0,p0,l0).kron(U(t1,p1,l1)).|00> + |11>)/(sqrt(2)) """
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0,1)
    # qc.x(0)
    qc.u(t0, p0, l0, 0)
    qc.u(t1, p1, l1, 1)
    return Statevector.from_instruction(qc)

def state_cepa(x0, y0, z0, x1, y1, z1):
    """Función para  crear el estado U(x0,y0,z0).kron(U(x1,y1,z1)).|00>"""
    qc = QuantumCircuit(2)
    qc.u(x0, y0, z0, 0)
    qc.u(x1, y1, z1, 1)
    return Statevector.from_instruction(qc)

def d_fs(x, y):
    """distancia flubini study """
    return np.arccos(np.real(x.inner(y)))

def cost(params):
    """Función de costo: 
    Muestra los parametro para cuando la distacia de flubini study es minima
     ¨"""
    
    
    t0, p0, l0, t1, p1, l1, x0, y0, z0, x1, y1, z1 = params
    entangled_state = state(t0, p0, l0, t1, p1, l1)
    separable_state = state_cepa(x0, y0, z0, x1, y1, z1)
    return d_fs(entangled_state, separable_state)

# Rango de parámetros
bounds = [
    (0, np.pi),      # t0
    (0, 2*np.pi),    # p0
    (0, 2*np.pi),    # l0
    (0, np.pi),      # t1
    (0, 2*np.pi),    # p1
    (0, 2*np.pi),    # l1
    (0, np.pi),      # x0
    (0, 2*np.pi),    # y0
    (0, 2*np.pi),    # z0
    (0, np.pi),      # x1
    (0, 2*np.pi),    # y1
    (0, 2*np.pi),    # z1
]

# Valores iniciales aleatorios
x0 = np.random.uniform([b[0] for b in bounds], [b[1] for b in bounds])

# Minimización
res = minimize(cost, x0, bounds=bounds, method='L-BFGS-B')

print("Distancia mínima encontrada (radianes):", res.fun)
print("Parámetros óptimos:", res.x)

# Si quieres convertir la distancia mínima a grados:
print("Distancia mínima (grados):", np.degrees(res.fun))


def state(t0, p0, l0, t1, p1, l1):
    """
    Crear el estado sx.ecr.U.U'|00>, este estado es el maximamente entrelazado con puertas nativas. 
    """
    qc = QuantumCircuit(2)
    qc.sx(0)
    qc.ecr(0, 1)
    qc.u(t0, p0, l0, 0)
    qc.u(t1, p1, l1, 1)
    return Statevector.from_instruction(qc)

# Estado separable con compuertas U en cada qubit
def state_cepa(x0, y0, z0, x1, y1, z1):
    qc = QuantumCircuit(2)
    qc.u(x0, y0, z0, 0)
    qc.u(x1, y1, z1, 1)
    return Statevector.from_instruction(qc)

x0 = np.random.uniform([b[0] for b in bounds], [b[1] for b in bounds])

# Minimización
res = minimize(cost, x0, bounds=bounds, method='L-BFGS-B')

print("Distancia mínima encontrada (radianes):", res.fun)
print("Parámetros óptimos:", res.x)

# Si quieres convertir la distancia mínima a grados:
print("Distancia mínima (grados):", np.degrees(res.fun))


(t0, p0, l0, t1, p1, l1)= res.x[6:12]
(x0, y0, z0, x1, y1, z1)= res.x[0:6]
params_00  = res.x[6:12]
theta = res.x
def gamma(k):
    """Función para encontarar el camino de estados que se encuentran en la geodesica"""
    theta =  res.fun
    vector = np.cos(k)*state_cepa(x0, y0, z0, x1, y1, z1) + np.sin(k)*(state(t0, p0, l0, t1, p1, l1) - np.cos(theta)*state_cepa(x0, y0, z0, x1, y1, z1) )/(np.sin(theta))
    
    return vector.data

    from qiskit.quantum_info import Operator,Statevector
############################################################
# Distancia de Fubini-Study
def d_fs(x, y):
    return np.arccos(np.real(x.inner(y)))

#############################################################
#  Operators

II = Operator.from_label("II")
XI = Operator.from_label("XI")
ZI = Operator.from_label("ZI")
SXI= Operator.from_label("XI").power(0.5)

IX = Operator.from_label("IX")
IZ = Operator.from_label("IZ")
ISX= Operator.from_label("IX").power(0.5)

XX = Operator.from_label("XX")
ZZ = Operator.from_label("ZZ")
ISX_ISX = (Operator.from_label("X").power(0.5)).tensor(Operator.from_label("X").power(0.5)) 

###################################################################
# Vectors 

initial = Statevector.from_label("00")

cir = QuantumCircuit(2)
cir.sx(0)
cir.ecr(0,1)

final = Statevector(cir)


###################################################################
gates_nativas = [II,XI,ZI,SXI,IX,IZ,ISX, XX, ZZ, ISX_ISX]
gates_str =     ["II","XI","ZI","SXI","IX","IZ","ISX", "XX", "ZZ", "ISX_ISX"]

uno = []
for i in gates_nativas:
    U = np.matmul(i.data, initial.data)
    uno.append(U)

"""Distacia en grados para una compuerta actuando en cada qubit |00>"""
s = 0
for i in uno:
    s += 1
    dis = d_fs(Statevector(i),final)
    print(  np.degrees( dis), " ", gates_str[s-1]  )    
############################################################################# 

# diatancia flubini staudy y su compuertas para le estado |00> y el estado finala con la ecr


dos = []
for i in gates_nativas:
    for j in uno:
        UU = np.matmul(i.data,j)   
        dos.append(UU)
            
            
gate_dos = []
for i in gates_str:
    for j in gates_str:
        gate_dos.append([i,j])      
        
l = 0
for i in dos:
    l += 1
    dis = d_fs(Statevector(i),final)
    print(  np.degrees( dis), "   ",  gate_dos[l-1])                 


tres = []
for i in gates_nativas:
    for j in dos:
        UUU = np.matmul(i.data,j)   
        tres.append(UUU)
            
            
gate_tres = []
for i in gates_str:
    for j in gates_str:
        for k in gates_str:
            gate_tres.append([i,j,k])      
        
l = 0
for i in tres:
    l += 1
    dis = d_fs(Statevector(i),final)
    print(  np.degrees( dis), "   ",  gate_tres[l-1])               



    # para el ghz

# theta(000,GHZ) = 0.7853
# theta(000, W) = 1.57079
# theta(GHZ,W) = 1.57079

GHZ = QuantumCircuit(3)
GHZ.h(0)
GHZ.cx(0,1)
GHZ.cx(1,2)

ghz = Statevector(GHZ)

initial_state = Statevector.from_label("000")

w = Statevector([0,1,1,0,1,0,0,0]/np.linalg.norm(np.array([0,1,1,0,1,0,0,0])))

dis_ghz = d_fs(ghz,initial_state)
dis_w   = d_fs(ghz,w)
print(f"la distancia en drados para el ghz es: {np.degrees(dis_ghz)},  para el estado w: {np.degrees(dis_w)}")

print(f"distancia entre los dos estados gaz y w :  {np.degrees(d_fs(w,ghz))}")
