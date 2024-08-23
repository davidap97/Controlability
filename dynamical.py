import sympy as sp


from sympy import Symbol, expand, simplify, trace
from sympy import Matrix as mx
from sympy.physics.quantum import TensorProduct as tp
from sympy.physics.quantum import tensor_product_simp as tp_simp

from sympy import *
from sympy import I, Float

#from sympy.physics.quantum import Commutator as com
from sympy.physics.quantum import Dagger as dag



threshold = 1.0e-15


# commutator
def km(A,B):
    return simplify(expand((A*B - B*A)))

# inner product
def scalar_product(A,B):
    Ad = dag(A)
    return trace(Ad*B)
    
# projection function
def proj(B,F):
    B = expand(B)
    B = simplify(B)
    F = expand(F)
    F = simplify(F)
    if F.norm()== 0.0:
        return F
    else:
        a = scalar_product(F,B)
        b = scalar_product(F,F)
        c = (a)/(b)
        return c*simplify(expand(F))
        
# Gramschmidt orthogonalization with input matrix C against matrices from set S
def GramSchmidt(C,S):
    Cgs = C
    for i,Si in enumerate(S):
        Cgs = Cgs - proj(Cgs,Si)
        Cgs = expand(Cgs)
        Cgs = simplify(Cgs)
    return Cgs

def get_new_elements(As_k):
    As_init = As_k[0]
    As_late = As_k[-1]
    
    # compute all non-zero commutators
    komms = []
    for i in range(len(As_init)):
        for l in range(len(As_late)):
            #print(i,l)
            km_As_il = km(As_init[i], As_late[l])
            km_As_il = expand(km_As_il)
            km_As_il = simplify(km_As_il)
            #print(km_As_il)
            if not (km_As_il).norm() < threshold:
                komms.append(km_As_il)
    return komms

# The objective is to filter all matrices that are linearly dependent
def get_lin_indep_comms(komms):
    lin_indep_comms = []
    for C in komms:
        C_gs = GramSchmidt(C, lin_indep_comms)
        C_gs = simplify(expand(C_gs))
        C_gs = simplify(C_gs)
        if not C_gs.norm() < threshold:
            lin_indep_comms.append(C_gs)
            
    return lin_indep_comms

# orthogonalisiere gegen die Basis
def get_new_elements_update(lin_indep_comms, As_k):
    new_As = []
    for i,C in enumerate(lin_indep_comms):
        C_gs = C
        for k in range(len(As_k)):
            #orthogonalize against As_k[k]
            #print(i,k)
           #print(C)
            C_gs = GramSchmidt(C_gs, As_k[k])
            #print(C_gs)
            C_gs = simplify(expand(C_gs))
        if not (C_gs).norm() < threshold:
            new_As.append(C_gs)

    return new_As

def get_dim_Lie(As_k):
    dim_Lie = 0
    for k in range(len(As_k)):
        print("depth", k, "dim", len(As_k[k]))
        dim_Lie = dim_Lie + len(As_k[k])
    return dim_Lie

def dynamicalLie(As):
    As_k = []
    As_k.append(As)
    k = 1
    new_elements = True
    while new_elements:
        new_lin_indep_elements = get_lin_indep_comms(
                        get_new_elements(As_k[0:k])
                    )
        new_lin_indep_elements = get_new_elements_update(
                    new_lin_indep_elements, As_k
                )
        new_lin_indep_elements = get_new_elements_update(
                    new_lin_indep_elements, As_k
                )
        new_lin_indep_elements = get_new_elements_update(
                    new_lin_indep_elements, As_k
                )
        As_k.append(
            get_lin_indep_comms(new_lin_indep_elements)
        )
        if len(As_k[k]) == 0 or k == 100:
            new_elements = False
            dim_Lie = get_dim_Lie(As_k)
            return dim_Lie
        else:
            print(k,len(As_k[k]))
            #pprint(As_k[k])
            for elem in As_k[k]:
                print(float(elem.norm()))
            k = k+1
            
            
            
# Pauli matrices
Me = mx([[1,0],[0,1]])
Mx = mx([[0,1],[1,0]])
My = mx([[0,-I],[I,0]])
Mz = mx([[1,0],[0,-1]])



w1 = 500
w2 = 530
w3 = 550


q1 = 11
q 2= 22
q3 = 31

# Define Two-qubit Hamiltonian
H00 = w1*tp(Mz,Me) + w2*tp(Me,Mz) + q1*tp(Mx,Mx) + q2*tp(Mz,Mz)

# Define Control
H01 = tp(Me,Mx)

# Define Basis
Bas = [H00,H01]

# Calculate Dynamical Lie Algebra
print(dynamicalLie(Bas))

