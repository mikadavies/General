# Python file / module (TBD) for the Continuum Mechanics part of Mechanics of Materials (MATS23101) at the University of Manchester.

import numpy as np
import math
from sympy import Symbol, diff

alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ" #allows 52 rank tensors?

class Tensor():
    def __init__(self, rank=2, max_index=3, values=[0]*9):
        self.rank = rank
        self.max = max_index
        self.tensor = self.create(values)
    
    def __getitem__(self, item):
        return self.tensor[item] # delegate to self.__getitem__
    
    def _nest(self, layers:int, length:int) -> list:
        """
        Nests a certain amount of lists (length) in a certain amount of layers (layers) 
        """
        return [[] if layers == 0 else self._nest(layers-1, length) for _ in range(length)]
    
    def create(self, values:list) -> list:
        """
        Creates a readable and usable tensor representation given a list of values
        """
        if self.rank == 1:
            self.tensor = values
        else:
            temp_arr = np.array_split(values, len(values)/self.max) #creates the dimensions
            temp_arr = [list(x) for x in temp_arr]
            ##print(temp_arr)
            shell = self._nest(self.rank, self.max) #empty tensor
            tensor = shell
            #to make an self-changing code dependant on the tensor rank
            #I used strings, as they can be modified at will
            access_code = "tensor"
            for _ in range(0, self.rank-1):
                access_code = access_code + "[" + alphabet[_] + "]" #determines which dimension to access
            access_command = "for a in range(0, self.max):"
            for _ in range(1, self.rank-1): 
                access_command = access_command +"\n" + "    "*_ +"for " + alphabet[_] + " in range(0, self.max):" #iterates through the dimensions regardless of nesting
            cmd = "it=0\n"+ access_command  + "\n" + "    "*(self.rank-1) + access_code + " = temp_arr[it]\n"+"    "*(self.rank-1)+"it+=1"
            ##print(cmd)
            exec(cmd)
            self.tensor = tensor
        return self.tensor
    
    def _rotation_angle_2D(self): #works only in 2 dimensions
        """
            Calculates the rotation angle of a 2D 2nd order tensor
        """
        return np.arctan((2*self.tensor[0][1])/(self.tensor[0][0]-self.tensor[1][1]))/2
    
    def _rotation_matrix_2D(self): # works only in 2 dimensions
        """
            Determines the 2D rotation matrix to reach the principal tensor
        """
        a = self._rotation_angle_2D()
        return [[np.cos(a), np.sin(a)], [-np.sin(a), np.cos(a)]]
    
    def rotate(self, R): # works in all dimensions
        """
            Rotates the tensor through a rotation matrix R
        """
        P = []
        for i in range(self.max):
            r = []
            for j in range(self.max):
                c = 0
                for k in range(self.max):
                    for l in range(self.max):
                        c += R[i][k]*self[k][l]*R[j][l]
                if math.isclose(c, int(c), abs_tol=0.0000001):
                    c = int(c)
                r.append(c)
            P.append(r)
            r = []
        return P
        
    def principal_tensor_2D(self): # works only in 2 dimensions
        """
            Determines a 2D tensor's principal values
        """
        R = self._rotation_matrix_2D()
        return(self.rotate(R))
    
    def invariants(self): # maximum of 3 dimensions for rank 2 tensors
        """
            Calculates the invariants of 1st and 2nd order tensors
        """
        if self.rank == 1:
            return np.sqrt(dot(self, self))
        elif self.rank == 2:
            i1 = self[0][0] + self[1][1] + self[2][2]
            i2 = self[0][0]*self[1][1] + self[1][1]*self[2][2] + self[0][0]*self[2][2] - self[0][1]*self[1][0] - self[1][2]*self[2][1] - self[0][2]*self[2][0]
            i3 = -self[0][2]*self[1][1]*self[2][0] + self[0][1]*self[1][2]*self[2][0] + self[0][2]*self[1][0]*self[2][1] - self[0][0]*self[1][2]*self[2][1] - self[0][1]*self[1][0]*self[2][2] + self[0][0]*self[1][1]*self[2][2]
            return (i1, i2, i3)
        else: print("ERROR: Tensor.invariants() currently only handles 1st and 2nd order tensors")


class DisplacementField():
    def __init__(self, rank=2, max_index=3, values=[0]*9, variables=["x1","x2"], constants=["A"]):
        self.rank = rank
        self.max = max_index
        self.tensor = self.create(values)
        self.vars = variables
        self.consts = constants
    
    def __getitem__(self, item):
        return self.tensor[item] # delegate to self.__getitem__
    
    def _nest(self, layers:int, length:int) -> list:
        """
        Nests a certain amount of lists (length) in a certain amount of layers (layers) 
        """
        return [[] if layers == 0 else self._nest(layers-1, length) for _ in range(length)]
    
    def create(self, values:list) -> list:
        """
        Creates a readable and usable tensor representation given a list of values
        """
        if self.rank == 1:
            self.tensor = values
        else:
            temp_arr = np.array_split(values, len(values)/self.max) #creates the dimensions
            temp_arr = [list(x) for x in temp_arr]
            ##print(temp_arr)
            shell = self._nest(self.rank, self.max) #empty tensor
            tensor = shell
            #to make an self-changing code dependant on the tensor rank
            #I used strings, as they can be modified at will
            access_code = "tensor"
            for _ in range(0, self.rank-1):
                access_code = access_code + "[" + alphabet[_] + "]" #determines which dimension to access
            access_command = "for a in range(0, self.max):"
            for _ in range(1, self.rank-1): 
                access_command = access_command +"\n" + "    "*_ +"for " + alphabet[_] + " in range(0, self.max):" #iterates through the dimensions regardless of nesting
            cmd = "it=0\n"+ access_command  + "\n" + "    "*(self.rank-1) + access_code + " = temp_arr[it]\n"+"    "*(self.rank-1)+"it+=1"
            ##print(cmd)
            exec(cmd)
            self.tensor = tensor
        return self.tensor
        
    def displacement_gradient(self):
        """
            Calculates the displacement gradient for a displacement field
        """
        
        variables = [Symbol(v) for v in self.vars]
        constants = [Symbol(c) for c in self.consts]
                
        gradient = []
        
        for i in range(self.max):
            for v in variables:
                partial = diff(self.tensor[i], v)
                gradient.append(partial)
        
        gradient = Tensor(2, len(self.vars), gradient)   
        return gradient
    
    def strain_tensor(self):
        """
            Calculates the strain tensor for a given displacement field
        """
        
        variables = [Symbol(v) for v in self.vars]
        constants = [Symbol(c) for c in self.consts]
        
        gradient = self.displacement_gradient()
        epsilon = []
        for i in range(gradient.max):
            for j in range(gradient.max):
                comp = (1/2)*(gradient[i][j]+gradient[j][i])
                epsilon.append(comp)
        
        epsilon = Tensor(2, gradient.max, epsilon)
        return epsilon
    
    def rotation_tensor(self):
        """
            Calculates the rotation tensor for a given displacement field
        """
        gradient = self.displacement_gradient()
        omega = []
        for i in range(gradient.max):
            for j in range(gradient.max):
                comp = (1/2)*(gradient[i][j]-gradient[j][i])
                omega.append(comp)
        
        omega = Tensor(2, gradient.max, omega)
        return omega
    


def dot(a:Tensor, b:Tensor): #works in all dimensions
    """
        Calculates the dot product between two first order tensors of equal length
    """
    if a.rank == b.rank == 1:
        if a.max == b.max:
            c = 0
            for i in range(a.max):
                c += a[i] * b[i]
            return c
        else: print("ERROR: dot() only accepts vectors with an equal number of components")
    else: print("ERROR: dot() only accepts rank 1 tensors (vectors)")

def double_contraction(A:Tensor, B:Tensor): # works in all dimensions
    """
        Calculates the double contraction between two second order tensors of equal length
    """
    if A.rank == B.rank == 2:
        if A.max == B.max:
            c = 0
            for i in range(A.max):
                for j in range(B.max):
                    c += A[i][j] * B[i][j]
            return c
        else: print("ERROR: double_contraction() only accepts matrices with an equal number of components")
    else: print("ERROR: double_contraction() only accepts rank 2 tensors (matrices)")
    
def matrix_vector_mult(A:Tensor, b:Tensor): # works in all dimensions
    if A.rank == 2 and b.rank == 1:
        if A.max == b.max:
            product = []
            for i in range(A.max):
                c = 0
                for j in range(A.max):
                    c += A[i][j]*b[j]
                product.append(c)
                product = Tensor(1, A.max, product)
            return product
        else: print("ERROR: matrix_vector_mult() matrix and vector need to be of equal dimensions")
    else: print("ERROR: matrix_vector_mult() only accepts a matrix and vector")        
    
def strain_compatibility(e, variables, constants):
        """
            Tests if strain is valid through teh strain compatibility equation
        """
        variables = [Symbol(v) for v in variables]
        constants = [Symbol(c) for c in constants]
        
        c = 0
        if e.max != 2:
            for i in range(e.max):
                for j in range(e.max):
                    for k in range(e.max):
                        for l in range(e.max):
                            c+= diff(diff(e[i][j], variables[k]), variables[l]) + diff(diff(e[k][l], variables[i]), variables[j]) - diff(diff(e[i][l], variables[j]), variables[k]) - diff(diff(e[j][k], variables[i]), variables[l])
        else:
            c = diff(diff(e[0][0], variables[1]), variables[1])+diff(diff(e[1][1], variables[0]), variables[0])-2*diff(diff(e[0][1], variables[0]), variables[1])        
        if c == 0: 
            print("Compatibility Equation Satisfied!")
            return True
        else: 
            print("Strain Compatibility Equation Not Satisfied: " + str(c))
            return False

def internal_traction_cauchy(s:Tensor, n:Tensor):
    return matrix_vector_mult(s, n)

def normal_stress_cauchy(t:Tensor, n:Tensor):
    return dot(t, n)

def shear_stress_cauchy(t:Tensor, sN, n=0):
    if n == 0:
        return np.sqrt(dot(t, t)-sN**2)
    elif sN == 0:
        return np.sqrt(dot(t, t)-normal_stress_cauchy(t, n)**2)

def hydrostatic_stress(s:Tensor):
    h = (1/3)*s.invariants()[0]
    H = Tensor(2, 3, [0]*9)
    for i in range(H.max):
        H[i][i] = h
    return H

def deviatoric_stress(s:Tensor):
    H = hydrostatic_stress(s)
    S = Tensor(2, 3, [0]*9)           
    for i in range(S.max):
        for j in range(S.max):
            S[i][j] = s[i][j] - H[i][j]
    return S

def von_mises_stress(s, S=0):
    if S == 0:
        S = deviatoric_stress(s)
    c = 0
    for i in range(S.max):
        for j in range(S.max):
            c += S[i][j] * S[j][i]     
    return np.sqrt((3/2)*c)

def principal_normal_stresses(s:Tensor):
    # principal normal stresses using eigenvectors and eigenvalues 
    e_val, e_vec = np.linalg.eig(s)
    p3, p2, p1 = np.sort(e_val)
    return (p1, p2, p3)

def max_normal_stress(s:Tensor):
    return principal_normal_stresses(s)[0] 

def max_shear_stress(s:Tensor):
    p = principal_normal_stresses(s)
    return (1/2)*(p[0]-p[-1]) 

def normal_stress_mohr(s:Tensor, n:Tensor):
    p = principal_normal_stresses(s)
    N = 0
    for i in range(n.max):
        N += p[i]+n[i]**2  
    return N

def shear_stress_mohr(s:Tensor, n:Tensor):
    p = principal_normal_stresses(s)
    N = normal_stress_mohr(s,n)
    S = 0
    for i in range(n.max):
        S += p[i]**2+n[i]**2  
    S -= N**2
    return np.sqrt(S)

def balanced_forces(s:Tensor, variables=["x1", "x2", "x3"], constants=["A"]):
    vars = [Symbol(v) for v in variables]
    consts = [Symbol(c) for c in constants]
    tot = 0
    for i in range(s.max):
        for j in range(s.max):
              tot += diff(s[j][i], vars[j])
    if tot == 0:
        print("Forces are balanced!")
        return True
    else:
        print("Forces are not balanced: " + str(tot))
        return False
    
def balanced_moments(s:Tensor):
    c = 0
    for i in range(s.max):
        for j in range(s.max):
            if s[i][j] != s[j][i]:
                c += 1 
                print("Moments are not balanced")
                return False
    if c == 0:
        print("Moments are balanced!")
        return True
    
def stress_equilibrium(s:Tensor, variables=["x1", "x2", "x3"], constants=["A"]):
    if balanced_forces(s, variables, constants) == True and balanced_moments(s) == True:
        print("Stress field in equilibrium!")
        return True
    else:
        print("Stress field not in equilibrium")
        return False

def work_general_constitutive_model(s:Tensor, e:Tensor):
    return (1/2)*double_contraction(s, e)
