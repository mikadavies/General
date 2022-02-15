# Continuum Mechanics Module

## Tensors
### New class - Tensor
Takes inputs tensor rank/order, number of dimensions, array of values. For example:
```
T = Tensor(2, 3, [11, 12, 13, 21, 22, 23, 31, 32, 33])
```
Creates a 2nd order tensor, T, with 3 dimensions (i, j = 1,2,3): <br>
![equation](https://latex.codecogs.com/gif.latex?T%3D%5Cbegin%7Bbmatrix%7D11%2612%2613%5C%5C21%2622%2623%5C%5C31%2632%2633%5Cend%7Bbmatrix%7D)
<br> 
The ```Tensor``` class was designed to be as similar to index notation as possible, therefore accessing any element in a tensor object can be done using the same through the ```T[i][j]``` format (or ```T[i]``` or ```T[i][j][k]``` etc...). Keep in mind Python uses 0 indexing so for i=1, j=1 in index notation, here it will be ```T[0][0]```.


#### Tensor methods
- ```Tensor.rank``` shows the tensor rank
- ```Tensor.max``` shows the number of dimensions
- ```Tensor.tensor``` displays the tensor as a list of lists, for T: [[11,12,13],[21,22,23],[31,32,33]]
- ```Tensor.create(list)```:
  - changes the tensor to use the values in ```list```
  - if the new tensor has different dimensions/rank, these should be changed using ```Tensor.max``` and ```Tensor.rank```
- ```Tensor._rotation_angle_2D()``` returns the rotation angle for a 2D tensor to become a principal tensor
- ```Tensor._rotation_matrix_2D()``` returns the rotation matrix for a 2D tensor to become a principal tensor
- ```Tensor.rotate(R)``` rotates a 2nd order tensor by the matrix ```R```
- ```Tensor.principal_tensor_2D()``` returns a 2D tensor's principal tensor
- ```Tensor.invariants()``` returns the invariants of the tensor (currently works for 1st order, and 2nd order in 3D)

### New class - DisplacementField
Takes same inputs as ```Tensor``` and additional variables and constants lists. See below:
```
vars = ["x1", "x2"] 
consts = ["A"]
u = DisplacementField(1, 2, ["A*(3*x1-x2)", "A*x1*x2**2"], vars, consts)
```
Creates a displacement field defined as: <br>
![equation](https://latex.codecogs.com/svg.image?u_i%20=%20%5Cbegin%7Bbmatrix%7D%20A%5Cleft(3x_1-x_2%5Cright)%20%5C%5C%20Ax_1x_2%5E2%20%5Cend%7Bbmatrix%7D)
<br>
where x1, x2 are variables, and A is an undefined constant

#### DisplacementField methods
- ```u.displacement_gradient()``` returns the displacement gradient tensor for the displacement field
- ```u.strain_tensor()``` returns the strain tensor for the displacement field
- ```u.rotation_tensor()``` returns the rotation tensor for the displacement field

### Functions
- ```dot(a, b)``` returns the dot product of two vectors a, b
- ```double_contraction(A, B)``` returns the double contraction of two 2nd order tensors A, B
- ```matrix_vector_mult(A, b)``` returns the matrix-vector multiplication of matrix A and vector b
- ```strain_compatibility(e)``` checks if a strain tensor e is valid through the strain compatibility equation
- ```internal_traction_cauchy(s, n)``` uses Cauchy's law to determine the internal traction from a stress s and a plane with normal n
- ```normal_stress_cauchy(t, n)``` uses Cauchy's law to determine the normal stress based on a traction t and normal n
- ```shear_stress_cauchy(t, sN, n)``` uses Cauchy's law to determine the shear stress based on a traction t and normal stress sN; in case the normal isn't known it can be set to zero and the normal n can be specified
- ```hydrostatic_stress(s)``` determines the hydrostatic stress for a stress tensor s
- ```deviatoric_stress(s)``` determines the deviatoric stress for a stress tensor s
- ```von_mises_stress(s, S)``` determines the von Mises stress from a stress tensor s, or alternatively from the deviatoric stress S (set to 0 by default); if s is unknown, it can be set to 0 and S can be specified instead
- ```principal_normal_stresses(s)``` determines the principal normal stresses for a stress tensor s using eigenvalues
- ```max_normal_stress(s)``` determines the maximum normal stress for a stress tensor s
- ```max_shear_stress(s)``` determines the maximum shear stress for a stress tensor s
- ```normal_stress_mohr(s, n)``` uses Mohr's circle to determine the normal stress for a stress tensor s and normal n
- ```shear_stress_mohr(s, n)``` uses Mohr's circle to determine the shear stress for a stress tensor s and normal n
- ```balanced_forces(s, vars, consts)``` determines whether or not the forces are balanced for a stress tensor s, whose components are described by equations with variables declared in vars and unknown constants declared in consts
- ```balanced_moments(s)``` determines whether or not the moments are balanced for a stress tensor s
- ```stress_equilibrium(s, vars, consts)``` determines whether or not a stress field s is in equilibrium
- ```work_general_constitutive_model(s, e)``` computes work based on the general constitutive model using stress s and strain e
