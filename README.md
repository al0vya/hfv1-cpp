1D adaptive finite volume solver, incorporating Haar wavelet adaptivity, for the shallow equations, Harten, Lax and van Leer approximate Riemann solver

6 test cases available.

User inputs:

- Number of mother elements on the adaptive mesh
- Maximum refinement level dictating number of sub-elements per mother element
- Error threshold to control the local resolution of the adaptive mesh; smaller values include finer sub-elements on the adaptive mesh: recommended value is 10e-3 (Kesserwani et al. 2019)


