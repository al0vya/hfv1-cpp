# HFV1-CPU

## Model description

This is a successor to the base shallow water model solver, <a href="https://github.com/al0vya/FV1_cpp/blob/master/README.md">FV1-CPU</a> and is based on <a href="https://www.sciencedirect.com/science/article/pii/S0309170819301770">this paper</a>. It consists of two paired parts: the multiresolution analysis (MRA) of 'Haar' wavelets and a finite volume (FV1) scheme, and is called 'HFV1-CPU'. HFV1-CPU is an *adaptive* solver and adapts the mesh to the solution by way of MRA. For example, it uses fine *sub*-elements in areas with high localised variation in water height, discharge and topography, and coarser sub-elements otherwise. Please read the information in the link for FV1-CPU for a brief description on what a mesh is and the FV1 scheme. 

### Multiresolution analysis

MRA gives access to a set of 'details', which are compared against a user-input error threshold `epsilon` to decide whether finer or coarser sub-elements are used in the mesh. Smaller values of `epsilon` include finer sub-elements, whereas larger values of `epsilon` include coarser ones. The more fine sub-elements are used, the greater the accuracy of the solution but also the greater the computational cost: <a href="https://www.sciencedirect.com/science/article/pii/S0309170819301770">Kesserwani et al. 2019</a> recommended a value of 10<sup>-3</sup> to maintain a balance between accuracy and cost for flood modelling applications. The finest possible sub-element resolution is dictated by a user-input maximum refinement level `L`, which subdivides a *mother* element into 2<sup>L</sup> sub-elements. For example, if a single mother element is used to model a domain 100 m in length using a maximum refinement level of 9, the mesh can comprise between 1 and 2<sup>9</sup> sub-elements, and the finest possible resolution is 100 m / 2<sup>9</sup> = 0.195 m.

## Running the model

HFV1-CPU has the same 6 test cases as <a href="https://github.com/al0vya/FV1_cpp/blob/master/README.md">FV1-CPU</a>. After building and then running the executable, the user must select:

* Which test case to run
* How many mother elements are to comprise the mesh
* A value for the error threshold `epsilon`
* A value for the maximum refinement level `L`

## Speedups relative to FV1-CPU

Below are speedups of HFV1-CPU over FV1-CPU for a 40 s wet dam break simulation, run in release mode on a machine with an i5-8250U CPU.

![speedups](/HFV1_cpp/adaptivity_speedup.png)
