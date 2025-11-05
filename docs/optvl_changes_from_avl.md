# Changes in OptVL relative to AVL

Overall OptVL maintains the exact same behavior as the original AVL source codes with the following explicit exceptions.
Any other observed differences between OptVL and AVL beyond what is listed here is a bug and should be reported [here](https://github.com/joanibal/OptVL/issues).


## Mach number from geometry file is set 
AVL does not use the Mach number specified in the geometry file. I believe this is a bug and have set the initial Mach number used for an anlysis to the value in the geometry file

## X cg, Y cg, Z cg of mass file always set
For some reason if you try to implicitly load the mass file into AVL using `avl aircraft` the cg location parameters would not be updated. 
If you used the more traditional `mass aircraft.mass` and `mset 0` approach then AVL would set the cg parameters correctly.
OptVL always sets the all the data from the mass file including the cg location.

## The modal system matrix
Before writing the system matrix used for the modal analysis to screen, AVL multiplies rows and columns associates with u,w,p,r,x, and z by -1. 
I believe this is done to change the frame from the geometry-axis to the body-axis, which is usually the frame of the equations of motion.
However, the eigenvectors found correspond to the unmodified matrix in AVL and OptVL.
I have chosen to return the matrix that fits the eigenvectors returned and thus does not have the sign changes.
One can get the modified matrix by using the `in_body_axis` option when getting the matrix from OptVL like this`get_system_matrix(in_body_axis=True)

