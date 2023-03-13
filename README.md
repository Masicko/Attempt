# Attempt

The idea is to solve complex symmetric (not hermitian) banded matrix using sliding window solving the system locally. I am not sure about a numerical stability, but it may work well. Maybe at least as a initial guess for an iterative solver(?).

## Method
Lets suppose we have a complex symmetric matrix $A$, complex right-hand side $b$ and we search the solution $x$ of the equation system

```math
Ax = b.
```

Since A is symmetrical, we can apply symmetrical changes to the matrix such that 
```math
(PDP^T) (P^{-T}x) = Pb,
```

where $P$ is a lower triagonal matrix and $D$ is a diagonal matrix.