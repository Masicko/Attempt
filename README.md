# Attempt

The idea is to solve complex symmetric (not hermitian) banded matrix using sliding window allowing to solve the system locally, i.e. with minimal memory usage. I am not sure about a numerical stability, but it may work well. Maybe at least as a initial guess for an iterative solver(?).

## Method

### General idea
Lets suppose we have a banded complex symmetric square matrix $A$ with dimensions $n$x$n$ and bandwidth $d$, complex right-hand side $b$ and we search the solution $x$ of the equation system

```math
Ax = b.
```

Since A is symmetrical, we can apply symmetrical changes to the matrix such that 
```math

\underbrace{(LAL^T)}_{D} \underbrace{(L^{-T}x)}_{y} = \underbrace{Lb}_{\tilde{b}},
```

where $L$ is a lower triangular (banded) matrix and $D$ is a diagonal matrix. Solving the diagonal problem and obtaing $x$ is then an easy task providing
```math
\begin{align}
\tilde{b} &= Lb,
\\
y &= D^{-1} \tilde{b},
\\
x &= L^T y.
\end{align}
```

Storing the whole $L$ would take a lot of memory. However, since we have a banded matrix we can compute entries of $D$, $y$, and $x$ locally using only a part of $L$.

### Algorithm
First, let us denote a band-radius $r$ = $\frac{d-1}{2}$. We construct $L$ only locally and we will store only $r+1$ columns (of length $r+1$) at once which yields $(r+1)^2$ memory requirements. Next, we need to store a changing values of $A$ depending on $LU$ factorization, but a window $U$ of dimensions $d$ x $d$ will be sufficient. We will proceed in the following algorithm.


#### Preparation

- matrix $U$ of dimensions $d$ x $d$ and copy $U = A[1:d, 1:d]$
- vectors $l_i$ of length $r+1$ for $i \in 1 \dots r+1$ (serving as "compressed" columns of $L$)
- vector $y$ of length $r + 1$

#### Main loop
for $k$ in $1, 2, \dots n$
- $m = \text{max}(k-r, 1)$
- (we have $l_{m}$, $l_{m+1}, \dots l_{k-1}$ and $y[m : k-1]$ stored)
- compute $l_k$ while changing values in $U$
- dispose $l_m$ (if $m >= 1$) and add $l_k$ to matrix $l$ (can be done via cyclic modulo indexing)

- $\tilde{b}[k] = L[m : k, k] x$ where $L[m : k, k]$ is implicitely constructed using $l_i$ vectors
- $d$ = $U[1,1]$

- dispose $y[m]$ (if $m >= 1$) and add $y[k] = \frac{\tilde{b}[k]}{d}$


- if $k >= r+1$ compute $x[k - r] = l_k \cdot y[k:k +r]$ (because $l_k = L^T[k, k : k + r]$)

- dispose first column and row from U and add the $k+1$-th row and column from A (this can be time consuming)

#### Last entries
Compute the remainng $x[n-r+1 : n]$ as

for $i$ in $n-r+1, \dots n$

- x[i] = $L^T[i, i : m] \cdot y[i : m]$ for $m = \text{min}(i +r, n)$
