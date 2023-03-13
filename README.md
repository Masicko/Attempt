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

/underbrace{(LAL^T)}{D} \underbrace{(L^{-T}x)}{y} = \underbrace{Lb}{\tilde{b}},
```

where $L$ is a lower triagonal matrix and $D$ is a diagonal matrix. Solving the diagonal problem and obtaing $x$ is then an easy task providing
```math
\begin{align}
D y &= \tilde{b}
\\
x =& L^T y.
\end{align}
```

Storing the whole $P$ would take a lot of memory. However, since we have a banded matrix we can compute entries of $D$, $y$, and $x$ locally using only a part of $L$. 

### Algorithm
First, let us denote $r$ = $\frac{d-1}{2}$. Let us denote the $k$-th column of $L$ as $l_k$ which has effectively maximum length of $r+1$ providing the banded structure of $A$. If $L_k$ denotes a identity matrix except the $k$-th column, where $l_k$ resides, we can write
```math
L = L_1 L_2 L_3 \dots L_n.
```
We will construct $L$ only locally and we will store only $r+1$ columns (of length $d$) at once yielding $d^2$ memory requirements. Next, we need to store a changing values of $A$ depending on $LU$ factorization, but a window $U$ of dimensions $d$x$d$ will be sufficient. We will proceed in the following algorithm.
//
Prepare:

- matrix $U$ of dimensions $d$ x $d$ and copy $U = A[1:d, 1:d]$
- matrix $l$ of dimensions $r+1$ x $r+1$
- array $y$ of length $r + 1$

for $k$ in $1, 2, \dots n$

- we have stored $l_{m}$, $l_{m+1}$, ... $l_{k-1}$, where $m = \text{max}(k-r, 1)$
- compute $l_k$ while changing values in $U$
- dispose $l_m$ and insert $l_k$ to matrix $l$ (rearanging columns can be time consuming)

- $d$ = $U[k]$

- $\tilde{b}[k] = L[k-r : k, k] x$ where $L[k-r : k, k]$ is implicitely constructed using $l_i$ vectors
- $y[k] = \frac{\tilde{b}[k]}{d}$
- dispose $y[k-r]$ and add $y[k]$

- if $k >= r+1$ compute $x[k - r] = L^T[k, k : k + r] y[k:k +r]$

- dispose first column and row from U and add the $k+1$-th row and column from A (this can be time consuming)

Compute the remainng $x[n-r+1 : n]$ as

for $i$ in $n-r+1, \dots n$

- x[i] = $L[i, i : m] y[i : m]$ for $m = \text{min}(i +r, n)$
