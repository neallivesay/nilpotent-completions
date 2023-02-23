# Explicit upper nilpotent completions

This project contains an implementation (written in Python3 and NumPy) of the algorithm specified in the following research article:

> TBA.

## Mathematical background

Let $r$ and $n$ be positive integers with $r<n$. Define $N_r$ to be the $n\times n$ matrix with ones on the $r$th subdiagonal and zeroes elsewhere.

For example, the $5\times 5$ matrix $N_5$ is shown below:
\[
N_5 =
\begin{bmatrix}
0&0&0&0&0\\
0&0&0&0&0\\
1&0&0&0&0\\
0&1&0&0&0\\
0&0&1&0&0\\
\end{bmatrix}.
\]

Note that $N_r$ is nilpotent, i.e., there exists a positive integer $k$ such that $N_r^k=0$.

Let $\lambda$ be a partition of $n$ with at most $r$ parts. By a theorem due to Mark Krupnik and Sasha Leibman, there exists a strictly upper triangular matrix $X$ such that $N_r+X$ is nilpotent of type $\lambda$ (i.e., the sizes of the Jordan blocks in $N_r+X$ correspond to the parts in $\lambda$).





## References




## License

This project is licensed under the terms of the MIT License.
