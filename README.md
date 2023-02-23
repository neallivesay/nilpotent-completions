# Explicit upper nilpotent completions

This project contains an implementation (written in Python3 and NumPy) of the algorithm specified in the following research article:

> TBA.

## Mathematical background

This project assumes that the reader is familiar with the following facts from linear algebra:
1. Any square matrix over the complex numbers is similar to a matrix in [jcf][].
2. Any $n\times n$ [nil][] can be associated to a [par][] of $n$ (by mapping the matrix to the partition with parts corresponding to the sizes of the Jordan blocks).

[jcf]: https://en.wikipedia.org/wiki/Jordan_normal_form	"Jordan canonical form"
[nil]: https://en.wikipedia.org/wiki/Nilpotent_matrix	"nilpotent matrix"
[par]: https://en.wikipedia.org/wiki/Partition_(number_theory)	"partition"

Let $A$ be an $n\times n$ matrix with complex entries. Moreover, suppose that $A$ is nilpotent; i.e., suppose that there exists a positive integer $k$ such that $A^k=0$. There is a correspondence between nilpotent $n\times n$ matrices and partitions of $n$.

Let $r$ and $n$ be positive integers with $r<n$. Define $N_r$ to be the $n\times n$ matrix with ones on the $r$th subdiagonal and zeroes elsewhere.

Let $\lambda$ be a partition of $n$ with at most $r$ parts. By a theorem due to Mark Krupnik and Sasha Leibman [1], there exists a strictly upper triangular matrix $X$ such that $N_r+X$ is nilpotent of type $\lambda$ (i.e., the sizes of the Jordan blocks in $N_r+X$ correspond to the parts in $\lambda$).


## References

[1] M. Krupnik and A. Leibman. Jordan structures of strictly lower triangular completions of nilpotent matrices. Integr. Equat. Oper. Th. 23 (1995), 459--471.

## License

This project is licensed under the terms of the MIT License.
