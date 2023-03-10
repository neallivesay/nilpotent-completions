# Upper nilpotent matrix completions

This project contains an implementation (written in Python3 and NumPy) of the matrix completion algorithm specified by Livesay, Sage, and Nguyen in [[LSN23]](#1).

## Matrix algebra background

This project assumes that the reader is familiar with the following facts from linear algebra:
1. Any square matrix over the complex numbers is similar to a matrix in [Jordan canonical form](https://en.wikipedia.org/wiki/Jordan_normal_form).
2. The collection of Jordan block sizes of an $n\times n$ matrix forms a [partition](https://en.wikipedia.org/wiki/Partition_(number_theory)) of $n$.
3. In particular, any $n\times n$ [nilpotent matrix](https://en.wikipedia.org/wiki/Nilpotent_matrix) can be associated to a partition of $n$.

If $A$ is a nilpotent matrix and $\lambda$ is a partition, we say that $A$ has _type_ $\lambda$ if the parts of $\lambda$ are the Jordan block sizes of $A$.

Let $r$ and $n$ be positive integers such that $r$ is less than $n$. Define $N_r$ to be the $n\times n$ matrix with ones on the $r$ th subdiagonal and zeroes elsewhere. It is easy to verify that $N_r$ is nilpotent.

Let $\lambda$ be a partition of $n$ with at most $r$ parts. By a theorem due to Mark Krupnik and Sasha Leibman [[KL95]](#2), there exists a strictly upper triangular matrix $X$ such that $N_r+X$ is nilpotent of type $\lambda$. We call the resulting matrix an _upper nilpotent matrix completion_ of $N_r$.

The algorithm defined in [[LSN23]](#1) produces explicit constructions of upper nilpotent completions for $N_r$ in all cases that they exists (i.e., when $\lambda$ has at most $r$ parts). This algorithm is implemented in this project.

## Motivation

An upper nilpotent completion of $N_r$ gives rise to a _homogeneous Coxeter connection_ (and a solution to an important special case of the _Deligne–Simpson problem_). See [[LSN23]](#1) or [[KLMNS22]](#3) for details.

## How to use

All of the source code is contained in one file: `nilp-completion.py`.

The source file `nilp-completion.py` contains three functions. The first function, `nilpotent_completion()`, is the implementation of the upper nilpotent completion algorithm defined in [[LSN23]](#1). The latter two functions, `partition_conjugate()` and `nilpotency_type()`, can be used for testing. The docstrings contains instructions for use and examples.

## References

<a id="1">[LSN23]</a>
N. Livesay, D. S. Sage, and B. Nguyen (2023). "Explicit constructions of connections on the projective line with a maximally ramified irregular singularity". <a href="https://arxiv.org/abs/2303.06581">arXiv:2303.06581</a>.

<a id="2">[KL95]</a>
M. Krupnik and A. Leibman (1995). "Jordan structures of strictly lower triangular completions of nilpotent matrices". Integral Equations and Operator Theory **23**, 459–471.

<a id="3">[KLMNS22]</a>
M. C. Kulkarni, N. Livesay, J. P. Matherne, B. Nguyen, and D. S. Sage (2022). "The Deligne–Simpson problem for connections on $G_m$ with a maximally ramified singularity". Advances in Mathematics **408**, 108596.

## License

This project is licensed under the terms of the MIT License.
