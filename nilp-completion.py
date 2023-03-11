r"""
Explicit upper nilpotent completions

The function `nilpotent_completion()` (defined herein) is an
implementation of Algorithm 1 specified in [1]_. See the article and/or
docstrings for details and examples. The function `nilpotency_type()`
may be used for testing.

References
----------
.. [1] N. Livesay, D. S. Sage, and B. Nguyen. Explicit constructions of
       connections on the projective line with a maximally ramified
       irregular singularity (2023).

"""

from numpy import matrix
from numpy import linalg  # for `matrix_rank()`


def nilpotent_completion(n, r, partition):
    r""" Returns an upper nilpotent completion of type `partition`

    Returns matrix of the form `N_r + X` with the following properties:
        - `N_r` and `X` are `n` by `n` matrices
        - `N_r` is the matrix with 1s on the rth subdiagonal
            (and 0s elsewhere)
        - `X` is strictly upper triangular
        - `N_r + X` is nilpotent of type `partition`; i.e., `partition`
            gives the sizes of the Jordan blocks for `N_r + X`.

    Parameters
    ----------
    n : int
        Must satisfy `1 < n`.

    r : int
        Must satisfy `1 <= r < n`.

    partition : iterable of ints
        An iterable of positive integers in monotone decreasing order
        which sum to `n`.

    Returns
    -------
    gamma : numpy matrix
        Returns `n` by `n` matrix that is nilpotent of type `partition`

    Notes
    -----
    Implementation is based on Algorithm 1 in [1]_.

    References
    ----------
    .. [1] N. Livesay, D. S. Sage, and B. Nguyen. Explicit
           constructions of connections on the projective line with a
           maximally ramified irregular singularity (2023).

    Examples
    --------
    The example below returns a `5` by `5` matrix of the form `N_2+X`
    which is nilpotent of type `[3, 2]`. Since `N_2` already has type
    `[3, 2]`, it simply returns `N_2`.

    >>> nilpotent_completion(5, 2, [3, 2])
    matrix([[0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0]])

    The remaining examples are less trivial; the upper-triangle
    contains some 1s.

    >>> nilpotent_completion(5, 2, [4, 1])
    matrix([[0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0]])

    >>> nilpotent_completion(5, 2, [5])
    matrix([[0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1],
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0]])
    """

    # raise Exception if any issues with input `r`, `n`, or `partition`
    if type(r) != int or type(n) != int:
        raise TypeError('inputs `n` and `r` must be integers')
    if not 1 <= r < n:
        raise ValueError('inputs `n` and `r` should satisfy 1<=r<n')
    if not hasattr(partition, '__iter__'):
        raise TypeError('input `partition` must be iterable')
    if any(type(partition[i]) != int for i in range(len(partition))):
        raise TypeError('input `partition` should be an iter of ints')
    if sum(partition) != n:
        raise ValueError('values in `partition` should sum to n')
    if any(partition[i] < partition[i+1] for i in
           range(len(partition)-1)):
        raise ValueError('input `partition` must be ' +
                         'monotone decreasing')
    if partition[-1] <= 0:
        raise ValueError('values in `partition` must be positive')
    if len(partition) > r:
        raise ValueError('input `partition` must have len at most `r`')

    # Initialization of variables:
    # - `gamma` is an `n` by `n` numpy matrix, initialized to `N_r`
    # i.e., the matrix with 1s on `r`th subdiagonal and 0s elsewhere.
    # - `tall` is a list of all the parts in `partition` larger than
    #   `ceiling{n/r}` possibly with extra copies of `ceiling{n/r}`
    # - `short` is a list of all the parts in `partition` smaller than
    #   `floor{n/r}`, possibly with extra copies of `floor{n/r}`
    # - `stock` and `cion` are integer "pointers"
    gamma = matrix([[(1 if row-col == r else 0) for col in range(n)]
                    for row in range(n)])
    ceiling = n//r + (1 if n % r != 0 else 0)
    tall = []
    short = []
    parts_equal_to_ceiling = 0
    parts_equal_to_floor = 0

    for part in partition:
        if part > ceiling:
            tall.append(part)
        elif part == n//r:
            parts_equal_to_floor += 1
        elif n % r != 0 and part == n//r+1:
            parts_equal_to_ceiling += 1
        elif part < n//r:
            short.append(part)
    # append some number of copies of `floor{n/r}` and `ceiling{n/r}`
    # to `short` and `tall`
    short += [n//r] * max(0, parts_equal_to_floor - (r-n % r))
    tall += [n//r+1] * max(0, parts_equal_to_ceiling - n % r)
    # initialize `stock` and `cion` pointers
    stock = min(parts_equal_to_floor, r-n % r)
    cion = stock + 1

    # helper functions
    def e_ij(i, j):
        """returns `n` by `n` matrix with 1 in `(i+1, j+1)`-entry and
        0s elsewhere (standard basis vector for `n` by `n` matrices)"""
        return matrix([[(1 if row == j and col == i else 0)
                        for col in range(n)] for row in range(n)])

    def height(i):
        """returns height of column `i` in initial embeding of `N_r`"""
        return n//r if i < r-n % r else ceiling

    def ordinal(i, m):
        """ returns the ordinal of `m`th vertex from top in column `i`
        in initial embedding of `N_r` """
        return n % r+i+(m-1)*r if i < r-n % r else i-(r-n % r)+(m-1)*r

    # two more helpful variables
    stock_height = height(stock)
    stock_top_ordinal = ordinal(stock, 1)

    # execution of primary loops (Loop 1 and Loop 2)
    while short:  # Loop 1
        if tall[0]-stock_height == height(cion)-short[0]+1 \
           and height(cion) < height(cion+1):  # Case 1a
            gamma += e_ij(ordinal(cion+1, height(cion+1)-short.pop(0)),
                          stock_top_ordinal)
            tall.pop(0)
            stock, cion = cion, cion+2
            stock_height = height(stock)
            stock_top_ordinal = ordinal(stock, 1)
        elif tall[0]-stock_height > height(cion)-short[0]:  # Case 1b
            m = height(cion) - short.pop(0)
            gamma += e_ij(ordinal(cion, m), stock_top_ordinal)
            stock_top_ordinal = ordinal(cion, 1)
            stock_height += m
            cion += 1
        elif tall[0]-stock_height == height(cion)-short[0]:  # Case 1c
            gamma += e_ij(ordinal(cion, height(cion)-short.pop(0)),
                          stock_top_ordinal)
            tall.pop(0)
            stock, cion = cion+1, cion+2
            stock_height = height(stock)
            stock_top_ordinal = ordinal(stock, 1)
        else:  # Case 1d (tall[0]-stock_height < height(cion)-short[0])
            m = tall.pop(0) - stock_height
            gamma += e_ij(ordinal(cion, m), stock_top_ordinal)
            stock = cion
            stock_top_ordinal = ordinal(cion, m+1)
            stock_height = height(cion) - m
            cion = stock + 1
    while tall:  # Loop 2
        while tall[0]-stock_height > ceiling:  # Little Loop
            gamma += e_ij(ordinal(cion, height(cion)),
                          stock_top_ordinal)
            stock_top_ordinal = ordinal(cion, 1)
            stock_height += height(cion)
            cion += 1
        if tall[0]-stock_height > height(cion):
            if height(cion+1) == ceiling:  # Case 2a
                gamma += e_ij(ordinal(cion+1, height(cion+1)),
                              stock_top_ordinal)
                stock, cion = cion, cion+2
                stock_height = height(stock)
                stock_top_ordinal = ordinal(stock, 1)
            else:  # Case 2b
                gamma += e_ij(ordinal(cion, height(cion)),
                              stock_top_ordinal)
                gamma += e_ij(ordinal(cion+1, 1), ordinal(cion, 1))
                stock_top_ordinal = ordinal(cion+1, 2)
                stock, cion = cion+1, cion+2
                stock_height = height(stock) - 1
        elif tall[0]-stock_height == height(cion):  # Case 2c
            gamma += e_ij(ordinal(cion, height(cion)),
                          stock_top_ordinal)
            stock, cion = cion+1, cion+2
            stock_top_ordinal = ordinal(stock, 1)
            stock_height = height(stock)
        else:  # Case 2d (tall[0]-stock_height < height(cion))
            m = tall[0] - stock_height
            gamma += e_ij(ordinal(cion, m), stock_top_ordinal)
            stock_top_ordinal = ordinal(cion, m+1)
            stock, cion = cion, cion+1
            stock_height = height(stock) - m
        tall.pop(0)

    return gamma

# useful functions for testing


def partition_conjugate(partition):
    r""" The conjugate (aka dual) of an integer partition.

    Parameters
    ----------
    partition : iterable of ints
        An iterable of positive integers in monotone decreasing order
        which sum to `n`.

    Returns
    -------
    conjugate : iterable of ints
        The conjugate partition of `partition`.

    Examples
    --------

    >>> partition_conjugate([5])
    [1, 1, 1, 1, 1]

    >>> partition_conjugate([1, 1, 1, 1, 1])
    [5]

    >>> partition_conjugate([5, 4, 4, 2, 1, 1])
    [6, 4, 3, 3, 1]

    >>> partition_conjugate([3, 3, 3, 2, 2])
    [5, 5, 3]
    """
    # step 1
    conjugate = [len(partition)] * partition[-1]
    # step 2
    i = len(partition) - 2
    while i >= 0 and partition[i] == partition[-1]:
        i -= 1
    if i >= 0:
        conjugate += [i+1] * (partition[i]-partition[-1])
    # steps 3+
    while i >= 0:
        j = i - 1
        while j >= 0 and partition[j] == partition[i]:
            j -= 1
        if j >= 0:
            conjugate += [j+1] * (partition[j] - partition[i])
        i = j
    return conjugate


def nilpotency_type(X):
    r"""
    Returns Jordan block sizes (aka type) of a nilpotent square matrix.

    The matrix is expected to be square. The "nilpotency type" of a
    matrix is the partition with parts given by the sizes of its Jordan
    blocks. An assertion error is raised if matrix is found to be
    non-nilpotent.

    Parameters
    ----------
    X : numpy matrix
        Must be square; i.e., must be `n` by `n` for some integer `n>0`

    Returns
    -------
    iterable of integers
        Returns a "partition"; i.e. an iterable of integers in
        monotone decreasing order which sum to `n`, where `n` is the
        number of rows (equiv. columns) in `X`. The entries (aka parts)
        in this partition are precisely the sizes of the Jordan blocks
        in `X`. This partition is also known as the "type" of `X`.

    Notes
    -----
    Suppose `X` is nilpotent of order k; i.e., `X^k == 0`. The type of
    `X` is given by the conjugate (aka dual) of the partition
    `(nullity(X), nullity(X^2)-nullity(X), nullity(X^3)-nullity(X^2),
        ..., nullity(X^k)-nullity(X^{k-1})`.

    Examples
    --------
    Each of the examples below is in Jordan canonical form. Thus it is
    straight-forward to verify nilpotency (since they are strictly
    upper-triangular). It is also easy to read the nilpotency type
    (the tuple consisting of the sizes of the Jordan blocks).

    >>> nilpotency_type(matrix([\
    ... [0,1],\
    ... [0,0]]))
    [2]

    >>> nilpotency_type(matrix([\
    ... [0,1,0,0,0],\
    ... [0,0,1,0,0],\
    ... [0,0,0,0,0],\
    ... [0,0,0,0,1],\
    ... [0,0,0,0,0]]))
    [3, 2]

    >>> nilpotency_type(matrix([\
    ... [0,1,0,0,0],\
    ... [0,0,1,0,0],\
    ... [0,0,0,1,0],\
    ... [0,0,0,0,1],\
    ... [0,0,0,0,0]]))
    [5]

    >>> nilpotency_type(matrix([\
    ... [0,0,0,0,0],\
    ... [0,0,0,0,0],\
    ... [0,0,0,0,0],\
    ... [0,0,0,0,0],\
    ... [0,0,0,0,0]]))
    [1, 1, 1, 1, 1]
    """
    n = len(X)
    X_to_power = X.copy()
    # by Rank-Nullity Theorem, `nullity(X) == n - rank(X)`
    previous_nullity = n - linalg.matrix_rank(X_to_power)
    nullities = [previous_nullity]
    power = 1
    while power <= n and previous_nullity < n:
        power += 1
        X_to_power *= X
        current_nullity = n - linalg.matrix_rank(X_to_power)
        nullities.append(current_nullity - previous_nullity)
        previous_nullity = current_nullity
    # if the `n`th power of an `n` by `n` matrix is not zero,
    # then the matrix is not nilpotent; raise an error in this case
    if power > n:
        raise ValueError('\n'+repr(X)+' is not nilpotent')
    return partition_conjugate(nullities)


# if this file is run directly from command line, then doctests are run
if __name__ == "__main__":
    import doctest
    doctest.testmod()
