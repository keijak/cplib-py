def ternary_search(f, xmin, xmax, epsilon=1e-6):
    """Ternary search.

    Args:
      f: An objective function. Must be convex downward.
      xmin: The lower bound of the range to search.
      xmax: The upper bound of the range to search.
      epsilon: The epsilon value for judging convergence.
    """
    l, r = xmin, xmax
    while r - l > epsilon:
        d = (r - l) / 3.0
        m1 = l + d
        m2 = l + 2.0 * d
        if f(m1) < f(m2):
            r = m2
        else:
            l = m1
    return r
