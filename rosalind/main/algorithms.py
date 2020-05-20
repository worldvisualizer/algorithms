# import statements
from collections import defaultdict


def fibonacci(n):
    """
        Parameters
        ----------
        n : int
            number of generations

        Returns
        -------
        int
            n-th fibonacci number
    """
    fibs = []
    for i in range(n + 1):
        if i == 0:
            fibs.append(0)
        elif i == 1 or i == 2:
            fibs.append(1)
        else:
            fibs.append(fibs[i-2] + fibs[i-1])
    return fibs[n]


def fibonacci_custom(n, multi):
    """
        Parameters
        ----------
        n : int
            number of generations
        multi : int
            number of multiplication factor for each gen.
            e.g.) number of pairs each pair generates
                  at every generations
                  (1 rabbit pair -> births 3 pairs)

        Returns
        -------
        int
            n-th fibonacci number
    """
    fibs = []
    for i in range(n + 1):
        if i == 0:
            fibs.append(0)
        elif i == 1 or i == 2:
            fibs.append(1)
        else:
            fibs.append(fibs[i-2] * multi + fibs[i-1])
    return fibs[n]


def hamming_distance(string1, string2):
    """
        Parameters
        ----------
        string1 : str
        string2 : str

        Returns
        -------
        int
            hamming distance of two strings
    """
    dist_counter = 0
    for n in range(len(string1)):
        if string1[n] != string2[n]:
            dist_counter += 1
    return dist_counter


if __name__ == '__main__':
    