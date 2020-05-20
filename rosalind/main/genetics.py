# import statements
from collections import defaultdict


def mendelian_genetics(k, m, n):
    """
        Parameters
        ----------
        k : int
            homozygous dominant
        m : int
            heterozygous
        n : int
            homozygous recessive

        Returns
        -------
        float
            probability of offspring showing dominant phenotype
    """
    dom_dom = k * (k-1)
    dom_het = k * m * 2
    het_het = m * (m-1) * 0.75
    het_res = m * n * 2 * 0.5 # never forget
    dom_res = k * n * 2
    total = (k + m + n) * (k + m + n - 1)
    return (dom_dom + dom_het + het_het + het_res + dom_res) / total



if __name__ == '__main__':
    print(mendelian_genetics(2,2,2))
    print(mendelian_genetics(19,21,24))