#!/usr/bin/env/ python

__author__ = "AdeeB NooR"
__copyright__ = "Copyright 2010, The GenomeDB project"
__credits__ = ["Adeeb Noor", "Rob Knight"]
__license__ = "GPL"
_version__ = "1.0-dev"
__maintainer__ = "Adeeb Noor"
__email__ = "adeeb.noor@colorado.edu"
__status__ = "Development"


def fact(x): 
    """ Performs a Factorial of x    
    x -- a number 
    """

    if x==0:
        return 1
    else:
        res = 1
        for i in range(1,x+1):
            res *= i
        return res

 
def n_choose_k(n,k):
    """ Performs a combination of n,k numbers 
	x, n are numbers 
    """

    if(n<0 or k < 0 or (n-k) < 0):
        return 0
	if ((n-k) < 0):
	   return 1

    return float(fact(n)/(fact(n-k)*fact(k)))

def hypergeometric_pmf(first_cell,total,total_first_row,total_first_column):
    """  Performs a Probability mass function (PMF) 
         K - first cell at the table , N - total number of table
	 m - total of first row , n - total of first colum 
    """

    k = first_cell
    N = total
    m = total_first_row
    n = total_first_column
    if((N-m) <= 0 or (n-k) <= 0):
        return 1
    return float((n_choose_k(m,k) * n_choose_k(N-m,n-k)/n_choose_k(N,n)))

def hypergeometric_cdf(values,size_population,number_items,number_samples):
    """ Performs a cummulative distribution function (CDF) 
        x - values of matrix , M - size of population
	K - numbder of items to be counted . N - number of samples
    """

    x = values
    M = size_population
    K = number_items
    N = number_samples
    cdf_cal = 0
    i = 0
    while(i<=x):
        cdf_cal +=  float(n_choose_k(K,i))*n_choose_k((M-K),(N-i))/float(n_choose_k(M,N))
        i = i + 1
    return float(cdf_cal)

def hypergeometric_sf(values,size_population,number_items,number_samples):
    """ Performs a Survival function(SF)
    """

    return float(1-(hypergeometric_cdf(values,size_population,number_items,number_samples)))

def fisher_exact_test(input_array):
    """ Performs a Fisher exact test on a 2x2 contingency table 
        Returns a tuple of (odds ratio, two-tailed P-value).
        Examples: >>> fisher_exact_test([[2, 5], [10, 4]]) (0.16, 0.15882352941176558)
    """

    if(input_array[1][0] == 0 or input_array[0][1] == 0):
	Ratio = 1
    else:
       Ratio = input_array[0][0] * input_array[1][1] / float(input_array[1][0] * input_array[0][1])
    row1 = input_array[0][0] + input_array[0][1]
    row2 = input_array[1][0] + input_array[1][1]
    column = input_array[0][0] + input_array[1][0]
    mode = int(float((column + 1) * (row1 + 1)) / (row1 + row2 + 2))
    probability_Exact = hypergeometric_pmf(input_array[0][0], row1 + row2, row1, column)
    probability_Mode = hypergeometric_pmf(input_array[0][0], row1 + row2, row1, column)
    if input_array[0][0] == mode :
        return Ratio
    elif input_array[0][0] < mode :
        probability_lower = hypergeometric_cdf(input_array[0][0], row1 + row2, row1, column)
        
        """ Search for upper half
        """ 
        min = mode
        max = column
        G = -1
        while min != max :
            G = max if (max == min + 1 and G == min) else \
                    (max + min) / 2

            probability_G = hypergeometric_pmf(G, row1 + row2, row1, column)
            if probability_G <= probability_Exact and hypergeometric_pmf(G - 1, row1 + row2, row1, column) > probability_Exact :
                break
            elif probability_G < probability_Exact :
                max = G
            else :
                min = G

        if G == -1 and min == max :
            G = min

        return Ratio, probability_lower + hypergeometric_sf(G - 1, row1 + row2, row1, column)
    else :
        probability_upper = hypergeometric_sf(input_array[0][0] - 1, row1 + row2, row1, column);

        """ Search for lower half 
        """
        min = 0
        max = mode
        G = -1
        while min != max :
            G = max if (max == min + 1 and G == min) else \
                    (max + min) / 2
            probability_G = hypergeometric_pmf(G, row1 + row2, row1, column);

            if probability_G <= probability_Exact and hypergeometric_pmf(G + 1, row1 + row2, row1, column) > probability_Exact :
                break;
            elif probability_G <= probability_Exact  :
                min = G
            else :
                max = G

        if G == -1 and min == max :
            G = min

        return Ratio, probability_upper + hypergeometric_cdf(G, row1 + row2, row1, column)

