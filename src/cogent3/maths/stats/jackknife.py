import numpy as np

from cogent3.util.table import Table


__author__ = "Anuj Pahwa, Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Anuj Pahwa", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"


def index_gen(length):
    """returns a callable

    Parameters
    ----------
    length : int
        length of series of integers to generate

    Returns
    -------
    callable

    Notes
    -----
    When invoked with an int, returns all indices except that provided.
    The result can be used with a numpy.take(data, indices).

    >>> gen_series = index_gen(4)
    >>> gen_series(0)
    [1, 2, 3]
    >>> gen_series(1)
    [0, 2, 3]
    """
    data = tuple(range(length))

    def gen(i):
        temp = list(data)
        temp.pop(i)
        return temp

    return gen


class JackknifeStats(object):
    """Computes the jackknife statistic for a particular statistical function
    as outlined by 'Tukey's Jackknife Method' Biometry by Sokal/Rohlf."""

    def __init__(self, length, calc_stat, gen_index=index_gen):
        """

        Parameters
        ----------
        length : int
            The length of the data set (since data is not passed to this class).
        calc_stat : callable
            A callback function that computes the required statistic of a defined dataset.
        gen_index
            A callback function that generates a list of indices that are used to sub-sample the dataset.
        """

        super(JackknifeStats, self).__init__()
        self.n = length
        self.calc_stat = calc_stat
        self.gen_index = gen_index(self.n)
        self._subset_statistics = None
        self._pseudovalues = None
        self._jackknifed_stat = None
        self._sample_statistic = None
        self._standard_error = None

    def jackknife(self):
        """Computes the jackknife statistics and standard error"""
        n = self.n
        n_minus_1 = n - 1

        # compute the statistic in question on the whole data set
        self._sample_statistic = self.calc_stat(list(range(self.n)))
        n_sample_statistic = n * self._sample_statistic
        # compute the jackknife statistic for the data by removing an element
        # in each iteration and computing the statistic.
        subset_statistics = []
        pseudovalues = []
        for index in range(self.n):
            stat = self.calc_stat(self.gen_index(index))
            subset_statistics.append(stat)
            pseudovalue = n_sample_statistic - n_minus_1 * stat
            pseudovalues.append(pseudovalue)

        self._pseudovalues = np.array(pseudovalues)
        self._subset_statistics = np.array(subset_statistics)
        self._jackknifed_stat = self._pseudovalues.mean(axis=0)

        # Compute the approximate standard error of the jackknifed estimate
        # of the statistic
        variance = np.square(self._pseudovalues - self._jackknifed_stat).sum(axis=0)
        variance_norm = np.divide(variance, n * n_minus_1)
        self._standard_error = np.sqrt(variance_norm)

    @property
    def sample_stat(self):
        if self._sample_statistic is None:
            self.jackknife()
        return self._sample_statistic

    @property
    def jackknifed_stat(self):
        if self._jackknifed_stat is None:
            self.jackknife()
        return self._jackknifed_stat

    @property
    def standard_error(self):
        if self._standard_error is None:
            self.jackknife()
        return self._standard_error

    @property
    def sub_sample_stats(self):
        """Return a table of the sub-sample statistics"""

        # if the statistics haven't been run yet.
        if self._subset_statistics is None:
            self.jackknife()

        # generate table
        title = "Subsample Stats"
        rows = []
        for index in range(self.n):
            row = [index]
            subset_statistics = self._subset_statistics[index]
            try:
                for value in subset_statistics:
                    row.append(value)
            except TypeError:
                row.append(subset_statistics)
            rows.append(row)

        header = ["i"]
        subset_stats = self._subset_statistics[0]

        try:
            num_datasets = len(subset_stats)
            for i in range(num_datasets):
                header.append(f"Stat_{i}-i")
        except TypeError:
            header.append("Stat-i")

        return Table(data=rows, header=header, title=title)

    @property
    def pseudovalues(self):
        """Return a table of the Pseudovalues"""

        # if the statistics haven't been run yet.
        if self._pseudovalues is None:
            self.jackknife()

        # detailed table
        title = "Pseudovalues"
        rows = []
        for index in range(self.n):
            row = [index]
            pseudovalues = self._pseudovalues[index]
            try:
                for value in pseudovalues:
                    row.append(value)
            except TypeError:
                row.append(pseudovalues)
            rows.append(row)

        header = ["i"]
        pseudovalues = self._pseudovalues[0]

        try:
            num_datasets = len(pseudovalues)
            for i in range(num_datasets):
                header.append(f"Pseudovalue_{i}-i")
        except TypeError:
            header.append("Pseudovalue-i")

        return Table(data=rows, header=header, title=title)

    @property
    def summary_stats(self):
        """Return a summary table with the statistic value(s) calculated for the
        the full data-set, the jackknife statistics and standard errors."""

        # if the statistics haven't been run yet.
        if self._jackknifed_stat is None:
            self.jackknife()

        header = ["Sample Stat", "Jackknife Stat", "Standard Error"]
        title = "Summary Statistics"
        rows = np.vstack(
            (self._sample_statistic, self._jackknifed_stat, self._standard_error)
        )
        rows = rows.transpose()
        return Table(header=header, data=rows, title=title)
