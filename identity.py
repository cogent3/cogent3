import warnings

def identity():
    """
    This function is marked for removal in cogent3 version 3.1.0.
    A DeprecationWarning is raised to inform users that it will be removed.
    """
    warnings.warn("This function will be removed in cogent3 version 3.1.0", DeprecationWarning)
