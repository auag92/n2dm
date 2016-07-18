import numpy as np

def check_parameter_array(numelements, name, a):
    if numelements > 1:
        ao = np.array(a)

        shape = (numelements, numelements)
        assert(ao.shape == shape)

        an = np.zeros(shape)

        for i in range(numelements):
            for j in range(i, numelements):
                an[i,j] = ao[i,j]
                an[j,i] = ao[i,j]

    else:
        if isinstance(a, (int, float)):
            an = np.array([[a]])
        else:
            an = a.copy()

    return an

def _ChkLJarray(data, n, var, pot="LennardJones"):
    "Check the parameters of the LennardJones factory function."
    # Is data a scalar or an array?
    try:
        l = len(data)
    except TypeError:
        if n == 1:
            return [data]
        else:
            raise TypeError("%s: %s parameter cannot be a scalar (%d elements)."
                            % (pot, var, n))
    # Get the shape
    try:
        shp = data.shape
    except AttributeError:
        data = np.array(data)
        shp = data.shape
    if shp != (n,n) and shp != (n*n,):
        raise TypeError("%s: %s parameter should be a %dx%d array"
                        % (pot, var, n, n))
    if shp == (n,n):
        for i in range(n):
            for j in range(i+1,n):
                if data[i,j] != 0 and data[i,j] != data[j,i]:
                    raise ValueError("%s: %s matrix must be a symmetric or a lower triangular matrix" % (pot, var))
    return data.flat[:]

