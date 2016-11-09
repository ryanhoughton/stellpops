import numpy as np

def nearest(x, dx=None, ceil=False, floor=False):
    """
    Purpose: to round a number to the nearest multiple of the dx you specify
             Particularly useful for plots 
             e.g. nearest(234,50) = 250 but nearest(234,50,/floor) = 200

             if dx==None, round to nearest Order of Mag
    """

    # sanity check
    if (ceil==True & floor==True): raise "Can't CEIL and FLOOR, just one or t'other"

    
    if dx==None:
        dx = 10.0**np.floor(np.log10(np.fabs(x)))
        if (~np.isfinite(np.log10(dx))): dx=10.0 # was rounding zero value

    near = float(x) / float(dx)
    if ceil:
        result = np.ceil(near)*dx
    elif floor:
        result = np.floor(near)*dx
    else:
        result = round(near)*dx

    return result
