from Numeric import array, matrixmultiply, transpose
from LinearAlgebra import inverse

def savitzky_golay(window_size=None,order=2):
    if window_size is None:
        window_size = order + 2

    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window size is too small for the polynomial")

    # A second order polynomial has 3 coefficients
    order_range = range(order+1)
    half_window = (window_size-1)//2
    B = array(
        [ [k**i for i in order_range] for k in range(-half_window, half_window+1)] )

    #           -1
    # [  T     ]      T
    # [ B  * B ]  *  B
    M = matrixmultiply(
           inverse(  matrixmultiply(transpose(B), B)),
           transpose(B)
           )
    return M

def savitzky_golay_weights(window_size=None, order=2, derivative=0):
    # The weights are in the first row
    # The weights for the 1st derivatives are in the second, etc.
    return savitzky_golay(window_size, order)[derivative]
