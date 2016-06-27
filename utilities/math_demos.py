'''
Demonstration of how to execute some basic operations:
    - Root finding
    - Interpolation
    - Integration 
    
Shows how the functions in the old "ClimateUtilities.py" can be replaced
with numpy and scipy routines.

Dependencies
------------
numpy, matplotlib, scipy

License
-------
BSD 3-clause (see https://www.w3.org/Consortium/Legal/2008/03-bsd-license.html)

'''

import string
import numpy as np
import matplotlib.pyplot as plt

#==============================================
#---Section: Math utilities-------------------
#==============================================

if __name__ == '__main__':
    
    # Usage Examples:
    # ===============
    
    # 1. Finding roots of a function (replaces the old "newtSolv")
    # ------------------------------------------------------------
    from scipy.optimize import newton
    
    def show_root(initial_guess, root):
        '''Show initial guess and resulting root value'''
        print('With initial guess {0}, we get the root {1}.'.format(
            initial_guess, root) )
        
    # Example 1.1: Function without parameters
    print('--- Example 1 ---')
    def g(x):
        return x*x - 1.
    
    initial_guess = 2.
    root = newton(g, initial_guess)
    show_root(initial_guess, root)

    # Example 1.2: Function with parameters
    print('--- Example 2 ---')
    def g(x, a, b):
        return a*x*x - b
    
    initial_guess = 2.
    a, b = 1, 2
    root = newton(g, initial_guess, args=(a, b))
    show_root(initial_guess, root)
    
    initial_guess = 1.
    root = newton(g, initial_guess, args=(a, b))
    show_root(initial_guess, root)
        
    # 2. Interpolation (replaces the old "interp" and "polint")
    # ---------------------------------------------------------    
    from scipy.interpolate import interp1d
    
    # Generate some data
    x = np.linspace(0, 10, num=11, endpoint=True)
    y = np.cos(-x**2/9.0)
    xnew = np.linspace(0, 10, num=41, endpoint=True)
    
    # Make a "linear" and a "cubic" interpolation object
    f = interp1d(x, y)
    f2 = interp1d(x, y, kind='cubic')
    
    # Calculate and plot the interpolated data
    plt.plot(x, y, 'o', label='data')
    plt.plot(xnew, f(xnew), '-', label='linear')
    plt.plot(xnew, f2(xnew),'--', label='cubic')
    plt.legend()
    plt.show()


    # 3. Romberg Integration (replaces the old "romberg")
    # ---------------------------------------------------
    from scipy.integrate import romberg
    
    # Define a function    
    def f(x):
        return x**2        
    
    # Integrate the function, from -1 to 2:
    a, b = -1, 2
    integral = romberg(f, a, b)
    
    print('The integral of f(x) between {0} and {1} is: {2}'.format(a, b, integral))