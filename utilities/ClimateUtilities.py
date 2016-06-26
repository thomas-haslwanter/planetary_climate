'''
Dependencies
------------
numpy, matplotlib

License
-------
BSD 3-clause (see https://www.w3.org/Consortium/Legal/2008/03-bsd-license.html)

ToDo
----
*Make and import a PathFinder module, which the instructor
 will customize for the site.  This module will help the
 students find locations of datasets and chapter scripts.
 Could provide "links" to find script that produced a given
 figure in the text.  Alternately, we could just
 define directory strings like WorkbookDatasets here in
 ClimateUtilities.

'''

import string
import numpy as np

#==============================================
#---Section: Math utilities-------------------
#==============================================

class Dummy:
    ''' A dummy class useful for passing parameters.
    Just make an instance, then add new members
    at will.
    '''
    pass

def polint(xa,ya,x):
    ''' Polynomial interpolation and extrapolation (adapted
    from Numerical Recipes.
    
    It's used in Romberg extrapolation, but could be useful
    for polynomial OLR fits and so forth as well. Also
    needs online documentation
    '''
    n = len(xa)
    if not (len(xa) == len(ya)):
            print("Input x and y arrays must be same length")
            return "Error"
        #Set up auxiliary arrays
    c = np.zeros(n, dtype=np.float)
    d = np.zeros(n, dtype=np.float)
    c[:] = ya[:]
    d[:] = ya[:]
    #Find closest table entry
    ns = 0
    diff = abs(xa[0]-x)
    for i in range(n):
            difft = abs(xa[i]-x)
            if difft < diff:
                diff = difft
                ns = i
    y=ya[ns]
    for m in range(1,n): 
        for i in range(n-m):
            ho=xa[i]-x
            hp=xa[i+m]-x
            w=c[i+1]-d[i]
            c[i] = ho*w/(ho-hp)
            d[i] = hp*w/(ho-hp)
        if 2*ns < (n-m):
            dy = c[ns]
        else:
            ns -= 1
            dy = d[ns]
        y += dy
        #You can also return dy as an error estimate. Here 
        #to keep things simple, we just return y.
    
    return y

class interp:
    '''
    Class for doing polynomial interpolation
    from a table, using polint
    
    Usage:
    ------
        Let xa be a list of independent variable
        values and ya be a list of the corresponding
        dependent variable values. Then, to create a function
        f (actually a callable object, techically) that interpolates
        or extrapolates to any value x, create f using
                    f = interp(xa,ya)
        Then you can get the value you want by invoking f(x)
        for your desired x.
        
        By default, the interpolator does fourth-order interpolation
        using the four nearest neighbors. You can change this by
        using an optional third argument to the creator. For
        example
        
                    f = interp(xa,ya,8)
        will use the 8 nearest neighbors (if they are available)
    
    '''
    
    def __init__(self, xa, ya, n=4):
        self.xa = np.array(xa)
        self.ya = np.array(ya)
        self.n = n
        
    def __call__(self,x):
        #Find the closes index to x
        if self.xa[0] < self.xa[-1]:
            i = np.searchsorted(self.xa,x)
        else:
            i = np.searchsorted(-self.xa,-x)
            
        i1 = max(i-self.n,0)
        i2 = min(i+self.n,len(self.xa))
        
        return polint(self.xa[i1:i2],self.ya[i1:i2],x)


class BetterTrap:
    ''' Quadrature (definite integral) by Romberg extrapolation.
    **ToDo: Add documentation and help string
    
    Before developing a general quadrature class, we'll
    implement a class which efficiently carries out trapezoidal rule
    integration with iterative refinement
    '''
    
    def __init__(self,f,params,interval,nstart):
        self.f = f
        self.n = nstart
        self.interval = interval
        self.params = params
        self.integral = self.dumbTrap(nstart)
        
    def dumbTrap(self,n):
        a = self.interval[0]
        b = self.interval[1]
        dx = (b-a)/n
        sum = dx*(self.f(a,self.params)+self.f(b,self.params))/2.
        for i in range(1,n):
            x = a+i*dx
            sum = sum + self.f(x,self.params)*dx
            
        return sum
    
    def refine(self):
        ''' Compute the sum of f(x) at the
        midpoints between the existing intervals.
        To get the refinement of the trapezoidal
        rule sum we just add this to half the
        previous result
        '''
        
        sum = 0.
        a = self.interval[0]
        b = self.interval[1]
        dx = (b-a)/self.n
        
        #Remember: n is the number of subintervals, not the number of endpoints.
        #Therefore we have one midpoint per subinterval. Keeping that
        #in mind helps us get the range of i right in the following loop
        
        for i in range(self.n):
            sum = sum + self.f(a+(i+.5)*dx,self.params)*(dx/2.)
            
        #The old trapezoidal sum was multiplied by the old dx. To get its
        #correct contribution to the refined sum, we must multiply it by .5,
        #because the new dx is half the old dx
        self.integral = .5*self.integral + sum
        
        #Update the number of intervals
        self.n = 2*self.n

class romberg:
    ''' Here I define a class called
    romberg, which assists in carrying out evaluation of
    integrals using romberg extrapolation. It assumes polint has
    been imported
    '''
    
    def __init__(self,f,nstart=4):
        self.nstart = nstart
        self.trap = None
        
        #-------------------------------------------------
        #This snippit of code allows the user to leave the parameter argument
        #out of the definition of f if it isn't needed
        self.fin = f
        
        #Find the number of arguments of f and append a
        #parameter argument if there isn't any.
        nargs = f.__code__.co_argcount
        if nargs == 2:
            self.f = f
        elif nargs ==1:
            def f1(x,param):
                return self.fin(x)
            self.f = f1
        else:
            name = f.__name__
            print('Error: %s has wrong number of arguments'%name)
        #-----------------------------------------------------
        
        #We keep lists of all our results, for doing Romberg extrapolation.
        #These are re-initialized after each call
        self.nList = []
        self.integralList = []
        
    def refine(self):
        self.trap.refine()
        self.integralList.append(self.trap.integral)
        self.nList.append(self.trap.n)
        dx = [1./(n*n) for n in self.nList]
        return polint(dx,self.integralList,0.)
    
    def __call__(self,interval,params=None,tolerance=1.e-6):
        ''' Use a __call__ method to return the result. The
        __call__ method takes the interval of integration
        as its mandatory first argument,takes an optional
        parameter argument as its second argument, and
        an optional keyword argument specifying the accuracy
        desired.
        **ToDo: Introduce trick to allow parameter argument of
        integrand to be optional, as in Integrator.  Also, make
        tolerance into a keyword argument
        '''
        
        self.nList = []
        self.integralList = []
        #Make a trapezoidal rule integrator
        self.trap = BetterTrap(self.f,params,interval,self.nstart)
        self.nList.append(self.nstart)
        self.integralList.append(self.trap.integral)
        #
        #Refine initial evaluation until 
        oldval = self.refine()
        newval = self.refine()
        while abs(oldval-newval)>tolerance:
            oldval,newval = newval,self.refine()
            
        return newval
        
class newtSolve:
    '''
    Newton method solver for function of 1 variable
    A class implementing Newton's method for solving f(x) = 0.

    Usage: solver = newtSolve(f), where f is a function with
    calling sequence f(x,params). Values of x such that
    f(x,params) = 0 are
    then found by invoking solver(guess), where guess
    is the initial guess.  The solver returns the string
    'No Convergence' if convergence fails. The argument
    params allows parameters to be passed to the function.
    It can be left out of the function definition if you don't
    need it. Note that params can be any Python object at all
    (including,e.g.,lists, functions or more complex user-defined objects)

    Optionally, one can specify the derivative function
    in the creator,e.g. solver = newtSolve(f,fp).
    If the derivative function isn't specified, the solver
    computes the derivative approximately using a centered
    difference. Note that in either case you can access
    the derivative function by invoking solver.deriv(x)
    As for f, fp can be optionally defined with a parameter
    argument if you need it. The same parameter object is
    passed to f and fp. 

    Use solver.setParams(value) to set the parameter object
    Alternately, the parameter argument can be passed as
    an optional second argument in the solver call. (see
    example below).

    Adjustable constants:
     eps         Increment for computing numerical approximation to
                 the derivative of f
     tolerance   Accuracy criterion for ending the iteration
                 (an approximation to the error in the root)
     nmax        maximum number of iterations

    e.g. to change the maximum number of iterations for an instance
    of the class, set solver.nmax = 10 .
    '''

    def __init__(self, f, fprime=None):
        self.fin = f
        #Find the number of arguments of f and append a
        #parameter argument if there isn't any.
        nargs = f.__code__.co_argcount
        if nargs == 2:
            self.f = f
        elif nargs ==1:
            def f1(x,param):
                return self.fin(x)
            self.f = f1
        else:
            name = f.__name__
            print('Error: %s has wrong number of arguments'%name)
        self.eps = 1.e-6
        def deriv(x,params):
            return (self.f(x+self.eps,params)- self.f(x-self.eps,params))/(2.*self.eps)
        if fprime == None:
            self.deriv = deriv 
        else:
            #A derivative function was explicitly specified
            #Check if it has a parameter argument
            nargs = fprime.__code__.co_argcount
            if nargs == 2:
                self.deriv = fprime #Has a parameter argument
            elif nargs == 1:
                self.fprimein = fprime
                def fprime1(x,param):
                    return self.fprimein(x)
                self.deriv = fprime1
            else:
                name = fprime.__name__
                print('Error: %s has wrong number of arguments'%name)
        self.tolerance = 1.e-6
        self.nmax = 100
        self.params = None
        
    def __call__(self,xGuess,params = None):
        if not (params == None):
            self.setParams(params)
        x = xGuess
        for i in range(self.nmax):
            dx = (self.f(x,self.params)/self.deriv(x,self.params))
            x = x - dx
            if abs(dx) < self.tolerance:
                return x
        return 'No Convergence'
    
    def setParams(self,params):
        #**ToDo: Check if f1 has a parameter argument
        #defined, and complain if it doesn't
        self.params = params
        
    def scan(self,interval,n=10):
        '''
        Finds initial guesses to roots in a specified
        interval, subdivided into n subintervals.
        e.g. if the instance is called "solver"
        solver.scan([0.,1.],100) generates a list
        of guesses between 0. and 1., with a resolution
        of .01. The larger n is, the less is the chance that
        a root will be missed, but the longer the search
        will take.  If n isn't specified, the default value is 10
        
        ToDo:
        Replace this with a bisection search, allowing user
        to specify the maximum number of distinct guesses that
        need to be found.
        '''
        guessList = []
        dx = (interval[1]-interval[0])/(n-1)
        flast = self.f(interval[0],self.params)
        for x in [interval[0]+ i*dx for i in range(1,n)]:
            fnow = self.f(x,self.params)
            if ((fnow >= 0.)&(flast <=0.)) or ((fnow <= 0.)&(flast >=0.)):
                guessList.append(x)
            flast = fnow
        return guessList

if __name__ == '__main__':
    
    # ---------------- Usage Examples for "newSolve" -------------------------
    
    def show_root(initial_guess, func):
        '''Show initial guess and resulting root value'''
        print('With initial guess {0}, we get the first root {1}.'.format(
            initial_guess, func(initial_guess)) )
        
    # Example 1: Function without parameters
    print('--- Example 1 ---')
    def g(x):
        return x*x - 1.
    
    roots = newtSolve(g)
    show_root(initial_guess=2., func=roots)

    # Example 2: Function with parameters
    print('--- Example 2 ---')
    def g(x,constants):
        return constants.a*x*x - constants.b
    
    roots = newtSolve(g)
    constants = Dummy()
    constants.a = 1.
    constants.b = 2.
    roots.setParams(constants)
    
    show_root(initial_guess=2., func=roots)
    show_root(initial_guess=1., func=roots)

    # Example 2a:
    print('--- Example 2a ---')
    # Instead of using roots.setParam(...) we could do
    print( roots(2.,constants) )
    print( roots(1.) )    # the parameters are remembered
    constants.a = 3.
    print( roots(1.,constants) )   # We changed the constants

    # Example 3: using scan to find initial guesses
    print('--- Example 3 ---')
    def g(x):
        return x*x - 1.
    
    roots = newtSolve(g)
    guesses = roots.scan([-2.,2.], 100)
    for guess in guesses:
        print( roots(guess) )
        
