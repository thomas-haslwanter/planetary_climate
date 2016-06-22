'''
Functions to calculate saturation vapor pressure

All units are mks units

Change Log
----------
06/15/2016: cleaned up, and ported to Python 3 [ThH]                

License
-------
BSD 3-clause (see https://www.w3.org/Consortium/Legal/2008/03-bsd-license.html)

'''

import numpy as np
import pandas as pd

def satvp_H2O(T, mode = 'general'):
    '''
    Saturation vapor pressure (SVP) computation used in the GFDL climate model. 
    
    See smithsonian meteorological tables page 350.
    Original source: GFDL climate model, circa 1995
    
    Parameters
    ----------
        T : Temperatur [Kelvin]
        
        type : string
            'general' ... blends over from "water" saturation ( T_Celsius>0 )
                         to "ice" saturation ( T_Celsius<-20 ) as the temperature
                         falls below 0C.
            'water' ... SVP over liquid water
            'heymsfield' ... alternate formula for SVP over liquid water
            'ice' ..... SVP over liquid ice, valid between -153C and 0C

    
    Returns
    -------
        pressure : saturation vapor pressure [Pascal]
    '''
    
    # Make sure that user does not use the function with Celsius
    if np.any(T<=0):
        print('Inputs to "satvp_H2O" are in [Kelvin], and have to be >0!')
        raise ValueError
    
    if mode == 'water':
        esbasw = 1013246.0
        tbasw =  373.16
        
        aa  = -7.90298 * (tbasw/T-1)
        b   =  5.02808 * np.log10(tbasw/T)
        c   = -1.3816e-07 * (  10.**( ((1-T/tbasw)*11.344)-1 )  )
        d   =  8.1328e-03 * (  10.**( ((tbasw/T-1)*(-3.49149))-1)  )
        
        e   = np.log10(esbasw)
        es_H2O  = 10.**(aa+b+c+d+e)
        
        return es_H2O * 0.1  #Convert to Pascals

    elif mode == 'ice':
        esbasi = 6107.1
        tbasi =  273.16
        
        aa  = -9.09718 * (tbasi/T-1.0)
        b   = -3.56654 * np.log10(tbasi/T)
        c   =  0.876793* (1.0-T/tbasi)
        e   = np.log10(esbasi)
        es_ice = 10.**(aa+b+c+e)
    
        return es_ice * 0.1  #Convert to Pascals
        
    elif mode == 'general':
        T_Celsius = T - 273.16
        
        # Write the function such that it also works with arrays
        svp_out = np.nan * np.ones_like(T)
        
        # Set the temperature ranges
        warm = T_Celsius > 0
        medium = np.logical_and(-20 <= T_Celsius, T_Celsius <= 0)
        cold = T_Celsius < -20
        
        # Return the corresponding values
        svp_out[warm] = satvp_H2O(T[warm], mode='water')
        svp_out[cold] = satvp_H2O(T[cold], mode='ice')
        # linear transition between "ice" and "water"
        svp_out[medium] = satvp_H2O(T[medium], mode='water') + \
                T_Celsius[medium]/20 * \
                   (satvp_H2O(T[medium],'water') - satvp_H2O(T[medium],'ice'))
        
        return svp_out
    
    elif mode == 'heymsfield':
        ts = 373.16
        sr = 3.0057166
        
        # Vapor pressure over water. Heymsfield formula
        ar = ts/T
        br = 7.90298 * (ar-1.)
        cr = 5.02808 * np.log10(ar);
        dw = (1.3816E-07) * (10.**(11.344*(1.-1./ar))-1.)
        er = 8.1328E-03 * ((10.**(-(3.49149*(ar-1.))) )-1.)
        vp = 10.**(cr-dw+er+sr-br)
        vp = vp * 1.0e02
        
        return(vp)        

class satvp:
    '''
    This class provides functions for the calculation of SVP from
    a simplified form of Clausius-Clapeyron. The use of the function is
    simplified, by setting up an object that stores the thermodynamic data
    needed, so it doesn't have to be re-entered each time. 
    Because of the __call__ method, once the object is created, it can be
    invoked like a regular function.
    
    Examples
    --------
    To set up a function e(T) that approximates the saturation
    vapor presure for a substance which has a latent heat of
    2.5e6 J/kg, a molecular weight of 18 and has vapor pressure
    3589. Pa at a temperature of 300K, create the function using:
    
        T0 = 300.
        e0 = 3589.
        svp = satvp(T0, e0, MolecularWeight=18., LatentHeat=2.5e6)
    
    and afterward you can invoke it simply as "svp(T)", where T
    is whatever temperature you want to evaluate it for.
    
    Alternately, satvp can be called with a gas object
    as the first argument, e.g.
    
        from gases import gas_props
        svp = satvps(gas_props.loc['CO2'])
    
    If no other arguments are given, the latent heat of sublimation
    will be used when "svp(T)" is called for temperatures below the triple
    point, and the latent heat of vaporization will be used for
    temperatures above the triple point. To allow you to force
    one or the other latent heats to be used, satvps_function takes
    an optional second argument when the first argument is a gas
    object.  Thus,
    
        svp = satvps(gas_props.loc['CO2'], 'ice')
        
    will always use the latent heat of sublimation, regardless of T,
    while
    
        svp = satvps(gas_props.loc['CO2'], 'liquid')
        
    will always use the latent heat of vaporization.
    '''
    
    def __init__(self, Gas_or_T0, e0_or_iceFlag=None, MolecularWeight=None,
                 LatentHeat=None):
        
        #Check if the first argument is a gas object. If not, assume
        #that the arguments give T0, e0, etc. as numbers
        self.iceFlag = e0_or_iceFlag
        if isinstance(Gas_or_T0, pd.Series):
            self.gas = Gas_or_T0
            self.M = Gas_or_T0.MolecularWeight
            self.T0 = Gas_or_T0.TriplePointT
            self.e0 = Gas_or_T0.TriplePointP
            
            if self.iceFlag == 'ice':
                self.L = Gas_or_T0.L_sublimation
            elif self.iceFlag == 'liquid':
                self.L = Gas_or_T0.L_vaporization
            else:
                self.iceFlag = 'switch'
                
            #self.M = Gas_or_T0.MolecularWeight
        else:
            self.L = LatentHeat
            self.M = MolecularWeight
            self.T0 = Gas_or_T0
            self.e0 = e0_or_iceFlag
            
    def __call__(self, T):
        
        #Decide which latent heat to use
        if self.iceFlag == 'switch':
            if T < self.gas.TriplePointT:
                L = self.gas.L_sublimation
            else:
                L = self.gas.L_vaporization
        else:
            L = self.L
            
        return satvp_simple(T, self.T0, self.e0, self.M, L)
            
    def satvp_simple(T, T0, e0, MolecularWeight, LatentHeat):
        ''' Saturation vapor pressure for any substance, computed using the
        simplified form of Clausius-Clapeyron assuming the perfect gas law and
        constant latent heat
        
        Parameters
        ----------
            T : Temperatur [Kelvin]
            T0 :
            e0 :
            MolecularWeight : 
            LatentHeat : 
        
        Returns
        -------
            pressure : saturation vapor pressure [Pascal]
        '''
        
        Rv = Rstar/MolecularWeight 
        
        return e0*math.exp(-(LatentHeat/Rv)*(1./T - 1./T0))

            


class MoistAdiabat:
    '''
    MoistAdiabat is a class which creates a callable object
    used to compute the moist adiabat for a mixture consisting
    of a condensible gas and a noncondensing gas.  The gases
    are specified as gas objects. By default, the condensible
    is water vapor and the noncondensible is modern Earth Air,
    if the gases are not specified.

    Usage:
          To create a function m that computes the moist
          adiabat for the gas Condensible mixed with the gas
          Noncondensible, do
                m = phys.MoistAdiabat(Condensible,Noncondensible)
          For example, to do a mixture of condensible CO2 in
          noncondensing N2, do
                m = phys.MoistAdiabat(phys.CO2,phys.N2)
          Once you have created the function, you give it
          the surface partial pressure of the noncondensible
          and the surface temperature when you call it, and it
          returns arrays consisting of pressure, temperature,
          molar concentration of the condensible, and mass
          specific concentration of the condensible. For example:
                p,T,molarCon,massCon = m(1.e5,300.)
          for a surface noncondensible pressure of 1.e5 Pascal and
          surface temperture of 300K.  The values returned
          are arrays. The pressure returned is total pressure at
          each level (condensible plus noncondensible).  By default,
          the compution chooses the pressure values on which to return
          the results.  For some purposes, you might want the results
          specified on a list of pressures of your own choosing.  The
          computation allows for this, by offering an interpolation
          option which returns the result interpolated to a pressure
          grid of your own choice, which is specified as an optional
          third argument to the function. Thus, to get the
          pressure values on a list consisting of [1000.,5000.,10000.] Pa,
          you would do:
                p,T,molarCon,massCon = m(1.e5,300.,[1000.,5000.,10000.])
          The calculation is still done at high resolution to preserve
          accuracy, but the results are afterward intepolated to the grid
          you want using polynomial interpolation. For your convenience,
          the pressure returned on the left hand side is a copy of
          the pressure list you specified as input.
          
   **ToDo:
      - Add help strings and documentation
      - The way the help strings for gas objects are
           set up makes the argument help box for the
           creator useless.  Fix this somehow
      - Add controls on resolution, top of atmosphere, etc.
           Do we want this to return molar or mass concentration?
           Maybe do both, but have result stored as an attribute
    '''
    
    def __init__(self, condensible, noncon):
        self.condensible = condensible
        self.noncon = noncon
        
        #Set up saturation vapor pressure function
        self.satvp = satvp(condensible)
        
        #Set up thermodynamic constants
        self.eps = condensible.MolecularWeight/noncon.MolecularWeight
        self.L = condensible.L_vaporization
        self.Ra = noncon.R
        self.Rc = condensible.R
        self.cpa = noncon.cp
        self.cpc = condensible.cp
        
        def slope(logpa,logT):
            '''Set up derivative function for integrator '''
            
            pa = math.exp(logpa)
            T = math.exp(logT)
            qsat = self.eps*(self.satvp(T)/pa)
            num = (1. + (self.L/(self.Ra*T))*qsat)*self.Ra
            den = self.cpa + (self.cpc + (self.L/(self.Rc*T) - 1.)*(self.L/T))*qsat
            
            return num/den
        
        self.slope = slope
        self.ptop = 1000. #Default top of atmosphere
        self.step = -.05 #Default step size for integration
        
    def __call__(self,ps,Ts,pgrid = None):
        
        #Initial conditions
        step = self.step  #Step size for integration
        ptop = self.ptop #Where to stop integratoin
        
        logpa = math.log(ps)
        logT = math.log(Ts)
        ad = integrator(self.slope,logpa,logT,step )
        
        #Initialize lists to save results
        pL = [math.exp(logpa) + self.satvp(math.exp(logT))]
        molarConL = [self.satvp(math.exp(logT))/pL[0]]
        TL = [math.exp(logT)]
        
        #Integration loop
        p = 1.e30 #Dummy initial value, to get started
        while p > ptop:
            ans = next(ad)
            pa = math.exp(ans[0])
            T = math.exp(ans[1])
            p = pa+self.satvp(T)
            pL.append(p)
            molarConL.append(self.satvp(T)/p)
            TL.append(T)
                
        pL = np.array(pL)
        TL = np.array(TL)
        molarConL = np.array(molarConL)
        
        #Now compute mass specific concentration
        Mc = self.condensible.MolecularWeight
        Mnc = self.noncon.MolecularWeight
        Mbar = molarConL*Mc +(1.-molarConL)*Mnc
        qL = (Mc/Mbar)*molarConL
        
        #The else clause below interpolates to a
        #specified pressure array pgrid, if desired.
        # interp is a class defined in ClimateUtilities
        #which creates a callable object which acts like
        #an interpolation function for the listed data give
        #as arguments.
        if pgrid == None:
            return pL,TL,molarConL,qL
        else:
            T1 = interp(pL,TL)
            mc1 = interp(pL,molarConL)
            q1 = interp(pL,qL)
            T = np.array([T1(pp) for pp in pgrid])
            mc = np.array([mc1(pp) for pp in pgrid])
            q = np.array([q1(pp) for pp in pgrid])
            
            return np.array(pgrid),T, mc, q            

if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    
    '''
    # For handling units, I use the package "pint"
    # This is a bit of an overkill here, but "pint" is a cool package for
    # working with different units
    import pint
    
    # Select the temperature range
    units = pint.UnitRegistry()
    Q_ = units.Quantity
    temp = Q_(np.arange(-30, 50, dtype=np.float), units.degC)
    T = temp.to(units.K).magnitude    
    
    # Example of "satvp_H2O" -----------------------------    
    plt.plot(temp.magnitude, satvp_H2O(T), label='general')
    plt.hold(True)
    plt.plot(temp.magnitude, satvp_H2O(T, mode='heymsfield'), label='heymsfield')
    plt.legend()
    
    plt.xlabel('Temperature [C]')
    plt.title('Saturation Vapor Pressure H2O')
    plt.ylabel('Pressure [Pa]')
    plt.show()
        
    # Example of "satvp" -----------------------------    
    
    '''
    # Example of "MoistAdiabat" -----------------------------    
    import gases
    
    water = gases.gas_props.loc['H2O']        
    air = gases.gas_props.loc['air']    
    
    ma = MoistAdiabat(condensible=water, noncon=air)    
    p,T,molarCon,massCon = ma(1.e5,300.)