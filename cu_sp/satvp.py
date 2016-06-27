'''
Functions to calculate saturation vapor pressure

All units are mks units

Change Log
----------
2016/06/15: cleaned up, and ported to Python 3 [ThH]                
2016/06/26: rewrote "moistAdiabat", eliminating the "integrator [ThH]

License
-------
BSD 3-clause (see https://www.w3.org/Consortium/Legal/2008/03-bsd-license.html)

'''

import numpy as np
import pandas as pd
import phys
from scipy.integrate import odeint

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
        p0 = 3589.
        MolecularWeight = 18.
        LatentHeat = 2.5e6
        properties = (T0, p0, MolecularWeight, LatentHeat)
        svp = satvp(properties)
    
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
    
    def __init__(self, properties, iceFlag='switch'):
        ''' Set the parameters of the object, so you can call it as a function 
        afterwards.
        If the parameters are passed as tuple, it has to have the following 
        sequence:
            - T0
            - p0
            - MolecularWeight
            - LatentHeat
        '''
        
        #Check if the first argument is a gas object. If not, assume
        #that the arguments give T0, p0, etc. as numbers
        if isinstance(properties, pd.Series):
            self.iceFlag = iceFlag
            
            self.T0  = properties.TriplePointT
            self.p0  = properties.TriplePointP
            self.M   = properties.MolecularWeight
            self.gas = properties
            
            if self.iceFlag == 'ice':
                self.L = properties.L_sublimation
            elif self.iceFlag == 'liquid':
                self.L = properties.L_vaporization
            else:
                self.iceFlag = 'switch'
                
        else:
            self.iceFlag = None
            
            self.T0 = properties[0]
            self.p0 = properties[1]
            self.M  = properties[2]
            self.L  = properties[3]
            
    def __call__(self, T):
        ''' Saturation vapor pressure for any substance, computed using the
        simplified form of Clausius-Clapeyron assuming the perfect gas law and
        constant latent heat
        
        Parameters
        ----------
            T : Temperatur [Kelvin]
        
        Returns
        -------
            pressure : saturation vapor pressure [Pascal]        
        '''
        
        #Decide which latent heat to use
            
        if self.iceFlag == 'switch':
            if not type(T) == np.ndarray:
                if T < self.gas.TriplePointT:
                    L = self.gas.L_sublimation
                else:
                    L = self.gas.L_vaporization
            else:
                L = np.ones_like(T) * self.gas.L_sublimation
                L[T>=self.gas.TriplePointT] = self.gas.L_vaporization
        else:
            L = self.L
            
        # Clausius-Clapeyron equation
        Rv = phys.Rstar/self.M
        p = self.p0 * np.exp( -(L/Rv)*(1./T - 1./self.T0) )
        return p

class MoistAdiabat:
    
    def __init__(self, condensible, noncon):
        '''Set up the function parameters'''        
        self.condensible = condensible
        self.noncon = noncon
        
        #Set up saturation vapor pressure function
        self.satvp = satvp(condensible)
        
        #Set up thermodynamic constants
        self.eps = condensible.MolecularWeight / noncon.MolecularWeight
        self.L   = condensible.L_vaporization
        self.Rc  = condensible.R
        self.cpc = condensible.cp
        self.Ra  = noncon.R
        self.cpa = noncon.cp
        
        self.ptop  = 100. #Default top of atmosphere
        
    def slope(self, log_T, log_pa):
        '''Derivative function defining the moist adiabat'''
        
        pa = np.exp(log_pa)
        T  = np.exp(log_T)
        r_sat = self.eps * self.satvp(T)/pa
        
        num = self.Ra * ( 1. + self.L*r_sat/(self.Ra*T) ) 
        den = self.cpa + (self.cpc + (self.L/(self.Rc*T)-1.) *self.L/T )*r_sat
        
        return num/den
        
    def __call__(self, p0, T0):
        '''Call to the resulting function'''
        
        # Initial conditions
        ln_p0, ln_T0 = np.log([p0, T0])
        
        p_air = np.logspace(np.log10(p0), np.log10(self.ptop), 101)
        ln_p_air = np.log(p_air)
        
        # Solve the differntial equation
        ln_T = odeint(self.slope, ln_T0, ln_p_air, full_output=True)
        T = np.exp(ln_T[0].ravel())
        
        # Total pressure is p_air + p_condensible
        p_c = self.satvp(T)
        p_total = p_air + p_c
    
        molarConL = p_c/p_total
        
        # Now compute mass specific concentration
        M_c  = self.condensible.MolecularWeight
        M_nc = self.noncon.MolecularWeight
        
        M_bar = molarConL*M_c + (1.-molarConL)*M_nc
        qL = (M_c/M_bar) * molarConL
        
        return p_total, T, molarConL, qL
    
if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    
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
        
    # Example of "MoistAdiabat" -----------------------------    
    import gases
    
    water = gases.gas_props.loc['H2O']        
    air = gases.gas_props.loc['air']    
    
    ma = MoistAdiabat(condensible=water, noncon=air)    
    p, T, molarCon, massCon = ma(1.e5,300.)

    plt.hold(True)
    plt.semilogy(T, p)
    plt.xlabel('Temperature [K]')        
    plt.ylabel('Pressure [Pa]')    
    plt.title('Moist Adiabat')
    plt.gca().invert_yaxis()
    plt.show()    

    # Example of "satvp" -----------------------------    
    # [To be done]
