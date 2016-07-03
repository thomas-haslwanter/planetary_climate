'''
Planetary database
Source for planetary data, and some of the data on
the moons, is http://nssdc.gsfc.nasa.gov/planetary/factsheet/

A "planet" contains basic planetary data:

- name .........Name of the planet
- a ........... Mean radius of planet [m]
- g ........... Surface gravitational acceleration [m/s**2]
- L ........... Annual mean solar constant (current) [W/m**2]
- albedo.... ...Bond albedo (fraction)
- mass ........ Mass of the planet [kg]
- rsm ......... Semi-major axis of orbit about Sun [m]
- year ........ Sidereal length of year [sec]
- eccentricity  Eccentricity [unitless]
- day ......... Mean tropical length of day [sec]
- obliquity ... Obliquity to orbit [degrees]
- Lequinox .... Longitude of equinox [degrees]
- Tsbar ....... Mean surface temperature [K]
- Tsmax ....... Maximum surface temperature [K]
- Tsmin ....... Minimum surface temperature [K]
- Type ........ [planet/none/moon]
- Around ...... Center of orbit

For gas giants, "surface" quantities are given at the 1 bar level

Todo
----

- Why do we have an entry "Lequinox", when no values are given?
- "Moon" is the only one where "Tsmin" is given
-  check Triton: obliquity to ecliptic;
                seasons are influenced by the inclination
                of Triton's orbit? (About 20 deg to Neptune's equator)

License
-------
BSD 3-clause (see https://www.w3.org/Consortium/Legal/2008/03-bsd-license.html)
'''

import pandas as pd
from io import StringIO

def get_properties():
    '''Properties of planets and some moons
    
    Returns
    -------
        props : pandas DataFrame
                contains important properties of the planets
                indexed by planet-name
                
        unit_dict: Dictionary
               Contains the corresponding (mks) units of the properties
                    
    Notes
    -----
        For the table:
            - "mass" is given in multiples of earth-mass
            - "day" is given in [h]
            - "years" is give in earth-days
            - for planets, the "rsm", "year", "L" and "eccentricity" are the same as for the corresponding planet.
        and the units are afterwards converted to SI
    '''
    
    data_string="""
    name,   a,          g,      albedo, L,      mass,   rsm,        year,       eccentricity,   day,        obliquity,  Lequinox,   Tsbar,  Tsmax,  Tsmin,  Type,   Around
    Earth,  6.371e6,    9.798,  0.306,  1367.6, 1.0,    149.60e9,   365.256,    0.0167,         24,         23.45,      None,       288,    None,   None,   planet, Sun
    Mercury,2.4397e6,   3.70,   0.119,  9126.6, 0.0553, 57.91e9,    87.969,     0.2056,         4222.6,     0.01,       None,       440.0,  725.0,  None,   planet, Sun
    Venus,  6.0518e6,   8.87,   0.750,  2613.9, 0.815,  108.21e9,   224.701,    0.0067,         2802.0,     177.36,     None,       737.0,  737.0,  None,   planet, Sun
    Mars,   3.390e6,    3.71,   0.250,  589.2,  0.107,  227.92e9,   686.98,     0.0935,         24.6597,    25.19,      None,       210.0,  295.0,  None,   planet, Sun
    Jupiter,69.911e6,   24.79,  0.343,  50.5,   317.8,  778.57e9,   4332.0,     0.0489,         9.9259,     3.13,       None,       165.0,  None,   None,   planet, Sun
    Saturn, 58.232e6,   10.44,  0.342,  14.90,  95.2,   1433.0e9,   10759.0,    0.0565,         10.656,     26.73,      None,       134.0,  None,   None,   planet, Sun
    Uranus, 25.362e6,   8.87,   0.300,  3.71,   14.5,   2872.46e9,  30685.4,    0.0457,         17.24,      97.77,      None,       76.0,   None,   None,   planet, Sun
    Neptune,26.624e6,   11.15,  0.290,  1.51,   17.2,   4495.06e9,  60189.0,    0.0113,         16.11,      28.32,      None,       72.0,   None,   None,   planet, Sun
    Pluto,  1.195e6,    0.58,   0.5,    0.89,   0.00218,5906.0e9,   90465.0,    0.2488,         153.2820,   122.53,     None,       50.0,   None,   None,   none,   Sun
    Moon,   1.737e6,    1.62,   0.11,   1367.6, 0.0123, 149.60e9,   365.256,    0.0167,         28.0,       None,       None,       None,   400.0,  100.0,  moon,   Earth
    Titan,  2.575e6,    1.35,   0.21,   14.90,  0.0225, 1433.0e9,   10759.0,    0.0565,         15.9452,    26.73       None,       95.0,   None,   None,   moon,   Saturn
    Europa, 1.560e6,    1.31,   0.67,   50.5,   0.008,  778.57e9,   4332.0,     0.0489,         3.551,      3.13,       None,       103.0,  125.0,  None,   moon,   Jupiter
    Triton, 1.3534e6,   0.78,   0.76,   1.51,   0.00359,4495.06e9,  60189.0,    0.0113,         5.877,      156.0,      None,       34.5,   None,   None,   moon,   Neptune """
    
    # Read in the data from that string
    props = pd.read_csv(StringIO(data_string), sep='[, ]+', engine='python', na_values='None')
    
    # Allow easy access to the data through the "formula"
    props.index = props['name']
    
    # Adjust to SI-units
    props.mass *= 5.9722e24  # mass of earth [kg]
    props.day *= 3600        # convert from [h] to [sec]
    props.year *= 24*3600    # convert from [earth-days] to [sec]
    
    # Create units dictionary
    units = ['none', 'm', 'm/s**2', 'fraction', 'W/m**2', 'kg', 'm', 'sec', 'None', 'sec', 'deg', 'deg', 'K', 'K', 'K', 'none', 'none']
    unit_dict = dict(zip(props.columns, units))
    
    return props, unit_dict

if __name__ == '__main__':
    
    # Get properties and units
    props, units = get_properties()
    # A few examples of how to access the planet-properties
    
    # List of properties stored
    print('The following informations are stored: -------')    
    print(props.columns)    
    
    # Show the properties of Earth
    print('The Properties of Earth: -------')
    print(props.loc['Earth'])    
    
    # Show the "mass"-s:
    print('Mass: ---------------')
    print(props.mass)        
        
    # Print all the "Critical*" Parameters    
    print('"T" Parameters: ---------------------')
    print(props.filter(regex='T.*'))    
    
    # Show the units    
    print('Units -----------------')
    for key in units:
        print('{0}: [{1}]'.format(key, units[key]))    

    input('Hit any key to continue ...')
