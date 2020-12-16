import sys

def get_mass_of_element(label: str, average=False):
    """Returns the isotope specific mass of an element given it's symbol.
    https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii2&isotype=some


    Args:
        label (str): A standard element symbol. The normal symbol (e.g. Cl) will return the mass
        of the isotopically dominant element. I'll just add elements as I need them. Should probably
        add an option to average the isotopes and get that mass.
    """
    masses = {"H":     [1.00782503223,  1.007975],
              "D":     [2.01410177812,  1.007975],
              "T":     [3.0160492779,   1.007975],
              "HE3":   [3.0160293201,   4.002602],
              "HE":    [4.00260325413,  4.002602],
              "LI6":   [6.0151228874,   6.9675],
              "LI":    [7.0160034366,   6.9675],
              "BE":    [9.012183065,    9.0121831],
              "B10":   [10.01293695,    10.8135],
              "B":     [11.00930536,    10.8135],
              "C":     [12.0000000,     12.0106],
              "C13":   [13.00335483507, 12.0106],
              "C14":   [14.0032419884,  12.0106],
              "N":     [14.00307400443, 14.006855],
              "N15":   [15.00010889888, 14.006855],
              "O":     [15.99491461957, 15.9994],
              "O17":   [16.99913175650, 15.9994],
              "O18":   [17.99915961286, 15.9994],
              "F":     [18.99840316273, 18.998403163],
              "NE":    [19.9924401762,  20.1797],
              "NE21":  [20.993846685,   20.1797],
              "NE22":  [21.991385114,   20.1797]
              }
    
    if average:
        try:
            mass = masses[label.upper()][1]
            return mass
        except KeyError:
            print("We don't have a mass for that element label. Consider adding it yourself.")
            sys.exit(1)
    else:
        try:
            mass = masses[label.upper()][0]
            return mass
        except KeyError:
            print("We don't have a mass for that element label. Consider adding it yourself.")
            sys.exit(1)