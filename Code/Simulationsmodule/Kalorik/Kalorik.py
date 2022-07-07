import matplotlib.pyplot as plt
from numpy import power
import cantera as ct

from Code.Simulationsmodule.Kalorik.Kraftstoffe import get_Kraftstoff


def main():
    pass


class ThermodynamischeDaten(object):

    def __init__(self):
        self.__konventionelleKraftstoffe = ["Benzin E5", "Benzin E10", "E85", "Diesel", "Biodiesel", "Schweroel",
                                            "Pflanzenoel"]
        self.__Wasserstoff = "Wasserstoff"
        self.__Ethanol = "Ethanol"
        self.__Methanol = "Methanol"
        self.__CNG = "CNG"
        self.__LPG = "LPG"

        ct.add_directory('../Code/Simulationsmodule/Kalorik/')

        self.__air = "O2:0.21,N2:0.78,Ar:0.01"

    def get_u(self, Kraftstoff="Benzin E5", Temperatur=300, LambdaVG=1, Druck=1e5, Verbrennungszustand=0):

        u = 0
        if Kraftstoff in self.__konventionelleKraftstoffe:
            u = 0.1445 * (1356.8 + ((489.6 + (46.4 / power(LambdaVG, 0.93))) * (power(Temperatur - 273.15, 1)) * (
                    10 ** -2)) + ((7.768 + (3.36 / power(LambdaVG, 0.8))) * (power(Temperatur - 273.15, 2)) * (
                    10 ** -4)) - ((0.0975 + (0.0485 / power(LambdaVG, 0.75))) * (power(Temperatur - 273.15, 3)) * (
                    10 ** -6)))
        """elif Kraftstoff == self.__Wasserstoff:
            gas = ct.Solution('h2.yaml')
            gas.set_equivalence_ratio(phi=1 / LambdaVG, fuel="H2:1", oxidizer=self.__air, basis='mass')
            gas.TP = Temperatur, Druck
            u = gas.int_energy_mass
        elif Kraftstoff == self.__Ethanol:
            gas = ct.Solution('ethanol.yaml')
            gas.set_equivalence_ratio(phi=1 / LambdaVG, fuel="C2H5OH:1", oxidizer=self.__air, basis='mass')
            gas.TP = Temperatur, Druck
            u = gas.int_energy_mass
        elif Kraftstoff == self.__Methanol:
            gas = ct.Solution('methanol.yaml')
            gas.set_equivalence_ratio(phi=1 / LambdaVG, fuel="CH3OH:1", oxidizer=self.__air, basis='mass')
            gas.TP = Temperatur, Druck
            u = gas.int_energy_mass
        elif Kraftstoff == self.__CNG:
            gas = ct.Solution('gri30.yaml')
            gas.set_equivalence_ratio(phi=1 / LambdaVG, fuel="CH4:1", oxidizer=self.__air, basis='mass')
            gas.TP = Temperatur, Druck
            u = gas.int_energy_mass
        elif Kraftstoff == self.__LPG:
            gas = ct.Solution('gri30.yaml')
            gas.set_equivalence_ratio(phi=1 / LambdaVG, fuel="C3H8:0.75,C4H10:0.25", oxidizer=self.__air, basis='mass')
            gas.TP = Temperatur, Druck
            u = gas.int_energy_mass"""

        return u

    def get_du_dt(self, Kraftstoff="Benzin E5", Temperatur=300, LambdaVG=1, Druck=1e5, Verbrennungszustand=0):
        du_dt = 0

        if Kraftstoff in self.__konventionelleKraftstoffe:
            du_dt = 0.00001445 * ((3.36 / power(LambdaVG, 0.8)) + 7.768) * (2 * Temperatur - 546.3) - 0.0000004335 * (
                    (0.0485 / power(LambdaVG, 0.75)) + 0.0975) * power(Temperatur - 273.15, 2) + (
                            0.067048 / power(LambdaVG, 0.93)) + 0.707472
            du_dt *= 1000
        """elif Kraftstoff == self.__Wasserstoff:
            gas = ct.Solution('h2.yaml')
            print((1e-10+LambdaVG*abs(1-Verbrennungszustand)))
            gas.set_equivalence_ratio(phi=1 / (1e-10+LambdaVG*abs(1-Verbrennungszustand)), fuel="H2:1", oxidizer=self.__air, basis='mass')
            gas.TP = Temperatur, Druck
            # gas.equilibrate("TP") # nach experimenten stellte sich herraus, dass der Gleichgewichtszustand für cv nicht berechnet werden muss
            gas.equilibrate("TP")
            du_dt = gas.cv_mass
        elif Kraftstoff == self.__Ethanol:
            gas = ct.Solution('ethanol.yaml')
            gas.set_equivalence_ratio(phi=1 / LambdaVG, fuel="C2H5OH:1", oxidizer=self.__air, basis='mass')
            gas.TP = Temperatur, Druck
            gas.equilibrate("TP")
            du_dt = gas.cv_mass
        elif Kraftstoff == self.__Methanol:
            gas = ct.Solution('methanol.yaml')
            gas.set_equivalence_ratio(phi=1 / LambdaVG, fuel="CH3OH:1", oxidizer=self.__air, basis='mass')
            gas.TP = Temperatur, Druck
            gas.equilibrate("TP")
            du_dt = gas.cv_mass
        elif Kraftstoff == self.__CNG:
            gas = ct.Solution('gri30.yaml')
            gas.set_equivalence_ratio(phi=1 / LambdaVG, fuel="CH4:1", oxidizer=self.__air, basis='mass')
            gas.TP = Temperatur, Druck
            gas.equilibrate("TP")
            du_dt = gas.cv_mass
        elif Kraftstoff == self.__LPG:
            gas = ct.Solution('gri30.yaml')
            gas.set_equivalence_ratio(phi=1 / LambdaVG, fuel="C3H8:0.75,C4H10:0.25", oxidizer=self.__air, basis='mass')
            gas.TP = Temperatur, Druck
            gas.equilibrate("TP")
            du_dt = gas.cv_mass"""
        return du_dt


if __name__ == "__main__":
    main()
