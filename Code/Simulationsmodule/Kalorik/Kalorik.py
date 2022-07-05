import matplotlib.pyplot as plt
from numpy import power
import cantera as ct

from Code.Simulationsmodule.Kalorik.Kraftstoffe import get_Kraftstoff


def main():
    gas = ct.Solution('Code/Simulationsmodule/Kalorik/H2_reaction_v1a.yaml')
    gas()


class ThermodynamischeDaten(object):

    def __init__(self):
        self.__konventionelleKraftstoffe = ["Benzin E5", "Benzin E10", "E85", "Diesel", "Biodiesel", "Schweroel",
                                            "Pflanzenoel"]
        self.__Wasserstoff = "Wasserstoff"

    def get_u(self, Kraftstoff="Benzin E5", Temperatur=300, LambdaPL=1):

        u = 0

        if Kraftstoff in self.__konventionelleKraftstoffe:
            u = 0.1445 * (1356.8 + ((489.6 + (46.4 / power(LambdaPL, 0.93))) * (power(Temperatur - 273.15, 1)) * (
                    10 ** -2)) + ((7.768 + (3.36 / power(LambdaPL, 0.8))) * (power(Temperatur - 273.15, 2)) * (
                    10 ** -4)) - ((0.0975 + (0.0485 / power(LambdaPL, 0.75))) * (power(Temperatur - 273.15, 3)) * (
                    10 ** -6)))
        elif Kraftstoff == self.__Wasserstoff:
            u = -1

        return u

    def get_du_dt(self, Kraftstoff="Benzin E5", Temperatur=300, LambdaPL=1):
        du_dt = 0
        if Kraftstoff in self.__konventionelleKraftstoffe:
            du_dt = 0.00001445 * ((3.36 / power(LambdaPL, 0.8)) + 7.768) * (2 * Temperatur - 546.3) - 0.0000004335 * (
                    (0.0485 / power(LambdaPL, 0.75)) + 0.0975) * power(Temperatur - 273.15, 2) + (
                            0.067048 / power(LambdaPL, 0.93)) + 0.707472
            du_dt *= 1000
        elif Kraftstoff == self.__Wasserstoff:
            du_dt = -1
        return du_dt


if __name__ == "__main__":
    main()
