from numpy import power
import cantera as ct


def main():
    pass


class ThermodynamischeDaten(object):

    def __init__(self):
        self.__konventionelleKraftstoffe = ["Benzin","Diesel"]
        self.__Wasserstoff = "Wasserstoff"
        self.__CNG = "CNG"

        ct.add_directory('../Code/Simulationsmodule/Kalorik/')
        ct.add_directory('..')

        self.__air = "O2:0.21,N2:0.78,Ar:0.01"

    def get_thermo_data(self, Kraftstoff="Benzin", Temperatur=300, LambdaVG=1, Druck=1e5):
        ret = {}

        def get_data():
            ret["cv"] = gas.cv_mass
            ret["cp"] = gas.cp_mass
            ret["kappa"] = gas.cp_mass / gas.cv_mass
            ret["u"] = gas.int_energy_mass

        if Kraftstoff in self.__konventionelleKraftstoffe:
            du_dt = 0.00001445 * ((3.36 / power(LambdaVG, 0.8)) + 7.768) * (2 * Temperatur - 546.3) - 0.0000004335 * (
                    (0.0485 / power(LambdaVG, 0.75)) + 0.0975) * power(Temperatur - 273.15, 2) + (
                            0.067048 / power(LambdaVG, 0.93)) + 0.707472
            ret["cv"] = du_dt * 1000
            ret["cp"] = du_dt + 287
            ret["kappa"] = ret["cp"] / ret["cv"]
            u = 0.1445 * (
                    1356.8 + ((489.6 + (46.4 / power(LambdaVG, 0.93))) * (power(Temperatur - 273.15, 1)) * (
                    10 ** -2)) + ((7.768 + (3.36 / power(LambdaVG, 0.8))) * (power(Temperatur - 273.15, 2)) * (
                    10 ** -4)) - ((0.0975 + (0.0485 / power(LambdaVG, 0.75))) * (power(Temperatur - 273.15, 3)) * (
                    10 ** -6)))
            ret["u"] = u * 1000
        elif Kraftstoff == self.__Wasserstoff:
            gas = ct.Solution('h2.yaml')
            gas.set_equivalence_ratio(phi=1 / LambdaVG, fuel="H2:1", oxidizer=self.__air, basis='mass')
            gas.TP = Temperatur, Druck
            gas.equilibrate("TP")
            get_data()
        elif Kraftstoff == self.__CNG:
            #GRIMECH12 aus Website
            gas = ct.Solution('grimech12.yaml')
            gas.set_equivalence_ratio(phi=1 / LambdaVG, fuel="CH4:1", oxidizer=self.__air, basis='mass')
            gas.TP = Temperatur, Druck
            gas.equilibrate("TP")
            get_data()

        return ret


if __name__ == "__main__":
    main()
