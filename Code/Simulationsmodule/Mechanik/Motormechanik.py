import scipy
from numpy import sin, cos, power, sqrt, pi


class Motormechanik(object):

    def __init__(self,Kurbelradius, Drehzahl, LambdaPL, Kompressionsvolumen, Hubraum, Pleuellaenge):
        from Code.Simulationsmodule.Prozessrechnung.Prozessrechnung import Realprozessrechnung
        self.__Kurbelradius = Kurbelradius
        self.__Drehzahl = Drehzahl
        self.__lambdaPL = LambdaPL
        self.__Kompressionsvolumen = Kompressionsvolumen
        self.__Hubraum = Hubraum
        self.__Pleuellaenge = Pleuellaenge

    """def __del__(self):
        del self.__RPR"""

    def Kolbengeschwindigkeit(self, phi):
        rad = phi * pi / 180
        ret = self.__Kurbelradius * ((self.__Drehzahl / 60) * 2 * pi) * (sin(rad) + (
                (cos(rad) * (self.__lambdaPL * sin(rad))) / sqrt(
            1 - power(self.__lambdaPL * sin(rad), 2))))
        return ret

    def Hubvolumen(self, phi):
        rad = phi * pi / 180
        ret = self.__Kompressionsvolumen + (
                    self.__Hubraum / (2 * self.__Kurbelradius)) * (
                      self.__Kurbelradius * (1 - cos(rad)) + self.__Pleuellaenge * (
                      1 - sqrt(1 - power(self.__lambdaPL * sin(rad), 2))))
        return ret

    def dV(self, phi):
        ret = scipy.misc.derivative(self.Hubvolumen, phi, dx=1e-6)
        return ret
