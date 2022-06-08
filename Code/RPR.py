import scipy.misc
from numpy import cos, sin, exp, power, pi, arange, zeros, sqrt
from time import time
from scipy.integrate import odeint, simps


# Alle Einheiten sind insofern nicht anders angegeben SI Einheiten!


class Realprozessrechnung(object):

    def __init__(self, Kurbelwinkelaufloesung=1):

        self.Tmax = 0
        self.execTime = 0
        self.Kreisprozessarbeit = 0
        self.Wirkunsgrad = 0
        self.pmax = 0
        self.IndizierterMitteldruck = 0
        self.Drehmoment = 0
        self.Leistung = 0

        self.Genauigkeit = Kurbelwinkelaufloesung
        self.ASP_pro_Umdrehung = 0.5  # 1: 2-Takt 0.5: 4-Takt
        self.T0 = 300  # Starttemperatur 300K
        self.p0 = 100000  # Startdruck 1bar
        self.phiES = 180 + 40  # +40  # 40° nach UT 180+40
        self.phiAOE = 540 - 60  # -60  # 60° vor UT 540-60
        self.ZZP = 360 - 10  # 20° vor OT Zuenden 360-20
        self.ZV = 0
        self.Drehzahl = 5800  # in U/min
        self.cv = 718
        self.R = 287
        self.Bohrung = 0.081
        self.Hub = 0.0864
        self.Kurbelradius = self.Hub / 2
        self.Pleuellaenge = 0.144
        self.epsilon = 10.3
        self.lambdaVerbrennung = 1
        self.Lmin = 14.7
        self.RGA = 0.04
        self.Hu = 41000000
        self.m_vibe = 1
        self.phiBA = self.ZZP + self.ZV
        self.phiBD = 50
        self.spezEnthalpieBB = 0
        self.isLuftansaugend = False
        self.Zylinderanzahl = 4

        #self.AnzahlStuetzstellen = int((self.phiAOE - self.phiES) / self.Genauigkeit) + 1
        #self.lambdaPL = self.Kurbelradius / self.Pleuellaenge
        #self.Hubraum = pi * 0.25 * self.Bohrung * self.Bohrung * self.Hub
        #self.Kompressionsvolumen = self.Hubraum / (self.epsilon - 1)
        #self.Vmin = self.Kompressionsvolumen
        #self.Vmax = self.Kompressionsvolumen + self.Hubraum
        #self.Luftdichte = self.p0 / (self.R * self.T0)
        #self.Luftmasse = self.Hubraum * self.Luftdichte
        #self.Brennstoffmasse = self.Luftmasse / (self.lambdaVerbrennung * self.Lmin)
        #self.RGA_Masse = (self.RGA / (1 - self.RGA)) * (self.Luftmasse + self.Brennstoffmasse)
        #self.m = self.Luftmasse + self.Brennstoffmasse + self.RGA_Masse
        #self.Qmax = self.Brennstoffmasse * self.Hu
        #self.omega = (self.Drehzahl / 60) * 2 * pi

        self.phiKW = arange(self.phiES, self.phiAOE + self.Genauigkeit, self.Genauigkeit)
        self.T = zeros(self.get_AnzahlStuetzstellen())
        self.deltaT = zeros(self.get_AnzahlStuetzstellen())
        self.deltaV = zeros(self.get_AnzahlStuetzstellen())
        self.deltaW = zeros(self.get_AnzahlStuetzstellen())
        self.deltaQb = zeros(self.get_AnzahlStuetzstellen())
        self.deltaQw = zeros(self.get_AnzahlStuetzstellen())
        self.deltaU = zeros(self.get_AnzahlStuetzstellen())
        self.V = zeros(self.get_AnzahlStuetzstellen())
        self.p = zeros(self.get_AnzahlStuetzstellen())
        self.lambdaVG = zeros(self.get_AnzahlStuetzstellen())
        self.Brennstoffmasse_phi = zeros(self.get_AnzahlStuetzstellen())
        self.V_darstellung = zeros(self.get_AnzahlStuetzstellen())
        self.p_darstellung = zeros(self.get_AnzahlStuetzstellen())

    def get_Genauigkeit(self):
        return self.Genauigkeit

    def set_Genauigkeit(self,value):
        self.Genauigkeit = value

    def get_ASP_pro_Umdrehung(self):
        return self.ASP_pro_Umdrehung

    def set_ASP_pro_Umdrehung(self,value):
        self.ASP_pro_Umdrehung = value

    def get_p0(self):
        return self.p0

    def set_p0(self,value):
        self.p0 = value

    def get_T0(self):
        return self.T0

    def set_T0(self,value):
        self.T0 = value

    def get_phiES(self):
        return self.phiES

    def set_phiES(self, value):
        self.phiES = value
        
    def get_phiAOE(self):
        return self.phiAOE
        
    def set_phiAOE(self,value):
        self.phiAOE = value
        
    def get_ZZP(self):
        return self.ZZP
        
    def set_ZZP(self,value):
        self.ZZP = value
    
    def get_ZV(self):
        return self.ZV
        
    def set_ZV(self,value):
        self.ZV = value

    def get_Drehzahl(self):
        return self.Drehzahl

    def set_Drehzahl(self,value):
        self.Drehzahl = value

    def get_cv(self):
        return self.cv

    def set_cv(self,value):
        self.cv = value

    def get_R(self):
        return self.R

    def set_R(self,value):
        self.R = value

    def get_Bohrung(self):
        return self.Bohrung

    def set_Bohrung(self,value):
        self.Bohrung = value

    def get_Hub(self):
        return self.Hub

    def set_Hub(self,value):
        self.Hub = value

    def get_Pleuellaenge(self):
        return self.Pleuellaenge

    def set_Pleuellaenge(self,value):
        self.Pleuellaenge = value

    def get_epsilon(self):
        return self.epsilon

    def set_epsilon(self,value):
        self.epsilon = value

    def get_lambdaVerbrennung(self):
        return self.lambdaVerbrennung

    def set_lambdaVerbrennung(self,value):
        self.lambdaVerbrennung = value

    def get_Lmin(self):
        return self.Lmin

    def set_Lmin(self,value):
        self.Lmin = value

    def get_RGA(self):
        return self.RGA

    def set_RGA(self,value):
        self.RGA = value

    def get_Hu(self):
        return self.Hu

    def set_Hu(self,value):
        self.Hu = value

    def get_m_vibe(self):
        return self.m_vibe

    def set_m_vibe(self,value):
        self.m_vibe = value

    def get_phiBA(self):
        return self.get_ZZP() + self.get_ZV()

    def get_Kurbelradius(self):
        return self.get_Hub() / 2

    def get_phiBD(self):
        return self.phiBD

    def set_phiBD(self,value):
        self.phiBD = value

    def get_spezEnthalpieBB(self):
        return self.spezEnthalpieBB

    def set_spezEnthalpieBB(self,value):
        self.spezEnthalpieBB = value

    def get_isLuftansaugend(self):
        return self.isLuftansaugend

    def set_isLuftansaugend(self,value):
        self.isLuftansaugend = value

    def get_Zylinderanzahl(self):
        return self.Zylinderanzahl

    def set_Zylinderanzahl(self,value):
        self.Zylinderanzahl = value

    def get_AnzahlStuetzstellen(self):
        return int((self.get_phiAOE() - self.get_phiES()) / self.get_Genauigkeit()) + 1

    def get_lambdaPL(self):
        return self.get_Kurbelradius() / self.get_Pleuellaenge()

    def get_Hubraum(self):
        return pi * 0.25 * self.get_Bohrung() * self.get_Bohrung() * self.get_Hub()

    def get_Kompressionsvolumen(self):
        return self.get_Hubraum() / (self.get_epsilon() - 1)

    def get_Vmin(self):
        return self.get_Kompressionsvolumen()

    def get_Vmax(self):
        return self.get_Kompressionsvolumen() + self.get_Hubraum()

    def get_Luftdichte(self):
        return self.get_p0() / (self.get_R() * self.get_T0())

    def get_Luftmasse(self):
        return self.get_Hubraum() * self.get_Luftdichte()

    def get_Brennstoffmasse(self):
        return self.get_Luftmasse() / (self.get_lambdaVerbrennung() * self.get_Lmin())

    def get_RGA_Masse(self):
        return (self.get_RGA() / (1 - self.get_RGA())) * (self.get_Luftmasse() + self.get_Brennstoffmasse())

    def get_m(self):
        return self.get_Luftmasse() + self.get_Brennstoffmasse() + self.get_RGA_Masse()

    def get_Qmax(self):
        return self.get_Brennstoffmasse() * self.get_Hu()

    def get_omega(self):
        return (self.get_Drehzahl() / 60) * 2 * pi

    def get_Kreisprozessarbeit(self):
        return simps(self.p, self.V)

    def get_Wirkungsgrad(self):
        return self.get_Kreisprozessarbeit() / self.get_Qmax()

    def get_pmax(self):
        return max(self.p)

    def get_T_array(self):
        return self.T

    def get_Tmax(self):
        return max(self.get_T_array())[0]

    def get_IndizierterMitteldruck(self):
        return self.get_Kreisprozessarbeit() / (100000 * (self.get_Vmax() - self.get_Vmin()))

    def get_Drehmoment(self):
        return self.get_ASP_pro_Umdrehung() * self.get_IndizierterMitteldruck() * 10 ** 5 * self.get_Hubraum() * self.get_Zylinderanzahl() / (
                2 * pi)

    def get_Leistung(self):
        return self.get_Drehmoment() * self.get_omega() * 1.36 / 1000

    def get_p_array(self):
        return self.p

    def get_phi_q_50(self):
        return 0

    def Hubvolumen(self, phi):
        rad = phi * pi / 180
        ret = self.get_Kompressionsvolumen() + (self.get_Hubraum() / (2 * self.get_Kurbelradius())) * (
                self.get_Kurbelradius() * (1 - cos(rad)) + self.get_Pleuellaenge() * (
                1 - sqrt(1 - power(self.get_lambdaPL() * sin(rad), 2))))
        return ret

    def dV(self, phi):
        ret = scipy.misc.derivative(self.Hubvolumen, phi, dx=1e-6)
        return ret

    def dQb(self, phi):
        if phi < self.get_phiBA():
            ret = 0
        else:
            ret = (self.get_Qmax() / self.get_phiBD()) * 6.908 * (self.get_m_vibe() + 1) * power((phi - self.get_phiBA()) / self.get_phiBD(),
                                                                               self.get_m_vibe()) * exp(
                -6.908 * power((phi - self.get_phiBA()) / self.get_phiBD(), self.get_m_vibe() + 1))
        return ret

    def dQw(self, phi, Temp):
        return 0

    def dmBB(self, phi):
        return 0

    def dU(self, phi, Temp):
        return self.dQb(phi) - self.dQw(phi, Temp) - (((self.get_m() * self.get_R() * Temp) / (self.Hubvolumen(phi))) * self.dV(
            phi)) - self.get_spezEnthalpieBB() * self.dmBB(phi)

    def dT(self, Temp, phi):
        return (1 / (self.get_m() * self.get_cv())) * self.dU(phi, Temp)

    def Druck(self, phi, Temp):
        return self.get_m() * self.get_R() * Temp / self.Hubvolumen(phi)

    def solve(self):

        StartZeit = time()

        T = odeint(self.dT, self.T0, self.phiKW)
        self.T = T

        for i in range(self.get_AnzahlStuetzstellen()):
            self.V[i] = self.Hubvolumen(self.phiKW[i])
            self.deltaV[i] = self.dV(self.phiKW[i])
            self.deltaQb[i] = self.dQb(self.phiKW[i])
            self.p[i] = self.Druck(self.phiKW[i], T[i])
            self.deltaU[i] = self.dU(self.phiKW[i], T[i])
            self.deltaQw[i] = self.dQw(self.phiKW[i], T[i])

            self.V_darstellung[i] = self.V[i] * 1000000
            self.p_darstellung[i] = self.p[i] / 100000

        self.Kreisprozessarbeit = self.get_Kreisprozessarbeit()
        self.Wirkunsgrad = self.get_Wirkungsgrad()
        self.pmax = self.get_pmax()
        self.IndizierterMitteldruck = self.get_IndizierterMitteldruck()
        self.Drehmoment = self.get_Drehmoment()
        self.Leistung = self.get_Leistung()

        EndZeit = time()

        self.execTime = EndZeit - StartZeit

        return T, EndZeit - StartZeit

    def printErgebnisse(self):
        print("Wk : {:.0f} J".format(self.get_Kreisprozessarbeit()))
        print("pmi : {:.2f} bar".format(self.get_IndizierterMitteldruck()))
        print("Wirkungsgrad : {:.2f} %".format(self.get_Wirkungsgrad() * 100))
        print("M : {:.0f} Nm bei ".format(self.get_Drehmoment()) + str(self.get_Drehzahl()) + " U/min")
        print("P : {:.0f} PS bei ".format(self.get_Leistung()) + str(self.get_Drehzahl()) + " U/min")
        print("pmax: {:.0f} bar".format(self.get_pmax() / 10 ** 5))
        print("Tmax: {:.0f} K".format(self.get_Tmax()))
        print("Berechnungszeit : {:.2f} s".format(self.execTime))
