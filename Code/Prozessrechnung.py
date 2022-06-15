import scipy.misc
from numpy import cos, sin, exp, power, pi, arange, zeros, sqrt
from time import time
from scipy.integrate import odeint, simps


# Alle Einheiten sind insofern nicht anders angegeben SI Einheiten!

class Realprozessrechnung(object):

    def __init__(self, Kurbelwinkelaufloesung=1, Zylinderanzahl=4, isLuftansaugend=False, spezEnthalpieBB=0, phiBD=50,
                 m_vibe=1, Hu=41000000, RGA=0.04, Lmin=14.7, lambdaVerbrennung=1, epsilon=10.3, Pleuellaenge=0.144,
                 Hub=0.0864, Bohrung=0.081, R=287, cv=718, Drehzahl=5800, ASP=0.5, T0=300,
                 p0=100000, phiES=220, phiAOE=480, ZZP=352):

        self.__Tmax = 0
        self.__execTime = 0
        self.__Kreisprozessarbeit = 0
        self.__Wirkunsgrad = 0
        self.__pmax = 0
        self.__IndizierterMitteldruck = 0
        self.__Drehmoment = 0
        self.__Leistung = 0

        self.__Genauigkeit = Kurbelwinkelaufloesung
        self.__ASP_pro_Umdrehung = ASP  # 1: 2-Takt 0.5: 4-Takt
        self.__T0 = T0  # Starttemperatur 300K
        self.__p0 = p0  # Startdruck 1bar
        self.__phiES = phiES  # 40° nach UT 180+40
        self.__phiAOE = phiAOE  # 60° vor UT 540-60
        self.__ZZP = ZZP  # 20° vor OT Zuenden 360-20
        self.__ZV = 0
        self.__Drehzahl = Drehzahl  # in U/min
        self.__cv = cv
        self.__R = R
        self.__Bohrung = Bohrung
        self.__Hub = Hub
        self.__Kurbelradius = self.__Hub / 2
        self.__Pleuellaenge = Pleuellaenge
        self.__epsilon = epsilon
        self.__lambdaVerbrennung = lambdaVerbrennung
        self.__Lmin = Lmin
        self.__RGA = RGA
        self.__Hu = Hu
        self.__m_vibe = m_vibe
        self.__phiBA = self.__ZZP + self.__ZV
        self.__phiBD = phiBD
        self.__spezEnthalpieBB = spezEnthalpieBB
        self.__isLuftansaugend = isLuftansaugend
        self.__Zylinderanzahl = Zylinderanzahl

        self.__phiKW = arange(self.__phiES, self.__phiAOE + self.__Genauigkeit, self.__Genauigkeit)
        self.__T = zeros(self.get_AnzahlStuetzstellen())
        self.__deltaT = zeros(self.get_AnzahlStuetzstellen())
        self.__deltaV = zeros(self.get_AnzahlStuetzstellen())
        self.__deltaW = zeros(self.get_AnzahlStuetzstellen())
        self.__deltaQb = zeros(self.get_AnzahlStuetzstellen())
        self.__deltaQw = zeros(self.get_AnzahlStuetzstellen())
        self.__deltaU = zeros(self.get_AnzahlStuetzstellen())
        self.__V = zeros(self.get_AnzahlStuetzstellen())
        self.__p = zeros(self.get_AnzahlStuetzstellen())
        self.__lambdaVG = zeros(self.get_AnzahlStuetzstellen())
        self.__Brennstoffmasse_phi = zeros(self.get_AnzahlStuetzstellen())
        self.__V_darstellung = zeros(self.get_AnzahlStuetzstellen())
        self.__p_darstellung = zeros(self.get_AnzahlStuetzstellen())
        self.__QB = zeros(self.get_AnzahlStuetzstellen())
        self.__normierter_Summenbrennverlauf = zeros(self.get_AnzahlStuetzstellen())

    def initArrays(self):
        self.__phiKW = arange(self.__phiES, self.__phiAOE + self.__Genauigkeit, self.__Genauigkeit)
        self.__T = zeros(self.get_AnzahlStuetzstellen())
        self.__deltaT = zeros(self.get_AnzahlStuetzstellen())
        self.__deltaV = zeros(self.get_AnzahlStuetzstellen())
        self.__deltaW = zeros(self.get_AnzahlStuetzstellen())
        self.__deltaQb = zeros(self.get_AnzahlStuetzstellen())
        self.__deltaQw = zeros(self.get_AnzahlStuetzstellen())
        self.__deltaU = zeros(self.get_AnzahlStuetzstellen())
        self.__V = zeros(self.get_AnzahlStuetzstellen())
        self.__p = zeros(self.get_AnzahlStuetzstellen())
        self.__lambdaVG = zeros(self.get_AnzahlStuetzstellen())
        self.__Brennstoffmasse_phi = zeros(self.get_AnzahlStuetzstellen())
        self.__V_darstellung = zeros(self.get_AnzahlStuetzstellen())
        self.__p_darstellung = zeros(self.get_AnzahlStuetzstellen())
        self.__QB = zeros(self.get_AnzahlStuetzstellen())
        self.__normierter_Summenbrennverlauf = zeros(self.get_AnzahlStuetzstellen())

    def get_Genauigkeit(self):
        return self.__Genauigkeit

    def set_Genauigkeit(self, value):
        self.__Genauigkeit = value

    def get_ASP_pro_Umdrehung(self):
        return self.__ASP_pro_Umdrehung

    def set_ASP_pro_Umdrehung(self, value):
        self.__ASP_pro_Umdrehung = value

    def get_p0(self):
        return self.__p0

    def set_p0(self, value):
        self.__p0 = value

    def get_T0(self):
        return self.__T0

    def set_T0(self, value):
        self.__T0 = value

    def get_phiES(self):
        return self.__phiES

    def set_phiES(self, value):
        self.__phiES = value

    def get_phiAOE(self):
        return self.__phiAOE

    def set_phiAOE(self, value):
        self.__phiAOE = value

    def get_ZZP(self):
        return self.__ZZP

    def set_ZZP(self, value):
        self.__ZZP = value

    def get_ZV(self):
        return self.__ZV

    def set_ZV(self, value):
        self.__ZV = value

    def get_Drehzahl(self):
        return self.__Drehzahl

    def set_Drehzahl(self, value):
        self.__Drehzahl = value

    def get_cv(self):
        return self.__cv

    def set_cv(self, value):
        self.__cv = value

    def get_R(self):
        return self.__R

    def set_R(self, value):
        self.__R = value

    def get_Bohrung(self):
        return self.__Bohrung

    def set_Bohrung(self, value):
        self.__Bohrung = value

    def get_Hub(self):
        return self.__Hub

    def set_Hub(self, value):
        self.__Hub = value

    def get_Pleuellaenge(self):
        return self.__Pleuellaenge

    def set_Pleuellaenge(self, value):
        self.__Pleuellaenge = value

    def get_epsilon(self):
        return self.__epsilon

    def set_epsilon(self, value):
        self.__epsilon = value

    def get_lambdaVerbrennung(self):
        return self.__lambdaVerbrennung

    def set_lambdaVerbrennung(self, value):
        self.__lambdaVerbrennung = value

    def get_Lmin(self):
        return self.__Lmin

    def set_Lmin(self, value):
        self.__Lmin = value

    def get_RGA(self):
        return self.__RGA

    def set_RGA(self, value):
        self.__RGA = value

    def get_Hu(self):
        return self.__Hu

    def set_Hu(self, value):
        self.__Hu = value

    def get_m_vibe(self):
        return self.__m_vibe

    def set_m_vibe(self, value):
        self.__m_vibe = value

    def get_phiBA(self):
        return self.get_ZZP() + self.get_ZV()

    def get_Kurbelradius(self):
        return self.get_Hub() / 2

    def get_phiBD(self):
        return self.__phiBD

    def set_phiBD(self, value):
        self.__phiBD = value

    def get_spezEnthalpieBB(self):
        return self.__spezEnthalpieBB

    def set_spezEnthalpieBB(self, value):
        self.__spezEnthalpieBB = value

    def get_isLuftansaugend(self):
        return self.__isLuftansaugend

    def set_isLuftansaugend(self, value):
        self.__isLuftansaugend = value

    def get_Zylinderanzahl(self):
        return self.__Zylinderanzahl

    def set_Zylinderanzahl(self, value):
        self.__Zylinderanzahl = value

    def get_AnzahlStuetzstellen(self) -> int:
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
        return simps(self.__p, self.__V)

    def get_Wirkungsgrad(self):
        return self.get_Kreisprozessarbeit() / self.get_Qmax()

    def get_pmax(self):
        return max(self.__p)

    def get_T_array(self):
        return self.__T

    def get_Tmax(self):
        return max(self.get_T_array())[0]

    def get_IndizierterMitteldruck(self):
        return self.get_Kreisprozessarbeit() / (100000 * (self.get_Vmax() - self.get_Vmin()))

    def get_Drehmoment(self):
        return self.get_ASP_pro_Umdrehung() * self.get_IndizierterMitteldruck() * 10 ** 5 * self.get_Hubraum() * \
               self.get_Zylinderanzahl() / (2 * pi)

    def get_Leistung(self):
        return self.get_Drehmoment() * self.get_omega() * 1.36 / 1000

    def get_p_array(self):
        return self.__p

    def get_phi_q_50(self):
        ret = 0

        for i in self.__phiKW:
            if self.get_normierter_Summenbrennverlauf(i) >= 0.5:
                ret = i
                break
        return ret - 360

    def get_phi_KW(self):
        return self.__phiKW

    def get_V_Array(self):
        return self.__V

    def get_V_Darstellung_Array(self):
        return self.__V_darstellung

    def get_p_Darstellung_Array(self):
        return self.__p_darstellung

    def get_dV_Array(self):
        return self.__deltaV

    def get_dQb_Array(self):
        return self.__deltaQb

    def get_dQw_Array(self):
        return self.__deltaQw

    def get_dW_Array(self):
        return self.__deltaW

    def get_dU_Array(self):
        return self.__deltaU

    def get_QB_Array(self):
        return self.__QB

    def get_lambdaVG_Array(self):
        return self.__lambdaVG

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
            ret = (self.get_Qmax() / self.get_phiBD()) * 6.908 * (self.get_m_vibe() + 1) * power(
                (phi - self.get_phiBA()) / self.get_phiBD(),
                self.get_m_vibe()) * exp(
                -6.908 * power((phi - self.get_phiBA()) / self.get_phiBD(), self.get_m_vibe() + 1))
        return ret

    def Qb(self, phi):
        ret = scipy.integrate.quad(self.dQb, self.__phiBA, phi)[0]
        return ret

    def get_normierter_Summenbrennverlauf(self, phi):
        ret = self.Qb(phi) / self.get_Qmax()
        return ret

    def get_normierter_Summenbrennverlauf_Array(self):
        return self.__normierter_Summenbrennverlauf

    def dQw(self, phi, Temp):
        return 0

    def dmBB(self, phi):
        return 0

    def dU(self, phi, Temp):
        return self.dQb(phi) - self.dQw(phi, Temp) - (
                ((self.get_m() * self.get_R() * Temp) / (self.Hubvolumen(phi))) * self.dV(
            phi)) - self.get_spezEnthalpieBB() * self.dmBB(phi)

    def dT(self, Temp, phi):
        return (1 / (self.get_m() * self.get_cv())) * self.dU(phi, Temp)

    def dT_Justi(self, Temp, phi):
        pass

    def Druck(self, phi, Temp):
        return self.get_m() * self.get_R() * Temp / self.Hubvolumen(phi)

    def solve(self):

        StartZeit = time()

        self.initArrays()

        T = odeint(self.dT, self.__T0, self.__phiKW)
        self.__T = T

        for i in range(self.get_AnzahlStuetzstellen()):
            self.__p[i] = self.Druck(self.__phiKW[i], T[i])
            self.__deltaU[i] = self.dU(self.__phiKW[i], T[i])
            self.__deltaQw[i] = self.dQw(self.__phiKW[i], T[i])

        self.__V = [self.Hubvolumen(i) for i in self.__phiKW]
        self.__deltaV = [self.dV(i) for i in self.__phiKW]
        self.__deltaQb = [self.dQb(i) for i in self.__phiKW]
        self.__V_darstellung = [i * 1000000 for i in self.__V]
        self.__QB = [self.Qb(i) for i in self.__phiKW]
        self.__normierter_Summenbrennverlauf = [self.get_normierter_Summenbrennverlauf(i) for i in self.__phiKW]
        self.__p_darstellung = [i / 100000 for i in self.__p]

        self.__Kreisprozessarbeit = self.get_Kreisprozessarbeit()
        self.__Wirkunsgrad = self.get_Wirkungsgrad()
        self.__pmax = self.get_pmax()
        self.__IndizierterMitteldruck = self.get_IndizierterMitteldruck()
        self.__Drehmoment = self.get_Drehmoment()
        self.__Leistung = self.get_Leistung()

        EndZeit = time()

        self.__execTime = EndZeit - StartZeit

        return T, EndZeit - StartZeit

    def writeTXT(self, dateiname="Ergebnisse_TXT"):
        file = open(dateiname + ".txt", "w")
        file.write("Winkel Volumen Druck Temperatur \n")
        for i in range(len(self.get_phi_KW())):
            string = f"{self.__phiKW[i]} {self.__V_darstellung[i]} {self.__p_darstellung[i]} {self.__T[i][0]}\n"
            file.write(string)
        file.close()

    def writeCSV(self, dateiname="Ergebnisse_CSV"):
        import csv
        file = open(dateiname + ".csv", "w")
        writer = csv.writer(file)
        writer.writerow("Winkel Volumen Druck Temperatur \n")
        writer.writerows([self.__phiKW, self.__V_darstellung, self.__p_darstellung, self.__T])
        file.close()

    def writeXLSX(self, dateiname="Ergebnisse"):
        try:
            import xlsxwriter
            workbook = xlsxwriter.Workbook(dateiname + ".xlsx")
            worksheet = workbook.add_worksheet()
            worksheet.write("A1", "Winkel")
            worksheet.write("B1", "Volumen")
            worksheet.write("C1", "Druck")
            worksheet.write("D1", "Temperatur")
            for j in range(len(self.__phiKW)):
                worksheet.write("A" + str(j + 2), self.__phiKW[j])
                worksheet.write("B" + str(j + 2), self.__V_darstellung[j])
                worksheet.write("C" + str(j + 2), self.__p_darstellung[j])
                worksheet.write("D" + str(j + 2), self.__T[j])
            chart = workbook.add_chart({'type': 'line'})
            chart.add_series({'categories':f'=Sheet1!A2:A{len(self.__phiKW)+2}','values': f'=Sheet1!C2:C{len(self.__phiKW)+2}','smooth':'True','name':'Druck'})
            chart.set_x_axis({'name':'Kurbelwinkel in °'})
            chart.set_y_axis({'name':'Druck in bar'})
            chart.set_title({'name':'Druck-Kurbelwinkel Diagramm'})
            worksheet.insert_chart('G3', chart)
            workbook.close()
        except:

            from tkinter import messagebox

            messagebox.showwarning("Dateifehler", "Bitte die Exceldatei schließen.")
            self.writeXLSX(dateiname)

    def printErgebnisse(self):
        print("Wk : {:.0f} J".format(self.get_Kreisprozessarbeit()))
        print("pmi : {:.2f} bar".format(self.get_IndizierterMitteldruck()))
        print("Wirkungsgrad : {:.2f} %".format(self.get_Wirkungsgrad() * 100))
        print("M : {:.0f} Nm bei ".format(self.get_Drehmoment()) + str(self.get_Drehzahl()) + " U/min")
        print("P : {:.0f} PS bei ".format(self.get_Leistung()) + str(self.get_Drehzahl()) + " U/min")
        print("pmax: {:.0f} bar".format(self.get_pmax() / 10 ** 5))
        print("Tmax: {:.0f} K".format(self.get_Tmax()))
        print("phi_q_50 liegt bei {:.1f} °KWnOT".format(self.get_phi_q_50()))
        print("Berechnungszeit : {:.2f} s".format(self.__execTime))
