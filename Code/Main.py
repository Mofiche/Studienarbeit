from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import odeint
import time

# Alle Einheiten sind insofern nicht anders angegeben SI Einheiten!

Genauigkeit = 0.01  # Aufloesung pro °KW

T0 = 300  # Starttemperatur 300K
p0 = 100000  # Startdruck 1bar
phiES = 180 + 40  # 40° nach UT 180+40
phiAOE = 540 - 60  # 60° vor UT 540-60
ZZP = 360 - 20  # 20° vor OT Zuenden 360-20
cv = 718
R = 287
b = 0.081
s = 0.0864
r = s / 2
l = 0.144
epsilon = 10.3
lambdaVerbrennung = 1
Lmin = 14.7
RGA = 0.05
Hu = 41000000
m_vibe = 1
phiBA = ZZP
phiBD = 40

AnzahlStuetzstellen = int((phiAOE - phiES) / Genauigkeit)+1
lambdaPL = r / l
Vh = np.pi * 0.25 * b * b * s
Vc = Vh / (epsilon - 1)
Vmin = Vc
Vmax = Vc + Vh
rhoL = p0 / (R * T0)
mL = Vh * rhoL
mB = mL / (lambdaVerbrennung * Lmin)
mRGA = (RGA / (1 - RGA)) * (mL + mB)
m = mL + mB + mRGA
Qmax = mB * Hu

# phiKW = np.linspace(phiES, phiAOE, AnzahlStuetzstellen)
phiKW = np.arange(phiES, phiAOE+Genauigkeit, Genauigkeit)
deltaT = np.zeros(AnzahlStuetzstellen)
deltaV = np.zeros(AnzahlStuetzstellen)
deltaW = np.zeros(AnzahlStuetzstellen)
deltaQb = np.zeros(AnzahlStuetzstellen)
V = np.zeros(AnzahlStuetzstellen)
p = np.zeros(AnzahlStuetzstellen)


def Hubvolumen(phi):
    rad = phi * np.pi / 180
    return Vc + (Vh / (2 * r)) * (r * (1 - np.cos(rad)) + l * (1 - np.sqrt(1 - np.power(lambdaPL * np.sin(rad), 2))))


def dV(phi):
    rad = phi * np.pi / 180
    return Vh * (np.sin(rad) * 0.5 + 0.25 * lambdaPL * (
            (np.sin(2 * rad)) / (np.sqrt(1 - np.power(lambdaPL * np.sin(rad), 2)))))


def dQb(phi):
    if phi < phiBA:
        ret = 0
    else:
        ret = (Qmax / phiBD) * 6.908 * (m_vibe + 1) * np.power((phi - phiBA) / phiBD, m_vibe) * np.exp(
            -6.908 * np.power((phi - phiBA) / phiBD, m_vibe + 1))
    # return ret
    return 0


def dT(Temp, phi):
    return (1 / (m * cv)) * (dQb(phi) - (((m * R * Temp) / (Hubvolumen(phi))) * dV(phi)))


def Druck(phi, Temp):
    return m * R * Temp / Hubvolumen(phi)

def Volumenarbeit(Druck,Volumen):
    W = 0
    for i in range(len(Druck)-1):
        deltaW = (Druck[i]+(Druck[i+1]-Druck[i])/2)*(Volumen[i+1]-Volumen[i])
        W += deltaW
    return W


StartZeit = time.time()

T = odeint(dT, T0, phiKW)

for i in range(AnzahlStuetzstellen):
    V[i] = Hubvolumen(phiKW[i])
    deltaV[i] = dV(phiKW[i])
    deltaQb[i] = dQb(phiKW[i])
    p[i] = Druck(phiKW[i], T[i])


EndZeit = time.time()

print(Volumenarbeit(p,V))

print("Berechnungszeit : {:.2f} s".format(EndZeit - StartZeit))

plt.plot(phiKW, V, label="V")
# plt.plot(phiKW, p, label="p")
plt.plot(phiKW, T, label="T")
# plt.plot(V,p, label="p-V)
plt.xlabel("°KW")
plt.legend()
plt.show()
