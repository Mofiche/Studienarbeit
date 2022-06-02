from matplotlib import pyplot as plt
import numpy as np
from Code import Numerik

Num = Numerik.RungeKutta4

# Alle Einheiten sind insofern nicht anders angegeben SI Einheiten!

T0 = 300  # Starttemperatur 300K
p0 = 100000  # Startdruck 1bar
phiES = 180 + 40  # 40° nach UT 180+40
phiAOE = 540 - 60  # 60° vor UT 540-60
cv = 718
R = 287
b = 0.081
s = 0.0864
r = s / 2
l = 0.144
epsilon = 10
lambdaVerbrennung = 1
Lmin = 14.7
RGA = 0.05

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

phi = []
deltaT = []
deltaV = []
V = []
T = []
p = []

for i in range(phiAOE - phiES + 1):
    phi.append(i + phiES)


def Hubvolumen(phi):
    rad = phi * np.pi / 180
    return Vc + ((Vh) / (2 * r)) * (r * (1 - np.cos(rad)) + l * (1 - np.sqrt(1 - np.power(lambdaPL * np.sin(rad), 2))))


def dV(phi):
    rad = phi * np.pi / 180
    return Vh * (np.sin(rad) * 0.5 + 0.25 * lambdaPL * (
                (np.sin(2 * rad)) / (np.sqrt(1 - np.power(lambdaPL * np.sin(rad), 2)))))


def dQb(phi):
    return 0


def dT(phi, T):
    return (1 / (m * cv)) * (dQb(phi) - (((m * R * T) / (Hubvolumen(phi))) * dV(phi)))


def Druck(phi, T):
    return m * R * T / Hubvolumen(phi)

for i in phi:
    V.append(Hubvolumen((i)))

phiTK, TRK ,dTRK = Num.solve_DGL_RK4(Num,dT,phiES,phiAOE,T0)

#for i in range(len(TRK)):
   # p.append((1/100000)*Druck(phi[i],TRK[i])) #in bar

plt.plot(phi,V, label="V")
plt.plot(phiTK,dTRK,label="dT")
plt.plot(phiTK,TRK,label="T")
plt.xlabel("°KW")
plt.legend()
plt.show()


print(rhoL)

