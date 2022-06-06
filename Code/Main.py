from matplotlib import pyplot as plt
from Code.RPR import Realprozessrechnung

Model = Realprozessrechnung(Kurbelwinkelaufloesung=1)

T, execTime = Model.solve()

Model.printErgebnisse()

y = []
x = []
min = 20
max = 200
schrittweite = 1
for i in range(min, max, schrittweite):
    x.append(i)
    Model.set_epsilon(i/10)
    Model.solve()
    y.append(Model.get_Tmax())
    print("Fortschritt : {:.2f} %".format(100 * (i - min) / (max - min)))

plt.plot(x, y, label="epsilon-eta")

# plt.plot(phiKW, V, label="V")
# plt.plot(phiKW,deltaV,label="dV")
# plt.plot(phiKW, p_darstellung, label="p")
# plt.plot(Prozessrechnung.phiKW, T, label="T")
# plt.plot(Prozessrechnung.V_darstellung, Prozessrechnung.p_darstellung, label="p-V")
# plt.plot(phiKW, deltaQb, label="dQb")
# plt.plot(phiKW, deltaU, label="dU")
# plt.plot(phiKW,deltaQw,label="dQw")

# plt.xlabel("°KW")
plt.legend()
plt.savefig("pyplot.png", dpi=1000)
plt.show()
