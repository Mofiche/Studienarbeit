from matplotlib import pyplot as plt
from Code.RPR import Realprozessrechnung

Model = Realprozessrechnung(Kurbelwinkelaufloesung=1)

T, execTime = Model.solve()

Model.printErgebnisse()

"""
y = []
x = []
min = 1
max = 1000
schrittweite = 10
for i in range(min, max+1, schrittweite):
    x.append(i)
    Model.set_m_vibe(i/100)
    Model.solve()
    y.append(Model.get_phi_q_50())
    print("Fortschritt : {:.2f} %".format(100 * (i - min) / (max - min)),Model.get_Leistung())
    #plt.plot(Model.V_darstellung,Model.p_darstellung,label=str(i))


plt.plot(x, y, label="")
"""



# plt.plot(phiKW, V, label="V")
# plt.plot(phiKW,deltaV,label="dV")
# plt.plot(phiKW, p_darstellung, label="p")
# plt.plot(Prozessrechnung.phiKW, T, label="T")
# plt.plot(Prozessrechnung.V_darstellung, Prozessrechnung.p_darstellung, label="p-V")
# plt.plot(phiKW, deltaQb, label="dQb")
# plt.plot(phiKW, deltaU, label="dU")
# plt.plot(phiKW,deltaQw,label="dQw")
# plt.plot(Model.get_phi_KW(),Model.get_normierter_Summenbrennverlauf_Array(),label="QB")

# plt.xlabel("°KW")
plt.legend()
plt.savefig("pyplot.png", dpi=1000)
plt.show()
