def main(args=None):
    from matplotlib import pyplot as plt
    from Code.Simulationsmodule.Prozessrechnung.Prozessrechnung import Realprozessrechnung  # , Kreisprozessrechnung

    Model = Realprozessrechnung(Kurbelwinkelaufloesung=1, Kraftstoff="Benzin E5",isLuftansaugend=True)

    T, execTime = Model.solve(modus="stat")

    Model.printErgebnisse()


    """y = []
    x = []
    min = 8
    max = 50
    schrittweite = 1
    for i in range(min, max + 1, schrittweite):
        x.append(i/10)
        Model.set_lambdaVerbrennung(i/10)
        _, zeit = Model.solve()
        #x.append(Model.get_Leistung())
        y.append(Model.get_Leistung())
        print("Fortschritt : {:.2f} %".format(100 * (i - min) / (max - min)))
        #plt.plot(Model.get_V_Darstellung_Array(),Model.get_p_Darstellung_Array(),label=str(i))
    plt.plot(x, y, label="")"""


    # plt.plot(Model.get_phi_KW(), Model.get_V_Darstellung_Array(), label="V")
    # plt.plot(__phiKW,__deltaV,label="dV")
    #plt.plot(Model.get_phi_KW(), Model.get_p_Darstellung_Array(), label="p")
    #plt.plot(Model.get_phi_KW(), Model.get_T_array(), label="T")
    # plt.plot(Model.get_V_Darstellung_Array(), Model.get_p_Darstellung_Array(), label="p-V")
    # plt.plot(Model.get_phi_KW(),Model.get_lambdaVG_Array())
    plt.plot(Model.get_phi_KW(), Model.get_dQb_Array(), label="dQb")
    plt.plot(Model.get_phi_KW(), Model.get_dU_Array(), label="dU")
    # plt.plot(Model.get_phi_KW(),__deltaQw,label="dQw")
    #plt.plot(Model.get_phi_KW(),Model.get_normierter_Summenbrennverlauf_Array(),label="__QB")

    # plt.xlabel("°KW")
    plt.legend()
    plt.savefig("Output\pyplot.png", dpi=1000)
    plt.show()


if __name__ == '__main__':
    main()
