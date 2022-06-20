def main(args=None):
    from matplotlib import pyplot as plt
    from Code.Prozessrechnung import Realprozessrechnung  # , Kreisprozessrechnung

    Model = Realprozessrechnung(Kurbelwinkelaufloesung=1, Kraftstoff="E100")
    # Seiliger = Kreisprozessrechnung()

    T, execTime = Model.solve(modus="stat")

    Model.printErgebnisse()

    """
    y = []
    x = []
    min = 1
    max = 5
    schrittweite = 1
    for i in range(min, max + 1, schrittweite):
        #x.append(i)
        Model.set_m_vibe(i)
        _, zeit = Model.solve()
        x.append(180+i)
        y.append(Model.get_Leistung())
        print("Fortschritt : {:.2f} %".format(100 * (i - min) / (max - min)), Model.get_Leistung())
        plt.plot(Model.get_V_Darstellung_Array(),Model.get_p_Darstellung_Array(),label=str(i))
    #plt.plot(x, y, label="")
    """

    # plt.plot(__phiKW, __V, label="V")
    # plt.plot(__phiKW,__deltaV,label="dV")
    # plt.plot(__phiKW, __p_darstellung, label="p")
    # plt.plot(Prozessrechnung.__phiKW, __T, label="T")
    plt.plot(Model.get_V_Darstellung_Array(), Model.get_p_Darstellung_Array(), label="p-V")
    #plt.plot(Model.get_phi_KW(),Model.get_lambdaVG_Array())
    # plt.plot(__phiKW, __deltaQb, label="dQb")
    # plt.plot(__phiKW, __deltaU, label="dU")
    # plt.plot(__phiKW,__deltaQw,label="dQw")
    # plt.plot(Model.get_phi_KW(),Model.get_normierter_Summenbrennverlauf_Array(),label="__QB")

    # plt.xlabel("°KW")
    plt.legend()
    plt.savefig("pyplot.png", dpi=1000)
    plt.show()


if __name__ == '__main__':
    main()
