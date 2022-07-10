import shutil


def main(args=None):
    from matplotlib import pyplot as plt
    from Code.Simulationsmodule.Prozessrechnung.Prozessrechnung import Realprozessrechnung  # , Kreisprozessrechnung
    import os

    #OUTPUT-Ordner löschen

    try:
        shutil.rmtree("Output\Bilder")
        os.makedirs("Output\Bilder")
    except OSError as e:
        print(e)
    else:
        print("Verzeichnis erfolgreich geleert")

    def plot(x, y, Name):
        plt.plot(x, y, label=Name)
        plt.legend()
        plt.savefig(f"Output\Bilder\{Name}.png", dpi=250)
        if Name == "p-V":
            plt.show()
        plt.close()

    Model = Realprozessrechnung(Kurbelwinkelaufloesung=1, Kraftstoff="Benzin", isLuftansaugend=False)

    #T, execTime = Model.solve(modus="stat")

    #Model.printErgebnisse()

    y = []
    x = []
    min = 10
    max = 100
    schrittweite = 25
    #for i in range(min, max + 1, schrittweite):
    for krst in ["Benzin","Diesel","Wasserstoff","CNG"]:
        x.append(krst)
        Model = Realprozessrechnung(Kurbelwinkelaufloesung=1, Kraftstoff=krst, isLuftansaugend=False)

        _, zeit = Model.solve()
        y.append(Model.get_Wirkungsgrad())
        #print("Fortschritt : {:.2f} %".format(100 * (i - min) / (max - min)))
        Model.printErgebnisse()
        plt.plot(Model.get_phi_KW(),Model.get_T_array(),label=str(krst))
        #plt.plot(Model.get_V_Darstellung_Array(),Model.get_p_Darstellung_Array(),label=str(krst))
    #plt.plot(x, y, label="")
    plt.legend()
    plt.show()
    exit()

    phi = Model.get_phi_KW()

    plot(phi, Model.get_p_Darstellung_Array(), "p-phi")
    plot(phi, Model.get_T_array(), "T-phi")
    plot(Model.get_V_Darstellung_Array(), Model.get_p_Darstellung_Array(), "p-V")
    plot(Model.get_V_Darstellung_Array(), Model.get_T_array(), "T-V")
    plot(phi, Model.get_dQb_Array(), "dQb-phi")
    plot(phi, Model.get_dU_Array(), "dU-phi")
    plot(phi, Model.get_dQw_Array(), "dQw-phi")
    plot(phi, Model.get_normierter_Summenbrennverlauf_Array(), "Summenbrennverlauf-phi")
    plot(phi,Model.get_alpha_array(),"alpha-phi")

    # plt.xlabel("°KW")
    """plt.legend()
    plt.savefig("Output\Bilder\pyplot.png", dpi=1000)
    plt.show()"""


if __name__ == '__main__':
    main()
