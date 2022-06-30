import matplotlib.pyplot as plt

def main():
    Thermo = ThermodynamischeDaten()
    print(Thermo.get_Anzahl())


class ThermodynamischeDaten(object):

    def __init__(self):

        self.__data = []
        self.__Anzahl = 0
        self.x = []

        with open("thermo30.dat", "r") as file:
            for x in file:
                if x == 'END\n':
                    break
                if x.find('4\n') == 79:
                    self.__Anzahl +=1
                self.__data.append(x)

        self.__data = self.__data[5:]

    def get_data(self,Verbindung = "02"):
        return Verbindung

    def get_Anzahl(self):
        return self.__Anzahl


if __name__ == "__main__":
    main()
