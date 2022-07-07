import warnings

R_Molar = 8.31446261815324

Benzin_E5 = {
    "Name": "Benzin E5",
    "Hu": 42200000,
    "Dichte": 755,
    "Mindestluftbedarf": 14.7,
    "Molare Masse": 93.8
}

Benzin_E10 = {
    "Name": "Benzin E10",
    "Hu": 42000000,
    "Dichte": 758.5,
    "Mindestluftbedarf": 14.27,
    "Molare Masse": 93.1
}

E85 = {
    "Name": "E85",
    "Hu": 28800000,
    "Dichte": 785,
    "Mindestluftbedarf": 9.96,
    "Molare Masse": 53.61
}

Ethanol = {
    "Name": "Ethanol",
    "Hu": 26800000,
    "Dichte": 789,
    "Mindestluftbedarf": 9.07,
    "Molare Masse": 46.07
}

Diesel = {
    "Name": "Diesel",
    "Hu": 42800000,
    "Dichte": 840,
    "Mindestluftbedarf": 14.6,
    "Molare Masse": 170
}

Biodiesel = {
    "Name": "Biodiesel",
    "Hu": 37100000,
    "Dichte": 880,
    "Mindestluftbedarf": 12.5,
    "Molare Masse": 296
}

Schweroel = {
    "Name": "Schweroel",
    "Hu": 41300000,
    "Dichte": 950,
    "Mindestluftbedarf": 14.6,
    "Molare Masse": 198
}

Wasserstoff = {
    "Name": "Wasserstoff",
    "Hu": 120000000,
    "Dichte": 0.01,
    "Mindestluftbedarf": 34.2,
    "Molare Masse": 2.02
}

LPG = {
    "Name": "LPG",
    "Hu": 45800000,
    "Dichte": 2.25,
    "Mindestluftbedarf": 15.5,
    "Molare Masse": 44.1
}

CNG = {
    "Name": "CNG",
    "Hu": 47700000,
    "Dichte": 0.83,
    "Mindestluftbedarf": 17.2,
    "Molare Masse": 16.04
}

Pflanzenoel = {
    "Name": "Pflanzenoel",
    "Hu": 37600000,
    "Dichte": 920,
    "Mindestluftbedarf": 12.7,
    "Molare Masse": 338.57
}

Methanol = {
    "Name": "Methanol",
    "Hu": 19700000,
    "Dichte": 795,
    "Mindestluftbedarf": 6.46,
    "Molare Masse": 32.04
}

Kraftstoffe = {
    "Benzin E5": Benzin_E5,
    "Benzin E10": Benzin_E10,
    "E85": E85,
    "Ethanol":Ethanol,
    "Diesel":Diesel,
    "Biodiesel":Biodiesel,
    "Schweroel":Schweroel,
    "Wasserstoff":Wasserstoff,
    "LPG":LPG,
    "CNG":CNG,
    "Pflanzenoel":Pflanzenoel,
    "Methanol":Methanol
}


def get_Kraftstoff(Name):
    try:
        ret = Kraftstoffe[Name]
        ret["Gaskonstante"] = R_Molar / ret["Molare Masse"]
        return ret
    except KeyError:
        warnings.warn("Kraftstoff nicht gefunden")
        return None

