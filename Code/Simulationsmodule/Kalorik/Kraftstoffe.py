import warnings

R_Molar = 8.31446261815324

Benzin = {
    "Name": "Benzin",
    "Hu": 42200000,
    "Dichte": 755,
    "Mindestluftbedarf": 14.7,
    "Molare Masse": 93.8
}
Diesel = {
    "Name": "Diesel",
    "Hu": 42800000,
    "Dichte": 840,
    "Mindestluftbedarf": 14.6,
    "Molare Masse": 170
}
Wasserstoff = {
    "Name": "Wasserstoff",
    "Hu": 120000000,
    "Dichte": 0.01,
    "Mindestluftbedarf": 34.2,
    "Molare Masse": 2.02
}
CNG = {
    "Name": "CNG",
    "Hu": 47700000,
    "Dichte": 0.83,
    "Mindestluftbedarf": 17.2,
    "Molare Masse": 16.04
}
Kraftstoffe = {
    "Benzin": Benzin,
    "Diesel":Diesel,
    "Wasserstoff":Wasserstoff,
    "CNG":CNG,

}

def get_Kraftstoff(Name):
    try:
        ret = Kraftstoffe[Name]
        ret["Gaskonstante"] = R_Molar / ret["Molare Masse"]
        return ret
    except KeyError:
        warnings.warn("Kraftstoff nicht gefunden")
        return None

