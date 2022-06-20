Moritz = {
    "Name":"Moritz",
    "Wohnung":"Zwickau",
    "Geschlecht":"m",
    "Alter":21
}

Julia = {
    "Name" : "Julia",
    "Wohnung":"Dresden",
    "Geschlecht":"w",
    "Alter":21
}

Schatzi = {
    "Moritz":Moritz,
    "Julia":Julia
}

Person = Schatzi["Moritz"]
Person["Doppel"] = 2*Person["Alter"]
print(Person["Doppel"])