import math
import numpy as np
import matplotlib.pyplot as plt

temp = []
basse = []
moy = []
haute = []
basseOrMoy = []
chauffe = []
constante = []

# Retourne pour une valeur temp l appartenance a la classe basse temperature
def basse_temp(temp):
    # SI La temperature est inferieur a 10
    if(temp <= 10):
        # ALORS la temperature est basse
        return 1
    # SI la temperature est superieur a 20
    if(temp >= 20):
        # ALORS la temperature n'est pas basse
        return 0
    # SI la temperature est compris entre 10 et 20
    # ALORS son appartenance a temperature basse est comprise entre 0 et 1
    return  float(20 - temp)/float(20 - 10)

# Retourne pour une valeur temp l'appartenance a la classe haute temperature
def haute_temp(temp):
    if(temp <= 20):
        return 0
    if(temp >= 30):
        return 1
    return (float(temp - 20) / float(30 - 20))

# Retourne pour une valeur temp l'appartenance a la classe temperature moyenne
def moy_temp(temp):
    a = 10
    m = 20
    b = 30
    return max(min(((float(temp - a) / float(m - a))),(float(b - temp)/float(b - m))),0)

def const(temp):
    return temp

def fuzzy_chauffe(temp):
    if(temp <= 7):
        return 0
    if(temp >= 10):
        return 1
    return  float(temp - 7)/float(10 - 7)

def fuzzy_min(ens1, ens2):
    min_t = []
    for i in range(0,len(ens1)):
        min_t.append(min(ens1[i],ens2[i]))
    return min_t

def fuzzy_max(ens1, ens2):
    max_t = []
    for i in range(0,len(ens1)):
        max_t.append(max(ens1[i],ens2[i]))
    return max_t


# Pour 400 valeur de 0 a 40
for i in range(0,4000):
    # De 0 a 40C
    temp.append(i/100)
    # basse contient les degres d'appartenance de la classe basse
    # de chaque temperature
    basse.append(basse_temp(i/100))
    moy.append(moy_temp(i/100))
    haute.append(haute_temp(i/100))
    chauffe.append(fuzzy_chauffe(i/100))
    basseOrMoy.append(max(basse_temp(i/100),moy_temp(i/100)))

for i in range(0,4000):
    constante.append(basse_temp(12))



temperature = 16
print("Pour 16C : ")
print("Appartenance basse :", basse_temp(temperature))
print("Appartenance moyenne :", moy_temp(temperature))
print("Appartenance haute :", haute_temp(temperature))

#plt.plot(temp, constante, 'b')
plt.plot(temp, chauffe, 'r')
#plt.plot(temp, fuzzy_min(constante,chauffe), 'r')
#plt.plot(temp, fuzzy_max(basse,fuzzy_max(haute,moy)), 'r')

#plt.plot(temp, fuzzy_min(constante,chauffe), 'r')
plt.grid()

plt.show()
