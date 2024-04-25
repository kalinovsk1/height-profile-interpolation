import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df_everest = pd.read_csv('data/MountEverest.csv')
df_gdansk = pd.read_csv('data/SpacerniakGdansk.csv')
df_kolorado = pd.read_csv('data/WielkiKanionKolorado.csv')

data_everest = []
data_gdansk = []
data_kolorado = []

for row in df_everest.values:
    data_everest.append(row)

for row in df_gdansk.values:
    data_gdansk.append(row)

for row in df_kolorado.values:
    data_kolorado.append(row)


def stworz_rowne_wezly(data, ilosc_wezlow):
    indeksy = np.linspace(0, len(data) - 1, num=ilosc_wezlow, dtype=int)
    wezly = [data[i] for i in indeksy]
    return indeksy, wezly


def stworz_random_wezly(data, ilosc_wezlow):
    zakres = np.arange(1, len(data)-1)
    indeksy = np.random.choice(zakres, size=ilosc_wezlow-2, replace=False)
    indeksy = np.append(indeksy, [0, len(data)-1])
    indeksy = sorted(indeksy)
    wezly = [data[i] for i in indeksy]
    return indeksy, wezly


def bazowa_funkcja(i, x, wezly):
    licznik = 1
    mianownik = 1
    for j in range(len(wezly)):
        if j != i:
            licznik *= x - wezly[j][0]
            mianownik *= wezly[i][0] - wezly[j][0]
    wynik = licznik/mianownik
    return wynik


def Langrange(x, wezly):
    wynik = 0
    for i in range(len(wezly)):
        wynik += wezly[i][1] * bazowa_funkcja(i, x, wezly)
    return wynik


def intepolowane_wartosci_Lagrange(ilosc_wezlow, data, typ_wezlow):
    if typ_wezlow == "rowne":
        wezly_indeks, wezly = stworz_rowne_wezly(data, ilosc_wezlow)
    else:
        wezly_indeks, wezly = stworz_random_wezly(data, ilosc_wezlow)
    lagrange = []
    for i in range(len(data)):
        if i not in wezly_indeks:
            lagrange.append([data[i][0], Langrange(data[i][0], wezly)])
        else:
            lagrange.append(data[i])
    return lagrange, wezly


def wykres(data, interpolacja_data, wezly_interpolacja, tytul, miejsce, typ_wezlow):
    data_odleglosc = [row[0] for row in data]
    data_wysokosc = [row[1] for row in data]
    interpolacja_data_odleglosc = [row[0] for row in interpolacja_data]
    interpolacja_data_wysokosc = [row[1] for row in interpolacja_data]
    wezly_data_odleglosc = [row[0] for row in wezly_interpolacja]
    wezly_data_wysokosc = [row[1] for row in wezly_interpolacja]
    plt.plot(data_odleglosc, data_wysokosc, label='rzeczywisty profil')
    plt.plot(interpolacja_data_odleglosc, interpolacja_data_wysokosc, 'r--', label='interpolowany profil')
    plt.plot(wezly_data_odleglosc, wezly_data_wysokosc, 'k.', label='węzły')
    plt.title(f'{tytul} dla {len(wezly_interpolacja)} {typ_wezlow} węzłów - {miejsce}')
    plt.xlabel("odległość")
    plt.ylabel("wysokość")
    plt.legend()
    plt.grid()
    plt.xlim([0, data[len(data) - 1][0]])
    plt.savefig(f'results/{tytul} {miejsce} {len(wezly_interpolacja)}')
    plt.show()


def wielomian(a, b, c, d, x, x_num):
    return a + b*(x - x_num) + c*(x - x_num)**2 + d*(x - x_num)**3


def stworz_uklad_rownan(ilosc_wezlow, data, typ_wezlow):
    if typ_wezlow == "rowne":
        wezly_indeks, wezly = stworz_rowne_wezly(data, ilosc_wezlow)
    else:
        wezly_indeks, wezly = stworz_random_wezly(data, ilosc_wezlow)

    wielkosc_macierzy = (len(wezly)-1) * 4
    A = np.zeros((wielkosc_macierzy, wielkosc_macierzy))
    a_index = 0
    param_start_index = 0
    S_numer = 0
    for i in range(wielkosc_macierzy//2):
        if i % 2 == 0:
            for j in range(wielkosc_macierzy):
                if j != a_index:
                    A[i][j] = 0
                else:
                    A[i][j] = 1
            a_index += 4
            S_numer += 1
        else:
            for j in range(0, param_start_index):
                A[i][j] = 0
            A[i][param_start_index] = 1
            A[i][param_start_index + 1] = wezly[S_numer][0] - wezly[S_numer-1][0]
            A[i][param_start_index + 2] = (wezly[S_numer][0] - wezly[S_numer - 1][0])**2
            A[i][param_start_index + 3] = (wezly[S_numer][0] - wezly[S_numer - 1][0])**3
            for j in range(param_start_index+4, wielkosc_macierzy):
                A[i][j] = 0
            param_start_index += 4
    pierw_pochodna_index = 1
    druga_pochodna_index = 2
    S_numer = 1
    for i in range(wielkosc_macierzy//2, wielkosc_macierzy-2):
        if i % 2 == 0:
            A[i][pierw_pochodna_index] = 1
            A[i][pierw_pochodna_index + 1] = 2 * (wezly[S_numer][0] - wezly[S_numer-1][0])
            A[i][pierw_pochodna_index + 2] = 3 * (wezly[S_numer][0] - wezly[S_numer - 1][0])**2
            A[i][pierw_pochodna_index + 4] = -1
            pierw_pochodna_index += 4
        else:
            A[i][druga_pochodna_index] = 2
            A[i][druga_pochodna_index + 1] = 6 * (wezly[S_numer][0] - wezly[S_numer-1][0])
            A[i][druga_pochodna_index + 4] = -2
            S_numer += 1
            druga_pochodna_index += 4
    A[wielkosc_macierzy-2][2] = 2
    A[wielkosc_macierzy-1][wielkosc_macierzy-2] = 2
    A[wielkosc_macierzy - 1][wielkosc_macierzy - 1] = 6 * (wezly[len(wezly)-1][0] - wezly[len(wezly)-2][0])

    b = np.zeros((wielkosc_macierzy, 1))
    b_index = 0
    for i in range(wielkosc_macierzy//2):
        if i % 2 == 0:
            b[i][0] = wezly[b_index][1]
            b_index += 1
        else:
            b[i][0] = wezly[b_index][1]

    return A, b, wezly_indeks, wezly


def funkcje_sklejane(ilosc_wezlow, data, typ_wezlow):
    A, b, wezly_indeks, wezly = stworz_uklad_rownan(ilosc_wezlow, data, typ_wezlow)
    x = np.linalg.solve(A, b)
    sklejane = []
    nr_wspolczynniki = -4
    nr_funkcji = -1
    licznik_funkcji = ilosc_wezlow-1
    for i in range(len(data)):
        if i in wezly_indeks and licznik_funkcji > 0:
            nr_wspolczynniki += 4
            nr_funkcji += 1
            licznik_funkcji -= 1
            sklejane.append(data[i])
        else:
            wartosc = wielomian(x[nr_wspolczynniki][0], x[nr_wspolczynniki+1][0], x[nr_wspolczynniki+2][0],
                                x[nr_wspolczynniki+3][0], data[i][0], wezly[nr_funkcji][0])
            sklejane.append([data[i][0], wartosc])
    return sklejane, wezly


wezly = [3, 5, 15, 30]

for i in range(len(wezly)):
    lagrange_everest, wezly_everest = intepolowane_wartosci_Lagrange(wezly[i], data_everest, "rowne")
    lagrange_gdansk, wezly_gdansk = intepolowane_wartosci_Lagrange(wezly[i], data_gdansk, "rowne")
    lagrange_kolorado, wezly_kolorado = intepolowane_wartosci_Lagrange(wezly[i], data_kolorado, "rowne")

    wykres(data_everest, lagrange_everest, wezly_everest, "Lagrange", "Mount Everest", "równoodległych")
    wykres(data_gdansk, lagrange_gdansk, wezly_gdansk, "Lagrange", "Spacerniak Gdańsk", "równoodległych")
    wykres(data_kolorado, lagrange_kolorado, wezly_kolorado, "Lagrange", "Kanion Kolorado", "równoodległych")

    sklejane_everest, wezly_everest = funkcje_sklejane(wezly[i], data_everest, "rowne")
    sklejane_gdansk, wezly_gdansk = funkcje_sklejane(wezly[i], data_gdansk, "rowne")
    sklejane_kolorado, wezly_kolorado = funkcje_sklejane(wezly[i], data_kolorado, "rowne")

    wykres(data_everest, sklejane_everest, wezly_everest, "funkcje sklejane", "Mount Everest", "równoodległych")
    wykres(data_gdansk, sklejane_gdansk, wezly_gdansk, "funkcje sklejane", "Spacerniak Gdańsk", "równoodległych")
    wykres(data_kolorado, sklejane_kolorado, wezly_kolorado, "funkcje sklejane", "Kanion Kolorado", "równoodległych")

wezly = [5, 15]

for i in range(len(wezly)):
    lagrange_everest, wezly_everest = intepolowane_wartosci_Lagrange(wezly[i], data_everest, "random")
    lagrange_gdansk, wezly_gdansk = intepolowane_wartosci_Lagrange(wezly[i], data_gdansk, "random")
    lagrange_kolorado, wezly_kolorado = intepolowane_wartosci_Lagrange(wezly[i], data_kolorado, "random")

    wykres(data_everest, lagrange_everest, wezly_everest, "Lagrange", "Mount Everest", "losowych")
    wykres(data_gdansk, lagrange_gdansk, wezly_gdansk, "Lagrange", "Spacerniak Gdańsk", "losowych")
    wykres(data_kolorado, lagrange_kolorado, wezly_kolorado, "Lagrange", "Kanion Kolorado", "losowych")

    sklejane_everest, wezly_everest = funkcje_sklejane(wezly[i], data_everest, "random")
    sklejane_gdansk, wezly_gdansk = funkcje_sklejane(wezly[i], data_gdansk, "random")
    sklejane_kolorado, wezly_kolorado = funkcje_sklejane(wezly[i], data_kolorado, "random")

    wykres(data_everest, sklejane_everest, wezly_everest, "funkcje sklejane", "Mount Everest", "losowych")
    wykres(data_gdansk, sklejane_gdansk, wezly_gdansk, "funkcje sklejane", "Spacerniak Gdańsk", "losowych")
    wykres(data_kolorado, sklejane_kolorado, wezly_kolorado, "funkcje sklejane", "Kanion Kolorado", "losowych")
