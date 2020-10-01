"""Obliczanie zbrojenia zgodnie z PN-EN 1992-1-1 wg algorytmów z książki M. Knauffa."""

from vector.si import *

# Beton konstrukcyjny klasy C20/25 (B25)
f_ck = 20 * MPa
f_cd = f_ck / 1.4
f_cm = (f_ck / MPa + 8) * MPa
f_ctm = 0.3 * (f_ck / MPa) ** (2 / 3) * MPa
E_cm = 22 * (0.1 * f_cm / MPa) ** 0.3 * GPa
fi_t_t0 = 3.2
E_c_eff = E_cm / (1 + fi_t_t0)

# Stal zbrojeniowa klasy C (A-IIIN) – gatunek B500SP EPSTAL
f_yk = 500 * MPa
f_yd = f_yk / 1.15
E_s = 200 * GPa

c_min_dur = 15 * mm  # dla klasy konstrukcji S4 i klasy ekspozycji XC1 (tab. 4.4N normy)
delta_c_dev = 5 * mm  # odchyłka otulenia 0 ÷ 10 mm zgodnie z (p. 4.4.1.3 normy)

b = 60 * cm  # szerokość standardowej płyty


def oblicz(fi_1, fi_r, h, g_k, q_k, l_eff, n=None):

    # Otulina i wysokość użyteczna
    c_min_b = fi_1
    c_min = max(c_min_b, c_min_dur, 10 * mm)
    c_nom = c_min + delta_c_dev
    a_1 = c_nom + fi_r + 0.5 * fi_1
    d = h - a_1

    # Statyka
    p = max(g_k * 1.35 + q_k * 1.5 * 0.7, g_k * 1.35 * 0.85 + q_k * 1.5)
    p_q = g_k + q_k * 0.3
    M_Ed = 0.125 * p * l_eff ** 2
    M_Ed_q = 0.125 * p_q * l_eff ** 2
    V_Ed = 0.5 * p * l_eff

    # Zginanie
    mi = M_Ed / (b * d ** 2 * f_cd)
    if mi > 0.371:
        raise ValueError("Przekroczno wartość graniczną mi.")
    omega = 0.9731 - (0.9469 - 1.946 * mi) ** 0.5
    A_s1_req = max(
        omega * b * d * (f_cd / f_yd), 0.26 * f_ctm / f_yk * b * d, 0.0013 * b * d,
    )
    A_s = 3.14159 * (0.5 * fi_1) ** 2
    if not n:
        n = int(A_s1_req // A_s) + (A_s1_req % A_s > 0)
    A_s1_prov = n * A_s
    if A_s1_prov > 0.5 * f_cd / f_yd * b * d or A_s1_prov > 0.04 * b * d:
        raise ValueError("Przekroczono graniczny stopień zbrojenia przekroju.")
    ro_1 = A_s1_prov / (b * d)
    M_Rd = A_s1_prov * f_yd * (d - 0.5138 * ((A_s1_prov * f_yd) / (b * f_cd)))
    print(M_Ed, M_Rd)

    # Ścinanie
    k = min(1 + (20 / (d / mm)) ** 0.5, 2.0)
    v_Rd_c = (0.18 / 1.4) * k * (100 * ro_1 * (f_ck / MPa)) ** (1 / 3) * MPa
    v_min = 0.035 * (k ** 3 * f_ck / MPa) ** 0.5 * MPa
    v_Rd_c = max(v_Rd_c, v_min)
    V_Rd_c = v_Rd_c * b * d
    print(V_Ed, V_Rd_c)

    # Ugięcia
    alfa_e = E_s / E_c_eff
    y = (alfa_e - 1) * A_s1_prov * (d - h / 2) / (b * h + (alfa_e - 1) * A_s1_prov)
    J = b * h ** 3 / 12
    J_I = J + b * h * y ** 2 + (alfa_e - 1) * A_s1_prov * (d - h / 2 - y) ** 2
    alfa_I = 5 / 48 * (M_Ed_q * l_eff ** 2) / (E_c_eff * J_I)
    M_cr = f_ctm * J_I / (h / 2 - y)
    kd = ((2 * ro_1 * alfa_e + (ro_1 * alfa_e) ** 2) ** 0.5 - ro_1 * alfa_e) * d
    J_II = (1 / 3) * b * kd ** 3 + alfa_e * A_s1_prov * (d - kd) ** 2
    J_e = (M_cr / M_Ed_q) ** 3 * J_I + (1 - (M_cr / M_Ed_q) ** 3) * J_II
    alfa_II = (J_I / J_e) * alfa_I
    dzeta = 1 - 0.5 * (M_cr / M_Ed_q) ** 2
    alfa_0 = max(l_eff / 250, 20 * mm)
    alfa = (dzeta * alfa_II) + ((1 - dzeta) * alfa_I) - alfa_0
    print(alfa, l_eff / 250)
