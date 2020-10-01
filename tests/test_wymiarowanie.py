from vector.wymiarowanie import oblicz
from vector.si import *


def test_wymiarowanie():
    oblicz(
        fi_1=10 * mm,
        fi_r=6 * mm,
        h=0.15 * m,
        g_k=(25 * 0.15 + 2.5) * kPa,
        q_k=2 * kPa,
        l_eff=5 * m,
    )
