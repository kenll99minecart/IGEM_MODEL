import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint



kad = 0.0001
kb1 = 0.5
kd1 = 0.5
kpt = 0.0015
kd2 = 0.5
kb2 = 0.5
kph = 0.05
kd3 = 0.5
kb3 = 0.5
tf = 8000

kap2 = 0.01
kap = 0.01
kpt2 = 102.1
kbLZ = 0.1
kdLZ = 0.1
kbSH3 = 0.4
kdSH3 = 0.01

EnvZPi = 0
EnvZPRi = 0
EnvZRPi = 0
EnvZRi = 0
OmpRPi = 0
EnvZi = 0.17
OmpRi = 6
Sci = 12
ScEnvZPi = 0
ScEnvZPRi = 0
ScEnvZRPi = 0
ScEnvZRi = 0
ScOmpRPi = 0
ScEnvZi = 0
ScOmpRi = 0
RFPi = 0
GFPi = 0

KC = 20e-3
KF = 1e-3
KF4 = 20e-3
kG = 1
kR = 1
kC = 0.01
kF = 0.01
kdG = 0.001
kdR = 0.001


def wrap(kad, kb1, kd1, kpt, kd2, kb2, kph, kd3, kb3, kap2, kap, kpt2, kbLZ, kdLZ, kbSH3, kdSH3, KC, kG, kC):
    def model(z, t):
        envz = z[0]
        envzp = z[1]
        envzpr = z[2]
        envzrp = z[3]
        envzr = z[4]
        ompr = z[5]
        omprp = z[6]
        sc = z[7]
        scenvz = z[8]
        scenvzp = z[9]
        scenvzpr = z[10]
        scenvzrp = z[11]
        scenvzr = z[12]
        scompr = z[13]
        scomprp = z[14]
        envzd = z[15]

        denvzdt = -kap * envz + kad * envzp + kd2 * envzrp - kb2 * envz * omprp + kd3 * envzr - kb3 * envz * ompr - \
            kbSH3 * sc * envz + kdSH3 * scenvz - kbSH3 * scompr * envz + kdSH3 * scenvzr - kbSH3 * scomprp * envz + \
            kdSH3 * scenvzrp - kbSH3 * envz + kdSH3 * envzd
        denvzddt = kbSH3 * envz - kdSH3 * envzd
        denvzpdt = kap * envz - kad * envzp - kb1 * envzp * (ompr ** 2) + kd1 * envzpr - kbSH3 * sc * envzp + kdSH3 * \
            scenvzp - kbSH3 * scompr * envzp + kdSH3 * scenvzpr
        denvzprdt = kb1 * envzp * (ompr ** 2) - kd1 * envzpr - kpt * envzpr
        denvzrpdt = kpt * envzpr - kd2 * envzrp + kb2 * envz * omprp - kph * envzrp
        denvzrdt = kph * envzrp + kb3 * envz * ompr - kd3 * envzr
        domprdt = -kb1 * envzp * (ompr ** 2) + kd1 * envzpr + kd3 * envzr - kb3 * envz * ompr - kbLZ * sc * ompr + kdLZ * \
                scompr - kbLZ * scenvz * ompr + kdLZ * scenvzr - kbLZ * scenvzp * ompr + kdLZ * scenvzpr
        domprpdt = kd2 * envzrp - kb2 * envz * omprp - kbLZ * sc * omprp + kdLZ * scomprp - kbLZ * scenvz * omprp + kdLZ * \
                scenvzrp
        dscdt = -kbSH3 * sc * envz - kbSH3 * sc * envzp - kbLZ * sc * ompr - kbLZ * sc * omprp + kdSH3 * scenvz + kdSH3 * \
                scenvzp + kdLZ * scompr + kdLZ * scomprp + kG * kC * omprp ** 2 / (omprp ** 2 + KC ** 2)
        dscenvzdt = -kap2 * scenvz + kad * scenvzp + kbSH3 * sc * envz - kdSH3 * scenvz - kbLZ * scenvz * ompr + kdLZ * \
                    scenvzr - kbLZ * scenvz * omprp + kdLZ * scenvzrp
        dscenvzpdt = kap2 * scenvz - kad * scenvzp + kbSH3 * sc * envzp - kdSH3 * scenvzp - kbLZ * scenvzp * ompr + kdLZ * \
                    scenvzpr
        dscenvzprdt = kbLZ * scenvzp * ompr - kdLZ * scenvzpr + kbSH3 * scompr * envzp - kdSH3 * scenvzpr - kpt2 * scenvzpr
        dscenvzrpdt = kbLZ * scenvz * omprp - kdLZ * scenvzrp + kbSH3 * scomprp * envz - kdSH3 * scenvzrp + kpt2 * \
                    scenvzpr - kph * scenvzrp
        dscenvzrdt = kbLZ * scenvz * ompr - kdLZ * scenvzr + kbSH3 * scompr * envz - kdSH3 * scenvzr + kph * scenvzrp
        dscomprdt = kbLZ * sc * ompr - kdLZ * scompr - kbSH3 * scompr * envz + kdSH3 * scenvzr - kbSH3 * scompr * envzp + \
                kdSH3 * scenvzpr
        dscomprpdt = kbLZ * sc * omprp - kdLZ * scomprp - kbSH3 * scomprp * envz + kdSH3 * scenvzrp

        dzdt = [denvzdt, denvzddt, denvzpdt, denvzprdt, denvzrpdt, denvzrdt, domprdt, domprpdt, dscdt, dscenvzdt,
                dscenvzpdt, dscenvzprdt, dscenvzrpdt, dscenvzrdt, dscomprdt, dscomprpdt]
        return dzdt

    z0 = [0.17, 0, 0, 0, 0, 6, 0, 12, 0, 0, 0, 0, 0, 0, 0, 0]
    t = np.linspace(0, 2000, 10000)
    return odeint(model, z0, t)


def main(kad, kb1, kd1, kpt, kd2, kb2, kph, kd3, kb3, kap2, kap, kpt2, kbLZ, kdLZ, kbSH3, kdSH3, KC, kG, kC):
    z = wrap(kad, kb1, kd1, kpt, kd2, kb2, kph, kd3, kb3, kap2, kap, kpt2, kbLZ, kdLZ, kbSH3, kdSH3, KC, kG, kC)
    if True:
        t = np.linspace(0, 2000, 10000)
        plt.plot(t, z[:, 0], 'b-', label='envz')
        plt.plot(t, z[:, 1], 'r-', label='envzp')
        plt.plot(t, z[:, 2], 'b-.', label='envzpr')
        plt.plot(t, z[:, 3], 'r-.', label='envzrp')
        plt.plot(t, z[:, 4], 'b--', label='envzr')
        plt.plot(t, z[:, 5], 'g-', label='ompr')
        plt.plot(t, z[:, 6], 'b-', label='omprp')
        plt.plot(t, z[:, 7], 'r-', label='sc')
        plt.plot(t, z[:, 8], 'b-.', label='scenvz')
        plt.plot(t, z[:, 9], 'r-.', label='scenvzp')
        plt.plot(t, z[:, 10], 'b--', label='scenvzpr')
        plt.plot(t, z[:, 11], 'g-', label='scenvzrp')
        plt.plot(t, z[:, 12], 'b-', label='scenvzr')
        plt.plot(t, z[:, 13], 'r-', label='scompr')
        plt.plot(t, z[:, 14], 'b-.', label='scomprp')
        plt.plot(t, z[:, 15], 'r-.', label='envzd')
        plt.ylabel('concentration')
        plt.xlabel('time')
        plt.legend(loc='best')
        plt.show()
    return z[:, 15][-1]


if __name__ == '__main__':
    r = main(kad, kb1, kd1, kpt, kd2, kb2, kph, kd3, kb3, kap2, kap, kpt2, kbLZ, kdLZ, kbSH3, kdSH3, KC, kG, kC)



