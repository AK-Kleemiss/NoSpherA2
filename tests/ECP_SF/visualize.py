import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import scipy.special as sp

fPI = 4 * np.pi


class ECP_Gaussian:
    def __init__(
        self,
        other=None,
        radial_exp=0,
        angular_momentum=0,
        exponent=1.0,
        coefficient=1.0,
    ):
        if other is not None:
            # Copy properties from the other instance
            self.n = other.n
            self.l = other.l
            self.e = other.e
            self.c = other.c
        else:
            # Assign properties from arguments
            self.n = radial_exp - 2
            self.l = angular_momentum
            self.e = exponent
            self.c = coefficient


class ECP:
    def __init__(self, Z=0):
        self.Z = Z
        self.ecp_gaussians = []
        self.N = 0
        self.L = 0
        self.l_starts = [0]

    def add_ECP_Gaussian(self, other):
        self.ecp_gaussians.append(other)
        self.N += 1
        if other.l > self.L:
            self.L = other.l
            self.l_starts.append(self.N - 1)
        elif other.l < self.L:
            raise ValueError("l is decreasing")
        else:
            pass

    def evaluate(self, r, l):
        res = 0
        counter = 0
        while self.l_starts[l] + counter < self.N:
            i = self.l_starts[l] + counter
            if l + 1 >= len(self.l_starts):
                pass
            elif i >= self.l_starts[l + 1]:
                break
            p = self.ecp_gaussians[i].n
            if p <= -1:
                p = 20 - p
            res += self.ecp_gaussians[i].c * np.exp(-self.ecp_gaussians[i].e * r**2) * r ** self.ecp_gaussians[i].n
            counter += 1
        return res


Rb_def2 = ECP(37)
Rb_def2.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=5.036551, coefficient=89.500198))
Rb_def2.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=1.970849, coefficient=0.493761))
Rb_def2.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=3.843114, coefficient=12.3169))
Rb_def2.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=4.258341, coefficient=58.568974))
Rb_def2.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=1.470709, coefficient=0.431791))
Rb_def2.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=3.843114, coefficient=12.3169))
Rb_def2.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=3.023127, coefficient=26.224898))
Rb_def2.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=0.650383, coefficient=0.962839))
Rb_def2.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=3.843114, coefficient=12.3169))
Rb_def2.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=3, exponent=3.843114, coefficient=-12.3169))

lanl_ECP = ECP(28)
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=0, angular_momentum=0, exponent=54.1980682, coefficient=3.0))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=1, angular_momentum=0, exponent=32.9053558, coefficient=27.3430642))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=13.6744890, coefficient=118.8028847))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=3.0341152, coefficient=43.4354876))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=0, angular_momentum=1, exponent=54.256334, coefficient=5.0))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=1, angular_momentum=1, exponent=26.0095593, coefficient=25.0504252))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=28.2012995, coefficient=92.6157463))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=9.4341061, coefficient=95.8249016))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=2.5321764, coefficient=26.2684983))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=0, angular_momentum=2, exponent=87.6328721, coefficient=3.0))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=1, angular_momentum=2, exponent=61.7373377, coefficient=22.5533557))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=32.4385104, coefficient=178.1241988))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=8.7537199, coefficient=76.9924162))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=1.6633189, coefficient=9.481827))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=1, angular_momentum=3, exponent=213.6143969, coefficient=-28.0))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=3, exponent=41.058538, coefficient=-134.9268852))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=3, exponent=8.708653, coefficient=-41.9271913))
lanl_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=3, exponent=2.6074661, coefficient=-5.933642))

e_br = [1.42, 3.56, 6.35, 10.4, 20.7, 32.4, 76.5, 142.1, 341.871, 1000.205, 220.0000]
d_br0 = [4.202568, 13.194570, 10.795897, 26.648345, 33.130248, 45.359127, 35.346799, -31.766999, -633.120024, -47285358, 3.0]
d_br1 = [3.611982, 10.720320, 5.828031, 34.390992, 61.018209, 115.349011, 257.006175, 161.631632, -144.676033, -9413723, 5.0]
d_br2 = [3.484330, -0.211214, 9.199904, -3.542208, 43.713385, 77.076843, 245.662800, 509.898585, -330.399360, -21536241, 7.0]
n_br = [2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0]
KG_ECP = ECP(35)
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=1.42, coefficient=4.202568))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=3.56, coefficient=13.194570))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=6.35, coefficient=10.795897))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=10.4, coefficient=26.648345))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=20.7, coefficient=33.130248))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=32.4, coefficient=45.359127))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=76.5, coefficient=35.346799))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=142.1, coefficient=-31.766999))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=0, exponent=341.871, coefficient=-633.120024))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=1, angular_momentum=0, exponent=1000.205, coefficient=-47285358))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=0, angular_momentum=0, exponent=220.0000, coefficient=3.0))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=1.42, coefficient=3.611982))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=3.56, coefficient=10.720320))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=6.35, coefficient=5.828031))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=10.4, coefficient=34.390992))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=20.7, coefficient=61.018209))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=32.4, coefficient=115.349011))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=76.5, coefficient=257.006175))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=142.1, coefficient=161.631632))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=1, exponent=341.871, coefficient=-144.676033))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=1, angular_momentum=1, exponent=1000.205, coefficient=-9413723))

KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=1.42, coefficient=3.484330))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=3.56, coefficient=-0.211214))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=6.35, coefficient=9.199904))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=10.4, coefficient=-3.542208))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=20.7, coefficient=43.713385))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=32.4, coefficient=77.076843))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=76.5, coefficient=245.662800))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=142.1, coefficient=509.898585))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=2, exponent=341.871, coefficient=-330.399360))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=1, angular_momentum=2, exponent=1000.205, coefficient=-21536241))
KG_ECP.add_ECP_Gaussian(ECP_Gaussian(radial_exp=2, angular_momentum=3, exponent=1.42, coefficient=213.6143969))

l_facs = [0.28209479**2, 0.48860251**2, 0.63078313**2, 0.74635267**2]


def ECP_func(r, Z=28, _ECP=Rb_def2):
    res = 0
    if Z < 0:
        raise ValueError("Z is negative")

    res = _ECP.evaluate(r, _ECP.L)

    for l in range(_ECP.L):
        res += _ECP.evaluate(r, l) * l_facs[l]

    return res / fPI


def get_numbers_from_filename(filename):
    # Use a regular expression to find all groups of digits in the filename
    numbers = re.findall(r"\d+", filename)
    str = ""
    str = str.join(numbers)

    # Convert the numbers to integers and return them
    return str


def Slm(l, m, theta, phi):
    if m == 0:
        return np.real(sp.sph_harm(m, l, phi, theta))
    elif m > 0:
        if l == 1:
            return np.real(sp.sph_harm(-m, l, phi, theta) - sp.sph_harm(m, l, phi, theta)) / np.sqrt(2)
        else:
            if m % 2 == 1:
                return np.real(sp.sph_harm(-m, l, phi, theta) - sp.sph_harm(m, l, phi, theta)) / np.sqrt(2)
            else:
                return np.real(sp.sph_harm(-m, l, phi, theta) + sp.sph_harm(m, l, phi, theta)) / np.sqrt(2)
    elif m < 0:
        if l == 1:
            return np.imag(sp.sph_harm(-m, l, phi, theta) + sp.sph_harm(m, l, phi, theta)) / np.sqrt(2)
        else:
            if m % 2 == 1:
                return np.imag(sp.sph_harm(-m, l, phi, theta) + sp.sph_harm(m, l, phi, theta)) / np.sqrt(2)
            else:
                return np.imag(sp.sph_harm(-m, l, phi, theta) - sp.sph_harm(m, l, phi, theta)) / np.sqrt(2)


def Slm_self(l, m, theta, phi):
    if l == 0:
        if m != 0:
            return 0.0
        else:
            return np.sqrt(1 / 4 / np.pi) * sp.lpmv(m, l, np.cos(theta))
    elif l == 1:
        if m == 0:
            return np.sqrt(3 / 4 / np.pi) * sp.lpmv(m, l, np.cos(theta))
        elif m == 1:
            return -np.sqrt(3 / 4 / np.pi) * np.cos(phi) * sp.lpmv(m, l, np.cos(theta))
        elif m == -1:
            return -np.sqrt(3 / 4 / np.pi) * np.sin(phi) * sp.lpmv(m, l, np.cos(theta)) * 2
        else:
            return 0.0
    elif l == 2:
        if m == 0:
            return np.sqrt(5 / 4 / np.pi) * sp.lpmv(m, l, np.cos(theta))
        elif m == 1:
            return -np.sqrt(5 / 12 / np.pi) * np.cos(phi) * sp.lpmv(m, l, np.cos(theta))
        elif m == -1:
            return -np.sqrt(5 / 12 / np.pi) * np.sin(phi) * sp.lpmv(m, l, np.cos(theta))
        elif m == 2:
            return np.sqrt(5 / 48 / np.pi) * np.cos(2 * phi) * sp.lpmv(m, l, np.cos(theta))
        elif m == -2:
            return np.sqrt(5 / 48 / np.pi) * np.sin(2 * phi) * sp.lpmv(m, l, np.cos(theta))
        else:
            return 0.0


def check_SLM():
    for l in range(4):
        for m in range(-l, l + 1, 1):
            print("l: ", l, "m: ", m)
            # Check if Slm and Slm_self are the same
            for theta in np.linspace(0, np.pi, 1):
                for phi in np.linspace(0, 2 * np.pi, 1):
                    print("theta: ", theta, "phi: ", phi)
                    print(f"Scipy:     {Slm(l,m,theta,phi):14.8f}")
                    # print(f"selfbuild: {Slm_self(l,m,theta,phi):14.8f}")


def br_test():
    y = np.zeros_like(x)
    z = np.zeros_like(x)
    a = np.zeros_like(x)
    y_lanl = np.zeros_like(x)
    z_lanl = np.zeros_like(x)
    a_lanl = np.zeros_like(x)
    c_lanl = np.zeros_like(x)
    lanl_density = np.zeros_like(x)
    KG_dens = np.zeros_like(x)

    for i in range(11):
        y += 0.28209479**2 * d_br0[i] * np.exp(-e_br[i] * x**2) * x ** (n_br[i] - 2)
        z += 0.48860251**2 * d_br1[i] * np.exp(-e_br[i] * x**2) * x ** (n_br[i] - 2)
        a += 0.63078313**2 * d_br2[i] * np.exp(-e_br[i] * x**2) * x ** (n_br[i] - 2)
        if n_br[i] == 2:
            KG_dens += (
                0.28209479**2 * 4 * d_br0[i] * e_br[i] * np.exp(-e_br[i] * x**2) * (-1 + e_br[i] * x**2)
                + 0.48860251**2 * 4 * d_br1[i] * e_br[i] * np.exp(-e_br[i] * x**2) * (-1 + e_br[i] * x**2)
                + 0.63078313**2 * 4 * d_br2[i] * e_br[i] * np.exp(-e_br[i] * x**2) * (-1 + e_br[i] * x**2)
            )
        elif n_br[i] == 1:
            KG_dens += (
                0.28209479**2 * d_br0[i] * np.exp(-e_br[i] * x**2) * (1 + 4 * e_br[i] ** 2 * x**4) / x**3
                + 0.48860251**2 * d_br1[i] * np.exp(-e_br[i] * x**2) * (1 + 4 * e_br[i] ** 2 * x**4) / x**3
                + 0.63078313**2 * d_br2[i] * np.exp(-e_br[i] * x**2) * (1 + 4 * e_br[i] ** 2 * x**4) / x**3
            )
        elif n_br[i] == 0:  #  used to be 2
            KG_dens += (
                0.28209479**2 * 4 * d_br2[i] * np.exp(-e_br[i] * x**2) * (1 + e_br[i] * x**2 + e_br[i] ** 2 * x**4) / x**4
                + 0.48860251**2 * 4 * d_br1[i] * np.exp(-e_br[i] * x**2) * (1 + e_br[i] * x**2 + e_br[i] ** 2 * x**4) / x**4
                + 0.63078313**2 * 4 * d_br0[i] * np.exp(-e_br[i] * x**2) * (1 + e_br[i] * x**2 + e_br[i] ** 2 * x**4) / x**4
            )

    for i in range(4):
        y_lanl += 0.28209479**2 * lanl_d1[i] * np.exp(-lanl_e1[i] * x**2) * x ** (lanl_n1[i])
        c_lanl += 0.74635267**2 * lanl_d4[i] * np.exp(-lanl_e4[i] * x**2) * x ** (lanl_n4[i])

        if lanl_n1[i] == 0:
            lanl_density -= 0.28209479**2 * 4 * lanl_d1[i] * lanl_e1[i] * np.exp(-lanl_e1[i] * x**2) * (-1 + lanl_e1[i] * x**2)
        elif lanl_n1[i] == 1:
            lanl_density -= 0.28209479**2 * lanl_d1[i] * np.exp(-lanl_e1[i] * x**2) * (1 - 8 * lanl_e1[i] * x**2 + 4 * lanl_e1[i] ** 2 * x**4) / x
        elif lanl_n1[i] == 2:
            lanl_density -= 0.28209479**2 * 4 * lanl_d1[i] * np.exp(-lanl_e1[i] * x**2) * (1 - 3 * lanl_e1[i] * x**2 + lanl_e1[i] ** 2 * x**4)

        if lanl_n4[i] == 0:
            lanl_density -= 0.74635267**2 * 4 * lanl_d4[i] * lanl_e4[i] * np.exp(-lanl_e4[i] * x**2) * (-1 + lanl_e4[i] * x**2)
        elif lanl_n4[i] == 1:
            lanl_density -= 0.74635267**2 * lanl_d4[i] * np.exp(-lanl_e4[i] * x**2) * (1 - 8 * lanl_e4[i] * x**2 + 4 * lanl_e4[i] ** 2 * x**4) / x
        elif lanl_n4[i] == 2:
            lanl_density -= 0.74635267**2 * 4 * lanl_d4[i] * np.exp(-lanl_e4[i] * x**2) * (1 - 3 * lanl_e4[i] * x**2 + lanl_e4[i] ** 2 * x**4)
    for i in range(5):
        z_lanl += 0.48860251**2 * lanl_d2[i] * np.exp(-lanl_e2[i] * x**2) * x ** (lanl_n2[i])
        a_lanl += 0.63078313**2 * lanl_d3[i] * np.exp(-lanl_e3[i] * x**2) * x ** (lanl_n3[i])

        if lanl_n2[i] == 0:
            lanl_density -= 0.48860251**2 * 4 * lanl_d2[i] * lanl_e2[i] * np.exp(-lanl_e2[i] * x**2) * (-1 + lanl_e2[i] * x**2)
        elif lanl_n2[i] == 1:
            lanl_density -= 0.48860251**2 * lanl_d2[i] * np.exp(-lanl_e2[i] * x**2) * (1 - 8 * lanl_e2[i] * x**2 + 4 * lanl_e2[i] ** 2 * x**4) / x
        elif lanl_n2[i] == 2:
            lanl_density -= 0.48860251**2 * 4 * lanl_d2[i] * np.exp(-lanl_e2[i] * x**2) * (1 - 3 * lanl_e2[i] * x**2 + lanl_e2[i] ** 2 * x**4)

        if lanl_n3[i] == 0:
            lanl_density -= 0.63078313**2 * 4 * lanl_d3[i] * lanl_e3[i] * np.exp(-lanl_e3[i] * x**2) * (-1 + lanl_e3[i] * x**2)
        elif lanl_n3[i] == 1:
            lanl_density -= 0.63078313**2 * lanl_d3[i] * np.exp(-lanl_e3[i] * x**2) * (1 - 8 * lanl_e3[i] * x**2 + 4 * lanl_e3[i] ** 2 * x**4) / x
        elif lanl_n3[i] == 2:
            lanl_density -= 0.63078313**2 * 4 * lanl_d3[i] * np.exp(-lanl_e3[i] * x**2) * (1 - 3 * lanl_e3[i] * x**2 + lanl_e3[i] ** 2 * x**4)

    b = a + z + y + 28 / x
    y += 28 / x
    z += 28 / x
    a += 28 / x
    KG_dens += 28 / x**3

    KG_radial_dens = KG_dens * x**2 * fPI
    lanl_radial_dens = lanl_density * x**2 * fPI

    fig, axs = plt.subplots(4, 1, figsize=(14, 8))

    axs[0].plot(x, y_lanl, linestyle="-", label="l=0")
    axs[0].plot(x, z_lanl, linestyle="-", label="l=1")
    axs[0].plot(x, a_lanl, linestyle="-", label="l=2")
    axs[0].plot(x, y_lanl + z_lanl + a_lanl - c_lanl, linestyle="-", label="sum")
    axs[0].plot(x, -c_lanl, linestyle="-", label="(-) l=3")

    axs[0].set_xlabel("r")
    axs[0].set_ylabel(r"$U_l(r)$ LANL")
    axs[0].legend()
    # axs[0].set_ylim(0,15.0)
    axs[0].set_xlim(0, 2.5)

    axs[1].plot(x, y, linestyle="-", label="l=0")
    axs[1].plot(x, z, linestyle="-", label="l=1")
    axs[1].plot(x, a, linestyle="-", label="l=2")
    axs[1].plot(x, b, linestyle="-", label="sum")

    axs[1].set_xlabel("r")
    axs[1].set_ylabel(r"$(U_l(r) - \frac{N_C}{r})$ KG")
    axs[1].legend()
    axs[1].set_ylim(0.0, 150.0)
    axs[1].set_xlim(0, 2.5)

    axs[2].plot(x, KG_dens, linestyle="-", label="KG")
    axs[2].plot(x, lanl_density, linestyle="-", label="LANL")

    axs[2].set_ylim(-50, 300.0)
    axs[2].set_xlim(0, 2.5)
    axs[2].legend()

    axs[3].plot(x, KG_radial_dens, linestyle="-", label="KG")
    axs[3].plot(x, lanl_radial_dens, linestyle="-", label="LANL")

    axs[3].set_ylim(-20, 100.0)
    axs[3].set_xlim(0, 2.5)
    axs[3].legend()

    fig.savefig("Br_test.png", dpi=300, bbox_inches="tight")
    exit(0)


# check_SLM()
if __name__ == "__main__":
    dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(dir)
    lf = os.path.join(dir, "core_dens.dat")
    # Read the file
    # 0 = k
    # 1 = Thakakr FF
    # 2 = Gaussian FF
    # 3 = Thakkar Core Density
    # 4 = Gaussian Core Density
    # 5 = Thakkar Density
    # 6 = ORCA Density
    # 7 = ORCA Core Density
    # 8 = r
    # 9 = ECP Density
    # 10 = QZVP Density
    # 11 = QZVP Core Density
    df = pd.read_csv(lf, delim_whitespace=True, header=None)

    upper = 13.0
    nr = 800000

    x = np.linspace(0.000001, upper, nr)
    dx = x[1]
    y = np.zeros_like(x)
    z = np.zeros_like(x)

    # br_test()

    def calculate_ecp_func(x):
        return ECP_func(x, 28, Rb_def2)

    import multiprocessing
    from itertools import repeat

    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        res = pool.starmap(ECP_func, zip(x))  # x ist der Vector für den was gemacht werden soll
        for i, r in enumerate(res):
            y[i] = r

    print("We have y")
    z = y * x**2 * fPI

    integral = sum(x**2 * y * (x[1] - x[0]) * fPI)
    print("ECP Integral: ", integral)

    integral2 = sum(df[8][1] * df[5] * (df[8]) ** 2 * fPI)
    print("Thakkar Integral:", integral2)
    integral2 = sum(df[8][1] * df[3] * (df[8]) ** 2 * fPI)
    print("Thakkar Core Integral:", integral2)
    integral2 = sum(df[8][1] * df[6] * (df[8]) ** 2 * fPI)
    print("ORCA Integral:", integral2)
    integral2 = sum(df[8][1] * df[7] * (df[8]) ** 2 * fPI)
    print("ORCA Core Integral:", integral2)

    integral2 = sum(df[8][1] * df[9] * (df[8]) ** 2 * fPI)
    print("ORCA ECP Integral:", integral2)

    integral2 = sum(df[8][1] * (df[6] - df[7]) * (df[8]) ** 2 * fPI)
    print("ORCA - ORCA Core Integral:", integral2)

    integral2 = sum(df[8][1] * df[10] * (df[8]) ** 2 * fPI)
    print("ORCA QZVP Integral:", integral2)
    integral2 = sum(df[8][1] * df[11] * (df[8]) ** 2 * fPI)
    print("ORCA QZVP Core Integral:", integral2)

    integral2 = sum(df[8][1] * (df[10] - df[11]) * (df[8]) ** 2 * fPI)
    print("ORCA - ORCA Core Integral:", integral2)

    # Create a figure with 3x2 subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    # axs[0][0].set_title('ECP SF')
    # axs[0][0].plot(df[0], df[1], label='Thakkar')
    # axs[0][0].plot(df[0], df[2], label='Gaussian')

    axs[0][0].set_title("Radial dens diff")
    temp = df[9] + df[3]  # ECP + Thakkar Core
    axs[0][0].plot(
        df[8],
        (df[8]) ** 2 * (temp - df[5]) * fPI,
        label="ECP + Thakkar Core - Thakkar full",
    )
    axs[0][0].plot(
        df[8],
        (df[8]) ** 2 * (temp - df[6]) * fPI,
        label="ECP + Thakkar Core - ORCA full",
    )
    axs[0][0].plot(df[8], (df[8]) ** 2 * (df[5] - df[6]) * fPI, label="Thakkar - ORCA full")
    axs[0][0].plot(df[8], (df[8]) ** 2 * (df[5] - df[10]) * fPI, label="Thakkar - ORCA QZVP")
    axs[0][0].set_xlim(-0.1, 3)
    axs[0][0].legend()

    axs[0][1].set_title("ECP Density")
    axs[0][1].plot(x, y, label="Gaussian python")
    axs[0][1].plot(df[8], df[3], label="Thakkar")
    # axs[0][1].plot(df[0]*dx_log, df[4], label='Gaussian')
    axs[0][1].plot(df[8], df[5], label="Thakkar full")
    axs[0][1].plot(df[8], df[6], label="ORCA")
    axs[0][1].plot(df[8], df[9], label="ECP")
    axs[0][1].plot(df[8], df[9] + df[3], label="ECP + Thakkar Core")
    axs[0][1].set_ylim(-5, 25)
    axs[0][1].set_xlim(-0.1, 3)
    axs[0][1].legend()

    axs[1][0].set_title("Valence Density")
    # axs[1][0].plot(x, y, label='Gaussian python')
    axs[1][0].plot(df[8], df[5] - df[3], label="Thakkar Valence")
    # axs[1][0].plot(df[0]*dx_log, df[5] - df[4], label='Gaussian Valence')
    axs[1][0].plot(df[8], df[6] - df[7], label="ORCA Valence")
    axs[1][0].plot(df[8], df[9], label="ORCA ECP Valence")
    axs[1][0].plot(df[8], df[6] - df[5], label="Thakkar vs ORCA")
    axs[1][0].plot(df[8], df[10] - df[11], label="QZVP Valence")

    axs[1][0].set_ylim(-2, 15)
    axs[1][0].set_xlim(-0.1, 2)
    axs[1][0].legend()

    axs[1][1].set_title("ECP radial Density")
    axs[1][1].plot(x, z, label="Gaussian python")
    axs[1][1].plot(df[8], (df[8]) ** 2 * df[3] * fPI, label="Thakkar")
    # axs[1][1].plot(df[0]*dx_log, (df[0]*dx_log)**2 * df[4] * fPI, label='Gaussian')
    axs[1][1].plot(df[8], (df[8]) ** 2 * df[5] * fPI, label="Thakkar full")
    axs[1][1].plot(df[8], (df[8]) ** 2 * df[6] * fPI, label="ORCA")
    axs[1][1].plot(df[8], (df[8]) ** 2 * df[9] * fPI, label="ECP")
    axs[1][1].plot(df[8], (df[8]) ** 2 * (df[9] + df[3]) * fPI, label="ECP +  Thakkar Core")
    axs[1][1].set_xlim(-0.1, 3)
    axs[1][1].set_ylim(-3, 60)
    axs[1][1].legend()

    # Show the figure
    # plt.show()
    # exit(0)
    fig.savefig("ECP_SF.png", dpi=300, bbox_inches="tight")
    exit(0)
