import numpy as np
from scipy.integrate import quad


# fit a function that sums to an integral of 0 over the whole range, that maps the two functions onto each other
def func(_x, params):
    return np.array([a * np.clip(np.exp(-b * _x**2), 0, 1) for a, b in zip(params[::2], params[1::2])]).sum(axis=0)


def int_func(_x, params):
    return func(_x, params) * _x**2 * 4 * np.pi


def diff(params, data, x):
    x2 = np.roll(x, 1)
    x2[0] = 0
    return (func(x, params) - data) * 4 * np.pi * x**-0.5


def diff1(params, data, x):
    x2 = np.roll(x, 1)
    x2[0] = 0
    return (func(x, params) - data) * 4 * np.pi * x**2


# integrate the differences
def diff3(params, data, x):
    x2 = np.roll(x, 1)
    x2[0] = 0
    return func(x, params) - data


def diff2(params, data, x):
    x2 = np.roll(x, 1)
    x2[0] = 0
    return (func(x, params) - data) * 4 * np.pi * x


# The callback function
def callback(xk, axs=None, _x=None, steps=0, data=None, element="H"):
    axs[0].plot(_x, func(_x, xk), label=f"fit: {to_minimize(xk, data, _x):.6e}")
    axs[1].plot(_x, int_func(_x, xk), label=f"fit: {to_minimize(xk, data, _x):.6e}")
    axs[2].plot(_x, func(_x, xk) * 4 * np.pi * _x, label=f"fit: {to_minimize(xk, data, _x):.6e}")
    plt.title(f"Step {steps}")
    axs[0].legend(ncol=3)
    axs[1].legend(ncol=3)
    axs[2].legend(ncol=3)
    # if there is more than 12 plots remove the first one
    if len(axs[0].lines) >= 12:
        axs[0].lines[2].remove()
        axs[1].lines[2].remove()
        axs[2].lines[2].remove()
        # rescale y-axis using automatic scaling
        axs[0].relim()
        axs[0].autoscale_view()
        axs[1].relim()
        axs[1].autoscale_view()
        axs[2].relim()
        axs[2].autoscale_view()
    plt.savefig(f"{element}_fit.png")


# The new function to minimize
def to_minimize(params, data, x):
    integral = quad(int_func, 0, np.inf, args=params, limit=1000, epsabs=1e-5, epsrel=1e-3)[0] ** 2 * 1e6
    # *x**0
    diff2_res = sum((diff3(params, data, x)) ** 2) * 1e-3
    # *x
    diff1_res = sum((diff2(params, data, x)) ** 2) * 1e-1
    # *x**2
    diff_res = sum((diff1(params, data[1000:], x[1000:])) ** 2) * 1e0
    # *x**-1
    diff_1_res = sum((diff(params, data[:5000], x[:5000])) ** 2) * 1e-9
    res = integral + diff_res + diff2_res + diff1_res + diff_1_res
    return res


class CallbackWrapper:
    def __init__(self, func, callback, figs, axs, data, x_outer, element="H"):
        self.func = func
        self.callback = callback
        self.best_result = None
        self.last_result = None
        self.figs = figs
        self.axs = axs
        self.data = data
        self.x_outer = x_outer
        self.steps = 0
        self.element = element

    def __call__(self, x, convergence=None, return_values=False):
        result = self.func(x, self.data, self.x_outer)
        # check for last result in the file and compare it to the new result
        # if self.last_result is None:
        if os.path.exists(f"{self.element}_best_result.txt"):
            with open(f"{self.element}_best_result.txt") as f:
                lines = f.readline()
                self.last_result = float(lines.split(":")[-1])
        else:
            self.last_result = np.inf
        if result < self.last_result:
            if self.best_result is None or result < self.best_result:
                import pandas as pd

                if convergence is not None:
                    print(f"Convergence: {convergence}")
                self.best_result = result
                self.steps += 1
                # sort results by b
                df = pd.DataFrame({"a": x[::2], "b": x[1::2]}).sort_values(by="b")
                # overwrite a file with the best result
                with open(f"{self.element}_best_result.txt", "w") as f:
                    f.write(f"best result: {result:.6e}\n")
                    f.write("a_vec =[                        b_vec = [\n")
                    # write df values separated by ,
                    for i, row in df.iterrows():
                        f.write(f"{row['a']:28.16f}, {row['b']:28.16f},\n")
                    f.write("       ]                                ]\n")
                self.callback(x, axs=self.axs, _x=self.x_outer, steps=self.steps, data=self.data, element=self.element)


if __name__ == "__main__":
    # include argparser for the Element symbol
    import argparse
    from scipy.optimize import minimize, basinhopping, differential_evolution
    import os

    # plotting
    import matplotlib.pyplot as plt
    import pandas as pd

    # argparse for the element symbol
    parser = argparse.ArgumentParser(description="Fit the difference between two atomic densities")
    parser.add_argument("element", type=str, help="The element symbol")
    args = parser.parse_args()

    print("Working on Element: ", args.element)

    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Rb atomic densities as pandas frame from the dat files and make them to numpy arrays
    def2_dens = pd.read_csv(f"{args.element}_def2.dat", sep="\s+", header=0).to_numpy()
    jorge_val = pd.read_csv(f"{args.element}_jorge_val.dat", sep="\s+", header=0).to_numpy()

    x = def2_dens[5000:, 0]
    def2 = def2_dens[5000:, 1]
    jorge = jorge_val[5000:, 1]
    data = jorge - def2

    # read a_vec and b_vec from the file if it exists, otherwise use template
    if not os.path.exists(f"{args.element}_best_result.txt"):
        a_vec = [
            -0.0009059044761585,
            0.0049165720445790,
            0.0026112054093839,
            0.5129121276932708,
            -16.5546029905644438,
            -0.0559565805267245,
            69.3065597961078481,
            -150.0372414878611664,
            224.5900807367371499,
            -64.1524732997385456,
            236.7772067557729940,
            235.1495391730983613,
            356.4627129901924150,
            521.1718946738141085,
            789.1739921585068487,
        ]
        b_vec = [
            0.1143892643492925,
            0.3677731762352405,
            0.6040692473928437,
            1.1197799383593285,
            6.4771516377879959,
            7.9348712618780306,
            21.3568053863003442,
            76.1922571598153837,
            173.7562999945938884,
            510.4032544955319395,
            6164.4564757970965729,
            26976.5326768802333390,
            107120.6178712070250185,
            719310.4340740576153621,
            5932074.4300301708281040,
        ]
    else:
        with open(f"{args.element}_best_result.txt") as f:
            lines = f.readlines()
            a_vec = []
            b_vec = []
            for line in lines[2:-1]:
                a, b, junk = line.split(",")
                a_vec.append(float(a))
                b_vec.append(float(b))
    guess = []
    for a, b in zip(a_vec, b_vec):
        guess += [a, b]
    # set the bounds for the parameters
    bounds = []
    for i in range(len(a_vec)):
        if a_vec[i] < 0:
            bounds += [
                (2 * a_vec[i], 0.5 * a_vec[i]),
                (b_vec[i] * 0.5, b_vec[i] * 2),
            ]
        else:
            bounds += [
                (0.5 * a_vec[i], 2 * a_vec[i]),
                (b_vec[i] * 0.5, b_vec[i] * 2),
            ]
    print("bounds: ", bounds)
    # minimize the error function
    integral = quad(int_func, 0, np.inf, args=guess, limit=1000, epsabs=1e-5, epsrel=1e-3)[0] ** 2
    diff_res = sum((diff1(guess, data[1000:], x[1000:])) ** 2)
    diff2_res = sum((diff3(guess, data, x)) ** 2)
    diff1_res = sum((diff2(guess, data, x)) ** 2)
    diff_1_res = sum((diff(guess, data[:5000], x[:5000])) ** 2)
    figs, axs = plt.subplots(3)
    callback_wrapper = CallbackWrapper(to_minimize, callback, figs, axs, data, x, element=args.element)
    callback_wrapper.best_result = to_minimize(guess, data, x)
    print(
        f"components: {integral:.6e} + {diff2_res:.6e} + {diff1_res:.6e} + {diff_res:.6e} + {diff_1_res:.6e}, total: {callback_wrapper.best_result:.6e}"
    )
    fig_width, fig_height = 12, 7  # these are in inches
    figs.set_size_inches(fig_width, fig_height)
    plt.tight_layout()
    axs[0].plot(x, data, label="difference")
    axs[1].plot(x, data * 4 * np.pi * x**2, label="difference")
    axs[2].plot(x, data * 4 * np.pi * x, label="difference")
    axs[0].set_xscale("log")
    axs[1].set_xscale("log")
    axs[2].set_xscale("log")
    axs[0].set_xlim(1e-4, x[-1])
    axs[1].set_xlim(1e-4, x[-1])
    axs[2].set_xlim(1e-4, x[-1])
    # set figure legend fontsize to tiny
    plt.rcParams.update({"legend.fontsize": "x-small"})
    axs[0].plot(x, func(x, guess), label=f"last best: {callback_wrapper.best_result:.4e}")
    axs[1].plot(x, int_func(x, guess), label=f"last best: {callback_wrapper.best_result:.4e}")
    axs[2].plot(
        x,
        func(x, guess) * 4 * np.pi * x,
        label=f"last best: {callback_wrapper.best_result:.4e}",
    )
    plt.title("Start guess")
    # make the legend in 3 columns
    axs[0].legend(ncol=3)
    axs[1].legend(ncol=3)
    axs[2].legend(ncol=3)
    plt.savefig(f"{args.element}_fit.png")
    # minimizer_kwargs = {
    #    "method": "L-BFGS-B",
    #    "args": (data, x),
    #    "bounds": bounds,
    #    "jac": "2-point",
    #    "callback": callback_wrapper,
    #    "options": {
    #        "ftol": 1e-10,
    #        "gtol": 1e-7,
    #        "maxfun": 50000,
    #        "disp": True,
    #        "maxiter": 1000,
    #        "maxls": 100,
    #    },
    # }
    # popt = basinhopping(
    #    to_minimize,
    #    guess,
    #    minimizer_kwargs=minimizer_kwargs,
    #    callback=callback_wrapper,
    #    niter=200,
    #    T=10,
    # ).x
    popt = minimize(
        to_minimize,
        guess,
        args=(data, x),
        method="L-BFGS-B",
        bounds=bounds,
        jac="3-point",
        callback=callback_wrapper,
        options={
            "ftol": 1e-5,
            "gtol": 1e-5,
            "maxfun": 200000,
            "disp": True,
            "maxiter": 1000,
            "maxls": 100,
        },
    ).x
    guess = popt
    a_vec = guess[::2]
    b_vec = guess[1::2]
    # set the bounds for the parameters
    bounds = []
    for i in range(len(a_vec)):
        if a_vec[i] < 0:
            bounds += [
                (2 * a_vec[i], 0.5 * a_vec[i]),
                (b_vec[i] * 0.5, b_vec[i] * 2),
            ]
        else:
            bounds += [
                (0.5 * a_vec[i], 2 * a_vec[i]),
                (b_vec[i] * 0.5, b_vec[i] * 2),
            ]
    popt = differential_evolution(
        to_minimize,
        bounds,
        x0=guess,
        args=(data, x),
        callback=callback_wrapper,
        strategy="best1bin",
        workers=10,
        updating="deferred",
        maxiter=100000,
        popsize=15,
        disp=True,
        tol=1e-10,
        polish=True,
    ).x
    # print out the optimized parameters pairwise:
    print("a_vec =[                        b_vec = [")
    for i in range(len(a_vec)):
        print(f"{popt[2 * i]:28.16f}, {popt[2 * i + 1]:28.16f},")
    print("       ]                                ]")
    integral = quad(int_func, 0, np.inf, args=popt, limit=1000, epsabs=1e-5, epsrel=1e-3)[0] ** 2
    diff_res = sum((diff1(popt, data[2000:], x[2000:])) ** 2)
    diff2_res = sum((diff3(popt, data, x)) ** 2)
    diff1_res = sum((diff2(popt, data, x)) ** 2)
    diff_1_res = sum((diff(popt, data[:5000], x[:5000])) ** 2)
    print(
        f"components: {integral:.6e} + {diff_res:.6e} + {diff2_res:.6e} + {diff1_res:.6e} + {diff_1_res:.6e}, total: {callback_wrapper.best_result:.6e}"
    )
