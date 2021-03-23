"""
Simulation causal time series
Date: Oct 2020
Author: Karim Assaad, karimassaad3@gmail.com
"""

import pandas as pd
import numpy as np
import random


def abso(x):
    return np.sqrt(np.power(x,2))


def identity(x):
    return x


functions_set = {0: abso, 1: np.tanh, 2:np.sin, 3: np.cos}


def uniform_with_gap(min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    while True:
        r = random.uniform(min_value, max_value)
        if min_gap>r or max_gap<r:
            break
    return r


def fork_generator(N=1000, non_linear=True, min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    summary = np.diag(np.diag(np.ones([3,3])))
    summary[0,1] = 2
    summary[1,0] = 1
    summary[0,2] = 2
    summary[2,0] = 1
    temporal = dict()
    temporal[0] = [(0, -1)]
    temporal[1] = [(1, -1), (0, -1)]
    temporal[2] = [(2, -1), (0, -2)]

    N=N+2
    ax = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bx = 0.1
    ay = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    by = 0.1
    aw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bw = 0.1

    axy = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    axw = uniform_with_gap(min_value,max_value, min_gap, max_gap)

    if non_linear:
        f = functions_set[random.randint(0, len(functions_set)-1)]
        g = functions_set[random.randint(0, len(functions_set)-1)]
    else:
        f = identity
        g = identity

    print("Fork: 1 <- 0 -> 2")
    print(f, g, ax, bx, ay,by, aw, bw, axy, axw)
    epsx = np.random.randn(N) ** 3
    epsy = np.random.randn(N) ** 3
    epsw = np.random.randn(N) ** 3

    x = np.zeros([N])
    y = np.zeros([N])
    w = np.zeros([N])

    x[0] = bx * epsx[0]
    x[1] = ax * x[0] + bx * epsx[1]
    y[1] = ay * y[0] + axy * x[0] + by * epsy[1]
    x[2] = ax * x[1] + bx * epsx[2]
    y[2] = ay * y[1] + axy * x[1] + by * epsy[2]
    w[2] = aw * w[1] + axw * x[0] + bw * epsw[2]
    for i in range(3, N):
        x[i] = ax * x[i - 1] + bx * epsx[i]
        y[i] = ay * y[i - 1] + axy * f(x[i - 1]) + by * epsy[i]
        w[i] = aw * w[i - 1] + axw * g(x[i - 2]) + bw * epsw[i]

    x = pd.DataFrame(x, columns=["V1"])
    y = pd.DataFrame(y, columns=["V2"])
    w = pd.DataFrame(w, columns=["V3"])

    series = pd.concat([x, y, w], axis=1, sort=False)
    series = series.drop(series.index[[0,1]])
    series = series.reset_index(drop=True)
    series.index.names = ['time_index']
    return series, summary, temporal


def v_structure_generator(N=1000, non_linear=True, min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    summary = np.diag(np.diag(np.ones([3,3])))
    summary[0,2] = 2
    summary[2,0] = 1
    summary[1, 2] = 2
    summary[2,1] = 1
    temporal = dict()
    temporal[0] = [(0, -1)]
    temporal[1] = [(1, -1)]
    temporal[2] = [(2, -1), (0, -2), (1, -1)]

    N=N+2
    ax = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bx = 0.1
    ay = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    by = 0.1
    aw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bw = 0.1

    axw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ayw = uniform_with_gap(min_value,max_value, min_gap, max_gap)

    if non_linear:
        f = functions_set[random.randint(0, len(functions_set)-1)]
        g = functions_set[random.randint(0, len(functions_set)-1)]
    else:
        f = identity
        g = identity

    print("V-structure: 0 -> 2 <- 1")
    print(f, g, ax, bx, ay,by, aw, bw, axw, ayw)
    epsx = np.random.randn(N) ** 3
    epsy = np.random.randn(N) ** 3
    epsw = np.random.randn(N) ** 3

    x = np.zeros([N])
    y = np.zeros([N])
    w = np.zeros([N])

    x[0] = bx * epsx[0]
    y[0] = by * epsy[0]
    x[1] = ax * x[0] + bx * epsx[1]
    y[1] = ay * y[0] + by * epsy[1]
    x[2] = ax * x[1] + bx * epsx[2]
    y[2] = ay * y[1] + by * epsy[2]
    w[2] = aw * w[1] + axw * x[0] + ayw * y[1] + bw * epsw[2]
    for i in range(3, N):
        x[i] = ax * x[i - 1] + bx * epsx[i]
        y[i] = ay * y[i - 1] + by * epsy[i]
        w[i] = aw * w[i - 1] + axw * f(x[i - 2]) + ayw * g(y[i - 1]) + bw * epsw[i]

    x = pd.DataFrame(x, columns=["V1"])
    y = pd.DataFrame(y, columns=["V2"])
    w = pd.DataFrame(w, columns=["V3"])

    series = pd.concat([x, y, w], axis=1, sort=False)
    series = series.drop(series.index[[0,1]])
    series = series.reset_index(drop=True)
    series.index.names = ['time_index']
    return series, summary, temporal


def mediator_generator(N=1000, non_linear=True, min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    summary = np.diag(np.diag(np.ones([3, 3])))
    summary[0, 1] = 2
    summary[1, 0] = 1
    summary[0, 2] = 2
    summary[2, 0] = 1
    summary[1, 2] = 2
    summary[2, 1] = 1

    temporal = dict()
    temporal[0] = [(0, -1), (0, -1)]
    temporal[1] = [(1, -1), (0, -1)]
    temporal[2] = [(2, -1), (0, -2), (1, -1)]

    N=N+2
    ax = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bx = 0.1
    ay = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    by = 0.1
    aw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bw = 0.1

    axy = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ayx = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    axw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ayw = uniform_with_gap(min_value,max_value, min_gap, max_gap)

    if non_linear:
        f = functions_set[random.randint(0, len(functions_set)-1)]
        g = functions_set[random.randint(0, len(functions_set)-1)]
        h = functions_set[random.randint(0, len(functions_set) - 1)]
    else:
        f = identity
        g = identity
        h = identity

    print("Mediator: 1 <- 0 -> 2 <- 1")
    print(f, g, h, ax, bx, ay,by, aw, bw, ayx, axy, axw)
    epsx = np.random.randn(N) ** 3
    epsy = np.random.randn(N) ** 3
    epsw = np.random.randn(N) ** 3

    x = np.zeros([N])
    y = np.zeros([N])
    w = np.zeros([N])

    x[0] = bx * epsx[0]
    x[1] = ax * x[0] + bx * epsx[1]
    y[1] = ay * y[0] + axy * x[0] + by * epsy[1]
    x[2] = ax * x[1] + ayx * y[1] + bx * epsx[2]
    y[2] = ay * y[1] + axy * x[1] + by * epsy[2]
    w[2] = aw * w[1] + axw * x[0] + ayw * y[1] + bw * epsw[2]
    for i in range(3, N):
        x[i] = ax * x[i - 1] + bx * epsx[i]
        y[i] = ay * y[i - 1] + axy * f(x[i - 1]) + by * epsy[i]
        w[i] = aw * w[i - 1] + axw * g(x[i - 2]) + ayw * h(y[i - 1]) + bw * epsw[i]

    x = pd.DataFrame(x, columns=["V1"])
    y = pd.DataFrame(y, columns=["V2"])
    w = pd.DataFrame(w, columns=["V3"])

    series = pd.concat([x, y, w], axis=1, sort=False)
    series = series.drop(series.index[[0,1]])
    series = series.reset_index(drop=True)
    series.index.names = ['time_index']
    return series, summary, temporal


def diamond_generator(N=1000, non_linear=True, min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    summary = np.diag(np.diag(np.ones([4,4])))
    summary[0,1] = 2
    summary[1,0] = 1
    summary[0,2] = 2
    summary[2,0] = 1
    summary[1, 3] = 2
    summary[3, 1] = 1
    summary[2, 3] = 2
    summary[3, 2] = 1
    temporal = dict()
    temporal[0] = [(0, -1)]
    temporal[1] = [(1, -1), (0, -1)]
    temporal[2] = [(2, -1), (0, -2)]
    temporal[3] = [(3, -1), (1, -1), (2, -1)]

    N=N+3
    ax = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bx = 0.1
    ay = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    by = 0.1
    aw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bw = 0.1
    az = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bz = 0.1

    axy = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    axw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ayz = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    awz = uniform_with_gap(min_value,max_value, min_gap, max_gap)

    if non_linear:
        f = functions_set[random.randint(0, len(functions_set)-1)]
        g = functions_set[random.randint(0, len(functions_set)-1)]
        h = functions_set[random.randint(0, len(functions_set)-1)]
        k = functions_set[random.randint(0, len(functions_set)-1)]
    else:
        f = identity
        g = identity
        h = identity
        k = identity

    print("Diamond: 3 <- 1 <- 0 -> 2 -> 3")
    print(f, g, h, k, ax, bx, ay,by, aw, bw, az, bz, axy, axw, ayz, awz)
    epsx = np.random.randn(N) ** 3
    epsy = np.random.randn(N) ** 3
    epsw = np.random.randn(N) ** 3
    epsz = np.random.randn(N) ** 3

    x = np.zeros([N])
    y = np.zeros([N])
    w = np.zeros([N])
    z = np.zeros([N])

    x[0] = bx * epsx[0]
    x[1] = ax * x[0] + bx * epsx[1]
    y[1] = ay * y[0] + axy * x[0] + by * epsy[1]
    x[2] = ax * x[1] + bx * epsx[2]
    y[2] = ay * y[1] + axy * x[1] + by * epsy[2]
    w[2] = aw * w[1] + axw * x[0] + bw * epsw[2]
    x[3] = ax * x[2] + bx * epsx[3]
    y[3] = ay * y[2] + axy * x[2] + by * epsy[3]
    w[3] = aw * w[2] + axw * x[1] + bw * epsw[3]
    z[3] = az * z[2] + ayz * h(y[2]) + awz * k(w[2]) + bz * epsz[3]
    for i in range(3, N):
        x[i] = ax * x[i - 1] + bx * epsx[i]
        y[i] = ay * y[i - 1] + axy * f(x[i - 1]) + by * epsy[i]
        w[i] = aw * w[i - 1] + axw * g(x[i - 2]) + bw * epsw[i]
        z[i] = az * z[i - 1] + ayz * h(y[i - 1]) + awz * k(w[i - 1]) + bz * epsz[i]

    x = pd.DataFrame(x, columns=["V1"])
    y = pd.DataFrame(y, columns=["V2"])
    w = pd.DataFrame(w, columns=["V3"])
    z = pd.DataFrame(z, columns=["V4"])

    series = pd.concat([x, y, w, z], axis=1, sort=False)
    series = series.drop(series.index[[0,1,2]])
    series = series.reset_index(drop=True)
    series.index.names = ['time_index']
    return series, summary, temporal


def seven_ts2h_generator(N=1000, non_linear=True, min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    summary = np.array([[1, 1, 0, 0, 0, 2, 0],
                      [2, 1, 1, 0, 0, 0, 2],
                      [0, 2, 1, 1, 0, 0, 0],
                      [0, 0, 2, 1, 2, 0, 0],
                      [0, 0, 0, 1, 1, 2, 0],
                      [2, 0, 0, 0, 1, 1, 2],
                      [0, 2, 0, 0, 0, 1, 1]])

    temporal = dict()
    temporal[0] = [(0, -1), (1, -1), (5, -1)]
    temporal[1] = [(1, -1), (0, -1), (2, -1)]
    temporal[2] = [(2, -1), (3, -1)]
    temporal[3] = [(3, -1)]
    temporal[4] = [(4, -1), (3, -1)]
    temporal[5] = [(5, -1), (4, -1), (6, -1)]
    temporal[6] = [(6, -1), (5, -1), (1, -1)]

    N=N+3
    at1 = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bt1 = 0.1
    at2 = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bt2 = 0.1
    aa = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ba = 0.1
    ab = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bb = 0.1
    ac = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bc = 0.1
    ad = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bd = 0.1
    ae = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    be = 0.1
    af = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bf = 0.1
    ah = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bh = 0.1

    atb = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ata = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    abe = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    acf = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ach = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ada = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ate = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    atd = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    afb = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ahd = uniform_with_gap(min_value,max_value, min_gap, max_gap)

    if non_linear:
        f0 = functions_set[random.randint(0, len(functions_set)-1)]
        f1 = functions_set[random.randint(0, len(functions_set)-1)]
        f2 = functions_set[random.randint(0, len(functions_set)-1)]
        f3 = functions_set[random.randint(0, len(functions_set)-1)]
        f4 = functions_set[random.randint(0, len(functions_set)-1)]
        f5 = functions_set[random.randint(0, len(functions_set)-1)]
        f6 = functions_set[random.randint(0, len(functions_set)-1)]
        f7 = functions_set[random.randint(0, len(functions_set)-1)]
        f8 = functions_set[random.randint(0, len(functions_set)-1)]
        f9 = functions_set[random.randint(0, len(functions_set)-1)]
    else:
        f0 = identity
        f1 = identity
        f2 = identity
        f3 = identity
        f4 = identity
        f5 = identity
        f6 = identity
        f7 = identity
        f8 = identity
        f9 = identity

    print("complex structure with two hidden causes")
    print(f8, f6, f2, f4, f9, f7, f0, f1, f3, f5, at1, bt1, at2, bt2, ae, be, ab, bb, af, bf, ac, bc, ah, bh, ad, bd,
          aa, ba, ate, atb, atd, ate, abe, afb, acf, ach, ahd, ada)

    epst1 = np.random.randn(N) ** 3
    epst2 = np.random.randn(N) ** 3
    epsa = np.random.randn(N) ** 3
    epsb = np.random.randn(N) ** 3
    epsc = np.random.randn(N) ** 3
    epsd = np.random.randn(N) ** 3
    epse = np.random.randn(N) ** 3
    epsf = np.random.randn(N) ** 3
    epsh = np.random.randn(N) ** 3

    t1 = np.zeros([N])
    t2 = np.zeros([N])
    a = np.zeros([N])
    b = np.zeros([N])
    c = np.zeros([N])
    d = np.zeros([N])
    e = np.zeros([N])
    f = np.zeros([N])
    h = np.zeros([N])

    t1[0] = bt1 * epst1[0]
    t2[0] = bt2 * epst2[0]
    c[0] = bc * epsc[0]
    f[0] = bf * epsf[0]
    h[0] = bh * epsh[0]
    d[0] = bd * epsd[0]
    a[0] = ba * epsa[0]
    b[0] = bb * epsb[0]
    e[0] = be * epse[0]

    t1[1] = at1 * t1[0] + bt1 * epst1[1]
    t2[1] = at2 * t2[0] + bt2 * epst2[1]
    c[1] = ac * c[0] + bc * epsc[1]
    f[1] = af * f[0] + acf * f0(c[0]) + bf * epsf[1]
    h[1] = ah * h[0] + ach * f1(c[0]) + bh * epsh[1]
    d[1] = ad * d[0] + atd * f2(t2[0]) + ahd * f3(h[0]) + bd * epsd[1]
    a[1] = aa * a[0] + ata * f4(t1[0]) + ada * f5(d[0]) + ba * epsa[1]
    b[1] = ab * b[0] + atb * f6(t1[0]) + afb * f7(f[0]) + bb * epsb[1]
    e[1] = ae * e[0] + ate * f8(t2[0]) + abe * f9(b[0]) + be * epse[1]

    for i in range(2, N):
        t1[i] = at1 * t1[i - 1] + bt1 * epst1[i]
        t2[i] = at2 * t2[i - 1] + bt2 * epst2[i]
        c[i] = ac * c[i - 1] + bc * epsc[i]

        f[i] = af * f[i - 1] + acf * f0(c[i - 1]) + bf * epsf[i]
        h[i] = ah * h[i - 1] + ach * f1(c[i - 1]) + bh * epsh[i]
        d[i] = ad * d[i - 1] + atd * f2(t2[i - 1]) + ahd * f3(h[i - 1]) + bd * epsd[i]

        a[i] = aa * a[i - 1] + ata * f4(t1[i - 1]) + ada * f5(d[i - 1]) + ba * epsa[i]
        b[i] = ab * b[i - 1] + atb * f6(t1[i - 1]) + afb * f7(f[i - 1]) + bb * epsb[i]

        e[i] = ae * e[i - 1] + ate * f8(t2[i - 1]) + abe * f9(b[i - 1]) + be * epse[i]

    a = pd.DataFrame(a, columns=["A"])
    b = pd.DataFrame(b, columns=["B"])
    f = pd.DataFrame(f, columns=["F"])
    c = pd.DataFrame(c, columns=["C"])
    h = pd.DataFrame(h, columns=["H"])
    d = pd.DataFrame(d, columns=["D"])
    e = pd.DataFrame(e, columns=["E"])

    series = pd.concat([e, b, f, c, h, d, a], axis=1, sort=False)
    series = series.drop(series.index[[0,1,2]])
    series = series.reset_index(drop=True)
    series.index.names = ['time_index']
    return series, summary, temporal


def seven_ts0h_generator(N=1000, non_linear=True, min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    summary = np.array([[1, 1, 0, 0, 0, 0, 0],
                      [2, 1, 1, 0, 0, 0, 0],
                      [0, 2, 1, 1, 0, 0, 0],
                      [0, 0, 2, 1, 2, 0, 0],
                      [0, 0, 0, 1, 1, 2, 0],
                      [0, 0, 0, 0, 1, 1, 2],
                      [0, 0, 0, 0, 0, 1, 1]])

    temporal = dict()
    temporal[0] = [(0, -1), (1, -1), (5, -1)]
    temporal[1] = [(1, -1), (0, -1), (2, -1)]
    temporal[2] = [(2, -1), (3, -1)]
    temporal[3] = [(3, -1)]
    temporal[4] = [(4, -1), (3, -1)]
    temporal[5] = [(5, -1), (4, -1), (6, -1)]
    temporal[6] = [(6, -1), (5, -1), (1, -1)]

    N=N+3
    aa = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ba = 0.1
    ab = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bb = 0.1
    ac = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bc = 0.1
    ad = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bd = 0.1
    ae = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    be = 0.1
    af = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bf = 0.1
    ah = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bh = 0.1

    abe = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    acf = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ach = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ada = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    afb = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ahd = uniform_with_gap(min_value,max_value, min_gap, max_gap)

    if non_linear:
        f0 = functions_set[random.randint(0, len(functions_set)-1)]
        f1 = functions_set[random.randint(0, len(functions_set)-1)]
        f3 = functions_set[random.randint(0, len(functions_set)-1)]
        f5 = functions_set[random.randint(0, len(functions_set)-1)]
        f7 = functions_set[random.randint(0, len(functions_set)-1)]
        f9 = functions_set[random.randint(0, len(functions_set)-1)]
    else:
        f0 = identity
        f1 = identity
        f3 = identity
        f5 = identity
        f7 = identity
        f9 = identity

    print("7 time series structure with no hidden causes")
    print(f9, f7, f0, f1, f3, f5, ae, be, ab, bb, af, bf, ac, bc, ah, bh, ad, bd, aa, ba, abe, afb, acf, ach, ahd, ada)
    epsa = np.random.randn(N) ** 3
    epsb = np.random.randn(N) ** 3
    epsc = np.random.randn(N) ** 3
    epsd = np.random.randn(N) ** 3
    epse = np.random.randn(N) ** 3
    epsf = np.random.randn(N) ** 3
    epsh = np.random.randn(N) ** 3

    a = np.zeros([N])
    b = np.zeros([N])
    c = np.zeros([N])
    d = np.zeros([N])
    e = np.zeros([N])
    f = np.zeros([N])
    h = np.zeros([N])

    c[0] = bc * epsc[0]
    f[0] = bf * epsf[0]
    h[0] = bh * epsh[0]
    d[0] = bd * epsd[0]
    a[0] = ba * epsa[0]
    b[0] = bb * epsb[0]
    e[0] = be * epse[0]

    c[1] = ac * c[0] + bc * epsc[1]
    f[1] = af * f[0] + acf * f0(c[0]) + bf * epsf[1]
    h[1] = ah * h[0] + ach * f1(c[0]) + bh * epsh[1]
    d[1] = ad * d[0] + ahd * f3(h[0]) + bd * epsd[1]
    a[1] = aa * a[0] + ada * f5(d[0]) + ba * epsa[1]
    b[1] = ab * b[0] + afb * f7(f[0]) + bb * epsb[1]
    e[1] = ae * e[0] + abe * f9(b[0]) + be * epse[1]

    for i in range(2, N):
        c[i] = ac * c[i - 1] + bc * epsc[i]

        f[i] = af * f[i - 1] + acf * f0(c[i - 1]) + bf * epsf[i]
        h[i] = ah * h[i - 1] + ach * f1(c[i - 1]) + bh * epsh[i]
        d[i] = ad * d[i - 1] + ahd * f3(h[i - 1]) + bd * epsd[i]

        a[i] = aa * a[i - 1] + ada * f5(d[i - 1]) + ba * epsa[i]
        b[i] = ab * b[i - 1] + afb * f7(f[i - 1]) + bb * epsb[i]

        e[i] = ae * e[i - 1] + abe * f9(b[i - 1]) + be * epse[i]

    a = pd.DataFrame(a, columns=["A"])
    b = pd.DataFrame(b, columns=["B"])
    f = pd.DataFrame(f, columns=["F"])
    c = pd.DataFrame(c, columns=["C"])
    h = pd.DataFrame(h, columns=["H"])
    d = pd.DataFrame(d, columns=["D"])
    e = pd.DataFrame(e, columns=["E"])

    series = pd.concat([e, b, f, c, h, d, a], axis=1, sort=False)
    series = series.drop(series.index[[0, 1, 2]])
    series = series.reset_index(drop=True)
    series.index.names = ['time_index']
    return series, summary, temporal


################################

def fork_generator_diff_sampling_rate(N=1000, non_linear=True, min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    summary = np.diag(np.diag(np.ones([3,3])))
    summary[0,1] = 2
    summary[1,0] = 1
    summary[0,2] = 2
    summary[2,0] = 1
    temporal = dict()
    temporal[0] = [(0, -1)]
    temporal[1] = [(1, -1), (0, -1)]
    temporal[2] = [(2, -1), (0, -2)]

    N=N+2
    ax = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bx = 0.1
    ay = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    by = 0.1
    aw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bw = 0.1

    axy = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    axw = uniform_with_gap(min_value,max_value, min_gap, max_gap)

    if non_linear:
        f = functions_set[random.randint(0, len(functions_set)-1)]
        g = functions_set[random.randint(0, len(functions_set)-1)]
    else:
        f = identity
        g = identity

    print("Fork: 1 <- 0 -> 2")
    print(f, g, ax, bx, ay,by, aw, bw, axy, axw)
    epsx = np.random.randn(N) ** 3
    epsy = np.random.randn(N) ** 3
    epsw = np.random.randn(N) ** 3

    x = np.zeros([N])
    y = np.zeros([N])
    w = np.zeros([N])

    x[0] = bx * epsx[0]
    x[1] = ax * x[0] + bx * epsx[1]
    y[1] = ay * y[0] + axy * x[0] + by * epsy[1]
    w[1] = np.nan

    for i in range(2, N):
        x[i] = ax * x[i - 1] + bx * epsx[i]
        y[i] = ay * y[i - 1] + axy * f(x[i - 1]) + by * epsy[i]
        if i % 2 == 0:
            w[i] = aw * w[i - 2] + axw * g(x[i - 2]) + bw * epsw[i]
        else:
            w[i] = np.nan

    x = pd.DataFrame(x, columns=["V1"])
    y = pd.DataFrame(y, columns=["V2"])
    w = pd.DataFrame(w, columns=["V3"])

    series = pd.concat([x, y, w], axis=1, sort=False)
    series = series.drop(series.index[[0,1]])
    series = series.reset_index(drop=True)
    series.index.names = ['time_index']
    return series, summary, temporal


def v_structure_generator_diff_sampling_rate(N=1000, non_linear=True, min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    summary = np.diag(np.diag(np.ones([3,3])))
    summary[0,2] = 2
    summary[2,0] = 1
    summary[1, 2] = 2
    summary[2,1] = 1
    temporal = dict()
    temporal[0] = [(0, -1)]
    temporal[1] = [(1, -1)]
    temporal[2] = [(2, -1), (0, -2), (1, -1)]

    N=N+2
    ax = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bx = 0.1
    ay = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    by = 0.1
    aw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bw = 0.1

    axw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ayw = uniform_with_gap(min_value,max_value, min_gap, max_gap)

    if non_linear:
        f = functions_set[random.randint(0, len(functions_set)-1)]
        g = functions_set[random.randint(0, len(functions_set)-1)]
    else:
        f = identity
        g = identity

    print("V-structure: 0 -> 2 <- 1")
    print(f, g, ax, bx, ay,by, aw, bw, axw, ayw)
    epsx = np.random.randn(N) ** 3
    epsy = np.random.randn(N) ** 3
    epsw = np.random.randn(N) ** 3

    x = np.zeros([N])
    y = np.zeros([N])
    w = np.zeros([N])

    x[0] = bx * epsx[0]
    y[0] = by * epsy[0]
    w[0] = np.nan
    x[1] = ax * x[0] + bx * epsx[1]
    y[1] = np.nan
    x[2] = ax * x[1] + bx * epsx[2]
    y[2] = ay * y[0] + by * epsy[2]
    w[2] = np.nan
    for i in range(3, N):
        x[i] = ax * x[i - 1] + bx * epsx[i]
        if i % 2 == 0:
            y[i] = ay * y[i - 2] + by * epsy[i]
            w[i] = np.nan
        else:
            y[i] = np.nan
            w[i] = aw * w[i - 2] + axw * f(x[i - 2]) + ayw * g(y[i - 1]) + bw * epsw[i]

    x = pd.DataFrame(x, columns=["V1"])
    y = pd.DataFrame(y, columns=["V2"])
    w = pd.DataFrame(w, columns=["V3"])

    series = pd.concat([x, y, w], axis=1, sort=False)
    series = series.drop(series.index[[0,1]])
    series = series.reset_index(drop=True)
    series.index.names = ['time_index']
    return series, summary, temporal


def mediator_generator_diff_sampling_rate(N=1000, non_linear=True, min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    summary = np.diag(np.diag(np.ones([3,3])))
    summary[0,1] = 2
    summary[1,0] = 1
    summary[0,2] = 2
    summary[2,0] = 1
    summary[1,2] = 2
    summary[2,1] = 1

    temporal = dict()
    temporal[0] = [(0, -1)]
    temporal[1] = [(1, -1), (0, -1)]
    temporal[2] = [(2, -1), (0, -2), (1, -1)]

    N=N+2
    ax = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bx = 0.1
    ay = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    by = 0.1
    aw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bw = 0.1

    axy = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    axw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ayw = uniform_with_gap(min_value,max_value, min_gap, max_gap)

    if non_linear:
        f = functions_set[random.randint(0, len(functions_set)-1)]
        g = functions_set[random.randint(0, len(functions_set)-1)]
        h = functions_set[random.randint(0, len(functions_set)-1)]
    else:
        f = identity
        g = identity
        h = identity

    print("Mediator: 0 -> 2 <- 1 <-0")
    print(f, g, ax, bx, ay,by, aw, bw, axy, axw)
    epsx = np.random.randn(N) ** 3
    epsy = np.random.randn(N) ** 3
    epsw = np.random.randn(N) ** 3

    x = np.zeros([N])
    y = np.zeros([N])
    w = np.zeros([N])

    x[0] = bx * epsx[0]
    x[1] = ax * x[0] + bx * epsx[1]
    y[1] = ay * y[0] + axy * x[0] + by * epsy[1]
    w[1] = np.nan
    for i in range(2, N):
        x[i] = ax * x[i - 1] + bx * epsx[i]
        y[i] = ay * y[i - 1] + axy * f(x[i - 1]) + by * epsy[i]
        if i % 2 == 0:
            w[i] = aw * w[i - 2] + axw * g(x[i - 2]) + ayw * h(y[i - 1]) + bw * epsw[i]
        else:
            w[i] = np.nan

    x = pd.DataFrame(x, columns=["V1"])
    y = pd.DataFrame(y, columns=["V2"])
    w = pd.DataFrame(w, columns=["V3"])

    series = pd.concat([x, y, w], axis=1, sort=False)
    series = series.drop(series.index[[0,1]])
    series = series.reset_index(drop=True)
    series.index.names = ['time_index']
    return series, summary, temporal


def diamond_generator_diff_sampling_rate(N=1000, non_linear=True, min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    # summary = np.zeros([4,4])
    summary = np.diag(np.diag(np.ones([4,4])))
    summary[0,1] = 2
    summary[1,0] = 1
    summary[0,2] = 2
    summary[2,0] = 1
    summary[1, 3] = 2
    summary[3, 1] = 1
    summary[2, 3] = 2
    summary[3, 2] = 1
    temporal = dict()
    temporal[0] = [(0, -1)]
    temporal[1] = [(1, -1), (0, -1)]
    temporal[2] = [(2, -1), (0, -2)]
    temporal[3] = [(3, -1), (1, -1), (2, -1)]

    N=N+3
    ax = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bx = 0.1
    ay = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    by = 0.1
    aw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bw = 0.1
    az = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bz = 0.1

    axy = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    axw = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ayz = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    awz = uniform_with_gap(min_value,max_value, min_gap, max_gap)

    if non_linear:
        f = functions_set[random.randint(0, len(functions_set)-1)]
        g = functions_set[random.randint(0, len(functions_set)-1)]
        h = functions_set[random.randint(0, len(functions_set)-1)]
        k = functions_set[random.randint(0, len(functions_set)-1)]
    else:
        f = identity
        g = identity
        h = identity
        k = identity

    print("Diamond: 3 <- 1 <- 0 -> 2 -> 3")
    print(f, g, h, k, ax, bx, ay,by, aw, bw, az, bz, axy, axw, ayz, awz)
    epsx = np.random.randn(N) ** 3
    epsy = np.random.randn(N) ** 3
    epsw = np.random.randn(N) ** 3
    epsz = np.random.randn(N) ** 3

    x = np.zeros([N])
    y = np.zeros([N])
    w = np.zeros([N])
    z = np.zeros([N])

    x[0] = bx * epsx[0]
    w[0] = np.nan
    x[1] = ax * x[0] + bx * epsx[1]
    # y[1] = ay * y[0] + axy * x[0] + by * epsy[1]
    y[1] = np.nan
    x[2] = ax * x[1] + bx * epsx[2]
    y[2] = ay * y[2] + axy * x[1] + by * epsy[2]
    # w[2] = aw * w[1] + axw * x[0] + bw * epsw[2]
    w[2] = np.nan
    x[3] = ax * x[2] + bx * epsx[3]
    # y[3] = ay * y[1] + axy * x[2] + by * epsy[3]
    y[3] = np.nan
    w[3] = aw * w[1] + axw * x[1] + bw * epsw[3]
    z[3] = az * z[1] + ayz * h(y[2]) + awz * k(w[1]) + bz * epsz[3]
    for i in range(4, N):
        x[i] = ax * x[i - 1] + bx * epsx[i]
        # y[i] = ay * y[i - 1] + axy * f(x[i - 1]) + by * epsy[i]
        # w[i] = aw * w[i - 1] + axw * g(x[i - 2]) + bw * epsw[i]
        if i % 2 == 0:
            y[i] = ay * y[i - 2] + axy * f(x[i - 1]) + by * epsy[i]
            w[i] = np.nan
        else:
            y[i] = np.nan
            w[i] = aw * w[i - 2] + axw * g(x[i - 2]) + bw * epsw[i]
        z[i] = az * z[i - 2] + ayz * h(y[i - 2]) + awz * k(w[i - 1]) + bz * epsz[i]

    x = pd.DataFrame(x, columns=["V1"])
    y = pd.DataFrame(y, columns=["V2"])
    w = pd.DataFrame(w, columns=["V3"])
    z = pd.DataFrame(z, columns=["V4"])

    series = pd.concat([x, y, w, z], axis=1, sort=False)
    series = series.drop(series.index[[0,1,2,3,4]])
    series = series.reset_index(drop=True)
    series.index.names = ['time_index']
    return series, summary, temporal


def seven_ts0h_generator_diff_sampling_rate(N=1000, non_linear=True, min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    summary = np.array([[1, 1, 0, 0, 0, 0, 0],
                      [2, 1, 1, 0, 0, 0, 0],
                      [0, 2, 1, 1, 0, 0, 0],
                      [0, 0, 2, 1, 2, 0, 0],
                      [0, 0, 0, 1, 1, 2, 0],
                      [0, 0, 0, 0, 1, 1, 2],
                      [0, 0, 0, 0, 0, 1, 1]])

    temporal = dict()
    temporal[0] = [(0, -1), (1, -1), (5, -1)]
    temporal[1] = [(1, -1), (0, -1), (2, -1)]
    temporal[2] = [(2, -1), (3, -1)]
    temporal[3] = [(3, -1)]
    temporal[4] = [(4, -1), (3, -1)]
    temporal[5] = [(5, -1), (4, -1), (6, -1)]
    temporal[6] = [(6, -1), (5, -1), (1, -1)]

    N=N+3
    # at1 = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    # bt1 = 0.1
    # at2 = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    # bt2 = 0.1
    aa = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ba = 0.1
    ab = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bb = 0.1
    ac = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bc = 0.1
    ad = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bd = 0.1
    ae = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    be = 0.1
    af = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bf = 0.1
    ah = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bh = 0.1

    # atb = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    # ata = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    abe = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    acf = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ach = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ada = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    # ate = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    # atd = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    afb = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ahd = uniform_with_gap(min_value,max_value, min_gap, max_gap)

    if non_linear:
        f0 = functions_set[random.randint(0, len(functions_set)-1)]
        f1 = functions_set[random.randint(0, len(functions_set)-1)]
        # f2 = functions_set[random.randint(0, len(functions_set)-1)]
        f3 = functions_set[random.randint(0, len(functions_set)-1)]
        # f4 = functions_set[random.randint(0, len(functions_set)-1)]
        f5 = functions_set[random.randint(0, len(functions_set)-1)]
        # f6 = functions_set[random.randint(0, len(functions_set)-1)]
        f7 = functions_set[random.randint(0, len(functions_set)-1)]
        # f8 = functions_set[random.randint(0, len(functions_set)-1)]
        f9 = functions_set[random.randint(0, len(functions_set)-1)]
    else:
        f0 = identity
        f1 = identity
        # f2 = identity
        f3 = identity
        # f4 = identity
        f5 = identity
        # f6 = identity
        f7 = identity
        # f8 = identity
        f9 = identity

    print("7 time series structure with no hidden causes")
    print(f9, f7, f0, f1, f3, f5, ae, be, ab, bb, af, bf, ac, bc, ah, bh, ad, bd, aa, ba, abe, afb, acf, ach, ahd, ada)
    # epst1 = np.random.randn(N) ** 3
    # epst2 = np.random.randn(N) ** 3
    epsa = np.random.randn(N) ** 3
    epsb = np.random.randn(N) ** 3
    epsc = np.random.randn(N) ** 3
    epsd = np.random.randn(N) ** 3
    epse = np.random.randn(N) ** 3
    epsf = np.random.randn(N) ** 3
    epsh = np.random.randn(N) ** 3

    # t1 = np.zeros([N])
    # t2 = np.zeros([N])
    a = np.zeros([N])
    b = np.zeros([N])
    c = np.zeros([N])
    d = np.zeros([N])
    e = np.zeros([N])
    f = np.zeros([N])
    h = np.zeros([N])

    # t1[0] = bt1 * epst1[0]
    # t2[0] = bt2 * epst2[0]
    c[0] = bc * epsc[0]
    f[0] = bf * epsf[0]
    h[0] = bh * epsh[0]
    d[0] = bd * epsd[0]
    a[0] = ba * epsa[0]
    b[0] = bb * epsb[0]
    e[0] = be * epse[0]

    # t1[1] = at1 * t1[0] + bt1 * epst1[1]
    # t2[1] = at2 * t2[0] + bt2 * epst2[1]
    c[1] = ac * c[0] + bc * epsc[1]
    f[1] = af * f[0] + acf * f0(c[0]) + bf * epsf[1]
    h[1] = ah * h[0] + ach * f1(c[0]) + bh * epsh[1]
    d[1] = ad * d[0] + ahd * f3(h[0]) + bd * epsd[1]
    a[1] = aa * a[0] + ada * f5(d[0]) + ba * epsa[1]
    b[1] = ab * b[0] + afb * f7(f[0]) + bb * epsb[1]
    e[1] = ae * e[0] + abe * f9(b[0]) + be * epse[1]

    for i in range(2, N):
        # t1[i] = at1 * t1[i - 1] + bt1 * epst1[i]
        # t2[i] = at2 * t2[i - 1] + bt2 * epst2[i]
        c[i] = ac * c[i - 1] + bc * epsc[i]

        f[i] = af * f[i - 1] + acf * f0(c[i - 1]) + bf * epsf[i]
        h[i] = ah * h[i - 1] + ach * f1(c[i - 1]) + bh * epsh[i]
        d[i] = ad * d[i - 1] + ahd * f3(h[i - 1]) + bd * epsd[i]

        if i % 2 == 0:
            a[i] = aa * a[i - 2] + ada * f5(d[i - 1]) + ba * epsa[i]
            b[i] = ab * b[i - 2] + afb * f7(f[i - 1]) + bb * epsb[i]
            e[i] = np.nan
        else:
            a[i] = np.nan
            b[i] = np.nan
            e[i] = ae * e[i - 2] + abe * f9(b[i - 1]) + be * epse[i]

    a = pd.DataFrame(a, columns=["A"])
    b = pd.DataFrame(b, columns=["B"])
    f = pd.DataFrame(f, columns=["F"])
    c = pd.DataFrame(c, columns=["C"])
    h = pd.DataFrame(h, columns=["H"])
    d = pd.DataFrame(d, columns=["D"])
    e = pd.DataFrame(e, columns=["E"])

    series = pd.concat([e, b, f, c, h, d, a], axis=1, sort=False)
    series = series.drop(series.index[[0, 1, 2]])
    series = series.reset_index(drop=True)
    series.index.names = ['time_index']
    return series, summary, temporal


def seven_ts2h_generator_diff_sampling_rate(N=1000, non_linear=True, min_value=-1, max_value=1, min_gap=-0.1, max_gap=0.1):
    summary = np.array([[1, 1, 0, 0, 0, 2, 0],
                      [2, 1, 1, 0, 0, 0, 2],
                      [0, 2, 1, 1, 0, 0, 0],
                      [0, 0, 2, 1, 2, 0, 0],
                      [0, 0, 0, 1, 1, 2, 0],
                      [2, 0, 0, 0, 1, 1, 2],
                      [0, 2, 0, 0, 0, 1, 1]])

    temporal = dict()
    temporal[0] = [(0, -1), (1, -1), (5, -1)]
    temporal[1] = [(1, -1), (0, -1), (2, -1)]
    temporal[2] = [(2, -1), (3, -1)]
    temporal[3] = [(3, -1)]
    temporal[4] = [(4, -1), (3, -1)]
    temporal[5] = [(5, -1), (4, -1), (6, -1)]
    temporal[6] = [(6, -1), (5, -1), (1, -1)]

    N=N+3
    at1 = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bt1 = 0.1
    at2 = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bt2 = 0.1
    aa = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ba = 0.1
    ab = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bb = 0.1
    ac = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bc = 0.1
    ad = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bd = 0.1
    ae = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    be = 0.1
    af = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bf = 0.1
    ah = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    bh = 0.1

    atb = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ata = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    abe = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    acf = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ach = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ada = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ate = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    atd = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    afb = uniform_with_gap(min_value,max_value, min_gap, max_gap)
    ahd = uniform_with_gap(min_value,max_value, min_gap, max_gap)

    if non_linear:
        f0 = functions_set[random.randint(0, len(functions_set)-1)]
        f1 = functions_set[random.randint(0, len(functions_set)-1)]
        f2 = functions_set[random.randint(0, len(functions_set)-1)]
        f3 = functions_set[random.randint(0, len(functions_set)-1)]
        f4 = functions_set[random.randint(0, len(functions_set)-1)]
        f5 = functions_set[random.randint(0, len(functions_set)-1)]
        f6 = functions_set[random.randint(0, len(functions_set)-1)]
        f7 = functions_set[random.randint(0, len(functions_set)-1)]
        f8 = functions_set[random.randint(0, len(functions_set)-1)]
        f9 = functions_set[random.randint(0, len(functions_set)-1)]
    else:
        f0 = identity
        f1 = identity
        f2 = identity
        f3 = identity
        f4 = identity
        f5 = identity
        f6 = identity
        f7 = identity
        f8 = identity
        f9 = identity

    print("complex structure with two hidden causes")
    print(f8, f6, f2, f4, f9, f7, f0, f1, f3, f5, at1, bt1, at2, bt2, ae, be, ab, bb, af, bf, ac, bc, ah, bh, ad, bd,
          aa, ba, ate, atb, atd, ate, abe, afb, acf, ach, ahd, ada)

    epst1 = np.random.randn(N) ** 3
    epst2 = np.random.randn(N) ** 3
    epsa = np.random.randn(N) ** 3
    epsb = np.random.randn(N) ** 3
    epsc = np.random.randn(N) ** 3
    epsd = np.random.randn(N) ** 3
    epse = np.random.randn(N) ** 3
    epsf = np.random.randn(N) ** 3
    epsh = np.random.randn(N) ** 3

    t1 = np.zeros([N])
    t2 = np.zeros([N])
    a = np.zeros([N])
    b = np.zeros([N])
    c = np.zeros([N])
    d = np.zeros([N])
    e = np.zeros([N])
    f = np.zeros([N])
    h = np.zeros([N])

    t1[0] = bt1 * epst1[0]
    t2[0] = bt2 * epst2[0]
    c[0] = bc * epsc[0]
    f[0] = bf * epsf[0]
    h[0] = bh * epsh[0]
    d[0] = bd * epsd[0]
    a[0] = ba * epsa[0]
    b[0] = bb * epsb[0]
    e[0] = be * epse[0]

    t1[1] = at1 * t1[0] + bt1 * epst1[1]
    t2[1] = at2 * t2[0] + bt2 * epst2[1]
    c[1] = ac * c[0] + bc * epsc[1]
    f[1] = af * f[0] + acf * f0(c[0]) + bf * epsf[1]
    h[1] = ah * h[0] + ach * f1(c[0]) + bh * epsh[1]
    d[1] = ad * d[0] + atd * f2(t2[0]) + ahd * f3(h[0]) + bd * epsd[1]
    a[1] = aa * a[0] + ata * f4(t1[0]) + ada * f5(d[0]) + ba * epsa[1]
    b[1] = ab * b[0] + atb * f6(t1[0]) + afb * f7(f[0]) + bb * epsb[1]
    e[1] = ae * e[0] + ate * f8(t2[0]) + abe * f9(b[0]) + be * epse[1]

    for i in range(2, N):
        t1[i] = at1 * t1[i - 1] + bt1 * epst1[i]
        t2[i] = at2 * t2[i - 1] + bt2 * epst2[i]
        c[i] = ac * c[i - 1] + bc * epsc[i]

        f[i] = af * f[i - 1] + acf * f0(c[i - 1]) + bf * epsf[i]
        h[i] = ah * h[i - 1] + ach * f1(c[i - 1]) + bh * epsh[i]
        d[i] = ad * d[i - 1] + atd * f2(t2[i - 1]) + ahd * f3(h[i - 1]) + bd * epsd[i]

        if i % 2 == 0:
            a[i] = aa * a[i - 2] + ada * f5(d[i - 1]) + ba * epsa[i]
            b[i] = ab * b[i - 2] + afb * f7(f[i - 1]) + bb * epsb[i]
            e[i] = np.nan
        else:
            a[i] = np.nan
            b[i] = np.nan
            e[i] = ae * e[i - 2] + abe * f9(b[i - 1]) + be * epse[i]

    a = pd.DataFrame(a, columns=["A"])
    b = pd.DataFrame(b, columns=["B"])
    f = pd.DataFrame(f, columns=["F"])
    c = pd.DataFrame(c, columns=["C"])
    h = pd.DataFrame(h, columns=["H"])
    d = pd.DataFrame(d, columns=["D"])
    e = pd.DataFrame(e, columns=["E"])

    series = pd.concat([e, b, f, c, h, d, a], axis=1, sort=False)
    series = series.drop(series.index[[0,1,2]])
    series = series.reset_index(drop=True)
    series.index.names = ['time_index']
    return series, summary, temporal


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 5:
        print(len(sys.argv))
        structure = sys.argv[1]
        n_samples = int(sys.argv[2])
        runs = int(sys.argv[3])
        non_linear = bool(int(sys.argv[4]))
        print('Argument List:', str(sys.argv))
    else:
        print('Missing arguments so will take default arguments')
        structure = "7ts2h"
        n_samples = 4000
        runs = 10
        non_linear = True
        print('Default Argument List:', str(structure), str(n_samples), str(runs))

    save_data = True

    get_data = {"fork": fork_generator, "v_structure": v_structure_generator, "mediator": mediator_generator,
                "diamond": diamond_generator, "7ts0h": seven_ts0h_generator, "7ts2h": seven_ts2h_generator,
                "dsr_fork": fork_generator_diff_sampling_rate,
                "dsr_v_structure": v_structure_generator_diff_sampling_rate,
                "dsr_mediator": mediator_generator_diff_sampling_rate,
                "dsr_diamond": diamond_generator_diff_sampling_rate,
                "dsr_7ts0h": seven_ts0h_generator_diff_sampling_rate,
                "dsr_7ts2h": seven_ts2h_generator_diff_sampling_rate}
    get_seed_constant = {"fork": 0.1, "v_structure": 0.2, "diamond": 0.3, "mediator": 0.7, "7ts0h": 0.5, "7ts2h": 0.5,
                         "dsr_fork": 0.1, "dsr_v_structure": 0.2, "dsr_mediator": 0.7, "dsr_diamond": 0.3,
                         "dsr_7ts0h": 0.5, "dsr_7ts2h": 0.5}

    for i in range(runs):
        random.seed(a=get_seed_constant[structure]*(i+1)/runs)
        np.random.seed((i+1))
        print("############################## Run "+str(i)+" ##############################")
        data, ground_truth_summary, ground_truth_temporal = get_data[structure](N=n_samples, non_linear=non_linear)
        print(data.iloc[:1000])

        if save_data:
            data.to_csv("../data/simulated_ts_data/" + "/" + str(structure) + "/data_" + str(i) + ".csv")
