"""
数值 ODE 考前 30 分钟填空练习

做法：
1. 搜索 TODO，手填关键公式。
2. 填完运行：python ode_fill_blank_practice.py
3. 看到“全部关键练习通过”才算完成。

文件末尾有注释答案区。当前文件故意不能正确通过，填完后可直接运行。
"""

import numpy as np


def euler_practice(f, t0, y0, T, h):
    ts = np.arange(t0, T + h / 2, h)
    ys = np.zeros(len(ts))
    ys[0] = y0
    for n in range(len(ts) - 1):
        # TODO 1：写 Euler 更新公式 y[n+1] = ...
        ys[n + 1] = np.nan
    return ts, ys


def rk4_practice(f, t0, y0, T, h):
    ts = np.arange(t0, T + h / 2, h)
    ys = np.zeros(len(ts))
    ys[0] = y0
    for n in range(len(ts) - 1):
        t, y = ts[n], ys[n]
        # TODO 2：补全 k1,k2,k3,k4
        k1 = np.nan
        k2 = np.nan
        k3 = np.nan
        k4 = np.nan
        ys[n + 1] = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return ts, ys


def rkf45_accept_reject_practice(t, y, h, y5, err, tol, ts, ys, hs):
    """只练最关键的接受/拒绝逻辑，不要求重抄 RKF45 全部系数。"""
    # TODO 3：
    # err <= tol 时：更新 t,y，并 append 到 ts,ys,hs。
    # err > tol 时：不许更新 t,y，也不许 append。
    pass
    return t, y


def second_order_system_practice(t, w):
    """题目：y''+y=0。令 w[0]=y, w[1]=y'。"""
    # TODO 4：返回 [w1', w2'] = [y', -y]
    return np.array([np.nan, np.nan])


def observed_order_practice(error_h, error_h2):
    """已知 E(h) 和 E(h/2)，返回观测收敛阶。"""
    # TODO 5：写 p=log(E(h)/E(h/2))/log(2)
    return np.nan


def check_answers():
    f = lambda t, y: y

    _, y_eu = euler_practice(f, 0.0, 1.0, 0.2, 0.1)
    assert np.allclose(y_eu, [1.0, 1.1, 1.21]), "Euler TODO 仍不正确"

    _, y_rk = rk4_practice(f, 0.0, 1.0, 0.1, 0.1)
    assert abs(y_rk[-1] - np.exp(0.1)) < 1e-7, "RK4 四个 k 仍不正确"

    ts, ys, hs = [0.0], [1.0], []
    t, y = rkf45_accept_reject_practice(0.0, 1.0, 0.1, 1.1, 1e-7, 1e-6,
                                         ts, ys, hs)
    assert (t, y) == (0.1, 1.1) and hs == [0.1], "RKF45 接受逻辑不正确"
    old = (t, y, list(ts), list(ys), list(hs))
    t, y = rkf45_accept_reject_practice(t, y, 0.2, 1.4, 1e-3, 1e-6,
                                         ts, ys, hs)
    assert (t, y, ts, ys, hs) == old, "RKF45 拒绝步不应更新或保存"

    assert np.allclose(second_order_system_practice(0.0, np.array([2.0, 3.0])),
                       [3.0, -2.0]), "二阶化系统 TODO 仍不正确"
    assert abs(observed_order_practice(0.04, 0.01) - 2.0) < 1e-12, "收敛阶 TODO 不正确"

    print("全部关键练习通过：Euler、RK4、RKF45 接受/拒绝、二阶系统、收敛阶。")


if __name__ == "__main__":
    check_answers()


# ============================================================================
# 答案区：默认注释。实在卡住再看。
# ============================================================================
#
# TODO 1:
# ys[n + 1] = ys[n] + h * f(ts[n], ys[n])
#
# TODO 2:
# k1 = f(t, y)
# k2 = f(t + h/2, y + h*k1/2)
# k3 = f(t + h/2, y + h*k2/2)
# k4 = f(t + h, y + h*k3)
#
# TODO 3:
# if err <= tol:
#     t += h
#     y = y5
#     ts.append(t)
#     ys.append(y)
#     hs.append(h)
#
# TODO 4:
# return np.array([w[1], -w[0]])
#
# TODO 5:
# return np.log(error_h / error_h2) / np.log(2)
