import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.io

# Constantes
g = 9.81
l1 = 2
l2 = 2.5
l2e = 1.8
m1 = 450
m2 = 650
F1 = -0.5 * m1 * g
F2 = -0.5 * m2 * g
uIz = 2.7
r = 0.3
xd_p = 22.2 # 80 km/h = 22.2 m/s

# Condições Iniciais
theta1 = 0
theta2 = 0
theta1_p = 0.4
theta2_p = -0.1
condicoes_iniciais = np.array([theta1, theta2, theta1_p, theta2_p])

# Parâmetros de simulação
H = [0.01, 0.005, 0.001]
# H = [0.05, 0.01, 0.005, 0.001]
tempo = 60

######
# EDOs


def fa(estado):
    a0 = math.pow(l1, 2)*l2*r * \
        (m2*math.cos(2*estado[0] - 2*estado[1]) - 2*m1 - m2)
    a1 = math.pow(l1, 2)*l2*r*m2*math.sin(2*estado[0] - 2*estado[1])
    a2 = 2*l1*math.pow(l2, 2)*r*m2*math.sin(estado[0] - estado[1])
    a3 = -2*l2*uIz*xd_p
    a4 = -2*l1*uIz*xd_p*math.cos(estado[0] - estado[1])
    a5 = -r*l1*(l2e*F2*math.sin(estado[0] - 2*estado[1]) +
                2*math.sin(estado[0])*(F1*l2 + (l2e*F2)/2))

    return (a1*math.pow(estado[2], 2) + a2*math.pow(estado[3], 2) + a3*estado[2] + a4*estado[3] + a5)/a0


def fb(estado, aceleracao):
    b0 = math.pow(l2, 2)*r*m2
    b1 = -l1*l2*r*m2*math.cos(estado[0] - estado[1])
    b2 = l1*l2*r*m2*math.sin(estado[0] - estado[1])
    b3 = -uIz*xd_p
    b4 = l2e*math.sin(estado[1])*r*F2

    return (b1*aceleracao[0] + b2*math.pow(estado[2], 2) + b3*estado[3] + b4)/b0


def edo(y, aceleracao):
    dy = np.array([y[2], y[3], fa(y), fb(y, aceleracao)])
    return dy


##########
# Soluções
def euler(c_i, fa, fb, t, h):
    passos = math.floor(t / h)
    estados = np.array([[0.0]*4]*passos)
    aceleracoes = np.array([[0.0]*2]*passos)

    estados[0] = c_i
    aceleracoes[0][0] = fa(c_i)
    aceleracoes[0][1] = fb(c_i, aceleracoes[0])

    for n in range(1, passos):
        dy = edo(estados[n-1], aceleracoes[n-1])
        estados[n] = estados[n-1] + h*dy
        aceleracoes[n][0] = dy[2]  # theta1_pp
        aceleracoes[n][1] = dy[3]  # theta2_pp

    return estados, aceleracoes


def RK2(c_i, fa, fb, t, h):
    passos = math.floor(t / h)
    estados = np.array([[0.0]*4]*passos)
    aceleracoes = np.array([[0.0]*2]*passos)

    estados[0] = c_i
    aceleracoes[0][0] = fa(c_i)
    aceleracoes[0][1] = fb(c_i, aceleracoes[0])

    for n in range(1, passos):
        k1 = edo(estados[n-1], aceleracoes[n-1])
        k2 = edo(estados[n-1] + (h/2)*k1, aceleracoes[n-1])
        estados[n] = estados[n-1] + h*k2
        aceleracoes[n][0] = k2[2]  # theta1_pp
        aceleracoes[n][1] = k2[3]  # theta2_pp

    return estados, aceleracoes

def RK4(c_i, fa, fb, t, h):
    passos = math.floor(t / h)
    estados = np.array([[0.0]*4]*passos)
    aceleracoes = np.array([[0.0]*2]*passos)

    estados[0] = c_i
    aceleracoes[0][0] = fa(c_i)
    aceleracoes[0][1] = fb(c_i, aceleracoes[0])

    for n in range(1, passos):
        k1 = edo(estados[n-1], aceleracoes[n-1])
        k2 = edo(estados[n-1] + (h/2)*k1, aceleracoes[n-1])
        k3 = edo(estados[n-1] + (h/2)*k2, aceleracoes[n-1])
        k4 = edo(estados[n-1] + h*k3, aceleracoes[n-1])
        estados[n] = estados[n-1] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        aceleracoes[n][0] = k4[2]  # theta1_pp
        aceleracoes[n][1] = k4[3]  # theta2_pp

    return estados, aceleracoes


######################
# Plotagem de gráficos
def graficos(e, a, h, titulo):
    plt.suptitle(titulo, fontsize=16)

    ax = plt.subplot(321)
    plt.plot(np.arange(0.0, tempo, h), e[:, 0])
    ax.set(xlabel='tempo (s)', ylabel=r'$\theta_1$ (rad)',
           title=r'$\theta_1$ com passo h =' + str(h))
    ax.grid()

    ax = plt.subplot(322)
    plt. plot(np.arange(0.0, tempo, h), e[:, 1])
    ax.set(xlabel='tempo (s)', ylabel=r'$\theta_2$ (rad)',
           title=r'$\theta_2$ com passo h =' + str(h))
    ax.grid()

    ax = plt.subplot(323)
    plt. plot(np.arange(0.0, tempo, h), e[:, 2])
    ax.set(xlabel='tempo (s)', ylabel=r'$\dot \theta_1$ (rad/s)',
           title=r'$\dot \theta_1$ com passo h =' + str(h))
    ax.grid()

    ax = plt.subplot(324)
    plt. plot(np.arange(0.0, tempo, h), e[:, 3])
    ax.set(xlabel='tempo (s)', ylabel=r'$\dot \theta_2$ (rad/s)',
           title=r'$\dot \theta_2$ com passo h =' + str(h))
    ax.grid()

    ax = plt.subplot(325)
    plt. plot(np.arange(0.0, tempo, h), a[:, 0])
    ax.set(xlabel='tempo (s)', ylabel=r'$\ddot \theta_1$ (rad/s²)',
           title=r'$\ddot \theta_1$ com passo h =' + str(h))
    ax.grid()

    ax = plt.subplot(326)
    plt. plot(np.arange(0.0, tempo, h), a[:, 1])
    ax.set(xlabel='tempo (s)', ylabel=r'$\ddot \theta_2$ (rad/s²)',
           title=r'$\ddot \theta_2$ com passo h =' + str(h))
    ax.grid()

    mng = plt.get_current_fig_manager()
    mng.full_screen_toggle()
    plt.tight_layout()
    plt.show()


######
# MAIN

# Metodo de euler
for h in H:
    e, a = euler(condicoes_iniciais, fa, fb, tempo, h)
    graficos(e, a, h, "Método de Euler")

# Metodo de euler modificado (RK2)
for h in H:
    e, a = RK2(condicoes_iniciais, fa, fb, tempo, h)
    graficos(e, a, h, "Método de Euler modificado (RK2)")

# Metodo de Runge-Kutta 4 (RK4)
for h in H:
    e, a = RK4(condicoes_iniciais, fa, fb, tempo, h)
    graficos(e, a, h, "Método de Runge-Kutta 4 (RK4)")

# scipy.io.savemat('estados-EP1-RK4-0_001.mat', dict(
#     theta=[e[:, 0], e[:, 2]], t=np.linspace(0, tempo, tempo/h)))