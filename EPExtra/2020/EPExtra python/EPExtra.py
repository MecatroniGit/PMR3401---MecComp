import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

plt.close('all')

# Declarando as constantes
r = 26                 #metros
Vmax = 500             #kV
Imax = 200             #Ampere
fi = 0                 #rad
frequencia = 60        #Hz
rc = 0.02              #metros
miar = 1.2566e-6       #H/m
misolo = 2 * 1.2567e-6 #H/m
sigmaar = 1e-10        #S/m
sigmasolo = 1e-2       #S/m

tamanho = 0.5            #tamanho do elemento triangular

Afio = -(miar*Imax)/(2*math.pi*rc)
Ae = np.power(tamanho, 2)/2 # Área do elemento


# Discretização do domínio

# Inicializando a matriz dos elementos e nós.
# A matriz dos nos contem as coordenadas de cada.
# A matriz dos elementos contém os índices dos nós
nos = [[0, 0]]          # [x, y, V, A]
elementos = [[0, 0, 0]] # [no1, no2, no3]

# Malhar apenas metade da figura por conta da simetria
# [x1, y1, x2, y2, x3, y3] -> nos locais
# Numero da linha = numero do elemento

#origem no centro da figura
x = 0
y = 0

while x < r:
    while y < np.sqrt(r**2 - x**2):
        # Pontas dos quadrados
        p1 = [x, y]
        p2 = [x+tamanho, y]
        p3 = [x, y+tamanho]
        p4 = [x+tamanho, y+tamanho]
        p5 = [x, -y]
        p6 = [x+tamanho, -y]
        p7 = [x, -y-tamanho]
        p8 = [x+tamanho, -y-tamanho]

        # Incluindo os pontos na lista de nós e pegando seus índices
        if not p1 in nos: nos.append(p1)
        i1 = nos.index(p1)
        if not p2 in nos: nos.append(p2)
        i2 = nos.index(p2)
        if not p3 in nos: nos.append(p3)
        i3 = nos.index(p3)
        if not p4 in nos: nos.append(p4)
        i4 = nos.index(p4)
        if not p5 in nos: nos.append(p5)
        i5 = nos.index(p5)
        if not p6 in nos: nos.append(p6)
        i6 = nos.index(p6)
        if not p7 in nos: nos.append(p7)
        i7 = nos.index(p7)
        if not p8 in nos: nos.append(p8)
        i8 = nos.index(p8)

        #triangulo de baixo acima do eixo x
        elementos.append([i1, i2, i4])

        #triangulo de cima acima do eixo x
        elementos.append([i1, i4, i3])

        #triangulo de cima abaixo do eixo x
        elementos.append([i5, i8, i6])

        #triangulo de baixo abaixo do eixo x
        elementos.append([i5, i7, i8])

        y += tamanho

    y = 0
    x += tamanho


nos = np.array(nos)
elementos = np.array(elementos[1:])

# Criando colunas para os potenciais
nos = np.insert(nos, 2, values=-1, axis=1) # Coluna do potencial elétrico
nos = np.insert(nos, 3, values=-1, axis=1) # Coluna do potencial magnético

# Aplicando condições de contorno
for n in nos:
    raio_no = np.sqrt(n[0]**2 + n[1]**2)
    # Caso o nó esteja na borda do domínio:
    if (raio_no > r - 1):
        n[2] = 0 # Potencial elétrico nulo
        n[3] = 0 # Potencial magnético nulo

    # Caso o nó esteja no carro
    if (n[0] <= 1.5 and n[1] <= 2 and n[1] >= 0):
        n[2] = 0 # Potencial elétrico nulo
    
    # Caso o nó esteja em um dos fios:
    if ((n[0] == 4 and n[1] == 14) or (n[0] == 6 and n[1] == 10)):
        n[2] = Vmax # Potencial elétrico máximo
        n[3] = Afio # Campo magnético


# Construção da matriz global KV
KV = np.zeros(shape = (len(nos), len(nos)))
for e in elementos:
    # Cada elemento no formato [no1, no2, no3]
    b1 = nos[e[1]][1] - nos[e[2]][1]
    b2 = nos[e[2]][1] - nos[e[0]][1]
    b3 = nos[e[0]][1] - nos[e[1]][1]

    c1 = nos[e[2]][0] - nos[e[1]][0]
    c2 = nos[e[0]][0] - nos[e[2]][0]
    c3 = nos[e[1]][0] - nos[e[0]][0]

    if (nos[e[0]][1] + nos[e[1]][1] + nos[e[2]][1]) >= 0: sigma = sigmaar
    else: sigma = sigmasolo

    KV[e[0]][e[0]] += sigma * (b1*b1 + c1*c1)/(4*Ae)
    KV[e[0]][e[1]] += sigma * (b1*b2 + c1*c2)/(4*Ae)
    KV[e[0]][e[2]] += sigma * (b1*b3 + c1*c3)/(4*Ae)

    KV[e[1]][e[0]] += sigma * (b1*b2 + c1*c2)/(4*Ae)
    KV[e[1]][e[1]] += sigma * (b2*b2 + c2*c2)/(4*Ae)
    KV[e[1]][e[2]] += sigma * (b2*b3 + c2*c3)/(4*Ae)

    KV[e[2]][e[0]] += sigma * (b1*b3 + c1*c3)/(4*Ae)
    KV[e[2]][e[1]] += sigma * (b2*b3 + c2*c3)/(4*Ae)
    KV[e[2]][e[2]] += sigma * (b3*b3 + c3*c3)/(4*Ae)


# Matriz de carregamento e mudando matriz K inserindo Vmax
carregamento_V = np.zeros(len(nos))
# V = Vmax nos pontos (4, 10) e (6, 14)
for i in range(len(nos)):
    # Caso o nó tenha um carregamento
    if (nos[i][2] != -1): 
        for j in range(len(KV)):
            carregamento_V[j] -= KV[j][i] * nos[i][2]
            KV[j][i] = 0 #colunas = 0
            KV[i][j] = 0 #linhas = 0

        KV[i][i] = 1
        carregamento_V[i] = nos[i][2]

# Solucionando o sistema:
V = np.linalg.solve(KV, carregamento_V)
for i in range(len(nos)):
    nos[i][2] = V[i]

# Construção da matriz global KA
KA = np.zeros(shape = (len(nos), len(nos)))
for e in elementos:
    # Cada elemento no formato [no1, no2, no3]
    b1 = nos[e[1]][1] - nos[e[2]][1]
    b2 = nos[e[2]][1] - nos[e[0]][1]
    b3 = nos[e[0]][1] - nos[e[1]][1]

    c1 = nos[e[2]][0] - nos[e[1]][0]
    c2 = nos[e[0]][0] - nos[e[2]][0]
    c3 = nos[e[1]][0] - nos[e[0]][0]

    if (nos[e[0]][1] + nos[e[1]][1] + nos[e[2]][1]) >= 0: mi = miar
    else: mi = misolo

    KA[e[0]][e[0]] += mi * (b1*b1 + c1*c1)/(4*Ae)
    KA[e[0]][e[1]] += mi * (b1*b2 + c1*c2)/(4*Ae)
    KA[e[0]][e[2]] += mi * (b1*b3 + c1*c3)/(4*Ae)

    KA[e[1]][e[0]] += mi * (b1*b2 + c1*c2)/(4*Ae)
    KA[e[1]][e[1]] += mi * (b2*b2 + c2*c2)/(4*Ae)
    KA[e[1]][e[2]] += mi * (b2*b3 + c2*c3)/(4*Ae)

    KA[e[2]][e[0]] += mi * (b1*b3 + c1*c3)/(4*Ae)
    KA[e[2]][e[1]] += mi * (b2*b3 + c2*c3)/(4*Ae)
    KA[e[2]][e[2]] += mi * (b3*b3 + c3*c3)/(4*Ae)

# Construindo a matriz de carregamentos A
carregamento_A = np.zeros(len(nos))
for i in range(len(nos)):
    # Caso o nó pertença ao fio
    if (nos[i][0] == 4 and nos[i][1] == 14) or (nos[i][0] == 6 and nos[i][1] == 10):
        carregamento_A[i] = Afio

    # Caso o nó pertença à borda
    if (nos[i][3] == 0):
        for j in range(len(KA)):
            carregamento_A[j] -= KA[j][i] * nos[i][3]
            KA[i][j] = 0
            KA[j][i] = 0
        KA[i][i] = 1
        carregamento_A[i] = nos[i][3]

# Solucionando o sistema
A = np.linalg.solve(KA, carregamento_A)
for i in range(len(nos)):
    nos[i][3] = A[i]


# Criando as colunas do elemento que serão usadas para armazenar
#   a densidade de fluxo Magnético e elétrico
elementos = np.insert(elementos, 3, values=0, axis=1) # Coluna do fluxo elétrico x
elementos = np.insert(elementos, 4, values=0, axis=1) # Coluna do fluxo elétrico y
elementos = np.insert(elementos, 4, values=0, axis=1) # Coluna do fluxo magnético x
elementos = np.insert(elementos, 5, values=0, axis=1) # Coluna do fluxo magnético y
elementos = np.insert(elementos, 6, values=0, axis=1) # Coluna do campo magnético y
elementos = np.insert(elementos, 7, values=0, axis=1) # Coluna do campo magnético y


# Calculando os vetores gradientes de V e A
for e in elementos:
    b1 = nos[e[1]][1] - nos[e[2]][1]
    b2 = nos[e[2]][1] - nos[e[0]][1]
    b3 = nos[e[0]][1] - nos[e[1]][1]

    c1 = nos[e[2]][0] - nos[e[1]][0]
    c2 = nos[e[0]][0] - nos[e[2]][0]
    c3 = nos[e[1]][0] - nos[e[0]][0]

    # Fluxo elétrico E
    e[3] = -(b1*nos[e[0]][2] + b2*nos[e[1]][2] + b3*nos[e[2]][2])/(2*Ae)
    e[4] = -(c1*nos[e[0]][2] + c2*nos[e[1]][2] + c3*nos[e[2]][2])/(2*Ae)
    
    # Fluxo magnético B
    e[5] = (c1*nos[e[0]][3] + c2*nos[e[1]][3] + c3*nos[e[2]][3])/(2*Ae)
    e[6] = -(b1*nos[e[0]][3] + b2*nos[e[1]][3] + b3*nos[e[2]][3])/(2*Ae)
    
    # Campo magnético H
    if (nos[e[0]][0] + nos[e[1]][0] + nos[e[2]][0]) >= 0:
        e[7] = e[5]/miar
        e[8] = e[6]/miar
    else:
        e[7] = e[5]/misolo
        e[8] = e[6]/misolo



# PLOTS

# Plotagem da discretização
# Duplicando e espelhando o sistema no eixo y:
# elementos = np.append(elementos, [[e[0]+len(nos), e[1]+len(nos), e[2]+len(nos), e[3], e[4], e[5], e[6], e[7], e[8]] for e in elementos], axis=0)
nos2 = np.append(nos, [[-n[0], n[1], 0, 0] for n in nos], axis=0)
fig = plt.figure(1)
line_collection1 = mc.LineCollection([[nos2[e[0]][:2], nos2[e[1]][:2]] for e in elementos], linewidths=0.6)
line_collection2 = mc.LineCollection([[nos2[e[0]][:2], nos2[e[2]][:2]] for e in elementos], linewidths=0.6)
line_collection3 = mc.LineCollection([[nos2[e[1]][:2], nos2[e[2]][:2]] for e in elementos], linewidths=0.6)
ax = fig.add_subplot(111, title="Discretização do Sistema")
ax.axis('equal')
ax.add_collection(line_collection1)
ax.add_collection(line_collection2)
ax.add_collection(line_collection3)
ax.scatter(nos2[:, 0], nos2[:, 1], alpha=0.5, s=3)

# Plotando V
fig = plt.figure(2)
X, Y = np.meshgrid(np.arange(-r-1, r+1, tamanho), np.arange(-r-1, r+1, tamanho))
Z = np.zeros(np.shape(X))
for i in range(len(nos)):
    Z[int((nos[i][1]+r)/tamanho)][int((nos[i][0]+r)/tamanho)] = V[i]
    Z[int((nos[i][1]+r)/tamanho)][-int((nos[i][0]+r)/tamanho)] = V[i]
ax1 = fig.add_subplot(121, title='Potencial Elétrico V (kV)')
ax1.axis('equal')
CS = ax1.contour(X, Y, Z, levels=15)
ax1.clabel(CS, inline=1, fmt='%1.0f', fontsize=6)
ax1.set_xlabel("x (m)")
ax1.set_ylabel("y (m)")
ax2 = fig.add_subplot(122, projection='3d', title='Potencial Elétrico V (kV)')
ax2.plot_surface(X, Y, Z, linewidth=5, cmap='viridis', antialiased=True)
ax2.set_xlabel("x (m)")
ax2.set_ylabel("y (m)")
ax2.set_zlabel("V (kV)")

# Plotando A
fig = plt.figure(3)
Z = np.zeros(np.shape(X))
for i in range(len(nos)):
    Z[int((nos[i][1]+r)/tamanho)][int((nos[i][0]+r)/tamanho)] = A[i]
    Z[int((nos[i][1]+r)/tamanho)][-int((nos[i][0]+r)/tamanho)] = A[i]
ax1 = fig.add_subplot(121, title='Potencial Vetor Magnético A (Wb/m)')
ax1.axis('equal')
# cax = ax.matshow(Z)
# fig.colorbar(cax)
CS = ax1.contour(X, Y, Z, levels=15)
ax1.clabel(CS, inline=1, fmt='%1.0f', fontsize=6)
ax1.set_xlabel("x (m)")
ax1.set_ylabel("y (m)")
ax2 = fig.add_subplot(122, projection='3d', title='Potencial Magnético A (Wb/m)')
ax2.plot_surface(X, Y, Z, linewidth=5, cmap='viridis', antialiased=True)
ax2.set_xlabel("x (m)")
ax2.set_ylabel("y (m)")
ax2.set_zlabel("V (kV)")

# Plotando o fluxo Elétrico E
fig = plt.figure(4)
ax1 = fig.add_subplot(111, title='Fluxo elétrico E')
ax1.axis('equal')
X = [(nos[e[0]][0] + nos[e[1]][0] + nos[e[2]][0])/3 for e in elementos]
Y = [(nos[e[0]][1] + nos[e[1]][1] + nos[e[2]][1])/3 for e in elementos]
X = np.append(X, np.negative(X))
Y = np.append(Y, Y)
dX = np.append(elementos[:, 3], np.negative(elementos[:, 3]))
dY = np.append(elementos[:, 4], elementos[:, 4])
M = np.sqrt(np.power(dX, 2) + np.power(dY, 2))
Q = ax1.quiver(X[::3], Y[::3], dX[::3], dY[::3], M[::3], units='x', linewidth=0.2, scale=15, cmap=cm.viridis)
qk = ax1.quiverkey(Q, 0.9, 0.9, 15, r'$15 kV$', labelpos='E', coordinates='figure')

# Plotando o fluxo Magnético B
fig = plt.figure(5)
ax2 = fig.add_subplot(111, title='Fluxo magnético B')
ax2.axis('equal')
dX = np.append(elementos[:, 5], np.negative(elementos[:, 5]))
dY = np.append(elementos[:, 6], elementos[:, 6])
M = np.sqrt(np.power(dX, 2) + np.power(dY, 2))
Q = ax2.quiver(X[::3], Y[::3], dX[::3], dY[::3], M[::3], units='x', linewidth=0.2, scale=100, cmap=cm.viridis)
qk = ax2.quiverkey(Q, 0.9, 0.9, 100, r'$100 T$', labelpos='E', coordinates='figure')

# Plotando o Campo Magnético H
fig = plt.figure(6)
ax2 = fig.add_subplot(111, title='Campo magnético H')
ax2.axis('equal')
dX = np.append(elementos[:, 7], np.negative(elementos[:, 7]))
dY = np.append(elementos[:, 8], elementos[:, 8])
M = np.sqrt(np.power(dX, 2) + np.power(dY, 2))
Q = ax2.quiver(X[::3], Y[::3], dX[::3], dY[::3], M[::3], units='x', linewidth=0.2, scale=10e7, cmap=cm.viridis)
qk = ax2.quiverkey(Q, 0.9, 0.9, 10e7, r'$10^7 H$', labelpos='E', coordinates='figure')

plt.show()