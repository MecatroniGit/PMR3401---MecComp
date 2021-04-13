import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import sys

V = 100 / 3.6 #Velocidade do vento em m/s
h = 3 #metros
d = 5 * h
L = 2 * h
H = 8 * h
comprimento_estrutura = 60

#Parâmetros de simulação
dx = dy = h/8
epsilon = 0.01

linhas = int(H//dy) + 1
colunas = int((2*d+L) //dx) + 1

#propriedades do ar
ro = 1.25 #kg/m^3
gama = 1.4
k = 0.026
cp = 1002
K = (ro * cp) / k

matriz = np.zeros(shape = (linhas, colunas))

def main():
  erro_max = 1
  contador = 0

  while erro_max > epsilon:
    erro_max = 0
    contador += 1

    for i in range(len(matriz) - 1):
      for j in range(len(matriz[0])):

        altura = 24 - i * dy
        comprimento = j * dx
        modulo = np.sqrt(np.power(18 - j * dx, 2) + np.power(altura - 3, 2))

        #considera corrente dentro do galpão igual a zero
        if comprimento >= 15 and comprimento <= 21 and altura <= 3:
          psi = 0

        elif comprimento >= 15 and comprimento <= 21 and altura >= 3  and altura <= 6 and modulo <= 3:
          psi = 0

        elif altura <= 3 and comprimento + dx > 15 and comprimento - dx < 21:

          #condicao da borda esquerda do galpão
          if 15 - comprimento < dx and 15 - comprimento > 0:
            b = round(15 % dx, 1)

            psi = (b / (b + 1)) * ((matriz[i][j-1] / (b + 1)) + ((matriz[i-1][j] + matriz[i+1][j]) / 2))

          #condicao da borda direita do galpao
          elif comprimento - 21 < dx and comprimento - 21 > 0:
            b = round(21 % dx, 1)
            
            psi = (b / (b + 1)) * ((matriz[i][j+1] / (b + 1)) + ((matriz[i-1][j] + matriz[i+1][j]) / 2))
      
        #condicao para borda do telhado
        elif comprimento >= 15 and comprimento <= 21 and altura < 6 + dy and altura >= 3:
          if comprimento < 18 and np.sqrt(np.power(18 - (j + 1) * dx, 2) + np.power(21 - i * dy, 2)) < 3 and np.sqrt(np.power(18 - j * dx, 2) + np.power(21 - (i + 1) * dy, 2)) < 3:
            #condicao para borda do teto esquerdo do galpão com uso de a e b    
            a = (21 - i * dy - np.sqrt(9 - np.power(j * dx - 18, 2))) / dy
            b = (abs(18 - j * dx) - 3 * np.cos(np.arcsin((21 - i * dy) / 3))) / dx
            
            psi = ((a * b) / (a + b)) * ((matriz[i][j-1] / (b + 1)) + (matriz[i-1][j] / (a + 1)))

          elif comprimento > 18 and np.sqrt(np.power(18 - (j - 1) * dx, 2) + np.power(21 - i * dy, 2)) < 3 and np.sqrt(np.power(18 - j * dx, 2) + np.power(21 - (i + 1) * dy, 2)) < 3:
            #condicao para a borda do teto direito do galpao com uso de a e b
            a = (21 - i * dy - np.sqrt(9 - np.power(j * dx - 18, 2))) / dy
            b = (abs(18 - j * dx) - 3 * np.cos(np.arcsin((21 - i * dy) / 3))) / dx
            
            psi = ((a * b) / (a + b)) * ((matriz[i][j+1] / (b + 1)) + (matriz[i-1][j] / (a + 1)))
            
          elif comprimento < 18 and np.sqrt(np.power(18 - (j + 1) * dx, 2) + np.power(21 - i * dy, 2)) < 3: 
            #condicao para a borda do teto esquerdo do galpão com uso de b
            b = abs((abs(18 - j * dx) - 3 * np.cos(np.arcsin((21 - i * dy) / 3))) / dx)

            psi = b * ((matriz[i-1][j] / (2 * (b + 1))) + (matriz[i+1][j] / (2 * (b + 1))) +
             (matriz[i][j-1] / ((b + 1) * (b + 1))))

          elif comprimento > 18 and np.sqrt(np.power(18 - (j - 1) * dx, 2) + np.power(21 - i * dy, 2)) < 3:
            #condicao para a borda do teto direito do galpão com uso de b
            b = (abs(18 - j * dx) - 3 * np.cos(np.arcsin((21 - i * dy) / 3))) / dx

            psi = b * ((matriz[i-1][j] / (2 * (b + 1))) + (matriz[i+1][j] / (2 * (b + 1))) + (matriz[i][j+1] / ((b + 1) * (b + 1))))

          elif np.sqrt(np.power(18 - j * dx, 2) + np.power(21 - (i + 1) * dy, 2)) < 3:
            #condicao para a borda do teto na parte superior com uso de a
            a = abs((21 - i * dy - np.sqrt(9 - np.power(j * dx - 18, 2))) / dy)

            psi = a * ((matriz[i][j+1] / (2 * a + 2)) + (matriz[i][j-1] / (2 * a + 2)) + (matriz[i-1][j] / ((a + 1) * (a + 1))))

          #atualiza os pontos dentro da condicao inicial, mas que não se encaixam em nenhum dos casos anteriores
          else: 
            psi = (matriz[i][j+1] + matriz[i][j-1] + matriz[i-1][j] + matriz[i+1][j]) / 4

        elif (i == 0 or j == 0) and j != (len(matriz[0]) - 1):
          #canto superior esquerdo
          if i == 0 and j == 0:
            psi = (2 * matriz[i][j+1] + 2 * matriz[i+1][j] + 2 * V * dy) / 4
          
          #limite esquerdo
          elif j == 0:
            psi = (2 * matriz[i][j+1] + matriz[i-1][j] + matriz[i+1][j]) / 4

          #limite superior sem contar canto superior esquerdo
          elif i == 0 and j != (len(matriz[0]) - 1):
            psi = (matriz[i][j+1] + matriz[i][j-1] + 2 * matriz[i+1][j] + 2 * V * dy) / 4
          

        elif i ==0 or j == (len(matriz[0]) - 1):
          #canto superior direito
          if i == 0 and j == (len(matriz[0]) - 1):
            psi = (2 * matriz[i][j-1] + 2 * matriz[i+1][j] + 2 * V * dy) / 4

          #limite direito
          elif j == (len(matriz[0]) - 1):
            psi = (2 * matriz[i][j-1] + matriz[i-1][j] + matriz[i+1][j]) / 4

        #para pontos i, j generalizados internos
        else:
          psi = (matriz[i][j+1] + matriz[i][j-1] + matriz[i-1][j] + matriz[i+1][j]) / 4

        # Sobrerrelaxação:
        fator = 1.85
        psi_novo = fator*psi + (1-fator)*matriz[i][j]
        erro_atual = abs((psi_novo - matriz[i][j])/psi_novo) if psi_novo != 0 else 0
        erro_max = max(erro_atual, erro_max) #atualiza com maior erro   
        matriz[i][j] = psi_novo

  print("Passo:", dx)
  print("Epsilon:", epsilon)
  print("Iterações:", contador)
  print("Média de psi da borda superior:", np.mean(matriz[0]))

  parte_1a(matriz)
  velocidadex, velocidadey = parte_1b(matriz)
  # dpressaodominio = parte_1c(velocidadex, velocidadey)
  # pressaoteto = parte_1d(dpressaodominio)
  # parte_1e(pressaoteto)
  malha = parte_2a(velocidadex, velocidadey)
  parte_2b(malha)


def parte_1a(matriz):
  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_title("Função de corrente do escoamento 2D")
  cax = ax.matshow(matriz, cmap='inferno')
  fig.colorbar(cax)
  plt.show()

  nx = int((2 * d + L)/dx) + 1
  ny = int((H)/dy) + 1
  x = np.linspace(0, 2*d+L, nx)
  y = np.linspace(0, H, ny)
  X, Y = np.meshgrid(x, y)
  fig1 = plt.figure()
  ax1 = fig1.gca(projection='3d')
  ax1.set_title("Função de corrente do escoamento 3D")
  s = ax1.plot_surface(X, Y, matriz, cmap=cm.inferno, linewidth=5, antialiased=True)
  ax1.set_xlabel("x (m)")
  ax1.set_ylabel("y (m)")
  ax1.set_zlabel("Psi (m²/s)")
  plt.show()

def parte_1b(matriz):
  velocidadex = np.zeros(shape = (linhas, colunas)) #matriz de velocidades na direção x
  velocidadey = np.zeros(shape = (linhas, colunas)) #matriz de velocidades na direção y

  #geração da matriz de velocidades em y, vetor v
  for i in range(len(matriz)): #percorre todas as linhas da matriz
    for j in range(len(matriz[0])): #percorre todas as colunas da matriz
      altura = 24 - i * dy
      comprimento = j * dx
      modulo = np.sqrt(np.power(18 - j * dx, 2) + np.power(altura - 3, 2))
      

      if comprimento >= 15 and comprimento <= 21 and altura <= 3:#considera corrente dentro do galpão igual a zero
        velocidadey[i][j] = 0
      
      elif comprimento >= 15 and comprimento <= 21 and altura >= 3  and altura <= 6 and modulo <= 3:
        velocidadey[i][j] = 0
      
      elif j == 0: #condição do limite esquerdo
        velocidadey[i][j] = 0
          
      elif j == len(matriz[0]) - 1: #condição do limite direito
        velocidadey[i][j] = 0
          
      elif 15 - comprimento < dx and 15 - comprimento > 0 and altura <= 3: #condição da borda do galpão esquerdo e suas proximidades
        velocidadey[i][j] = (-3 * matriz[i][j] + 4 * matriz[i][j-1] - matriz[i][j-2]) / (2 * dx)
          
      elif comprimento - 21 < dx and comprimento - 21 > 0 and altura <= 3: #condição da borda do galpão direito e suas proximidades
        velocidadey[i][j] = (matriz[i][j+2] - 4 * matriz[i][j+1] + 3 * matriz[i][j]) / (2 * dx)

      elif comprimento < 18 and comprimento >= 15 and np.sqrt(np.power(18 - (j + 1) * dx, 2) + np.power(21 - i * dy, 2)) < 3 and altura >= 3:
      #condição telhado esquerdo
        b = (abs(18 - j * dx) - 3 * np.cos(np.arcsin((21 - i * dy) / 3))) / dx

        velocidadey[i][j] = (-1 / (b * (b + 1) * dx)) * (-b * b * matriz[i][j-1] + (b * b - 1) * matriz[i][j])

      elif comprimento > 18 and comprimento <= 21 and np.sqrt(np.power(18 - (j - 1) * dx, 2) + np.power(21 - i * dy, 2)) < 3 and altura >= 3:
      #condição do telhado direito
        b = (abs(18 - j * dx) - 3 * np.cos(np.arcsin((21 - i * dy) / 3))) / dx

        velocidadey[i][j] = (-1 / (b * (b + 1)* dx)) * (b * b * matriz[i][j+1] + (1 - b * b) * matriz[i][j])
          
      else: #pontos internos generalizados
        velocidadey[i][j] = (matriz[i][j-1] - matriz[i][j+1]) / (2 * dx)

  #geração da matriz de velocidade em x, vetor u
  for i in range(len(matriz)):
    for j in range(len(matriz[0])):
      altura = 24 - i * dy
      comprimento = j * dx
      modulo = np.sqrt(np.power(18 - j * dx, 2) + np.power(altura - 3, 2))

      if comprimento >= 15 and comprimento <= 21 and altura < 3:#considera corrente dentro do galpão igual a zero
        velocidadex[i][j] = 0
          
      elif comprimento >= 15 and comprimento <= 21 and altura >= 3  and altura <= 6 and modulo <= 3:
        velocidadex[i][j] = 0
      
      elif i == 0: #condição do limite superior
        velocidadex[i][j] = V

      elif i == len(matriz) - 1: #condição do limite inferior
        velocidadex[i][j] = 0

      elif comprimento >= 15 and comprimento <= 21 and np.sqrt(np.power(21 - (i + 1) * dy, 2) + np.power(18 - j * dx, 2)) < 3:
        #condição da borda do galpão superior e suas proximidades
        a = (21 - i * dy - np.sqrt(9 - np.power(j * dx - 18, 2))) / dy

        velocidadex[i][j] = (1/(a * (a + 1) * dy)) * ((1 - a * a) * matriz[i][j] + a * a * matriz[i-1][j])

      else: #pontos internos generalizados
        velocidadex[i][j] = (matriz[i-1][j] - matriz[i+1][j]) / (2 * dy)


  max_vel = np.amax(np.sqrt(np.power(velocidadex, 2) + np.power(velocidadey, 2)))
  print("Velocidade máxima:", max_vel)
  x = np.linspace(0, 36, (36/dx) + 1.0)
  y = np.linspace(24, 0, (24/dx) + 1.0)

  fig3, ax3 = plt.subplots()
  ax3.set_title("Vetores de velocidade absoluta do escoamento")
  M = np.hypot(velocidadex, velocidadey)
  Q = ax3.quiver(x, y, velocidadex, velocidadey, M, units='x', width=0.08, scale=20, cmap=cm.inferno)
  qk = ax3.quiverkey(Q, 0.9, 0.9, 20, r'$20 \frac{m}{s}$', labelpos='E', coordinates='figure')

  plt.show()

  return velocidadex, velocidadey

def parte_1c(velocidadex, velocidadey):
  x = np.linspace(0, 36, (36/dx) + 1.0)
  y = np.linspace(24, 0, (24/dx) + 1.0)

  dpressaodominio = np.zeros(shape = (linhas, colunas)) #matriz da vriação de pressão no domínio

  for i in range(len(matriz)):
    for j in range(len(matriz[0])):          
      dpressaodominio[i][j] = ((ro * (gama - 1)) / (2 * gama)) * (-np.power(velocidadex[i][j], 2) - np.power(velocidadey[i][j], 2))

  min_pres = np.amin(dpressaodominio)
  print("Pressão mínima:", min_pres)

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_title("Variação de pressão no domínio")
  cax = ax.matshow(dpressaodominio, interpolation='nearest', cmap=cm.inferno)
  fig.colorbar(cax)
  plt.show()
  
  nx = int((2 * d + L)/dx) + 1
  ny = int((H)/dy) + 1
  x = np.linspace(2*d+L, 0, nx)
  y = np.linspace(H, 0, ny)
  X, Y = np.meshgrid(x, y)
  ax1 = fig.gca(projection="3d")
  ax1.set_title("Variação de pressão p(x, y) - p_atm no domínio")
  s = ax1.plot_surface(X, Y, dpressaodominio, cmap=cm.inferno, linewidth=5, antialiased=True)
  ax1.set_xlabel("x (m)")
  ax1.set_ylabel("y (m)")
  ax1.set_zlabel("Pressão (pascal)")
  plt.show()

  return dpressaodominio

def parte_1d(dpressaodominio):
  l = []
  c = []
  for j in range(round(15//dx), round(21//dx) + 1):
    for i in range(round(18//dy - 1), round(21//dy) + 1):
      altura = 24 - i * dy
      
      if altura - np.sqrt(9 - np.power(j * dx - 18, 2)) - 3 <= dy and altura - np.sqrt(9 - np.power(j * dx - 18, 2)) - 3 > 0:
        l.append(i)
        c.append(j)

  pressaoteto = []
  for i in range(len(l)):      
    pressaoteto.append(dpressaodominio[l[i]][c[i]])

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_title("Variação da pressão na borda do telhado do galpão")
  ax.plot(pressaoteto)
  ax.set_xlabel("x (m)")
  ax.set_ylabel("Pressão (N/m²)")
  plt.show()

  return pressaoteto

def parte_1e(pressaoteto):
  comprimento = 60 #comprimento do telhado
  forca_total = 0 #armazena a somatória das pressões sobre o telhado

  for i in range(len(pressaoteto)):
    forca_total = forca_total + pressaoteto[i] *dx * comprimento

  print("Força total:", forca_total)

def parte_2a(velocidadex, velocidadey):
  malha = np.zeros(shape=(linhas, colunas))
  contador = 0
  erro_max = 1

  while erro_max > epsilon:
    contador += 1
    erro_max = 0
    
    for i in range(len(malha)): #percorre as linhas da matriz
      for j in range(len(malha[0])): #percorre as colunas da matriz
        altura = 24 - i * dy
        comprimento = j * dx
        modulo = np.sqrt(np.power(18 - j * dx, 2) + np.power(altura - 3, 2))

        #configura temperatura dentro do galpão             
        if comprimento >= 15 and comprimento <= 21 and altura <= 3:#considera corrente dentro do galpão igual a zero
          temperatura = 40
    
        elif comprimento >= 15 and comprimento <= 21 and altura >= 3  and altura <= 6 and modulo <= 3:
          temperatura = 40

        #condicao do contorno do sistema
        elif i == 0 or j == 0 or i == len(malha) - 1 or j == len(malha[0]) - 1:

          #condicao do canto esquerdo superior
          if i==0 and j==0:
            temperatura = 20
              
          #condicao do canto direito superior
          elif i == 0 and j == len(malha[0]) - 1:
            temperatura = (malha[i][j-1] + malha[i+1][j]) / 2

          #condicao do canto esquerdo inferior
          elif i == len(malha) - 1 and j == 0:
            temperatura = 20

          #condicao do canto direito inferior
          elif i == len(malha) - 1 and j == len(malha[0]) - 1:
            temperatura = (malha[i][j-1] + malha[i-1][j]) / 2

          #condicao do limite esquerdo
          elif j == 0:
            temperatura = 20

          #condicao do limite superior
          elif i == 0:
            temperatura = (malha[i][j+1] + malha[i][j-1] + 2 * malha[i+1][j] + K * velocidadex[i][j] * dx * malha[i][j-1]) / (K * dx * velocidadex[i][j] + 4)

          #condicao do limite direito
          elif j == len(malha[0]) - 1:
            if velocidadey[i][j] >= 0:
              temperatura = (2 * malha[i][j-1] + malha[i-1][j] + malha[i+1][j] + K * dx * velocidadey[i][j] * malha[i+1][j]) / (K * dx * velocidadey[i][j] + 4)

            elif velocidadey[i][j] < 0:
              temperatura = (2 * malha[i][j-1] + malha[i-1][j] + malha[i+1][j] - K * dx * velocidadey[i][j] * malha[i+1][j]) / (-K * dx * velocidadey[i][j] + 4)

          #condicao do limite inferior
          elif i == len(malha) - 1:
            temperatura = (malha[i][j+1] + malha[i][j-1] + 2 * malha[i-1][j] + K * velocidadex[i][j] * dx * malha[i][j-1]) / (K * dx * velocidadex[i][j] + 4)

        elif altura <= 3 and comprimento + dx > 15 and comprimento - dx < 21:
          #condicao da borda esquerda do galpao
          if 15 - comprimento < dx and 15 - comprimento > 0:
            b = round(15 % dx, 1)

            if velocidadey[i][j] >= 0:
              temperatura = ((80/(b*(b+1))) + (2*malha[i][j-1]/(b+1)) + malha[i-1][j] + malha[i+1][j] + K*dx*(velocidadex[i][j]*malha[i][j-1] + velocidadey[i][j]*malha[i+1][j])) / ((2/b) + 2 + K*dx*(velocidadex[i][j] + velocidadey[i][j]))
                              
            elif velocidadey[i][j] < 0:
              temperatura = ((80/(b*(b+1))) + (2*malha[i][j-1]/(b+1)) + malha[i-1][j] + malha[i+1][j] + K*dx*(velocidadex[i][j]*malha[i][j-1] - velocidadey[i][j]*malha[i-1][j])) / ((2/b) + 2 + K*dx*(velocidadex[i][j] - velocidadey[i][j]))

          #condição da borda direita do galpao
          elif comprimento - 21 < dx and comprimento - 21 > 0:
            b = round(21 % dx, 1)

            if velocidadey[i][j] >= 0:
              temperatura = ((80/(b*(b+1))) + (2*malha[i][j+1]/(b+1)) + malha[i-1][j] + malha[i+1][j] + K*dx*(velocidadex[i][j]*malha[i][j-1] + velocidadey[i][j]*malha[i+1][j])) / ((2/b) + 2 + K*dx*(velocidadex[i][j] + velocidadey[i][j]))

            elif velocidadey[i][j] < 0:
              temperatura = ((80/(b*(b+1))) + (2*malha[i][j+1]/(b+1)) + malha[i-1][j] + malha[i+1][j] + K*dx*(velocidadex[i][j]*malha[i][j-1] - vrlocidadey[i][j]*malha[i-1][j])) / ((2/b) + 2 + K*dx*(velocidadex[i][j] - velocidadey[i][j]))

        #condicao da borda do telhado
        elif comprimento >= 15 and comprimento <= 21 and altura < 6 + dy and altura >= 3:
          
          #condicao para borda do teto esquerda com uso de a e b
          if comprimento < 18 and np.sqrt(np.power(18 - (j + 1) * dx, 2) + np.power(21 - i * dy, 2)) < 3 and np.sqrt(np.power(18 - j * dx, 2) + np.power(21 - (i + 1) * dy, 2)) < 3:
            a = (21 - i * dy - np.sqrt(9 - np.power(j * dx - 18, 2))) / dy
            b = (abs(18 - j * dx) - 3 * np.cos(np.arcsin((21 - i * dy) / 3))) / dx

            if velocidadey[i][j] >= 0:
              temperatura = ((80/(b*(b+1))) + (80/(a*(a+1))) + (2*malha[i][j-1]/(b+1)) + (2*malha[i-1][j]/(a+1)) + K*dx*(velocidadex[i][j]*malha[i][j-1] + velocidadey[i][j]*malha[i+1][j])) / ((2/b) + (2/a) + K*dx*(velocidadex[i][j] + velocidadey[i][j]))

            elif velocidadey[i][j] < 0:
              temperatura = ((80/(b*(b+1))) + (80/(a*(a+1))) + (2*malha[i][j-1]/(b+1)) + (2*malha[i-1][j]/(a+1)) + K*dx*(velocidadex[i][j]*malha[i][j-1] - velocidadey[i][j]*malha[i-1][j])) / ((2/b) + (2/a) + K*dx*(velocidadex[i][j] - velocidadey[i][j]))

          #condicao para borda do teto direita com uso de a e b
          elif comprimento > 18 and np.sqrt(np.power(18 - (j - 1) * dx, 2) + np.power(21 - i * dy, 2)) < 3 and np.sqrt(np.power(18 - j * dx, 2) + np.power(21 - (i + 1) * dy, 2)) < 3:
            a = (21 - i * dy - np.sqrt(9 - np.power(j * dx - 18, 2))) / dy
            b = (abs(18 - j * dx) - 3 * np.cos(np.arcsin((21 - i * dy) / 3))) / dx

            if velocidadey[i][j] >= 0:
              temperatura = ((80/(b*(b+1))) + (80/(a*(a+1))) + (2*malha[i][j+1]/(b+1)) + (2*malha[i-1][j]/(a+1)) + K*dx*(velocidadex[i][j]*malha[i][j-1] + velocidadey[i][j]*malha[i+1][j])) / ((2/b) + (2/a) + K*dx*(velocidadex[i][j] + velocidadey[i][j]))

            elif velocidadey[i][j] < 0:
              temperatura = ((80/(b*(b+1))) + (80/(a*(a+1))) + (2*malha[i][j+1]/(b+1)) + (2*malha[i+1][j]/(a+1)) + K*dx*(velocidadex[i][j]*malha[i][j-1] - velocidadey[i][j]*malha[i-1][j])) / ((2/b) + (2/a) + K*dx*(velocidadex[i][j] - velocidadey[i][j]))


          #condicao para borda inferior do teto esquerda com uso de b
          elif comprimento < 18 and np.sqrt(np.power(18 - (j + 1) * dx, 2) + np.power(21 - i * dy, 2)) < 3:
            b = abs((abs(18 - j * dx) - 3 * np.cos(np.arcsin((21 - i * dy) / 3))) / dx)

            if velocidadey[i][j] >= 0:
              temperatura = ((80/(b*(b+1))) + (2*malha[i][j-1]/(b+1)) + malha[i-1][j] + malha[i+1][j] + K*dx*(velocidadex[i][j]*malha[i][j-1] + velocidadey[i][j]*malha[i+1][j])) / ((2/b) + 2 + K*dx*(velocidadex[i][j] + velocidadey[i][j]))

            elif velocidadey[i][j] < 0:
              temperatura = ((80/(b*(b+1))) + (2*malha[i][j-1]/(b+1)) + malha[i-1][j] + malha[i+1][j] + K*dx*(velocidadex[i][j]*malha[i][j-1] - velocidadey[i][j]*malha[i-1][j])) / ((2/b) + 2 + K*dx*(velocidadex[i][j] - velocidadey[i][j]))

          #condicao para borda inferior do teto direita com uso de b
          elif comprimento > 18 and np.sqrt(np.power(18 - (j - 1) * dx, 2) + np.power(21 - i * dy, 2)) < 3:
            b = abs((abs(18 - j * dx) - 3 * np.cos(np.arcsin((21 - i * dy) / 3))) / dx)

            if velocidadey[i][j] >= 0:
              temperatura = ((80/(b*(b+1))) + (2*malha[i][j+1]/(b+1)) + malha[i-1][j] + malha[i+1][j] + K*dx*(velocidadex[i][j]*malha[i][j-1] + velocidadey[i][j]*malha[i+1][j])) / ((2/b) + 2 + K*dx*(velocidadex[i][j] + velocidadey[i][j]))

            elif velocidadey[i][j] < 0:
              temperatura = ((80/(b*(b+1))) + (2*malha[i][j+1]/(b+1)) + malha[i-1][j] + malha[i+1][j] + K*dx*(velocidadex[i][j]*malha[i][j-1] - velocidadey[i][j]*malha[i-1][j])) / ((2/b) + 2 + K*dx*(velocidadex[i][j] - velocidadey[i][j]))

          #condicao para borda superior do teto com uso de a
          elif np.sqrt(np.power(18 - j * dx, 2) + np.power(21 - (i + 1) * dy, 2)) < 3:
            a = (21 - i * dy - np.sqrt(9 - np.power(j * dx - 18, 2))) / dy

            if velocidadey[i][j] >= 0:
              temperatura = ((80/(a*(a+1))) + (2*malha[i-1][j]/(a+1)) + malha[i][j-1] + malha[i][j+1] + K*dx*(velocidadex[i][j]*malha[i][j-1] + velocidadey[i][j]*malha[i+1][j])) / ((2/a) + 2 + K*dx*(velocidadex[i][j] + velocidadey[i][j]))
                
            elif velocidadey[i][j] < 0:
              temperatura = ((80/(a*(a+1))) + (2*malha[i-1][j]/(a+1)) + malha[i][j-1] + malha[i][j+1] + K*dx*(velocidadex[i][j]*malha[i][j-1] - velocidadey[i][j]*malha[i-1][j])) / ((2/a) + 2 + K*dx*(velocidadex[i][j] - velocidadey[i][j]))

          #atualiza os pontos dentro da condicao inicial, mas que não se encaixam em nenhum dos casos anteriores
          else:
            if velocidadey[i][j] >= 0:
              temperatura = (K*dx*(velocidadex[i][j]*malha[i][j-1] + velocidadey[i][j]*malha[i+1][j]) + malha[i][j+1] + malha[i][j-1] + malha[i-1][j] + malha[i+1][j]) / (K*dx*(velocidadex[i][j] + velocidadey[i][j]) + 4)

            elif velocidadey[i][j] < 0:
              temperatura = (K*dx*(velocidadex[i][j]*malha[i][j-1] - velocidadey[i][j]*malha[i-1][j]) + malha[i][j+1] + malha[i][j-1] + malha[i-1][j] + malha[i+1][j]) / (K*dx*(velocidadex[i][j] - velocidadey[i][j]) + 4)

        #atualiza pontos internos genéricos
        else:
          if velocidadey[i][j] >= 0:
            temperatura = (K*dx*(velocidadex[i][j]*malha[i][j-1] + velocidadey[i][j]*malha[i+1][j]) + malha[i][j+1] + malha[i][j-1] + malha[i-1][j] + malha[i+1][j]) / (K*dx*(velocidadex[i][j] + velocidadey[i][j]) + 4)

          elif velocidadey[i][j] < 0:
            temperatura = (K*dx*(velocidadex[i][j]*malha[i][j-1] - velocidadey[i][j]*malha[i-1][j]) + malha[i][j+1] + malha[i][j-1] + malha[i-1][j] + malha[i+1][j]) / (K*dx*(velocidadex[i][j] - velocidadey[i][j]) + 4)

        # Sobrerrelaxação:
        fator = 1.15
        temp_nova = fator*temperatura + (1-fator)*malha[i][j]
        erro_atual = abs((temp_nova - malha[i][j])/temp_nova) if temp_nova != 0 else 0
        erro_max = max(erro_atual, erro_max) #atualiza com maior erro   
        malha[i][j] = temp_nova

  print("Iterações:", contador)

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_title("Temperatura do domínio (°C)")
  cax = ax.matshow(malha, interpolation='nearest', cmap=cm.inferno)
  fig.colorbar(cax)

  plt.show()

  return(malha)

def parte_2b(malha):
  l = []
  c = []

  #armazena pontos da borda do prédio

  for j in range(round(15//dx) - 1, round(21//dx) + 2):
    for i in range(round(18//dy - 1), round(24//dy)):
      altura = 24 - i * dy
      comprimento = j * dx

      if altura < 3 and (comprimento == 14.625 or comprimento == 21.375):
        l.append(i)
        c.append(j)                

      elif altura - np.sqrt(9 - np.power(j * dx - 18, 2)) - 3 <= dy and altura - np.sqrt(9 - np.power(j * dx - 18, 2)) - 3 > 0:
        l.append(i)
        c.append(j)

  q = 0 #armazena calor total trocado no prédio
  
  for i in range(len(l)):
  #calcula calor total trocado em cada ponto
    derivadax = 0
    derivaday = 0
    calor = 0
    altura = 24 - l[i] * dy
    comprimento = c[i] * dx

    if altura < 3 and comprimento < 18:
      derivadax = -(3 * malha[l[i]][c[i]] - 4 * malha[l[i]][c[i]-1] + malha[l[i]][c[i]-2]) / (2 * dx)

    elif altura < 3 and comprimento > 18:
      derivadax = (-malha[l[i]][c[i]+2] + 4 * malha[l[i]][c[i]+1] - 3 * malha[l[i]][c[i]]) / (2 * dx)
        
    elif altura >= 3 and comprimento <18:
      derivadax = -(3 * malha[l[i]][c[i]] - 4 * malha[l[i]][c[i]-1] + malha[l[i]][c[i]-2]) / (2 * dx)
      derivaday = (-malha[l[i]-2][c[i]] + 4 * malha[l[i]-1][c[i]] - 3 * malha[l[i]][c[i]]) / (2 * dy)

    elif altura >= 3 and comprimento > 18:
      derivadax = (-malha[l[i]][c[i]+2] + 4 * malha[l[i]][c[i]+1] - 3 * malha[l[i]][c[i]]) / (2 * dx)
      derivaday = (-malha[l[i]-2][c[i]] + 4 * malha[l[i]-1][c[i]] - 3 * malha[l[i]][c[i]]) / (2 * dy)

    else:
      derivadax = (malha[l[i]][c[i]+1] - malha[l[i]][c[i]-1]) / (2 * dx)
      derivaday = (-malha[l[i]-2][c[i]] + 4 * malha[l[i]-1][c[i]] - 3 * malha[l[i]][c[i]]) / (2 * dy)

    calor = k * derivadax * dx * comprimento_estrutura + k * derivaday * dy * comprimento_estrutura

    q = q + calor

  print(q)


main()