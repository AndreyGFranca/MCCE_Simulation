import random, math, pylab, os, cmath, time
import numpy as np

def show_conf(L, sigma, title, fname, delxy):
    pylab.axes()
    for [x, y] in L:
        for ix in range(-1, 2):
            for iy in range(-1, 2):
                cir = pylab.Circle((x + ix, y + iy), radius=sigma,  fc='r')
                pylab.plot([x + ix], [y + iy], marker='o', markersize=0.5, color="black")
                pylab.gca().add_patch(cir)
    pylab.axis('scaled')
    #pylab.xlabel('eixo x')
    #pylab.ylabel('eixo y')
    #pylab.title(title)
    major_ticks = np.arange(0, 1.0, delxy)
    pylab.grid(True)
    pylab.xticks(major_ticks)
    pylab.yticks(major_ticks)
    pylab.axis([0.0, 1.0, 0.0, 1.0])
    pylab.savefig(fname)

    pylab.close()

####################################################################

#Lendo os dados do arquivo gerado pelo codigo em C e armazenando em uma lista L
#print L
####################################################################

#Parametros de entrada

N = 8 ** 2 #numero de discos, tem que ser um numero cuja raiz quadrada seja inteira
eta = 0.3 #densidade, ou fator de preenchimento = soma das areas dos discos dividido pela area do quadrado
sigma = math.sqrt(eta / (N * math.pi)) #raio do disco em funcao de N e da densidade eta
N_sqrt = math.sqrt(N)
Q = 20 # numero de iteracoes
delxy = 1.0/ (2.0 * N_sqrt)
two_delxy = 2.0 * delxy

##################################################################

#grafico da configuraco final
titulo1 = '$N=$'+str(N)+', $\eta =$'+str(eta)
nome1 = 'inicial'+'_N_'+str(N) + '_eta_'+str(eta) + '.png'

#grafico da configuraco final
titulo1 = '$N=$'+str(N)+', $\eta =$'+str(eta)
nome1 = 'inicial'+'_N_'+str(N) + '_eta_'+str(eta) + '.png'
for i in range(1000):
    dados = open('resultados/' +str(i), "r")
    L = []
    for linha in dados:
        a, b = linha.split(',')
        L.append([float(a), float(b)])
    show_conf(L, sigma, "figures/"+str(i),"resultados/"+ str(i), two_delxy)

