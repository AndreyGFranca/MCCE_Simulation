import random, math, pylab, os, cmath, time
import matplotlib.pyplot as plt

def show_conf(L, sigma, title, fname):
    pylab.axes()
    for [x, y] in L:
        for ix in range(-1, 2):
            for iy in range(-1, 2):
                cir = pylab.Circle((x + ix, y + iy), radius=sigma,  fc='r')
                pylab.gca().add_patch(cir)
    pylab.axis('scaled')
    pylab.xlabel('eixo x')
    pylab.ylabel('eixo y')
    pylab.title(title)
    pylab.axis([0.0, 1.0, 0.0, 1.0])
    pylab.savefig(fname)
    pylab.close()


def grafico(tempo,Psi,erro,title,nome):
    #plt.plot(tempo, Psi,'-b', lw = 3)
    plt.errorbar(tempo, Psi, yerr = erro, fmt='-o',lw = 2, ms = 8)
    plt.title(title,fontsize=18)
    plt.xticks(color='k', size=14) #k indica cor preta
    plt.yticks(color='k', size=14)
    plt.xlabel('Tempo $t$ (s)', {'color': 'k','fontsize': 16}) #k indica cor preta
    plt.ylabel('Orientacao global $|\Psi_6|$', {'color': 'k','fontsize': 16})
    plt.savefig(nome)
    plt.show()
    plt.close()



####################################################################

#Lendo os dados do arquivo gerado pelo codigo em C e armazenando em uma lista L
#print L
####################################################################

#Parametros de entrada
Nsqrt = 32
N = Nsqrt ** 2 #numero de discos, tem que ser um numero cuja raiz quadrada seja inteira
eta = 0.6 #densidade, ou fator de preenchimento = soma das areas dos discos dividido pela area do quadrado
sigma = math.sqrt(eta / (N * math.pi)) #raio do disco em funcao de N e da densidade eta
Q = 20 # numero de iteracoes
Psi_av = []#list(0.0 for i in xrange(0,Q)) #lista com Q zeros
stdev = []#list(0.0 for i in xrange(0,Q))
tempo = [] #comecando em zero
samples = 5
##################################################################


titulo3 = '$S=$'+str(samples)+', $N=$'+str(Nsqrt)+'$^2$, $\eta =$'+str(eta)
nome3 = 'Psi_av'+'_S_'+str(samples)+'_N_'+str(N)+'_Q_'+str(Q) + '_eta_'+str(eta) + '.png'

psi_dados = open("resultados.csv", "r")
for linha in psi_dados:
    a, b, c = linha.split(',')
    tempo.append(float(a))
    Psi_av.append(float(b))
    stdev.append(float(c))


grafico(tempo,Psi_av,stdev,titulo3,nome3)
