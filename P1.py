"""
===========
Clusters
===========
"""

import math
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import os
import json
import numpy as np
from matplotlib.ticker import ( AutoLocator, AutoMinorLocator)

"Constantes do Material"

Sigma_CTE = 0;              "J/m² -> energia de superfície - Interface cristal - líquido super resfriado"
Q_CTE  =    0;              "K/m -> taxa aquecimento"
Tm_CTE =    0;              "Kelvin Temp de fusão"
Tg_CTE =    0;              "Kelvin Tem transição vítrea"
Hm_CTE =    0;
Cs_CTE =    0;
D_CTE =     0;              "metros -> d atômico"

"Parâmetros empíricos equação VFT"
A_CTE =     0;           
B_CTE =     0;              "Kelvin"
To_CTE =    0;              "Kelvin"


"Constantes"
NAv_CTE = 6.0221e23;        "mol^-1 -> Avogadro"
R_CTE = 8.3144621;          "J/mol.K -> cte dos gases"
Kb_CTE = 1.3806488e-23;     "m².kg/s²K -> boltsmann"
H_CTE = 6.626e-34;          "J.s -> Plank"

"Importa os parâmetros do vidro"
if os.path.isfile('P1.json'):
    with open('P1.json', "r") as f:
        my_dict = json.load(f)
    locals().update(my_dict)
else:
    print("Não foi possível carregar as variáveis, o software será fechado")
    os._exit(1)

"Cria o diretório de exportação"
if not os.path.exists("Export"):
    os.mkdir("Export")

"Calcula o Volume molar"
Vm_CTE = D_CTE**3 * NAv_CTE;    "volume molar do cristal"



"""
===========
BLOCO DE EQUAÇÕES (TIPO LAMBDA)
===========
"""

# VISCOSIDADE
"Eq VFT - Viscosidade - Pa.s"
funcVFT =               lambda Ty: pow(10, A_CTE + (B_CTE/(Ty-To_CTE)))

# NUCLEAÇÃO
"∆Gv - Russel"
dGv =                   lambda C1, T, C2: C1 - T*C2
"I Nucleação - Aproximada"
Inucleacao =            lambda T: 3698966000 * math.exp(-(T - 460.9398)**2/(2*12.36206**2))

# CRESCIMENTO
"Eq Variação energia livre -> dGibbs"
funcGibbs =             lambda Ty: Hm_CTE * (Tm_CTE - Ty) * Ty / (Tm_CTE**2)
"Eq Sitios pref de crescimento -> fy"
funcSitiosPref =        lambda dGibbs: Cs_CTE*D_CTE*dGibbs / (4 * math.pi * Sigma_CTE * Vm_CTE)
"Eq difusão efetiva -> zy"
funcDifusaoEfetiva =    lambda Ty, ny: Kb_CTE * Ty / (D_CTE * ny)
"Eq taxa de crescimento de cristais -> uy"
funcCrescCrist =        lambda zy, fy, dGibbs, Ty: (zy/D_CTE)*fy*(1-math.exp((-dGibbs/(R_CTE*Ty))))

# TTT
X_TTT =                 lambda I, u, t: (math.pi/3) * I * u**3 * t**4
t_TTT =                 lambda I, u, X: ((3 * X) / ( math.pi * I * u**3 ))**(1/4)


"""
===========
BLOCO DE FUNÇÕES DE PLOTAGEM
===========
"""

# Função de plotagem das imagens - 2 parâmetros opcionais apenas para quando for plotar 2 curvas juntas
def plotExport(n, xPlot, yPlot, typePlot, xLeg, yLeg, titulo, secondXPlot=0, secondYPlot=0, thirdXPlot=0, thirdYPlot=0, label1=None, label2=None, label3=None, xLimLeft=None, xLimRight=None, yLimTop=None, yLimBottom=None):
    print("Salvo:   " + str(n) + " - " + titulo + ".png")
    fig, ax = plt.subplots()
    if(typePlot == "Linear"):
        ax.plot(xPlot, yPlot, label=label1)
        # Checa se há uma segunda curva a ser plotada e se ela compartilha o eixo Y
        if (secondXPlot != 0):
            if (secondYPlot != 0):
                ax.plot(secondXPlot, secondYPlot, label=label2)
            else:
                ax.plot(secondXPlot, yPlot, label=label2)
    elif(typePlot == "Log"):
        ax.semilogy(xPlot, yPlot, label=label1)
        ax.get_yaxis().set_major_locator(LogLocator(numticks=12))
        # Checa se há uma segunda curva a ser plotada e se ela compartilha o eixo Y
        if secondXPlot != 0:
            if (secondYPlot != 0):
                ax.plot(secondXPlot, secondYPlot, label=label2)
            else:
                ax.plot(secondXPlot, yPlot, label=label2)
        # Checa se há uma terceira curva a ser plotada e se ela compartilha o eixo Y
        if thirdXPlot != 0:
            if (thirdYPlot != 0):
                ax.plot(thirdXPlot, thirdYPlot, label=label3)
            else:
                ax.plot(thirdXPlot, thirdYPlot, label=label3)
        
    else:
        ax.semilogx(xPlot, yPlot, label=label1)
        # Checa se há uma segunda curva a ser plotada e se ela compartilha o eixo Y
        if secondXPlot != 0:
            if (secondYPlot != 0):
                ax.semilogx(secondXPlot, secondYPlot, label=label2)
            else:
                ax.semilogx(secondXPlot, yPlot, label=label2)
        # Checa se há uma terceira curva a ser plotada e se ela compartilha o eixo Y
        if thirdXPlot != 0:
            if (thirdYPlot != 0):
                ax.semilogx(thirdXPlot, thirdYPlot, label=label3)
            else:
                ax.semilogx(thirdXPlot, thirdYPlot, label=label3)
                
    if xLimLeft != None:
        plt.xlim(left=xLimLeft)
    if xLimRight != None:
        plt.xlim(right=xLimRight)
    if yLimTop != None:
        plt.ylim(top=yLimTop)
    if yLimBottom != None:
        plt.ylim(bottom=yLimBottom)
        
    ax.set(xlabel=xLeg, ylabel=yLeg, title=titulo)
    ax.grid()
    ax.legend(prop={'size': 18})

    plt.gcf().set_size_inches(8, 10)
    fig.savefig("Export/" + str(n) + " - " + titulo + ".png", dpi = 300)
    plt.close(fig)
    #plt.show()


def plot_2_different_scales(n, xPlot, yPlot, yPlot2):
    print("Salvo:   " + str(n) + ".png")
    fig, ax1 = plt.subplots()
    
    color = 'tab:red'
    ax1.set_xlabel("Temperatura (K)")
    ax1.set_ylabel("Taxa de Nucleação (1/m³.s)", color=color)
    ax1.semilogy(xPlot, yPlot, color=color, label="Taxa de Nucleação")
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_ylim([1E-20,1E10])
    
    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel("Cresc. Cristais - U (m/s)", color=color)
    ax2.semilogy(xPlot, yPlot2, color=color, label="Cresc. Cristais")
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_ylim([1E-20,1E10])
    
    ax1.grid()
    ax1.set(title="Nucleação e Crescimento de Cristais vs Temperatura")
    fig.tight_layout()
    plt.gcf().set_size_inches(9, 5)
    fig.savefig("Export/" + str(n) + ".png", dpi = 300)
    plt.close(fig)
    
def plot_3_curves_in_the_same_image(n, xPlot, yPlot, yPlot2, yPlot3):
    print("Salvo:   " + str(n) + ".png")
    fig, ax = plt.subplots()
    fig.subplots_adjust(bottom=0.2, right=0.85)
    
    newax = fig.add_axes(ax.get_position())
    newax.patch.set_visible(False)
    newax.yaxis.set_label_position('right')
    newax.yaxis.set_ticks_position('right')
    newax.spines['bottom'].set_position(('outward', 50))
    
    ax.semilogx(yPlot, xPlot, 'r-', label="Nucleação (1/m³.s)")
    ax.semilogx(yPlot2, xPlot,'b-', label="Crescimento (m/s)")
    ax.set_xlabel('Nucleação (1/m³.s)', color='red')
    ax.set_xlabel('Crescimento (m/s)', color='blue')
    ax.set_ylabel('Temperatura (K)', color='black')
    ax.set_xlim([1E-20,1E40])
    ax.set_ylim([600,1300])
    ax.grid()
    ax.legend(loc="upper left", prop={'size': 8})
    
    newax.semilogx(yPlot3, xPlot, 'g-', label="TTT")
    newax.set_xlim([1E-20,1E20])
    newax.set_ylim([600,1300])
    newax.legend(loc="upper right", prop={'size': 8})
    newax.set_xlabel('Tempo (s)', color='green')

    plt.gcf().set_size_inches(7, 10)
    fig.savefig("Export/" + str(n) + ".png", dpi = 300)
    plt.close(fig)
    
      
def plot_CCT_Resfriamento(n, xPlot, yPlot):
    print("Salvo:   " + str(n) + ".png")
    
    fig, ax = plt.subplots()
    ax.semilogx(xPlot, yPlot)
        
    ax.set(xlabel="TTT (s)", ylabel="Temperatura (K)", title="Diagrama CCT de Resfriamento")
    ax.grid()
    ax.legend(prop={'size': 18})
    
    taxas = [7E-4, 1E-3, 1E-2, 1E-1, 1E0, 1, 10]
    for Q in taxas:
        CCT = []
        CCT_Tempo = []
        Tempo = 0
        Temperature = Tm_CTE
        while Temperature > Tg_CTE-127:
            Tempo += 1
            CCT_Tempo.append(Tempo)
            CCT.append(Temperature)
            Temperature -= Q
        labelTaxa = str(Q) + " K/s"
        ax.semilogx(CCT_Tempo, CCT, label=labelTaxa)

    ax.legend(prop={'size': 18})
    plt.gcf().set_size_inches(8, 10)
    fig.savefig("Export/" + str(n) + ".png", dpi = 300)
    plt.close(fig)
    #plt.show()
    
def plot_CCT_Aquecimento(n, xPlot, yPlot):
    print("Salvo:   " + str(n) + ".png")
    
    fig, ax = plt.subplots()
    ax.semilogx(xPlot, yPlot)
        
    ax.set(xlabel="TTT (s)", ylabel="Temperatura (K)", title="Diagrama CCT de Aquecimento")
    ax.grid()
    ax.legend(prop={'size': 18})
    
    taxas = [2.5E-4, 1E-3, 1E-2, 1E-1, 1E0, 1, 10]
    for Q in taxas:
        CCT = []
        CCT_Tempo = []
        Tempo = 0
        Temperature = Tg_CTE-127
        while Temperature < Tm_CTE:
            Tempo += 1
            CCT_Tempo.append(Tempo)
            CCT.append(Temperature)
            Temperature += Q
        labelTaxa = str(Q) + " K/s"
        ax.semilogx(CCT_Tempo, CCT, label=labelTaxa)

    ax.legend(loc="upper right", prop={'size': 8})
    plt.gcf().set_size_inches(8, 10)
    fig.savefig("Export/" + str(n) + ".png", dpi = 300)
    plt.close(fig)
    #plt.show()
    
    

"""
===========
PROGRAMA PRINCIPAL
===========
"""

"Calcula uma lista de 'tempo'"
temp_y = list(range(Tg_CTE-100, Tm_CTE, 1))
 
#Define as listas da parte 1
dGibbs = [];  
dGibbs2 = [];

fy = [];
ny = [];
zy = [];
uy = [];
TTT = [];
TTT6= [];

I_nucleacao = [];


# Calcula valores em função da TEMPERATURA
for T in temp_y:    
    """
    ===========
    VISCOSIDADE
    ===========
    """
    "Calcula a viscosidade no tempo y - Pa.s"
    ny_step = funcVFT(T)
    ny.append(ny_step)
    
    
    """
    ===========
    NUCLEAÇÃO
    ===========
    """
    "Calcula a variação de energia livre de Gibbs no tempo y ATRAVÉS DO MÉTODO X"
    dG_step = dGv(53368, T, 39.37)
    dGibbs.append(dG_step)
    
    I_nucleacao_step = Inucleacao(T-273)
    I_nucleacao.append(I_nucleacao_step)
    
    
    """
    ===========
    CRESCIMENTO
    ===========
    """
    "Calcula ΔGibbs ATRAVÉS DO MÉTODO X"
    dGibbs_step = funcGibbs(T)
    dGibbs2.append(dGibbs_step)
    
    "Calcula o número de sítios prefereniais de crescimento"
    fy_step = funcSitiosPref(dGibbs_step)
    fy.append(fy_step)
    
    "Calcula a difusão efetiva"
    zy_step = funcDifusaoEfetiva(T, ny_step)
    zy.append(zy_step)
    
    "Calcula a taxa de crescimento de cristais"
    uy_step = funcCrescCrist(zy_step, fy_step, dGibbs_step, T)
    uy.append(uy_step)  
    
    
    """
    ===========
    TTT
    ===========
    """
    if I_nucleacao_step < 1E-30 or uy_step == 0:
        TTT_step = float('inf')
        TTT_step6 = float('inf')
    else:
        TTT_step = t_TTT(I_nucleacao_step, uy_step, 1E-3)
        TTT_step6 = t_TTT(I_nucleacao_step, uy_step, 1E-6)
        #TTT_step = math.log(TTT_step)
        
    TTT.append(TTT_step)
    TTT6.append(TTT_step6)
    
# Plota as curvas comuns
plotExport(1, temp_y, dGibbs,   "Linear",   "Temperatura (K)",  "ΔGibbs (J/mol)", "ΔGibbs vs Temperatura")
plotExport(1, temp_y, dGibbs2,  "Linear",   "Temperatura (K)",  "ΔGibbs (J/mol)", "ΔGibbs2 vs Temperatura")
plotExport(2, temp_y, ny,       "Log",      "Temperatura (K)",  "Viscosidade (Pa.s)",     "Viscosidade vs Temperatura")
plotExport(3, temp_y, uy,       "Log",      "Temperatura (K)",  "Cresc. Cristais - U (m/s)",      "Taxa de Crescimento de Cristais vs Temperatura")
plotExport(4, temp_y[62:153], I_nucleacao[62:153], "Log", "Temperatura (K)", "Taxa de Nucleação (1/m³.s)",     "Taxa de Nucleação")
plotExport(5, TTT[40:680], temp_y[40:680],  "LogX", "Tempo (s)",    "Temperatura (K)",     "Curva TTT", secondXPlot=TTT6[40:680], yLimBottom=Tg_CTE-200, yLimTop=Tm_CTE, xLimLeft=1, label1="α = 1E-3", label2="α = 1E-6")
plotExport(6, TTT, temp_y,      "LogX",     "Tempo (s)",            "Temperatura (K)",     "Curva TTT Com Crist", uy, temp_y, I_nucleacao, temp_y, xLimLeft=1E-20)
plotExport(7, temp_y, uy,       "Log",      "Temperatura (K)",  "Cresc. Cristais - U (m/s)",      "Taxa de Crescimento de Cristais vs Temperatura", temp_y, I_nucleacao)

# Plota nucleação e crescimento
plot_2_different_scales(8, temp_y, I_nucleacao, uy)

# Plota nucleação, crescimento e TTT
plot_3_curves_in_the_same_image(9, temp_y, I_nucleacao, uy, TTT)

# Plota CCT
plot_CCT_Resfriamento(10, TTT, temp_y)
plot_CCT_Aquecimento(11, TTT, temp_y)