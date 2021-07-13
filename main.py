'''
Main - Desempenho de aeronaves
'''
# ========================
import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere

from planeio import planeio_class
from cruzeiro import cruzeiro_class
from ascendente import ascendente_class
from manobras import manobras_class
from decolagem import decolagem_class
# ========================
plt.close('all')
sealevel = Atmosphere(0)
rho0 = sealevel.density[0]
g = 9.81
beta = 9296  # m
X = 62
# ========================

gerais = {'T0': 2*89000, 'n': 0.85 + 0.15*X/99, 'CD0':0.025 + X/20000,
          'k': 0.042, 'S': 92.53, 'c': 0.64/3600, 'CLmax': 1.78}

pesos = {'BOW': 28667*g, 'MTOW': 53000*g, 'payload_max': 16150*g,  
         'useful_fuel': 12970*g }

planeio = {'V': np.linspace(20,300), 'h': [1524], 'm': 28667 + 16150 + 4000,
           'h2': 3048, 'h1': 0}

ascendente = {'h': [0], 'V': np.linspace(0,400,500, endpoint=True),
              'm': 53000}

cruzeiro = {'h': [10972.8], 'V': np.linspace(20,400,500, endpoint=True),
            'm': 42000, 'W1': 42000*g, 'W2': (42000 - 10000)*g}

manobras = {'h': [914.4], 'V': np.linspace(80,400,500, endpoint=True), 
            'm': 37400, 'fc': np.linspace(1.0001,10,1000, endpoint = True)}

decolagem = {'h': 661, 'T': 8 + 273.15, 'P':  102600}

# ========================

# ----------------------------- QUESTÃO 1 ------------------------------------
pl = planeio_class(gerais, planeio, g, beta, rho0)

# Item A)

xmax_CL = pl.max_alcance_CL()
[xmax_V, Vxmax] = pl.max_alcance_V()

# Item B)

CL_pl = pl.CL_funcao(Vxmax)
CD_pl = pl.polar(CL_pl)
E_pl  = pl.Eficiencia(CL_pl, CD_pl)
hponto_pl = pl.hponto_funcao(Vxmax, E_pl)

# ----------------------------- QUESTÃO 2 ------------------------------------
asc = ascendente_class(gerais, ascendente, g, beta, rho0)

# Item A)

T_asc = asc.empuxo(ascendente['h'])
V_max_hponto = asc.vel_max_razao_subida(asc.W, T_asc, ascendente['h'])
D_asc = asc.arrasto_analise(ascendente['h'], V_max_hponto, asc.W)[2]
hponto_asc = asc.razao_subida(ascendente['h'], T_asc, D_asc, V_max_hponto, asc.W)

# Item B)

[teto_abs, teto_servico, teto_operacional] = asc.tetos(asc.W, tol = 0.1)

# ----------------------------- QUESTÃO 3 ------------------------------------
cr = cruzeiro_class(gerais, pesos, cruzeiro, g, beta, rho0)

# Item A)
[V1_cr,V2_cr] = cr.velocidades_cruzeiro(cruzeiro['h'], cruzeiro['V'], cruzeiro['m']*9.81)

# Item B)
x_v_CL = cr.max_alcance_V_CL(np.array(V2_cr), cruzeiro['h'], cruzeiro['W1'], cruzeiro['W2'])

# Item C)
#cr.diag_cargaPaga(10972.8, 0.8)

# ----------------------------- QUESTÃO 4 ------------------------------------
mbr = manobras_class(gerais, manobras, g, beta, rho0)

# Item A)

T_mbr = mbr.empuxo_manobras([1]*len(manobras['h']))
[rtt, CLtt] = mbr.raio_minimo(T_mbr, CLmax = 1.78)


# Item B)
omega_ft = mbr.menor_tempo(T_mbr, CLmax = 1.78)[1]
t = 2*np.pi/omega_ft[0]  # [s]

# Item C)
Vft = mbr.menor_tempo(T_mbr, CLmax = 1.78)[3]
rft = mbr.menor_tempo(T_mbr, CLmax = 1.78)[2]

CLft = mbr.menor_tempo(T_mbr, CLmax = 1.78)[0]  # Check CL

'''
# EXTRA: Diagrama Desempenho em Curva
fc_max_mbr = mbr.fc_maximo(manobras['fc'], manobras['V'])
fc_mbr = np.linspace(1.0001, fc_max_mbr, 100, endpoint= True)  # Não alterar

[D_mbr, label1] = mbr.arrasto_curva(manobras['fc'], manobras['V'])
[V1_mbr,V2_mbr] = cr.velocidades_cruzeiro(manobras['h'], manobras['V'], mbr.W)
[V1_manobras, V2_manobras] = mbr.velocidades_manobras_solver(V1_mbr, V2_mbr, fc_mbr)
omega1_manobras = mbr.cinematica(V1_manobras, fc_mbr)[0]
omega2_manobras = mbr.cinematica(V2_manobras, fc_mbr)[0]


mbr.plot_manobra(V1_manobras, V2_manobras, omega1_manobras, omega2_manobras, 
                 manobras['V'], fc_max_mbr, fc_mbr, T_mbr, CLmax = 1.78)
'''
# ----------------------------- QUESTÃO 5 ------------------------------------
dec = decolagem_class(decolagem, gerais, pesos)
[V_s, V_R, V_Lo, V2] = dec.velocidades()

# Item A)

# Item B)

