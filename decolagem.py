'''
Funções para decolagem
Abner Micael de Paula Souza - 10788676
'''
# ====================================
import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
import atmos

plt.close('all')

# ====================================

class decolagem_class():
    
    def __init__(self, decolagem, gerais, pesos):
        
        self.h = decolagem['h']
        self.T = decolagem['T']
        self.P = decolagem['P']
        
        self.rho = atmos.calculate('rho', Tv = self.T, p = self.P)
        
        self.CD0 = gerais['CD0']
        self.k = gerais['k']
        self.S = gerais['S']
        self.c = gerais['c']
        self.T0 = gerais['T0']
        self.n = gerais['n']
        self.CLmax = gerais['CLmax']
        
        self.W = pesos['MTOW']
        return
    
    def velocidades(self):
        
        # Velocidade de estol
        V_s =  np.sqrt(2*self.W/(self.CLmax*self.S*self.rho)) # TODO: ALTERAR
        
        # Velocidade de rotação
        V_R = 1.1*V_s
        
        # Velocidade de 'liftoff'
        V_Lo = 1.15*V_s
        
        # Velocidade V2 (um motor inoperante)
        V2 = 1.25*V_s
        
        return V_s, V_R, V_Lo, V2