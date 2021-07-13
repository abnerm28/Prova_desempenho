'''
Funções para voo em planeio
Autor: Abner Micael de Paula Souza - 10788676
'''
# ============================================
import numpy as np
import matplotlib.pyplot as plt

from ambiance import Atmosphere

plt.close('all')
# ============================================
class planeio_class():
    
    def __init__(self, gerais, planeio, g, beta, rho0):
        
        self.CD0 = gerais['CD0']
        self.k = gerais['k']
        self.S = gerais['S']
        self.c = gerais['c']
        self.T0 = gerais['T0']
        self.n = gerais['n']
        
        self.V = planeio['V']
        self.m = planeio['m']
        self.h = planeio['h']
        self.h1 = planeio['h1']
        self.h2 = planeio['h2']
        
        self.g = g
        self.beta = beta
        self.rho0 = rho0
        
        self.W = self.m*self.g 
        self.rho = Atmosphere(self.h).density[0]
        
        return
    
    def polar(self, CL):
        CD = self.CD0 + self.k*(CL**2)
        return CD
    
    def CL_funcao(self, V):
        return self.W/(0.5*self.rho*V**2*self.S)
    
    def Eficiencia(self,CL,CD):
        return CL/CD
    
    def gamma_funcao(self, E):
        gamma =  -1/E    
        return gamma
    
    def gamma_min_funcao(self):
        CL = np.sqrt(self.CD0/self.k)
        CD = 2*self.CD0
        Emax = CL/CD
        gamma_min = -1/Emax 
        return gamma_min, Emax
    
    def hponto_funcao(self, V, E):
        hponto = -V/E
        return hponto
        
    def hponto_min_funcao(self, V):
        CL = np.sqrt(3*self.CD0/self.k)
        CD = 4*self.CD0
        Emax = CL/CD
        hponto_min = -V/Emax
        return hponto_min
        
    def max_alcance_CL(self):
    
        E = np.sqrt(self.CD0/self.k)/(2*self.CD0)
        x = E*(self.h2 - self.h1)
        return x
    
    def max_alcance_V(self):
    
        a1 = (self.rho0**2*self.CD0*self.S**2*np.exp(-self.h1/self.beta))/(4*self.W**2*self.k) 
        a2 = (self.rho0**2*self.CD0*self.S**2*np.exp(-self.h2/self.beta))/(4*self.W**2*self.k) 
        Vx = (3/(a1*a2))**(1/8)
        
        x = (self.beta*self.rho0*self.S*(a1 - a2)*Vx**6)/(2*self.W*self.k*(1 + a1*a2*Vx**8))
        
        return x , Vx
    
    def plot_planeio(self):
    
        Emax = self.gamma_min_funcao()[1]
        CL = self.CL_funcao(self.V)
        CD = self.polar(CL)
        
        E = self.Eficiencia(CL, CD)
        hponto = self.hponto_funcao(self.V, E)
        
        plt.figure()
        plt.plot(self.V, hponto, 'k')
        plt.plot(self.V, -self.V/Emax, 'r', label = '$\\gamma_{min}$')
        plt.xlabel('Velocidade [m/s]')
        plt.ylabel('Razão de descida $\\dot{h}\:[m/s]$')
        plt.grid(True)
        plt.legend(loc = 'best', framealpha = 1)
        return

# ============================================

