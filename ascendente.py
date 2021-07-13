'''
Funções para voo ascendente
Autor: Abner Micael de Paula Souza - 10788676
'''
# ============================================
import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere

plt.close('all')
sealevel = Atmosphere(0)
rho0 = sealevel.density[0]
# ============================================
class ascendente_class():
    
    def __init__(self, gerais, ascendente, g, beta, rho0):
        
        self.CD0 = gerais['CD0']
        self.k = gerais['k']
        self.S = gerais['S']
        self.c = gerais['c']
        self.T0 = gerais['T0']
        self.n = gerais['n']
        
        self.g = g
        self.beta = beta
        self.rho0 = rho0
        
        self.h = ascendente['h']
        self.V = ascendente['V']
        self.m = ascendente['m']
        
        self.W = self.m*self.g 
        self.rho = Atmosphere(self.h).density[0]
        
        return
    
    def sigma_funcao(self, h):
        sigma = Atmosphere(h).density[0]/self.rho0
        return sigma
    
    def arrasto_analise(self, h, V, W):
        '''
        Diferentemente da função polar(), essa função quebra
        o arrasto em duas parcelas: parasita + induzido.
        '''
        D_parasita_resp, D_induzido_resp, D_total_resp = [],[],[]
    
        for i in h:
            D_parasita = 0.5*self.sigma_funcao(i)*self.rho0*(V**2)*self.S*self.CD0
            D_induzido = (2*self.k*W**2)/(self.sigma_funcao(i)*self.rho0*self.S*V**2)
            D_total = D_parasita + D_induzido
            
            D_parasita_resp.append(D_parasita)
            D_induzido_resp.append(D_induzido)
            D_total_resp.append(D_total)
        
        return D_parasita_resp, D_induzido_resp, D_total_resp
    
    def empuxo(self, h):
        T_resp = []
        
        for i in h:
            T = self.T0*self.sigma_funcao(i)**self.n
            T_resp.append(T)
            
        return T_resp
    
    def gamma_funcao(self, T, D, W, h, V):
        gamma_resp = []
        
        for i in np.arange(0, len(h)):
            gamma = ([T[i]]*len(V) - D[i])/W
            gamma_resp.append(gamma)
        return gamma_resp
    
    def razao_subida(self, h, T, D, V, W):
        hponto_resp = []
        for i in np.arange(0, len(h)):
    
            hponto = (([T[i]]*len(V))*V - D[i]*V)/W
            hponto_resp.append(hponto)
        return hponto_resp
    
    def vel_max_razao_subida(self, W, T, h):
        
        V1_resp = np.zeros(len(h))
        
        for i in np.arange(len(h)):
            
            rho = Atmosphere(h[i]).density[0]
            V1 = np.sqrt((T[i]/self.S)/(3*rho*self.CD0)*(1 + np.sqrt(1 + 12*self.CD0*self.k*(W/T[i])**2)))
            V1_resp[i] = V1
            
        return V1_resp
    
    def tetos(self, W, tol):
        
        h = np.linspace(0,15000,1000, endpoint = True)
        D_resp = np.zeros(len(h))
        
        T = self.empuxo(h)
        Vmax = self.vel_max_razao_subida(W, T, h)
        
        for j in np.arange(len(h)):
            D = self.arrasto_analise([h[j]], Vmax[j], W)[2]
            D_resp[j] = D[0]
        
        hponto_max = (T*Vmax - D_resp*Vmax)/W
    
        for i in np.arange(len(hponto_max)):
            if(abs(hponto_max[i]) <= tol):
                teto_absoluto =h[i]
                
            elif(abs(hponto_max[i] - 0.508) <= tol):
                teto_servico = h[i]
            
            elif(abs(hponto_max[i] - 1.524) <= tol):
                teto_operacional = h[i]
                
        return teto_absoluto, teto_servico, teto_operacional
    
    def plot_ascendente(self, V, h, W):
    
        T = self.empuxo(h)
        D = self.arrasto_analise(h, V, W)[2]
        hponto = self.razao_subida(h, T, D, V, W)
        
        plt.figure()
        for i in np.arange(0,len(h)):
            plt.plot(V,hponto[i], label = 'h = {:.2f}'.format(h[i]))
        plt.xlabel('Velocidade [m/s]')
        plt.ylabel('$\\dot{h}\:\:[m/s]$')
        plt.grid(True)
        plt.ylim(bottom = 0)
        plt.legend(loc = 'best', framealpha = 1)
        plt.savefig('razaoSubida_V.svg')
        plt.show()
        return


# ============================================
