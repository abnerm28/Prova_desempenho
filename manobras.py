'''
Funções para voo em manobra
Autor: Abner Micael de Paula Souza - 10788676
'''
# ============================================
import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere
from scipy.optimize import fsolve

plt.close('all')
sealevel = Atmosphere(0)
rho0 = sealevel.density[0]
# ============================================

class manobras_class():
    
    def __init__(self, gerais, manobras, g, beta, rho0):
        
        self.CD0 = gerais['CD0']
        self.k = gerais['k']
        self.S = gerais['S']
        self.c = gerais['c']
        self.T0 = gerais['T0']
        self.n = gerais['n']
        
        self.g = g
        self.beta = beta
        self.rho0 = rho0
        
        self.h = manobras['h']
        self.V = manobras['V']
        self.m = manobras['m']
        self.fc = manobras['fc']
        
        self.W = self.m*self.g
        return
    
    def CL_funcao(self, fc, V):
        CL = (2*fc*self.W)/(Atmosphere(self.h).density[0]*V**2*self.S)
        return CL
    
    def sigma_funcao(self):
        sigma = Atmosphere(self.h).density[0]/self.rho0
        return sigma
    
    def arrasto_curva(self, fc, V):
        D_curva_resp = []
        label = []
      
        for i in fc:
            CL = self.CL_funcao(i, V)
            D_curva = 0.5*Atmosphere(self.h).density[0]*V**2*self.S*(self.CD0 + self.k*CL**2) 
            D_curva_resp.append(D_curva)
            label.append(i)
            
        return D_curva_resp, label
    
    def raio_minimo(self, T, CLmax):
        ''' 
        Retorna o raio mínimo para uma curva coordenada,
        verificando se a condição de estol é atendida'''
        
        rtt_resp = []
        CLtt_resp = []
        
        Em = np.sqrt(self.CD0/self.k)/(2*self.CD0)
        
        for i in T:
            Vtt = 2*(self.k*(self.W/self.S)/(self.rho0*self.sigma_funcao()*(i/self.W)))**0.5
            ntt = np.sqrt(2 - 1/(Em**2*(i/self.W)**2))
            rtt = Vtt**2/(self.g*np.sqrt(ntt**2 - 1))
            CLtt = (i/self.W)*ntt/(2*self.k)
            
            rtt_resp.append(rtt)
            CLtt_resp.append(CLtt)
        
        for i in np.arange(len(CLtt_resp)):
            if (CLtt_resp[i] > CLmax):
                
                CLtt_resp[i] = CLmax
                E_estol = CLmax/(self.CD0 + self.k*CLmax**2)
                n_estol = (T[i]/self.W)*E_estol
                V_estol = (2*n_estol*(self.W/self.S)/(self.rho0*self.sigma_funcao()*CLmax))**0.5
                rtt_resp[i] = V_estol**2/(self.g*np.sqrt(n_estol**2 - 1))
                
        return  rtt_resp, CLtt_resp
    
    def menor_tempo(self, T, CLmax):
        Em = np.sqrt(self.CD0/self.k)/(2*self.CD0)
        Vft = (2*(self.W/self.S)/(self.rho0*self.sigma_funcao()))**(1/2)*(self.k/self.CD0)**(1/4)
        
        CLft_resp = []
        omegaft_resp, rft_resp = [], []
        
        for i in T:    
            nft = np.sqrt(2*i/self.W*Em - 1)
            CLft = nft*np.sqrt(self.CD0/self.k)
            omegaft = self.g/Vft*np.sqrt(nft**2 - 1)
            rft = Vft**2/(self.g*np.sqrt(nft**2 - 1))
            
            CLft_resp.append(CLft)
            omegaft_resp.append(omegaft)
            rft_resp.append(rft)
        
        for i in np.arange(len(CLft_resp)):
            if (CLft_resp[i] > CLmax):
                
                CLft_resp[i] = CLmax
                E_estol = CLmax/(self.CD0 + self.k*CLmax**2)
                n_estol = (T[i]/self.W)*E_estol
                V_estol = (2*n_estol*(self.W/self.S)/(self.rho0*self.sigma_funcao()*CLmax))**0.5
                omegaft_resp[i] = self.g/V_estol*np.sqrt(n_estol**2 - 1)
                rft_resp[i] = V_estol**2/(self.g*np.sqrt(n_estol**2 - 1))
                
                Vft = V_estol
            
        return  CLft_resp, omegaft_resp, rft_resp, Vft
    
    def estol(self, T, CLmax, cond, fc_mbr):
        
        if(cond == 'analise'):
            E_estol = CLmax/(self.CD0 + self.k*CLmax**2)
            n_estol = (T[0]/self.W)*E_estol
            V_estol = (2*n_estol*(self.W/self.S)/(self.rho0*self.sigma_funcao()*CLmax))**0.5
            omega_estol = self.g/V_estol*np.sqrt(n_estol**2 - 1)
            r_estol = V_estol**2/(self.g*np.sqrt(n_estol**2 - 1))
            
            return V_estol, omega_estol
        
        elif(cond == 'plot'):
            E_estol = CLmax/(self.CD0 + self.k*CLmax**2)
            n_estol = fc_mbr
            V_estol = (2*n_estol*(self.W/self.S)/(self.rho0*self.sigma_funcao()*CLmax))**0.5
            omega_estol = self.g/V_estol*np.sqrt(n_estol**2 - 1)
            r_estol = V_estol**2/(self.g*np.sqrt(n_estol**2 - 1))     
            
            return V_estol, omega_estol
            
    def cinematica(self, V, fc):
        ''' 
        A velocidade (V) de entrada é aquela em que T = D, para um
        dado fator de carga (fc). '''
        
        omega_resp, raio_resp = [], []
        
        for i in np.arange(len(fc)):
            omega = self.g/V[i]*np.sqrt(fc[i]**2 - 1)
            raio = V[i]**2/self.g*(np.sqrt(fc[i]**2 - 1))**(-1)
            omega_resp.append(omega)
            raio_resp.append(raio)
                
        return omega_resp, raio_resp
    
    def fc_maximo(self, fc, V):
        
        T_aux = self.empuxo_manobras(fc)
        D_aux = self.arrasto_curva(fc, V)[0]
        
        i = 0
        while(max([T_aux[i]]*len(V) - D_aux[i]) > 0):
            i += 1   
        
        fc_max = [fc[i + 1]]
        
        return fc_max
    
    def empuxo_manobras(self, fc):
        T_resp = []
            
        for i in fc:
            T = self.T0*self.sigma_funcao()**self.n
            T_resp.append(T)
            
        return T_resp
    
    def velocidades_manobras_eq(self, x, fc):
        ''' 
        Monta o sistema de equações a serem resolvidas '''
        
        V = x
        CL = self.CL_funcao(fc, V)
        
        T = self.T0*self.sigma_funcao()**self.n
        D = 0.5*Atmosphere(self.h).density[0]*V**2*self.S*(self.CD0 + self.k*CL**2)
        
        eq = T-D
    
        return eq
    
    def velocidades_manobras_solver(self, V1, V2, fc):
        ''' 
        Calcula quais as possíveis velocidades de manobra (os extremos) para
        uma determinada altitude (h) '''
        
        V1_man_resp, V2_man_resp = [], []
        
        V1_init = V1
        V2_init = V2
        
        for i in fc:
            V1_man = fsolve(self.velocidades_manobras_eq, V1_init , args = (i))
            V2_man = fsolve(self.velocidades_manobras_eq, V2_init , args = (i))
    
            V1_man_resp.append(V1_man)
            V2_man_resp.append(V2_man)
            
            V1_init = V1_man
            V2_init = V2_man
            
        return V1_man_resp, V2_man_resp
    
    def plot_manobra(self, V1_manobras, V2_manobras, omega1_manobras,
                     omega2_manobras, V, fc_max, fc, T, CLmax):
    
        [Vestol, omega_estol] = self.estol(T, CLmax, 'plot', fc)
        
        plt.figure()
        plt.xlabel('Velocidade linear [m/s]')
        plt.ylabel('Velocidade angular [rad/s]')
        plt.grid(True)
        plt.plot(V1_manobras, omega1_manobras, 'k')
        plt.plot(V2_manobras, omega2_manobras, 'k')
        plt.plot(V, self.cinematica(V, fc_max*len(V))[0], label = '$n_m$ = {:.2f}'.format(fc_max[0]))
        plt.plot(Vestol, omega_estol, 'r', label = 'Estol')
        plt.legend(loc = 'best', framealpha = 1)
        plt.ylim(top = 0.5)
        plt.savefig('desempenho_em_curva.svg')
        return
# ============================================