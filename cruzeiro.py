'''
Funções para voo em cruzeiro
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

class cruzeiro_class():
    
    def __init__(self, gerais, pesos, cruzeiro, g, beta, rho0):
        
        self.CD0 = gerais['CD0']
        self.k = gerais['k']
        self.S = gerais['S']
        self.c = gerais['c']
        self.T0 = gerais['T0']
        self.n = gerais['n']
        
        self.g = g
        self.beta = beta
        self.rho0 = rho0
  
        self.h = cruzeiro['h']
        self.V = cruzeiro['V']
        self.W1 = cruzeiro['W1']
        self.W2 = cruzeiro['W2']
        self.m =  cruzeiro['m']
        self.W = self.m*self.g
        
        self.BOW = pesos['BOW']
        self.MTOW = pesos['MTOW']
        self.payload_max = pesos['payload_max']
        self.useful_fuel = pesos['useful_fuel']
        return
    
    
    def polar(self, CL):
        CD = self.CD0 + self.k*(CL**2)
        return CD
    
    def CL_funcao(self, V, h, W):
        CL = (2*W)/(self.sigma_funcao(h)*self.rho0*self.S*V**2)
        return CL
    
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
            
    def arrasto_min(self, h, W):
        V = np.sqrt(2*W/(self.sigma_funcao(h)*self.rho0*self.S)*np.sqrt(self.k/self.CD0))
        Dmin = 2*W*np.sqrt(self.CD0*self.k)
        return V, Dmin
    
    def empuxo(self, h):
        T_resp = []
        
        for i in h:
            T = self.T0*self.sigma_funcao(i)**self.n
            T_resp.append(T)
            
        return T_resp
    
    def velocidades_cruzeiro(self, h, V, W):
        Em = np.sqrt(self.CD0/self.k)/(2*self.CD0)
        V1_resp = []
        V2_resp = []
        
        for i in h:
            V1 = np.sqrt((self.T0*self.sigma_funcao(i)**self.n/(self.sigma_funcao(i)*self.rho0*self.S*self.CD0))*(1 - np.sqrt(1 - W**2/(Em**2*(self.T0*self.sigma_funcao(i)**self.n)**2))))
            V2 = np.sqrt((self.T0*self.sigma_funcao(i)**self.n/(self.sigma_funcao(i)*self.rho0*self.S*self.CD0))*(1 + np.sqrt(1 - W**2/(Em**2*(self.T0*self.sigma_funcao(i)**self.n)**2))))
        
            V1_resp.append(V1)
            V2_resp.append(V2)
            
        return V1_resp, V2_resp       
    
    def max_alcance_V_CL(self, V1, h1, W1, W2):
       
        CL1 = self.CL_funcao(V1, h1, W1)
        E1 = CL1/self.polar(CL1)
        zeta = self.zeta_funcao(W1, W2)
        x = (2*V1*E1)*np.log(1/(1 - zeta))/self.c
        
        return x
    
    def max_alcance_V_h(self, V1, h, W1, W2):
        
        CL1 = self.CL_funcao(V1, h, W1) 
        E1 = CL1/self.polar(CL1)
        Em = np.sqrt(self.CD0/self.k)/(2*self.CD0)
        
        zeta = self.zeta_funcao(W1, W2)
        x = (2*V1*Em/self.c)*np.arctan(E1*zeta/(2*Em*(1 - self.k*E1*CL1*zeta))) 
        
        return x
    def zeta_funcao(self, W1, W2):
        
        z = (W1 - W2)/W1
        
        return z
    
    # TODO: arrumar
    def plot_cruzeiro(self, h, V, D_total, D_induzido, D_parasita, V_Dmin, Dmin, V1, V2, T):
    
        if(len(h) == 1):
            plt.figure()
            plt.plot(V, D_total[0], 'k', label = 'Arrasto total')
            plt.plot(V, D_induzido[0], 'r', label = 'Arrasto induzido')
            plt.plot(V, D_parasita[0], 'b', label = 'Arrasto parasita')
            plt.plot(V_Dmin, Dmin, 'og', label = 'Arrasto mínimo')
            plt.xlabel('Velocidade [m/s]')
            plt.ylabel('$Arrasto\:\:[N]$')
            plt.grid(True)
            plt.legend(loc = 'best', framealpha = 1)
            plt.savefig('Figure_1.svg')
            
            
            
            plt.figure()
            plt.plot(V, D_total[0], 'k', label = 'Arrasto total')
            plt.plot(V, [T]*len(V), 'r', label = 'Empuxo')
            plt.plot([V1, V2], [T, T], 'og', label = 'Velocidades de cruzeiro')
            plt.xlabel('Velocidade [m/s]')
            plt.ylabel('$D,T\:\:[N]$')
            plt.grid(True)
            plt.legend(loc = 'best', framealpha = 1)
            plt.savefig('TD_vs_V.svg')
        
        else:
            
            plt.figure()
            plt.plot(V1, h, 'k')
            plt.plot(V2, h, 'r')
            #plt.plot(50.4,12910,'og', label = 'Teto absoluto')
            plt.xlabel('Velocidade [m/s]')
            plt.ylabel('h [m]')
            plt.grid(True)
            #plt.legend(loc = 'best', framealpha = 1)
            plt.savefig('h_vs_M.svg')
        
        return

    def diag_cargaPaga(self, h, M):
        
        V = M*Atmosphere(h).speed_of_sound[0]
        
        # Ponto A:
        '''
        A aeronave não possui combustível nesse ponto, não havendo variação de peso.
        '''
        W1_A = self.BOW + self.payload_max
        W2_A = W1_A
        x_A = self.max_alcance_V_h(V, h , W1_A, W2_A)
        
        # Ponto B:
        '''
        Nesse ponto, adiciona-se combustível até atingir-se o MTOW da aeronave.
        '''
        W1_B = self.MTOW
        W2_B = W1_B - (self.MTOW - (self.BOW + self.payload_max))
        x_B = self.max_alcance_V_h(V,h, W1_B, W2_B)
        
        # Ponto C:
        '''
        Nesse ponto, partiu-se com o máximo combustível permitido.
        '''
        W1_C = self.MTOW
        W2_C = W1_C - self.useful_fuel
        x_C = self.max_alcance_V_h(V, h, W1_C, W2_C)
        
        # Ponto D:
        '''
        Nesse ponto, partiu-se sem carga paga, com o máximo de combustível.
        '''
        W1_D = self.BOW + self.useful_fuel
        W2_D = self.BOW
        x_D = self.max_alcance_V_h(V, h, W1_D, W2_D)
        
        # figura
        plt.figure()
        
        plt.plot([x_A/1000, x_B/1000], [self.payload_max/self.g, self.payload_max/self.g], 'o-k', linewidth = 2)
        plt.plot([x_B/1000, x_C/1000], [self.payload_max/self.g,(self.MTOW - (self.BOW + self.useful_fuel))/self.g], 'o-k', linewidth = 2)
        plt.plot([x_C/1000, x_D/1000], [(self.MTOW - (self.BOW + self.useful_fuel))/self.g, 0], 'o-k', linewidth = 2)
        #plt.plot(x_D/1000, 0)
        
        plt.xlabel("Alcance [km]")
        plt.ylabel("Carga paga [kg]")
        plt.savefig("diag_cargaPaga.svg")
        
        return
# ============================================

