## this class defines the "memory terms" in QUAPI IF
## to initialize, need time-step of t and t'' (delt for all except at the 2 end points),
## need jw spectral density of bath oscillators
import numpy as np
from scipy.integrate import quad

class g:

    ### redefine params esp t, dt, k
    def __init__(self,t,dt,k,jw,s,temp):
        self.t = t
        self.dt = dt
        self.k = k
        self.jw = jw
        self.s = s
        self.temp = temp

    def im_gplus_t(self,t):
        a = quad(lambda w: self.jw(w)*(w**-2)*((1+np.exp(w/self.temp))**(-1))*(np.sin(w*t)-w*t),0,np.inf,epsabs=1.49e-15)
        # if a[1]/a[0] >= 0.05:
        #     print("----- Large error in imEta t -----")

        return a[0]

    def re_gplus_t(self,t):

        a = quad(lambda w: self.jw(w)*(w**-2)*((1+np.exp(w/self.temp))**(-1))*(1-np.cos(w*t)),0,np.inf,epsabs=1.49e-15)
        #print('----t----',t)
        #print('------ eta 1------',eta1)
        # if eta1[1]/eta1[0] >= 0.05:
        #     print("----- Large error in reEta t -----")

        return a[0]

    def im_gminus_t(self,t):
        a = quad(lambda w: self.jw(w)*(w**-2)*((1+np.exp(-w/self.temp))**(-1))*(np.sin(w*t)-w*t),0,np.inf,epsabs=1.49e-15)
        # if a[1]/a[0] >= 0.05:
        #     print("----- Large error in imEta t -----")

        return a[0]

    def re_gminus_t(self,t):

        a = quad(lambda w: self.jw(w)*(w**-2)*((1+np.exp(-w/self.temp))**(-1))*(1-np.cos(w*t)),0,np.inf,epsabs=1.49e-15)
        #print('----t----',t)
        #print('------ eta 1------',eta1)
        # if eta1[1]/eta1[0] >= 0.05:
        #     print("----- Large error in reEta t -----")

        return a[0]


    def gplus_t(self,t):
        #print(self.reEta_t(t))
        return 1*self.re_gplus_t(t) + self.im_gplus_t(t)*1J


    def gminus_t(self,t):
        #print(self.reEta_t(t))
        return 1*self.re_gminus_t(t) + self.im_gminus_t(t)*1J


    def gplus_arr (self):
        garr = np.zeros(self.k,dtype=complex)

        for i in np.arange(self.k):
            garr[i] = self.gplus_t(i*self.dt)

        return garr

    def gplus_arr_half (self):
        garr_half = np.zeros(self.k,dtype=complex)

        for i in np.arange(self.k):
            garr_half[i] = self.gplus_t((i+0.5)*self.dt)
        return garr_half


    def gminus_arr (self):
        garr = np.zeros(self.k,dtype=complex)

        for i in np.arange(self.k):
            garr[i] = self.gminus_t(i*self.dt)

        return garr

    def gminus_arr_half (self):
        garr_half = np.zeros(self.k,dtype=complex)

        for i in np.arange(self.k):
            garr_half[i] = self.gminus_t((i+0.5)*self.dt)
        return garr_half
