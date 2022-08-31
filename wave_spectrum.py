'''
DISCLAIMER: This is a Python implementation of a MATLAB script. 
The original description is posted below:

Function to evaluate different type of wave spectra.

Use:  S = wavespec(SpecType,Par,W,PlotFlag)

Input:
SpecType	- Spectrum type
Par		- Spectrum parameters
W			- Column vector of wave frequencies [rad/s]
PlotFlag	- 1 to plot the spectrum, 0 for no plot

Output:
S			- Column vector of wave spectrum values [m^2 s]; evaluated at W(k) 
    
SpecType and Par =[p1,p2,p3,...pn]:
SpecType =1  Bretschneitder (p1=A,p2=B)
SpecType =2  Pierson-Moskowitz   (p1=Vwind20) 
SpecType =3, ITTC-Modified Pierson-Moskowitz (p1=Hs,p2=T0)
SpecType =4, ITTC-Modified Pierson-Moskowitz (p1=Hs,p2=T1)
SpecType =5, ITTC-Modified Pierson-Moskowitz (p1=Hs,p2=Tz)
SpecType =6, JONSWAP (p1=Vwind10,p2=Fetch)
SpecType =7, JONSWAP (p1=Hs,p2=w0,p3=gamma)
SpecType =8, Torsethaugen (p1=Hs,p2=w0) 

Bretschneither:
p1 = A
p2 = B
S(w)=A*w^(-5)*exp(-B/(w^4));
Reference [6]

Pierson-Moskowitz:
p1 = Vwind20 --Average wind speed @20m above sea level [m/s]
A=8.1e-3*9.81^2;  
B=0.74*(9.81/Vwind20)^4;
S(w)=A*w^(-5)*exp(-B/(w^4)); [m^2 s]
Reference [5,6]


ITTC-Modified Pierson-Moskowitz (Hs,T0):
p1 = Hs  --Significant wave height (Hs = 4 m0^1/2) [m]
p2 = T0  --Modal period      (T0 = 2 pi /w0)       [s]
A=487*Hs^2/T0^4;
B=1949/T0^4;
S(w)=A*w^(-5)*exp(-B/(w^4)); [m^2 s]
Reference [1]


ITTC-Modified Pierson-Moskowitz (Hs,T1):
p1 = Hs  --Significant wave height (Hs = 4 m0^1/2) [m]
p2 = T1   --Average wave period (T1 = 2 pi m0/m1)     [s]
A=173*Hs^2/T1^4;
B=691/T1^4;
S(w)=A*w^(-5)*exp(-B/(w^4)); [m^2 s]
Reference [1]


ITTC-Modified Pierson-Moskowitz (Hs,Tz):
p1 = Hs --Significant wave height (Hs = 4 m0^1/2)             [m]
p2 = Tz --Average zero-crossing period (T = 2 pi (m0/m1)^1/2) [s]
A=123*Hs^2/Tz^4;
B=495/Tz^4;
S(w)=A*w^(-5)*exp(-B/(w^4)); [m^2 s]
Reference [1]


JONSWAP (Vwind10,Fetch):
p1 = Vwind10 --wind speed @ 10m over sea surface [m/sec]
p2 = fetch --distance to georaphical boundary [m]
g=9.81;
xtilde= g*fetch/(Vwind10^2);
f0=3.5*(g/Vwind10)*xtilde^-0.33;
w0=2*pi*f0;
alpha=0.076*xtilde^-0.22;
gamma =3.3;
sigma=0.07  if w<w0, sigma=0.09 otherwise;
S(w)=S1*S2 -  [m^2 sec]  
with,
S1=alpha*g^2*(W^-5)*exp(-(5/4)*(w0/w)^4);
S2=gamma^(exp(-(w-w0)^2/(2*(sigma*w0)^2)));   
Reference [2]


JONSWAP (Hs,w0, gamma):
p1 = Hs - Significant wave height (Hs = 4 sqrt(m0)) [m]
p2 = w0 - Modal Freq. [rad/sec] (Recomended 1.25<w0*sqrt(Hc)<1.75)
p3 = gamma - Peakedness factor (Recommended between 1 and 5; usually 3.3, set to zero to use DNV formula)
alpha=0.2*Hc^2*w0^4/g^2;
g=9.81 [m/s^2]
sigma=0.07  if w<w0, sigma=0.09 otherwise;
S(w)=S1*S2  [m^2 s]  
with,
S1=alpha*g^2*(W^-5)*exp(-(5/4)*(w0/w)^4);
S2=gamma^(exp(-(w-w0)^2/(2*(sigma*w0)^2))); 
Reference [3]


Torsethaugen (Hs,w0):
The Torsethaugen spectrum is an empirical two peaked spectrum for swell 
and developing sea based on experimental data from the North Sea. For 
small peak frequencies, i.e.  0 < wmax <= 0.6 only one peak in the 
spectrum appears. Returns the spectral density function S of the 
Torsethaugen spectrum for the frequencies in the vector  W [rad/s].

p1 = Hs  --significant wave height [m] 
p2 = w0  -- Modal (peak) frequency [rad/s]
S  - vector of power spectral densities [m^2s]
Reference [4]


References: 

[1] A.R.M.J LLoyd (1998) "Seakeeping: Ship Behaviour in Rough Wheather."
    Published by A.R.M.J LLoyd, Gosport UK.ISBN 0-9532634-01

[2] Ochi, M.K. (1998) "Ocean Waves, The stochastic Approach"
    Cambridge Ocean Technology Series, Vol 6,
    Cambridge University Press.

[3] Sørensen, A.J. (2005) "Marine Cybernetics: Modelling and Control"
    Lecture Notes for TMR4241 Marine Control Systems, NTNU, 2005.

[4] K.Torsethaugen (1996): "Model for a Doubly Peaked Wave Spectra"
     Sintef report no.: STF22 A96204 prepared for Norsk Hydro.

[5] T.I. Fossen (2002) "Marine Control Systems" Marine Cybernetics. 

[6] Lewis E.V. "Principles of Naval Architecture volume  III
    Motions in Waves and Controllability." SNAME, 1989.

Created by: Tristan Perez in 2005  
Revisions: 2007-03-09 minor fixes Øyvind Smogeli
            Adapted for MSS V3.0 March 2005. Grouped all spectra in one
            function. Added Spec_type=1,4,5,6.
           2009-09-11 fixed scaling of JONSWAP sepctrum
________________________________________________________________

MSS GNC is a Matlab toolbox for guidance, navigation and control.
The toolbox is part of the Marine Systems Simulator (MSS).

Copyright (C) 2008 Thor I. Fossen and Tristan Perez

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: contact@marinecontrol.org
URL:    <http://www.marinecontrol.org>
'''
import math
import numpy as np
import matplotlib.pyplot as plt

def wave_spectrum(SpecType, Par, W, PlotFlag):
    S = []
    m = len(W)
    '''
    [m,n] = W.shape
    if n > m: 
        print("Error: W must be a column vector")
        return 
    '''
    def bretschneither(Par, W, S):
        A = Par[0]
        B = Par[1]
        for k in range(0, m):
            Sa = A*W[k]**(-5)*np.exp(-B/(W[k]**4))
            S.append(Sa)
        return S

    def piersonmoskowitz(Par, W, S):
        Vwind20 = Par[0]
        A = 0.0081*9.81**2
        B = 0.74*(9.81/Vwind20)**4
        for k in range(0, m):
            Sa = A*W[k]**(-5)*np.exp(-B/(W[k]**4))
            S.append(Sa)
        return S

    def ittcmodifiedpmt0(Par, W, S):
        Hs = Par[0]
        T0 = Par[1]
        A = 487*Hs**2/(T0**4)
        B = 1949/(T0**4)
        for k in range(0, m):
            Sa = A*W[k]**(-5)*np.exp(-B/(W[k]**4))
            S.append(Sa)
        return S

    def ittcmodifiedpmt1(Par, W, S):
        Hs = Par[0]
        T1 = Par[1]
        A = 173*Hs**2/(T1**4)
        B = 691/(T1**4)
        for k in range(0, m):
            Sa = A*W[k]**(-5)*np.exp(-B/(W[k]**4))
            S.append(Sa)
        return S

    def ittcmodifiedpmtz(Par, W, S):
        Hs = Par[0]
        Tz = Par[1]
        A = 123*Hs**2/(Tz**4)
        B = 495/(Tz**4)
        for k in range(0, m):
            Sa = A*W[k]**(-5)*np.exp(-B/(W[k]**4))
            S.append(Sa)
        return S

    def jonswapwind(Par, W, S):
        Vw10 = Par[0]
        fetch = Par[1]
        g = 9.81
        xtilde = g*fetch/(Vw10**2)
        f0 = 3.5*(g/Vw10)*xtilde**(-0.33)
        w0 = 2*np.pi*f0
        alpha = 0.076*xtilde**(-0.22)
        for k in range(0, m):
            if W[k] < w0:
                sigma = 0.07
            else:
                sigma = 0.09
            S1 = alpha*g**2*(W[k]**(-5))*np.exp(-(5/4)*(w0/W[k])**4)
            S2 = 3.3**(np.exp(-(W[k]-w0)**2/(2*(sigma*w0)**2)))
            S.append(S1*S2)
        return S

    def jonswapgamma(Par, W, S):
        Hs = Par[0]
        w0 = Par[1]
        gamma = Par[2]
        g = 9.81
        alpha = 0.2*Hs**2*w0**4/(g**2)
        if gamma < 1 or gamma > 7:
            if gamma != 0:
                print("Warning: gamma value in wave_spectrum function outside validity range, using DNV formula")
            k=2*np.pi/(w0*np.sqrt(Hs))
            if k <= 3.6:
                gamma = 5
            elif k <= 5:
                gamma = np.exp(5.75-1.15*k)
            else:
                gamma = 1
        for k in range(0, m):
            if W[k] < w0:
                sigma = 0.07
            else:
                sigma = 0.09
            S1 = alpha*g**2*(W[k]**(-5))*np.exp(-(5/4)*(w0/W[k])**4)
            S2 = 3.3**(np.exp(-(W[k]-w0)**2/(2*(sigma*w0)**2)))
            Conv_factor = 1-0.287*np.log(gamma)
            S.append(Conv_factor*S1*S2)
        return S

    def torset_spec(Hs, wo, omg):
        '''
        The Torsethaugen spectrum is an empirical two peaked spectrum for swell and 
        developing sea based on experimental data from the North Sea. For small peak
        frequencies, i.e.  0 < wmax <= 0.6 only one peak in the spectrum appears.
        Returns the spectral density function S of the Torsethaugen spectrum for 
        the frequencies:  0 < w < wmax (rad/s).

        Ouputs:
        S     	- vector of power spectral densities (m^2s)

        Inputs:
        Hs    	- significant wave height (m) - mean of the ones third highest waves
        wo    	- peak frequency (rad/s)
        omg  	- vector of frequencies at which to calculate S

        Ref: K.Torsethaugen (1996): "Model for a Doubly Peaked Wave Spectra"
            Sintef report no.: STF22 A96204 prepared for Norsk Hydro.

        Author:     G. Kleiven, Norsk Hydro 
        Date:       2000-06-15
        Revisions:  2001-07-06,Svein I. Sagatun, Norsk Hydro - minor revisions
                    2001-10-14,Thor I. Fossen - IO compatible with the GNC toolbox
                    2005-03-12 Øyvind Smogeli - Revised to comply with MSS, added output consistency test
                    2007-10-08 Øyvind Smogeli - Bug fix for scaling of spectrum magnitude

        ---------------------------------------------------------------------------
        source code: Norsk Hydro
        ---------------------------------------------------------------------------
        '''
        Tp = 2*np.pi/wo
        f2pii = 2*np.pi
        fwtp = f2pii/Tp
        Nfrq = len(omg)

        if Hs > 0:
            af = 6.6
            ae = 2.0
            au = 25
            a10 = 0.7
            a1 = 0.5
            kg = 35.0
            kg0 = 2.5
            kg1 = 1.0
            r = 0.857
            k0 = 0.5
            k00 = 3.2
            m0 = 4.0
            b1 = 2.0
            a20 = 0.6
            a2 = 0.3
            a3 = 6.0
            s0 = 0.08
            s1 = 3.0
            b2 = 0.7
            b3 = 3.0
            sigma_a2 = 2*0.07**2
            sigma_b2 = 2*0.09**2
            tf = af*Hs**(1/3)
            if Tp < tf:
                tl = ae*(Hs**0.5)
                eps1 = (tf-Tp)/(tf-tl)
                rpw = (1-a10)*np.exp(-((eps1/a1)**2))+a10
                hsw = rpw*Hs
                hss = np.sqrt(1-rpw**2)*Hs
                tpw = Tp
                tps = tf+b1
                sp = ((f2pii/9.81)*hsw/(Tp**2))
                gammaw = kg*(1+kg0*np.exp(-Hs/kg1))*(sp**r)
                gammas = 1
                nw = k0*np.sqrt(Hs)+k00
                mw = m0
                ns = nw
                ms = mw
                if ms < 1:
                    ms = 1
                g_argw = (nw-1)/mw
                g_args = (ns-1)/ms
                if g_args < 0: 
                    g0s = 1/((1/ms)*math.gamma(g_args)*((ns/ms)**(-g_args)))
                else:
                    g0s = 1/((1/ms)*math.gamma(g_args)/((ns/ms)**(g_args)))
                if g_argw < 0:
                    g0w = 1/((1/mw)*math.gamma(g_argw)*((nw/mw)**(-g_argw)))
                else:
                    g0w = 1/((1/mw)*math.gamma(g_argw)/((nw/mw)**(g_argw)))
                a1m = 4.1
                b1m = 2.0*(mw**0.28)-5.3
                c1m = -1.45*(mw**0.1)+0.96
                a2m = 2.2/(mw**3.3)+0.57
                b2m = -0.58*mw**0.37+0.53
                c2m = -1.04/(mw**1.9)+0.94
                if c1m < 0:
                    f1w = a1m/((nw-b1m)**(-c1m))
                else:
                    f1w = a1m *(nw-b1m)**c1m
                if b2m < 0:
                    f2w = a2m/(nw**(-b2m))+c2m
                else:
                    f2w = a2m*nw**b2m+c2m
                b1m = 2.0*(ms**0.28)-5.3
                c1m = -1.45*(ms**0.1)+0.96
                a2m = 2.2/(ms**3.3)+0.57
                b2m = -0.58*ms**0.37+0.53
                c2m = -1.04/(ms**1.9)+0.94
                if c1m < 0:
                    f1s = a1m/((ns-b1m)**(-c1m))
                else:
                    f1s = a1m *(ns-b1m)**c1m
                if b2m < 0:
                    f2s = a2m/(ns**(-b2m))+c2m
                else:
                    f2s = a2m*ns**b2m+c2m
                agammaw = (1+f1w*np.log(gammaw)**(f2w))/gammaw
                agammas = (1+f1s*np.log(gammas)**(f2s))/gammas
            else:
                tu = au
                epsu = (Tp-tf)/(tu-tf)
                rps = (1-a20)*np.exp(-(epsu/a2)**2)+a20
                hss = rps*Hs
                hsw = np.sqrt(1-rps**2)*Hs
                tps = Tp
                ns = k0*np.sqrt(Hs)+k00
                ms = m0
                nw = ns
                mw = m0*(1-b2*np.exp(-Hs/b3))
                s4 = s0*(1-np.exp(-Hs/b3))
                g_argw = (nw-1)/mw
                g_args = (ns-1)/ms
                if g_args < 0: 
                    g0s = 1/((1/ms)*math.gamma(g_args)*((ns/ms)**(-g_args)))
                else:
                    g0s = 1/((1/ms)*math.gamma(g_args)/((ns/ms)**(g_args)))
                if g_argw < 0:
                    g0w = 1/((1/mw)*math.gamma(g_argw)*((nw/mw)**(-g_argw)))
                else:
                    g0w = 1/((1/mw)*math.gamma(g_argw)/((nw/mw)**(g_argw)))
                tpw = ((g0w*hsw**2)/(16*s4*(0.4**nw)))**(1/(nw-1))
                sf = ((f2pii/9.81)*Hs/(tf**2))
                gammaw = 1
                gamma_f = kg*(1+kg0*np.exp(-Hs/kg1))*sf**r
                gammas = gamma_f*(1+a3*epsu)
                a1m = 4.1
                b1m = 2.0*(mw**0.28)-5.3
                c1m = -1.45*(mw**0.1)+0.96
                a2m = 2.2/(mw**3.3)+0.57
                b2m = -0.58*mw**0.37+0.53
                c2m = -1.04/(mw**1.9)+0.94
                if c1m < 0:
                    f1w = a1m/((nw-b1m)**(-c1m))
                else:
                    f1w = a1m *(nw-b1m)**c1m
                if b2m < 0:
                    f2w = a2m/(nw**(-b2m))+c2m
                else:
                    f2w = a2m*nw**b2m+c2m
                b1m = 2.0*(ms**0.28)-5.3
                c1m = -1.45*(ms**0.1)+0.96
                a2m = 2.2/(ms**3.3)+0.57
                b2m = -0.58*ms**0.37+0.53
                c2m = -1.04/(ms**1.9)+0.94
                if c1m < 0:
                    f1s = a1m/((ns-b1m)**(-c1m))
                else:
                    f1s = a1m *(ns-b1m)**c1m
                if b2m < 0:
                    f2s = a2m/(ns**(-b2m))+c2m
                else:
                    f2s = a2m*ns**b2m+c2m
                agammaw = (1+f1w*np.log(gammaw)**(f2w))/gammaw
                agammas = (1+f1s*np.log(gammas)**(f2s))/gammas
            fdenorm_s = (tps*(hss**2))/16
            fdenorm_w = (tpw*(hsw**2))/16

            f = omg/f2pii
            fnw = f*tpw
            ind = max([i for i,v in enumerate(fnw) if v < 1])
            ftest1 = []
            ftest1[0:ind] = np.exp(-(((fnw[0:ind]-1)**2)/sigma_a2))
            ftest1[ind:Nfrq] = np.exp(-(((fnw[ind:Nfrq]-1)**2)/sigma_b2))
            #gamma_wf = gammaw**ftest1
            gamma_wf = list(map(lambda x:pow(gammaw,x),ftest1))
            gamma_ws_1 = fnw**(-nw)
            gamma_ws_2 = np.exp(-(nw/mw)*fnw**(-mw))
            gamma_ws = gamma_ws_1*gamma_ws_2
            sw = g0w*agammaw*gamma_ws*gamma_wf*fdenorm_w/(2*np.pi)

            fns = f*tps
            ins = max([i for i,v in enumerate(fns) if v < 1])
            ftest2 = []
            ftest2[0:ins] = np.exp(-(((fns[0:ins]-1)**2)/sigma_a2))
            ftest2[ins:Nfrq] = np.exp(-(((fns[ins:Nfrq]-1)**2)/sigma_b2))
            gamma_sf = list(map(lambda x:pow(gammas,x),ftest2))
            #gamma_sf = gammas**ftest2
            gamma_ss_1 = fns**(-ns)
            gamma_ss_2 = np.exp(-(ns/ms)*fns**(-ms))
            gamma_ss = gamma_ss_1*gamma_ss_2
            ss = g0s*agammas*gamma_ss*gamma_sf*fdenorm_s/(2*np.pi)

            S = sw + ss
        else:
            S = np.zeros(len(omg)).tolist()
        
        if sum(S.imag) != 0:
            S = np.zeros(len(omg)).tolist()
            print("Torsethaugen spectrum input outside validity range in wave_spectrum, complex output set to zero")



    def torsethaugen(Par, W, S):
        Hs = Par[0]
        w0 = Par[1]
        N = len(W)
        wmax = max(W)
        S = torset_spec(Hs,w0,W)
        return S

   # options = {1 : bretschneither(Par, W, S),
   #             2 : piersonmoskowitz(Par, W, S),
    #            3 : ittcmodifiedpmt0(Par, W, S),
     #           4 : ittcmodifiedpmt1(Par, W, S),
      #          5 : ittcmodifiedpmtz(Par, W, S),
       #         6 : jonswapwind(Par, W, S),
        #        8 : torsethaugen(Par, W, S)}
   # S = options[SpecType]
    if SpecType == 1:
        return bretschneither(Par, W, S)
    if SpecType == 2:
        return piersonmoskowitz(Par, W, S)
    if SpecType == 3:
        return ittcmodifiedpmt0(Par, W, S)
    if SpecType == 4:
        return ittcmodifiedpmt1(Par, W, S)
    if SpecType == 5:
        return ittcmodifiedpmtz(Par, W, S)
    if SpecType == 6:
        return jonswapwind(Par, W, S)
    if SpecType == 7:
        return jonswapgamma(Par, W, S)
    if SpecType == 8:
        return torsethaugen(Par, W, S)
    
