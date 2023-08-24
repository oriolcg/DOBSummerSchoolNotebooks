import numpy as np
from scipy import optimize
from typing import List, Tuple
import copy

#----- Start Constants -----
g = 9.81
rad2deg = 90 / np.pi
H0Start = 10.0
LBedStart = 10.0
#----- End Constants -----


#----- Start Functions -----

def catenary_xz(W_pm, EA, V0, H0, L_susp):
    
    W = W_pm*L_susp*g
        
    s = L_susp    
    z_t1 = 1 + (V0/H0)**2
    z_t2 = 1 + ((V0 - W*s/L_susp)/H0)**2
    z = W*s/EA * (s/2/L_susp - V0/W) - H0*L_susp/W*(np.sqrt(z_t1) - np.sqrt(z_t2))            
    x_susp = H0*s/EA + H0*L_susp/W*(np.arcsinh(V0/H0) - np.arcsinh((V0 - W*s/L_susp)/H0))
    
    return x_susp, z

def catenary_end_xz(W_pm, L_tot, EA, V0, H0, L_susp):
    
    W = W_pm*L_susp*g
        
    s = L_susp    
    z_t1 = 1 + (V0/H0)**2
    z_t2 = 1 + ((V0 - W*s/L_susp)/H0)**2
    z = W*s/EA * (s/2/L_susp - V0/W) - H0*L_susp/W*(np.sqrt(z_t1) - np.sqrt(z_t2))            
    x_susp = H0*s/EA + H0*L_susp/W*(np.arcsinh(V0/H0) - np.arcsinh((V0 - W*s/L_susp)/H0))
    x_bed = L_tot - L_susp
    x = x_bed + x_susp
        
    return x, z

#----- End Functions -----


#----- Start moorSeg -----
class moorSeg:
    
    # Prevent adding attributes outside of __init__
    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "%r is a frozen class" % self )
        object.__setattr__(self, key, value)

    
    def __init__(self, LTot=100, Wpm=376, EA=1861000000):
        self.EA = EA
        self.Wpm = Wpm
        self.LTot = abs(LTot)        
        self.__LBed = self.LTot
        self.__LSusp = 0.0        
        self.xe = self.LTot
        self.ze = 0.0

        self.__isfrozen = True

    
    def end_xz(self, V0, H0, LBed):
        LSusp = self.LTot - LBed        
        x, z = catenary_end_xz( self.Wpm, self.LTot, self.EA,
            V0, H0, LSusp )        
        return x,z

    
    def setLSusp(self, LSusp):
        self.__LSusp = LSusp
        self.__LBed = self.LTot - LSusp

    
    def setLBed(self, LBed):
        self.__LBed = LBed
        self.__LSusp = self.LTot - LBed

    
    def retLBed(self):
        return self.__LBed

    
    def retLSusp(self):
        return self.__LSusp


    def retLTot(self):
        return self.LTot

#----- End moorSeg -----


#----- Start moorLine2D -----
class moorLine2D:

    # Prevent adding attributes outside of __init__
    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "%r is a frozen class" % self )
        object.__setattr__(self, key, value)


    def __init__(self, seg=[moorSeg()]):

        self.nSeg = len(seg)
        self.seg = copy.deepcopy(seg)
        
        # Calculated        
        self.xe = 0.0
        self.ze = 0.0
        self.V0 = 0.0
        self.Ve = 0.0
        self.H0 = 0.0

        self.__isfrozen = True


    def retLBed(self):
        return sum( [iSeg.retLBed() for iSeg in self.seg] )

    def retLSusp(self):
        return sum( [iSeg.retLSusp() for iSeg in self.seg] )

    def retLTot(self):
        return sum( [iSeg.LTot for iSeg in self.seg] )


    def setEnd_xz(self, V0, H0, LBed):        

        self.V0 = V0
        self.H0 = H0

        if(LBed < 0.0):
            LBed = 0.0

        V0r = V0

        LBedRun = LBed
        for iSeg in self.seg:            
            if(LBedRun >= iSeg.LTot):
                iSeg.setLSusp(0.0)                
                LBedRun -= iSeg.retLBed()
                iSeg.xe = iSeg.retLBed()
                iSeg.ze = 0.0
                # print(iSeg.xe, iSeg.ze)
            
            else:                
                iSeg.setLBed(LBedRun)                
                LBedRun=0.0 #disabling it

                iSeg.xe, iSeg.ze = iSeg.end_xz(V0r, H0, iSeg.retLBed())
                V0r -= iSeg.Wpm * iSeg.retLSusp() * g                
                # print(iSeg.xe, iSeg.ze)

        self.Ve = V0r
        self.xe = sum( [iSeg.xe for iSeg in self.seg] )
        self.ze = sum( [iSeg.ze for iSeg in self.seg] )
        # print(self.xe, self.ze)
        # print()

        return        


    def objFnc(sol, *data):
        H0 = sol[0]
        LBed = sol[1]
        V0, xtarg, ztarg, lmoor  = data
        lmoor.setEnd_xz(V0, H0, LBed)

        return np.array( [abs(ztarg-lmoor.ze),  abs(xtarg-lmoor.xe)] )

    
    def solveEnd(self, V0, xtarg, ztarg,
        H0Start=H0Start, LBedStart=LBedStart):        

        data = V0, xtarg, ztarg, self
        root = optimize.fsolve(moorLine2D.objFnc, 
            np.array([H0Start, LBedStart]),
            args=data)

        self.setEnd_xz(V0, root[0], root[1])


    def plotLine2D(self, ds):
        LTot = self.retLTot()
        LBed = self.retLBed()
        LSusp  = self.retLSusp()

        # Taking only suspended portion
        # Excluding last point
        s = np.arange(LBed, LTot, ds)
        x = 0.0 * s
        z = 0.0 * s
        segTyp = [0 for si in s]
        x[0] = LBed

        segN = 0
        iSeg = self.seg[segN]
        Lr = iSeg.retLTot()
        V0r = self.V0
        H0 = self.H0

        # Excluding first point
        # First point is already set 
        for i in range(1, len(s)):
            si = s[i]
            while si > Lr:                
                segN += 1
                iSeg = self.seg[segN]
                Lr += iSeg.retLTot()

            dx, dz = catenary_xz( iSeg.Wpm, iSeg.EA,
                V0r, H0, ds )        
            x[i] = x[i-1] + dx
            z[i] = z[i-1] + dz
            segTyp[i] = segN
            V0r -= iSeg.Wpm * ds * g

        x = np.hstack([0.0, x, self.xe])
        z = np.hstack([0.0, z, self.ze])
        segTyp = [0] + segTyp + [self.nSeg-1]

        return x, z, segTyp

#----- End moorLine2D -----


#----- Start moorLine3D -----
class moorLine3D:

    # Prevent adding attributes outside of __init__
    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "%r is a frozen class" % self )
        object.__setattr__(self, key, value)


    def __init__(self, xan, yan, zan=0.0, seg=[moorSeg()]):
        self.xan = xan
        self.yan = yan
        self.zan = zan
        self.l2d = moorLine2D( seg )

        self.xflRef = 0.0
        self.yflRef = 0.0
        self.zflRef = 0.0
        self.H0Start = H0Start
        self.LBedStart = LBedStart

        self.xfl = 0.0
        self.yfl = 0.0
        self.zfl = 0.0
        self.nx = 1.0
        self.ny = 0.0

        self.__isfrozen = True

    
    def solveFl_xyz(self, V0, xfl, yfl, zfl, 
        H0Start=H0Start, LBedStart=LBedStart):

        self.xfl = xfl
        self.yfl = yfl
        self.zfl = zfl        

        xtarg = np.sqrt( (self.xfl - self.xan)**2 + 
            (self.yfl - self.yan)**2 )
        ztarg = self.zfl - self.zan

        self.nx = (self.xfl - self.xan) / xtarg
        self.ny = (self.yfl - self.yan) / xtarg

        self.l2d.solveEnd( V0, xtarg, ztarg, 
            H0Start=H0Start, LBedStart=LBedStart)    

    
    def retFlForce(self):
        Fx = -self.l2d.H0 * self.nx
        Fy = -self.l2d.H0 * self.ny
        Fz = self.l2d.Ve

        return Fx, Fy, Fz


    def retLBed(self):
        return self.l2d.retLBed()

    def retLSusp(self):
        return self.l2d.retLSusp()

    def retLTot(self):
        return self.l2d.retLTot()


    def setFlRef_xyz(self, xflRef, yflRef, zflRef):
        self.xflRef = xflRef
        self.yflRef = yflRef
        self.zflRef = zflRef


    def plotLine3D(self, ds):

        x2d, z2d, segTyp = self.l2d.plotLine2D(ds)

        x = [self.xan + xi*self.nx for xi in x2d]
        y = [self.yan + xi*self.ny for xi in x2d]
        z = [self.zan + zi for zi in z2d]

        return x, y, z, segTyp


#----- End moorLine3D -----