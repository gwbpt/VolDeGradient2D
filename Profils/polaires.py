#!/usr/bin/env python

# Polaire type 2 Re(Cl) = Re0 * sqrt(Cl)
 
from __future__ import print_function, division
from collections import OrderedDict

from math import sqrt, log10

import numpy as np
from scipy import interpolate, optimize

#POLAR_DIR = '../Profils/Polaires' # relatif a l'appelant !

def debug(msg):
    print("Debug : " + msg)
    quit()

g  =  9.81   # m/s2 earth gravity
ro =  1.225  # air kg/m3 
mu = 15.6e-6 # air viscosity m2/s

two_g_ro  = 2.0 * g / ro
ro_mu = ro / mu

def speed2Cl(speed, wingLoad): # kg/m2
    # mass * g = lift = 0.5 * Cl * ro * Area * speed**2
    # Cl = 2.0*g/ro * (mass/Area) / speed**2
    Cl = two_g_ro * wingLoad / speed**2
    return Cl # Cz

def Cl2Speed(Cl, wingLoad): # kg/m2
    # mass * g = lift = 0.5 * Cl * ro * Area * speed**2
    # speed = sqrt( g * (mass/Area) /(0.5 * ro * Cl))
    speed = sqrt( two_g_ro * wingLoad / Cl )
    return speed # m/s

def speed2Reynolds(speed, wingChord): # m/s, m
    reynolds = speed * wingChord * ro_mu
    return reynolds
    
def speedReynolds(Cl, wingLoad=5.0, wingChord=0.20): # 5.0kg/m2 = 50gr/dm2
    speed = Cl2Speed(Cl, wingLoad)
    reynolds = speed2Reynolds(speed, wingChord=wingChord)
    #print("wingLoad :%5.1f kg/m2, wingChord :%4.2f m, Cl :%4.2f => %5.1f m/s(%3.0fkm/h), Re:%6.0fk"%(wingLoad, wingChord, Cl, speed, speed*3.6, reynolds/1000.0))
    return speed, reynolds
    
if 0:
    # rc light glider
    wingLoad  = 2.6
    wingChord = 0.20
    speed     = 6.5 # m/s
    Cl        = speed2Cl(speed, wingLoad)
    speed, reynolds = speedReynolds(Cl, wingLoad=wingLoad, wingChord=wingChord)
    
    # model air graft
    wingLoad  = 4.0
    wingChord = 0.15
    speed     = 9.0 # m/s
    Cl        = speed2Cl(speed, wingLoad)
    speed, reynolds = speedReynolds(Cl, wingLoad=wingLoad, wingChord=wingChord)
    
    # small aircraft
    wingLoad  = 50.0
    wingChord = 1.5
    speed     = 50.0 # m/s
    Cl        = speed2Cl(speed, wingLoad)
    speed, reynolds = speedReynolds(Cl, wingLoad=wingLoad, wingChord=wingChord)
    
    # line aircraft
    wingLoad  = 500.0
    wingChord = 3.5
    speed     = 230.0 # m/s
    Cl        = speed2Cl(speed, wingLoad)
    speed, reynolds = speedReynolds(Cl, wingLoad=wingLoad, wingChord=wingChord)
    
    quit()

def addToDict(dic, alpha, reId, datas, nbRe=5):
    if alpha in dic :
        lst = dic[alpha]
    else :
        lst = [None] * nbRe
        dic[alpha] = lst
    lst[reId] = datas
        
def readPolarsFiles(fullFileNameWoReynolds, reynolds=(50e3,100e3,200e3,500e3,1000e3)):
    #log10Res = np.array([log10(Re) for Re in reynolds]) # usefull for interpolation
    tables = list()
    alphaGlobalRange = [+100.0, -100.0]
    alphaCommonRange = [-100.0, +100.0]
    allAlphas = dict()
    nbRe = len(reynolds)
    for reId, reyn in enumerate(reynolds):
        fileName = "%s-%.0f.txt"%(fullFileNameWoReynolds, reyn)
        print("------------------ read file : %s"%fileName)
        table = None
        with open(fileName, 'r') as f:
            for line in f:
                line = line.strip()
                if table is None :
                    if line[:6] == "------":
                        table = list()
                        print()
                    else:
                        print(line)
                else:
                    #print(line)
                    fields = line.split()
                    floats = [float(e) for e in fields]
                    alpha = floats[0]
                    addToDict(allAlphas, alpha, reId, floats[1:], nbRe=nbRe)
                    table.append(floats)
            almin, almax = table[0][0], table[-1][0]
            print("%3d lines read from %6.3f to %6.3f deg"%(len(table), almin, almax))
            
            if almin < alphaGlobalRange[0] : alphaGlobalRange[0] = almin 
            if almax > alphaGlobalRange[1] : alphaGlobalRange[1] = almax
            
            if almin > alphaCommonRange[0] : alphaCommonRange[0] = almin 
            if almax < alphaCommonRange[1] : alphaCommonRange[1] = almax
            
            tables.append((reyn, np.array(table)))
    print("alphaGlobalRange: ", alphaGlobalRange)
    print("alphaCommonRange: ", alphaCommonRange)
    
    npRe    = np.array(reynolds)
    log10Re = np.log10(npRe)
    
    # inter-extrapolate1d over alphaGlobalRange
    minal, maxal = alphaGlobalRange
    globalAlphas = np.arange(minal, maxal+1e-6, 0.5)
    nGlobAl = globalAlphas.size
    print("nGlobAl :", nGlobAl)
    
    fctCls = list()
    fctCds = list()
    fctCms = list()
    
    shape = (nGlobAl, nbRe)   
    globCls = np.zeros(shape)
    globCds = np.zeros(shape)
    globCms = np.zeros(shape)
    
    for iRe, (reyn, table) in enumerate(tables):            
        alphas = table[:, 0]
        Cls    = table[:, 1]
        Cds    = table[:, 2] 
        Cms    = table[:, 4]
        
        fctCl = interpolate.interp1d(alphas, Cls, fill_value='extrapolate')
        fctCd = interpolate.interp1d(alphas, Cds, fill_value='extrapolate')
        fctCm = interpolate.interp1d(alphas, Cms, fill_value='extrapolate')
        
        fctCls.append(fctCl)
        fctCds.append(fctCd)
        fctCms.append(fctCm)
        
        globCls[:, iRe] = fctCl(globalAlphas)
        globCds[:, iRe] = fctCd(globalAlphas)
        globCms[:, iRe] = fctCm(globalAlphas)
        
    # interpolate2d
    globClFct = interpolate.interp2d(log10Re, globalAlphas, globCls)
    globCdFct = interpolate.interp2d(log10Re, globalAlphas, globCds)
    globCmFct = interpolate.interp2d(log10Re, globalAlphas, globCms)
    
    globalTable = log10Re, globalAlphas, globCls, globCds, globCms, globClFct, globCdFct, globCmFct
    
    return tables, globalTable

class Profils(object):
    def __init__(self, profilsDir=""):
        self.profilsDir = profilsDir
        self.profils    = OrderedDict()
        
    def add(self, name="xf-n0009sm-il"):
        p = Profil(name, self.profilsDir)
        self.profils[name] = p
        return p
        
class Profil(object):
    def __init__(self, name, profilsDir):
        print("Profil('%s')"%name)
        self.name      = name
        self.synthetic = False
        
        if   name == "ApproximNACA0009":
            self.synthetic = True
            
            self.ClMax  = 0.9
            self.ClMin  =-0.9
            self.Cl0    = 0.0
            self.dCl_da = 1.8/16
            
            self.Cd0        = 0.012
            self.alpha0     = 0.0 # deg
            self.Cd_alpha4  = (0.10-self.Cd0) / 9.0**4
            self.alphaRange = -9, 9 # deg integer
            
            self.Cm0   = 0.0
            
            reynolds = 100e3
            self.tables = [reynolds, np.array(self.polarGen(reynolds))]
            
        elif name == "ApproximRG14":
            self.synthetic = True
            
            self.ClMax  = 1.0
            self.ClMin  =-0.5
            self.Cl0    = 0.015
            self.dCl_da = 1.5/15
            
            self.Cd0        = 0.015
            self.alpha0     = 2.5 # deg
            self.Cd_alpha4  = (0.10-self.Cd0) / 12.0**4
            self.alphaRange = -9, 14
            
            self.Cm0   = 0.04
            
            reynolds = 100e3
            self.tables = [reynolds, np.array(self.polarGen(reynolds))]
            
        else:
            fileName = profilsDir + '/' + name
            self.tables, globalTable = readPolarsFiles(fileName)
            self.log10Re, self.globAls, self.globCls, self.globCds, self.globCms, self.globClFct, self.globCdFct, self.globCmFct = globalTable
            
    def ClCdCmSynth(self, alpha, reynolds=100.000):
        Cl = self.Cl0 + self.dCl_da * alpha # linear approximation
        if   Cl > self.ClMax : 
            Cl = self.ClMax
        elif Cl < self.ClMin :
            Cl = self.ClMin
        Cd = self.Cd0 + self.Cd_alpha4 * (alpha-self.alpha0)**4 # 
        return Cl, Cd, self.Cm0

    def polarGen(self, reynolds=100e3):
        aMin, aMax = self.alphaRange
        n = aMax - aMin + 1
        table = list()
        for i in range(n):
            alpha = aMin + i
            Cl, Cd, Cm = self.ClCdCmSynth(alpha, reynolds)
            print("%2d:%6.1f,%6.2f,%6.3f"%(i, alpha, Cl, Cd))
            table.append([alpha, Cl, Cd, 0.0, Cm, 0.0, 0.0])      
        return table
        
    def ClCdCm(self, alpha, Re):
        if self.synthetic : return self.ClCdCmSynth(alpha, Re)
        
        logRe = log10(Re)
        Cl = self.globClFct(logRe, alpha)
        Cd = self.globCdFct(logRe, alpha)
        Cm = self.globCmFct(logRe, alpha)
        return Cl, Cd, Cm  
        
    def polarAtCteLift(self, wingLoad=2.5, wingChord=0.20): # ie 25gr/dm2, 0.2m
        # a given Cl give a Speed such as Cl * V**2 = Cte 
        # which give a Renolds Re0 * V which give a Cd
        
        # mass * g = lift = 0.5 * Cl * ro * Area * speed**2
        # speed = sqrt( g * (mass/Area) /(0.5 * ro * Cl))
        print("polarAtCteLift(wingLoad=%5.2fkg/m2, wingChord=%5.2f m)"%(wingLoad, wingChord))
        n = 10
        polarCteLift = np.zeros((n,4))
        for i in range(n):
            Cl = 0.1 * (i+1)
            speed, reynolds = speedReynolds(Cl, wingLoad=wingLoad, wingChord=wingChord)
            logRe = log10(reynolds)
            #print("Cl:%8.3f, %8.2f m/s, Re=%8.0f"%(Cl, speed, reynolds))
            
            #ClsAtRe = self.globClFct([reynolds], self.globAls)
            #print("ClsAtRe.size :", ClsAtRe.size)
            
            def ClFctAtRe(alpha) : return self.globClFct(logRe, alpha) - Cl
            #print("type(ClFctAtRe):", type(ClFctAtRe))
            res = optimize.fsolve(ClFctAtRe, 0.0)
            alpha0 = res[0]
            Cl0    = self.globClFct(logRe, alpha0)
            Cd0    = self.globCdFct(logRe, alpha0)
            Cm0    = self.globCmFct(logRe, alpha0)
            print("%8.2f m/s, %8.3f deg, %8.0f Re : Cl:%6.2f, Cd:%7.4f, Cm:%8.3f, Cl/Cd:%6.1f, desc:%6.3f m/s, pow:%6.3f"%(speed, alpha0, reynolds, Cl0, Cd0, Cm0, Cl0/Cd0, speed*Cd0/Cl0, Cd0*speed**2))
            polarCteLift[i, 0] = alpha0
            polarCteLift[i, 1] = Cl0
            polarCteLift[i, 2] = Cd0
            polarCteLift[i, 3] = Cm0
            
        return polarCteLift    
          
    def polarsVsWingload(self, wingChord=0.20): # 0.2m
        self.polarWingload = list()
        for wingLoad in (2.0, 5.0, 10.0):
            polar = self.polarAtCteLift(wingLoad=wingLoad, wingChord=wingChord)
            self.polarWingload.append((wingLoad, polar))
    
    def plotExtrapolatedPolar(self):     
        print("plotExtrapolatedPolar %s"%(self.name))
        
        Reys = 10.0**self.log10Re

        alRange = -10.0,14.001
        cdRange =   0.0, 0.0801
        clRange =  -0.7, 1.3001
        cl_cd_Range = -50.0, 100.1
        
        alpha_ticks = np.arange(alRange[0], alRange[1], 2.0)
        Cd_ticks    = np.arange(0.0, cdRange[1], 0.01)
        Cl_ticks    = np.arange(clRange[0], clRange[1], 0.1)
        Cl_Cd_ticks = np.arange(-50.0, 100.0, 10.0) # finesse
        
        fig2 = plt.figure(figsize=(8, 8))
        fig2.suptitle('Extrapolated Polar of %s'%self.name, fontsize=16)
        
        axClvCd = fig2.add_subplot(321)
        for i, reyn in enumerate(Reys):
            axClvCd.plot(self.globCds[:,i], self.globCls[:,i], '.-', label="re=%3.0fK"%(reyn*1e-3))
        axClvCd.set_xlim(cdRange)
        axClvCd.set_ylim(clRange)
        axClvCd.set_xticks(Cd_ticks)                                                       
        axClvCd.set_yticks(Cl_ticks)                                                       
        axClvCd.legend()       
        plt.grid()
        plt.title("Cl v Cd")

        axClvCdPct = fig2.add_subplot(323)
        for wingload, polar in self.polarWingload:
            axClvCdPct.plot(polar[:,2], polar[:,1], '.-', label="wingload=%4.1f kg/m2"%(wingload))
        axClvCdPct.set_xlim(cdRange)
        axClvCdPct.set_ylim(clRange)
        axClvCdPct.set_xticks(Cd_ticks)                                                       
        axClvCdPct.set_yticks(Cl_ticks)                                                       
        axClvCdPct.legend()
        plt.grid()
        plt.title("Cl v Cd at Cte Portance")

        axClvAl = fig2.add_subplot(322)
        for i, reyn in enumerate(Reys):
            axClvAl.plot(self.globAls, self.globCls[:,i], '.-', label="re=%3.0fK"%(reyn*1e-3))
        axClvAl.set_xlim(alRange)
        axClvAl.set_ylim(clRange)
        axClvAl.set_xticks(alpha_ticks)                                                       
        axClvAl.set_yticks(Cl_ticks)                                                       
        axClvAl.legend()
        plt.grid()
        plt.title("Cl v alpha")

        axCl_CdvAl = fig2.add_subplot(326)
        for i, reyn in enumerate(Reys):
            Cl_Cds = self.globCls[:,i] / self.globCds[:,i]
            axCl_CdvAl.plot(self.globAls, Cl_Cds, '.-', label="re=%3.0fK"%(reyn*1e-3))
        axCl_CdvAl.set_xlim(alRange)
        axCl_CdvAl.set_ylim(cl_cd_Range)
        axCl_CdvAl.set_xticks(alpha_ticks)                                                       
        axCl_CdvAl.set_yticks(Cl_Cd_ticks)                                                       
        axCl_CdvAl.legend()
        plt.grid()
        plt.title("Cl/Cd vs alpha")

        axCdvAl = fig2.add_subplot(324)
        for i, reyn in enumerate(Reys):
            axCdvAl.plot(self.globAls, self.globCds[:,i], '.-', label="re=%3.0fK"%(reyn*1e-3))
        axCdvAl.set_xlim(alRange)
        axCdvAl.set_ylim(cdRange)
        axCdvAl.set_xticks(alpha_ticks)                                                       
        axCdvAl.set_yticks(Cd_ticks)                                                       
        axCdvAl.legend()
        plt.grid()
        plt.title("Cd v alpha")

        #plt.show()

    def plotPolaire(self):     
        print("polaire %s"%(self.name))
        
        al_min, al_max  = +100.0, -100.0
        cl_min, cl_max =   +10.0,  -10.0
        cm_min, cm_max =  +100.0, -100.0
        cd_max = 0.0
        cl_cd_min, cl_cd_max =  +1000.0,  -1000.0
        
        for reyn, table in self.tables:            
            if table is not None:
                alphas = table[:, 0]
                a_min, a_max = alphas[0], alphas[-1]
                if a_min < al_min : al_min = a_min           
                if a_max > al_max : al_max = a_max
                  
                Cls    = table[:, 1]
                clmin = np.amin(Cls)
                if clmin < cl_min : cl_min = clmin  
                clmax = np.amax(Cls)
                if clmax > cl_max : cl_max = clmax  
                
                Cms    = table[:, 4]
                cmmin = np.amin(Cms)
                if cmmin < cm_min : cm_min = cmmin  
                cmmax = np.amax(Cms)
                if cmmax > cm_max : cm_max = cmmax  
                
                Cds    = table[:, 2]
                cdmax = np.amax(Cds)
                if cdmax > cd_max : cd_max = cdmax
                
                Cl_Cds = Cls / Cds
                clcdmin = np.amin(Cl_Cds)
                if clcdmin < cl_cd_min : cl_cd_min = clcdmin  
                clcdmax = np.amax(Cl_Cds)
                if clcdmax > cl_max : cl_cd_max = clcdmax  
                
        fmt = "alpha range :%6.2f ..%6.2f ; CL range :%8.4f ..%8.4f ; CM range :%8.4f ..%8.4f ; CD max :%6.5f"
        print(fmt%(al_min, al_max, cl_min, cl_max, cm_min, cm_max, cd_max)) 
           
        al_min_r = np.around(al_min - 0.8, 0)
        al_max_r = np.around(al_max + 0.8, 0)
        alRange  = al_min_r, al_max_r
        
        cd_max_r = np.around(cd_max + 0.008, 2)
        cdRange = 0.0, cd_max
        
        cl_min_r = np.around(cl_min - 0.08, 1)
        cl_max_r = np.around(cl_max + 0.08, 1)
        clRange = cl_min_r, cl_max_r
              
        cm_min_r = np.around(cm_min - 0.08, 1)
        cm_max_r = np.around(cm_max + 0.08, 1)
        cmRange  = cm_min_r, cm_max_r
              
        cl_cd_min_r = np.around(cl_cd_min - 6.0, -1)
        cl_cd_max_r = np.around(cl_cd_max + 6.0, -1)
        cl_cd_Range = cl_cd_min_r, cl_cd_max_r
              
        amint = ((al_min_r-1)//2)*2.0
        amaxt = ((al_max_r+1)//2)*2.0 + 0.01
        
        alpha_ticks = np.arange(amint, amaxt, 2.0)
        Cd_ticks    = np.arange(0, cd_max_r, 0.01)
        Cl_ticks    = np.arange(cl_min_r, cl_max_r+0.01, 0.1)
        Cm_ticks    = np.arange(cm_min_r, cm_max_r+0.01, 0.1)
        Cl_Cd_ticks = np.arange(cl_cd_min_r, cl_cd_max_r, 10.0) # finesse

        print(fmt%(al_min_r, al_max_r, cl_min_r, cl_max_r, cm_min_r, cm_max_r, cd_max_r)) 

        if 0:                                               
            print("alpha_ticks :", alpha_ticks)                                          
            print("Cd_ticks    :", Cd_ticks   )                                          
            print("Cl_ticks    :", Cl_ticks   )                                             
            print("Cl_Cd_ticks    :", Cl_Cd_ticks   )                                             

        fig1 = plt.figure(figsize=(8, 8))
        fig1.suptitle('Polar graphs of %s'%self.name, fontsize=16)
        
        axClvCd = fig1.add_subplot(321)
        for reyn, table in self.tables:
            Cls    = table[:, 1]
            Cds    = table[:, 2] 
            axClvCd.plot(Cds, Cls, '.-', label="re=%6.0f"%reyn)
        axClvCd.legend() # shadow=True, fancybox=False, loc="upper right", bbox_to_anchor=[0, 1], ncol=1, title="Legend"
        #axClvCd.get_legend().get_title().set_color("red")
        
        axClvCd.set_xlim(cdRange)
        axClvCd.set_ylim(clRange)
        axClvCd.set_xticks(Cd_ticks)                                                       
        #axClvCd.set_xticks(minor_ticks, minor=True)                                           
        axClvCd.set_yticks(Cl_ticks)                                                       
        #axClvCd.set_yticks(minor_ticks, minor=True)                                        
        plt.grid()
        plt.title("Cl v Cd")

        axClvCdPct = fig1.add_subplot(323)
        for reyn, table in self.tables:
            Cls    = table[:, 1]
            Cds    = table[:, 2] 
            axClvCdPct.plot(Cds, Cls, '.-', label="re=%6.0f"%reyn)
        axClvCdPct.legend()
        
        axClvCdPct.set_xlim(cdRange)
        axClvCdPct.set_ylim(clRange)
        axClvCdPct.set_xticks(Cd_ticks)                                                       
        #axClvCdPct.set_xticks(minor_ticks, minor=True)                                           
        axClvCdPct.set_yticks(Cl_ticks)                                                       
        #axClvCdPct.set_yticks(minor_ticks, minor=True)                                        
        plt.grid()
        plt.title("Cl v Cd at Cte Portance")

        axClvAl = fig1.add_subplot(322)
        for reyn, table in self.tables:
            alphas = table[:, 0]
            Cls    = table[:, 1]
            axClvAl.plot(alphas, Cls, '.-', label="re=%6.0f"%reyn)
        axClvAl.legend()
        axClvAl.set_xlim(alRange)
        axClvAl.set_ylim(clRange)
        axClvAl.set_xticks(alpha_ticks)                                                       
        axClvAl.set_yticks(Cl_ticks)                                                       
        plt.grid()
        plt.title("Cl v alpha")

        axCl_CdvAl = fig1.add_subplot(326)
        for reyn, table in self.tables:
            alphas = table[:, 0]
            Cls    = table[:, 1]
            Cds    = table[:, 2] 
            Cl_Cds = Cls / Cds
            axCl_CdvAl.plot(alphas, Cl_Cds, '.-', label="re=%6.0f"%reyn)
        axCl_CdvAl.legend()
        axCl_CdvAl.set_xlim(alRange)
        axCl_CdvAl.set_ylim(cl_cd_Range)
        axCl_CdvAl.set_xticks(alpha_ticks)                                                       
        axCl_CdvAl.set_yticks(Cl_Cd_ticks)                                                       
        plt.grid()
        plt.title("Cl/Cd vs alpha")

        axCdvAl = fig1.add_subplot(324)
        for reyn, table in self.tables:
            alphas = table[:, 0]
            Cds    = table[:, 2] 
            axCdvAl.plot(alphas, Cds, '.-', label="re=%6.0f"%reyn)
        axCdvAl.legend()
        axCdvAl.set_xlim(alRange)
        axCdvAl.set_ylim(cdRange)
        axCdvAl.set_xticks(alpha_ticks)                                                       
        axCdvAl.set_yticks(Cd_ticks)                                                       
        plt.grid()
        plt.title("Cd v alpha")

        axCmvAl = fig1.add_subplot(325)
        for reyn, table in self.tables:
            alphas = table[:, 0]
            Cms    = table[:, 4]
            axCmvAl.plot(alphas, Cms, '.-', label="re=%6.0f"%reyn)
        axCmvAl.legend()
        axCmvAl.set_xlim(alRange)
        axCmvAl.set_xticks(alpha_ticks)                                                       
        if 0:
            axCmvAl.set_ylim(cmRange)
            axCmvAl.set_yticks(Cm_ticks)                                                       
        plt.grid()
        plt.title("Cm v alpha")

        #plt.show()

#--------------------------------
if __name__ == "__main__":

    POLAR_DIR = '../Profils/Polaires'
    
    profils = Profils(profilsDir=POLAR_DIR)
    p = profils.add("xf-rg14-il") # ApproximNACA0009 ApproximRG14 xf-n0009sm-il xf-rg14-il xf-rg15-il
    
    Re = 100e3
    for i in range(-10, 11):
        alpha = 1.0 * i
        Cl, Cd, Cm = p.ClCdCm(alpha, Re)
        print("%6.1f deg, Cl:%6.2f, Cd%6.3f, Cm%6.3f"%(alpha, Cl, Cd, Cm))
    
    import matplotlib.pyplot as plt
    
    p.plotPolaire()
    p.polarsVsWingload(wingChord=0.20)
    p.plotExtrapolatedPolar()
    
    plt.show()
    

