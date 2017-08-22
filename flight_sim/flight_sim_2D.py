#!/usr/bin/env python

# Polaire type 2 Re(Cl) = Re0 * sqrt(Cl)
# remplaces les angles par un repere de rotation
# v 63 : introduction girouette 
from __future__ import print_function, division

from pprint import pprint

if 0:
    import inspect
    for e in dir(inspect):
        if "is" in e : print(e)
    quit()
from inspect import * # getmembers, ismemberdescriptor, isfunction, ismethod, isclass, isbuiltin
IS_LIST = (isabstract,
isawaitable,
isbuiltin,
isclass,
iscode,
iscoroutine,
iscoroutinefunction,
isdatadescriptor,
isframe,
isfunction,
isgenerator,
isgeneratorfunction,
isgetsetdescriptor,
ismemberdescriptor,
ismethod,
ismethoddescriptor,
ismodule,
isroutine,
istraceback,)

from collections import OrderedDict

import sys

#sys.path.insert(0, "/mnt/LinuxData/users/baptige/PYTHON/FlightSimulators")
sys.path.insert(0, "..")
import Profils.polaires as polaires
print(dir(polaires.Profil))
    
POLAR_DIR = '../Profils/Polaires' # relatif a l'appelant !
    
import math

import numpy as np

#print(dir(np))

#POLAR_DIR = '../Profils/Polaires'

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
    
#-------------------------------------------------------------------------------

def strCenteredOver(s, n):
    oversize = n - len(s)
    nSpaceBefore = nSpaceAfter = 0
    if oversize <= 0 : return s[:n]
    
    nSpaceAfter  = oversize//2
    nSpaceBefore = oversize - nSpaceAfter
    return "."*nSpaceBefore + s + "."*nSpaceAfter

def getSubAttr(obj, attrName):
    if '.' not in attrName :
        return getattr(obj, attrName)
        
    names = attrName.split('.')
    subObj = obj
    for name in names:
        subObj = getattr(subObj, name)
        
    return subObj

def str2FctParams(s):
    s = s.strip()
    fn, s = s.split('(')
    fctName = fn.strip()
    ps = s.split(')')
    prms = ps[0].split(',')
    params = [p.strip() for p in prms]
    return fctName, params

def appendToDict(dic, key, datas):
    if key in dic :
        lst = dic[key]
    else :
        lst = list()
        dic[key] = lst
    lst.append(datas)
        
        
RAD2DEG = 180.0/np.pi
DEG2RAD = np.pi/180.0

def rotate(v, a):
    x  , y   = v
    cos, sin = math.cos(a), math.sin(a)
    return x*cos - y*sin, y*cos + x*sin

def pv(v1, v2): return v1[0]*v2[1] - v1[1]*v2[0]
    
if 0:
    b0  = (-2.0, 0.0); f0  = (0.0, -1.0)
    b90 = ( 0.0,-2.0); f90 = (1.0,  0.0)
    print("b0 , f0  :", pv(b0 , f0 ))
    print("b90, f90 :", pv(b90, f90))
    quit()
    
vecDim, rotDim = 2, 1 # 2, 1 = simul2d; 3, 3 = simul3d 
vecShape = (vecDim,) # numpy
rotShape = (rotDim,)
mRotShap = (vecDim, vecDim)

def mRot(a): # rad
    cosa, sina = math.cos(a), math.sin(a)
    return np.array(((cosa, sina),(-sina, cosa)))

def microRot(da): # rad
    #cosa, sina = math.cos(a), math.sin(a)
    return np.array(((1.0, da),(-da, 1.0)))
    
def dmRot_da(a): # rad
    dcos_da, dsin_da = -math.sin(a), math.cos(a)
    return np.array(((dcos_da, dsin_da),(-dsin_da, dcos_da)))
    
def mRotSpeed(da): # rad
    #dcos_da, dsin_da = -math.sin(a), math.cos(a) ???????????????
    return np.array(((0.0, da),(-da, 0.0)))
    
def norm(v):
    vl = v.reshape((1,-1))
    #nv2 = v[0]**2 + v[1]**2
    nv2 = vl @ vl.T
    if nv2.shape != (1,1) : print("v.shape :", v.shape, nv2.shape)
    return np.sqrt(nv2)
        
if 0:
    a = np.pi/6
    v = np.array((np.cos(a), np.sin(a)))
    print("norm(v) :", norm(v))

    mra = mRot(0.0)
    pmra = mra
    mra[:,:] = mRot(a)
    print("a :%6.3f :\n"%a, mra)
    
    b = a
    mrb = np.zeros(mRotShap)
    mrb[:,:] = mRot(b)
    print("b :%6.3f :\n"%a, mrb)
    
    mra[:,:] = mra @ mrb # mra @= mrb not yet suported
    print("pmra :\n", pmra)
    
    mra[:,:] = mra @ mrb # mra @= mrb not yet suported
    print("pmra :\n", pmra)
    
    quit()

gravity = np.zeros(vecShape)
shoulderHeight = 1.6 # m 

def infoArray(objName, arrayName, arrayNames2Info):     
    if '[' in arrayName:
        name, idxs = arrayName.split('[')
        idxs = idxs[:-1].strip()
    else:
        name = arrayName
        idxs = ':' # all
    varName = objName + '_' + name 
    #print("arrayName : %s[%s]"%(name, idxs))
    try:
        firstIdx, nbFloats = arrayNames2Info[name] # find in second dict
    except KeyError :
        print("'%s' not found in arrayNames2Info: '%s' !"%(name, arrayNames2Info))
        return None
    if idxs != ':' :
        # only [:] or [n] understood
        #print("idxs : '%s'"%idxs)
        idx = int(idxs)
        assert idx < nbFloats
        firstIdx += idx
        nbFloats = 1
        varName += '_%d'%idx 
    #print("\t%-16s => "%arrayName, varName , firstIdx, nbFloats)
    return varName, firstIdx, nbFloats

#----------------------------------------------------            
class Recorder(object):

    def __init__(self, simulation):
        self.objDict     = simulation.objDict
        self.t           = simulation.t
        self.floatNb     = self.t.size
        self.arraylst    = list()
        self.objsInfo    = OrderedDict()
        self.tPltLst     = OrderedDict()
        self.varLstXY    = list()
        self.fctLstXY    = list()
        
    def addToLst(self, objName, attribMethodNames):
        obj = self.objDict[objName]
        arrayNames2Info = OrderedDict()
        for memberName in attribMethodNames:
            subObj = getSubAttr(obj, memberName) # , default)
            #if memberName == "wing.lift" : print("*** addToLst :", memberName, subObj)
            typMember = type(subObj)
            #print("type(member) :", typMember)
            assert typMember is np.ndarray
            self.arraylst.append(subObj)
            nbFloats = subObj.size
            arrayNames2Info[memberName] = (self.floatNb, nbFloats)
            self.floatNb += nbFloats
        self.objsInfo[objName] = arrayNames2Info

    def start(self, nbRecord=10):
        print("Recorder.start")
        self.genHeaders()
        self.nbRecord = nbRecord
        m = 1
        for e in self.arraylst:
            m += e.size
        #print("m : %2d"%m)
        self.records = np.zeros((self.nbRecord, m))
        self.recordIdx = 0
        self.record()
        
    def record(self):
        if self.recordIdx >= self.nbRecord :
            print("Recorder full !") 
            return # !!!!!!
        #print("recording %3d"%self.recordIdx) 
        record = self.records[self.recordIdx, :]
        record[0] = self.t
        m = 1
        #print("record lst :", self.lst)
        for e in self.arraylst:
            nbFloat = e.size
            #print("record e :\n", e)
            record[m:m+nbFloat] = e
            m += nbFloat
        self.recordIdx += 1

    def genHeaders(self):
        l1 = " i  : t (s) " ; l2 = " i  : t (s) "
        for objName, arrayNames2Info in self.objsInfo.items():
            s = "" ; ns = 0
            for arrayName, (offset, size) in arrayNames2Info.items():
                n = size*8
                #if size >= 2 : n += 2 # for parenthesis
                s  += ";" + strCenteredOver(arrayName, n)
                ns += n + 1
            l2 += s
            l1 += ";" + strCenteredOver(objName, ns-1)
        self.headers = [l1, l2] 
        
    def rec2Str(self, i):
        record = self.records[i]
        ss = ["%6.3f"%record[0]]
        m = 1
        for e in self.arraylst:
            nbFloat = e.size
            if nbFloat == 1 :
                f = record[m]
                m += 1
                ss.append("%7.3f"%f)
            else:
                e = record[m:m+nbFloat]
                m += nbFloat
                fs = ["%7.3f"%f for f in e]
                ss.append(",".join(fs))
        line = "; ".join(ss)
        return "%4d: %s"%(i, line)
    
    def replay(self, n=4):
        print("replay :")
        
        nRec = self.recordIdx # len(self.records)
         
        for line in self.headers:
            print(line)
        if nRec <= 2*n : 
            range0 = (range(nRec))
            for i in range0 : print(self.rec2Str(i))
        else:
            ranges = (range(0,n), range(nRec-n, nRec))
        
            for i in ranges[0]: print(self.rec2Str(i))
            print("....")
            for i in ranges[1]: print(self.rec2Str(i))      
    
    # ex: addPlotVsTime(('posX', 'm', (0, 4.0)), ('moment', 'Nm', None), ('energy', 'J', None))
    def addPlotVsTime(self, plotList):
        #print("addPlotVsTime :", plotList)
        for name, unitY, rangeY in plotList:
            #print("name, unitY, rangeY :", name, unitY, rangeY) 
            varsLst = list()
            fctsLst = list()
            self.tPltLst[name] = (unitY, rangeY, varsLst, fctsLst)
         
    # ex: addToPlot('alphas', ('Plane1', ('stabAlpha' , 'wingAlpha' )))
    def addToPlot(self, pltName, objNameArraysNames):
        #print("Recorder.addToPlot '%s' :"%pltName)
        try:
            unitY, rangeY, varsLst, fctsLst = self.tPltLst[pltName]
        except KeyError :
            print("'%s' not found in PlotVsTime: '%s' !"%(pltName, self.tPltLst))
            return
            
        #print("objNameArraysNames :", objNameArraysNames)    
        objName, arraysNames = objNameArraysNames
        obj = self.objDict[objName]
        #print("    from '%s'"%objName)
        try:
            arrayNames2Info = self.objsInfo[objName] # find in main dict
        except KeyError :
            print("'%s' not found in objsInfo: '%s' !"%(objName, self.objsInfo))
            return
            
        #print("arraysNames :", arraysNames)
        for arrayName in arraysNames:
            if '(' in arrayName:
                fctName, params = str2FctParams(arrayName)
                #print("Function : %s(%s)"%(fctName, ",".join(params)))
                fct = getSubAttr(obj, fctName)
                paramsNameFirstIdxNbFloats = list()
                for param in params:
                    #print("param :", param)
                    nameFirstIdxNbFloats = infoArray(objName, param, arrayNames2Info)
                    if nameFirstIdxNbFloats is None: break
                    paramsNameFirstIdxNbFloats.append(nameFirstIdxNbFloats)
                if len(paramsNameFirstIdxNbFloats) >= 1 :
                    fctsLst.append((fct, paramsNameFirstIdxNbFloats))
            else:   
                nameFirstIdxNbFloats = infoArray(objName, arrayName, arrayNames2Info)
                if nameFirstIdxNbFloats is not None:
                    varsLst.append(nameFirstIdxNbFloats)
            
    def addToPlotXY(self, objNameArraysNames):
        print("addToPlotXY : objNameArraysNames :", objNameArraysNames)    
        objName, arraysNames = objNameArraysNames
        obj = self.objDict[objName]
        #print("    from '%s'"%objName)
        try:
            arrayNames2Info = self.objsInfo[objName] # find in main dict
        except KeyError :
            print("'%s' not found in objsInfo: '%s' !"%(objName, self.objsInfo))
            return
            
        #print("arraysNames :", arraysNames)
        for arrayName in arraysNames:
            if '(' in arrayName:
                fctName, params = str2FctParams(arrayName)
                #print("Function : %s(%s)"%(fctName, ",".join(params)))
                fct = getSubAttr(obj, fctName)
                paramsNameFirstIdxNbFloats = list()
                for param in params:
                    #print("param :", param)
                    nameFirstIdxNbFloats = infoArray(objName, param, arrayNames2Info)
                    if nameFirstIdxNbFloats is None: break
                    paramsNameFirstIdxNbFloats.append(nameFirstIdxNbFloats)
                if len(paramsNameFirstIdxNbFloats) >= 1 :
                    self.fctLstXY.append((fct, paramsNameFirstIdxNbFloats))
            else:   
                nameFirstIdxNbFloats = infoArray(objName, arrayName, arrayNames2Info)
                if nameFirstIdxNbFloats is not None:
                    self.varLstXY.append(nameFirstIdxNbFloats)
            
    def plot(self):
        print("Recorder.plot()")
        import matplotlib.pyplot as plt

        nbSubplot = len(self.tPltLst)
        
        rcds = self.records[:self.recordIdx, :]
        ts = rcds[:, 0]
        
        minTime   = 0.0
        maxTime   = np.amax(ts)
        timeRange = (minTime, maxTime)
        timeTicks = np.arange(minTime, maxTime+1e-9, maxTime/20)

        fig = plt.figure(figsize=(20, 11))
        
        # default left=0.125,right=0.9,bottom=0.1,top=0.9,wspace=0.2,hspace=0.2
        fig.subplots_adjust(top=0.96,bottom=0.03, left=0.04,right=0.98, hspace=0.4)
        
        #tight_layout(pad=1.08, h_pad=None, w_pad=None, rect=None)       
        i = 0
        for name, (unitY, rangeY, varsLst, fctsLst) in self.tPltLst.items():
            i += 1
            ax = fig.add_subplot(nbSubplot, 1, i)
            ax.set_xlim(timeRange)
            ax.set_xticks(timeTicks)                                                       
            if rangeY is not None:
                ax.set_ylim(rangeY)
                minY, maxY = rangeY
                deltaY = maxY - minY
                ax.set_yticks(np.arange(minY, maxY+deltaY*1e-6, deltaY/8))
            ax.grid()
                
            for varName, firstIdx, nbFloats in varsLst:
                ax.plot(ts, rcds[:, firstIdx:firstIdx+nbFloats], '-')
                
            for fct, paramsNameFirstIdxNbFloats in fctsLst:
                params = list()
                for paramNameFirstIdxNbFloats in paramsNameFirstIdxNbFloats: 
                    varName, firstIdx, nbFloats = paramNameFirstIdxNbFloats
                    param = rcds[:, firstIdx:firstIdx+nbFloats]
                    params.append(param)
                    print("fct :%s ; param shape"%fct.__name__, param.shape)
                calculated = fct(*params)
                print("calculated.shape", calculated.shape)
                ax.plot(ts, calculated, '-')
            plt.title(name)
        
        if 1: # XY graph
            fig = plt.figure(figsize=(15, 10))
            
            varName, firstIdx, nbFloats = self.varLstXY[0] # pos
            xys  = rcds[:, firstIdx:firstIdx+nbFloats]
            
            varName, firstIdx, nbFloats = self.varLstXY[1] # wingForce
            wfxys = rcds[:, firstIdx:firstIdx+nbFloats] # rapport de bras de levier
            
            varName, firstIdx, nbFloats = self.varLstXY[2] # stabForce
            sfxys = rcds[:, firstIdx:firstIdx+nbFloats]
            
            if 0:
                varName, firstIdx, nbFloats = self.varLstXY[3]
                dxys = rcds[:, firstIdx:firstIdx+nbFloats]
            
                varName, firstIdx, nbFloats = self.varLstXY[4] # accel
                axys = rcds[:, firstIdx:firstIdx+nbFloats]
            
            plt.plot(xys[:,0], xys[:, 1], '.-b')
            plt.quiver(xys[:,0], xys[:, 1], wfxys[:,0], wfxys[:, 1], color='r', scale=1000, width=0.002) #
            plt.quiver(xys[:,0], xys[:, 1], sfxys[:,0], sfxys[:, 1], color='m', scale= 200, width=0.002) #
            #plt.quiver(xys[:,0], xys[:, 1],  dxys[:,0],  dxys[:, 1], color='r', scale=10, width=0.002) #
            #plt.quiver(xys[:,0], xys[:, 1],  axys[:,0],  axys[:, 1], color='c', scale=100, width=0.002) #
            plt.axis('equal')
            plt.grid()
        plt.show()

class Simulation(object):
    def __init__(self, wind=None):
        self.t        = np.zeros((1,)) # time in sec
        self.objDict  = dict()
        self.wind     = wind
        self.recorder = Recorder(self)
        
    def add(self, obj, name):
        obj.simul = self    
        self.objDict[name] = obj
        
    #def addToRecorderList(self, objName, floatArraysNames): self.recorder.addToLst(objName, floatArraysNames)
    def addToRecorderList(self, *args) : self.recorder.addToLst     (*args)    
    def addPlotVsTime    (self, *args) : self.recorder.addPlotVsTime(*args)     
    def addToPlot        (self, *args) : self.recorder.addToPlot    (*args)
    def addToPlotXY      (self, *args) : self.recorder.addToPlotXY  (*args)
    def plot             (self, *args) : self.recorder.plot         (*args)
        
    def info(self):
        #print(self.atm.info())
        print(self.plane.info())
    
    def step(self, t, dt):
        for name, obj in self.objDict.items():
            obj.step(t, dt)

    def start(self, end_t=1.0, dt=10e-3, underSampling=1):
        self.dt = dt
        self.underSampling = underSampling
        self.end_t = end_t - 1e-3*self.dt
        for name, obj in self.objDict.items():
            obj.launch()
        self.end_t = end_t - 1e-3*self.dt
        nbRecord = int(self.end_t/self.dt)//self.underSampling + 2
        self.recorder.start(nbRecord=nbRecord)
        self.samplingCnt = 0
        self.run()
        
    def run(self):
        self.end = 0
        while self.end == 0 :
            if self.t >= self.end_t :
                print("normal stop at %6.3f s"%self.t)
                self.end = 1 # normal stop by time
            self.step(self.t, self.dt)
            self.t += self.dt
            
            self.samplingCnt += 1
            if self.samplingCnt >= self.underSampling : 
                self.recorder.record()
                self.samplingCnt = 0

        print("End simulation")
        
    def info(self):
        print("Simulation info")

class MechanicPoint(object):

    def __init__(self, mass=1.0, initialPos=(0.0,10.0), initialSpeed=(0.0, 0.0)):
        self.simul  = None 
        self.mass         = mass         # kg
        self.initialPos   = initialPos   # m       
        self.initialSpeed = initialSpeed # m/s
        self.accel  = np.zeros(vecShape) # m/s2
        self.speed  = np.zeros(vecShape) # m/s
        self.pos    = np.zeros(vecShape) # m
        self.force  = np.zeros(vecShape) # Newton 
        
    def info(self):
        return "mass:%6.3f kg"%(self.mass)
            
    def energyKin(self, speed): # speed could be np.array
        #print("MechanicPoint.energyKin : speed.shape :", speed.shape)
        print("speed.shape :", speed.shape)
        n, m = speed.shape
        spd = np.zeros((n,1))
        spd[:,0] = speed[:,0]**2 + speed[:,1]**2
        print("spd.shape :", spd.shape)
        return 0.5 * self.mass * spd
        
    def energyPot(self, altitude): # altitude could be np.array
        return self.mass * g * altitude
           
    def energyTot(self, speed, altitude):
        ep = self.energyPot(altitude)
        print("ep.shape :", ep.shape)
        
        ek = self.energyKin(speed) 
        print("ek.shape :", ek.shape)
        
        et = ek + ep
        print("et.shape :", et.shape)
        
        return et
    
    def launch(self):
        self.pos  [:] = self.initialPos
        self.speed[:] = self.initialSpeed
        self.force[:] = 0.0, self.mass * -g # suppose objet stand on ground
        self.accel[:] = gravity[:] + self.force/self.mass
        
    def step(self, t, dt=0.1):
        self.accel[:] = gravity[:] + self.force/self.mass
        self.speed += self.accel * dt
        self.pos   += self.speed * dt
        
def CalculateMassInertiaCenterOfGravity(massDist):
    # moment of inertia / Mass*dist**2
    mass, balance,  inertia = 0.0, 0.0, 0.0
    for m, d in massDist:
        mass    += m
        balance += m * d   
        inertia += m * d**2
    centerOfGravity = balance / mass    
    print("mass : %5.3f kg, GC: %5.3f m, inertia:%6.4f kgm2"%(mass, centerOfGravity, inertia))
    return mass, centerOfGravity, inertia
    
class MechanicObject(MechanicPoint):

    def __init__(self, massDist=((0.2,-0.6),(0.2,0.0),(0.6,0.2)), gCfixed=False,\
                 initialPos=(0.0,0.0), initialAngle=0.0, initialSpeed=(0.0, 0.0), initialRot=0.0):
                 
        self.massDist = massDist # (mass, dist) from Center
        mass, centerOfGravity, inertia = CalculateMassInertiaCenterOfGravity(self.massDist)
        
        self.gCfixed = gCfixed # MassCenter no speed : girouette
        if self.gCfixed : initialSpeed = 0.0
        MechanicPoint.__init__(self, mass=mass, initialPos=initialPos, initialSpeed=initialSpeed)
        
        self.Jz           = inertia      # kg.m2
        self.initialAngle = initialAngle # rd
        self.initialRot   = initialRot   # rd/s
        
        self.d2a_dt2 = np.zeros(rotShape) # rd/s2
        self.dij_dt  = np.zeros(mRotShap) # /s
        self.repere  = mRot(initialAngle)
        self.repereU = self.repere[0,:]
        self.da_dt   = np.zeros(rotShape) # rd/s
        #self.angle   = np.zeros(rotShape) # rd
        self.torque  = np.zeros(rotShape) # N*m
    
    def angle(self, vec):
        a = np.arctan2(vec[:,1], vec[:,0])
        print("angle() :", a[0:5])
        return a*RAD2DEG

    def info(self):
        return MechanicPoint.info(self) + "; RotInertia:%6.3f kg*m2"%(self.Jz)

    def absPosSpeedAtInternPos(self, posInObject):
        pos   = self.pos   + posInObject @ self.repere
        speed = self.speed + posInObject @ self.dij_dt
        return pos, speed
        
    def launch(self):
        MechanicPoint.launch(self)
        self.repere[:,:] = mRot(self.initialAngle)
        if 1:
            u = self.repere[0,:]
            rs = np.array((u,))
            a = self.angle(rs)
            print("initialAngle : %6.1f deg"%(a))
        self.da_dt[0]    = self.initialRot
        
    def step(self, t, dt=0.1):
        if not self.gCfixed : MechanicPoint.step(self, t, dt=dt)
        
        self.d2a_dt2[0]  = self.torque/self.Jz
        self.da_dt[0]   += self.d2a_dt2 * dt
        self.dij_dt[:,:] = mRotSpeed(self.da_dt) @ self.repere
        self.repere[:,:] = self.repere @ mRot(self.da_dt * dt)
        
class Wind(object):
    def __init__(self, windGradient=None, v0=-7.0, h0=2.5, z0=0.1): #25km at 2.5m
        self.windGradient = windGradient
        self.v0 = v0
        self.h0 = h0
        self.z0 = z0
        
    def step(self, t, dt=0.1):
        pass    
        
    def speed(self, pos):
        h = pos[1]
        vec = np.zeros(vecShape)
        if not self.windGradient is None: # linear
            windSpeed = self.windGradient * h
        else: # exponential
            if h <= self.z0 : 
                windSpeed = 0.0
            else:
                windSpeed = self.v0 *(np.log(h/self.z0)/np.log(self.h0/self.z0))
        vec[0] = windSpeed
        #print("Wind.speed :", vec)
        return vec

if 0:
    w0 = Wind() #v0, h0, z0)
    w1 = Wind(windGradient=0.5)
    print("Wind profile:")
    for i in range(2, 100):
        h = 0.5 * i
        speed0 = w0.speed((0.0, h))
        speed1 = w1.speed((0.0, h))
        print("%5.1f m :%5.1f m/s :%5.1f m/s"%(h, speed0[0], speed1[0]))
    quit()

class Pilote(object):
    def __init__(self, durationCmds=((0.1,(0.0,)),(0.7,(-0.3,)),(0.45,(-0.7,))), gain=1.0):
        self.ramp = 2.0/0.5 # -100% to +100% in 0.5 sec
        self.gain = gain
        self.elevCmd = np.zeros((1,)) # rd  elevator neg to go up
        self.elevCmd[:] = durationCmds[0][1]
        self.durationCmds = durationCmds # [(duration, values), (duration, values) ...] 
        self.timeOfnextChange = 0.0
        
    def readCmds(self):
        if self.idx >= len(self.durationCmds):
            duration, cmds = 1e12, (0.0,) # zero forever
        else:
            duration, cmds = self.durationCmds[self.idx]
            self.idx += 1
        self.timeOfnextChange += duration
        self.cmds = cmds
        print("cmds :", self.cmds)
        
    def launch(self):
        self.idx = 0
        self.readCmds()
    
    def step(self, t, dt, alt, repere, speed):
        dmax = self.ramp * dt
        if t >= self.timeOfnextChange: self.readCmds()
        delta = self.cmds[0] - self.elevCmd[0]
        if   delta < -dmax : delta = -dmax
        elif delta >  dmax : delta =  dmax
        self.elevCmd[0] += delta
        
class AeroSurf(object):
    def __init__(self, name, span, chord, pos, angle, profil):
        self.name   = name
        self.span   = span  # m
        self.chord  = chord # m
        self.area   = self.span * self.chord # m2
        self.pos    = pos   # relative to plane GC
        self.angle  = angle # relative to plane main axe (calage aile)
        self.profil = profil
        self.kWing  = 0.5 * ro * self.area
        
        self.alpha  = np.zeros(1,)
        self.lift   = np.zeros(1,)
        self.drag   = np.zeros(1,)
        self.moment = np.zeros(1,)

    def aeroForce(self, airSpeed, alphaPlane, incidenceCmd): # m/s, rd, rd
        alphaDeg = (alphaPlane + self.angle + incidenceCmd) * RAD2DEG
        '''
        if   alphaDeg >  180.0 : alphaDeg -= 360.0
        elif alphaDeg < -180.0 : alphaDeg += 360.0
        '''
        self.alpha[:] = alphaDeg
        
        Re = speed2Reynolds(airSpeed, self.chord)
        Cl, Cd, Cm = self.profil.ClCdCm(self.alpha, Re)
        # Cd *= 0.8 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        kWingSquareSpeed = self.kWing * airSpeed**2
        
        self.lift[:]   = Cl * kWingSquareSpeed
        self.drag[:]   = Cd * kWingSquareSpeed
        self.moment[:] = Cm * kWingSquareSpeed * self.chord    
        
class Plane(MechanicObject):

    def __init__(self, massDist=((0.2,-0.6),(0.2,0.0),(0.6,0.2)), girouette=False, \
                initialPos=(0.0,shoulderHeight), initialAngle=0.0, initialSpeed=(10.0, 0.0), initialRot=0.0,\
                wingProfilName="xf-rg14-il", wingSpan=1.8, wingChord=0.19, wingPos=(-0.010, 0.0),\
                stabSpan=0.45, stabChord=0.10, stabPos=(-0.650, 0.0), stabDeg=-2.0, wingDeg=5.0,\
                GCajust=0.0, wind=None, pilote=None):
                
        MechanicObject.__init__(self, massDist=massDist, gCfixed=girouette,\
                initialPos=initialPos, initialAngle=initialAngle,\
                initialSpeed=initialSpeed, initialRot=initialRot)             
                
        self.profils = polaires.Profils(profilsDir=POLAR_DIR)
        self.pilote  = pilote
        self.wind    = wind
        self.GCajust = np.array((GCajust, 0.0)) # x positif = centrage avant
        
        wingProfil   = self.profils.add(wingProfilName)
        wingAngle    = wingDeg * DEG2RAD # rd angle with frame
        wingPos      = np.array(wingPos) - self.GCajust
        self.wing    = AeroSurf("wing", wingSpan, wingChord, wingPos, wingAngle, wingProfil)
        
        tailProfil   = self.profils.add("xf-n0009sm-il")
        stabPos      = np.array(stabPos) - self.GCajust # -0.650, 0.20 = rear stab T shape
        stabAngle    = stabDeg * DEG2RAD # rd callage  
        self.stab    = AeroSurf("stab", stabSpan, stabChord, stabPos, stabAngle, tailProfil)
        self.elevCmd2StabIncidence = -0.1/1.0 # rd/100% up: cmd pos and incidence neg
        
        stabVolume = (self.stab.area * -self.stab.pos[0]) / (self.wing.area * self.wing.chord)
        print("stabVolume : 0.300 < %5.3f < 0.900 "%stabVolume)
        
        self.wingForce  = np.zeros(vecShape)
        self.stabForce  = np.zeros(vecShape)
        
        self.wingTorque = np.zeros(rotShape)
        self.stabTorque = np.zeros(rotShape)
                
        self.wingAirSpeed = np.zeros((1,))
        self.stabAirSpeed = np.zeros((1,))
        
        self.alphaPlane   = np.zeros(rotShape)
        self.alphaPlAtTail= np.zeros(rotShape)
                
    def info(self):
        fmt = "mass:%6.3fkg; wingArea:%6.1fdm2; stabArea:%6.2fdm2; FAI_area:%6.2fdm2"
        return fmt%(self.mass, self.wingArea*1e2, self.stabArea*1e2, (self.wingArea+self.stabArea)*1e2)
        
    # q = 0.5* ro * V**2 ; Fz = q * S * Cz ; Fx = q * S * Cx

    def aeroForceOf(self, aeroSurf, incidenceCmd=0.0):
        absolutePos, absSpeed = self.absPosSpeedAtInternPos(aeroSurf.pos)
        if self.wind :
            windSpeed = self.wind.speed(absolutePos)
            if 3.0 < self.simul.t < 3.1 : print("%6.3f s : absolutePos :%6.3f m; windSpeed :%6.3f m/s"%(self.simul.t, absolutePos[1], windSpeed[0]))
            airSpeed = windSpeed - absSpeed
        else:
            airSpeed = - absSpeed
            
        airSpeedNorm = np.sqrt(airSpeed.dot(airSpeed)) # norm() not exist
        airSpeed_u   = airSpeed/airSpeedNorm
        airSpeed_v   = np.array((airSpeed_u[1], -airSpeed_u[0])) # perpendicular
        
        #airRepere[0,:] = airSpeed_u
        #airRepere[1,:] = airSpeed_v
        
        alphaPlane = np.arcsin(pv(-airSpeed_u, self.repere[0,:]))
        
        aeroSurf.aeroForce(airSpeedNorm, alphaPlane, incidenceCmd)
        
        force = airSpeed_u * aeroSurf.drag + airSpeed_v * aeroSurf.lift 
        brasDeLevier = aeroSurf.pos @ self.repere
            
        torque = pv(brasDeLevier, force)
        return force, torque, airSpeedNorm, alphaPlane*RAD2DEG
            
    def aeroForces(self):
        self.wingForce[:], self.wingTorque[:], self.wingAirSpeed[:], self.alphaPlane[:] = self.aeroForceOf(self.wing)
        
        incidenceCmd = 0.0 if self.pilote is None else self.pilote.elevCmd[0]*self.elevCmd2StabIncidence 
        self.stabForce[:], self.stabTorque[:], self.stabAirSpeed[:], self.alphaPlAtTail[:] = self.aeroForceOf(self.stab, incidenceCmd)
        
        self.force[:]  = self.wingForce + self.stabForce
        self.torque[:] = self.wing.moment + self.stabTorque + self.wingTorque  # + stab.moment # neglectable
        
    def g_number(self, lift): # lift could be np.array ot time
        # projection de la force sur l'axe vertical local
        #print("MechanicPoint.g_number")
        gravityReactionForce = g *self.mass
        return lift / gravityReactionForce
        
    def launch(self):
        if self.pilote is not None : self.pilote.launch()
        MechanicObject.launch(self)
        self.aeroForces()
    
    def step(self, t, dt=0.1):
        if self.pilote is not None : 
            self.pilote.step(t, dt=dt, alt=self.pos[1], repere=self.repere, speed=self.speed[1])
        MechanicObject.step(self, t, dt=dt)
        self.aeroForces()
                
#==============================================================================

gravity[1] = -g 
windGradient =  0.00  # m/s / m
GCajust      = 0.010
stabDeg      = 0.0
wingDeg      = 2.0
girouette    = False
wingSpan     = 1.8
wingChord    = 0.19
initialPos   = (0.0,shoulderHeight)
if 1:
    sim  = Simulation()
    sim.addPlotVsTime( (('pos'      , 'm'   , (-20.0, 60.0)), # ( 0.0, 4.0)),
                        ('angle'    , 'deg' , ( -8.0,  8.0)), #(-np.pi, np.pi)), # ,(-3.14, 3.14)
                        ('airSpeed' , 'm/s' , (8.0,40.0)), # ( 6.0,11.0)), 
                        #('rotSpeed' , 'rd/s', (-1.5, 2.5)),
                        ('cmds'     , ''    , (-1.0, 1.0)), 
                        ('alphas'   , 'deg' , (-8.0, 8.0)), # (-6.0, 6.0)), 
                        ('lifts'    , 'N '  , (-40.0,120.0)), # (-1.0, 1.0)), 
                        ('moment'   , 'Nm'  , (-4.0, 8.0)), # (-1.0, 1.0)), 
                        ('g_number' , ''    , (-4.0,20.0)), # (-1.0, 1.0)), 
                        ('energy'   , 'J'   , (-200.0,1000.0)),
                        )) # ( 0.0,50.0))) )
    if 1:
        end_t = 6.0
        wingProfilName = "xf-rg14-il"
        if 0: # giroutte
            girouette    = True
            wingSpan     = 0.100
            wingChord    = 0.020
            windGradient = -1.0  # 10m/s at 10m  m/s / m
            initialPos   = (0.0, 10.0)
            initialSpeed = (0.0, 0.0)
            initialAngle = 8.0 * DEG2RAD
            durationCmds = ((4.0,(0.0,)), (4.0,(-0.8,)), (4.0,(0.0,)), (4.0,(0.8,)), (4.0,(0.0,)))
            end_t = 20.0
        elif 0: # giroutte avec ailes
            GCajust      = -0.020
            wingProfilName="ApproximNACA0009"
            girouette    = True
            windGradient = -1.0  # 10m/s at 10m  m/s / m
            initialPos   = (0.0, 10.0)
            initialSpeed = (0.0, 0.0)
            initialAngle = 8.0 * DEG2RAD
            durationCmds = ((4.0,(0.0,)), (4.0,(-0.5,)), (4.0,(0.0,)), (4.0,(0.5,)), (4.0,(0.0,)))
            end_t = 20.0
        elif 0: # little slope
            initialSpeed = (10.0, 0.5)
            initialAngle = -1.5 * DEG2RAD
            durationCmds = ((8.0,(0.76,)), (0.49,(-0.4,)), (9.0,(0.3,)))
        elif 0: # little slope centrage arriere
            GCajust      = -0.020
            initialSpeed = (10.0, 0.0)
            initialAngle = 0.0 * DEG2RAD
            durationCmds = ((8.0,(0.59,)), (0.49,(-0.4,)), (9.0,(0.3,)))
        elif 1: # big slope centrage arriere
            wingProfilName="ApproximNACA0009"
            GCajust      = -0.050
            stabDeg      = 0.0
            wingDeg      = 1.0
            initialSpeed = (10.0, -4.0)
            initialAngle = -0.40 
            durationCmds = ((5.0,(0.0,)), (0.5,(-0.20,)),(50.0,(0.0,)))
            end_t = 12.0
        elif 0: # looping in gradian
            windGradient = -0.50  # m/s / m
            GCajust      = 0.0
            initialSpeed = (30.0, 0.0)
            initialAngle = -1.2 * DEG2RAD # * np.pi
            #durationCmds= ((0.5,(-0.1,)), (4.29,(0.75,)), (1.5,(-0.05,)), (5.65,(0.75,)), (9.0,(0.0,)))
            durationCmds = ((0.5,(-0.1,)), (4.49,(0.61,)), (1.4,(-0.064,)), (4.80,(0.72,)), (9.0,(-0.01,)))
            end_t = 12.0
        elif 0: # oval looping in gradian 18G 30 m/s after first loop 38 m/s
            windGradient = -0.50  # m/s / m
            initialSpeed = (30.0, 0.0)
            initialAngle = -1.2 * DEG2RAD # * np.pi
            durationCmds = ((0.5,(-0.1,)), (0.4,(0.33,)), (2.0,(-0.1,)), (1.6,(1.1,)), (1.1,(-0.30,)), (0.77,(1.2,)), (9.0,(-0.1,)))
            end_t = 7.0
        elif 0: # cool oval looping in gradian 11G 25 m/s after first loop 29.7 m/s
            windGradient = -0.50  # m/s / m
            initialSpeed = (25.0, 0.0)
            initialAngle = -1.2 * DEG2RAD # * np.pi
            durationCmds = ((0.5,(-0.1,)), (0.5,(0.4,)), (2.0,(-0.0,)), (2.7,(1.0,)), (1.3,(-0.30,)), (1.07,(0.99,)), (9.0,(-0.11,)))
            end_t = 10.0
        elif 0: # cool oval looping in gradian centrage moins avant 11G 25 m/s after first loop 30.2 m/s maxHeight=50m => wind=25m/s ie 90km/h 
            GCajust      = 0.000
            windGradient = -0.50  # m/s / m
            initialSpeed = (25.0, 0.0)
            initialAngle = -1.2 * DEG2RAD # * np.pi
            durationCmds = ((0.5,(-0.1,)), (0.5,(0.4,)), (2.0,(-0.0,)), (2.7,(0.96,)), (1.3,(-0.30,)), (1.07,(0.91,)), (9.0,(-0.11,)))
            end_t = 10.0
        elif 0: #  cool oval looping in log gradient centrage moins avant 10G 25 m/s after first loop 21.5 m/s maxHeight=36 m =>  
            GCajust      = 0.000
            windGradient = None # gradiant logarithmique 25 km/h a 2.5 m
            initialSpeed = (25.0, 0.0)
            initialAngle = -1.2 * DEG2RAD # * np.pi
            durationCmds = ((0.5,(-0.1,)), (0.4,(0.3,)), (0.9,(-0.1,)), (1.6,(0.99,)), (1.1,(-0.30,)), (1.10,(1.074,)), (9.0,(-0.072,)))
            end_t = 7.0
        elif 0: # gradiant 0.25 cool oval looping in gradient centrage moins avant 11G 25 m/s after first loop 22.7 m/s maxHeight=30.5m => wind=7.6m/s ie 27 km/h 
            GCajust      = 0.000
            windGradient = -0.25  # m/s / m
            initialSpeed = (25.0, 0.0)
            initialAngle = -1.2 * DEG2RAD # * np.pi
            durationCmds = ((0.5,(-0.1,)), (0.5,(0.3,)), (0.5,(-0.0,)), (2.6,(0.99,)), (1.0,(-0.30,)), (1.25,(1.2,)), (9.0,(-0.04,)))
            end_t = 8.0
        elif 0: # cool looping
            initialSpeed = (30.0, 0.0)
            initialAngle = -1.2 * DEG2RAD # * np.pi
            durationCmds = ((0.5,(-0.1,)), (4.6,(0.75,)), (9.0,(0.0,)))
        elif 1: # looping
            initialSpeed = (30.0, 0.0)
            initialAngle = -1.2 * DEG2RAD # * np.pi
            durationCmds = ((0.5,(-0.1,)), (1.6,(0.80,)), (1.6,(0.40,)), (1.6,(0.84,)), (9.0,(0.0,)))
        elif 0: # square looping
            initialSpeed = (40.0, 0.0)
            initialAngle = -1.2 * DEG2RAD # * np.pi
            durationCmds = ((1.0,(-0.13,)), (0.55,(0.95,)), (0.8,(-0.10,)), (0.8,(0.94,)), (9.0,(-0.33,)))
        elif 0: # passage dos
            initialSpeed = (30.0, 0.0)
            initialAngle = -1.2 * DEG2RAD # * np.pi
            durationCmds = ((0.5,(-0.1,)), (1.6,(0.80,)), (0.5,(0.40,)), (0.6,(-0.0,)), (9.0,(-0.71,)))
        elif 0: # looping no gravity  New ALGO !!!!!!
            gravity[1] = 0.0 
            GCajust    = -0.020
            wingProfilName="ApproximNACA0009"
            initialSpeed, initialAngle = (10.0, 0.0), np.pi*0 # (0.0, 10.0), np.pi/2 #(7.0, 7.0), np.pi/4
            durationCmds = ((0.25,(0.0,)),(13.7,(1.0,)),(9.0,(0.0,)))
            end_t = 10.0
        else: # looping no gravity
            gravity[1] = 0.0 
            wingProfilName="ApproximNACA0009"
            initialSpeed, initialAngle = (10.0, 0.0), np.pi*0 # (0.0, 10.0), np.pi/2 #(7.0, 7.0), np.pi/4
            durationCmds = ((0.1,(0.0,)),(13.7,(1.0,)),(9.0,(0.0,)))
            end_t = 15.0
            
        wind  = Wind(windGradient=windGradient)
        pi2 = Pilote(durationCmds=durationCmds)
        #pl2 = Plane(wind=wind, initialSpeed=(10.0, 0.5), initialAngle=0.10, pilote=pi2)
        pl2 = Plane(wind=wind, girouette=girouette, wingProfilName=wingProfilName, GCajust=GCajust,\
                    wingSpan=wingSpan, wingChord=wingChord,\
                    stabDeg=stabDeg, wingDeg=wingDeg, initialPos=initialPos,initialSpeed=initialSpeed, initialAngle=initialAngle, pilote=pi2)
        sim.add(pl2, 'Plane2') 
        sim.addToRecorderList('Plane2',('pos','repereU','speed','da_dt','stabTorque','wingTorque',\
                              'stab.alpha','stab.lift','wing.alpha','wing.lift',\
                              'pilote.elevCmd','wingAirSpeed','stabForce','wingForce',\
                              'wing.moment','force','alphaPlane','alphaPlAtTail'))
        
        sim.addToPlot('pos'     , ('Plane2', ('pos',)))
        sim.addToPlot('angle'   , ('Plane2', ('angle(repereU)',)))
        sim.addToPlot('airSpeed', ('Plane2', ('wingAirSpeed',)))
        #sim.addToPlot('rotSpeed', ('Plane2', ('da_dt',)))
        sim.addToPlot('cmds'    , ('Plane2', ('pilote.elevCmd',)))
        sim.addToPlot('alphas'  , ('Plane2', ('stab.alpha', 'wing.alpha', 'alphaPlAtTail'))) # ,'alphaPlane'
        sim.addToPlot('moment'  , ('Plane2', ('stabTorque', 'wingTorque', 'wing.moment')))
        sim.addToPlot('energy'  , ('Plane2', ('energyKin(speed)','energyPot(pos[1])','energyTot(speed,pos[1])')))
        sim.addToPlot('lifts'   , ('Plane2', ('wing.lift', 'stab.lift')))
        sim.addToPlot('g_number', ('Plane2', ('g_number(wing.lift)',)))
        
        sim.addToPlotXY(('Plane2', ('pos','wingForce','stabForce')))
        
    if 0:       
        pl1 = Plane(initialSpeed=(8.15,-0.2), initialAngle=-0.024)
        sim.add(pl1, 'Plane1') 
        sim.addToRecorderList('Plane1',('pos','angle','speed','stabTorque','wingTorque','stabAlpha','wingAlpha','wingAirSpeed'))
        
        sim.addToPlot('pos'      , ('Plane1', ('pos[1]',)))
        sim.addToPlot('angle'    , ('Plane1', ('angle' ,)))
        sim.addToPlot('airSpeed' , ('Plane1', ('wingAirSpeed',)))
        sim.addToPlot('alphas'   , ('Plane1', ('wingAlpha', )))
        sim.addToPlot('moment'   , ('Plane1', ('stabTorque', 'wingTorque')))
        sim.addToPlot('energy'   , ('Plane1', ('energyKin(speed)','energyPot(pos[1])')))
        
    if 0:
        print("inspect :")
        if 0:
            for isSomething in IS_LIST:
                print("%s :"%(isSomething.__name__))
                pprint(getmembers(sim, isSomething), width=120)
            name = 'objDict' # 'wind'
            print("getattr_static(sim, '%s') :"%name, getattr_static(sim, name))
        else:
            print("isfunction :")      
            pprint(getmembers(sim, isfunction)) 
            print("ismethod :")      
            pprint(getmembers(sim, ismethod)) 
            print("ismemberdescriptor :")      
            pprint(getmembers(sim, ismemberdescriptor)) # 
        quit()
        
    sim.start(end_t=end_t, dt=0.010, underSampling=1)
    sim.recorder.replay(n=20)
    sim.plot()
else:
    sim = Simulation()
    mp1 = MechanicPoint()
    sim.add(mp1, 'MechPt1') 
    sim.addToRecorderList('MechPt1',('accel','speed','pos','eKinet','ePoten','eTotal'))
    
    mp2 = MechanicPoint(initialPos=(0.0, 0.0), initialSpeed=(0.0, 10.0))
    sim.add(mp2, 'MechPt2') 
    sim.addToRecorderList('MechPt2',('accel','speed','pos','eKinet','ePoten','eTotal'))
    
    sim.start(end_t=1.0, dt=0.1)
    
    sim.recorder.replay() 
    
    sim.plot()  
    
#--------------------------------------   

