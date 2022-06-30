# -*- coding: utf-8 -*-
import rlcompleter, readline
readline.parse_and_bind('tab:complete')
import sys, os, glob
import numpy as np

"diphasique_statistiques_0.967020.txt"
"diphasique_statistiques_1.030141.txt"

def buildStringOfColumns(fic):
   st = ""
   f = open(fic,"r")
   lines = f.readlines()
   for val,key in [(int(s.split()[2])-1,s.split()[4])  for s in lines if "# colonne" in s]:
      st += key+" "
      pass
   f.close()
   return st

def getTempsIntegrationFile(fic):
   f = open(fic,"r")
   lines = f.readlines()
   tintegration = lines[0].split()[2]
   f.close()
   return float(tintegration)

def getTempsIntegrationList(fics):
   lt = []
   for fic in fics:
      lt.append(getTempsIntegrationFile(fic))
      pass
   return np.array(lt)

def getFloatBtwIters(tab, ltintegration, it,it2):
   dt = ltintegration[it2]-ltintegration[it]
   res = (ltintegration[it2]*tab[it2,:] - ltintegration[it]*tab[it,:])/dt
   return res


class StatReader:
    def __init__(self, head, fold, finstatfile = None, inistatfile = None, diphasique=True, master_plot=False):
	self.fic_stats = glob.glob("*phasique_statistiques_*.txt")
        self.fic_stats.sort(key=lambda item: float(item.strip("monodiphasique_statistiques_.txt")))
	if (finstatfile == None):
	    finstatfile = self.fic_stats[-1]
	if (inistatfile != None):
            if (not os.path.isfile(inistatfile)):
                raise Exception("Missing input ini file : %s"%(inistatfile))
            pass
        #
        if (not os.path.isfile(finstatfile)):
            raise Exception("Missing input fin file : %s"%(finstatfile))
        #
        self.finstatfile = finstatfile
    	# self.dvar = dns.buildDicoColonnes(finstatfile)
        self.inistatfile = inistatfile
	self.loadStats()
        return 
    #
    #
    def loadStats(self):
	print "\tLoading stats... "
	fin = np.loadtxt(self.finstatfile)
        z = fin[:,0]
        if (self.inistatfile == None):
            print "\t\tReading file : ", self.finstatfile
            time = float(self.finstatfile.strip("monodiphasique_moyenne_spatiale_.txt"))
            tintegration =  getTempsIntegrationFile(self.finstatfile)
            print "\t\tTime list : ", time, "     (Integration over dt=", tintegration, ")"
	    mat = fin
	    dt = tintegration
	else:
	    fic_stats = [self.inistatfile, self.finstatfile]
            ltimes = np.array([float(f.strip("monodiphasique_moyenne_spatiale_.txt")) for f in fic_stats])
            tintegration = getTempsIntegrationList(fic_stats)
	    ini = np.loadtxt(self.inistatfile)
            # Compute mat from the difference btw the 2 times in resu :
	    mat = getFloatBtwIters(np.array([ini,fin]),tintegration,0,1)
	    mat[:,0] = z
            dt = tintegration[1]-tintegration[0]
            print "\t\tTime list : ", ltimes, "     (Integration over dt=", dt, ")"
	    print "\t\tMat is restricted to the average between those 2 times"
	    pass
        self.mat = mat
	self.z = z
        self.nz =len(self.z)
	#
	print "\t\tCell centers are btw ", z[0], " and ", z[-1], " (nz=%d)"%self.nz
        tab_dz = z[1:]-z[:-1]
        dz=tab_dz.mean()
        print "\t\tMesh step size along z (Delta z) : ", dz
        if (tab_dz.max()-tab_dz.min())/dz>1e-10:
            raise Exception("Maillage a pas variable!!!")
        #
	self.dz = dz
	self.tintegration = dt
	#
        return 
    #
    #
    def saveStats(self):
        st = buildStringOfColumns(self.finstatfile)
	np.savetxt("Stats.dt_ev", self.mat, header=st)
	return 
    #
diph = StatReader("DNS", "DIPH", inistatfile="diphasique_statistiques_14.919289.txt", finstatfile = "diphasique_statistiques_15.006300", master_plot=True)
diph.saveStats()
