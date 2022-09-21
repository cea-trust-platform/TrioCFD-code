# -*- coding: utf-8 -*-
import rlcompleter, readline
readline.parse_and_bind('tab:complete')
import sys, os, glob
import numpy as np
from optparse import OptionParser

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
   print fic, tintegration
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

def getDiffBtwIters(tab, dt, it,it2):
   #Un-weighted difference
   res = (tab[it2,:] - tab[it,:])/dt 
   return res
   
def process_cmd_line(argv):
    """Processes the passed command line arguments."""
    parser = OptionParser(usage="usage: %prog [options]")
    
    parser.add_option("-e", "--end_stat_file", dest="endfile", type="string", default=None,
                      help="File *statistiques*txt to read the final cumulated stats.")
    
    parser.add_option("-i", "--ini_stat_file", dest="inifile", type="string", default=None,
                      help="File *statistiques*txt to read the initial cumulated stats.")

    parser.add_option( "--hide", action="store_true", dest="hide", 
                      help="Do not show plots.")

    #parser.set_defaults(folds="")
    
    (options, args) = parser.parse_args(argv)
    return options


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
	#
	# Trying to read instantaneous fields if present : 
	self.unsteady = True
	self.fininstfile = (self.finstatfile).replace("statistiques", "moyenne_spatiale")
	if (self.inistatfile):
	    self.iniinstfile = (self.inistatfile).replace("statistiques", "moyenne_spatiale")
	else:
	    self.iniinstfile = None
	    pass
	
	if not (os.path.isfile(self.fininstfile) and os.path.isfile(self.iniinstfile)):
	    print "One (or 2) instantaneous fields are missing. We cannot create unsteady file"
	    self.unsteady = False 
	    pass
	#
	self.loadStats()
        return 
    #
    #
    def loadStats(self):
	print "\tLoading stats... "
	fin = np.loadtxt(self.finstatfile)
        z = fin[:,0]
	#
	# Trying to read instantaneous fields if present : 
	if (os.path.isfile(self.fininstfile)):
 	    self.finu = np.loadtxt(self.fininstfile)
        else:
	    print("Last instantanous field is missing, we cannot load and create Final_snapshot.dt_ev file")
            self.finu = None
	    pass
	#
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
        if self.unsteady:
	    fic_insts = [self.iniinstfile, self.fininstfile]
	    print "\t\tUnsteady. Time difference btw: ", fic_insts
	    #ltimes = np.array([float(f.strip("monodiphasique_statistiques_moyenne_spatiale_.txt")) for f in fic_insts])
	    #print ltimes
	    #return
            #tintegration = getTempsIntegrationList(fic_insts)
	    iniu = np.loadtxt(self.iniinstfile)
	    #finu = np.loadtxt(self.fininstfile)
            # Compute mat from the difference btw the 2 times in resu :
	    self.matu = getDiffBtwIters(np.array([iniu,self.finu]),dt,0,1)
            print "\t\tTime list : ", ltimes, "     (Time gap dt=", dt, ")"
	    pass
	#
        li = glob.glob("*acceleration.out")
        if (len(li)==1):
           matacc=np.loadtxt(li[0])
           tdeb=ltimes[0]
           tfin=ltimes[1]
           itdeb=np.argmin(abs(matacc[:,1]-tdeb))
           itmax=np.argmin(abs(matacc[:,1]-tfin))
           t=matacc[itdeb:itmax,1]
           acc=matacc[itdeb:itmax,7]
           mean_acc=acc.mean()
           print "The mean acceleration source term is : ", mean_acc
           f = open("mean_acc.txt", 'w')
           f.write("%g %g %g\n"%(tdeb,tfin,mean_acc))
           f.close()
           pass
        #
        return 
    #
    #
    def saveStats(self):
        st = buildStringOfColumns(self.finstatfile)
        if self.unsteady:
           stu = buildStringOfColumns(self.fininstfile)
           np.savetxt("Unsteady.dt_ev", self.matu, header=stu)
           pass
	#
	if (self.finu != None):
           stu = buildStringOfColumns(self.fininstfile)
	   np.savetxt("Final_snapshot.dt_ev", self.finu, header=stu)
	   pass
	np.savetxt("Stats.dt_ev", self.mat, header=st)
	return 
    #

options = process_cmd_line(sys.argv[1:])

diph = StatReader("DNS", "DIPH", finstatfile=options.endfile, inistatfile=options.inifile, master_plot=True)
diph.saveStats()

