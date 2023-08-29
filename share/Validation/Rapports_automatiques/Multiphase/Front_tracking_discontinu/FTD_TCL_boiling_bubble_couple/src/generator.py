import lib2to3.fixes.fix_metaclass
import math
import sys, os, shutil, re, copy, subprocess, time
import numpy as np
from pip._internal.utils.misc import captured_output

if len(sys.argv) != 1:
    raise Exception(f"usage: python {sys.argv[0]} []")

coarsening = 0

cwd = os.path.join(os.getcwd())
srcdir = cwd
root = os.path.join(os.getcwd(), "..")
import socket
host=socket.gethostname()[0:5]
runs = os.path.join(root, f"RUNS-{host}")
if not os.path.isdir(root): os.mkdir(root)

file_list = ["post_run", "clean", "visit_Re_Nu2.py"]


def meshing(Nx_regu, Ny_regu, per_regu=0.333):
    """
   Nx_regu: Number of grid = Num of node - 1
    per_regu : percentage of regular mesh
   """
    fac = 1.05  # factor progressive > 1

    Nx_prog = math.log(1. - (1. - per_regu) / per_regu * (1. - fac) * Nx_regu, fac)
    Ny_prog = math.log(1. - (1. - per_regu) / per_regu * (1. - fac) * Ny_regu, fac)
    # print("[DIR X] number of NODES in progeresive zone :", int(Nx_prog) + 1)
    # print("[DIR Y] number of NODES in progeresive zone :", int(Ny_prog) + 1)
    return int(Nx_prog), int(Ny_prog)


class Case:
    v = True  # Verbose mode
    running_jobs = []  # The list of currently running jobs
    # Filling the list from some template files:
    templates = {}
    for k in ["sub_file_template", "template_steady.data", "template_source.data"]:
        with open(os.path.join(srcdir, k), "r") as template_file:
            templates[k] = template_file.readlines()

    # Are we on a cluster with a msub launcher?
    launcher = shutil.which("ccc_msub")
    bash = False
    if launcher is None:
        launcher = "./"
        bash = True
    else:
        launcher += " -A gen7712 "
        pass

    def __init__(self, cases, repo, resolution):
        self.cases = cases
        self.directory = repo
        self.resolution = resolution  # in µm
        self.dx = float(resolution) * 1e-6
        self.dict = copy.deepcopy(cases.dict)
        # 1/3 is regular mesh:
        pre_regu = 1. / 3.
        Nx1 = int(self.cases.Lx * pre_regu / self.dx)
        Ny1 = int(self.cases.Ly * pre_regu / self.dx)
        Nx2, Ny2 = meshing(Nx1, Ny1)
        Ny3 = int(self.cases.ZZ / self.dx)
        # print("Mesh info:", Nx1, Ny1, Nx2, Ny2)
        # Here we store NODES numbers :
        self.dict["Nx1"] = Nx1 + 1
        self.dict["Ny1"] = Ny1 + 1
        self.dict["Nx2"] = Nx2 + 1
        self.dict["Ny2"] = Ny2 + 1
        self.dict["Ny3"] = Ny3 + 1
        self.dict["nprocx"] = 3
        if self.bash:
            self.dict["nprocy"] = 4
        else:
            self.dict["nprocy"] = math.floor((Ny1+Ny2)/60)
        self.dict["nprocs"] = self.dict["nprocx"]*self.dict["nprocy"]
        self.options = []
        return

    def createCase(self):
        ok = True
        os.makedirs(self.directory, exist_ok=True)
        # The pre_run in steady creates and fills MESH directory
        self.createSteady()
        self.createNOTCL("NOTCL")
        for nmeso in self.cases.nmeso.split("-"):
            self.createTCLcase(f"TCL{nmeso}", nmeso=nmeso)
        return ok

    def createNOTCL(self, c):
        calc_dir = self.createTCLcase(c)
        # Specificities of calcs without TCL:
        fic = os.path.join(calc_dir, "source.data")
        with open(f"{fic}.tmp", "w") as out:
            with open(fic, "r") as fjdd:
                for line in fjdd.readlines():
                    # line = re.sub(f'^.*{"TCL"}.*', f'{}', line)
                    if "TCL" not in line:
                        out.write(line)
        os.rename(f"{fic}.tmp", fic)
        return

    def createTCLcase(self, c, nmeso=None):
        self.options.append(c)
        calc_dir = self.createCalculation(os.path.join(c, "R0"))
        # Specificities of Std calcs:
        with open(os.path.join(calc_dir, "injection.txt"), "w") as inj:
            phase = 0
            bubble = f" {phase} x^2+(y-{self.cases.r0}*cos({self.cases.theta}*pi/180.0))^2-{self.cases.r0}^2"
            inj.write(f"@time@ {bubble}\n")
            inj.write(f"100.0000 {bubble}\n")  # injection should have a line at a very large time

        if nmeso is not None:
            fic = (os.path.join(calc_dir, "source.data"))
            cmd = f"sed -i 's/n_extend_meso 4/n_extend_meso {nmeso}/' {fic}"
            cmd = ["sed", "-i", "-e", f"\"s/n_extend_meso 4/n_extend_meso {nmeso}/\"", f"{fic}"]
            st = " ".join(cmd)
            try:
                os.system(st)
            except Exception as e:
                try:
                    # Je ne comprends pas que :
                    # 1. p n'ai pas une portee au dela du try
                    # 2. cette commande ne fonctionne pas. Elle dit que le sed sort en erreur
                    p= subprocess.run(cmd, check=True, text=True)
                except Exception as e2:
                    p = subprocess.run(cmd)
                    print(p.stdout)
                    print(p.stderr)
                    p.wait()
                    raise e2
                raise e
        return calc_dir

    def createSteady(self):
        calc_dir = self.createCalculation("STEADY")
        # Specificities of Steady:
        with open(os.path.join(calc_dir, ".not_reprise"), "w") as fic:
            fic.write('')

    def createCalculation(self, fold="STEADY"):
        calc_dir = os.path.join(self.directory, fold)
        os.makedirs(calc_dir, exist_ok=True)
        shutil.copy(os.path.join(srcdir, "pre_run"), calc_dir)

        with open(os.path.join(calc_dir, "sub_file"), "w") as sf:
            nproc = self.dict["nprocs"]
            for line in self.templates["sub_file_template"]:
                if fold == "STEADY":
                    line = re.sub(f'^jdd=.*$', f'jdd=PAR_steady', line)
                line = re.sub(f'@case@', f'{self.cases.name}/M{self.resolution}/{fold}', line)
                line = re.sub(f'@nproc@', f'{nproc}', line)
                sf.write(line)
        os.chmod(os.path.join(calc_dir, "sub_file"), 0o755)

        jdd = "source"
        if fold == "STEADY":
            jdd = "steady"
        with open(os.path.join(calc_dir, f"{jdd}.data"), "w") as fjdd:
            for line in self.templates[f"template_{jdd}.data"]:
                for key in self.dict.keys():
                    val = self.dict[key]
                    line = re.sub(f'@{key}@', f'{val}', line)
                fjdd.write(line)

        for f in file_list:
            shutil.copy(os.path.join(srcdir, f), calc_dir)

        return calc_dir

    def getNumberOfOtherRunning(self):
        """

        :return: The number of jobs that are running but NOT from subprocess
                     so that they are not in
        """
        n_by_subprocess = len(self.running_jobs)
        nall = 0
        for root, dirs, files in os.walk(self.directory):
            for file in files:
                if file == ".running":
                    nall +=1
        others = nall - n_by_subprocess
        return others

    def waitAnyJob(self, nmax_jobs=2):
        """
        Blocks until some job are finished if nmax is reached
        :param self:
        :param nmax_jobs: Maximum number of simultaneous jobs authorised
        :return:
        """
        if (len(self.running_jobs) >= nmax_jobs):
            # On attend la fin d'un des jobs:
            while True:
                time.sleep(6)
                for job in self.running_jobs:
                    if job.poll() is not None:
                        # This job is finished. Remove from list and return
                        ok = (job.returncode == 0)
                        self.running_jobs.remove(job)
                        return ok
        return True

    def waitThisJob(self, run, subdir):
        output, error = run.communicate()
        ret_code = run.wait()
        ok = (ret_code == 0)  # Attend la fin du run steady
        if not ok:
            if run.stdout is not None: print(run.stdout)
            if run.stderr is not None: print(run.stderr)
            print(f"output = {output}")
            raise Exception(f"run {run.args} failed in {self.directory}/{subdir}")
        self.running_jobs.remove(run)
        return ok

    def run(self, cmd="sub_file", suff="", launcher=launcher, subdir="", waitme=False):
        if self.v: print("launching: ", f"{launcher}{suff}{cmd}".split(" "))
        run = subprocess.Popen(f"{launcher}{suff}{cmd}".split(" "),
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True)
        if self.bash or launcher == "./":
            nprocs = self.dict["nprocs"]
            if self.v: print(f"Job {cmd} in {subdir} is running under #{run.pid} on {nprocs} procs.")
            self.running_jobs.append(run)
            jobsteady_ID = run.pid
            if waitme:
                ok = self.waitThisJob(run, subdir)
            else:
                others = self.getNumberOfOtherRunning()# Other running jobs
                ok = self.waitAnyJob(nmax_jobs=2-others)  # To check if too many jobs are running already.
        else:
            stdout_data, stderr_data = run.communicate()
            print(stdout_data)
            jobsteady_ID = stdout_data.split(" ")[3].strip() # Supprime le saut de ligne
            ok = True
            try:
                int(jobsteady_ID)
            except Exception as e:
                if (type(e) == ValueError):
                    return print(e), False
                else:
                    raise e
            if self.v: print(f"Job {cmd} in {subdir} is running under JobID # {jobsteady_ID}")
        return jobsteady_ID, ok

    def runCase(self):
        if os.getenv("project_directory") is None:
            raise Exception("TrioCFD environment is required but not found")
        cwd = os.getcwd()
        os.chdir(self.directory)
        # Running the steady:
        os.chdir("STEADY")
        _, pre_ok = self.run(cmd="pre_run", launcher="./", subdir="STEADY", waitme=True)
        jobsteady_ID, run_ok = self.run(subdir="STEADY", waitme=True)  # Running sub_file on steady.
        ok = pre_ok and run_ok
        if self.bash:
            suff = ""
        else:
            suff = f"-E \"--dependency=afterok:{jobsteady_ID}\" "

        os.chdir("..")
        # Running different options:
        for opt in self.options:
            dest = os.path.join(opt, "R0")
            os.chdir(dest)
            _, pre2_ok = self.run(cmd="pre_run", launcher="./", subdir=dest, waitme=True)
            _, cas_ok = self.run(subdir=dest, suff=suff)  # Running sub_file on steady.
            ok = ok and pre2_ok and cas_ok
            os.chdir("../..")
        print(os.getcwd())
        os.chdir(cwd)
        return ok


class Cases:
    def __init__(self, case, resolution, nmeso, Lx, Ly, ZZ, dT, theta, r0=0.0002):
        self.name = case
        self.M = resolution  # in µm
        self.nmeso = nmeso
        self.Lx = Lx
        self.Ly = Ly
        self.ZZ = ZZ
        self.dT = dT
        self.theta = theta
        self.r0 = r0  # Initial bubble radius
        self.dict = {}
        self.dict["rmax"] = Lx
        self.dict["zmax"] = Ly
        self.dict["zsol"] = ZZ
        self.dict["dT"] = dT
        self.dict["theta"] = theta
        self.dict["delta_th"] = 0.0005  # m
        self.dcases = {}
        return

    def createCases(self, root='.'):
        ok = True
        os.makedirs(os.path.join(root, self.name), exist_ok=True)
        for resol in self.M.split("-"):
            repo = os.path.join(root, self.name, f"M{resol}")
            cas = Case(self, repo, resol)
            ok = cas.createCase() and ok
            if not ok: return ok
            self.dcases[(self.name, f"M{resol}")] = cas
        return ok

    def runCases(self):
        ok = True
        for cas in self.dcases.values():
            ok = cas.runCase() and ok
        return ok

    def __str__(self):
        st = f"Member of Cases: {self.name}"
        st += "\nSubdirs:"
        for t in self.dcases.keys():
            dest = t[0]
            for i in range(1, len(t)):
                dest = os.path.join(dest, t[i])
            st += f" {dest}"
        return st


configurations = ["T2"]
unit = np.ones(len(configurations))
rs = [0.009]
zs = [0.012]
zssol = [0.001]
dTs = [8.5]
thetas = [50.0]
sMs = [4.9]  # in µm
hMs = [3.8]  # in µm
Qmicros = [30.5]
Ms = ["40"]  # Mesh sizes in µm
nmesos = ["1"]
dts = 1e-8 * unit
tmaxs = 50e-3 * unit
r0 = [0.0002]  # initial bubble radius
nb_pas_dt_max = 1000000000 * unit

d_all = {}
for idx, config in enumerate(configurations):
    r = rs[idx]
    z = zs[idx]
    zz = zssol[idx]
    theta = thetas[idx]
    dT = dTs[idx]
    M = Ms[idx]
    nmeso = nmesos[idx]
    print(f"CONFIGURATION: {config} {r} {z} {M} $dT $theta $sM $hM $Qmicro $eth $r0 $Ntot")
    cas = Cases(config, resolution=M, nmeso=nmeso, Lx=r, Ly=z, ZZ=zz, dT=dT, theta=theta)
    ok = cas.createCases(root=runs)
    d_all[config] = cas
    print("CASES CREATED:")
    print(f"{cas.dcases.keys()}")

print("RUNNING CASES")
for idx, config in enumerate(configurations):
    cas = d_all[config]
    print(cas)
    cas.runCases()
