from numpy import array, arange, outer, sin, loadtxt
from saxs_mod import compute_intensity
from scipy.optimize import minimize
import json
import os
localpath = os.path.dirname(os.path.realpath(__file__)) + '/'
# # (atomlabel,self.number,self.atomname,self.resname,\
#             self.chain,self.resnum,self.alternate,self.x,self.y,self.z,\
#             self.occ,self.Bfactor,self.elem,self.charge)
# PDB format, old python format style. Please revise if needed.
pdbfmt = "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"

# 3-Letter Amino Acid look-up table
# MSE is selenoMet
aa3 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',\
    'GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M',\
    'ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T',\
    'VAL':'V','TRP':'W','TYR':'Y','MSE':'M','GLD':'X','AuX':'X', 'UNK': '.'}

# zero form factor
FF0 = {
    "H":0.999953,"He":0.999872,"Li":2.99,"Be":3.99,"B":4.99,"C":5.9992,"N":6.9946,\
    "O":7.9994,"F":8.99,"Ne":9.999,"Na":10.9924,"Mg":11.9865,"Al":12.99,"Si":13.99,\
    "P":14.9993,"S":15.9998,"Cl":16.99,"Ar":17.99,"K":18.99,"Ca2+":18.0025,"Cr":23.99,\
    "Mn":24.99,"Fe2+":24.0006,"Co":26.99,"Ni":27.99,"Cu":28.99,"Zn2+":27.9996,"Se":33.99,\
    "Br":34.99,"I":52.99,"Ir":76.99,"Pt":77.99,"Au":78.9572,"Hg":79.99,"CH":6.99915,\
    "CH2":7.99911,"CH3":8.99906,"NH":7.99455,"NH2":8.99451,"NH3":9.99446,"OH":8.99935,\
    "OH2":9.9993,"SH":16.9998
}

FF0net = {
    "H":-0.720147,"He":-0.720228,"Li":1.591,"Be":2.591,"B":3.591,"C":0.50824,"N":6.16294,\
    "O":4.94998,"F":7.591,"Ne":6.993,"Na":7.9864,"Mg":8.9805,"Al":9.984,"Si":10.984,\
    "P":13.0855,"S":9.36656,"Cl":13.984,"Ar":16.591,"K":15.984,"Ca2+":14.9965,"Cr":20.984,\
    "Mn":21.984,"Fe2+":20.9946,"Co":23.984,"Ni":24.984,"Cu":25.984,"Zn2+":24.9936,"Se":30.9825,\
    "Br":31.984,"I":49.16,"Ir":70.35676,"Pt":71.35676,"Au":72.324,"Hg":73.35676,"CH":-0.211907,\
    "CH2":-0.932052,"CH3":-1.6522,"NH":5.44279,"NH2":4.72265,"NH3":4.0025,"OH":4.22983,\
    "OH2":3.50968,"SH":8.64641   
}

# Excluded volume
Vexc = {
    "H": 5.15, "C": 16.44, "N": 2.49, "O":9.13, "S": 19.86, "P": 5.73,\
    "Au":19.86, "Se": 9.0, "Na": 6.538, "K": 14.71, "CH": 21.59, "CH2": 26.74,\
    "CH3":31.89, "NH": 7.64, "NH2": 12.79, "NH3":17.94, "OH":14.28,"SH":25.1,\
    "Cl": 22.45, "Mg": 21.69, "Ca": 51.63
}

# VDW radius for atoms
# V = 4/3 * pi * r^3
# r^3 = V * 3/(4 * pi)  
rads = {
    "C":1.7,"N":1.55,"O":1.52,"F":1.47,"P":1.8,"S":1.8,"Cl":1.75,\
    "H":1.07,"Cu":1.4,"K":2.75,"Na":2.27, "Mg":1.73, "Co":2.0, "Au": 1.66,\
    "He": 1.4, "Hg": 1.55, "Ni":1.63, "Li":1.82, "Br":1.85, "I":1.98, "Pb":2.02
}
# Dictionary for number of hydrogen bound per atom type
with open(localpath+'pdbHbound_dict.json','rt') as f:
    Hdict = json.load(f)

class PDB:
    def __init__(self, filename=None):
        self.residues = None
        self.atoms    = None
        self.Natoms   = 0
        self.residues = dict()
        self.chains   = []
        self.r        = None
        self.Pr       = None
        self.Iq       = None
        self.beamprofile = None
        self.saxscalc = False # flag for SAXS calculation
        if filename is None:
            self.filename = 'File is not specified. Creating a new PDB object'
        else:
            self.read(filename)
        if self.atoms is not None:
            self.calcCofM() # compute the center of mass

    def read(self, filename):
        self.filename = filename
        fhd = open(filename, 'rt')
        dat = fhd.readlines()
        fhd.close()
        
        dat = [d for d in dat if d.startswith('ATOM') or d.startswith('HETATM')]
        # parse into Atom objects
        self.atoms = []
        for d in dat:
            if d.startswith('ATOM'):
                atm = Atom(hetatm=False, number=int(d[6:11]), name=d[12:16],\
                    alternate=d[16:17],resname=d[17:20],chain=d[21:22],\
                    resnum=int(d[22:26]),x=float(d[30:38]),y=float(d[38:46]),\
                    z=float(d[46:54]),occ=float(d[54:60]),Bfactor=float(d[60:66]),\
                    element=d[76:78],charge=d[78:80])        
            elif d.startswith('HETATM'):
                atm = Atom(hetatm=True, number=int(d[6:11]), name=d[12:16],\
                    alternate=d[16:17],resname=d[17:20],chain=d[21:22],\
                    resnum=int(d[22:26]),x=float(d[30:38]),y=float(d[38:46]),\
                    z=float(d[46:54]),occ=float(d[54:60]),Bfactor=float(d[60:66]),\
                    element=d[76:78],charge=d[78:80])   
            if atm.chain not in self.chains:
                self.chains.append(atm.chain)
            self.atoms.append(atm)

        self.Natoms = len(self.atoms)
        # get sequence from PDB
        self.parse_sequence()
        

    def parse_sequence(self):
        if self.atoms is not None:
            try:
                seqpair = [(s.resnum, s.resname, s.chain) for s in self.atoms\
                     if s.resname in aa3]
            except:
                print "Residue name is not recognized."
                
            for chain in self.chains: # do per chain
                # only get protein-only residues
                unique_res = []
                for res in seqpair: 
                    if res not in unique_res and res[2] is chain:
                        unique_res.append(res)
                # find min and max of sequence names
                resmin = 99999
                resmax = 1
                for res in unique_res:
                    if res[0]<=resmin:
                        resmin = res[0]
                    if res[0]>resmax:
                        resmax = res[0]
                res_range = range(1,resmax+1)
                residues = []
                counter = 0

                for i in res_range:
                    if unique_res[counter][0]!=i:
                        residues.append('-')
                    elif unique_res[counter][0]==i:
                        residues.append(aa3[unique_res[counter][1]])
                        counter += 1
                sequence = ''.join(residues)
                self.residues[chain] = sequence
    # temporary functions to return arrays for computation
    def calc_profiles(self,**kwargs):
        dr = 0.5
        if 'q' in kwargs:
            self.q = kwargs['q']
            dq = self.q[1]-self.q[0]
        elif 'qmax' in kwargs:
            if 'dq' in kwargs:
                dq = kwargs['dq']
                self.q = arange(0,kwargs['qmax'],dq)
            else:
                print "No dq was specified, using 0.001"
                dq = 0.001
                self.q = arange(0,kwargs['qmax'],dq)

        else:
            print "No q or qmax,dq was specified. Calculating to 0.5 by 0.001"
            dq = 0.001
            self.q = arange(0,0.5,dq)

        xyz = array([[e.x,e.y,e.z,e.ff0] for e in self.atoms])
        x   = xyz[:,0]
        y   = xyz[:,1]
        z   = xyz[:,2]
        fflist=xyz[:,3]

        if 'beamprofile' in kwargs:
            beamfn = kwargs['beamprofile']
            beamdat = readPDH(beamfn)
            beamdat = beamdat[:,:2] # get only first two columns
            # and normalize so that integral of beam is 0.5
            beamdat[:,0] *= 0.1
            dy = beamdat[1,0] - beamdat[0,0]
            beamdat[:,1] = 0.5*(beamdat[:,1]/(beamdat[:,1].sum() * dy))
            # convert unit to angstrom
            
            Iq,Nr = compute_intensity(self.q,x,y,z,fflist,beamdat)
        else:
            beamdat = array([[1,1],[1,1]]) # dummy input for fortran program
            Iq,Nr = compute_intensity(self.q,x,y,z,fflist,beamdat)

        r = arange(0,Nr*dr,dr)

        if self.q[0]==0.0:
            self.Iq = Iq/Iq[0]
        else:
            self.Iq = Iq/Iq[0]

        # then convert to P(r)
        qr = outer(self.q,r)
        Iq = Iq - Iq.min()
        Iq = Iq/Iq.max()

        Pr = Iq.dot(qr * sin(qr) * dq)

        # consider removing this if you have some smearing
        self.r = r
        self.Pr= Pr

        self.saxscalc = True

    def addAtom(self, **kwargs):        
        if self.atoms is None:
            self.atoms = []
        self.Natoms += 1
        newatom = Atom(number=self.Natoms,**kwargs)
        self.atoms.append(newatom)       
        
    def writeChain(self, filename, chain):
        sel = [s for s in self.atoms if s.chain==chain and s.het==False]
        with open(filename,'wt') as f:
            for s in sel:
                f.write(s.__str__() + "\n")
        print "Successfully written Chain {:s} into {:s}".format(chain, filename)

    def save(self, filename):
        with open(filename,'wt') as f:
            for s in self.atoms:
                f.write(s.__str__() + "\n")
        print "Successfully written PDB into {:s}".format(filename)

    def boundingBox(self, padding=10.0):
        if self.atoms is not None:
            coords = array([(a.x,a.y,a.z) for a in self.atoms])
            xmin,xmax = (coords[:,0].min(), coords[:,0].max())
            ymin,ymax = (coords[:,1].min(), coords[:,1].max())
            zmin,zmax = (coords[:,2].min(), coords[:,2].max())
            return ((xmin-padding,xmax+padding), (ymin-padding,ymax+padding),\
             (zmin-padding,zmax+padding))

    def fitSAXS(self, q_obs, Iq_obs, sd_obs, p0=[0.1,1e-4], beamprofile=None):
        # create a local function to return Chi^2
        def chisq(par,obs,calc,sd):
            s = par[0]
            b = par[1]
            Nobs = obs.size
            scalc = s*calc + b # computed profile
            chi2  = (obs - scalc)**2 / sd**2
            return chi2.sum()/Nobs
        
        if not self.saxscalc:
            # compute profiles if not yet
            if beamprofile is not None:
                self.calc_profiles(q=q_obs, beamprofile=beamprofile)
            else:
                self.calc_profiles(q=q_obs)

            calc = self.Iq
            opt = minimize(chisq, p0, args=(Iq_obs,calc,sd_obs,),method='L-BFGS-B')
        
        elif self.saxscalc:
            opt = minimize(chisq, p0, args=(Iq_obs,calc,sd_obs,),method='L-BFGS-B')

        print "Chi-sq value {:7.3f}".format(opt.fun)

        self.saxs_scale = opt.x[0]
        self.saxs_bg    = opt.x[1]

        return opt.fun

    def calcCofM(self):
        # computes the center of mass of atom
        xyz = [(a.x,a.y,a.z) for a in self.atoms]
        xyz = array(xyz)
        self.center = xyz.mean(axis=0)
        self.ch_centers = dict()
        for c in self.chains:
            wrkatoms = [(a.x,a.y,a.z) for a in self.atoms if a.chain==c]
            wrkatoms = array(wrkatoms)
            self.ch_centers[c] = wrkatoms.mean(axis=0)

    def translate(self,dx=0,dy=0,dz=0):
        # move entire PDB by dx,dy,dz
        for a in self.atoms:
            a.x += dx
            a.y += dy
            a.z += dz
        return True

    def nextChain(self):
        alphabets = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        takenchains = ''.join(self.chains)
        Nchains = len(takenchains)
        return alphabets[Nchains]

    def addPDB(self,PDBadd,x=0.0,y=0.0,z=0.0):
        # PDBadd is the PDB wished to be placed and appended to the current
        # PDB at position (x,y,z) by its center of mass

        # figure out how much to move
        dx = x-PDBadd.center[0]
        dy = y-PDBadd.center[1]
        dz = z-PDBadd.center[2]

        PDBadd.translate(dx,dy,dz)
        # figure out what chains exists in current PDB
        nc = self.nextChain()

        # and append a new chain to this structure
        i = 3
        for a in PDBadd.atoms:
            a.chain = nc
            a.number= self.Natoms + i
            i += 1

        for a in PDBadd.atoms:
            self.atoms.append(a)
        
        self.chains.append(nc)

        print "PDB successfully appended as chain %s" % (nc)

    def moveChain(self,chainID,xtarg,ytarg,ztarg):
        # assuming chain is on the origin !!!
        # get center of mass of chainID
        self.calcCofM() # compute center of mass
        cofm = self.ch_centers[chainID]
        dx = xtarg - cofm[0]
        dy = ytarg - cofm[1]
        dz = ztarg - cofm[2]

        for a in self.atoms:
            if (a.chain==chainID):
                a.x += dx
                a.y += dy
                a.z += dz


class Atom:
    def __init__(self,element='N',name='N',number=1,charge='',x=0,y=0,z=0,\
        occ=1.0, Bfactor=20.0, chain='A', alternate='',resname='ALA',\
        resnum=1,hetatm=False):
        self.het   = hetatm
        self.elem  = element.strip()
        self.atomname = name.strip()
        self.charge= charge.strip()
        self.chain = chain.strip()
        self.x     = x
        self.y     = y
        self.z     = z
        self.occ   = occ
        self.Bfactor = Bfactor
        self.number= number
        self.alternate = alternate
        self.resname=resname.strip()
        self.resnum =resnum
        self.boundH =0
        self.ff0  = None
        self.vexc = None
        fixname = self.elem
        # fix upper case
        if len(self.elem)==2:
            fixname = fixname[0] + fixname[1].lower()
        # look up radius   
        try:
            self.radius=rads[fixname]
        except KeyError:
            print "Element {:s} doesn't have a radius assigned".format(self.elem)
            self.radius = 1.6
        self.accountFF0() # compute net form factor

    def __str__(self):
        if self.het:
            atomlabel = 'HETATM'    
        else:
            atomlabel = 'ATOM'
        out = pdbfmt % (atomlabel,self.number,self.atomname,self.alternate,\
            self.resname,self.chain,self.resnum,'',self.x,self.y,self.z,\
            self.occ,self.Bfactor,self.elem,self.charge)
        return out

    def accountH(self):
        # right now only supports Protein
        # implement Nucleic Acid later
        if self.atomname not in ['C','N','O','P','S','H']:
            try:
                self.boundH = Hdict[self.resname][self.atomname]
                if self.elem=='C': # Carbon atom types
                    if self.boundH==1:
                        self.vexc = Vexc["CH"]
                    elif self.boundH==2:
                        self.vexc = Vexc["CH2"]
                    elif self.boundH==3:
                        self.vexc = Vexc["CH3"]
                    else:
                        self.vexc = Vexc[self.elem]
                if self.elem=='N': # Nitrogen atom types
                    if self.Hbound==1:
                        self.vexc = Vexc["NH"]
                    elif self.Hbound==2:
                        self.vexc = Vexc["NH2"]
                    elif self.Hbound==3:
                        self.vexc = Vexc["NH3"]
                    else:
                        self.vexc = Vexc[self.elem]
                if self.elem=='O': # Oxygen atom types
                    if self.Hbound==1:
                        self.vexc = Vexc["OH"]
                    else:
                        self.vexc = Vexc["O"]
                if self.elem=='S':
                    if self.Hbound==1:
                        self.vexc = Vexc["SH"]
                    else:
                        self.vexc = Vexc[self.elem]
            except:
                self.boundH=0
        elif self.atomname in Vexc:
            self.vexc  = Vexc[self.elem]
            self.boundH= 0
        else:
            self.boundH= 0
            self.vexc  = 0

    def accountFF0(self):
        self.accountH()
        if self.atomname not in ['C','N','O','P','S','H']:
            try:
                if self.atomname in FF0net:
                    # print "Atom found {:s}".format(self.atomname)
                    self.ff0  = FF0net[self.atomname]
                    self.boundH = 0
                else:
                    self.boundH = Hdict[self.resname][self.atomname]
                    if self.elem=='C':
                        if self.boundH==1:
                            self.ff0 = FF0net["CH1"]
                        elif self.boundH==2:
                            self.ff0 = FF0net["CH2"]
                        elif self.boundH==3:
                            self.ff0 = FF0net["CH3"]
                    if self.elem=='N':
                        if self.Hbound==1:
                            self.ff0 = FF0net["NH"]
                        elif self.Hbound==2:
                            self.ff0 = FF0net["NH2"]
                        elif self.Hbound==3:
                            self.ff0 = FF0net["NH3"]
                    if self.elem=='O':
                        if self.Hbound==2:
                            self.ff0 = FF0net["OH2"]
                        if self.Hbound==1:
                            self.ff0 = FF0net["OH"]
                    if self.elem=='S':
                        if self.Hbound==1:
                            self.ff0 = FF0net["SH"]
            except:
                try:
                    self.ff0 = FF0net[self.elem]
                    self.boundH=0
                except:
                    fixname = self.elem
                    fixelem = fixname[0] + fixname[1].lower()
                    self.ff0 = FF0net[fixelem]
                    self.boundH=0
        else:
            self.ff0 = FF0net[self.elem]

def readPDH(fn):
    try:
        with open(fn,'rt') as f:
            dat = f.readlines()
        beamdat = []
        if dat[1].startswith('SAXS'): #PDH file
            Npts = int(dat[2].split()[0])
            for i in range(5,Npts+1):
                tmp = [float(d) for d in dat[i].split()]
                beamdat.append(tmp)
            beamdat = array(beamdat)
        else:
            beamdat = loadtxt(fn)
    except:
        print "File %s not found."%(fn)
        beamdat = False

    return beamdat
