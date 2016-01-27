from numpy import array, arange, outer, sin
from saxs_mod import compute_intensity
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
    'VAL':'V','TRP':'W','TYR':'Y','MSE':'M','GLD':'X'}

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
    "CH3":31.89, "NH": 7.64, "NH2": 12.79, "NH3":17.94, "OH":14.28,"SH":25.1
}
# Dictionary for number of hydrogen bound per atom type
with open(localpath+'pdbHbound_dict.json','rt') as f:
    Hdict = json.load(f)

class PDB:
    def __init__(self, filename=None):
        self.residues = None
        self.atoms    = None
        self.residues = dict()
        self.chains   = []
        self.r        = None
        self.Pr       = None
        self.Iq       = None
        self.beamprofile = None
        if filename is None:
            self.filename = 'File is not specified.'
        else:
            self.read(filename)

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

        # get sequence from PDB
        self.parse_sequence()
        

    def parse_sequence(self):
        if self.atoms is not None:
            seqpair = [(s.resnum, s.resname, s.chain) for s in self.atoms\
                 if s.resname in aa3]
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
        Iq,Nr = compute_intensity(self.q,x,y,z,fflist)
        r = arange(0,Nr*dr,dr)

        if self.q[0]==0.0:
            self.Iq = Iq/Iq[0]
        else:
            print "Not supported yet"

        # then convert to P(r)
        qr = outer(self.q,r)
        Iq = Iq - Iq.min()
        Iq = Iq/Iq.max()
        Pr = Iq.dot(qr * sin(qr) * dq)
        self.r = r
        self.Pr= Pr
        
    def write(self, filename, chain):
        sel = [s for s in self.atoms if s.chain==chain and s.het==False]
        with open(filename,'wt') as f:
            for s in sel:
                f.write(s.__str__() + "\n")
        print "Successfully written Chain {:s} into {:s}".format(chain, filename)


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
        # Carbon atom types
        if self.atomname not in ['C','N','O','P','S','H']:
            try:
                self.boundH = Hdict[self.resname][self.atomname]
                if self.elem=='C':
                    if self.boundH==1:
                        self.vexc = Vexc["CH"]
                    elif self.boundH==2:
                        self.vexc = Vexc["CH2"]
                    elif self.boundH==3:
                        self.vexc = Vexc["CH3"]
                    else:
                        self.vexc = Vexc[self.elem]
                if self.elem=='N':
                    if self.Hbound==1:
                        self.vexc = Vexc["NH"]
                    elif self.Hbound==2:
                        self.vexc = Vexc["NH2"]
                    elif self.Hbound==3:
                        self.vexc = Vexc["NH3"]
                    else:
                        self.vexc = Vexc[self.elem]
                if self.elem=='O':
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

            