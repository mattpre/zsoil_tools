##    Reading ZSoil binary results
##    Copyright (C) 2015-2018  Matthias Preisig
##    Last modification: 2018/10/08
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU Affero General Public License as published
##    by the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU Affero General Public License for more details.
##
##    You should have received a copy of the GNU Affero General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import struct,re,numpy,math
from numpy import linalg as la
from struct import unpack
from struct import unpack_from
import sys

class ele_group:
    def __init__(self):
        self.type = ''
        self.nres_types = 0
        self.res_labels = []
        self.ncomp = []
        self.comp_labels = []

class ele_group_nl:
    def __init__(self):
        self.type = ''
        self.nres_types = 0
        self.res_labels = []
        self.ncomp = []
        self.comp_labels = []

class nodal_res:
    def __init__(self):
        self.type = ''
        self.nres_types = 0
        self.res_labels = []
        self.ncomp = []
        self.comp_labels = []

class mesh:
    def __init__(self):
        self.nElements = 0
        self.nNodes = 0

class ele_info:
    def __init__(self):
        self.inel = []
        self.mat = []
        self.EF = []
        self.LF = []
        self.ps = []    # property set after .eda file
        self.loc_syst = [] # rotation matrix for projecting results into ref_vect3D
        self.base = []
        self.area = [] # tributary areas for gauss points abs(dXdxi X dXdeta)
        self.type = []  # for contact: 0 for 2D, 1 for 3D, 2 for pile interface, 3 for tip interf.
        self.parent = []    # element number of beam or shell etc for contacts
        self.dir_nodes = [] # for beams, the numbers of nodes indicating local y-direction (2 per beam)
        self.nint = []
        self.nlayer = []

class nodal:
    def __init__(self):
        self.disp = []
        self.rot = []
        self.ppres = []
        self.pres_head = []
        self.r_disp = []
        self.r_rot = []
        self.r_ppres = []
        self.r_pres_head = []
        self.a_disp = []
        self.a_rot = []
        self.v_disp = []
        self.v_rot = []
        self.temp = []
        self.i_disp = []
        self.i_rot = []
        self.i_ppres = []
        self.i_pres_head = []
        
class vol:
    def __init__(self):
        self.stress = []
        self.strain = []
        self.pla_code = []
        self.str_level = []
        self.satur = []
        self.flu_veloc = []
        self.Et = []
##        self.DAMAGE = []
        # principal stresses are not automatically computed:
        self.princ = []
        self.invar = [] #I1, J2, p, q
        
class shell:
    def __init__(self):
        self.thick = []
        self.smforce = []
        self.smoment = []
        self.sqforce = []

class nlshell:
    def __init__(self):
        self.position = []
        self.stress = []
        self.strain = []
        self.pla_code = []
        self.str_level = []
        self.damage = []
        
class beam:
    def __init__(self):
        self.force = []
        self.moment = []
        self.dm = []
        self.edisp_tra = []
        self.edisp_rot = []

class nlbeam:
    def __init__(self):
        self.position = []
        self.stress = []
        self.strain = []
        self.pla_code = []
        self.str_level = []
        
class truss:
    def __init__(self):
        self.stress = []
        self.strain = []
        self.pla_code = []
        self.str_level = []
        self.force = []
        
class cnt:
    def __init__(self):
        self.stress = []
        self.tstress = []
        self.strain = []
        self.pla_code = []
        self.str_level = []

class pile:
    def __init__(self):
        self.beams = []
        self.cnt1D = []
        self.cnt0D = -1
        self.nBeams = 0
        self.volumics = []

class nail:
    def __init__(self):
        self.beams = []
        self.cnt1D = []
        self.cnt0D = -1
        self.nBeams = 0

class anchor:
    def __init__(self):
        self.trusses = []
        self.cnt1D = []
        self.cnt0D = -1
        self.nTrusses = 0
        
class mem:
    def __init__(self):
        self.smforce = []
        self.strain = []
        self.pla_code = []
        self.str_level = []

class property_set:
    def __init__(self):
        self.number = 0
        self.mat = 0
        self.matmodel = ''
        self.EF = 0
##        self.EF_active = []
        self.LF = 0
        self.group_labels = ['ELAS','GEOM','DENS','NONL','MAIN']
        self.data = dict()
        for gl in self.group_labels:
            self.data[gl] = []

class material:
    def __init__(self):
        self.label = ''
        self.number = 0
        self.ps = 0
        self.model = ''
        self.unit_weight = dict()
        self.elastic = dict()
        self.flow = dict()
        self.nonlinear = dict()
        self.initial = dict()        
        
class time_step:
    def __init__(self):
        #info from .his file:
        self.number = 0
        self.niter = 0
        self.conv_status = 1
        self.sf = 0.0
        self.plt_status = 0
        self.time = 0.0
        self.step_count = 0
        self.push_over = 0
        self.PSH = 0.0
        self.PSH_disp = 0.0
        self.solver = 0
        # from log file:
        self.Enorm = 0.0
        self.Fnorm = 0.0

        self.nModes = 0
        self.modes = []
        self.f0 = []
        self.effm = []
        self.pf = []

        # nodal results:
        self.nodal = nodal()

        # continuum results:
        self.vol = vol()

        # shell results:
        self.shell = shell()
        self.nlshell = nlshell()

        # beam results:
        self.beam = beam()
        self.nlbeam = nlbeam()

        # truss results:
        self.truss = truss()

        # contact results:
        self.cnt = cnt()

        # 1D contact results:
        self.cnt1D = cnt()

        # 0D contact results (e.g. pile tip interface):
        self.cnt0D = cnt()

        # membrane results:
        self.mem = mem()

class mode:
    def __init__(self,dim):
        self.nodal = nodal()
        self.nodal.disp = [[] for kd in range(dim)]

class zsoil_results:
    def __init__(self,pathname,problem_name,verbose_level=0):
        self.pathname = pathname
        self.problem_name = problem_name
        # verbose_level = 0: print every step read
        # verbose_level = 1: print every 100th step read
        # verbose_level = 2: suppress reading file messages
        # verbose_level = 3: suppress error messages
        self.verbose_level = verbose_level
        self.ele_group_labels = []
        self.ele_groups = []
        self.nEleGroups = 0
        self.ele_group_labels_nl = []
        self.ele_groups_nl = []
        self.nEleGroups_nl = 0
        self.nodal_res = []
        self.nNodalRes = 0
        self.nodalRead = False
        self.shellsRead = False
        self.shellsNLRead = False
        self.volumicsRead = False
        self.beamsRead = False
        self.beamsNLRead = False
        self.trussesRead = False
        self.contactsRead = False
        self.membranesRead = False

        self.property_sets = []

        self.jobtype = 0

        self.nSteps = 0
        self.steps = []
        self.out_steps = []
        self.out_diff_steps = []
        self.ref_vect2D = [0,1,0]
        self.nModes = 0

        self.LTF = []
        self.nLTF = 0
        self.EF = []
        self.nEF = 0
        self.materials = []

        self.accel_LTF = 0
        self.accel_mult = 0

        self.nElements = 0
        self.nNodes = 0
        self.nVolumics = 0
        self.nContacts = 0
        self.nBeams = 0
        self.nShells = 0
        self.nTrusses = 0
        self.nPiles = 0
        self.nNodalLinks = 0
        self.nMembranes = 0
        self.num_volumics = []
        self.num_shells = []
        self.num_beams = []
        self.num_trusses = []
        self.num_contacts = []
        self.num_membranes = []

        self.coords = []
        self.vol = ele_info()
        self.cnt = ele_info()
        self.beam = ele_info()
        self.shell = ele_info()
        self.truss = ele_info()
        self.mem = ele_info()

        # macro elements:
        self.piles = []
        self.nails = []
        self.anchors = []

        self.lists = []
        self.nodallinks = []

        self.logfile = 0
        
    def give_time_step(self,kt):
        if len(self.steps)<=kt:
            astep = time_step()
            astep.number = kt
            self.steps.append(astep)
        return self.steps[kt]

    def read_his(self):
        if self.verbose_level<2:
            print(self.pathname + '/' + self.problem_name + '.his')
        file = open(self.pathname + '/' + self.problem_name + '.his')

        lcount = 0
        while 1:
            line = file.readline()
            if not line: break
            v = line.split()
            if len(v)>6:
                s = self.give_time_step(lcount)
                s.type = int(v[0])
                s.niter = int(v[1])
                s.conv_status = int(v[2])
                s.sf = float(v[3])
                s.plt_status = int(v[4])
                s.time = float(v[5])
                s.step_count = int(v[6])
                s.push_over = int(v[7])
                s.PSH = float(v[9])
                s.PSH_disp = float(v[10])
                s.solver = int(v[14])
                lcount += 1
        self.nSteps = lcount

        file.close()

    def read_log(self):
        from zsoil_tools import logfile
        self.logfile = logfile.Logfile(self.pathname+'/'+self.problem_name+'.log')
        self.logfile.scanTimeSteps()
        for kt in range(len(self.logfile.tvect)):
            self.logfile.readTstep(kt)

    def read_dat(self):
        if self.verbose_level<2:
            print(self.pathname + '/' + self.problem_name + '.dat')
        file = open(self.pathname + '/' + self.problem_name + '.dat')
        lines = file.readlines()

        lcount = 0
        kl = 0
        while kl<len(lines):
            line = lines[kl]
            kl += 1
            try:
                ele_type = line.split()[1]
            except:
                ele_type = ''
            if not line: break
            if 'ELEM' in line[:4]:
                self.nElements = int(line.split()[1])
            elif ele_type=='CNNPL':
                v = line.split()
                inel = []
                for kk in range(2):
                    inel.append(int(v[3+kk]))
                self.num_contacts.append(int(v[0]))
                self.nContacts += 1
                self.cnt.inel.append(inel)
                self.cnt.ps.append(int(v[2]))
                self.cnt.type.append(3)
                self.cnt.parent.append(self.num_beams[-1])
            elif ele_type in ['SXQ4']:
                v = line.split()
                self.num_shells.append(int(v[0]))
                self.nShells += 1
                inel = []
                for kk in range(4):
                    inel.append(int(v[5+kk]))
                self.shell.inel.append(inel)
                self.shell.ps.append(int(v[2]))
            elif ele_type in ['SHQ4']:
                v = line.split()
                self.num_shells.append(int(v[0]))
                self.nShells += 1
                inel = []
                for kk in range(8):
                    inel.append(int(v[5+kk]))
                self.shell.inel.append(inel)
                self.shell.ps.append(int(v[2]))
            elif ele_type=='C_Q4':
                v = line.split()
                self.num_contacts.append(int(v[0]))
                self.nContacts += 1
                inel = []
                for kk in range(len(v)-3):
                    inel.append(int(v[3+kk]))
                self.cnt.inel.append(inel)
                self.cnt.ps.append(int(v[2]))
                self.cnt.type.append(1)
##                self.cnt.parent.append(self.num_shells[-1])
            elif ele_type=='C_L2':
                v = line.split()
                self.num_contacts.append(int(v[0]))
                self.nContacts += 1
                inel = []
                for kk in range(4*(dim-1)):
                    inel.append(int(v[3+kk]))
                self.cnt.inel.append(inel)
                self.cnt.ps.append(int(v[2]))
                self.cnt.type.append(0)
                if self.nBeams:
                    self.cnt.parent.append(self.num_beams[-1])
            elif ele_type=='CB_L2':
                v = line.split()
                inel = []
                for kk in range(4):
                    inel.append(int(v[3+kk]))
                self.num_contacts.append(int(v[0]))
                self.nContacts += 1
                self.cnt.inel.append(inel)
                self.cnt.ps.append(int(v[2]))
                self.cnt.type.append(2)
                self.cnt.parent.append(self.num_beams[-1])
            elif ele_type=='Q4' or ele_type=='Q4ES' or ele_type=='B8' or ele_type=='B8ES':
                v = line.split()
                self.num_volumics.append(int(v[0]))
                self.nVolumics += 1
                inel = []
                for kk in range(4*(dim-1)):
                    inel.append(int(v[5+kk]))
                self.vol.inel.append(inel)
                self.vol.ps.append(int(v[2]))
            elif ele_type=='M_Q4':
                v = line.split()
                self.num_membranes.append(int(v[0]))
                self.nMembranes += 1
                inel = []
                for kk in range(4):
                    inel.append(int(v[5+kk]))
                self.mem.inel.append(inel)
                self.mem.ps.append(int(v[2]))
            elif 'BEL2' in ele_type:
                v = line.split()
                self.num_beams.append(int(v[0]))
                self.nBeams += 1
                inel = []
                for kk in range(2):
                    inel.append(int(v[5+kk]))
                self.beam.dir_nodes.append([int(v[kk]) for kk in [9,10]])
                self.beam.inel.append(inel)
                self.beam.ps.append(int(v[2]))
                parent_ele = self.num_beams[-1]
            elif ele_type=='TRS2':
                v = line.split()
                self.num_trusses.append(int(v[0]))
                self.nTrusses += 1
                inel = []
                for kk in range(2):
                    inel.append(int(v[5+kk]))
                self.truss.inel.append(inel)
                self.truss.ps.append(int(v[2]))
            elif 'NODE' in line[:4]:
                self.nNodes = int(line.split()[1])
                crd = [[] for k in range(dim)]
                for k in range(self.nNodes):
                    line = lines[kl]
                    kl += 1
                    cpos = [[10,30],[30,50],[50,70]]
                    for kk in range(dim):
                        crd[kk].append(float(line[cpos[kk][0]:cpos[kk][1]]))
                    if 'L_REC' in line:
                        line = lines[kl]
                        kl += 1
                self.coords = crd
            elif 'JOB_TYPE' in line:
                self.jobtype = line.split()[1]
                if self.jobtype=='3D':
                    dim = 3
                else:
                    dim = 2
            elif 'PROP ' in line:
                self.nPS = int(line.split()[1])
                line = lines[kl]
                kl += 1
                for kps in range(self.nPS):
                    ps = property_set()
                    self.property_sets.append(ps)
                    ps.number = kps+1
                    v = line.split()
                    if not 'DUMMY' in line:
                        ps.mat = int(v[4])
                    ps.EF = int(v[2])-1
                    ps.LF = int(v[3])
                    ps.matmodel = v[1]

                    lastpos = kl-1
                    line = lines[kl]
                    kl += 1
                    group_header = line.split()[0]
                    while len(line)<2:                        
                        line = lines[kl]
                        kl += 1
                    while not line[1]==' ':
                        if not line[0]==' ' and not line[:5]=='HUMID' and not line[:5]=='CREEP' and not line[0]=='-':
                            kl = lastpos
                            break
                        lastpos = kl-1
                        line = lines[kl]
                        while len(line)<2:
                            line = lines[kl]
                            kl += 1
                        try:
                            v = line.split()
                            ps.data[group_header].extend([float(val) for val in v])
                        except:
                            group_header = line.split()[0]
                        kl += 1
                        while len(line)<2:
                            line = lines[kl]
                            kl += 1
            elif 'EXISTFUN' in line:
                self.nEF = int(line.split()[1])
                for kef in range(self.nEF):
                    line = lines[kl+1]
                    kl += 2
                    v = line.split()
                    self.EF.append([float(v[0]),float(v[1])])
            elif 'ACCEL_GLOB' in line:
                self.accel_LTF = int(line.split()[1])
                line = lines[kl]
                kl += 1
                self.accel_mult = float(line.split()[0])
            elif 'PILES' in line:
                self.nPiles = int(line.split()[1])
                for kp in range(self.nPiles):
                    line = lines[kl]
                    kl += 1
                    if len(self.piles)==kp:
                        aPile = pile()
                        self.piles.append(aPile)
                    else:
                        aPile = self.piles[kp]
                    v = line.split()
                    aPile.nBeams = int(v[1])
                    aPile.cnt0D = int(v[2])
                    for ke in range(aPile.nBeams):
                        line = lines[kl]
                        kl += 1
                        v = line.split()
                        aPile.beams.append(int(v[0]))
                        aPile.cnt1D.append(int(v[1]))
            elif 'NAILS' in line:
                self.nNails = int(line.split()[1])
                for kp in range(self.nNails):
                    line = lines[kl]
                    kl += 1
                    aNail = nail()
                    v = line.split()
                    aNail.nBeams = int(v[1])
##                    aNail.cnt0D = int(v[2])
                    for ke in range(aNail.nBeams):
                        line = lines[kl]
                        kl += 1
                        v = line.split()
                        aNail.beams.append(int(v[0]))
                        aNail.cnt1D.append(int(v[1]))
                    self.nails.append(aNail)
            elif 'ANCHOR_HEADS' in line:
                self.nAnchors = int(line.split()[1])
                for kp in range(self.nAnchors):
                    line = lines[kl]
                    kl += 1
                    anAnchor = anchor()
                    v = line.split()
                    anAnchor.nTrusses = int(v[1])
##                    anAnchor.cnt0D = int(v[2])
                    for ke in range(anAnchor.nTrusses):
                        line = lines[kl]
                        kl += 1
                        v = line.split()
                        anAnchor.trusses.append(int(v[0]))
                        anAnchor.cnt1D.append(int(v[1]))
                    self.anchors.append(anAnchor)
                    
            elif 'UNITS' in line:
                pass
            elif 'LINK' in line:
                self.nNodalLinks = int(line.split()[1])
                for klink in range(self.nNodalLinks):
                    line = lines[kl]
                    kl += 1
                    v = line.split()
                    self.nodallinks.append((int(v[0]),int(v[1]),int(v[3])))
            elif 'LIST' in line:
                nList = int(line.split()[1])
                labels = []
                for klist in range(nList):
                    line = lines[kl]
                    kl += 1
                    nVal = int(line.split()[2])
                    label = line[37:-1]
                    labels.append(label)
                    # rtype defines how to read the values
                    vals = []
                    if int(line.split()[3])==1:
                        for kv in range(int(numpy.ceil(nVal/10.))):
                            line = lines[kl]
                            kl += 1
                            for v in line.split():
                                vals.append(int(v))
                    elif int(line.split()[3])==2:
                        for kv in range(nVal):
                            line = lines[kl]
                            kl += 1
                            v = line.split()
                            vals.append((int(v[0]),int(v[1]),int(v[2])))
                    self.lists.append((vals,label))
                    if 'P' in label and 'N' in label and not '_' in label:
                        kp = int(label.split()[0][1:])-1
                        if len(self.piles)==kp:
                            aPile = pile()
                            self.piles.append(aPile)
                        else:
                            aPile = self.piles[kp]
                        kb = int(label.split()[1][1:])
                        eles = vals
                        if len(self.piles[kp].volumics)==0:
                            self.piles[kp].volumics.append([eles[0]])
                        else:
                            if not eles[0]==self.piles[kp].volumics[-1][-1]:
                                self.piles[kp].volumics[-1].append(eles[0])
                            self.piles[kp].volumics.append([eles[-1]])
                        

        file.close()

        for ele in range(self.nVolumics):
            self.vol.mat.append(self.property_sets[self.vol.ps[ele]-1].mat)
            self.vol.EF.append(self.property_sets[self.vol.ps[ele]-1].EF)
            self.vol.LF.append(self.property_sets[self.vol.ps[ele]-1].LF)
        for ele in range(self.nShells):
            self.shell.mat.append(self.property_sets[self.shell.ps[ele]-1].mat)
            self.shell.EF.append(self.property_sets[self.shell.ps[ele]-1].EF)
            self.shell.LF.append(self.property_sets[self.shell.ps[ele]-1].LF)
        for ele in range(self.nBeams):
            self.beam.mat.append(self.property_sets[self.beam.ps[ele]-1].mat)
            self.beam.EF.append(self.property_sets[self.beam.ps[ele]-1].EF)
            self.beam.LF.append(self.property_sets[self.beam.ps[ele]-1].LF)
        for ele in range(self.nTrusses):
            self.truss.mat.append(self.property_sets[self.truss.ps[ele]-1].mat)
            self.truss.EF.append(self.property_sets[self.truss.ps[ele]-1].EF)
            self.truss.LF.append(self.property_sets[self.truss.ps[ele]-1].LF)
        for ele in range(self.nContacts):
            self.cnt.mat.append(self.property_sets[self.cnt.ps[ele]-1].mat)
            self.cnt.EF.append(self.property_sets[self.cnt.ps[ele]-1].EF)
            self.cnt.LF.append(self.property_sets[self.cnt.ps[ele]-1].LF)
        for ele in range(self.nMembranes):
            self.mem.mat.append(self.property_sets[self.mem.ps[ele]-1].mat)
            self.mem.EF.append(self.property_sets[self.mem.ps[ele]-1].EF)
            self.mem.LF.append(self.property_sets[self.mem.ps[ele]-1].LF)

        

    def read_rcf(self):
        if self.verbose_level<2:
            print(self.pathname + '/' + self.problem_name + '.rcf')
        file = open(self.pathname + '/' + self.problem_name + '.rcf')

        line = 'dummy'
        lcount = 0
        while 1:
            if not line: break
            line = file.readline()
            lcount += 1

            if lcount==1:
                self.nGroups = int(line)

            # if all element groups are read, read nodal records:
            elif len(self.ele_groups)==self.nGroups:
                ngroup = nodal_res()
                self.nNodalRes = int(line)
                for kr in range(self.nNodalRes):
                    line = file.readline()
                    lcount += 1
                    ngroup.res_labels.append(line.split()[0])
                    ngroup.ncomp.append(int(line.split()[1]))
                    comps = []
                    if ngroup.ncomp[-1]>1:
                        for kc in range(ngroup.ncomp[-1]):
                            comps.append(line.split()[kc+2])
                    ngroup.comp_labels.append(comps)
                n = int(file.readline())
                while n>0:
                    self.nNodalRes += n
                    line = file.readline()
                    ngroup.res_labels.append(line.split()[0])
                    ngroup.ncomp.append(int(line.split()[1]))
                    comps = []
                    if ngroup.ncomp[-1]>1:
                        for kc in range(ngroup.ncomp[-1]):
                            comps.append(line.split()[kc+2])
                    ngroup.comp_labels.append(comps)
                    n = int(file.readline())
                self.nodal_res.append(ngroup)
                return 1

            # read element groups:
            else:
                # new element group:
                egroup = ele_group()
                self.ele_groups.append(egroup)

                egroup.type = line.split()[0]
                self.ele_group_labels.append(line.split()[0])
                self.nres_types = int(line.split()[1])
                # read all records for this group:
                for kt in range(self.nres_types):
                    line = file.readline()
                    lcount += 1
                    egroup.res_labels.append(line.split()[0])
                    # reading NS-CONTACT is not yet implemented:
                    if egroup.res_labels[-1]=='ADDOUT':
                        egroup.ncomp.append(0)
                    else:
                        egroup.ncomp.append(int(line.split()[1]))
                    comps = []
##                    if egroup.ncomp[-1]>1:
                    if len(line.split())>2:
                        for kc in range(egroup.ncomp[-1]):
                            comps.append(line.split()[kc+2])
                    egroup.comp_labels.append(comps)

        file.close()

    def read_lay(self):
        if self.verbose_level<2:
            print(self.pathname + '/' + self.problem_name + '.lay')
        file = open(self.pathname + '/' + self.problem_name + '.lay')

        line = 'dummy'
        lcount = 0
        while 1:
            line = file.readline()
            if not line: break
            lcount += 1

            if lcount==1:
                self.nGroups_nl = int(line)

            else:
                # new element group:
                egroup = ele_group_nl()
                self.ele_groups_nl.append(egroup)

                egroup.type = line.split()[0]
                self.ele_group_labels_nl.append(line.split()[0])
                self.nres_types_nl = int(line.split()[1])
                # read all records for this group:
                for kt in range(self.nres_types_nl):
                    line = file.readline()
                    lcount += 1
                    egroup.res_labels.append(line.split()[0])
                    egroup.ncomp.append(int(line.split()[1]))
                    comps = []
                    if len(line.split())>2:
                        for kc in range(egroup.ncomp[-1]):
                            comps.append(line.split()[kc+2])
                    egroup.comp_labels.append(comps)
##                return 1

        file.close()

    def read_s00(self,arg=[],res_type='displacements'):
        if '/v+' in arg or self.verbose_level==1:
            verbose = False
            printN = 100
        elif '/v' in arg:
            if len(arg)==2:
                verbose = True
                printN = 0
            elif arg[:2]=='/v':
                verbose = False
                printN = int(arg[2:])
        elif self.verbose_level<2:
            verbose = False
            printN = 1
        else:
            verbose = False
            printN = 1
            
        if res_type=='displacements':
            # read nodal results:
            if self.verbose_level<2:
                print('reading nodal results:')
            f = open(self.pathname + '/' + self.problem_name + '.s00', "rb")
        elif res_type=='reactions':
            # read nodal reactions:
            if self.verbose_level<2:
                print('reading nodal reactions:')
            f = open(self.pathname + '/' + self.problem_name + '.r00', "rb")
        elif res_type=='accelerations':
            # read nodal accelerations:
            if self.verbose_level<2:
                print('reading nodal accelerations:')
            f = open(self.pathname + '/' + self.problem_name + '.a00', "rb")
        elif res_type=='velocities':
            # read nodal velocities:
            if self.verbose_level<2:
                print('reading nodal velocities:')
            f = open(self.pathname + '/' + self.problem_name + '.v00', "rb")
        elif res_type=='init':
            # read initial displacements:
            if self.verbose_level<2:
                print('reading initial displacements:')
            f = open(self.pathname + '/' + self.problem_name + '.uro', "rb")

        kkt = -1
        
        ncomp = 0
        for rt in self.nodal_res[0].res_labels:
            if res_type in ['displacements','reactions','init']:
                if rt=='DISP_TRA':
                    ncomp += self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]
                elif rt=='DISP_ROT':
                    ncomp += self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]
                elif rt=='PPRESS':
                    ncomp += self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]
                elif rt=='PRES_HEAD':
                    ncomp += self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]
                elif rt=='TEMP':
                    ncomp += self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]
            elif res_type=='velocities':
                if rt=='VELOC-TRA':
                    ncomp += self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]
            elif res_type=='accelerations':
                if rt=='ACCEL-TRA':
                    ncomp += self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]

        if res_type in ['init']:
            read_steps = [0]
        else:
            read_steps = self.out_steps
        for kt in range(len(self.steps)):
            if kt in read_steps:
                nbytes = 0
                kkt += 1
                if not (verbose==True or self.verbose_level>=2):
                    if kkt%printN==0:
                        print('reading step ' + str(kt+1) + ' out of ' + str(self.nSteps))
                step = self.give_time_step(kt)
                disp = []
                rot = []
                ppres = []
                pres_head = []
                temp = []
        
                nodal_res = []
                for rt in self.nodal_res[0].res_labels:
                    if res_type in ['displacements','reactions','init']:
                        if rt=='DISP_TRA':
                            nodal_res.append(0)
                            for kcomp in range(self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]):
                                disp.append([])
                        elif rt=='DISP_ROT':
                            nodal_res.append(1)
                            for kcomp in range(self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]):
                                rot.append([])
                        elif rt=='PPRESS':
                            nodal_res.append(2)
                            for kcomp in range(self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]):
                                ppres.append([])
                        elif rt=='PRES_HEAD':
                            nodal_res.append(3)
                            for kcomp in range(self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]):
                                pres_head.append([])
                        elif rt=='TEMP':
                            nodal_res.append(4)
                            for kcomp in range(self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]):
                                temp.append([])
                    elif res_type=='velocities':
                        if rt=='VELOC-TRA':
                            nodal_res.append(0)
                            for kcomp in range(self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]):
                                disp.append([])
                    elif res_type=='accelerations':
                        if rt=='ACCEL-TRA':
                            nodal_res.append(0)
                            for kcomp in range(self.nodal_res[0].ncomp[self.nodal_res[0].res_labels.index(rt)]):
                                disp.append([])
                        
                for n in range(self.nNodes):
                    bytes_read = f.read(4*ncomp)
                    vals = unpack_from('f'*ncomp,bytes_read)
                    ind = 0
                    for kk,k in enumerate(nodal_res):
                        if k==0:
                            for kr in range(self.nodal_res[0].ncomp[kk]):
                                disp[kr].append(vals[ind])
                                ind += 1
                        elif k==1:
                            for kr in range(self.nodal_res[0].ncomp[kk]):
                                rot[kr].append(vals[ind])
                                ind += 1
                        elif k==2:
                            for kr in range(self.nodal_res[0].ncomp[kk]):
                                ppres[kr].append(vals[ind])
                                ind += 1
                        elif k==3:
                            for kr in range(self.nodal_res[0].ncomp[kk]):
                                pres_head[kr].append(vals[ind])
                                ind += 1
                        elif k==4:
                            for kr in range(self.nodal_res[0].ncomp[kk]):
                                temp[kr].append(vals[ind])
                                ind += 1
                if res_type=='displacements':
                    step.nodal.disp = disp
                    step.nodal.rot = rot
                    step.nodal.ppres = ppres
                    step.nodal.pres_head = pres_head
                    step.nodal.temp = temp
                elif res_type=='reactions':
                    step.nodal.r_disp = disp
                    step.nodal.r_rot = rot
                    step.nodal.r_ppres = ppres
                    step.nodal.r_pres_head = pres_head
                elif res_type=='accelerations':
                    step.nodal.a_disp = disp
##                    step.nodal.a_rot = rot
                elif res_type=='velocities':
                    step.nodal.v_disp = disp
##                    step.nodal.v_rot = rot
                elif res_type=='init':
                    step.nodal.i_disp = disp
                    step.nodal.i_rot = rot
                    step.nodal.i_ppres = ppres
                    step.nodal.i_pres_head = pres_head
            else:
                offset = ncomp*self.nNodes*4
                f.seek(offset,1)
        f.close()
        self.nodalRead = True

    def read_s01(self,arg=[]):
        if not 'VOLUMICS' in self.ele_group_labels:
            print('No volumics to be read.')
            return 0

        if '/v+' in arg:
            verbose = False
            printN = 100
        elif '/v' in arg:
            if len(arg)==2:
                verbose = True
                printN = 0
            elif arg[:2]=='/v':
                verbose = False
                printN = int(arg[2:])
        else:
            verbose = False
            printN = 1
        # read VOLUMICS:
        if self.verbose_level<2:
            print('reading volumic results:')
        gind = self.ele_group_labels.index('VOLUMICS')

        egroup = self.ele_groups[gind]
        nint = egroup.ncomp[egroup.res_labels.index('NINT')]
        if self.jobtype=='PLANESTRAIN':
            scomp = 4
            fcomp = 2
        elif self.jobtype=='3D':
            scomp = 6
            fcomp = 3
        

        f = open(self.pathname + '/' + self.problem_name + '.s01', "rb")
        k = 0

        kkt = -1
        byte = 'dummy'
        for kt in range(len(self.steps)):
            if kt in self.out_steps:
                kkt += 1
                if not (verbose==True or self.verbose_level>=2):
                    if kkt%printN==0:
                        print('reading step ' + str(kt+1) + ' out of ' + str(self.nSteps))
                step = self.give_time_step(kt)

                # initialize result vectors for step kt:
                vol_res = [[],[]]
                for rt in egroup.res_labels:
                    krt = egroup.res_labels.index(rt)
                    if rt=='NINT':
                        pass
                    elif rt=='STRESESS' or rt=='STRESSES':
                        stress = []
                        for c in egroup.comp_labels[krt]:
                            stress.append([])
                    elif rt=='STRAINS':
                        strain = []
                        for c in egroup.comp_labels[krt]:
                            strain.append([])
                    elif rt=='PLA_CODE':
                        pla_code = []
                    elif rt=='STR_LEVEL':
                        str_level = []
                    elif rt=='SATUR':
                        satur = []
                    elif rt=='FLU_VELOC':
                        flu_veloc = []
                        for c in egroup.comp_labels[krt]:
                            flu_veloc.append([])
                    else:  # general method for any result type (could replace lines above)
                        vol_res[0].append(rt)
                        if egroup.ncomp[krt]==1:
                            ares = []
                        else:
                            ares = [[] for c in egroup.comp_labels[krt]]
                        vol_res[1].append(ares)
                        

                # read results:
                ncomp = sum(egroup.ncomp)
                for n in range(self.nVolumics):
                    bytes_read = f.read(4*ncomp)
                    if bytes_read=='':break
                    vals = unpack_from('f'*ncomp,bytes_read)
                    ind = 0
                    for rt in egroup.res_labels:
                        krt = egroup.res_labels.index(rt)
                        if rt=='NINT':
                            ind += 1
                        elif rt=='STRESESS' or rt=='STRESSES':
                            for k in range(egroup.ncomp[krt]):
                                stress[k].append(vals[ind])
                                ind += 1
                        elif rt=='STRAINS':
                            for k in range(egroup.ncomp[krt]):
                                strain[k].append(vals[ind])
                                ind += 1
                        elif rt=='PLA_CODE':
                            pla_code.append(vals[ind])
                            ind += 1
                        elif rt=='STR_LEVEL':
                            str_level.append(vals[ind])
                            ind += 1
                        elif rt=='SATUR':
                            satur.append(vals[ind])
                            ind += 1
                        elif rt=='FLU_VELOC':
                            for k in range(egroup.ncomp[krt]):
                                flu_veloc[k].append(vals[ind])
                                ind += 1
                        else:  # general method for any result type (could replace lines above)
                            res_ind = vol_res[0].index(rt)
                            if egroup.ncomp[krt]==1:
                                vol_res[1][res_ind].append(vals[ind])
                                ind += 1
                            else:
                                for k in range(egroup.ncomp[krt]):
                                    vol_res[1][res_ind][k].append(vals[ind])
                                    ind += 1
                if 'STRAINS' in egroup.res_labels:
                    step.vol.strain = strain
                if 'STRESESS' in egroup.res_labels or 'STRESSES' in egroup.res_labels:
                    step.vol.stress = stress
                if 'PLA_CODE' in egroup.res_labels:
                    step.vol.pla_code = pla_code
                if 'STR_LEVEL' in egroup.res_labels:
                    step.vol.str_level = str_level
                if 'SATUR' in egroup.res_labels:
                    step.vol.satur = satur
                if 'FLU_VELOC' in egroup.res_labels:
                    step.vol.flu_veloc = flu_veloc
                for kk in range(len(vol_res[0])):
                    setattr(step.vol,vol_res[0][kk],vol_res[1][kk])
        #        res.steps.append(step)
            else:
                offset = sum(egroup.ncomp)*self.nVolumics*4
                f.seek(offset,1)
                    

        f.close()
        self.volumicsRead = True
        
    def read_s02(self,arg=[]):
        if not 'SHELLS' in self.ele_group_labels:
            print('No shells to be read.')
            return 0


        # compute local coordinate systems:
        for kele in range(self.nShells):
            inel = self.shell.inel[kele]
            crd = []
            for kk in range(len(inel)):
                crd.append([self.coords[0][inel[kk]-1],
                            self.coords[1][inel[kk]-1],
                            self.coords[2][inel[kk]-1]])
            
            exi = numpy.array([0.5*(crd[1][0]+crd[2][0])-0.5*(crd[0][0]+crd[3][0]),
                               0.5*(crd[1][1]+crd[2][1])-0.5*(crd[0][1]+crd[3][1]),
                               0.5*(crd[1][2]+crd[2][2])-0.5*(crd[0][2]+crd[3][2])])
            exi /= la.norm(exi)
            eta = numpy.array([0.5*(crd[2][0]+crd[3][0])-0.5*(crd[0][0]+crd[1][0]),
                               0.5*(crd[2][1]+crd[3][1])-0.5*(crd[0][1]+crd[1][1]),
                               0.5*(crd[2][2]+crd[3][2])-0.5*(crd[0][2]+crd[1][2])])
            eta /= la.norm(eta)
            e3 = numpy.cross(exi,eta)
            e3 /= la.norm(e3)
            
            a = exi+eta
            a /= la.norm(a)
            b = numpy.cross(e3,a)
            b /= la.norm(b)
            
            f = math.sqrt(2.)/2
            e1 = f*(a-b)
            e2 = f*(a+b)

            v = numpy.array(self.ref_vect2D)
            vp = v - numpy.inner(e3,v)*e3
            vpn = la.norm(vp)
            if abs(vpn)<1e-6:
                v = numpy.array([0,0,1])
                vp = v - numpy.inner(e3,v)*e3
                vpn = la.norm(vp)
            cosa = numpy.inner(e1,vp)/vpn
            sina = numpy.inner(e2,vp)/vpn
            
            self.shell.loc_syst.append(numpy.matrix([[cosa,-sina],[sina,cosa]]))
            self.shell.base.append([e1,e2,e3])
        
        if '/v+' in arg:
            verbose = False
            printN = 100
        elif '/v' in arg:
            if len(arg)==2:
                verbose = True
                printN = 0
            elif arg[:2]=='/v':
                verbose = False
                printN = int(arg[2:])
        else:
            verbose = False
            printN = 1
        # read SHELLS:
        if self.verbose_level<2:
            print('reading shell results:')
        gind = self.ele_group_labels.index('SHELLS')

        egroup = self.ele_groups[gind]
        nint = egroup.ncomp[egroup.res_labels.index('NINT')]        

        f = open(self.pathname + '/' + self.problem_name + '.s02', "rb")
        k = 0

        kkt = -1
        byte = 'dummy'
        for kt in range(len(self.steps)):
            if kt in self.out_steps:
                kkt += 1
                if not (verbose==True or self.verbose_level>=2):
                    if kkt%printN==0:
                        print('reading step ' + str(kt+1) + ' out of ' + str(self.nSteps))
                step = self.give_time_step(kt)

                # initialize result vectors for step kt:
                for rt in egroup.res_labels:
                    krt = egroup.res_labels.index(rt)
                    if rt=='THICK':
                        step.shell.thick = []
                    elif rt=='SMFORCE':
                        for c in egroup.comp_labels[krt]:
                            step.shell.smforce.append([])
                    elif rt=='SMOMENT':
                        for c in egroup.comp_labels[krt]:
                            step.shell.smoment.append([])
                    elif rt=='SQFORCE':
                        for c in egroup.comp_labels[krt]:
                            step.shell.sqforce.append([])

                # read results:
                ncomp = sum(egroup.ncomp)
                for n in range(self.nShells):
                    bytes_read = f.read(4*ncomp)
                    if bytes_read=='':break
                    vals = unpack_from('f'*ncomp,bytes_read)
                    ind = 0
                    for rt in egroup.res_labels:
                        krt = egroup.res_labels.index(rt)
                        if rt=='NINT':
                            self.shell.nint.append(int(vals[ind]))
                            ind += 1
                        elif rt=='NLAYER':
                            self.shell.nlayer.append(vals[ind])
                            ind += 1
                        elif rt=='THICK':
                            step.shell.thick.append(vals[ind])
                            ind += 1
                        elif rt=='SMFORCE':
                            f0 = numpy.matrix([[vals[ind],vals[ind+2]],[vals[ind+2],vals[ind+1]]])
                            ind += 3

                            rot = self.shell.loc_syst[n]
                            ff = rot.transpose()*f0
                            ff *= rot
                            step.shell.smforce[0].append(ff[0,0])
                            step.shell.smforce[1].append(ff[1,1])
                            step.shell.smforce[2].append(ff[0,1])
                        elif rt=='SMOMENT':
                            f0 = numpy.matrix([[vals[ind],vals[ind+2]],[vals[ind+2],vals[ind+1]]])
                            ind += 3

                            rot = self.shell.loc_syst[n]
                            ff = rot.transpose()*f0
                            ff *= rot                            
                            step.shell.smoment[0].append(ff[0,0])
                            step.shell.smoment[1].append(ff[1,1])
                            step.shell.smoment[2].append(ff[0,1])
                        elif rt=='SQFORCE':
                            f0 = numpy.matrix([[vals[ind]],[vals[ind+1]]])
                            ind += 3

                            rot = self.shell.loc_syst[n]
                            ff = rot.transpose()*f0
                            for k in range(egroup.ncomp[krt]):
                                step.shell.sqforce[k].append(ff[k,0])
                                ind += 1
##                            for k in range(egroup.ncomp[krt]):
##                                step.shell.sqforce[k].append(vals[ind])
##                                ind += 1
            elif kt==0:
                # read results:
                ncomp = sum(egroup.ncomp)
                for n in range(self.nShells):
                    bytes_read = f.read(4*ncomp)
                    if bytes_read=='':break
                    vals = unpack_from('f'*ncomp,bytes_read)
                    ind = 0
                    for rt in egroup.res_labels:
                        krt = egroup.res_labels.index(rt)
                        if rt=='NINT':
                            self.shell.nint.append(int(vals[ind]))
                            ind += 1
                        elif rt=='NLAYER':
                            self.shell.nlayer.append(vals[ind])
                            ind += 1
                        elif rt=='THICK':
                            ind += 1
                        elif rt=='SMFORCE':
                            ind += 3
                        elif rt=='SMOMENT':
                            ind += 3
                        elif rt=='SQFORCE':
                            ind += 3
                            for k in range(egroup.ncomp[krt]):
                                ind += 1                
            else:
                offset = sum(egroup.ncomp)*self.nShells*4
                f.seek(offset,1)
                

        f.close()
        self.shellsRead = True

    def read_L02(self,arg=[]):
        if not 'SHELLS' in self.ele_group_labels:
            if self.verbose_level<3:
                print('No shells to be read.')
            return 0

        
        if '/v+' in arg:
            verbose = False
            printN = 100
        elif '/v' in arg:
            if len(arg)==2:
                verbose = True
                printN = 0
            elif arg[:2]=='/v':
                verbose = False
                printN = int(arg[2:])
        else:
            verbose = False
            printN = 1
        # read SHELLS:
        if self.verbose_level<2:
            print('reading non-linear shell results:')
        gind = self.ele_group_labels_nl.index('SHELLS')

        for kt in range(len(self.steps)):
            step = self.give_time_step(kt)
            nint = 1#step.shell.nint[0]
            if len(self.shell.nlayer):
                nTotLayers = 0
                for ke in range(self.nShells):
                    nLayers = int(self.shell.nlayer[ke])
                    nTotLayers += nLayers*nint


        egroup = self.ele_groups_nl[gind]     

        f = open(self.pathname + '/' + self.problem_name + '.L02', "rb")
        k = 0

        kkt = -1
        byte = 'dummy'
        for kt in range(len(self.steps)):
            if kt in self.out_steps:
                kkt += 1
                if not (verbose==True or self.verbose_level>=2):
                    if kkt%printN==0:
                        print('reading step ' + str(kt+1) + ' out of ' + str(self.nSteps))
                step = self.give_time_step(kt)
                if len(self.shell.nlayer)==0 and self.nShells>0:
                    print('shell results (s02) have to be read first')
                nint = self.shell.nint[0]

                # initialize result vectors for step kt:
                for rt in egroup.res_labels:
                    krt = egroup.res_labels.index(rt)
                    if rt=='POSITION':
                        step.nlshell.position = []
                    elif rt=='STRESESS' or rt=='STRESSES':
                        step.nlshell.stress = []
                    elif rt=='STRAINS':
                        step.nlshell.strain = []
                    elif rt=='PLA_CODE':
                        step.nlshell.pla_code = []
                    elif rt=='STR_LEVEL':
                        step.nlshell.str_level = []
                    elif rt=='DAMAGE':
                        step.nlshell.damage = []

                # read results:
                for kgp in range(nint):
                    for ke in range(self.nShells):
                        nLayers = int(self.shell.nlayer[ke])
                        position = []
##                        position = [[] for c in egroup.comp_labels[egroup.res_labels.index('POSITION')]]
                        if 'STRESESS' in egroup.res_labels:
                            stress = [[] for c in egroup.comp_labels[egroup.res_labels.index('STRESESS')]]
                        else:
                            stress = [[] for c in egroup.comp_labels[egroup.res_labels.index('STRESSES')]]
                        strain = [[] for c in egroup.comp_labels[egroup.res_labels.index('STRAINS')]]
                        pla_code = []
                        str_level = []
                        damage = []
                        for kl in range(nLayers):
                            for rt in egroup.res_labels:
                                krt = egroup.res_labels.index(rt)
                                ncomp = len(egroup.comp_labels[krt])
                                if rt=='POSITION':
                                    byte = f.read(4)
                                    position.append(struct.unpack('f',byte)[0])
                                elif rt=='STRESESS' or rt=='STRESSES':
                                    for k in range(ncomp):
                                        byte = f.read(4)
##                                    print('stress',struct.unpack('f',byte)[0])
                                        stress[k].append(struct.unpack('f',byte)[0])
                                elif rt=='STRAINS':
                                    for k in range(ncomp):
                                        byte = f.read(4)
##                                    print('strain',struct.unpack('f',byte)[0])
                                        strain[k].append(struct.unpack('f',byte)[0])
                                elif rt=='PLA_CODE':
                                    byte = f.read(4)
##                                    print('pla_code',struct.unpack('f',byte)[0])
                                    pla_code.append(struct.unpack('f',byte)[0])
                                elif rt=='STR_LEVEL':
                                    byte = f.read(4)
##                                    print('str_level',struct.unpack('f',byte)[0])
                                    str_level.append(struct.unpack('f',byte)[0])
                                elif rt=='DAMAGE':
                                    byte = f.read(4)
                                    damage.append(struct.unpack('f',byte)[0])
                        step.nlshell.position.append(position)
                        step.nlshell.stress.append(stress)
                        step.nlshell.strain.append(strain)
                        step.nlshell.pla_code.append(pla_code)
                        step.nlshell.str_level.append(str_level)
                        step.nlshell.damage.append(damage)

            else:
                offset = sum(egroup.ncomp)*nTotLayers*4
                f.seek(offset,1)
                

        f.close()
        self.shellsNLRead = True
        
    def read_s03(self,arg=[]):
        if not 'TRUSSES' in self.ele_group_labels:
            if self.verbose_level<3:
                print('No trusses to be read.')
            return 0

        if '/v+' in arg:
            verbose = False
            printN = 100
        elif '/v' in arg:
            if len(arg)==2:
                verbose = True
                printN = 0
            elif arg[:2]=='/v':
                verbose = False
                printN = int(arg[2:])
        else:
            verbose = False
            printN = 1
        # read Trusses:
        if self.verbose_level<2:
            print('reading truss results:')
        gind = self.ele_group_labels.index('TRUSSES')

        egroup = self.ele_groups[gind]
        nint = egroup.ncomp[egroup.res_labels.index('NINT')]        

        f = open(self.pathname + '/' + self.problem_name + '.s03', "rb")
        k = 0

        kkt = -1
        byte = 'dummy'
        for kt in range(len(self.steps)):
            if kt in self.out_steps:
                kkt += 1
                if not (verbose==True or self.verbose_level>=2):
                    if kkt%printN==0:
                        print('reading step ' + str(kt+1) + ' out of ' + str(self.nSteps))
                step = self.give_time_step(kt)

                # initialize result vectors for step kt:
                for rt in egroup.res_labels:
                    krt = egroup.res_labels.index(rt)
                    if rt=='STRESESS' or rt=='STRESSES':
##                        for c in egroup.comp_labels[krt]:
##                            step.truss.stress.append([])
                        step.truss.stress = []
                    elif rt=='STRAINS':
##                        for c in egroup.comp_labels[krt]:
##                            step.truss.strain.append([])
                        step.truss.strain = []
                    elif rt=='PLA_CODE':
                        step.truss.pla_code = []
                    elif rt=='STR_LEVEL':
                        step.truss.str_level = []
                    elif rt=='FORCE':
                        step.truss.force = []

                # read results:
                for ke in range(self.nTrusses):
                    for rt in egroup.res_labels:
                        krt = egroup.res_labels.index(rt)
                        byte = f.read(4)
                        if rt=='STRESESS' or rt=='STRESSES':
                            step.truss.stress.append(struct.unpack('f',byte)[0])
                        elif rt=='STRAINS':
                            step.truss.strain.append(struct.unpack('f',byte)[0])
                        elif rt=='PLA_CODE':
                            step.truss.pla_code.append(struct.unpack('f',byte)[0])
                        elif rt=='STR_LEVEL':
                            step.truss.str_level.append(struct.unpack('f',byte)[0])
                        elif rt=='FORCE':
                            step.truss.force.append(struct.unpack('f',byte)[0])
        #        res.steps.append(step)
            else:
                offset = sum(egroup.ncomp)*self.nTrusses*4
                f.seek(offset,1)

        f.close()
        self.trussesRead = True
        
    def read_s04(self,arg=[]):
        if not 'BEAMS' in self.ele_group_labels:
            if self.verbose_level<3:
                print('No beams to be read.')
            return 0

        if '/v+' in arg:
            verbose = False
            printN = 100
        elif '/v' in arg:
            if len(arg)==2:
                verbose = True
                printN = 0
            elif arg[:2]=='/v':
                verbose = False
                printN = int(arg[2:])
        else:
            verbose = False
            printN = 1
        # read Beams:
        if self.verbose_level<2:
            print('reading beam results:')
        gind = self.ele_group_labels.index('BEAMS')

        egroup = self.ele_groups[gind]

        f = open(self.pathname + '/' + self.problem_name + '.s04', "rb")
        k = 0

        kkt = -1
        byte = 'dummy'
        for kt in range(len(self.steps)):
            if kt in self.out_steps:
                kkt += 1
                if not (verbose==True or self.verbose_level>=2):
                    if kkt%printN==0:
                        print('reading step ' + str(kt+1) + ' out of ' + str(self.nSteps))
                step = self.give_time_step(kt)

                # initialize result vectors for step kt:
                for rt in egroup.res_labels:
                    krt = egroup.res_labels.index(rt)
                    if rt=='NLAYER':
                        step.beam.nlayer = []
                    elif rt=='FORCE':
                        for c in egroup.comp_labels[krt]:
                            step.beam.force.append([])
                    elif rt=='MOMENT':
                        for c in egroup.comp_labels[krt]:
                            step.beam.moment.append([])
                        step.beam.dm = []
                    elif rt=='EDISP_TRA':
                        for c in egroup.comp_labels[krt]:
                            step.beam.edisp_tra.append([])
                    elif rt=='EDISP_ROT':
                        for c in egroup.comp_labels[krt]:
                            step.beam.edisp_rot.append([])

                # read results:
                nint = 1
                kgp = 0
                while kgp < nint:
                    for ke in range(self.nBeams):
                        inel = self.beam.inel[ke]
                        dx = self.coords[0][inel[1]-1]-self.coords[0][inel[0]-1]
                        dy = self.coords[1][inel[1]-1]-self.coords[1][inel[0]-1]
                        dl = 0.5*math.sqrt(dx**2+dy**2)
                        for rt in egroup.res_labels:
                            krt = egroup.res_labels.index(rt)
                            ncomp = len(egroup.comp_labels[krt])
                            if rt=='NINT':
                                byte = f.read(4)
                                nint = int(struct.unpack('f',byte)[0])
                                self.beam.nint.append(nint)
                            elif rt=='NLAYER':
                                byte = f.read(4)
                                self.beam.nlayer.append(struct.unpack('f',byte)[0])
                            elif rt=='FORCE':
                                for k in range(len(egroup.comp_labels[krt])):
                                    byte = f.read(4)
                                    step.beam.force[k].append(struct.unpack('f',byte)[0])
                            elif rt=='MOMENT':
                                for k in range(len(egroup.comp_labels[krt])):
                                    byte = f.read(4)
                                    step.beam.moment[k].append(struct.unpack('f',byte)[0])
                                step.beam.dm.append(step.beam.force[1][ke]*dl)
                            elif rt=='EDISP_TRA':
                                for k in range(len(egroup.comp_labels[krt])):
                                    byte = f.read(4)
                                    step.beam.edisp_tra[k].append(struct.unpack('f',byte)[0])
                            elif rt=='EDISP_ROT':
                                for k in range(len(egroup.comp_labels[krt])):
                                    byte = f.read(4)
                                    step.beam.edisp_rot[k].append(struct.unpack('f',byte)[0])
                    kgp += 1
        #        res.steps.append(step)
            elif kt==0:
                step = self.give_time_step(kt)

                # read results:
                nint = 1
                kgp = 0
                while kgp < nint:
                    for ke in range(self.nBeams):
                        for rt in egroup.res_labels:
                            krt = egroup.res_labels.index(rt)
                            ncomp = len(egroup.comp_labels[krt])
                            if rt=='NINT':
                                byte = f.read(4)
                                nint = int(struct.unpack('f',byte)[0])
                                self.beam.nint.append(nint)
                            elif rt=='NLAYER':
                                byte = f.read(4)
                                nlayer = int(struct.unpack('f',byte)[0])
                                self.beam.nlayer.append(nlayer)
                            elif rt=='FORCE':
                                for k in range(len(egroup.comp_labels[krt])):
                                    byte = f.read(4)
                            elif rt=='MOMENT':
                                for k in range(len(egroup.comp_labels[krt])):
                                    byte = f.read(4)
                            elif rt=='EDISP_TRA':
                                for k in range(len(egroup.comp_labels[krt])):
                                    byte = f.read(4)
                            elif rt=='EDISP_ROT':
                                for k in range(len(egroup.comp_labels[krt])):
                                    byte = f.read(4)
                    kgp += 1
                offset = sum(egroup.ncomp)*4*sum([self.beam.nint[ke] for ke in range(self.nBeams)])
            else:
##                step = self.give_time_step(kt)
##                # read results:
##                nint = 1
##                kgp = 0
##                while kgp < nint:
##                    for ke in range(self.nBeams):
##                        for rt in egroup.res_labels:
##                            krt = egroup.res_labels.index(rt)
##                            offset = 0
##                            if rt=='NINT':
##                                byte = f.read(4)
##                                nint = int(struct.unpack('f',byte)[0])
##                                step.beam.nint.append(nint)
##                            elif rt=='NLAYER':
##                                byte = f.read(4)
##                                step.beam.nlayer.append(struct.unpack('f',byte)[0])
##                            elif rt=='FORCE':
##                                for k in range(len(egroup.comp_labels[krt])):
##                                    offset += 4
##                            elif rt=='MOMENT':
##                                for k in range(len(egroup.comp_labels[krt])):
##                                    offset += 4
##                            elif rt=='EDISP_TRA':
##                                for k in range(len(egroup.comp_labels[krt])):
##                                    offset += 4
##                            elif rt=='EDISP_ROT':
##                                for k in range(len(egroup.comp_labels[krt])):
##                                    offset += 4
##                            f.seek(offset,1)
##                    kgp += 1
                offset = sum(egroup.ncomp)*self.nBeams*4
                f.seek(offset,1)


        f.close()
        self.beamsRead = True

    def read_L04(self,arg=[],minmax=False):
        if not 'BEAMS' in self.ele_group_labels_nl:
            if self.verbose_level<3:
                print('No beams to be read.')
            return 0

        maxmin = False
        if '/v+' in arg:
            verbose = False
            printN = 100
        elif '/v' in arg:
            if len(arg)==2:
                verbose = True
                printN = 0
            elif arg[:2]=='/v':
                verbose = False
                printN = int(arg[2:])
        else:
            verbose = False
            printN = 1
        # read Beams:
        if self.verbose_level<2:
            print('reading non-linear beam results:')
        gind = self.ele_group_labels_nl.index('BEAMS')

        for kt in self.out_steps:
            step = self.give_time_step(kt)
            if len(step.beam.force):
                nTotLayers = 0
                for ke in range(self.nBeams):
                    nint = self.beam.nint[ke]
                    nLayers = int(self.beam.nlayer[ke])
                    nTotLayers += nLayers*nint

        egroup = self.ele_groups_nl[gind]

        f = open(self.pathname + '/' + self.problem_name + '.L04', "rb")
        k = 0

        kkt = -1
        byte = 'dummy'
        for kt in range(len(self.steps)):
            if kt in self.out_steps:
                kkt += 1
                if not (verbose==True or self.verbose_level>=2):
                    if kkt%printN==0:
                        print('reading step ' + str(kt+1) + ' out of ' + str(self.nSteps))
                step = self.give_time_step(kt)
                if len(self.beam.nlayer)==0 and self.nBeams>0:
                    print('beam results (s04) have to be read first')
                nint = self.beam.nint[0]

                # initialize result vectors for step kt:
                for rt in egroup.res_labels:
                    krt = egroup.res_labels.index(rt)
                    if rt=='POSITION':
                        step.nlbeam.position = []
                    elif rt=='STRESESS' or rt=='STRESSES':
                        step.nlbeam.stress = []
                    elif rt=='STRAINS':
                        step.nlbeam.strain = []
                    elif rt=='PLA_CODE':
                        step.nlbeam.pla_code = []
                    elif rt=='STR_LEVEL':
                        step.nlbeam.str_level = []

                # read results:
                for kgp in range(nint):
                    for ke in range(self.nBeams):
                        position = [[] for c in egroup.comp_labels[egroup.res_labels.index('POSITION')]]
                        stress = []
                        strain = []
                        pla_code = []
                        str_level = []
                        nLayers = int(self.beam.nlayer[ke])
                        for kl in range(nLayers):
                            for rt in egroup.res_labels:
                                krt = egroup.res_labels.index(rt)
                                ncomp = len(egroup.comp_labels[krt])
                                if rt=='POSITION':
                                    for k in range(len(egroup.comp_labels[krt])):
                                        byte = f.read(4)
                                        position[k].append(struct.unpack('f',byte)[0])
                                elif rt=='STRESESS' or rt=='STRESSES':
                                    byte = f.read(4)
                                    stress.append(struct.unpack('f',byte)[0])
                                elif rt=='STRAINS':
                                    byte = f.read(4)
                                    strain.append(struct.unpack('f',byte)[0])
                                elif rt=='PLA_CODE':
                                    byte = f.read(4)
                                    pla_code.append(struct.unpack('f',byte)[0])
                                elif rt=='STR_LEVEL':
                                    byte = f.read(4)
                                    str_level.append(struct.unpack('f',byte)[0])
                        if minmax:
                            step.nlbeam.position.append(0)
                            temp = [0,0,0,0]
                            temp_ind = []
                            if len(stress)>0:
                                if len(stress)>30:
                                    ## assuming 30 concrete layers
                                    temp = [min(stress[:30]),
                                            max(stress[:30]),
                                            min(stress[30:]),
                                            max(stress[30:])]
                                    temp_ind = [stress.index(v) for v in temp]
                                else:
                                    temp = [min(stress[:30]),
                                            max(stress[:30]),
                                            0,0]
                                    temp_ind = [stress.index(temp[0]),
                                                stress.index(temp[1]),
                                                0,0]
                            else:
                                temp = [0,0,0,0]
                            step.nlbeam.stress.append(temp)
                            step.nlbeam.strain.append([strain[k] for k in temp_ind])
                            step.nlbeam.pla_code.append(0)
                            step.nlbeam.str_level.append(0)
                        else:
                            step.nlbeam.position.append(position)
                            step.nlbeam.stress.append(stress)
                            step.nlbeam.strain.append(strain)
                            step.nlbeam.pla_code.append(pla_code)
                            step.nlbeam.str_level.append(str_level)
            #        res.steps.append(step)
            else:
                offset = sum(egroup.ncomp)*nTotLayers*4
                f.seek(offset,1)

        f.close()
        self.beamsNLRead = True
        
    def read_s07(self,arg=[]):
        if not 'CONTACT' in self.ele_group_labels:
            if self.verbose_level<3:
                print('No contacts to be read.')
            return 0

        if '/v+' in arg:
            verbose = False
            printN = 100
        elif '/v' in arg:
            if len(arg)==2:
                verbose = True
                printN = 0
            elif arg[:2]=='/v':
                verbose = False
                printN = int(arg[2:])
        else:
            verbose = False
            printN = 1

        # compute local coordinate systems:
        Xi = [-1,1,1,-1]
        Eta = [-1,-1,1,1]
        dNdxi = []
        dNdeta = []
        for kgp in range(4):
            xi = Xi[kgp]
            eta = Eta[kgp]
            dNdxi.append([-0.25*(1-eta),0.25*(1-eta),0.25*(1+eta),-0.25*(1+eta)])
            dNdeta.append([-0.25*(1-xi),-0.25*(1+xi),0.25*(1+xi),0.25*(1-xi)])

        for kele in range(self.nContacts):
            if self.cnt.type[kele]==1:
                inel = self.cnt.inel[kele]
                crd = []
                for kk in range(4):
                    crd.append([self.coords[0][inel[kk]-1],
                                self.coords[1][inel[kk]-1],
                                self.coords[2][inel[kk]-1]])

                loc_syst = []
                base = []
                area = []
                for kgp in range(4):
                    dxdxi = sum([dNdxi[kgp][k]*crd[k][0] for k in range(4)])
                    dydxi = sum([dNdxi[kgp][k]*crd[k][1] for k in range(4)])
                    dzdxi = sum([dNdxi[kgp][k]*crd[k][2] for k in range(4)])
                    dxdeta = sum([dNdeta[kgp][k]*crd[k][0] for k in range(4)])
                    dydeta = sum([dNdeta[kgp][k]*crd[k][1] for k in range(4)])
                    dzdeta = sum([dNdeta[kgp][k]*crd[k][2] for k in range(4)])
                    dXdxi = numpy.array([dxdxi,dydxi,dzdxi])
                    dXdeta = numpy.array([dxdeta,dydeta,dzdeta])
                    exi = dXdxi/numpy.linalg.norm(dXdxi)
                    eta = dXdeta/numpy.linalg.norm(dXdeta)
                    area.append(numpy.linalg.norm(numpy.cross(dXdxi,dXdeta)))
                    
        
##                exi = numpy.array([0.5*(crd[1][0]+crd[2][0])-0.5*(crd[0][0]+crd[3][0]),
##                                   0.5*(crd[1][1]+crd[2][1])-0.5*(crd[0][1]+crd[3][1]),
##                                   0.5*(crd[1][2]+crd[2][2])-0.5*(crd[0][2]+crd[3][2])])
##                
##                exi /= la.norm(exi)
##                eta = numpy.array([0.5*(crd[2][0]+crd[3][0])-0.5*(crd[0][0]+crd[1][0]),
##                                   0.5*(crd[2][1]+crd[3][1])-0.5*(crd[0][1]+crd[1][1]),
##                                   0.5*(crd[2][2]+crd[3][2])-0.5*(crd[0][2]+crd[1][2])])
##                eta /= la.norm(eta)
                    e3 = numpy.cross(exi,eta)
                    e3 /= la.norm(e3)
                    
                    a = exi+eta
                    a /= la.norm(a)
                    b = numpy.cross(e3,a)
                    b /= la.norm(b)
                    
                    f = math.sqrt(2.)/2
                    e1 = f*(a-b)
                    e2 = f*(a+b)

                    v = numpy.array(self.ref_vect2D)
                    vp = v - numpy.inner(e3,v)*e3
                    vpn = la.norm(vp)
                    if abs(vpn)<1e-6:
                        v = numpy.array([0,0,1])
                        vp = v - numpy.inner(e3,v)*e3
                        vpn = la.norm(vp)
                    cosa = numpy.inner(e1,vp)/vpn
                    sina = numpy.inner(e2,vp)/vpn

                    loc_syst.append(numpy.matrix([[cosa,-sina],[sina,cosa]]))
                    base.append([e1,e2,e3])
                self.cnt.loc_syst.append(loc_syst)
                self.cnt.base.append(base)
                self.cnt.area.append(area)
##            elif self.cnt.type[kele]==1:
##                inel = self.cnt.inel[kele]
##                crd = []
##                for kk in range(2):
##                    crd.append([self.coords[0][inel[kk]-1],
##                                self.coords[1][inel[kk]-1]])
##                dx = crd[1][0]-crd[0][0]
##                dy = crd[1][1]-crd[0][1]
##                norm = math.sqrt(dx**2+dy**2)
##                self.cnt.loc_syst.append(numpy.matrix([[dx/norm,dy/norm],[-dy/norm,dx/norm]]))
            else:
                self.cnt.loc_syst.append(0)
                self.cnt.base.append(0)
            
        # read Contact:
        if self.verbose_level<2:
            print('reading contact results:')
        gind = self.ele_group_labels.index('CONTACT')

        egroup = self.ele_groups[gind]
        nint = egroup.ncomp[egroup.res_labels.index('NINT')]

        f = open(self.pathname + '/' + self.problem_name + '.s07', "rb")
        k = 0

        nGPdict = {0:2,1:4,2:2,3:1} # 0: 2D, 1: 3D plane contact, 2: pile/anchor, 3: pile tip
        byte = 'dummy'
        kkt = -1
        for kt in range(len(self.steps)):
            if kt in self.out_steps:
                kkt += 1
                if not (verbose==True or self.verbose_level>=2):
                    if kkt%printN==0:
                        print('reading step ' + str(kt+1) + ' out of ' + str(self.nSteps))
                step = self.give_time_step(kt)

                # initialize result vectors for step kt:
                for rt in egroup.res_labels:
                    krt = egroup.res_labels.index(rt)
                    if rt=='NINT':
                        step.cnt.nint = [0 for k in range(self.nContacts)]
                    elif rt=='STRESESS' or rt=='STRESSES':
                        for c in egroup.comp_labels[krt]:
                            step.cnt.stress.append([[0 for kgp in range(nGPdict[self.cnt.type[k]])] for k in range(self.nContacts)])
                    elif rt=='T_STRESSES':
                        for c in egroup.comp_labels[krt]:
                            step.cnt.tstress.append([[0 for kgp in range(nGPdict[self.cnt.type[k]])] for k in range(self.nContacts)])
                    elif rt=='STRAINS':
                        for c in egroup.comp_labels[krt]:
                            step.cnt.strain.append([[0 for kgp in range(nGPdict[self.cnt.type[k]])] for k in range(self.nContacts)])
                    elif rt=='PLA_CODE':
                        step.cnt.pla_code = [[0 for kgp in range(nGPdict[self.cnt.type[k]])] for k in range(self.nContacts)]
                    elif rt=='STR_LEVEL':
                        step.cnt.str_level = [[0 for kgp in range(nGPdict[self.cnt.type[k]])] for k in range(self.nContacts)]

                # read results for contact (plane and 1D):
                ncomp = sum(egroup.ncomp)
                for ke in range(self.nContacts):
                    nGP = nGPdict[self.cnt.type[ke]]
                    for kgp in range(nGP):
                        bytes_read = f.read(4*ncomp)
                        if bytes_read=='':break
                        vals = unpack_from('f'*ncomp,bytes_read)
                        ind = 0
                        for rt in egroup.res_labels:
                            krt = egroup.res_labels.index(rt)
                            if rt=='NINT':
                                nint = int(vals[ind])
                                if not nint==nGP:
                                    print('Check reading of contacts, nint=%i'%(nint))
                                ind += 1
                            elif rt=='STRESESS' or rt=='STRESSES':
                                if self.cnt.type[ke]==1:
                                    f0 = numpy.matrix([[vals[ind]],[vals[ind+1]]])

                                    rot = self.cnt.loc_syst[ke][kgp]
                                    ff = rot.transpose()*f0
                                    step.cnt.stress[0][ke][kgp] = ff[0,0]
                                    step.cnt.stress[1][ke][kgp] = ff[1,0]
##                                    step.cnt.stress[0][ke][kgp] = vals[ind+0]#ff[0,0]
##                                    step.cnt.stress[1][ke][kgp] = vals[ind+1]#ff[1,0]
                                    step.cnt.stress[2][ke][kgp] = vals[ind+2]

                                    ind += 3
                                elif self.cnt.type[ke] in [0]:
                                    step.cnt.stress[0][ke][kgp] = vals[ind]
                                    step.cnt.stress[1][ke][kgp] = vals[ind+1]
                                    ind += 2
                                elif self.cnt.type[ke] in [2,3] and self.jobtype=='3D': # pile interf in 3D, needs to be modified for 2D anchors
                                    step.cnt.stress[0][ke][kgp] = vals[ind]
                                    step.cnt.stress[1][ke][kgp] = vals[ind+1]
                                    step.cnt.stress[2][ke][kgp] = vals[ind+2]
                                    ind += 3
                                else:
                                    step.cnt.stress[0][ke][kgp] = vals[ind]
                                    step.cnt.stress[1][ke][kgp] = vals[ind+1]
                                    ind += 2
                            elif rt=='T_STRESSES':
                                if self.cnt.type[ke]==1:
                                    f0 = numpy.matrix([[vals[ind]],[vals[ind+1]]])

                                    rot = self.cnt.loc_syst[ke][kgp]
                                    ff = rot.transpose()*f0
                                    step.cnt.tstress[0][ke][kgp] = ff[0,0]
                                    step.cnt.tstress[1][ke][kgp] = ff[1,0]
                                    step.cnt.tstress[2][ke][kgp] = vals[ind+2]
                                    ind += 3
                                elif self.cnt.type[ke] in [0]:
                                    step.cnt.tstress[0][ke][kgp] = vals[ind]
                                    step.cnt.tstress[1][ke][kgp] = vals[ind+1]
                                    ind += 2
                                elif self.cnt.type[ke] in [2,3] and self.jobtype=='3D':
                                    step.cnt.tstress[0][ke][kgp] = vals[ind]
                                    step.cnt.tstress[1][ke][kgp] = vals[ind+1]
                                    step.cnt.tstress[2][ke][kgp] = vals[ind+2]
                                    ind += 3
                                else:
                                    step.cnt.tstress[0][ke][kgp] = vals[ind]
                                    step.cnt.tstress[1][ke][kgp] = vals[ind+1]
                                    ind += 2
                            elif rt=='STRAINS':
                                if self.cnt.type[ke]==1:
                                    f0 = numpy.matrix([[vals[ind]],[vals[ind+1]]])

                                    rot = self.cnt.loc_syst[ke][kgp]
                                    ff = rot.transpose()*f0
                                    step.cnt.strain[0][ke][kgp] = ff[0,0]
                                    step.cnt.strain[1][ke][kgp] = ff[1,0]
                                    step.cnt.strain[2][ke][kgp] = vals[ind+2]
                                    ind += 3
                                elif self.cnt.type[ke] in [0]:
                                    step.cnt.strain[0][ke][kgp] = vals[ind]
                                    step.cnt.strain[1][ke][kgp] = vals[ind+1]
                                    ind += 2
                                elif self.cnt.type[ke] in [2,3] and self.jobtype=='3D':
                                    step.cnt.strain[0][ke][kgp] = vals[ind]
                                    step.cnt.strain[1][ke][kgp] = vals[ind+1]
                                    step.cnt.strain[2][ke][kgp] = vals[ind+2]
                                    ind += 3
                                else:
                                    step.cnt.strain[0][ke][kgp] = vals[ind]
                                    step.cnt.strain[1][ke][kgp] = vals[ind+1]
                                    ind += 2
                            elif rt=='PLA_CODE':
                                step.cnt.pla_code[ke][kgp] = vals[ind]
                                ind += 1
                            elif rt=='STR_LEVEL':
                                step.cnt.str_level[ke][kgp] = vals[ind]
                                ind += 1
            else:
                offset = sum(egroup.ncomp)*2*self.cnt.type.count(0)*4
                offset += sum(egroup.ncomp)*4*self.cnt.type.count(1)*4
                offset += sum(egroup.ncomp)*2*self.cnt.type.count(2)*4
                offset += sum(egroup.ncomp)*1*self.cnt.type.count(3)*4
                f.seek(offset,1)
                
        f.close()
        self.contactsRead = True


    def read_s15(self,arg=[]):
        if not 'MEMBRANE' in self.ele_group_labels:
            if self.verbose_level<3:
                print('No membranes to be read.')
            return 0

        if '/v+' in arg:
            verbose = False
            printN = 100
        elif '/v' in arg:
            if len(arg)==2:
                verbose = True
                printN = 0
            elif arg[:2]=='/v':
                verbose = False
                printN = int(arg[2:])
        else:
            verbose = False
            printN = 1
        # read MEMBRANES:
        if self.verbose_level<2:
            print('reading membrane results:')
        gind = self.ele_group_labels.index('MEMBRANE')

        egroup = self.ele_groups[gind]
        nint = egroup.ncomp[egroup.res_labels.index('NINT')]
        scomp = 4
        fcomp = 2
        

        f = open(self.pathname + '/' + self.problem_name + '.s15', "rb")
        k = 0

        kkt = -1
        byte = 'dummy'
        for kt in range(len(self.steps)):
            if kt in self.out_steps:
                kkt += 1
                if verbose==False:
                    if kkt%printN==0:
                        print('reading step ' + str(kt+1) + ' out of ' + str(self.nSteps))
                step = self.give_time_step(kt)

                # initialize result vectors for step kt:
                vol_res = [[],[]]
                for rt in egroup.res_labels:
                    krt = egroup.res_labels.index(rt)
                    if rt=='NINT':
                        pass
                    elif rt=='SMFORCE':
                        smforce = []
                        for c in egroup.comp_labels[krt]:
                            smforce.append([])
                    elif rt=='STRAINS':
                        strain = []
                        for c in egroup.comp_labels[krt]:
                            strain.append([])
                    elif rt=='PLA_CODE':
                        pla_code = []
                    elif rt=='STR_LEVEL':
                        str_level = []
                    else:  # general method for any result type (could replace lines above)
                        vol_res[0].append(rt)
                        if egroup.ncomp[krt]==1:
                            ares = []
                        else:
                            ares = [[] for c in egroup.comp_labels[krt]]
                        vol_res[1].append(ares)
                        

                # read results:
                ncomp = sum(egroup.ncomp)
                for n in range(self.nMembranes):
                    bytes_read = f.read(4*ncomp)
                    if bytes_read=='':break
                    vals = unpack_from('f'*ncomp,bytes_read)
                    ind = 0
                    for rt in egroup.res_labels:
                        krt = egroup.res_labels.index(rt)
                        if rt=='NINT':
                            ind += 1
                        elif rt=='STRAINS':
                            for k in range(egroup.ncomp[krt]):
                                strain[k].append(vals[ind])
                                ind += 1
                        elif rt=='PLA_CODE':
                            pla_code.append(vals[ind])
                            ind += 1
                        elif rt=='STR_LEVEL':
                            str_level.append(vals[ind])
                            ind += 1
                        elif rt=='SMFORCE':
                            for k in range(egroup.ncomp[krt]):
                                smforce[k].append(vals[ind])
##                                if k==0:
##                                    if abs(vals[ind])>0.001:
##                                        print(n,vals[ind])
                                ind += 1
                        else:  # general method for any result type (could replace lines above)
                            res_ind = vol_res[0].index(rt)
                            if egroup.ncomp[krt]==1:
                                vol_res[1][res_ind].append(vals[ind])
                                ind += 1
                            else:
                                for k in range(egroup.ncomp[krt]):
                                    vol_res[1][res_ind][k].append(vals[ind])
                                    ind += 1
                if 'STRAINS' in egroup.res_labels:
                    step.mem.strain = strain
                if 'SMFORCE' in egroup.res_labels:
                    step.mem.smforce = smforce
                if 'PLA_CODE' in egroup.res_labels:
                    step.mem.pla_code = pla_code
                if 'STR_LEVEL' in egroup.res_labels:
                    step.mem.str_level = str_level
                for kk in range(len(vol_res[0])):
                    setattr(step.mem,vol_res[0][kk],vol_res[1][kk])
            else:
                offset = sum(egroup.ncomp)*self.nMembranes*4
                f.seek(offset,1)
                    

        f.close()
        self.membranesRead = True

    def read_eda(self):
        file = open(self.pathname + '/' + self.problem_name + '.eda')
        if self.verbose_level<2:
            print(self.pathname + '/' + self.problem_name + '.eda')

        for line in iter(lambda: file.readline(), ""):
            if 'MATERIAL PROPERTIES' in line:
                line = file.readline()
                line = file.readline()
                line = file.readline()
                line = file.readline()
                while not '--------------------' in line:
                    if 'Material label' in line:
                        mat = material()
                        self.materials.append(mat)
                        mat.label = (line.split(':')[1][:-1].lstrip()).rstrip()
                        mat.ps = int(file.readline().split(':')[1])
                        mat.number = int(file.readline().split(':')[1])
                        mat.model = (file.readline().split(':')[1][:-1].lstrip()).rstrip()
                        line = file.readline()
                    elif '[' in line:
                        if 'UNIT WEIGHT PARAMETERS' in line:
                            line = file.readline()
                            while not '[' in line and len(line)>4 and not line[:20]=='-'*20:
                                if ':' in line:
                                    v = line.split(':',1)
                                    mat.unit_weight[v[0]] = float(v[1])
                                line = file.readline()
                        elif 'ELASTIC PARAMETERS' in line:
                            line = file.readline()
                            while not '[' in line and len(line)>4 and not line[:20]=='-'*20:
                                if ':' in line:
                                    v = line.split(':',1)
                                    mat.elastic[(v[0].lstrip()).rstrip()] = float(v[1])
                                line = file.readline()
                        elif 'FLOW PARAMETERS' in line:
                            line = file.readline()
                            while not '[' in line and len(line)>4 and not line[:20]=='-'*20:
                                if ':' in line:
                                    v = line.split(':',1)
                                    mat.flow[(v[0].lstrip()).rstrip()] = float(v[1])
                                line = file.readline()
                        elif 'NONLINEAR PARAMETERS' in line:
                            line = file.readline()
                            while not '[' in line and len(line)>4 and not line[:20]=='-'*20:
                                if ':' in line:
                                    v = line.split(':',1)
                                    mat.nonlinear[(v[0].lstrip()).rstrip()] = float(v[1])
                                line = file.readline()
                        elif 'INITIAL STATE KO PARAMETERS' in line:
                            line = file.readline()
                            while not '[' in line and len(line)>4 and not line[:20]=='-'*20:
                                if ':' in line:
                                    v = line.split(':',1)
                                    mat.initial[(v[0].lstrip()).rstrip()] = float(v[1])
                                line = file.readline()
                        else:
                            line = file.readline()
                            while not '[' in line:
                                line = file.readline()
                    else:
                        line = file.readline()
        file.close()        

        
    def read_LTF(self):
        file = open(self.pathname + '/' + self.problem_name + '.dat')

        for line in iter(lambda: file.readline(), ""):
            if 'LOADTIME' in line:
                self.nLTF = int(line.split()[1])
                for kltf in range(self.nLTF):
                    line = file.readline()
                    n = int(line.split()[1])
                    tx = []
                    ltf = []
                    for kt in range(n):
                        line = file.readline()
                        v = line.split()
                        tx.append(float(v[0]))
                        ltf.append(float(v[1]))
                    self.LTF.append([tx,ltf])
        file.close()

    def read_eig(self):
        for step in self.steps:
            ## only 1 eigenvalue analysis possible per calculation
            if step.type==7:
                if self.jobtype=='3D':
                    dim = 3
                else:
                    dim = 2

                file = open(self.pathname + '/' + self.problem_name + '.eig')

                line = file.readline()
                v = line.split(';')
                step.nModes = len(v)-2
                line = file.readline()
                line = file.readline()
                v = line.split(';')
                step.f0 = [float(v[kk+1]) for kk in range(step.nModes)]
                for kd in range(dim):
                    line = file.readline()
                    v = line.split(';')
                    print(len(v),step.nModes,line)
                    step.effm.append([float(v[kk+1]) for kk in range(step.nModes)])
                for kd in range(dim):
                    line = file.readline()
                    v = line.split(';')
                    step.pf.append([float(v[kk+1]) for kk in range(step.nModes)])
                line = file.readline()

                DISP = []
                for line in file:
                    v = line.split(';')
                    DISP.append([float(v[kk+3]) for kk in range(step.nModes*dim)])
                ##    coords[0].append(float(v[1]))
                ##    coords[1].append(float(v[2]))
                file.close()
                for kf in range(step.nModes):
                    aMode = mode(dim)
                    for kd in range(dim):
                        aMode.nodal.disp[kd] = [DISP[kn][kf*dim+kd] for kn in range(self.nNodes)]
                    step.modes.append(aMode)


    def write_vtk(self,res,pname):

        step = res.steps[0]

        if len(res.shell.inel):
            ##############################
            ##       shell
            ##############################
                
            fo = open(pname + '_shell.vtk', 'w')
        
            # header:
            string = '# vtk DataFile Version 3.0\n'
            string += 'shells\n'
            string += 'ASCII\n\n'
            fo.write(string)
        
            # dataset:
            string = 'DATASET UNSTRUCTURED_GRID\n'
            fo.write(string)
        
            # nodes:
            string = 'POINTS ' + str(len(res.coords[0])) + ' float\n'
            fo.write(string)
            for kn in range(len(res.coords[0])):
                fo.write(str(res.coords[0][kn]) + ' ' + str(res.coords[1][kn])
                         + ' ' + str(res.coords[2][kn]) + '\n')
            fo.write('\n')
        
            # elements:
            string = 'CELLS ' + str(len(res.shell.inel)) + ' ' + str(len(res.shell.inel)*5) + '\n'
            fo.write(string)
            for inel in res.shell.inel:
                string = '4'
                for k in inel:
                    string += ' ' + str(int(k)-1)
                fo.write(string + '\n')
            fo.write('\n')
        
            # cell types:
            string = 'CELL_TYPES ' + str(len(res.shell.inel)) + '\n'
            fo.write(string)
            for inel in res.shell.inel:
                fo.write('9\n')
            fo.write('\n')
        
            # cell data:
            string = 'CELL_DATA ' + str(len(res.shell.inel)) + '\n'
            shell_comp = ['fx','fy','mx','my','tx','ty','f1','f2']
            for k,c in enumerate([step.shell.smforce[0],step.shell.smforce[1],step.shell.smoment[0],step.shell.smoment[1]]):
                string = 'SCALARS ' + shell_comp[k] + ' float 1\n'
                string += 'LOOKUP_TABLE default\n'
                fo.write(string)
                for l in c:
                    fo.write(str(l) + ' ')
                fo.write('\n')

            fo.close()



        ##############################
        ##       vol
        ##############################
            
        fo = open(pname + '_vol.vtk', 'w')
    
        # header:
        string = '# vtk DataFile Version 3.0\n'
        string += 'volumics\n'
        string += 'ASCII\n\n'
        fo.write(string)
    
        # dataset:
        string = 'DATASET UNSTRUCTURED_GRID\n'
        fo.write(string)
    
        # nodes:
        string = 'POINTS ' + str(len(res.coords[0])) + ' float\n'
        fo.write(string)
        for kn in range(len(res.coords[0])):
            fo.write(str(res.coords[0][kn]) + ' '
                     + str(res.coords[1][kn]) + ' '
                     + str(res.coords[2][kn]) + '\n')
        fo.write('\n')
    
        # elements:
        string = 'CELLS ' + str(len(res.vol.inel)) + ' ' + str(len(res.vol.inel)*9) + '\n'
        fo.write(string)
        for inel in res.vol.inel:
            string = '8'
            for k in inel:
                string += ' ' + str(int(k)-1)
            fo.write(string + '\n')
        fo.write('\n')
    
        # cell types:
        string = 'CELL_TYPES ' + str(len(res.vol.inel)) + '\n'
        fo.write(string)
        for inel in res.vol.inel:
            fo.write('12\n')
        fo.write('\n')
    
        # cell data:
        string = 'CELL_DATA ' + str(len(res.vol.inel)) + '\n'
        fo.write(string + '\n')
        shell_comp = ['sxx','syy','sxy','szz','sxz','syz']
        for k,c in enumerate(step.vol.stress):
            string = 'SCALARS ' + shell_comp[k] + ' float 1\n'
            string += 'LOOKUP_TABLE default\n'
            fo.write(string)
            for l in c:
                fo.write(str(l) + ' ')
            fo.write('\n')

        fo.close()

    def compute_princ_stresses(self,ele_type,steps=0,elist=0):
        if steps==0:
            steps = range(len(self.steps)-1)
        if elist==0:
            elist = range(len(self.nVolumics))
        if 'vol' in ele_type:
            eg = self.ele_groups[self.ele_group_labels.index('VOLUMICS')]
            eg.res_labels.append('PRINC')
            eg.ncomp.append(3)
            eg.comp_labels.append(['1','2','3'])
            for kt in steps:
                if len(steps)*len(elist)>1e6:
                    if steps.index(kt)%int(len(steps)/10)==0:
                        print('computation of principal stresses - %i %%'%(1.0*steps.index(kt)/len(steps)*100))
                step = self.steps[kt]
                sxx = step.vol.stress[0]
                syy = step.vol.stress[1]
                sxy = step.vol.stress[2]
                szz = step.vol.stress[3]
                if not '2D' in ele_type:
                    sxz = step.vol.stress[4]
                    syz = step.vol.stress[5]
                step.vol.princ = [[],[],[]]
                for kele in elist:
                    if '2D' in ele_type:
                        mat = numpy.array([[sxx[kele],sxy[kele],0],
                                           [sxy[kele],syy[kele],0],
                                           [0,0,szz[kele]]])
                    else:
                        mat = numpy.array([[sxx[kele],sxy[kele],sxz[kele]],
                                           [sxy[kele],syy[kele],syz[kele]],
                                           [sxz[kele],syz[kele],szz[kele]]])
                    ev = list(la.eigvals(mat))
                    ev.sort(reverse=True)
                    step.vol.princ[0].append(ev[0])
                    step.vol.princ[1].append(ev[1])
                    step.vol.princ[2].append(ev[2])

    def compute_invariants(self,ele_type='vol',res_type='stress',steps=0):
        if steps==0:
            steps = range(len(self.steps))
        if 'vol' in ele_type:
            for kt in steps:
                if len(steps)*self.nVolumics>1e6:
                    if steps.index(kt)%int(len(steps)/10)==0:
                        print('computation of principal stresses - %i %%'%(1.0*steps.index(kt)/len(steps)*100))
                step = self.steps[kt]
                if res_type=='stress':
                    tensor = step.vol.stress
                elif res_type=='strain':
                    tensor = step.vol.strain
                else:
                    print('Error in zsoil_results.compute_invariants(): res_type %s not valid'%(res_type))
                sxx = tensor[0]
                syy = tensor[1]
                sxy = tensor[2]
                szz = tensor[3]
                
                
                step.vol.invar = [[],[],[],[]]
                for kele in range(self.nVolumics):
                    if self.jobtype=='PLANESTRAIN' or self.jobtype=='AXISYMETRY':
                
                        sigma = numpy.array([[sxx[kele],sxy[kele],0.0],
                                             [sxy[kele],syy[kele],0.0],
                                             [0.0,0.0,szz[kele]]])
                        sigv = float(1)/3*numpy.trace(sigma)*numpy.array([[1,0,0],
                                                                          [0,1,0],
                                                                          [0,0,1]])
                        sigdev = sigma - sigv
                        I1 = numpy.trace(sigma)
                        J2 = float(1)/2*numpy.trace(numpy.dot(sigdev,sigdev))
                        p = I1/3.0
                        q = math.sqrt(3*J2)
                    step.vol.invar[0].append(I1)
                    step.vol.invar[1].append(J2)
                    step.vol.invar[2].append(p)
                    step.vol.invar[3].append(q)

    def compute_area(self,inel):
        
        x = []
        y = []
        for kk in range(len(inel)):
            x.append(self.coords[0][inel[kk]-1])
            y.append(self.coords[1][inel[kk]-1])
        if len(inel)==3:
            a = 0.5*abs((x[0]-x[2])*(y[1]-y[0])-(x[0]-x[1])*(y[2]-y[0]))
        elif len(inel)==4:
            a1 = 0.5*abs((x[0]-x[2])*(y[1]-y[0])-(x[0]-x[1])*(y[2]-y[0]))
            a2 = 0.5*abs((x[0]-x[3])*(y[2]-y[0])-(x[0]-x[2])*(y[3]-y[0]))
            a = a1+a2

        return a

    def compute_volume(self,inel):
        
        vv = numpy.array([[self.coords[kc][inel[kk]-1]-self.coords[kc][inel[0]-1] for kc in range(3)] for kk in [1,2,5]])
        a = numpy.linalg.det(vv)
        vv = numpy.array([[self.coords[kc][inel[kk]-1]-self.coords[kc][inel[0]-1] for kc in range(3)] for kk in [2,3,7]])
        a += numpy.linalg.det(vv)
        vv = numpy.array([[self.coords[kc][inel[kk]-1]-self.coords[kc][inel[0]-1] for kc in range(3)] for kk in [4,5,7]])
        a += numpy.linalg.det(vv)
        vv = numpy.array([[self.coords[kc][inel[kk]-1]-self.coords[kc][inel[0]-1] for kc in range(3)] for kk in [2,5,7]])
        a += numpy.linalg.det(vv)
        vv = numpy.array([[self.coords[kc][inel[kk]-1]-self.coords[kc][inel[2]-1] for kc in range(3)] for kk in [5,6,7]])
        a += numpy.linalg.det(vv)

        return a*0.5

    def give_beams(self,x0,y0,x1,y1):
        sys.path.append('D:/Mandats/python/utils')
        import tools

        blist = []
        for ke in range(self.nBeams):
            inel = self.beam.inel[ke]
            crd0 = (self.coords[0][inel[0]-1],self.coords[1][inel[0]-1])
            crd1 = (self.coords[0][inel[1]-1],self.coords[1][inel[1]-1])
            d0 = tools.dist(x0,y0,x1,y1,crd0[0],crd0[1])
            d1 = tools.dist(x0,y0,x1,y1,crd1[0],crd1[1])
            dx0 = (crd0[0]-x0)**2+(crd0[1]-y0)**2
            dx1 = (crd1[0]-x0)**2+(crd1[1]-y0)**2
            if abs(d0)<1.e-6 and abs(d1)<1.e-6:
                blist.append((ke,(0.5*(crd0[0]+crd1[0])-x0)**2+
                              (0.5*(crd0[1]+crd1[1])-y0)**2,
                              sorted([(inel[0],dx0),(inel[1],dx1)],key=lambda d:d[1])))
        blist.sort(key=lambda d:d[1])

        return [(val[0],[val[2][0][0],val[2][1][0]]) for val in blist]

    def compute_loc_syst(self,ele_type=['beams']):
        if 'beams' in ele_type:
            # compute local coordinate systems for all beams:
            for ke in range(self.nBeams):
                inel = self.beam.inel[ke]
                dn = self.beam.dir_nodes[ke]
                v0 = numpy.array([self.coords[kd][inel[0]-1] for kd in range(3)])
                v1 = numpy.array([self.coords[kd][inel[1]-1] for kd in range(3)])
                d0 = numpy.array([self.coords[kd][dn[0]-1] for kd in range(3)])
                d1 = numpy.array([self.coords[kd][dn[1]-1] for kd in range(3)])
                x = (v1-v0)
                x /= la.norm(x)
                y0 = d0-v0
                y0 /= la.norm(y0)
                y1 = d1-v1
                y1 /= la.norm(y1)
                z0 = numpy.cross(x,y0)
                z1 = numpy.cross(x,y1)
                self.beam.loc_syst.append([numpy.array([x,y0,z0]),
                                           numpy.array([x,y1,z1])])

    def get_LTF_at(self,kLTF,t):
        LTF = self.LTF[kLTF]
        v = LTF[1][0]
        if len(LTF[0])>1:
            for kk in range(1,len(LTF[0])):
                if t>=LTF[0][kk-1] and t<LTF[0][kk]:
                    t0 = LTF[0][kk-1]
                    t1 = LTF[0][kk]
                    v0 = LTF[1][kk-1]
                    v1 = LTF[1][kk]
                    v = v0 + (v1-v0)/(t1-t0)*(t-t0)
                elif t>=LTF[0][-1]:
                    v = LTF[1][-1]

        return v

    def get_tstr(self,t,t0=False):
            if t0:
                intpart = int(float('%1.2f'%(t)))
                tstr = str(intpart).rjust(3,'0')+'_'+('%1.0f'%(100*(t-intpart))).rjust(2,'0')
                intpart = int(float('%1.2f'%(t0)))
                tstr += '-'+str(intpart).rjust(3,'0')+'_'+('%1.0f'%(100*(t0-intpart))).rjust(2,'0')
            else:
                intpart = int(float('%1.2f'%(t)))
                tstr = str(intpart).rjust(3,'0')+'_'+('%1.0f'%(100*(t-intpart))).rjust(2,'0')

            return tstr
