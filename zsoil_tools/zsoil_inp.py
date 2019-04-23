##    Reading and writing ZSoil .inp-files
##    Copyright (C) 2018  Matthias Preisig
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


import struct,re,math
import numpy as np
from numpy import linalg as la
from struct import unpack
from struct import unpack_from
import vtk

face_nds = [[0,1,2,3],
            [0,1,5,4],
            [0,4,7,3],
            [1,2,6,5],
            [2,3,7,6],
            [4,5,6,7]]

def eformat(f, prec, exp_digits):
    s = "%.*e"%(prec, f)
    mantissa, exp = s.split('e')
    # add 1 to digits as 1 is taken by sign +/-
    return "%se%+0*d"%(mantissa, exp_digits+1, int(exp))

class mesh:
    def __init__(self):
        self.nElements = 0
        self.nNodes = 0

class ele_info:
    """Element data
    """
    def __init__(self):
        self.number = []
        self.inel = []
        self.mat = []
        self.EF = []
        self.LF = []
        self.rm1 = []
        self.rm2 = []
        self.loc_syst = []
        self.center = [[],[],[]]
        self.dir = []

class dummy:
    def __init__(self):
        self.rm1 = []
        self.rm2 = []

class time_step:
    def __init__(self):
        #info from .his file:
        self.number = 0
        self.time = 0.0
        self.step_count = 0
        self.nodal = dummy()
        self.vol = dummy()
        self.shell = dummy()

class Material:
    def __init__(self):
        self.number = 0
        self.name = str()
        self.type = str()
        self.cross_section = 0

class ExistFun:
    def __init__(self):
        self.number = 0
        self.name = ''
        self.nInstances = 0
        self.data = []

class LoadFun:
    def __init__(self):
        self.number = 0
        self.name = ''
        self.nSteps = 0
        self.data = []
        self.flags = ''
        self.options = []

    def at(self,t):
        if t<self.data[0][0]:
            v = self.data[0][1]
        elif t>=self.data[-1][0]:
            v = self.data[-1][1]
        else:
            for kd in range(len(self.data)-1):
                if t>=self.data[kd][0] and t<self.data[kd+1][0]:
                    break
            t0 = self.data[kd][0]
            v0 = self.data[kd][1]
            t1 = self.data[kd+1][0]
            v1 = self.data[kd+1][1]
            v = v0+(v1-v0)/(t1-t0)*(t-t0)

        return v*self.options[1]

class SurfaceLoad:
    def __init__(self):
        self.number = 0
        self.name = ''
        self.type = ''
        self.faces = []
        self.LF = 0
        self.data = []
        self.dir = -1
        self.ref_node = -1
        self.CG = [0,0,0]
        self.resultant = [0,0,0]

    def compute_resultant(self,mesh):
        if self.type=='UNI_LOAD':
            A = 0
            for kf,face in enumerate(self.faces):
                inel = face[1]
                v = []
                for kk in range(4):
                    v.append(np.array([mesh.coords[k][inel[kk]-1] for k in range(3)]))
                a1 = 0.5*la.norm(np.cross(v[1]-v[0],v[0]-v[2]))
                a2 = 0.5*la.norm(np.cross(v[2]-v[0],v[0]-v[3]))
                a = a1+a2
                A += a

                self.resultant[0] += self.data[0]*a
                self.resultant[1] += self.data[1]*a
                self.resultant[2] += self.data[2]*a

                xm = 0.25*sum([v[k][0] for k in range(4)])
                ym = 0.25*sum([v[k][1] for k in range(4)])
                zm = 0.25*sum([v[k][2] for k in range(4)])
                self.CG[0] += xm*a
                self.CG[1] += ym*a
                self.CG[2] += zm*a
            for k in range(3):
                self.CG[k] /= A
        elif self.type=='GRAD_LOAD':
            A = 0
            for kf,face in enumerate(self.faces):
                inel = face[1]
                v = []
                dx = []
                dy = []
                dz = []
                val = []
                rn = self.ref_node-1
                for kk in range(4):
                    dx = mesh.coords[0][inel[kk]-1]-mesh.coords[0][rn]
                    dy = mesh.coords[1][inel[kk]-1]-mesh.coords[1][rn]
                    dz = mesh.coords[2][inel[kk]-1]-mesh.coords[2][rn]
                    val.append(self.data[0]+dx*self.data[1]
                               +dy*self.data[2]
                               +dz*self.data[3])
                    v.append(np.array([mesh.coords[k][inel[kk]-1] for k in range(3)]))
                a1 = 0.5*la.norm(np.cross(v[1]-v[0],v[0]-v[2]))
                a2 = 0.5*la.norm(np.cross(v[2]-v[0],v[0]-v[3]))
                a = a1+a2
                A += a

                self.resultant[self.dir] += 0.25*sum(val)*a

                xm = 0.25*sum([v[k][0] for k in range(4)])
                ym = 0.25*sum([v[k][1] for k in range(4)])
                zm = 0.25*sum([v[k][2] for k in range(4)])
                self.CG += 0.25/9*v[0]*(4*val[0]+2*val[3]+2*val[1]+val[2])*a
                self.CG += 0.25/9*v[1]*(4*val[1]+2*val[0]+2*val[2]+val[3])*a
                self.CG += 0.25/9*v[2]*(4*val[2]+2*val[1]+2*val[3]+val[0])*a
                self.CG += 0.25/9*v[3]*(4*val[3]+2*val[2]+2*val[0]+val[1])*a
            for k in range(3):
                self.CG[k] /= A

class ReinfSet:
    def __init__(self):
        self.nLayers = 0
        self.nMaterials = 0
        self.unknown = 0
        self.name = str()
        self.layers = []
        self.materials = []

class ReinfLayer:
    def __init__(self):
        self.number = 0
        self.length_type = 0
        self.length_type_dict = {0:'Full length',
                                 1:'Left side (total)',
                                 2:'Left side (relative)',
                                 3:'Right side (total)',
                                 4:'Right side (relative)',
                                 5:'Partial (total)',
                                 6:'Partial (relative)'}
        self.dist_l = 0
        self.dist_r = 0
        self.ypos_type = 0
        self.ypos_type_dict = {0:'From bot.',
                               1:'From top',
                               2:'From center'}
        self.ydist = 0
        self.zoffset_l = 0
        self.zoffset_r = 0
        self.diam = 0
        self.nBars = 0
        self.prestress = 0
        self.material = 0
        self.enabled = True

class ReinfMember:
    def __init__(self):
        self.reSet = str()
        self.reSet_num = 0
        self.nBeams = 0
        self.beams = []

class LayeredBeamComponent:
    def __init__(self):
        self.number = 0
        self.name = str()
        self.type = 0
        self.type_dict = {0:'Elastic',
                          1:'Elastic-perfectly plastic',
                          2:'User defined'}
        self.E = 0
        self.nu = 0
        self.E0_setup = 0
        self.reg_soft = 0
        self.char_len = 0
        self.coupTC_soft = 0
        self.creep_type = 0
        self.creep_A = 0
        self.creep_B = 0
        self.ft_fun = str()
        self.fc_fun = str()

class CrossSection:
    def __init__(self):
        self.def_type = 0
        self.def_type_dict = {0:'Profiles',
                              1:'User',
                              2:'Values'}
        self.name = ''
        self.type = 0
        self.type_dict = {0:'Plane rectangle',
                          1:'Rectangular tube',
                          2:'Plane circular',
                          3:'Circular tube',
                          4:'Box',
                          5:'I-bisym',
                          6:'I-asym',
                          7:'T-section'}
        self.dimensions = []
        self.values = []
        
        
class NodalMass:
    def __init__(self):
        self.node = 0
        self.mass = 0
        self.name = str()

class zsoil_inp:
    """Main data structure
    """
    def __init__(self,pathname,problem_name):
        self.pathname = pathname
        self.problem_name = problem_name
        self.analysis_type = -1

        astep = time_step()
        astep.time = 0
        self.steps = [astep]
        self.out_steps = [0]

        self.nElements = 0
        self.nNodes = 0
        self.nVolumics = 0
        self.nVolumics2D = 0
        self.nContacts = 0
        self.nBeams = 0
        self.nShells = 0
        self.nTrusses = 0
        self.num_volumics = []
        self.num_shells = []
        self.num_beams = []
        self.num_trusses = []
        self.num_contacts = []

        self.nPoints = 0
        self.nLines = 0
        self.nPlines = 0
        self.nArcs = 0
        self.nSurfaces = []
        self.ptcrds = [[],[],[]]
        self.line = []
        self.arc = []
        self.pline = []
        self.surfaces = []
        
        self.nBCs = 0
        self.nSurfaceLoads = 0
        self.nNodalLoads = 0
        self.nBeamLoads = 0
        self.nLF = 0
        self.nEF = 0
        self.nMaterials = 0
        self.nReinfSets = 0
        self.nReinfMembers = 0
        self.nLayeredBeamComponents = 0
        self.nNodalMasses = 0

        self.coords = [[],[],[]]
        self.vol = ele_info()
        self.cnt = ele_info()
        self.beam = ele_info()
        self.shell = ele_info()
        self.truss = ele_info()

        self.BCs = [[],[],[]]
        self.EFs = dict()
        self.LFs = dict()
        self.materials = dict()
        self.layered_beam_components = dict()
        self.surfLoads = []
        self.nodalLoads = []
        self.beamLoads = []
        self.nodal_masses = []
        self.reinf_sets = []
        self.reinf_members = []
        self.reinf_set_names_list = []
        self.reinf_set_per_beam = []

        self.vtkVol = 0
        self.vtkShell = 0

    def read_inp(self,sections=['ing','i0g','ibg','ilg','ics','inb',
                                'pob','gob','inl','ibf','gsl','itg',
                                'brc'],
                 debug=False):
        ## Read inp file
        # 'ing': nodes
        # 'i0g': volumics
        # 'ibg': beams
        # 'ilg': shells
        # 'ics': contacts
        # 'inb': BC's (part 1)
        # 'pob': points
        # 'gob': geometrical objects (lines, polylines, arcs, circles)
        # 'inl': nodal loads
        # 'ibf': beam loads
        # 'gsl': surface loads (macro level)
        # 'itg': trusses
        # 'brc': reinforcement sets

        file = open(self.pathname + '/' + self.problem_name + '.inp')
        
        line = file.readline()
        line = file.readline()
        self.analysis_type = int(line.split()[0])
        self.nVolumics = int(line.split()[11])
        self.nVolumics2D = int(line.split()[9])
        self.nNodes = int(line.split()[5])
        self.nBCs = int(line.split()[6])
        self.nSurfaceLoads = int(line.split()[14])
        self.nNodalLoads = int(line.split()[7])
        self.nMaterials = int(line.split()[2])
        self.nEF = int(line.split()[3]) # not correct in v16.03
        self.nLF = int(line.split()[4]) # not correct in v16.03
        line = file.readline()
        self.nContacts = int(line.split()[1])
        self.nShells = int(line.split()[14])
        self.nBeams = int(line.split()[8])
        self.nBeamLoads = int(line.split()[8])
        self.nNodalMasses = int(line.split()[23])
        line = file.readline()
        self.nPoints = int(line.split()[0])
        self.nLines = int(line.split()[2])
        self.nArcs = int(line.split()[3])
        self.nPlines = int(line.split()[1])

        for line in iter(lambda: file.readline(), ""):
            if debug:
                if line[0]=='.':
                    print line,
            if '.ing' in line and 'ing' in sections:
                if debug:
                    print 'reading ing'
                for kn in range(self.nNodes):
                    line = file.readline()
                    v = line.split()
                    self.coords[0].append(float(v[1]))
                    self.coords[1].append(float(v[2]))
                    self.coords[2].append(float(v[3]))
            elif '.i0g' in line and 'i0g' in sections:
                if debug:
                    print 'reading i0g'
                for ke in range(self.nVolumics+self.nVolumics2D):
                    line = file.readline()
                    if 'B8' in line:
                        v = line.split()
                        inel = []
                        if v[2]=='B8':
                            for kk in range(8):
                                inel.append(int(v[kk+3]))
                            center = [0,0,0]
                            for kn in inel:
                                center[0] += self.coords[0][kn-1]/8
                                center[1] += self.coords[1][kn-1]/8
                                center[2] += self.coords[2][kn-1]/8
                            pos = 10
                        self.vol.inel.append(inel)
                        self.vol.number.append(int(v[1]))
                        self.vol.mat.append(int(v[pos+4]))
                        self.vol.EF.append(int(v[pos+7]))
                        self.vol.LF.append(int(v[pos+8]))
                        self.vol.rm1.append(int(v[pos+5]))
                        self.vol.rm2.append(int(v[pos+6]))
                        self.vol.center[0].append(center[0])
                        self.vol.center[1].append(center[1])
                        self.vol.center[2].append(center[2])
                        self.num_volumics.append(int(v[0]))
                    elif 'Q4' in line:
                        v = line.split()
                        inel = []
                        if v[2]=='Q4':
                            for kk in range(4):
                                inel.append(int(v[kk+3]))
                            center = [0,0,0]
                            for kn in inel:
                                center[0] += self.coords[0][kn-1]/4
                                center[1] += self.coords[1][kn-1]/4
                                center[2] += self.coords[2][kn-1]/4
                            pos = 5
                        self.vol.inel.append(inel)
                        self.vol.number.append(int(v[1]))
                        self.vol.mat.append(int(v[pos+4]))
                        self.vol.EF.append(int(v[pos+7]))
                        self.vol.LF.append(int(v[pos+8]))
                        self.vol.rm1.append(int(v[pos+5]))
                        self.vol.rm2.append(int(v[pos+6]))
                        self.vol.center[0].append(center[0])
                        self.vol.center[1].append(center[1])
                        self.vol.center[2].append(center[2])
                        self.num_volumics.append(int(v[0]))
                    elif 'T3' in line:
                        v = line.split()
                        inel = []
                        if v[2]=='T3':
                            for kk in range(3):
                                inel.append(int(v[kk+3]))
                            center = [0,0,0]
                            for kn in inel:
                                center[0] += self.coords[0][kn-1]/3
                                center[1] += self.coords[1][kn-1]/3
                                center[2] += self.coords[2][kn-1]/3
                            pos = 10
                        self.vol.inel.append(inel)
                        self.vol.number.append(int(v[1]))
                        self.num_volumics.append(int(v[0]))
##                        self.vol.mat.append(int(v[pos+4]))
##                        self.vol.EF.append(int(v[pos+7]))
##                        self.vol.LF.append(int(v[pos+8]))
##                        self.vol.rm1.append(int(v[pos+5]))
##                        self.vol.rm2.append(int(v[pos+6]))
##                        self.vol.center[0].append(center[0])
##                        self.vol.center[1].append(center[1])
##                        self.vol.center[2].append(center[2])
            elif '.ibg' in line and 'ibg' in sections:
                if debug:
                    print 'reading ibg'
                for ke in range(self.nBeams):
                    line = file.readline()
                    if 'BEL2' in line:
                        v = line.split()
                        inel = []
                        for kk in range(2):
                            inel.append(int(v[3+kk]))
                        self.beam.inel.append(inel)
                        self.beam.number.append(int(v[0]))
                        file.readline()
                        line = file.readline()
                        self.beam.dir.append(line)
                        line = file.readline()
                        v = line.split()
                        self.beam.mat.append(int(v[0]))
                        self.beam.EF.append(int(v[3]))
                        file.readline()
                        self.num_beams.append(int(v[0]))
            elif '.itg' in line and 'itg' in sections:
                if debug:
                    print 'reading itg'
                for ke in range(self.nTrusses):
                    line = file.readline()
                    if 'TRS2' in line:
                        v = line.split()
                        inel = []
                        for kk in range(2):
                            inel.append(int(v[3+kk]))
                        self.truss.inel.append(inel)
                        self.truss.number.append(int(v[0]))
                        line = file.readline()
                        v = line.split()
                        self.truss.mat.append(int(v[0]))
                        self.truss.EF.append(int(v[1]))
                        self.truss.LF.append(int(v[2]))
                        self.truss.rm1.append(int(v[3]))
                        self.truss.rm2.append(int(v[4]))
                        file.readline()
                        self.num_trusses.append(int(v[0]))
            elif '.ilg' in line and 'ilg' in sections:
                if debug:
                    print 'reading ilg'
                for ke in range(self.nShells):
                    line = file.readline()
                    if 'SXQ4' in line:
                        v = line.split()
                        inel = []
                        if v[2]=='SXQ4':
                            for kk in range(4):
                                inel.append(int(v[kk+3]))
                            center = [0,0,0]
                            for kn in inel:
                                center[0] += self.coords[0][kn-1]/4
                                center[1] += self.coords[1][kn-1]/4
                                center[2] += self.coords[2][kn-1]/4
                            pos = 6
                        self.shell.inel.append(inel)
                        self.shell.number.append(int(v[1]))
                        self.shell.mat.append(int(v[pos+4]))
                        self.shell.EF.append(int(v[pos+7]))
                        self.shell.LF.append(int(v[pos+8]))
                        self.shell.rm1.append(int(v[pos+5]))
                        self.shell.rm2.append(int(v[pos+6]))
                        self.shell.center[0].append(center[0])
                        self.shell.center[1].append(center[1])
                        self.shell.center[2].append(center[2])
                        self.num_shells.append(int(v[0]))
                if debug:
                    print 'read %i volumics'%(ke)
            elif '.ics' in line and 'ics' in sections:
                if debug:
                    print 'reading ics'
                self.cnt.doublesided = []
                self.cnt.side_pos = []
                self.cnt.side_neg = []
                for ke in range(self.nContacts):
                    line = file.readline()
                    if 'C_Q4' in line:
                        v = line.split()
                        self.cnt.number.append(int(v[3]))
                        if v[6]=='1':
                            self.cnt.doublesided.append(False)
                            if v[7]=='1':
                                self.cnt.side_pos.append(True)
                                self.cnt.side_neg.append(False)
                            elif v[7]=='2':
                                self.cnt.side_pos.append(False)
                                self.cnt.side_neg.append(True)
                        elif v[6]=='2':
                            self.cnt.doublesided.append(True)
                            self.cnt.side_pos.append(True)
                            self.cnt.side_neg.append(True)
                            
                        line = file.readline()
                        line = file.readline()
                        v = line.split()
                        volele = [int(v[0])]
                        volface = [int(v[1])]
                        
                        line = file.readline()
                        v = line.split()
                        inel = [[]]
                        for kk in range(8):
                            inel[0].append(int(v[kk]))
                        center = [0,0,0]
                        for kn in inel[0][:4]:
                            center[0] += self.coords[0][kn-1]/4
                            center[1] += self.coords[1][kn-1]/4
                            center[2] += self.coords[2][kn-1]/4
                        line = file.readline()
                        line = file.readline()
                        v = line.split()
                        mat = [int(v[1])]
                        EF = [int(v[2])]
                        LF = [int(v[3])]
                        if self.cnt.doublesided[-1]:
                            line = file.readline()
                            v = line.split()
                            volele.append(int(v[0]))
                            volface.append(int(v[1]))
                            line = file.readline()
                            v = line.split()
                            inel.append([])
                            for kk in range(8):
                                inel[1].append(int(v[kk]))
                            line = file.readline()
                            line = file.readline()
                            v = line.split()
                            mat.append(int(v[1]))
                            EF.append(int(v[2]))
                            LF.append(int(v[3]))
                        self.cnt.inel.append(inel)
                        self.cnt.mat.append(mat)
                        self.cnt.EF.append(EF)
                        self.cnt.LF.append(LF)
                        self.cnt.center[0].append(center[0])
                        self.cnt.center[1].append(center[1])
                        self.cnt.center[2].append(center[2])
                        self.num_contacts.append(int(v[1]))
                    else:
                        print 'Error: check reading of contacts'
            elif '.inb' in line and 'inb' in sections:
                if debug:
                    print 'reading inb'
                for kn in range(self.nBCs):
                    line = file.readline()
                    v = line.split()
                    if v[3]=='1':
                        self.BCs[0].append(int(v[1]))
                    if v[8]=='1':
                        self.BCs[1].append(int(v[1]))
                    if v[13]=='1':
                        self.BCs[2].append(int(v[1]))
            elif '.pob' in line and 'pob' in sections:
                if debug:
                    print 'reading pob'
                for kp in range(self.nPoints):
                    line = file.readline()
                    v = line.split()
                    self.ptcrds[0].append(float(v[1]))
                    self.ptcrds[1].append(float(v[2]))
                    self.ptcrds[2].append(float(v[3]))
            elif '.gob' in line and 'gob' in sections:
                if debug:
                    print 'reading gob'
                line = file.readline()
                n = int(line)
                for kp in range(n):
                    line = file.readline()
                    if 'GEOMLINE' in line:
                        line = file.readline()
                        line = file.readline()
                        v = line.split()
                        self.line.append([int(v[0]),int(v[1])])
                    elif 'GEOMARC' in line:
                        line = file.readline()
                        line = file.readline()
                        v = line.split()
                        anArc = [[int(v[0]),int(v[1])]]
                        line = file.readline()
                        v = line.split()
                        anArc.append([float(v[0]),float(v[1]),float(v[2])])
                        line = file.readline()
                        v = line.split()
                        anArc.append([float(v[0]),float(v[1]),float(v[2])])
                        line = file.readline()
                        v = line.split()
                        anArc.append([int(v[0]),int(v[1])])
                        self.arc.append(anArc)
                    elif 'GEOMPLINE' in line:
                        line = file.readline()
                        line = file.readline()
                        v = line.split()
                        aPline = [[int(v[0]),int(v[1])],[[],[],[]]]
                        line = file.readline()
                        npts = int(line)
                        for kpt in range(npts):
                            line = file.readline()
                            v = line.split()
                            for kc in range(3):
                                aPline[1][kc].append(float(v[kc]))
                        self.pline.append(aPline)
            elif '.inl' in line and 'inl' in sections:
                if debug:
                    print 'reading inl'
                for kp in range(self.nNodalLoads):
                    line = file.readline()
                    v = line.split()
                    nLoad = [int(v[1]),[float(val) for val in v[4:10]],int(v[10])]
                    line = file.readline()
                    nLoad.append(line[:-1])
                    self.nodalLoads.append(nLoad)
            elif '.ibf' in line and 'ibf' in sections:
                if debug:
                    print 'reading ibf'
                for kp in range(self.nBeams):
                    line = file.readline()
                    kb = int(line.split()[0])
                    nl = int(line.split()[1])
                    for kl in range(nl):
                        line = file.readline()
                        v = line.split()
                        bLoad = [kb,[float(val) for val in v[:6]],int(v[7])]
                        line = file.readline()
                        bLoad.append(line[:-1])
                        self.beamLoads.append(bLoad)
                self.nBeamLoads = len(self.beamLoads)
            elif '.gsl' in line and 'gsl' in sections:
                if debug:
                    print 'reading gsl'
                line = file.readline()
                while not line=='\n':
                    sl = SurfaceLoad()
                    if 'UNI_LOAD' in line:
                        sl.type = 'UNI_LOAD'
                        v = line.split()
                        nFaces = int(v[5])
                        sl.LF = int(v[3])
                        line = file.readline()
                        sl.name = line[:-1]
                        line = file.readline()
                        v = line.split()
                        sl.data = [float(val) for val in v]
                        for kkf in range(nFaces):
                            v = file.readline().split()
                            kele = int(v[0])
                            kf = int(v[1])
                            if kele in self.num_volumics:
                                inel0 = self.vol.inel[self.num_volumics.index(kele)]
                                inel = [inel0[kk] for kk in face_nds[kf-1]]
                            elif kele in self.num_shells:
                                inel = self.shell.inel[self.num_shells.index(kele)]
                            sl.faces.append([kele,inel])
                    elif 'GRAD_LOAD' in line:
                        sl.type = 'GRAD_LOAD'
                        v = line.split()
                        nFaces = int(v[5])
                        sl.LF = int(v[3])
                        line = file.readline()
                        sl.ref_node = int(line)
                        line = file.readline()
                        sl.name = line[:-1]
                        line = file.readline()
                        v = line.split()
                        sl.dir = int(v[0])
                        sl.data = [float(val) for val in v[1:]]
                        faces = []
                        for kkf in range(nFaces):
                            v = file.readline().split()
                            kele = int(v[0])
                            kf = int(v[1])
                            if kele in self.num_volumics:
                                inel0 = self.vol.inel[self.num_volumics.index(kele)]
                                inel = [inel0[kk] for kk in face_nds[kf-1]]
                            elif kele in self.num_shells:
                                inel = self.shell.inel[self.num_shells.index(kele)]
                            sl.faces.append([kele,inel])
                    self.surfLoads.append(sl)
                    line = file.readline()
                self.nSurfaceLoads = len(self.surfLoads)
            elif '.brc' in line:
                v = file.readline().split()
                self.nReinfSets = int(v[0])
                self.nReinfMembers = int(v[1])
                self.reinf_set_per_beam = [0 for ke in range(self.nBeams)]

                for kRSet in range(self.nReinfSets):
                    v = file.readline().split()
                    RS = ReinfSet()
                    RS.nLayers = int(v[0])
                    RS.nMaterials = int(v[1])
                    RS.unknown = int(v[2])
                    RS.name = file.readline()[:-1]
                    self.reinf_set_names_list.append(RS.name)
                    RS.materials = [int(v) for v in file.readline().split()]
                    for kl in range(RS.nLayers):
                        layer = ReinfLayer()
                        v = file.readline().split()
                        layer.number = kl+1
                        layer.length_type = int(v[2])
                        layer.dist_l = float(v[3])
                        layer.dist_r = float(v[4])
                        layer.ypos_type = int(v[5])
                        layer.ydist = float(v[6])
                        layer.zoffset_l = float(v[7])
                        layer.zoffset_r = float(v[8])
                        layer.diam = float(v[11])
                        layer.nBars = int(v[12])
                        layer.prestress = float(v[14])
                        layer.material = int(v[16])
                        if v[15]=='0':
                            layer.enabled = False
                        RS.layers.append(layer)
                    self.reinf_sets.append(RS)
                for kRMem in range(self.nReinfMembers):
                    v = file.readline().split()
                    RM = ReinfMember()
                    RM.nBeams = int(v[1])
                    if not v[0]=='0' and not v[2]=='0' and not v[3]=='0.100000':
                        print v
                    RM.reSet = file.readline()[:-1]
                    RM.reSet_num = self.reinf_set_names_list.index(RM.reSet)
                    v = file.readline().split()
                    for kk in range(RM.nBeams):
                        ke = int(v[kk])
                        RM.beams.append(ke)
                        self.reinf_set_per_beam[ke-1] = RM.reSet_num
                    self.reinf_members.append(RM)
            elif '.inm' in line:
                for km in range(self.nNodalMasses):
                    line = file.readline()
                    v = line.split()
                    m = NodalMass()
                    m.mass = float(v[4])
                    m.node = int(v[1])
                    m.number = km+1
                    line = file.readline()
                    line = file.readline()
                    m.name = line[:-1]
                    self.nodal_masses.append(m)
            elif 'NUM_MATERIALS' in line:
                if len(self.materials)==0:  # standard materials
                    self.nMaterials = int(line.split('=')[1])
                    for km in range(self.nMaterials):
                        line = file.readline()
                        mat = Material()
                        mat.number = int(line.split()[2])
                        line = file.readline()
                        mat.type = line[:-1]
                        line = file.readline()
                        mat.name = line[:-1]
                        line = file.readline()
                        mat.buttons = [int(v) for v in (line.split('=')[1]).split()]
                        while not 'DAMP->' in line:
                            if 'GEOM->' in line and mat.buttons[2]==1:
                                sect = CrossSection()
                                sect.def_type = int(line.split()[1])
                                if sect.def_type==1:
                                    line = file.readline()
                                    sect.dimensions = [float(v) for v in line.split()]
                                    line = file.readline()
                                    sect.values = [float(v) for v in line.split()]
				elif sect.def_type==0:
                                    line = file.readline()
                                    line = file.readline()
                                    sect.name = line[:-1]
                                    line = file.readline()
                                    sect.values = [float(v) for v in line.split()]
                                mat.cross_section = sect
                            line = file.readline()
                        self.materials[mat.number] = mat
                else:
                    self.nFiberMaterials = int(line.split('=')[1])
                    for km in range(self.nFiberMaterials):
                        while not 'DAMP->' in line:
                            line = file.readline()
                        
            elif 'LAYERED_BEAM_COMPONENTS' in line:
                self.nLayeredBeamComponents = int(line.split()[1])
                for kc in range(self.nLayeredBeamComponents):
                    line = file.readline()
                    comp = LayeredBeamComponent()
                    comp.number = int(line)
                    line = file.readline()
                    comp.name = line[:-1]
                    v = file.readline().split()
                    comp.type = int(v[3])
                    comp.E = float(v[0])
                    comp.nu = float(v[1])
                    comp.E0_setup = int(v[8])
                    comp.reg_soft = int(v[6])
                    comp.char_len = float(v[7])
                    comp.coupTC_soft = int(v[9])
                    comp.creep_type = int(v[10])
                    comp.creep_A = float(v[11])
                    comp.creep_B = float(v[12])
                    comp.ft_fun = file.readline()[:-1]
                    comp.fc_fun = file.readline()[:-1]
                    self.layered_beam_components[comp.number] = comp
            elif 'EXIST_FUNC' in line:
                if debug:
                    print 'reading EXIST_FUNC'
                self.nEF = int(line.split()[1])
                exfun0 = ExistFun()
                exfun0.number = 0
                exfun0.name = 'permanent'
                exfun0.nInstances = 1
                exfun0.data = [[0,1e36]]
                self.EFs[0] = exfun0
                for kEF in range(self.nEF):
                    exfun = ExistFun()
                    line = file.readline()
                    exfun.number = int(line.split()[0])
                    exfun.name = line.split(' ',1)[1][:-1]
                    line = file.readline()
                    exfun.nInstances = int(line)
                    v = file.readline().split()
                    for kl in range(exfun.nInstances):
                        exfun.data.append([float(v[kl*2]),float(v[kl*2+1])])
                    self.EFs[exfun.number] = exfun
            elif 'LOAD_FUN' in line:
                if debug:
                    print 'reading LOAD_FUN'
                self.nLF = int(line.split()[1])
                for kLF in range(self.nLF):
                    lfun = LoadFun()
                    line = file.readline()
                    lfun.number = int(line.split()[0])
                    lfun.nSteps = int(line.split()[1])
                    lfun.name = line.split(' ',2)[1][:-1]
                    line = file.readline()
                    lfun.flags = line[:-1]
                    v = file.readline().split()
                    lfun.options = [float(v[0]),float(v[1]),int(v[2])]
                    for kl in range(lfun.nSteps):
                        v = file.readline().split()
                        lfun.data.append([float(v[0]),float(v[1])])
                    self.LFs[lfun.number] = lfun
                # 0-Load function:
                lfun = LoadFun()
                lfun.number = 0
                lfun.data = [[0,1]]
                lfun.options = [0,1,0]
                self.LFs[0] = lfun
                
        file.close()

    def write_inp(self):
        f = open(self.pathname + '/' + self.problem_name + '.inp')
        of = open(self.pathname + '/' + self.problem_name + '_new.inp','w')

        for line in iter(lambda: f.readline(), ""):
            if '.i0g' in line:
                of.write(line)
                for ke in range(self.nVolumics+self.nVolumics2D):
                    line = f.readline()
                    if 'B8' in line:
                        v = line.split()
                        ke = int(v[1])-1
                        pos = 10
                        v[pos+4] = str(self.vol.mat[ke])
                        v[pos+7] = str(self.vol.EF[ke])
                        v[pos+8] = str(self.vol.LF[ke])
                        v[pos+5] = str(self.vol.rm1[ke])
                        v[pos+6] = str(self.vol.rm2[ke])
                        lmod = str()
                        for vv in v:
                            lmod += vv + ' '
                        of.write(lmod[:-1]+'\n')
                    else:
                        of.write(line)
            elif '.ilg' in line:
                of.write(line)
                for ke in range(self.nVolumics+self.nVolumics2D):
                    line = f.readline()
                    if 'SXQ4' in line:
                        v = line.split()
                        ke = int(v[1])-1
                        pos = 6
                        v[pos+4] = str(self.shell.mat[ke])
                        v[pos+7] = str(self.shell.EF[ke])
                        v[pos+8] = str(self.shell.LF[ke])
                        v[pos+5] = str(self.shell.rm1[ke])
                        v[pos+6] = str(self.shell.rm2[ke])
                        lmod = str()
                        for vv in v:
                            lmod += vv + ' '
                        of.write(lmod[:-1]+'\n')
                    else:
                        of.write(line)
            elif '.ics' in line:
                of.write(line)
                for ke in range(self.nContacts):
                    line = f.readline()
                    if 'C_Q4' in line:
                        of.write(line)
                        for kk in range(4):
                            line = f.readline()
                            of.write(line)
                        line = f.readline()
                        v = line.split()
                        v[1] = str(self.cnt.mat[ke][0])
                        v[2] = str(self.cnt.EF[ke][0])
                        v[3] = str(self.cnt.LF[ke][0])
                        lmod = str()
                        for vv in v:
                            lmod += vv + ' '
                        of.write(lmod[:-1]+'\n')
                        if self.cnt.doublesided[ke]:
                            for kk in range(3):
                                line = f.readline()
                                of.write(line)
                            line = f.readline()
                            v = line.split()
                            v[1] = str(self.cnt.mat[ke][1])
                            v[2] = str(self.cnt.EF[ke][1])
                            v[3] = str(self.cnt.LF[ke][1])
                            lmod = str()
                            for vv in v:
                                lmod += vv + ' '
                            of.write(lmod[:-1]+'\n')
                    else:
                        of.write(line)
            else:
                of.write(line)

        f.close()
        of.close()

    def write_new_inp(self):
        self.nVolumics = 0
        self.nVolumics2D = 0
        if len(self.vol.inel):
            if len(self.vol.inel[0])==4:
                self.nVolumics2D = len(self.vol.inel)
            elif len(self.vol.inel[0])==8:
                self.nVolumics = len(self.vol.inel)
                
        self.nNodes = len(self.coords[0])
        self.nPoints = len(self.ptcrds[0])
        self.nSurfaces = len(self.surfaces)
        self.nLines = len(self.line)

        f = open(self.pathname + '/' + self.problem_name + '.inp')
        of = open(self.pathname + '/' + self.problem_name + '_new.inp','w')

        of.write(f.readline())
        line = f.readline()
        v = line.split()
        v[11] = str(self.nVolumics)
        v[9] = str(self.nVolumics2D)
        v[5] = str(self.nNodes)
        v[6] = str(self.nBCs)
        v[14] = str(self.nSurfaceLoads)
        v[2] = str(self.nMaterials)
        v[3] = str(self.nEF)
        v[4] = str(self.nLF)
        string = str()
        for vv in v:
            string += vv + ' '
        of.write(string+'\n')
        line = f.readline()
        v = line.split()
        v[1] = str(self.nContacts)
        v[14] = str(self.nShells)
        string = str()
        for vv in v:
            string += vv + ' '
        of.write(string+'\n')
        
        for line in iter(lambda: f.readline(), ""):
            if '.ing' in line:
                of.write(line)
                for kn in range(self.nNodes):
                    string = '%i'%(kn+1)
                    for kk in range(3):
                        string += ' '+eformat(self.coords[kk][kn],12,3)
                    of.write(string+' 0\n')
                while not line=='\n':
                    line = f.readline()
                of.write('\n')
            elif '.i0g' in line:
                of.write(line)
                if not len(self.vol.mat):
                    self.vol.mat = [1 for k in range(self.nVolumics)]
                if not len(self.vol.rm1):
                    self.vol.rm1 = [0 for k in range(self.nVolumics)]
                if not len(self.vol.rm2):
                    self.vol.rm2 = [0 for k in range(self.nVolumics)]
                if not len(self.vol.EF):
                    self.vol.EF = [0 for k in range(self.nVolumics)]
                if not len(self.vol.LF):
                    self.vol.LF = [0 for k in range(self.nVolumics)]
                for ke in range(self.nVolumics):
                    string = str(ke+1)+' '+str(ke+1)+' B8'
                    for k in self.vol.inel[ke]:
                        string += ' '+str(k)
                    string += ' 1 1 1'
                    string += ' %i'%(self.vol.mat[ke])
                    string += ' %i'%(self.vol.rm1[ke])
                    string += ' %i'%(self.vol.rm2[ke])
                    string += ' %i'%(self.vol.EF[ke])
                    string += ' %i 0\n'%(self.vol.LF[ke])
                    of.write(string)
                for ke in range(self.nVolumics2D):
                    string = str(ke+1)+' '+str(ke+1)+' Q4'
                    for k in self.vol.inel[ke]:
                        string += ' '+str(k)
                    try:
                        string += ' '+str(self.vol.mat[ke])
                        string += ' '+str(self.vol.rm1[ke])
                        string += ' '+str(self.vol.rm2[ke])
                        string += ' '+str(self.vol.EF[ke])
                        string += ' '+str(self.vol.LF[ke])
                    except:
                        string += ' 1 1 1 0 0'
                    string += ' 0 0 0\n'
                    of.write(string)
                while not line=='\n':
                    line = f.readline()
                of.write('\n')
            elif '.ilg' in line:
                of.write(line)
##                for ke in range(self.nVolumics+self.nVolumics2D):
##                    line = f.readline()
##                    if 'SXQ4' in line:
##                        v = line.split()
##                        ke = int(v[1])-1
##                        pos = 6
##                        v[pos+4] = str(self.shell.mat[ke])
##                        v[pos+7] = str(self.shell.EF[ke])
##                        v[pos+8] = str(self.shell.LF[ke])
##                        v[pos+5] = str(self.shell.rm1[ke])
##                        v[pos+6] = str(self.shell.rm2[ke])
##                        lmod = str()
##                        for vv in v:
##                            lmod += vv + ' '
##                        of.write(lmod[:-1]+'\n')
##                    else:
##                        of.write(line)
            elif '.ics' in line:
                of.write(line)
##                for ke in range(self.nContacts):
##                    line = f.readline()
##                    if 'C_Q4' in line:
##                        of.write(line)
##                        for kk in range(4):
##                            line = f.readline()
##                            of.write(line)
##                        line = f.readline()
##                        v = line.split()
##                        v[1] = str(self.cnt.mat[ke][0])
##                        v[2] = str(self.cnt.EF[ke][0])
##                        v[3] = str(self.cnt.LF[ke][0])
##                        lmod = str()
##                        for vv in v:
##                            lmod += vv + ' '
##                        of.write(lmod[:-1]+'\n')
##                        if self.cnt.doublesided[ke]:
##                            for kk in range(3):
##                                line = f.readline()
##                                of.write(line)
##                            line = f.readline()
##                            v = line.split()
##                            v[1] = str(self.cnt.mat[ke][1])
##                            v[2] = str(self.cnt.EF[ke][1])
##                            v[3] = str(self.cnt.LF[ke][1])
##                            lmod = str()
##                            for vv in v:
##                                lmod += vv + ' '
##                            of.write(lmod[:-1]+'\n')
##                    else:
##                        of.write(line)
            elif '.pob' in line:
                of.write(line)
                
                for kn in range(self.nPoints):
                    string = '%i'%(kn+1)
                    for kk in range(3):
                        string += ' '+eformat(self.ptcrds[kk][kn],12,3)
                    of.write(string+' 0\n')
                while not line=='\n':
                    line = f.readline()
                of.write('\n')
            elif '.gob' in line:
                of.write(line)

                of.write('%i\n'%(self.nLines))
                for ke in range(self.nLines):
                    inel = self.line[ke]
                    of.write('%i GEOMLINE 0 %i\n'%(ke+1,len(inel)))
                    of.write('No name\n')
                    string = '%i'%(inel[0]+1)
                    for kk in range(len(inel)-1):
                        string += ' %i'%(inel[kk+1]+1)
                    of.write(string+'\n')
                while not line=='\n':
                    line = f.readline()
                of.write('\n')
            elif '.gos' in line:
                of.write(line)

                of.write('%i\n'%(self.nSurfaces))
                for ke in range(self.nSurfaces):
                    inel = self.surfaces[ke][0]
                    ll = self.surfaces[ke][1]
                    of.write('%i GEOMSURF 0 %i\n'%(ke+1,len(inel)))
                    of.write('No name\n')
                    string = '%i'%(inel[0]+1)
                    for kk in range(len(inel)-1):
                        string += ' %i'%(inel[kk+1]+1)
                    of.write(string+'\n')
                    of.write('%i\n'%(len(inel)))
                    string = '%i'%(ll[0]+1)
                    for kk in range(len(ll)-1):
                        string += ' %i'%(ll[kk+1]+1)
                    of.write(string+'\n')
                while not line=='\n':
                    line = f.readline()
                of.write('\n')
            else:
                of.write(line)

        f.close()
        of.close()

    def write_vtu(self,filename,pathname='.',steps=[],transform=False):
        if len(steps)==0:
            steps = self.get_all_steps()

        for kt in range(len(steps)):
            time = steps[kt]
        
            self.vtkVol = vtk.vtkUnstructuredGrid()
            points = vtk.vtkPoints()
            for kn in range(self.nNodes):
                if transform:
                    points.InsertNextPoint(transform(self.coords[0][kn],
                                                     self.coords[1][kn],
                                                     self.coords[2][kn]))
                else:
                    points.InsertNextPoint(self.coords[0][kn],
                                           self.coords[1][kn],
                                           self.coords[2][kn])
            
            names = ['mat','EF','LF','rm1','rm2']
            
            volumics = vtk.vtkCellArray()
            vol_arrays = []
            for ka in range(5):
                a = vtk.vtkIntArray()
                a.SetName(names[ka])
                vol_arrays.append(a)

            for ke in range(self.nVolumics):
                if self.exists(self.vol.EF[ke],time):
                    anEle = vtk.vtkHexahedron()
                    for kk in range(8):
                        anEle.GetPointIds().SetId(kk,self.vol.inel[ke][kk]-1)
                    volumics.InsertNextCell(anEle)
                    vol_arrays[0].InsertNextTuple1(self.vol.mat[ke])
                    vol_arrays[1].InsertNextTuple1(self.vol.EF[ke])
                    vol_arrays[2].InsertNextTuple1(self.vol.LF[ke])
                    vol_arrays[3].InsertNextTuple1(self.vol.rm1[ke])
                    vol_arrays[4].InsertNextTuple1(self.vol.rm2[ke])
            self.vtkVol.SetPoints(points)
            self.vtkVol.SetCells(vtk.VTK_HEXAHEDRON,volumics)
            
            for ke in range(self.nVolumics2D):
                if self.exists(self.vol.EF[ke],time):
                    anEle = vtk.vtkQuad()
                    for kk in range(4):
                        anEle.GetPointIds().SetId(kk,self.vol.inel[ke][kk]-1)
                    volumics.InsertNextCell(anEle)
                    vol_arrays[0].InsertNextTuple1(self.vol.mat[ke])
                    vol_arrays[1].InsertNextTuple1(self.vol.EF[ke])
                    vol_arrays[2].InsertNextTuple1(self.vol.LF[ke])
                    vol_arrays[3].InsertNextTuple1(self.vol.rm1[ke])
                    vol_arrays[4].InsertNextTuple1(self.vol.rm2[ke])
            self.vtkVol.SetPoints(points)
            self.vtkVol.SetCells(vtk.VTK_QUAD,volumics)
            
            data = self.vtkVol.GetCellData()
            for ka in range(5):
                data.AddArray(vol_arrays[ka])
            
            writer = vtk.vtkXMLUnstructuredGridWriter();
    ##        writer.SetDataModeToAscii()
            writer.SetDataModeToBinary()
            writer.SetFileName(pathname+'/'+filename+'_T%i_%i'%(time,(time-int(time))*100)+'.vtu')
            writer.SetInputData(self.vtkVol)
            writer.Write()

    def create_ug(self):
        
        ug = vtk.vtkUnstructuredGrid()
        points = vtk.vtkPoints()
        volumics = vtk.vtkCellArray()

        for kn in range(self.nNodes):
            points.InsertNextPoint(self.coords[0][kn],
                                   self.coords[1][kn],
                                   self.coords[2][kn])

        for ke in range(self.nVolumics):
            anEle = vtk.vtkHexahedron()
            for kk in range(8):
                anEle.GetPointIds().SetId(kk,self.vol.inel[ke][kk]-1)
            volumics.InsertNextCell(anEle)
        ug.SetPoints(points)
        ug.SetCells(vtk.VTK_HEXAHEDRON,volumics)

        return ug
            

    def get_all_steps(self):
        steps = set()
        for ef in self.EFs.iteritems():
            for d in ef[1].data:
                steps.add(d[0])
                steps.add(d[1])
                
        return list(sorted(list(steps)))

    def exists(self,EFnum,time):
        data = self.EFs[EFnum].data
        for d in data:
            if (time>d[0] and time<=d[1]) or (time==0 and d[0]==0):
                return True

        return False
