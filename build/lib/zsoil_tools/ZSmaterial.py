# -*- coding: cp1252 -*-
##    Reading ZSoil material data from .inp-file, writing data to .xlsx-file
##    Copyright (C) 2017  Matthias Preisig
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


class GlobalData:
    def __init__(self):
        self.units = []
        self.nMat = 0
        self.unit_dict = {'M':'[-]',
                          'K0y':'[-]',
                          'K0x':'[-]',
                          'gamma_D':'[F/L3]',
                          'KSR_0':'[-]',
                          'nu_ur':'[-]',
                          'OCR':'[-]',
                          'Eref_ur':'[P]',
                          'Eref_50':'[P]',
                          'e_0':'[-]',
                          'gamma_f':'[F/L3]',
                          'INT5':'[-]',
                          'INT4':'[-]',
                          'INT7':'[-]',
                          'INT6':'[-]',
                          'INT1':'[-]',
                          'INT3':'[-]',
                          'E_oed':'[P]',
                          'Eref_0':'[P]',
                          'INT9':'[-]',
                          'nu':'[-]',
                          'phi':'[D]',
                          'psi':'[D]',
                          'E':'[P]',
                          'D':'[-]',
                          'H':'[-]',
                          'INT2':'[-]',
                          'SS_ext':'[-]',
                          'pmin_co':'[P]',
                          'sigref_oed':'[P]',
                          'R_f':'[-]',
                          'c':'[P]',
                          'sig_L':'[P]',
                          'KNC_0':'[-]',
                          'm':'[-]',
                          'gam_07':'[-]',
                          'INT8':'[-]',
                          'f_t':'[P]',
                          'f_c':'[P]',
                          'sig_ref':'[P]',
                          'gamma':'[F/L3]',
                          'Ix':'[L4]',
                          'Iy':'[L4]',
                          'Iz':'[L4]',
                          'b':'[L]',
                          'h':'[L]',
                          'inp_type':'[-]',
                          'Ax':'[L2]',
                          'Ay':'[L2]',
                          'Az':'[L2]',
                          'Kt/Kn':'[-]',
                          'Kn':'[-]',
                          'inherit':'[-]',
                          'k_x':'[L/T]',
                          'k_y':'[L/T]',
                          'k_z':'[L/T]',
                          'S_rres':'[-]',
                          'alpha':'[-]',
                          'pen':'[-]',
                          'cutoff':'[-]',
                          'K_f':'[-]',
                          'Profile':'[-]',
                          'Database':'[-]',
                          'm_b':'[-]',
                          'a':'[-]',
                          's':'[-]'}

class ZSmaterial:
    def __init__(self):
        self.version = 0
        self.number = 0
        self.name = 0
        self.type = 0
        self.buttons = []
        self.groups = []
        
        self.elas = {}
        self.main = {}
        self.nonl = {}
        self.geom = {}
        self.dens = {}
        self.flow = {}
        self.creep = {}
        self.heat = {}
        self.humid = {}
        self.inis = {}
        self.stab = {}
        self.damp = {}

    def read_Data(self,data):
        continuum = ['Elastic','HS-small strain stiffness','Mohr-Coulomb']
        
        lines = data[0]
        self.groups.append('elas')
        if 'Elastic' in self.type or 'Mohr-Coulomb' in self.type \
           or self.type=='Multilaminate' or 'Hoek-Brown' in self.type:
            self.elas['E'] = float(lines[0].split()[1])
            self.elas['nu'] = float(lines[0].split()[2])
        elif (self.type=='HS-small strain stiffness' or
              self.type=='Densification model'):
            self.elas['Eref_ur'] = float(lines[0].split()[1])
            self.elas['sig_ref'] = float(lines[0].split()[4])
            self.elas['nu_ur'] = float(lines[0].split()[2])
            self.elas['m'] = float(lines[0].split()[3])
            self.elas['sig_L'] = float(lines[0].split()[5])
            self.elas['SS_ext'] = int(float(lines[0].split()[7]))
            self.elas['Eref_0'] = float(lines[0].split()[8])
            self.elas['gam_07'] = float(lines[0].split()[9])
        elif self.type=='Contact':
            self.elas['Kn'] = float(lines[0].split()[1])
            self.elas['Kt/Kn'] = float(lines[0].split()[2])

        lines = data[1]
        if self.type=='Elastic Beam':
            self.groups.append('geom')
            self.geom['inp_type'] = int(lines[0].split()[1])
            if self.geom['inp_type']==0:  # specification by profile DB
                self.geom['Database'] = lines[1]
                self.geom['Profile'] = lines[2]
                ind = 3
            elif self.geom['inp_type']==1:  # specification by geometry
                self.geom['b'] = float(lines[1].split()[1])
                self.geom['h'] = float(lines[1].split()[3])
                ind = 2
            elif self.geom['inp_type']==2:  # specification by values
                ind = 1
            else:
                print 'Error in read_Data'
            self.geom['Ix'] = float(lines[ind].split()[3])
            self.geom['Iy'] = float(lines[ind].split()[4])
            self.geom['Iz'] = float(lines[ind].split()[5])
            self.geom['Ax'] = float(lines[ind].split()[0])
            self.geom['Ay'] = float(lines[ind].split()[1])
            self.geom['Az'] = float(lines[ind].split()[2])

        lines = data[2]
        if self.type=='Beam':
            self.groups.append('main')
            self.main['flex_based'] = int(lines[0].split()[1])

        lines = data[3]
        if 'interface' not in self.type and 'hinge' not in self.type:
            self.groups.append('dens')
            self.dens['gamma'] = float(lines[0].split()[1])
            if self.type in continuum:
                self.dens['gamma_f'] = float(lines[0].split()[2])
                self.dens['e_0'] = float(lines[0].split()[3])
                self.dens['gamma_D'] = float(lines[0].split()[6])

        lines = data[4]
        if self.buttons[4]==1 and not self.type=='Seepage' and 'hinge' not in self.type:
            self.groups.append('flow')
            self.flow['k_x'] = float(lines[0].split()[1])
            self.flow['k_y'] = float(lines[0].split()[2])
            self.flow['k_z'] = float(lines[0].split()[3])
            self.flow['K_f'] = float(lines[0].split()[11])
            if self.type in continuum:
                self.flow['S_rres'] = float(lines[0].split()[12])
                self.flow['alpha'] = float(lines[0].split()[13])
                self.flow['pen'] = float(lines[0].split()[17])
                self.flow['cutoff'] = float(lines[0].split()[18])

        lines = data[5]
        if self.buttons[5]==1:
            self.groups.append('creep')
            self.creep['Ad'] = float(lines[0].split()[4])
            self.creep['Av'] = float(lines[0].split()[6])
            self.creep['Bd'] = float(lines[0].split()[5])
            self.creep['Bv'] = float(lines[0].split()[7])
            self.creep['EF1'] = int(lines[0].split()[2])
            self.creep['EF2'] = int(lines[0].split()[3])
            self.creep['a'] = float(lines[0].split()[8])
            self.creep['b'] = float(lines[0].split()[9])

        lines = data[6]
        if 'Mohr-Coulomb' in self.type:
            self.groups.append('nonl')
            self.nonl['c'] = float(lines[0].split()[4])
            self.nonl['phi'] = float(lines[0].split()[2])
            self.nonl['psi'] = float(lines[0].split()[3])
        elif 'Hoek-Brown' in self.type:
            self.groups.append('nonl')
            self.nonl['f_c'] = float(lines[1].split()[0])
##            self.nonl['psi'] = float(lines[0].split()[3])
            self.nonl['m'] = float(lines[1].split()[1])
            self.nonl['s'] = float(lines[1].split()[2])
            self.nonl['a'] = float(lines[1].split()[3])
        elif (self.type=='HS-small strain stiffness' or
              self.type=='Densification model'):
            self.groups.append('nonl')
            self.nonl['phi'] = float(lines[0].split()[1])
            self.nonl['psi'] = float(lines[0].split()[2])
            self.nonl['c'] = float(lines[0].split()[3])
            self.nonl['Eref_50'] = float(lines[0].split()[4])
            self.nonl['R_f'] = float(lines[0].split()[5])
            self.nonl['f_t'] = int(float(lines[0].split()[6]))
            self.nonl['D'] = float(lines[0].split()[7])
            self.nonl['INT1'] = int(float(lines[0].split()[8]))
            self.nonl['INT2'] = int(float(lines[0].split()[9]))
            self.nonl['INT3'] = int(float(lines[0].split()[10]))
            self.nonl['H'] = float(lines[1].split()[0])
            self.nonl['M'] = float(lines[1].split()[1])
            self.nonl['KNC_0'] = float(lines[1].split()[2])
            self.nonl['sigref_oed'] = float(lines[1].split()[3])
            self.nonl['E_oed'] = float(lines[1].split()[4])
            self.nonl['INT4'] = int(float(lines[2].split()[0]))
            self.nonl['INT5'] = int(float(lines[2].split()[1]))
            self.nonl['INT6'] = int(float(lines[3].split()[0]))
            self.nonl['pmin_co'] = int(float(lines[3].split()[1]))
            self.nonl['INT7'] = int(float(lines[3].split()[2]))
            self.nonl['OCR'] = int(float(lines[3].split()[3]))
            self.nonl['KSR_0'] = int(float(lines[3].split()[4]))
            self.nonl['INT8'] = int(float(lines[3].split()[5]))
            self.nonl['INT9'] = int(float(lines[3].split()[6]))
        elif self.type=='Contact':
            self.groups.append('nonl')
            self.nonl['inherit'] = float(lines[0].split()[5])
            if self.nonl['inherit']==0:
                self.nonl['phi'] = float(lines[0].split()[1])
                self.nonl['psi'] = float(lines[0].split()[2])
                self.nonl['c'] = float(lines[0].split()[3])

        lines = data[7]
        if self.buttons[7]==1:
            self.groups.append('heat')
            self.heat['dil'] = float(lines[0].split()[1])
            if self.type in continuum:
                self.heat['cond'] = float(lines[0].split()[2])
                self.heat['cap'] = float(lines[0].split()[3])

        lines = data[8]
        if self.buttons[8]==1:
            self.groups.append('humid')
            self.humid['dil'] = float(lines[0].split()[1])
            self.humid['a'] = float(lines[0].split()[2])
            self.humid['WI'] = float(lines[0].split()[3])
            self.humid['DI'] = float(lines[0].split()[4])

        lines = data[9]
        if self.buttons[9]==1:
            self.groups.append('inis')
            self.inis['K0x'] = float(lines[0].split()[1])
            self.inis['K0y'] = float(lines[0].split()[2])

        lines = data[12]
        if self.buttons[13]==1:
            self.groups.append('damp')
            self.damp['alpha'] = float(lines[0].split()[1])
            self.damp['beta'] = float(lines[0].split()[2])

    def read_mat(global_data,fname):    
        f = file(fname)

        materials = []
        for line in iter(lambda: f.readline(), ""):
            if 'NUM_MATERIALS=' in line:
                global_data.nMat = int(float(line.split('=')[1]))
                line = f.readline()
                for km in range(global_data.nMat):
                    mat = ZSmaterial.ZSmaterial()
                    materials.append(mat)
                    # readmat proc:
                    mat.version = line.split()[1]
                    mat.number = line.split()[2]
                    line = f.readline()
                    mat.type = line[:-1]
                    line = f.readline()
                    mat.name = line[:-1]
                    line = f.readline()
                    v = line.split('=')[1].split()
                    for val in v:
                        mat.buttons.append(int(float(val)))
                    data = [[] for k in range(13)]
                    line = f.readline()
                    while not 'MATERIAL' in line:
                        if 'ELAS' in line:
                            ind = 0
                            data[0].append(line)
                        elif 'GEOM' in line:
                            ind = 1
                            data[1].append(line)
                        elif 'MAIN' in line:
                            ind = 2
                            data[2].append(line)
                        elif 'DENS' in line:
                            ind = 3
                            data[3].append(line)
                        elif 'FLOW' in line:
                            ind = 4
                            data[4].append(line)
                        elif 'CREEP' in line:
                            ind = 5
                            data[5].append(line)
                        elif 'NONL' in line:
                            ind = 6
                            data[6].append(line)
                        elif 'HEAT' in line:
                            ind = 7
                            data[7].append(line)
                        elif 'HUMID' in line:
                            ind = 8
                            data[8].append(line)
                        elif 'INIS' in line:
                            ind = 9
                            data[9].append(line)
                        elif 'STAB' in line:
                            ind = 10
                            data[10].append(line)
                        elif 'DISC' in line:
                            ind = 11
                            data[11].append(line)
                        elif 'DAMP' in line:
                            ind = 12
                            data[12].append(line)
                        elif not line=='\n':
                            data[ind].append(line)
                        line = f.readline()
                    print mat.number,mat.name,mat.type
                    mat.read_Data(data)
                    # end readmat proc
            elif 'STANDARD' in line:
                line = f.readline()
                v = line.split()
                for vv in v:
                    global_data.units.append(vv)
                global_data.units.append(v[0][0]+'Pa')

        f.close()

        return materials

    def write_mat(self,path,prob):
        workbook = xlsxwriter.Workbook(path + '/mat_' + prob + '.xlsx')

        # format definitions:
        subscript = workbook.add_format({'font_script':2})
        superscript = workbook.add_format({'font_script':1})
        symbol = workbook.add_format()
        symbol.set_font_name('Symbol')
        header = workbook.add_format({'bg_color':'#C5BE97','border':1,'align':'vcenter'})

        groupnames = ['Elastic','Geometry','Main','Density','Flow',
                      'Creep','Nonlinear','Heat','Humidity','Stability','Damping',
                      'Initial state']
        groupnames_active = []
        gps = ['elas','geom','main','dens','flow',
               'creep','nonl','heat','humid','stab','damp','inis']

        gps_active0 = set()
        gps_active = []
        for km,mat in enumerate(materials):
            for gp in mat.groups:
                gps_active0.add(gp)
        for kgp,gp in enumerate(gps):
            if gp in gps_active0:
                gps_active.append(gp)
                groupnames_active.append(groupnames[kgp])
        nGroups = len(gps_active)

        headers = [set() for kg in range(nGroups)]
        for km,mat in enumerate(materials):
            for kg in range(nGroups):
                data = getattr(mat,gps_active[kg])
                for key in data.iterkeys():
                    headers[kg].add(key)
        for kg,h in enumerate(headers):
            headers[kg] = list(h)

        for kg in range(nGroups):
            worksheet = workbook.add_worksheet(groupnames_active[kg])

            # write general headers:
            worksheet.write(0,0,'Version',header)
            worksheet.write(0,1,'Number',header)
            worksheet.write(0,2,'Name',header)
            worksheet.write(0,3,'Type',header)
            worksheet.set_column(0,0,6.86)
            worksheet.set_column(1,1,7.29)
            worksheet.set_column(2,2,30)
            worksheet.set_column(3,3,21.86)

            # write group-specific headers:
            for kh,head in enumerate(headers[kg]):
                if '_' in head:
                    worksheet.write_rich_string(0,4+kh,head.split('_')[0],subscript,head.split('_')[1],header)
                else:
                    worksheet.write_rich_string(0,4+kh,head,header)

            # write units:
            for kh,head in enumerate(headers[kg]):
                astr = re.sub('L',gd.units[1],gd.unit_dict[head])
                astr = re.sub('F',gd.units[0],astr)
                astr = re.sub('D',u'°',astr)
                astr = re.sub('T',gd.units[3],astr)
                astr = re.sub('H',gd.units[4],astr)
                astr = re.sub('P',gd.units[5],astr)
                worksheet.write_string(1,4+kh,astr,header)

            for km,mat in enumerate(materials):
                name = unicode(mat.name,sys.stdin.encoding)
                worksheet.write(km+2,0, 'v'+mat.version)
                worksheet.write_number(km+2,1, int(mat.number))
                worksheet.write(km+2,2, name)
                worksheet.write(km+2,3, mat.type)

                data = getattr(mat,gps_active[kg])
                for key in data:
                    worksheet.write(km+2,4+headers[kg].index(key),data[key])

        # create summary continuum:

        header = workbook.add_format({'bg_color':'#C5BE97','border':1,'align':'center'})
        header.set_align('vcenter')
        cell = workbook.add_format({'border':1,'align':'center'})
        cell.set_align('vcenter')

        worksheet = workbook.add_worksheet('Summary')

        # headers:
        worksheet.merge_range(0,0,1,0,u'Matériau',header)
        worksheet.merge_range(0,1,1,1,'Type',header)
        worksheet.write_rich_string(0,2,symbol,'g',header)
        worksheet.write_rich_string(0,3,symbol,'g',subscript,'D',header)
        worksheet.write_rich_string(0,4,'E',subscript,'50',header)
        worksheet.write_rich_string(0,5,'E',subscript,'ur',header)
        worksheet.write_rich_string(0,6,'E',subscript,'0',header)
        worksheet.write_rich_string(0,7,symbol,'s',subscript,'h,ref',header)
        worksheet.write_rich_string(0,8,symbol,'f','\'',header)
        worksheet.write(0,9,'c\'',header)
        worksheet.set_column(0,0,12)
        worksheet.set_column(1,1,22)
        worksheet.set_column(2,12,8)
        ##worksheet.set_column(3,3,5)
        ##worksheet.set_column(3,3,5)
        ##worksheet.set_column(3,3,5)
        ##worksheet.set_column(3,3,5)
        ##worksheet.set_column(3,3,5)

        # units:
        worksheet.write_rich_string(1,2,'[kN/m',superscript,'3',']',header)
        worksheet.write_rich_string(1,3,'[kN/m',superscript,'3',']',header)
        worksheet.write(1,4,'[MPa]',header)
        worksheet.write(1,5,'[MPa]',header)
        worksheet.write(1,6,'[MPa]',header)
        worksheet.write(1,7,'[kPa]',header)
        worksheet.write(1,8,u'[°]',header)
        worksheet.write(1,9,'[kPa]',header)

        # write material data:
        kkm = -1
        for km,mat in enumerate(materials):
            if (mat.type=='HS-small strain stiffness'
                or mat.type=='Densification model'
                or 'Mohr' in mat.type
                or 'Hoek-Brown' in mat.type):
                print km,mat.type,mat.name
                kkm += 1
                worksheet.write(kkm+2,0, unicode(mat.name,sys.stdin.encoding),cell)
                worksheet.write(kkm+2,1, mat.type,cell)
                worksheet.write(kkm+2,2, mat.dens['gamma'],cell)
                try:
                    worksheet.write(kkm+2,3, mat.dens['gamma_D'],cell)
                except:
                    pass
                if (mat.type=='HS-small strain stiffness' or
                    mat.type=='Densification model'):
                    worksheet.write(kkm+2,4, mat.nonl['Eref_50']*1e-3,cell)
                    worksheet.write(kkm+2,5, mat.elas['Eref_ur']*1e-3,cell)
                    worksheet.write(kkm+2,6, mat.elas['Eref_0']*1e-3,cell)
                    worksheet.write(kkm+2,7, mat.elas['sig_ref'],cell)
                else:
                    worksheet.merge_range(kkm+2,4,kkm+2,6, mat.elas['E']*1e-3,cell)
                try:
                    worksheet.write(kkm+2,8, mat.nonl['phi'],cell)
                    worksheet.write(kkm+2,9, mat.nonl['c'],cell)
                except:
                    worksheet.write(kkm+2,8, '',cell)
                    worksheet.write(kkm+2,9, '',cell)
                

        workbook.close()
    
