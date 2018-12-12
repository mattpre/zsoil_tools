import os.path,time,re

class Logfile:
    def __init__(self,f):
        self.filepath = f
        temp = str(self.filepath.split('/')[-1])
        self.name = temp[:temp.rfind('.')]
        if len(self.name)>40:
            self.name = self.name[-40:]
        self.file = open(self.filepath)
        self.norms = []
        self.tvect = []
        self.sf = []
        self.tfilepos = []
        self.currentStep = 0
        self.currentIteration = 0
        self.convStatVect = []
        self.conv_status = -2
        self.conv_dict = {-1:'Converged',0:'Diverged',-2:'N/A',-3:'Killed'}
        self.last_completed_step = 0
        self.notes = ''
        self.date = 0
        self.length = 0

    def scanTimeSteps(self):
        f = self.file
        if len(self.tfilepos):
            f.seek(self.tfilepos[-1])
        for line in iter(lambda: f.readline(), ""):
            if 'Time               :' in line:
                tcur = float(line.split(':')[1])
                self.tvect.append(tcur)
                self.sf.append(-1)
                kt = len(self.tvect)-1
                self.norms.append(0)
                self.convStatVect.append(-2)
                self.tfilepos.append(f.tell())
            elif 'Safety Factor      :' in line:
                tcur = float(line.split(':')[1])
                self.tvect.append(self.tvect[-1])
                self.sf.append(tcur)
                kt = len(self.tvect)-1
                self.norms.append(0)
                self.convStatVect.append(-2)
                self.tfilepos.append(f.tell())
            elif 'Actual Load Factor :' in line:
                tcur = float(line.split(':')[1])
                self.tvect.append(0)
                self.sf.append(-1)
                kt = len(self.tvect)-1
                self.norms.append(0)
                self.convStatVect.append(-2)
                self.tfilepos.append(f.tell())
            elif ('End of computation' in line or
                  'Divergence or lack' in line):
                self.convStatVect[-1] = 0
        self.readTstep(len(self.tvect)-1)

        # get time stamp of last update:
        date = os.path.getmtime(self.filepath)
        if date<=self.date:
            # check if file length has changed:
            if self.length<f.tell():
                self.date = time.time()
        else:
            self.date = date
        self.length = f.tell()

        self.currentStep = len(self.tvect)-1

    def readTstep(self,kt):
        f = self.file
        f.seek(self.tfilepos[kt])
        
        HSnorm = []
        HSit = []
        fnorm = []
        Enorm = []
        line = f.readline()    # '------------'
        while 'ITER' not in line:
            line = f.readline()
        col_headers = line.split()

        try:
            fpos = col_headers.index('RHS_FORC')
        except:
            fpos = col_headers.index('RHS_FFLX')
        try:
            epos = col_headers.index('DEN_DEFO')
        except:
            epos = col_headers.index('DEN_FLOW')
        it = 0
        lc = 0
##        status = -2 # still running
        while len(line)>2:
            v = line.split()
            if v[0]==str(it+1):
                fnorm.append(float(v[fpos]))
                Enorm.append(float(v[epos]))
                it += 1
            elif 'HS-s' in line:
                HSnorm.append(float(line.split(':')[1]))
                HSit.append(it)
            line = f.readline()
            if '   -----' in line:
                status = -1 # step converged
                break
            elif '--------' in line:
                status = -3 # step killed
                break
            elif 'The total execution' in line:
                status = 0 # step diverged or stopped
                break
            lc += 1
        norms = []
        norms.append(fnorm)
        norms.append(Enorm)
        norms.append(HSnorm)
        norms.append(HSit)
        self.norms[kt] = norms
        self.currentIteration = len(fnorm)
