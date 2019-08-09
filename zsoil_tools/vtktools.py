#########################################################
##                                                     ##
##         Library for plotting zsoil results          ##
##              developed by M. Preisig                ##
##                 mpreisig@geomod.ch                  ##
##                      2015-2018                      ##
##                                                     ##
#########################################################

import numpy as np
import math
import vtk
from matplotlib import colors


res2import = {'NODAL':['DISP_TRA','DISP_ROT','PPRESS','PRES_HEAD'],
              'BEAMS':['FORCE','MOMENT'],
              'TRUSSES':['FORCE'],
              'SHELLS':['SMFORCE','SQFORCE','SMOMENT','THICK'],
              'VOLUMICS':['STRESESS','STRESSES','STRAINS','STR_LEVEL','PRINC','FLU_VELOC','SATUR'],
              'CONTACT':['STRESESS','STRESSES','STRAINS','PLA_CODE','STR_LEVEL']}

def create_cell_data(mesh,eg,res_group,step_res,res_labels,nEle,
                     vtk_constructor,vtk_cell_type,elementsRead):
    rdict = {'FORCE':'force',
             'STRESESS':'stress',
             'STRESSES':'stress',
             'STRAINS':'strain',
             'STR_LEVEL':'str_level',
             'THICK':'thick',
             'PRINC':'princ',
             'PLA_CODE':'pla_code',
             'PPRES':'ppres',
             'PRES_HEAD':'pres_head',
             'SATUR':'satur',
             'FLU_VELOC':'flu_veloc'
             }
    cells = vtk.vtkCellArray()
    EF = vtk.vtkIntArray()
    EF.SetNumberOfComponents(1)
    EF.SetName('EF')
    LF = vtk.vtkIntArray()
    LF.SetNumberOfComponents(1)
    LF.SetName('LF')
    mat = vtk.vtkIntArray()
    mat.SetNumberOfComponents(1)
    mat.SetName('mat')
    arrays = []
    if elementsRead:
        for rt in range(len(eg.res_labels)):
            if eg.res_labels[rt] in res_labels:
                anArr = vtk.vtkFloatArray()
                anArr.SetNumberOfComponents(eg.ncomp[rt])
                for kc,c in enumerate(eg.comp_labels[rt]):
                    anArr.SetComponentName(kc,c)
                anArr.SetName(eg.res_labels[rt])
                if rdict.has_key(eg.res_labels[rt]):
                    lab = rdict[eg.res_labels[rt]]
                else:
                    lab = eg.res_labels[rt].lower()
                arrays.append([anArr,lab,eg.ncomp[rt]])
    if type(step_res)==list:
        step0_res = step_res[0]
        step1_res = step_res[1]
        for ke in range(nEle):
            ele = vtk_constructor()
            ids = ele.GetPointIds()
            inel = res_group.inel[ke]
            for kk,kn in enumerate(inel):
                ids.SetId(kk,kn-1)
            cells.InsertNextCell(ele)
            EF.InsertNextTuple1(res_group.EF[ke])
            LF.InsertNextTuple1(res_group.LF[ke])
            mat.InsertNextTuple1(res_group.mat[ke])
            for anArr in arrays:
                aRes0 = getattr(step0_res,anArr[1])
                aRes1 = getattr(step1_res,anArr[1])
                try:
                    if anArr[2]==1:
                        anArr[0].InsertNextTuple1(aRes1[ke]-aRes0[ke])
                    else:
                        anArr[0].InsertNextTuple([aRes1[kcomp][ke]-aRes0[kcomp][ke] for kcomp in range(len(aRes1))])
                except:
                    print anArr
                    print('Results for %s have not been read for the requested step.'%(eg.type))
    else:
        for ke in range(nEle):
            ele = vtk_constructor()
            ids = ele.GetPointIds()
            inel = res_group.inel[ke]
            for kk,kn in enumerate(inel):
                ids.SetId(kk,kn-1)
            cells.InsertNextCell(ele)
            EF.InsertNextTuple1(res_group.EF[ke])
            LF.InsertNextTuple1(res_group.LF[ke])
            mat.InsertNextTuple1(res_group.mat[ke])
            for anArr in arrays:
                aRes = getattr(step_res,anArr[1])
                try:
                    if anArr[2]==1 and len(aRes)==nEle:
                        anArr[0].InsertNextTuple1(aRes[ke])
                    else:
                        anArr[0].InsertNextTuple([comp[ke] for comp in aRes])
                except:
                    print anArr
                    print('Results for %s have not been read for the requested step.'%(eg.type))
    mesh.SetCells(vtk_cell_type,cells)
    cdata = mesh.GetCellData()
    cdata.AddArray(mat)
    cdata.AddArray(EF)
    cdata.AddArray(LF)
    for anArr in arrays:
        cdata.AddArray(anArr[0])

    return cdata,cells

def create_point_data(mesh,nodal,step_nodal,res_labels,nNodes):
    rdict = {'DISP_TRA':'disp',
             'DISP_ROT':'rot',
             'PPRESS':'ppres',
             'PRES_HEAD':'pres_head'}
    pdata = mesh.GetPointData()
    arrays = []
    for rt in range(len(nodal.res_labels)):
        if nodal.res_labels[rt] in res_labels:
            anArr = vtk.vtkFloatArray()
            anArr.SetNumberOfComponents(nodal.ncomp[rt])
            for kc,c in enumerate(nodal.comp_labels[rt]):
                anArr.SetComponentName(kc,c)
            anArr.SetName(nodal.res_labels[rt])
            if rdict.has_key(nodal.res_labels[rt]):
                lab = rdict[nodal.res_labels[rt]]
            else:
                lab = nodal.res_labels[rt].lower()
            arrays.append([anArr,lab,nodal.ncomp[rt]])

    for anArr in arrays:
##        print anArr
        if type(step_nodal)==list:
            aRes0 = getattr(step_nodal[0],anArr[1])
            aRes1 = getattr(step_nodal[1],anArr[1])
            if len(aRes1)>0:
                for kn in range(nNodes):
                    anArr[0].InsertNextTuple([aRes1[kcomp][kn]-aRes0[kcomp][kn] for kcomp in range(len(aRes1))])
            pdata.AddArray(anArr[0])
        else:
            aRes = getattr(step_nodal,anArr[1])
            if len(aRes)>0:
                for kn in range(nNodes):
    ##                if anArr[2]==1:
    ##                    anArr[0].InsertNextTuple(aRes[kn])
    ##                else:
                    anArr[0].InsertNextTuple([comp[kn] for comp in aRes])
            pdata.AddArray(anArr[0])

    return pdata

def write_unstructured_grid(filename,mesh,cdata,nEle,EFs,time,verbose,
                            outline=False,cut=[]):
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetDataModeToBinary()
    writer.SetFileName(filename+'.vtu')
    extract = vtk.vtkExtractCells()
    extract.SetInputData(mesh)
    eleList = vtk.vtkIdList()
    EF = cdata.GetArray('EF')
    for k in range(nEle):
        if (EFs[EF.GetValue(k)][0]<time and EFs[EF.GetValue(k)][1]>=time) or (EFs[EF.GetValue(k)][0]==0 and time==0):
            a=eleList.InsertNextId(k)
    extract.SetCellList(eleList)
    grid = extract.GetOutputPort()
    
    if outline:
        gf = vtk.vtkGeometryFilter()
        gf.SetInputConnection(grid)
        cpd = vtk.vtkCleanPolyData()
        cpd.SetInputConnection(gf.GetOutputPort())
        cpd.Update()
        af = vtk.vtkAppendFilter()
        af.AddInputData(cpd.GetOutput())
        grid = af.GetOutputPort()
    elif len(cut):
        plane = vtk.vtkPlane()
        plane.SetOrigin(cut[0])
        plane.SetNormal(cut[1])
        cutter = vtk.vtkCutter()
        cutter.SetInputConnection(grid)
        cutter.SetCutFunction(plane)
        cutter.Update()
        cpd = vtk.vtkCleanPolyData()
        cpd.SetInputConnection(cutter.GetOutputPort())
        cpd.Update()
        af = vtk.vtkAppendFilter()
        af.AddInputData(cpd.GetOutput())
        grid = af.GetOutputPort()
        writer.SetFileName(filename+'_cut.vtu')

    writer.SetInputConnection(grid)
    writer.Write()
    if not verbose:
        print '%i elements written to %s'%(eleList.GetNumberOfIds(),filename)

def get_tstr(t,t0=False):
        if t0:
            intpart = int(float('%1.2f'%(t)))
            tstr = str(intpart).rjust(3,'0')+'_'+('%1.0f'%(100*(t-intpart))).rjust(2,'0')
            intpart = int(float('%1.2f'%(t0)))
            tstr += '-'+str(intpart).rjust(3,'0')+'_'+('%1.0f'%(100*(t0-intpart))).rjust(2,'0')
        else:
            intpart = int(float('%1.2f'%(t)))
            tstr = str(intpart).rjust(3,'0')+'_'+('%1.0f'%(100*(t-intpart))).rjust(2,'0')

        return tstr

def write_vtu(res,tsteps='all',verbose=True,
              beams=False,vol=False,shells=False,trusses=False,cnt=False,
              disp=True,outline=False,cut=[],refstep=False):

    dim = len(res.coords)
    points = vtk.vtkPoints()
    if dim==3:
        for kn in range(res.nNodes):
            points.InsertNextPoint(res.coords[0][kn],
                                   res.coords[1][kn],
                                   res.coords[2][kn])
    else:
        for kn in range(res.nNodes):
            points.InsertNextPoint(res.coords[0][kn],
                                   res.coords[1][kn],
                                   0)
    mesh = vtk.vtkUnstructuredGrid()
    mesh.SetPoints(points)

    if tsteps=='all':
        tsteps = res.out_steps

    pvdFile = open(res.problem_name+'.pvd','w')
    pvdFile.write('<?xml version="1.0"?>\n')
    pvdFile.write('<VTKFile type="Collection" version="0.1"')
    pvdFile.write(' byte_order="LittleEndian"')
    pvdFile.write(' compressor="vtkZLibDataCompressor">')
    pvdFile.write('<Collection>\n')

    for kt in tsteps:
        step = res.steps[kt]
        if not verbose:
            print 'writing step %i'%(kt)
        if refstep:
            tstr = get_tstr(step.time,refstep.time)
        else:
            tstr = get_tstr(step.time)

        if res.nodalRead:
            nodal_res = res.nodal_res[0]
            res_labels = res2import['NODAL']
            if refstep:
                STEPS = [refstep.nodal,step.nodal]
            else:
                STEPS = step.nodal
            pdata = create_point_data(mesh,nodal_res,STEPS,
                                      res_labels,res.nNodes)

        if beams:
            # write vtu for beams:
            if 'BEAMS' in res.ele_group_labels:
                ri = res.ele_group_labels.index('BEAMS')
                eg = res.ele_groups[ri]
                res_labels = res2import['BEAMS']
            else:
                res_labels = []
                eg = None
            if refstep:
                STEPS = [refstep.beam,step.beam]
            else:
                STEPS = step.beam
            cdata,cells = create_cell_data(mesh,eg,res.beam,STEPS,
                                           res_labels,res.nBeams,
                                           vtk.vtkLine,vtk.VTK_LINE,
                                           res.beamsRead)

            write_unstructured_grid(res.problem_name+'_'+tstr+'_beam',
                                    mesh,cdata,res.nBeams,res.EF,
                                    step.time,verbose)
            pvdFile.write('<DataSet timestep="%1.2f" group="" part="0" file="%s"/>\n'%
                          (step.time,res.problem_name+'_'+tstr+'_beam.vtu'))

        elif vol:
            # write vtu for volumics:
            if 'VOLUMICS' in res.ele_group_labels:
                ri = res.ele_group_labels.index('VOLUMICS')
                eg = res.ele_groups[ri]
                res_labels = res2import['VOLUMICS']
            else:
                res_labels = []
                eg = None
            if refstep:
                STEPS = [refstep.vol,step.vol]
            else:
                STEPS = step.vol
            if dim==3:
                cdata,cells = create_cell_data(mesh,eg,res.vol,STEPS,
                                               res_labels,res.nVolumics,
                                               vtk.vtkHexahedron,vtk.VTK_HEXAHEDRON,
                                               res.volumicsRead)
                write_unstructured_grid(res.problem_name+'_'+tstr+'_vol',
                                        mesh,cdata,res.nVolumics,res.EF,
                                        step.time,verbose,outline,cut)
            else:
                cdata,cells = create_cell_data(mesh,eg,res.vol,step.vol,
                                               res_labels,res.nVolumics,
                                               vtk.vtkQuad,vtk.VTK_QUAD,
                                               res.volumicsRead)
                write_unstructured_grid(res.problem_name+'_'+tstr+'_vol',
                                        mesh,cdata,res.nVolumics,res.EF,
                                        step.time,verbose)
            pvdFile.write('<DataSet timestep="%1.2f" group="" part="1" file="%s"/>\n'%
                          (step.time,res.problem_name+'_'+tstr+'_vol.vtu'))

        elif shells:
            # write vtu for shells:
            if 'SHELLS' in res.ele_group_labels:
                ri = res.ele_group_labels.index('SHELLS')
                eg = res.ele_groups[ri]
                res_labels = res2import['SHELLS']
            else:
                res_labels = []
                eg = None
            if refstep:
                STEPS = [refstep.shell,step.shell]
            else:
                STEPS = step.shell
            cdata,cells = create_cell_data(mesh,eg,res.shell,STEPS,
                                           res_labels,res.nShells,
                                           vtk.vtkQuad,vtk.VTK_QUAD,
                                           res.shellsRead)

            write_unstructured_grid(res.problem_name+'_'+tstr+'_shell',
                                    mesh,cdata,res.nShells,res.EF,
                                    step.time,verbose)
            pvdFile.write('<DataSet timestep="%1.2f" group="" part="2" file="%s"/>\n'%
                          (step.time,res.problem_name+'_'+tstr+'_shell.vtu'))

        elif cnt:
            # write vtu for trusses:
            if 'CONTACT' in res.ele_group_labels:
                ri = res.ele_group_labels.index('CONTACT')
                eg = res.ele_groups[ri]
                res_labels = res2import['CONTACT']
            else:
                res_labels = []
                eg = None
            if refstep:
                STEPS = [refstep.cnt,step.cnt]
            else:
                STEPS = step.cnt
            cdata,cells = create_cell_data(mesh,eg,res.cnt,STEPS,
                                           res_labels,res.nContacts,
                                           vtk.vtkQuad,vtk.VTK_QUAD,
                                           res.contactsRead)

            write_unstructured_grid(res.problem_name+'_'+tstr+'_cnt',
                                    mesh,cdata,res.nContacts,res.EF,
                                    step.time,verbose)
            pvdFile.write('<DataSet timestep="%1.2f" group="" part="3" file="%s"/>\n'%
                          (step.time,res.problem_name+'_'+tstr+'_cnt.vtu'))

        elif trusses:
            # write vtu for trusses:
            if 'TRUSSES' in res.ele_group_labels:
                ri = res.ele_group_labels.index('TRUSSES')
                eg = res.ele_groups[ri]
                res_labels = res2import['TRUSSES']
            else:
                res_labels = []
                eg = None
            if refstep:
                STEPS = [refstep.truss,step.truss]
            else:
                STEPS = step.truss
            cdata,cells = create_cell_data(mesh,eg,res.truss,STEPS,
                                           res_labels,res.nTrusses,
                                           vtk.vtkLine,vtk.VTK_LINE,
                                           res.trussesRead)

            write_unstructured_grid(res.problem_name+'_'+tstr+'_truss',
                                    mesh,cdata,res.nTrusses,res.EF,
                                    step.time,verbose)
            pvdFile.write('<DataSet timestep="%1.2f" group="" part="4" file="%s"/>\n'%
                          (step.time,res.problem_name+'_'+tstr+'_truss.vtu'))
    pvdFile.write('</Collection>\n</VTKFile>')
    pvdFile.close()


def project_on_plane(base,origin,pt):
    x1 = np.dot(base[0],[pt[k]-origin[k] for k in range(3)])
    x2 = np.dot(base[1],[pt[k]-origin[k] for k in range(3)])
    return (x1,x2)

def get_polylines(segments0):
    # first find separate polylines:
    segments = list(segments0)
    last = segments.pop(0)
    polylines = [[last]]
    while len(segments):
        flag = False
        for seg in segments[1:]:
            if pts_ident(last[0],seg[0]) or pts_ident(last[0],seg[1]):
                polylines[-1].append(segments.pop(segments.index(seg)))
                flag = True
                break
        if not flag:
            last = segments.pop(0)
            polylines.append([last])
    return polylines
                

def pts_ident(pt0,pt1):
    if len(pt0)==2:
        if abs((pt0[0]-pt1[0])*(pt0[1]-pt1[1]))<1e-12:
            return True        
    else:
        if abs((pt0[0]-pt1[0])*(pt0[1]-pt1[1])*(pt0[2]-pt1[2]))<1e-18:
            return True
    return False        
        

def get_section_diagram(mesh,plane):

    origin = plane.GetOrigin()
    normal = plane.GetNormal()
    if not abs(normal[1]-1)<1e-6:
        base = np.array([np.cross(normal,(0,1,0)),(0,1,0)])
    else:
        base = np.array([np.cross(normal,(1,0,0)),(1,0,0)])

    cutEdges = vtk.vtkCutter()
    cutEdges.SetInputData(mesh)
    cutEdges.SetCutFunction(plane)
    cutEdges.Update()
    output = cutEdges.GetOutput()
    bounds = output.GetBounds()
    lscale = ((bounds[1]-bounds[0])**2
              +(bounds[3]-bounds[2])**2
              +(bounds[5]-bounds[4])**2)**0.5
    cd = output.GetCellData()

    mbd = vtk.vtkMultiBlockDataSet()
    ncomp = 0
    for ka in range(cd.GetNumberOfArrays()):
        anArr = cd.GetArray(ka)
        ncomp += anArr.GetNumberOfComponents()
    mbd.SetNumberOfBlocks(ncomp)

    kkc = 0
    for ka in range(cd.GetNumberOfArrays()):
        anArr = cd.GetArray(ka)
        ncomp = anArr.GetNumberOfComponents()
        for kc in range(ncomp):
            arange = anArr.GetRange(kc)

            pd = vtk.vtkPolyData()
            points = vtk.vtkPoints()
            cells = vtk.vtkCellArray()

            try:
                scale = lscale/max([abs(arange[1]),abs(arange[0])])
            except:
                scale = lscale
            for ke in range(output.GetNumberOfCells()):
                Ids = output.GetCell(ke).GetPointIds()
                v = [anArr.GetTuple(Ids.GetId(kk))[kc] for kk in range(2)]
                quad = vtk.vtkQuad()
                pt0 = np.array(output.GetPoint(Ids.GetId(0)))
                pt1 = np.array(output.GetPoint(Ids.GetId(1)))
                n = np.cross(pt1-pt0,normal)/np.linalg.norm(pt1-pt0)
                points.InsertNextPoint(pt0)
                points.InsertNextPoint(pt1)
                points.InsertNextPoint(pt1+n*anArr.GetTuple(ke)[kc]*scale)
                points.InsertNextPoint(pt0+n*anArr.GetTuple(ke)[kc]*scale)
                quad.GetPointIds().SetId(0,ke*4+0)
                quad.GetPointIds().SetId(1,ke*4+1)
                quad.GetPointIds().SetId(2,ke*4+2)
                quad.GetPointIds().SetId(3,ke*4+3)
                cells.InsertNextCell(quad)
            pd.SetPoints(points)
            pd.SetPolys(cells)
            mbd.GetMetaData(kkc).Set(vtk.vtkCompositeDataSet.NAME(),anArr.GetName())
            mbd.SetBlock(kkc,pd)
            kkc += 1
        
    return mbd

def get_section(mesh,plane,origin=0,loc_syst=[],matlist=[],EFlist=[],LFlist=[],thlist=[],disp=False):

    if origin==0:
        origin = plane.GetOrigin()
    normal = plane.GetNormal()
    if not len(loc_syst):
        if not abs(normal[1]-1)<1e-6:
            base = np.array([np.cross(normal,(0,1,0)),(0,1,0)])
        else:
            base = np.array([np.cross(normal,(1,0,0)),(1,0,0)])
    else:
        base = [loc_syst[1],loc_syst[2]]

    cutEdges = vtk.vtkCutter()
    cutEdges.SetInputData(mesh)
    cutEdges.SetCutFunction(plane)
    cutEdges.SetGenerateTriangles(0)
    cutEdges.Update()
    
    lines = cutEdges.GetOutput()
    
    cdata = lines.GetCellData()

    points = lines.GetPoints()
    mat = cdata.GetArray('mat')
    EF = cdata.GetArray('EF')
    LF = cdata.GetArray('LF')
    th = cdata.GetArray('THICK')
    M = cdata.GetArray('SMOMENT')
    N = cdata.GetArray('SMFORCE')
    T = cdata.GetArray('SQFORCE')

    if disp:
        pdata = lines.GetPointData()
        D = pdata.GetArray('DISP_TRA')

    segments = []
    for kl in range(lines.GetNumberOfCells()):
        if len(matlist)==0 or mat.GetTuple1(kl) in matlist:
            if len(EFlist)==0 or EF.GetTuple1(kl) in EFlist:
                if len(LFlist)==0 or LF.GetTuple1(kl) in LFlist:
                    if len(thlist)==0 or th.GetTuple1(kl) in thlist:
                        line = lines.GetCell(kl)
                        id0 = line.GetPointId(0)
                        id1 = line.GetPointId(1)
                        pt0 = points.GetPoint(id0)
                        pt1 = points.GetPoint(id1)
                        if disp:
                            vals = [M.GetTuple(kl),N.GetTuple(kl),T.GetTuple(kl),
                                    [D.GetTuple(id0),D.GetTuple(id1)]]
                        else:
                            vals = [M.GetTuple(kl),N.GetTuple(kl),T.GetTuple(kl)]
                        segments.append([project_on_plane(base,origin,pt0),
                                         project_on_plane(base,origin,pt1),vals])
##                if id1>id0:
##        ##            v1 = M.GetTuple(id1)[0]
##                    ax.plot([v0,v0],[pt0[1],pt1[1]],'k')
        
    return segments

def get_section_vol(mesh,plane,loc_syst,celldata=False,
                    array='DISP_TRA',component=-1,mesh0=0,
                    matlist=[],EFlist=[],LFlist=[]):

    cut = vtk.vtkCutter()
    cut.SetInputData(mesh)
    cut.SetCutFunction(plane)
    cut.Update()
    
    if celldata:
        output0 = cut.GetOutput()
        cdata = output0.GetCellData()
        mat = cdata.GetArray('mat')
        EF = cdata.GetArray('EF')
        LF = cdata.GetArray('LF')
        extract = vtk.vtkExtractCells()
        extract.SetInputData(output0)
        eleList = vtk.vtkIdList()
        for ke in range(output0.GetNumberOfCells()):
            if len(matlist)==0 or mat.GetTuple1(ke) in matlist:
                if len(EFlist)==0 or EF.GetTuple1(ke) in EFlist:
                    if len(LFlist)==0 or LF.GetTuple1(ke) in LFlist:
                        a=eleList.InsertNextId(ke)
        extract.SetCellList(eleList)
        output1 = extract.GetOutputPort()
        geom = vtk.vtkGeometryFilter()
        geom.SetInputConnection(output1)
        geom.Update()
        output = geom.GetOutput()
        
        c2p = vtk.vtkCellDataToPointData()
        c2p.SetInputData(output)
        c2p.Update()
        pdata = c2p.GetOutput().GetPointData()
    else:
        output = cut.GetOutput()
        pdata = output.GetPointData()
    anArray = pdata.GetArray(array)
    if anArray is None:
        print 'Array '+array+' is not found.'
    val = []
    inel = []
    crd = [[],[]]


    orig = plane.GetOrigin()
    
    for kp in range(output.GetNumberOfPoints()):
        p = output.GetPoint(kp)
        crd_2D = project_on_plane(loc_syst,orig,p)
        crd[0].append(crd_2D[0])
        crd[1].append(crd_2D[1])
        if component==-1:
            val.append(anArray.GetTuple(kp))
        else:
            val.append(anArray.GetTuple(kp)[component])

    for kp in range(output.GetNumberOfPolys()):
        p = output.GetCell(kp)
        if p.GetNumberOfPoints()==3:
            inel.append([p.GetPointId(kk) for kk in range(3)])

    return val,crd,output

def get_curved_section_vol(mesh,polyline,ylim,dx,direction=np.array([0,1,0]),celldata=False,
                           array='DISP_TRA',component=-1,mesh0=0,
                           matlist=[],EFlist=[],LFlist=[]):

    locator = vtk.vtkCellLocator()
    locator.SetDataSet(mesh)
    locator.BuildLocator()

    if celldata:
        print('celldata not implemented!')
    else:
        volarr = mesh.GetPointData().GetArray(array)
    linearr = vtk.vtkFloatArray()
    linearr.SetName(array)
    linearr.SetNumberOfComponents(volarr.GetNumberOfComponents())
    crd_2d = vtk.vtkFloatArray()
    crd_2d.SetName('2Dcrds')
    crd_2d.SetNumberOfComponents(2)
    
    bounds = mesh.GetBounds()
    # right now only works for vertical raycasting!
    if abs(direction[1]-1)<1e-6:
        minmax = [bounds[2],bounds[3]]
    minmax = ylim

    output = vtk.vtkPolyData()
    lines = vtk.vtkCellArray()
    points = vtk.vtkPoints()
    unique = vtk.vtkPointLocator()
    unique.InitPointInsertion(points,bounds)

    linepts = []
    for kp in range(polyline.GetNumberOfPoints()-1):
        p0 = np.array(polyline.GetPoints().GetPoint(kp))
        p1 = np.array(polyline.GetPoints().GetPoint(kp+1))
        length = np.linalg.norm(p1-p0)
        nP = int(np.ceil(length/float(dx)))
        for kkp in range(nP):
            linepts.append(p0+kkp*(p1-p0)/nP)
    linepts.append(p1)

    x2d = 0
    for kp in range(len(linepts)):
        p = linepts[kp]
        if kp>0:
            x2d += np.linalg.norm(p-linepts[kp-1])
        ids = vtk.vtkIdList()
        ind = locator.FindCellsAlongLine([p[0],minmax[1],p[2]],
                                         [p[0],minmax[0],p[2]],0.001,ids)

        for ke in range(ids.GetNumberOfIds()):
            cell = mesh.GetCell(ids.GetId(ke))
            line = vtk.vtkLine()
            n = 0
            for kf in range(cell.GetNumberOfFaces()):
                face = cell.GetFace(kf)
                t = vtk.mutable(0)
                x = [0,0,0]
                pcrd = [0,0,0]
                subId = vtk.mutable(0)
                ind = face.IntersectWithLine([p[0],minmax[1],p[2]],
                                             [p[0],minmax[0],p[2]],0.001,t,x,pcrd,subId)
                if ind:
                    ind = vtk.mutable(0)
                    isnew = unique.InsertUniquePoint(x,ind)                        
                    line.GetPointIds().SetId(n,ind)
                    if isnew:
                        cp = [0,0,0]
                        subId = vtk.mutable(0)
                        pcrd = [0,0,0]
                        dist2 = vtk.mutable(0)
                        weights = [0,0,0,0]
                        face.EvaluatePosition(x,cp,subId,pcrd,dist2,weights)
                        tups = [volarr.GetTuple(face.GetPointId(kk)) for kk in range(4)]
                        linearr.InsertNextTuple3(sum([tups[kk][0]*weights[kk] for kk in range(4)]),
                                                 sum([tups[kk][1]*weights[kk] for kk in range(4)]),
                                                 sum([tups[kk][2]*weights[kk] for kk in range(4)]))
                        crd_2d.InsertNextTuple2(x2d,x[1])
                    if n==1:
                        lines.InsertNextCell(line)
                    elif n>1:
                        print('Intersection with %i points'%(n))
                    n += 1

    output.SetPolys(lines)
    output.SetPoints(unique.GetPoints())
    output.GetPointData().AddArray(linearr)
    output.GetPointData().AddArray(crd_2d)
##    w = vtk.vtkXMLPolyDataWriter()
##    w.SetFileName('lines.vtp')
##    w.SetInputData(output)
##    w.Write()
    
    pdata = output.GetPointData()

    anArray = pdata.GetArray(array)
    if anArray is None:
        print 'Array '+array+' is not found.'
    val = []
    crd = [[],[]]

    
    for kp in range(output.GetNumberOfPoints()):
        p = output.GetPoint(kp)
        crd[0].append(crd_2d.GetTuple2(kp)[0])
        crd[1].append(crd_2d.GetTuple2(kp)[1])
        if component==-1:
            val.append(anArray.GetTuple(kp))
        else:
            val.append(anArray.GetTuple(kp)[component])

    return val,crd,output
        
from matplotlib.patches import Polygon
def get_patches(output,loc_syst,orig,array='mat'):

    points = output.GetPoints()
    anArray = output.GetCellData().GetArray(array)
    patches = []
    cvect = []
    for kc in range(output.GetNumberOfCells()):
        cell = output.GetCell(kc)
        Ids = cell.GetPointIds()
        val = anArray.GetValue(kc)
        xy = np.zeros(3*2)
        xy.shape = (3,2)
        for k in range(3):
            pt = points.GetPoint(Ids.GetId(k))
            crd_2D = project_on_plane(loc_syst,orig,pt)
            xy[k][0] = crd_2D[0]
            xy[k][1] = crd_2D[1]
        poly = Polygon(np.array(xy))
        patches.append(poly)
        cvect.append(val)
##        LIM[0][0] = min(LIM[0][0],x[kpt])
##        LIM[0][1] = min(LIM[0][1],x[kpt])
##        LIM[0][0] = max(LIM[0][0],x[kpt])
##        LIM[0][1] = max(LIM[0][1],x[kpt])

    return patches,cvect

def extract_interfaces(output,array='mat'):
    # - find exterior edges for all int domains
    # - find layer interfaces

    anArray = output.GetCellData().GetArray(array)
    points = output.GetPoints()

    values = []
    Edges = []   # edge id: point 1 x 1e6 + point 2 (points ordered)
    for kc in range(output.GetNumberOfCells()):
        val = anArray.GetValue(kc)
        if val not in values:
            values.append(val)
            Edges.append([])
        ind = values.index(val)
        edges = Edges[ind]

        Ids0 = output.GetCell(kc).GetPointIds()
        Ids = [Ids0.GetId(k) for k in range(3)]
        Ids.append(Ids[0])
        for ke in range(3):
            edge = min(Ids[ke],Ids[ke+1])*1000000 + max(Ids[ke],Ids[ke+1])
            if edge in edges:
                edges.remove(edge)
            else:
                edges.append(edge)

    # find layer interfaces:
    interfaces = []
    interf_labels = []
    for kind1 in range(len(Edges)):
        for kind2 in range(kind1+1,len(Edges)):
            interf = []
            for edge in Edges[kind1]:
                if edge in Edges[kind2]:
                    interf.append([points.GetPoint(edge/1000000),
                                   points.GetPoint(edge%1000000)])
            if len(interf)>0:
                interf_labels.append([kind1,kind2])
                interfaces.append(interf)

    return interfaces

def create_diff_dataset(mesh0,mesh1):
    if not mesh0.GetNumberOfCells()==mesh1.GetNumberOfCells():
        print('Number of elements is not equal!')
        return 0
    if not mesh0.GetNumberOfNodes()==mesh1.GetNumberOfNodes():
        print('Number of nodes is not equal!')
        return 0

class pl_view:
    def __init__(self,grid):
        self.ren = vtk.vtkRenderer()
        self.renWin = vtk.vtkRenderWindow()
        self.renWin.AddRenderer(self.ren)
        self.mapper = vtk.vtkDataSetMapper()
        self.grid = grid
        self.mapper.SetInputData(grid)
        self.vrange = []

        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        p=self.actor.GetProperty()
        p.SetEdgeVisibility(True)
        p.SetLineWidth(2)

        self.ren.AddActor(self.actor)

        self.scalar_name = ''
        self.scalarBar = 0
        self.cmap = 0
        self.colors = [[41,28,166],
                       [75,11,244],
                       [60,138,255],
                       [61,167,254],
                       [63,190,252],
                       [69,215,245],
                       [83,232,232],
                       [95,220,194],
                       [88,226,143],
                       [81,238,77],
                       [143,251,64],
                       [187,251,117],
                       [216,254,99],
                       [255,255,0],
                       [241,231,35],
                       [239,216,80],
                       [238,186,77],
                       [242,139,64],
                       [254,71,67],
                       [233,6,1],
                       [193,80,4],
                       [167,1,34]]

    def set_array(self,aname,comp=0,vrange=[]):
        pdata = self.grid.GetPointData()
        vector = pdata.GetArray(aname)
        scalar = vtk.vtkFloatArray()
        compnames = ['X','Y','Z']
        cname = aname+' '+compnames[comp]
        scalar.SetName(cname)
        scalar.SetNumberOfTuples(vector.GetNumberOfTuples())
        scalar.SetNumberOfComponents(1)
        scalar.CopyComponent(0,vector,comp)
        pdata.AddArray(scalar)
        if len(vrange)==0:
            self.vrange = scalar.GetRange()
        else:
            self.vrange = vrange
        self.mapper.SetScalarRange(self.vrange)
        pdata.SetActiveScalars(cname)
        self.scalar_name = cname

    def set_view(self,origin_offset=[0,0,0],camera_offset=[1,2,1],
                 viewsize = [1253, 817],view={}):
        bbox = self.grid.GetBounds()
        dims = [bbox[1]-bbox[0],bbox[3]-bbox[2],bbox[5]-bbox[4]]
        if len(view)==0:
            cen = self.grid.GetCenter()
            CameraFocalPoint = [cen[0]+origin_offset[0],
                                cen[1]+origin_offset[1],
                                cen[2]+origin_offset[2]]
            CameraPosition = [cen[0]+origin_offset[0]+dims[0]*camera_offset[0],
                              cen[1]+origin_offset[1]+dims[1]*camera_offset[1],
                              cen[2]+origin_offset[2]+dims[2]*camera_offset[2]]
            CameraViewUp = [0, 1, 0]
        else:
            CameraFocalPoint = view['CameraFocalPoint']
            CameraPosition = view['CameraPosition']
            CameraViewUp = view['CameraViewUp']
        Background = [0.32, 0.34, 0.43]
        Background = [0.9,0.9,0.9]

        self.renWin.SetSize(viewsize)
        camera = self.ren.GetActiveCamera()
        camera.SetPosition(CameraPosition)
        camera.SetFocalPoint(CameraFocalPoint)
        camera.SetViewUp(CameraViewUp)
        self.ren.SetBackground(Background)

        light = vtk.vtkLight()
        LightPosition = [CameraFocalPoint[0]+(CameraPosition[0]-CameraFocalPoint[0])*1,
                         CameraFocalPoint[1]+(CameraPosition[1]-CameraFocalPoint[1])*2,
                         CameraFocalPoint[2]+(CameraPosition[2]-CameraFocalPoint[2])*1]
        print CameraPosition
        print CameraFocalPoint
        print LightPosition
        light.SetPosition(LightPosition)
        light.SetFocalPoint(CameraFocalPoint)
        light.SetIntensity(1.2)
        self.ren.AddLight(light)

    def set_colorbar(self,title=''):
        self.cmap = vtk.vtkDiscretizableColorTransferFunction()
        colors = list(reversed(self.colors))
        for kc,c in enumerate(colors):
            self.cmap.AddRGBPoint(self.vrange[0]*(float(kc)/(len(colors)-1))
                             +self.vrange[1]*(1.-float(kc)/(len(colors)-1)),
                             c[0]/255.,c[1]/255.,c[2]/255.)
        self.cmap.SetNumberOfValues(10)
        self.cmap.SetDiscretize(True)
        self.cmap.SetVectorModeToComponent()
        self.cmap.Build()

        self.mapper.SetLookupTable(self.cmap)
        self.mapper.InterpolateScalarsBeforeMappingOn()

        self.scalarBar = vtk.vtkScalarBarActor()
        self.scalarBar.SetLookupTable(self.cmap)
        if title=='':
            self.scalarBar.SetTitle(self.scalar_name)
        else:
            self.scalarBar.SetTitle(title)
        tp = self.scalarBar.GetTitleTextProperty()
        tp.SetItalic(False)
        tp.SetBold(False)
        tp.SetColor(0,0,0)
        tp.SetShadow(0)
        tp.SetFontSize(8)
        tp.SetFontFamily(1)
        self.scalarBar.SetTitleTextProperty(tp)
        self.scalarBar.SetLabelTextProperty(tp)
        self.scalarBar.SetWidth(0.15)
        self.scalarBar.SetHeight(0.8)
        self.scalarBar.SetNumberOfLabels(11)
        self.scalarBar.SetAnnotationLeaderPadding(16)

        self.ren.AddActor2D(self.scalarBar)

    def write(self,name=''):
        if name=='':
            name = self.scalar_name
        renderLarge = vtk.vtkWindowToImageFilter()
        renderLarge.SetInput(self.renWin)
        renderLarge.SetMagnification(4)

        writer = vtk.vtkPNGWriter()
        writer.SetInputConnection(renderLarge.GetOutputPort())
        writer.SetFileName(name+'.png')
        self.renWin.SetOffScreenRendering(True)
        self.renWin.Render()
        writer.Write()

    def add_text(self,string,pos=[20,20],size=64):
        textActor = vtk.vtkTextActor()
        textActor.GetTextProperty().SetFontSize(size)
        textActor.SetPosition(pos[0],pos[1])
        self.ren.AddActor2D(textActor)
        textActor.SetInput(string)
        textActor.GetTextProperty().SetColor(0.0,0.0,0.0)

    def plot_to_mpl(self,fname,bounds=[],figsize=[14,14],vrange=[],label=''):
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        import matplotlib as mpl
        from matplotlib import gridspec
        import matplotlib.image as mpimg
        
        if len(vrange)==0:
            vrange = self.vrange
        if label=='':
            label = self.scalar_name
            
        fig = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(1,2,width_ratios=[1,0.15])
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])
        col = [[v/255. for v in vv] for vv in self.colors]
        cmap = colors.ListedColormap(col)
        ticks = np.linspace(vrange[0],vrange[1],11)
        norm = mpl.colors.Normalize(vmin=min(ticks), vmax=max(ticks))
        cb = mpl.colorbar.ColorbarBase(ax1,cmap=cmap,norm=norm,boundaries=ticks)
        cb.set_ticks(ticks)
        cb.ax.tick_params(labelsize=16)
        cb.set_label(label,size=20)

        img=mpimg.imread(fname+'.png')
        ax0.imshow(img)
        if len(bounds)>0:
            ax0.set_xlim(bounds[:2])
            ax0.set_ylim(bounds[2:])
            ax0.axis('off')

        return fig

def get_lut(ncol=20,lut_type='mat',maxind=20,vrange=(0,1)):

    if lut_type=='mat':
        c = ['#ffdcd2',
             '#ffa4a4',
             '#f98568',
             '#da180e',
             '#ffffc6',
             '#def538',
             '#b0b000',
             '#878e2b',
             '#dbfdc6',
             '#8bf391',
             '#5ac960',
             '#658750',
             '#e0e4fe',
             '#bb9af1',
             '#548bcf',
             '#fdcbfe',
             '#e75ae3',
             '#ad5ab4',
             '#abe3e7',
             '#67b1ae']
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfTableValues(ncol*int(np.ceil(maxind/20.)))
        for k in range(ncol*int(np.ceil(maxind/20.))):
            cv = colors.hex2color(c[k%20])
            lut.SetTableValue(k,cv[0],cv[1],cv[2])
    elif lut_type=='maps':
        c = ['#4b0bf4',
             '#3c8aff',
             '#3da7fe',
             '#3fbefc',
             '#45d7f5',
             '#53e8e8',
             '#5fdcc2',
             '#58e28f',
             '#51ee4d',
             '#8ffb40',
             '#bbfb75',
             '#d8fe63',
             '#ffff00',
             '#f1e723',
             '#efd850',
             '#eeba4d',
             '#f28b40',
             '#fe4743',
             '#e90601',
             '#c15004']
        
        lut = vtk.vtkDiscretizableColorTransferFunction()
        lut.DiscretizeOn()
        lut.SetNumberOfValues(ncol)
        dv = vrange[1]-vrange[0]
        for k,cc in enumerate(c):
            cv = colors.hex2color(cc)
            lut.AddRGBPoint(k*dv/ncol+vrange[0],cv[0],cv[1],cv[2])

    return lut



def compute_volume(pts):
    
    x = []
    y = []
    for kk in range(len(pts)):
        x.append(pts[kk][0])
        y.append(pts[kk][1])
    if len(pts)==3:
        a = 0.5*abs((x[0]-x[2])*(y[1]-y[0])-(x[0]-x[1])*(y[2]-y[0]))
    elif len(pts)==4:
        a1 = 0.5*abs((x[0]-x[2])*(y[1]-y[0])-(x[0]-x[1])*(y[2]-y[0]))
        a2 = 0.5*abs((x[0]-x[3])*(y[2]-y[0])-(x[0]-x[2])*(y[3]-y[0]))
        a = a1+a2

    return a   
        
