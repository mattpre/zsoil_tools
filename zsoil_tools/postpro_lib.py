#########################################################
##                                                     ##
##         Library for plotting zsoil results          ##
##              developed by M. Preisig                ##
##                 mpreisig@geomod.ch                  ##
##                    2015 - 2018                      ##
##                                                     ##
#########################################################

import numpy as np
import os,math,cmath
import matplotlib.pyplot as plt
import vtk


def write_vtu(res,tsteps='all',verbose=True,
              beams=False,vol=False,shells=False,trusses=False,
              disp=True):

    if tsteps=='all':
        tsteps = res.out_steps

    for kt in tsteps:
        step = res.steps[kt]
        if not verbose:
            print 'writing step %i'%(kt)
        tstr = str(int(step.time)).rjust(3,'0')+'_'+str(int((step.time-int(step.time))*100)).rjust(2,'0')

        if vol:
            # write vtu for volumics:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.DebugOn()
            writer.SetDataModeToBinary()
            writer.SetFileName(res.problem_name+'_'+tstr+'_vol.vtu')
            pdata = res.vol.mesh.GetPointData()
            if step.nodal.disp.GetNumberOfTuples()>0:
                pdata.AddArray(step.nodal.disp)
            if step.nodal.ppres.GetNumberOfTuples()>0:
                pdata.AddArray(step.nodal.ppres)
            cdata = res.vol.mesh.GetCellData()
            if step.vol.strain.GetNumberOfTuples()>0:
                cdata.AddArray(step.vol.strain)
                cdata.AddArray(step.vol.stress)
                cdata.AddArray(step.vol.str_level)
            cdata.AddArray(res.vol.EF)
            cdata.AddArray(res.vol.LF)
            cdata.AddArray(res.vol.mat)
            extract = vtk.vtkExtractCells()
            extract.SetInputData(res.vol.mesh)
            eleList = vtk.vtkIdList()
            for k in range(res.nVolumics):
                if res.EF[res.vol.EF.GetValue(k)][0]<=step.time and res.EF[res.vol.EF.GetValue(k)][1]>step.time:
                    a=eleList.InsertNextId(k)
            if not verbose:
                print '%i volumics written'%(eleList.GetNumberOfIds())
            extract.SetCellList(eleList)
            grid = extract.GetOutputPort()
            writer.SetInputConnection(grid)
            writer.Write()

        if shells:
            # write vtu for shells:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetDataModeToBinary()
            writer.SetFileName(res.problem_name+'_'+tstr+'_shell.vtu')
            pdata = res.shell.mesh.GetPointData()
            if step.nodal.disp.GetNumberOfTuples()>0:
                pdata.AddArray(step.nodal.disp)
                pdata.AddArray(step.nodal.ppres)
            cdata = res.shell.mesh.GetCellData()
            cdata.AddArray(step.shell.smforce)
            cdata.AddArray(step.shell.smoment)
            cdata.AddArray(step.shell.sqforce)
            cdata.AddArray(res.shell.EF)
            cdata.AddArray(res.shell.LF)
            cdata.AddArray(res.shell.mat)
            cdata.AddArray(step.shell.thick)
            extract = vtk.vtkExtractCells()
            extract.SetInputData(res.shell.mesh)
            eleList = vtk.vtkIdList()
            for k in range(res.nShells):
                if res.EF[res.shell.EF.GetValue(k)][0]<=step.time and res.EF[res.shell.EF.GetValue(k)][1]>step.time:
                    a=eleList.InsertNextId(k)
            if not verbose:
                print '%i shells written'%(eleList.GetNumberOfIds())
            extract.SetCellList(eleList)
            grid = extract.GetOutputPort()
            writer.SetInputConnection(grid)
            writer.Write()

        if trusses:
            # write vtu for trusses:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetDataModeToBinary()
            writer.SetFileName(res.problem_name+'_'+tstr+'_truss.vtu')
            pdata = res.truss.mesh.GetPointData()
            if step.nodal.disp.GetNumberOfTuples()>0:
                pdata.AddArray(step.nodal.disp)
                pdata.AddArray(step.nodal.ppres)
            cdata = res.truss.mesh.GetCellData()
            cdata.AddArray(step.truss.force)
            cdata.AddArray(res.truss.EF)
            cdata.AddArray(res.truss.LF)
            cdata.AddArray(res.truss.mat)
            extract = vtk.vtkExtractCells()
            extract.SetInputData(res.truss.mesh)
            eleList = vtk.vtkIdList()
            for k in range(res.nTrusses):
                if res.EF[res.truss.EF.GetValue(k)][0]<=step.time and res.EF[res.truss.EF.GetValue(k)][1]>step.time:
                    a=eleList.InsertNextId(k)
            if not verbose:
                print '%i trusses written'%(eleList.GetNumberOfIds())
            extract.SetCellList(eleList)
            grid = extract.GetOutputPort()
            writer.SetInputConnection(grid)
            writer.Write()

        if beams:
            # write vtu for beams:
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetDataModeToBinary()
            writer.SetFileName(res.problem_name+'_'+tstr+'_beam.vtu')
            pdata = res.beam.mesh.GetPointData()
            if step.nodal.disp.GetNumberOfTuples()>0:
                pdata.AddArray(step.nodal.disp)
                pdata.AddArray(step.nodal.ppres)
            cdata = res.beam.mesh.GetCellData()
            cdata.AddArray(step.beam.force)
            cdata.AddArray(step.beam.moment)
            cdata.AddArray(res.beam.EF)
            cdata.AddArray(res.beam.LF)
            cdata.AddArray(res.beam.mat)
            extract = vtk.vtkExtractCells()
            extract.SetInputData(res.beam.mesh)
            eleList = vtk.vtkIdList()
            for k in range(res.nBeams):
                if res.EF[res.beam.EF.GetValue(k)][0]<=step.time and res.EF[res.beam.EF.GetValue(k)][1]>step.time:
                    a=eleList.InsertNextId(k)
            if not verbose:
                print '%i beams written'%(eleList.GetNumberOfIds())
            extract.SetCellList(eleList)
            grid = extract.GetOutputPort()
            writer.SetInputConnection(grid)
            writer.Write()

##def write_diff(mesh0,mesh1,name):
##    diff = vtk.vtkUnstructuredGrid()
##    diff.DeepCopy(mesh1)
####    diff.SetPoints(mesh1.GetPoints())
##    cells = diff.GetCells()
####    diff.SetCells(mesh1.GetCellType(0),cells)
##    cdata = diff.GetCellData()
##    
##    cdata0 = mesh0.GetCellData()
##    cdata1 = mesh1.GetCellData()
##    Ids1 = cdata1.GetArray('vtkOriginalCellIds')
##    tmp = cdata0.GetArray('vtkOriginalCellIds')
##    Ids0 = dict(enumerate([tmp.GetTuple1(k) for k in range(tmp.GetNumberOfTuples())]))
##    print Ids0
##    for k in range(cdata1.GetNumberOfArrays()):
####        anArray1 = cdata1.GetArray(k)
##        anArray0 = cdata0.GetArray(k)
##        print anArray0.GetName()
##        anArray = cdata.GetArray(k)
##        nc = anArray.GetNumberOfComponents()
##        count = 0
##        for kc in range(cells.GetNumberOfCells()):
##            tup = anArray.GetTuple(kc)
####            tup0 = anArray0.GetTuple(Ids0.index(Ids1.GetTuple1(kc)))
##            tup0 = anArray0.GetTuple(Ids0[Ids1.GetTuple1(kc)])
####            if Ids1.GetTuple1(kc) in ind:
####                tup0 = anArray0.GetTuple(count)
##            tup1 = tuple(map(lambda x,y:x-y,tup,tup0))
####            print tup1
##            anArray.SetTuple(kc,tup1)
####                count += 1
####        Id1 = Ids1.GetTuple1(kc)
####        Id0 = Ids0.index(Id1)
####        for k in range(cdata0.GetNumberOfArrays()):
####            anArray0 = cdata0.GetArray(k)
####            anArray0.SetName(anArray0.GetName()+'0')
####            cdata.AddArray(anArray0)
####            anArray1 = cdata1.GetArray(k)
####            anArray1.SetName(anArray1.GetName()+'1')
####            cdata.AddArray(anArray1)
##
##    writer = vtk.vtkXMLUnstructuredGridWriter()
##    writer.SetFileName(name)
##    writer.SetInputData(diff)
##    writer.Write()

def write_vtu_diff(res,step0,step1,name,verbose=True,
                   beams=False,vol=False,shells=False,trusses=False,
                   disp=True):

    if not verbose:
        print 'writing step %i'%(kt)
##    tstr = str(int(step1.time)).rjust(3,'0')+'_'+str(int((step1.time-int(step1.time))*100)).rjust(2,'0')
##    tstr += '-'+str(int(step0.time)).rjust(3,'0')+'_'+str(int((step0.time-int(step0.time))*100)).rjust(2,'0')

    if vol:
        # write vtu for volumics:
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.DebugOn()
        writer.SetDataModeToBinary()
        writer.SetFileName(res.problem_name+'_'+name+'_vol.vtu')
        pdata = res.vol.mesh.GetPointData()
        if step1.nodal.disp.GetNumberOfTuples()>0:
            anArray = vtk.vtkFloatArray()
            anArray.SetNumberOfComponents(3)
            anArray.SetName(step1.nodal.disp.GetName())
            for kt in range(step1.nodal.disp.GetNumberOfTuples()):
                tup0 = step0.nodal.disp.GetTuple(kt)
                tup1 = step1.nodal.disp.GetTuple(kt)
                anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
            pdata.AddArray(anArray)
            
            anArray = vtk.vtkFloatArray()
            anArray.SetNumberOfComponents(1)
            anArray.SetName(step1.nodal.ppres.GetName())
            for kt in range(step1.nodal.ppres.GetNumberOfTuples()):
                tup0 = step0.nodal.ppres.GetTuple(kt)
                tup1 = step1.nodal.ppres.GetTuple(kt)
                anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
            pdata.AddArray(anArray)

        for att in dir(step0.vol):
            if att not in ['__doc__', '__init__', '__module__']:
                cdata = res.vol.mesh.GetCellData()
                anArray = vtk.vtkFloatArray()
                att0 = getattr(step0.vol,att)
                att1 = getattr(step1.vol,att)
                if att1.GetNumberOfTuples()>0:
                    anArray.SetNumberOfComponents(att1.GetNumberOfComponents())
                    anArray.SetName(att1.GetName())
                    for kt in range(att1.GetNumberOfTuples()):
                        tup0 = att0.GetTuple(kt)
                        tup1 = att1.GetTuple(kt)
                        anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
                    cdata.AddArray(anArray)
        
##            cdata = res.vol.mesh.GetCellData()
##            anArray = vtk.vtkFloatArray()
##            anArray.SetNumberOfComponents(step1.vol.strain.GetNumberOfComponents())
##            anArray.SetName(step1.vol.strain.GetName())
##            for kt in range(step1.vol.strain.GetNumberOfTuples()):
##                tup0 = step0.vol.strain.GetTuple(kt)
##                tup1 = step1.vol.strain.GetTuple(kt)
##                anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
##            cdata.AddArray(anArray)
##            
##            anArray = vtk.vtkFloatArray()
##            anArray.SetNumberOfComponents(step1.vol.stress.GetNumberOfComponents())
##            anArray.SetName(step1.vol.stress.GetName())
##            for kt in range(step1.vol.stress.GetNumberOfTuples()):
##                tup0 = step0.vol.stress.GetTuple(kt)
##                tup1 = step1.vol.stress.GetTuple(kt)
##                anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
##            cdata.AddArray(anArray)
##            
##            anArray = vtk.vtkFloatArray()
##            anArray.SetNumberOfComponents(step1.vol.str_level.GetNumberOfComponents())
##            anArray.SetName(step1.vol.str_level.GetName())
##            for kt in range(step1.vol.str_level.GetNumberOfTuples()):
##                tup0 = step0.vol.str_level.GetTuple(kt)
##                tup1 = step1.vol.str_level.GetTuple(kt)
##                anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
##            cdata.AddArray(anArray)
        
        cdata.AddArray(res.vol.EF)
        cdata.AddArray(res.vol.LF)
        cdata.AddArray(res.vol.mat)
        extract = vtk.vtkExtractCells()
        extract.SetInputData(res.vol.mesh)
        eleList = vtk.vtkIdList()
        for k in range(res.nVolumics):
            if res.EF[res.vol.EF.GetValue(k)][0]<=step1.time and res.EF[res.vol.EF.GetValue(k)][1]>step1.time:
                a=eleList.InsertNextId(k)
        if not verbose:
            print '%i volumics written'%(eleList.GetNumberOfIds())
        extract.SetCellList(eleList)
        grid = extract.GetOutputPort()
        writer.SetInputConnection(grid)
        writer.Write()

    if shells:
        # write vtu for shells:
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetDataModeToBinary()
        writer.SetFileName(res.problem_name+'_'+name+'_shell.vtu')
        pdata = res.shell.mesh.GetPointData()
        if step1.nodal.disp.GetNumberOfTuples()>0:
            anArray = vtk.vtkFloatArray()
            anArray.SetNumberOfComponents(3)
            anArray.SetName(step1.nodal.disp.GetName())
            for kt in range(step1.nodal.disp.GetNumberOfTuples()):
                tup0 = step0.nodal.disp.GetTuple(kt)
                tup1 = step1.nodal.disp.GetTuple(kt)
                anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
            pdata.AddArray(anArray)
            
            anArray = vtk.vtkFloatArray()
            anArray.SetNumberOfComponents(1)
            anArray.SetName(step1.nodal.ppres.GetName())
            for kt in range(step1.nodal.ppres.GetNumberOfTuples()):
                tup0 = step0.nodal.ppres.GetTuple(kt)
                tup1 = step1.nodal.ppres.GetTuple(kt)
                anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
            pdata.AddArray(anArray)

        for att in dir(step0.shell):
            if att not in ['__doc__', '__init__', '__module__']:
                cdata = res.shell.mesh.GetCellData()
                anArray = vtk.vtkFloatArray()
                att0 = getattr(step0.shell,att)
                att1 = getattr(step1.shell,att)
                if att1.GetNumberOfTuples()>0:
                    anArray.SetNumberOfComponents(att1.GetNumberOfComponents())
                    anArray.SetName(att1.GetName())
                    for kt in range(att1.GetNumberOfTuples()):
                        tup0 = att0.GetTuple(kt)
                        tup1 = att1.GetTuple(kt)
                        anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
                    cdata.AddArray(anArray)

        cdata.AddArray(res.shell.EF)
        cdata.AddArray(res.shell.LF)
        cdata.AddArray(res.shell.mat)
        cdata.AddArray(step1.shell.thick)
        extract = vtk.vtkExtractCells()
        extract.SetInputData(res.shell.mesh)
        eleList = vtk.vtkIdList()
        for k in range(res.nShells):
            if res.EF[res.shell.EF.GetValue(k)][0]<=step1.time and res.EF[res.shell.EF.GetValue(k)][1]>step1.time:
                a=eleList.InsertNextId(k)
        if not verbose:
            print '%i shells written'%(eleList.GetNumberOfIds())
        extract.SetCellList(eleList)
        grid = extract.GetOutputPort()
        writer.SetInputConnection(grid)
        writer.Write()

    if trusses:
        # write vtu for trusses:
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetDataModeToBinary()
        writer.SetFileName(res.problem_name+'_'+name+'_truss.vtu')
        pdata = res.truss.mesh.GetPointData()
        if step1.nodal.disp.GetNumberOfTuples()>0:
            anArray = vtk.vtkFloatArray()
            anArray.SetNumberOfComponents(3)
            anArray.SetName(step1.nodal.disp.GetName())
            for kt in range(step1.nodal.disp.GetNumberOfTuples()):
                tup0 = step0.nodal.disp.GetTuple(kt)
                tup1 = step1.nodal.disp.GetTuple(kt)
                anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
            pdata.AddArray(anArray)
        cdata = res.truss.mesh.GetCellData()

        for att in dir(step0.truss):
            if att not in ['__doc__', '__init__', '__module__']:
                cdata = res.truss.mesh.GetCellData()
                anArray = vtk.vtkFloatArray()
                att0 = getattr(step0.truss,att)
                att1 = getattr(step1.truss,att)
                if att1.GetNumberOfTuples()>0:
                    anArray.SetNumberOfComponents(att1.GetNumberOfComponents())
                    anArray.SetName(att1.GetName())
                    for kt in range(att1.GetNumberOfTuples()):
                        tup0 = att0.GetTuple(kt)
                        tup1 = att1.GetTuple(kt)
                        anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
                    cdata.AddArray(anArray)

        cdata.AddArray(res.truss.EF)
        cdata.AddArray(res.truss.LF)
        cdata.AddArray(res.truss.mat)
        extract = vtk.vtkExtractCells()
        extract.SetInputData(res.truss.mesh)
        eleList = vtk.vtkIdList()
        for k in range(res.nTrusses):
            if res.EF[res.truss.EF.GetValue(k)][0]<=step1.time and res.EF[res.truss.EF.GetValue(k)][1]>step1.time:
                a=eleList.InsertNextId(k)
        if not verbose:
            print '%i trusses written'%(eleList.GetNumberOfIds())
        extract.SetCellList(eleList)
        grid = extract.GetOutputPort()
        writer.SetInputConnection(grid)
        writer.Write()

    if beams:
        # write vtu for beams:
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetDataModeToBinary()
        writer.SetFileName(res.problem_name+'_'+name+'_beam.vtu')
        pdata = res.beam.mesh.GetPointData()
        if step1.nodal.disp.GetNumberOfTuples()>0:
            anArray = vtk.vtkFloatArray()
            anArray.SetNumberOfComponents(3)
            anArray.SetName(step1.nodal.disp.GetName())
            for kt in range(step1.nodal.disp.GetNumberOfTuples()):
                tup0 = step0.nodal.disp.GetTuple(kt)
                tup1 = step1.nodal.disp.GetTuple(kt)
                anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
            pdata.AddArray(anArray)
        cdata = res.beam.mesh.GetCellData()

        for att in dir(step0.beam):
            if att not in ['__doc__', '__init__', '__module__','dm']:
                cdata = res.beam.mesh.GetCellData()
                anArray = vtk.vtkFloatArray()
                att0 = getattr(step0.beam,att)
                att1 = getattr(step1.beam,att)
                if att1.GetNumberOfTuples()>0:
                    anArray.SetNumberOfComponents(att1.GetNumberOfComponents())
                    anArray.SetName(att1.GetName())
                    for kt in range(att1.GetNumberOfTuples()):
                        tup0 = att0.GetTuple(kt)
                        tup1 = att1.GetTuple(kt)
                        anArray.InsertNextTuple(tuple(map(lambda x,y:x-y,tup1,tup0)))
                    cdata.AddArray(anArray)

        cdata.AddArray(res.beam.EF)
        cdata.AddArray(res.beam.LF)
        cdata.AddArray(res.beam.mat)
        extract = vtk.vtkExtractCells()
        extract.SetInputData(res.beam.mesh)
        eleList = vtk.vtkIdList()
        for k in range(res.nBeams):
            if res.EF[res.beam.EF.GetValue(k)][0]<=step1.time and res.EF[res.beam.EF.GetValue(k)][1]>step1.time:
                a=eleList.InsertNextId(k)
        if not verbose:
            print '%i beams written'%(eleList.GetNumberOfIds())
        extract.SetCellList(eleList)
        grid = extract.GetOutputPort()
        writer.SetInputConnection(grid)
        writer.Write()
        

def project_on_plane(base,origin,pt):
    x1 = np.dot(base[0],[pt[k]-origin[k] for k in range(3)])
    x2 = np.dot(base[1],[pt[k]-origin[k] for k in range(3)])
    return (x1,x2)

def project_on_polyplane(polyplane,p):
    polyline = polyplane.GetPolyLine()

    xpcoord = 0
    for kl in range(polyline.GetNumberOfPoints()-1):
        line = vtk.vtkLine()
        p0 = polyline.GetPoints().GetPoint(kl)
        p1 = polyline.GetPoints().GetPoint(kl+1)
        length = ((p0[0]-p1[0])**2+(p0[1]-p1[1])**2+(p0[2]-p1[2])**2)**0.5
        t = vtk.mutable(0)
        tt = [0,0,0]
        d = line.DistanceToLine((p[0],-p[2],0),p0,p1,t,tt)
##        if t>-1e-6 and t<1+1e-6 and d<1e-6:
        if d<1e-6:
            xpcoord += t*length
##            print(p,p0,p1,t)
            break
        xpcoord += length

    return (xpcoord,p[1])

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
        

def get_section(mesh,plane,origin=0,loc_syst=0,matlist=[],EFlist=[],LFlist=[]):

    if origin==0:
        origin = plane.GetOrigin()
    if loc_syst==0:
        normal = plane.GetNormal()
        if not abs(normal[1]-1)<1e-6:
            base = np.array([np.cross(normal,(0,1,0)),(0,1,0)])
        else:
            base = np.array([np.cross(normal,(1,0,0)),(1,0,0)])

    cutEdges = vtk.vtkCutter()
    cutEdges.SetInputConnection(mesh)
    cutEdges.SetCutFunction(plane)
    cutEdges.SetGenerateTriangles(0)
    cutStrips = vtk.vtkStripper()
    cutStrips.PassCellDataAsFieldDataOn()
    cutStrips.PassThroughCellIdsOn()
    cutStrips.SetInputConnection(cutEdges.GetOutputPort())
    cutStrips.Update()
    cutStrips1 = vtk.vtkStripper()
    cutStrips1.PassCellDataAsFieldDataOn()
    cutStrips1.PassThroughCellIdsOn()
    cutStrips1.SetInputConnection(cutStrips.GetOutputPort())
    cutStrips1.Update()
    cutPoly = vtk.vtkPolyData()
    cutPoly.SetPoints(cutStrips.GetOutput().GetPoints())
    cutPoly.SetPolys(cutStrips.GetOutput().GetLines())
    cdata = cutPoly.GetCellData()

    # Extract boundary from cutPoly
    cutBoundary = vtk.vtkFeatureEdges()
    cutBoundary.SetInputData(cutPoly)
    cutBoundary.Update()
        

    probe = vtk.vtkProbeFilter()
    probe.SetInputConnection(cutBoundary.GetOutputPort())
    probe.SetSourceConnection(mesh)
    probe.Update()

    aPd = vtk.vtkPolyData()
    probePoints = vtk.vtkPoints()
    lines = probe.GetOutput()
    for kl in range(lines.GetNumberOfCells()):
        line = lines.GetCell(kl)
        pts = line.GetPoints()
        pt0 = pts.GetPoint(0)
        pt1 = pts.GetPoint(1)
        probePoints.InsertNextPoint(0.5*(pt0[0]+pt1[0]),
                                    0.5*(pt0[1]+pt1[1]),
                                    0.5*(pt0[2]+pt1[2]))
    aPd.SetPoints(probePoints)

    probePts = vtk.vtkProbeFilter()
    probePts.SetInputData(aPd)
    probePts.SetSourceConnection(mesh)
    probePts.Update()
    
##            # Write cut line to .vtk file
##            polyWriter = vtk.vtkPolyDataWriter()
##            polyWriter.SetInputConnection(probe.GetOutputPort())
##        ##    polyWriter.SetInputConnection(cutBoundary.GetOutputPort())
##            polyWriter.SetFileName("cutPoly.vtk")
##            polyWriter.Write()

    points = lines.GetPoints()
    cdata = probePts.GetOutput().GetPointData()
    mat = cdata.GetArray('mat')
    EF = cdata.GetArray('EF')
    LF = cdata.GetArray('LF')
    M = cdata.GetArray('SMOMENT')
    N = cdata.GetArray('SMFORCE')
    T = cdata.GetArray('SQFORCE')

    segments = []
    for kl in range(lines.GetNumberOfCells()):
        if len(matlist)==0 or mat.GetTuple1(kl) in matlist:
            if len(EFlist)==0 or EF.GetTuple1(kl) in EFlist:
                if len(LFlist)==0 or LF.GetTuple1(kl) in LFlist:
                    line = lines.GetCell(kl)
                    id0 = line.GetPointId(0)
                    id1 = line.GetPointId(1)
                    pt0 = points.GetPoint(id0)
                    pt1 = points.GetPoint(id1)
                    vals = [M.GetTuple(kl),N.GetTuple(kl),T.GetTuple(kl)]
                    segments.append([project_on_plane(base,origin,pt0),
                                     project_on_plane(base,origin,pt1),vals])
##                if id1>id0:
##        ##            v1 = M.GetTuple(id1)[0]
##                    ax.plot([v0,v0],[pt0[1],pt1[1]],'k')
        
    return segments

def contourf(ax,val,crd,output,loc_syst,orig,levels=0):
    loc = locator(output,val,loc_syst,orig)

    # contour data:
    bounds = [[min(crd[0])-2,max(crd[0])+2],
              [min(crd[1])-2,max(crd[1])+2]]
    xnew = np.linspace(bounds[0][0],bounds[0][1],100)
    znew = np.linspace(bounds[1][0],bounds[1][1],100)
    X,Z = np.meshgrid(xnew,znew)
    V = np.zeros([len(X),len(X[0])])
    val1 = []
    for kx in range(len(X[0])):
        for kz in range(len(X)):
            V[kz][kx] = loc.interpolate(X[kz][kx],Z[kz][kx],0)
            if not np.isnan(V[kz][kx]):
                val1.append(V[kz][kx])
##            else:
##                val1.append(0)
    try:
        if levels==0:
            levels = np.linspace(min(val1),max(val1),10)
    except:
        pass
    CS = ax.contourf(X,Z,V,levels,extend='both')

    return CS

def contourf_curved_section(ax,val,crd,output,levels=0,alpha=0):
    loc = locator_curved(output,val,alpha)

    # contour data:
    bounds = [[min(crd[0]),max(crd[0])],
              [min(crd[1]),max(crd[1])]]
    arr_crd_2d = output.GetPointData().GetArray('2Dcrds')
    xy = [[arr_crd_2d.GetTuple(kk)[0] for kk in range(arr_crd_2d.GetNumberOfTuples())],
          [arr_crd_2d.GetTuple(kk)[1] for kk in range(arr_crd_2d.GetNumberOfTuples())]]
    bounds = [[min(xy[0]),max(xy[0])],
              [min(xy[1]),max(xy[1])]]
    xnew = np.linspace(bounds[0][0],bounds[0][1],100)
    znew = np.linspace(bounds[1][0],bounds[1][1],100)
    X,Z = np.meshgrid(xnew,znew)
    V = np.zeros([len(X),len(X[0])])
    val1 = []
    for kx in range(len(X[0])):
        for kz in range(len(X)):
            V[kz][kx] = loc.interpolate(X[kz][kx],Z[kz][kx],0)
##            print(X[kz][kx],Z[kz][kx],V[kz][kx])
            if not np.isnan(V[kz][kx]):
                val1.append(V[kz][kx])
##            else:
##                val1.append(0)
    try:
        if levels==0:
            levels = np.linspace(min(val1),max(val1),10)
    except:
        pass
    CS = ax.contourf(X,Z,V,levels,extend='both')

    return CS

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

class locator:
    def __init__(self,pd,val,loc_syst,orig):
        self.loc = vtk.vtkCellLocator()
        self.loc = vtk.vtkModifiedBSPTree()
        self.pd2 = vtk.vtkPolyData()
        self.pd2.DeepCopy(pd)
        
        pts = vtk.vtkPoints()
        for kp in range(pd.GetNumberOfPoints()):
            p = pd.GetPoint(kp)
            crd_2D = project_on_plane(loc_syst,orig,p)
            pts.InsertNextPoint(crd_2D[0],crd_2D[1],0)
        self.pd2.SetPoints(pts)
        self.loc.SetDataSet(self.pd2)
        self.loc.BuildLocator()
        self.values = val        

    def interpolate(self,x,y,z):
        aCell = vtk.vtkTriangle()
        pcoords = [0,0,0]
        weights = [0,0,0]
        
        cid = self.loc.FindCell((x,y,z))
        if cid>-1:
            aCell = self.pd2.GetCell(cid)
            pts = [(aCell.GetPoints().GetPoint(k)[0],
                    aCell.GetPoints().GetPoint(k)[1]) for k in range(3)]
            res = aCell.BarycentricCoords((x,y),pts[0],pts[1],pts[2],pcoords)
            v = [self.values[aCell.GetPointId(k)] for k in range(3)]
            return sum([pcoords[kk]*v[kk] for kk in range(3)])

class locator_curved:
    def __init__(self,pd,val,alpha=0):
        self.loc = vtk.vtkCellLocator()
        self.loc = vtk.vtkModifiedBSPTree()
        
        arr_crd_2d = pd.GetPointData().GetArray('2Dcrds')
        
        pts = vtk.vtkPoints()
        for kp in range(pd.GetNumberOfPoints()):
            p = pd.GetPoint(kp)
            crd_2D = arr_crd_2d.GetTuple(kp)
            pts.InsertNextPoint(crd_2D[0],crd_2D[1],0)
        pd.SetPoints(pts)

        del2d = vtk.vtkDelaunay2D()
        del2d.SetInputData(pd)
        del2d.SetTolerance(0)
        if alpha>0:
            del2d.SetAlpha(alpha)
        del2d.Update()
        self.pd2 = del2d.GetOutput()
##        w = vtk.vtkXMLPolyDataWriter()
##        w.SetFileName('del.vtp')
##        w.SetInputData(self.pd2)
##        w.Write()
##        sys.exit(0)
        
        self.loc.SetDataSet(self.pd2)
        self.loc.BuildLocator()
        self.values = val        

    def interpolate(self,x,y,z):
        aCell = vtk.vtkTriangle()
        pcoords = [0,0,0]
        weights = [0,0,0]
        
        cid = self.loc.FindCell((x,y,z))
        if cid>-1:
            aCell = self.pd2.GetCell(cid)
            if aCell.GetClassName()=='vtkTriangle':
                pts = [(aCell.GetPoints().GetPoint(k)[0],
                        aCell.GetPoints().GetPoint(k)[1]) for k in range(3)]
                res = aCell.BarycentricCoords((x,y),pts[0],pts[1],pts[2],pcoords)
                v = [self.values[aCell.GetPointId(k)] for k in range(3)]
                return sum([pcoords[kk]*v[kk] for kk in range(3)])
##        else:
##            return 0

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

def get_legend(lut,categories={},hfrac=0.8,vpad=0.12,hwratio=4,dpi=96,
               label='Label [Unit]',ncolors=30):

    from matplotlib.patches import Rectangle
    from PIL import Image
    
    fig = plt.figure(figsize=(4.5,hwratio*4.5),dpi=dpi)
    bbox_props = dict(boxstyle='round',fc="w", ec="k", lw=2)

    ax = fig.add_subplot(111)

    if len(categories.keys()):
        # Legend with categories (materials)
        keys = categories.keys()
        nval = len(keys)
        dh = min(hfrac/(ncolors*(1-vpad)+(ncolors+1)*vpad)*(1-vpad),1./(nval*(1-vpad)+(nval+1)*vpad)*(1-vpad))
        ds = dh/(1-vpad)*vpad
        for kc in range(nval):
            p = Rectangle((2*ds,1-(kc+1)*(dh+ds)),2-4*ds,dh,fc=lut.GetTableValue(keys[kc]%20-1),ec='k')
            ax.add_patch(p)
            ax.text(0.1,1-(kc+1)*(dh+ds)+dh/3,str(keys[kc])+': '+categories[keys[kc]],
                    va='center',size=int(600*dh))
    else:
        # Colormap legend
        nval = lut.GetNumberOfValues()
        vpad *= 5
        dh = min(hfrac/(ncolors*(1-vpad)+(ncolors+1)*vpad)*(1-vpad),1./(nval*(1-vpad)+(nval+1)*vpad)*(1-vpad))
        ds = dh/(1-vpad)*vpad
        vrange = lut.GetRange()
        for kc in range(nval):
            val = vrange[0] + (vrange[1]-vrange[0])/nval*kc
            col = [0,0,0]
            lut.GetColor(val,col)
            p = Rectangle((2*ds,1-(kc+1)*(dh+ds)),2-4*ds,dh,fc=col,
                          ec='None')
            ax.add_patch(p)
            ax.text(0.1,1-(kc+1)*(dh+ds)-ds/2,'%1.2e'%(val),
                    va='center',size=20)
    ax.annotate(label,size=24,
                xy=(0.5,1-hfrac-0.02),xycoords='axes fraction',ha='center',va='bottom')
    ax.annotate('Model created\nand computed\nwith ZSoil',size=24,
                xy=(0.5,0.01),xycoords='axes fraction',ha='center',va='bottom')
    ax.set_xlim(0,2)
    ax.set_ylim(0,1)
    ax.axis('off')
    fig.tight_layout()
    fig.savefig('legend.png')
    leg = Image.open('legend.png')

    return leg
