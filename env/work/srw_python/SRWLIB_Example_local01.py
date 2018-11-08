from srwlib import *
from uti_plot import *
import numpy

if not srwl_uti_proc_is_master(): exit()

####################################################
# LIGHT SOURCE
__num_x=4000
__num_y=4000

mesh = SRWLRadMesh(_eStart=0.11698,
                   _eFin  =0.11698,
                   _ne    =1,
                   _xStart=-0.05,
                   _xFin  =0.05,
                   _nx    =__num_x,
                   _yStart=-0.05,
                   _yFin  =0.05,
                   _ny    =__num_y,
                   _zStart=0.0)

stk = SRWLStokes()
stk.allocate(1,__num_x,__num_y)
stk.mesh = mesh

wfr = SRWLWfr()
wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
wfr.mesh = mesh

initial_mesh = deepcopy(wfr.mesh)

GsnBm = SRWLGsnBm()
GsnBm.x         = 0.0
GsnBm.y         = 0.0
GsnBm.z         = 0.0
GsnBm.xp        = 0.0
GsnBm.yp        = 0.0
GsnBm.avgPhotEn = 0.11698
GsnBm.pulseEn   = 0.001
GsnBm.repRate   = 1
GsnBm.polar     = 1
GsnBm.sigX      = 0.02
GsnBm.sigY      = 0.02
GsnBm.sigT      = 1e-14
GsnBm.mx        = 0
GsnBm.my        = 0

wfr.partBeam.partStatMom1.x = GsnBm.x
wfr.partBeam.partStatMom1.y = GsnBm.y
wfr.partBeam.partStatMom1.z = GsnBm.z
wfr.partBeam.partStatMom1.xp = GsnBm.xp
wfr.partBeam.partStatMom1.yp = GsnBm.yp

srwl.CalcElecFieldGaussian(wfr, GsnBm, [0.0])

mesh0 = deepcopy(wfr.mesh)
arI = array('f', [0]*mesh0.nx*mesh0.ny)
srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, mesh0.eStart, 0, 0)
arIx = array('f', [0]*mesh0.nx)
srwl.CalcIntFromElecField(arIx, wfr, 6, 0, 1, mesh0.eStart, 0, 0)
arIy = array('f', [0]*mesh0.ny)
srwl.CalcIntFromElecField(arIy, wfr, 6, 0, 2, mesh0.eStart, 0, 0)
#save ascii file with intensity
#srwl_uti_save_intens_ascii(arI, mesh0, <file_path>)
plotMesh0x = [1000*mesh0.xStart, 1000*mesh0.xFin, mesh0.nx]
plotMesh0y = [1000*mesh0.yStart, 1000*mesh0.yFin, mesh0.ny]
uti_plot2d1d (arI, plotMesh0x, plotMesh0y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity Before Propagation'])

####################################################
# BEAMLINE

srw_oe_array = []
srw_pp_array = []

drift_before_oe_0 = SRWLOptD(7.4)
pp_drift_before_oe_0 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_0)
srw_pp_array.append(pp_drift_before_oe_0)

oe_0=SRWLOptA(_shape='c',
               _ap_or_ob='a',
               _Dx=0.05,
               _Dy=0.05,
               _x=0.0,
               _y=0.0)

pp_oe_0 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_0)
srw_pp_array.append(pp_oe_0)

oe_1=SRWLOptL(_Fx=7.4, _Fy=7.4, _x=0.0, _y=0.0)

pp_oe_1 = [0,0,1.0,4,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_1)
srw_pp_array.append(pp_oe_1)

drift_before_oe_2 = SRWLOptD(0.3)
pp_drift_before_oe_2 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_2)
srw_pp_array.append(pp_drift_before_oe_2)

acceptance_slits_oe_2=SRWLOptA(_shape='r',
               _ap_or_ob='a',
               _Dx=0.1,
               _Dy=0.09998489885977781,
               _x=0.0,
               _y=0.0)

oe_2 = SRWLOptMirPl(_size_tang=0.1414,
                     _size_sag=0.1,
                     _ap_shape='r',
                     _sim_meth=2,
                     _treat_in_out=1,
                     _nvx=0,
                     _nvy=0.7071067811865476,
                     _nvz=-0.7071067811865475,
                     _tvx=0,
                     _tvy=-0.7071067811865475,
                     _x=0.0,
                     _y=0.0)
oe_2.set_dim_sim_meth(_size_tang=0.1414,
                      _size_sag=0.1,
                      _ap_shape='r',
                      _sim_meth=2,
                      _treat_in_out=1)
oe_2.set_orient(_nvx=0,
                 _nvy=0.7071067811865476,
                 _nvz=-0.7071067811865475,
                 _tvx=0,
                 _tvy=-0.7071067811865475,
                 _x=0.0,
                 _y=0.0)


pp_acceptance_slits_oe_2 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]
pp_oe_2 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(acceptance_slits_oe_2)
srw_pp_array.append(pp_acceptance_slits_oe_2)

srw_oe_array.append(oe_2)
srw_pp_array.append(pp_oe_2)

drift_before_oe_3 = SRWLOptD(0.1)
pp_drift_before_oe_3 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_3)
srw_pp_array.append(pp_drift_before_oe_3)

acceptance_slits_oe_3=SRWLOptA(_shape='r',
               _ap_or_ob='a',
               _Dx=0.06,
               _Dy=0.05999093931586669,
               _x=0.0,
               _y=0.0)

oe_3 = SRWLOptMirPl(_size_tang=0.08484,
                     _size_sag=0.06,
                     _ap_shape='r',
                     _sim_meth=2,
                     _treat_in_out=1,
                     _nvx=0,
                     _nvy=0.7071067811865476,
                     _nvz=-0.7071067811865475,
                     _tvx=0,
                     _tvy=-0.7071067811865475,
                     _x=0.0,
                     _y=0.0)
oe_3.set_dim_sim_meth(_size_tang=0.08484,
                      _size_sag=0.06,
                      _ap_shape='r',
                      _sim_meth=2,
                      _treat_in_out=1)
oe_3.set_orient(_nvx=0,
                 _nvy=0.7071067811865476,
                 _nvz=-0.7071067811865475,
                 _tvx=0,
                 _tvy=-0.7071067811865475,
                 _x=0.0,
                 _y=0.0)


pp_acceptance_slits_oe_3 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]
pp_oe_3 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(acceptance_slits_oe_3)
srw_pp_array.append(pp_acceptance_slits_oe_3)

srw_oe_array.append(oe_3)
srw_pp_array.append(pp_oe_3)

oe_4=SRWLOptA(_shape='r',
               _ap_or_ob='o',
               _Dx=0.05,
               _Dy=0.0026,
               _x=0.0,
               _y=0.0)

pp_oe_4 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_4)
srw_pp_array.append(pp_oe_4)

drift_before_oe_5 = SRWLOptD(4.168)
pp_drift_before_oe_5 = [0,0,1.0,1,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_5)
srw_pp_array.append(pp_drift_before_oe_5)

oe_5=SRWLOptA(_shape='r',
               _ap_or_ob='a',
               _Dx=0.026,
               _Dy=0.03,
               _x=0.0,
               _y=0.0)

pp_oe_5 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_5)
srw_pp_array.append(pp_oe_5)

drift_before_oe_6 = SRWLOptD(0.5)
pp_drift_before_oe_6 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_6)
srw_pp_array.append(pp_drift_before_oe_6)

oe_6=SRWLOptA(_shape='r',
               _ap_or_ob='a',
               _Dx=0.01,
               _Dy=0.03,
               _x=0.0,
               _y=0.0)

pp_oe_6 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(oe_6)
srw_pp_array.append(pp_oe_6)

drift_before_oe_7 = SRWLOptD(2.332)
pp_drift_before_oe_7 = [0,0,1.0,0,0,1.0,1.0,1.0,1.0,0,0.0,0.0]

srw_oe_array.append(drift_before_oe_7)
srw_pp_array.append(pp_drift_before_oe_7)



####################################################
# PROPAGATION

optBL = SRWLOptC(srw_oe_array, srw_pp_array)
srwl.PropagElecField(wfr, optBL)

mesh1 = deepcopy(wfr.mesh)
arI1 = array('f', [0]*mesh1.nx*mesh1.ny)
srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, mesh1.eStart, 0, 0)
arI1x = array('f', [0]*mesh1.nx)
srwl.CalcIntFromElecField(arI1x, wfr, 6, 0, 1, mesh1.eStart, 0, 0)
arI1y = array('f', [0]*mesh1.ny)
srwl.CalcIntFromElecField(arI1y, wfr, 6, 0, 2, mesh1.eStart, 0, 0)
#save ascii file with intensity
#srwl_uti_save_intens_ascii(arI1, mesh1, <file_path>)
plotMesh1x = [1000*mesh1.xStart, 1000*mesh1.xFin, mesh1.nx]
plotMesh1y = [1000*mesh1.yStart, 1000*mesh1.yFin, mesh1.ny]
uti_plot2d1d(arI1, plotMesh1x, plotMesh1y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity After Propagation'])
uti_plot_show()
