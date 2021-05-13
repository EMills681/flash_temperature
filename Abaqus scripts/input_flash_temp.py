# -*- coding: mbcs -*-
#### Open Abaqus - set the work directory to the location where these python files are stored (file -> set work directory...)
#### Run this script in Abaqus (file -> run script... -> input_flash_temp)

import geofun
import flashtfun
import math

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

from abaqus import getInput

#User inputs required parameters
while True:
    try:
        m = float(getInput('Enter gear module (mm):','2'))
        break
    except ValueError:
        print("Invalid input, please enter a number")

while True:
    try:
        z1 = int(getInput("Enter number of teeth: ",'30'))
        break
    except ValueError:
        print("Invalid input, please enter an integer")

while True:
    try:
        pressure_angle = float(getInput("Enter pressure angle (degrees): ",'20'))
        break
    except ValueError:
        print("Invalid input, please enter a number")

while True:
    try:
        b = float(getInput("Enter face width (mm): ",'17.4'))
        break
    except ValueError:
        print("Invalid input, please enter a number")

while True:
    try:
        T1 = float(getInput("Enter torque (Nm): ",'8.0'))
        break
    except ValueError:
        print("Invalid input, please enter a number")

while True:
    try:
        N1 = float(getInput("Enter rotational speed (rpm): ",'1000'))
        break
    except ValueError:
        print("Invalid input, please enter a number")

while True:
    try:
        mu = float(getInput("Enter coefficient of friction (experimental): ",'0.35'))
        break
    except ValueError:
        print("Invalid input, please enter a number")

#other inputs
## Gear parameters
mg = 1
psi = math.radians(pressure_angle)

## Step
step = 11            # number of partitians - DO NOT CHANGE WITHOUT ADJUSING CODE seq's!!!
start_temp = 20.0    # gear temp at start of simulation
surround_temp = 20.0 # surrounding temp

## Material properties (Material: POM)
E = 3590                    # Modulus of elasticity (N/mm2 = MPa)
v = 0.34                    # Poisson's ratio
hc = 0.230                  # Heat conductivity (N/s/C)
rhoM = 1420                 # Density (kg/m3)
cM = 1285                   # Specific heat per unit mass (J/kg/C)

## Air properties
rhoair = 1.204              # Density of air (kg/m3) @ 20 oC
vair = 1.825 * 10**-5       # Dynamic viscocity of air (kg/m/s) @ 20 oC
Cpair = 1007                # Specific heat capacity of air (J/kg/C) @ 20 oC
kair = 0.02514              # Thermal conductivity of air (W/m/C) @ 20 oC
h = flashtfun.heat_trans_co(m, z1, N1, rhoair, vair, Cpair, kair)

#Flash temp function values (flashtfun.py)
datum = flashtfun.datum_points(m, z1, psi, step)
temp = flashtfun.temps_and_time(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM, mu)[:,0]
maxtemp = flashtfun.tmax(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM, mu) + start_temp
time = min(flashtfun.temps_and_time(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM, mu)[:,1])
globalseed = 0.8 * min(flashtfun.Hertz_contact(m, z1, psi, b, step, mg, E, v, T1)) ##minimum Hertzian half width * 0.9

###Geometry coordinates###
spline1 = geofun.spline1(m, z1, psi)
spline2 = geofun.spline2(m, z1, psi)
arcs = geofun.dedendum_coords(m,z1, psi)
clearance = geofun.clearance(m)

########################################################################################

##sketch tooth shape
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].Spline(points=(geofun.array_to_tuple(spline1)))
mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0.0, 0.0)
    , direction=CLOCKWISE, point1=(spline2[-1,0], spline2[-1,1]), point2=(spline1[-1,0], spline1[-1,1]))
mdb.models['Model-1'].sketches['__profile__'].Spline(points=(geofun.array_to_tuple(spline2)))

## add clearance section
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spline1[0,0], spline1[0,1]), point2=(arcs[0,0], arcs[0,1]))
mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0.0, 0.0)
    , direction=CLOCKWISE, point1=(arcs[0,0], arcs[0,1]), point2=(arcs[1,0], arcs[1,1]))

mdb.models['Model-1'].sketches['__profile__'].Line(point1=(spline2[0,0], spline2[0,1]), point2=(arcs[2,0], arcs[2,1]))
mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0.0, 0.0)
    , direction=CLOCKWISE, point1=(arcs[3,0], arcs[3,1]), point2=(arcs[2,0], arcs[2,1]))

## add 0.1 bulk of gear #### INCREASE BETWEEN 0.1 and 1 if you want to show more of the gear
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(arcs[1,0], arcs[1,1]), point2=(0.9*arcs[1,0], 0.9*arcs[1,1]))
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0.9*arcs[3,0], 0.9*arcs[3,1]), point2=(arcs[3,0], arcs[3,1]))
mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0.0, 0.0)
    , direction=CLOCKWISE, point1=(0.9*arcs[3,0], 0.9*arcs[3,1]), point2=(0.9*arcs[1,0], 0.9*arcs[1,1]))

## create part from sketch
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseSolidExtrude(depth=b, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

## add fillet
mdb.models['Model-1'].parts['Part-1'].Round(edgeList=(
    mdb.models['Model-1'].parts['Part-1'].edges[3], ), radius= clearance)
mdb.models['Model-1'].parts['Part-1'].Round(edgeList=(
    mdb.models['Model-1'].parts['Part-1'].edges[17], ), radius= clearance)

###  Model Attributes ###
mdb.models['Model-1'].setValues(absoluteZero=-273.15, stefanBoltzmann=5.6e-11)

### MATERIAL (POM) ###
mdb.models['Model-1'].Material(name='POM')
mdb.models['Model-1'].materials['POM'].Density(table=((rhoM * 10**-12, ), ))
mdb.models['Model-1'].materials['POM'].Elastic(table=((E, v), ))
mdb.models['Model-1'].materials['POM'].Conductivity(table=((hc, ), ))
mdb.models['Model-1'].materials['POM'].SpecificHeat(table=((cM * 10**6, ), ))
mdb.models['Model-1'].HomogeneousSolidSection(material='POM', name='Section-1', 
    thickness=None)
mdb.models['Model-1'].parts['Part-1'].Set(cells=
    mdb.models['Model-1'].parts['Part-1'].cells.getSequenceFromMask(('[#1 ]', 
    ), ), name='Set-1')
mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Part-1'].sets['Set-1'], sectionName=
    'Section-1', thicknessAssignment=FROM_SECTION)

### DATUM Splits
for i in range(step-1):
    mdb.models['Model-1'].parts['Part-1'].DatumPointByCoordinate(coords=(
        datum[i,0], datum[i,1], 0.0))
    mdb.models['Model-1'].parts['Part-1'].DatumPointByCoordinate(coords=(
        datum[i,0], datum[i,1], b))

########### this will change if step changes!!############################
seq = ('#200', '#400', '#800', '#1000', '#2000', '#4000', '#8000', '#10000', '#20000', '#40000')
##########################################################################
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='Part-1-1', 
    part=mdb.models['Model-1'].parts['Part-1'])

for i in range(step-1):
    mdb.models['Model-1'].parts['Part-1'].PartitionFaceByShortestPath(faces=
        mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask((seq[i], 
        ), ), point1=mdb.models['Model-1'].parts['Part-1'].datums[2*i+5], point2=
        mdb.models['Model-1'].parts['Part-1'].datums[2*i+6])

### Initial Conditions
mdb.models['Model-1'].rootAssembly.regenerate()
mdb.models['Model-1'].rootAssembly.Set(cells=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].cells.getSequenceFromMask(
    ('[#1 ]', ), ), edges=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.getSequenceFromMask(
    ('[#ffffffff #fffffff ]', ), ), faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
    ('[#3fffff ]', ), ), name='Set-1', vertices=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].vertices.getSequenceFromMask(
    ('[#ffffffff #ff ]', ), ))
mdb.models['Model-1'].Temperature(createStepName='Initial', 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
    UNIFORM, magnitudes=(start_temp, ), name='Predefined Field-1', region=
    mdb.models['Model-1'].rootAssembly.sets['Set-1'])

### BC's for each temperature step **
set = ()
bc = ()

for i in range(step):
    set = set + ('Set-' + str(i+2),)
    bc = bc + ('BC-' + str(i+1),)

########### this will need to change if step changes!!############################
seq = ('#100', '#80', '#40', '#20', '#10', '#8', '#4', '#2', '#1', '#80000')
##########################################################################

##BC's / steps across the whole tooth face
mdb.models['Model-1'].HeatTransferStep(deltmx= maxtemp, initialInc= time, 
    maxInc= 2 * (time * 10 **4), minInc=time/2, name= 'Heating', previous= 'Initial', 
    timePeriod=1200)
mdb.models['Model-1'].steps['Heating'].setValues(maxNumInc=500)
for i in range(step-1):
    mdb.models['Model-1'].rootAssembly.Set(faces=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
        (seq[i], ), ), name= set[i])
    mdb.models['Model-1'].TemperatureBC(amplitude=UNSET, createStepName= 'Heating', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude= temp[i] + start_temp, name=
        bc[i], region=mdb.models['Model-1'].rootAssembly.sets[set[i]])

### Interactions
#tip
mdb.models['Model-1'].rootAssembly.Surface(name='Surf-1', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
    ('[#40000 ]', ), ))
mdb.models['Model-1'].FilmCondition(createStepName='Heating', definition=
    EMBEDDED_COEFF, filmCoeff= h[0], filmCoeffAmplitude='', name='Int-1', 
    sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='', 
    sinkTemperature= surround_temp, surface=
    mdb.models['Model-1'].rootAssembly.surfaces['Surf-1'])
#sides
mdb.models['Model-1'].rootAssembly.Surface(name='Surf-2', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
    ('[#300000 ]', ), ))
mdb.models['Model-1'].FilmCondition(createStepName='Heating', definition=
    EMBEDDED_COEFF, filmCoeff= h[1], filmCoeffAmplitude='', name='Int-2', 
    sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='', 
    sinkTemperature= surround_temp, surface=
    mdb.models['Model-1'].rootAssembly.surfaces['Surf-2'])
#non-meshing tooth
mdb.models['Model-1'].rootAssembly.Surface(name='Surf-3', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
    ('[#23e00 ]', ), ))
mdb.models['Model-1'].FilmCondition(createStepName='Heating', definition=
    EMBEDDED_COEFF, filmCoeff= h[2], filmCoeffAmplitude='', name='Int-3', 
    sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='', 
    sinkTemperature= surround_temp, surface=
    mdb.models['Model-1'].rootAssembly.surfaces['Surf-3'])

### CREATE A BASIC MESH
mdb.models['Model-1'].rootAssembly.setElementType(elemTypes=(ElemType(
    elemCode=DCC3D8, elemLibrary=STANDARD), ElemType(elemCode=DC3D6, 
    elemLibrary=STANDARD), ElemType(elemCode=DC3D4, elemLibrary=STANDARD)), 
    regions=(
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].cells.getSequenceFromMask(
    ('[#1 ]', ), ), ))
mdb.models['Model-1'].rootAssembly.seedPartInstance(deviationFactor=0.1, 
    minSizeFactor=0.1, regions=(
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], ), size=globalseed)
mdb.models['Model-1'].rootAssembly.generateMesh(regions=(
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], ))


print("Gear module: ")
print(m)
print("Number of teeth: ")
print(z1)
print("Pressure angle: ")
print(pressure_angle)
print("Face width: ")
print(b)
print("Torque: ")
print(T1)
print("Rotational speed: ")
print(N1)
print("Pressure angle: ")
print(mu)
print("Minimum Hertzian half width: ")
print(min(flashtfun.Hertz_contact(m, z1, psi, b, step, mg, E, v, T1)))
print("tip surface heat transfer coefficient: ")
print(h[0])
print("side surfaces heat transfer coefficient: ")
print(h[1])
print("non-meshing tooth surfaces heat transfer coefficient: ")
print(h[2])