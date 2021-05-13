# -*- coding: mbcs -*-
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

# mesh is not created in this script!

#User inputs required parameters
while True:
    try:
        m = float(getInput('Enter gear module:','2'))
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
        pressure_angle = float(getInput("Enter pressure angle: ",'20'))
        break
    except ValueError:
        print("Invalid input, please enter a number")

while True:
    try:
        b = float(getInput("Enter face width: ",'17.4'))
        break
    except ValueError:
        print("Invalid input, please enter a number")

while True:
    try:
        T1 = float(getInput("Enter pinion torque (Nm): ",'8.0'))
        break
    except ValueError:
        print("Invalid input, please enter a number")

while True:
    try:
        N1 = float(getInput("Enter pinion rotational speed (rpm): ",'1000'))
        break
    except ValueError:
        print("Invalid input, please enter a number")

while True:
    try:
        mu = float(getInput("Enter coefficient of friction (experimental): ",'0.35'))
        break
    except ValueError:
        print("Invalid input, please enter a number")

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

#Temperature values
datum = flashtfun.datum_points(m, z1, psi, step)
flux = flashtfun.heat_flux(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM, mu)
temp = flashtfun.temps_and_time(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM, mu)[:,0] + start_temp
maxtemp = flashtfun.tmax(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM, mu) + start_temp
time = flashtfun.temps_and_time(m, z1, T1, b, N1, psi, step, mg, E, v, hc, rhoM, cM, mu)[:,1]
rest_t = flashtfun.rotation_time(m, z1, psi, step, N1)

###Create Geometry ###
spline1 = geofun.spline1(m, z1, psi)
spline2 = geofun.spline2(m, z1, psi)
arcs = geofun.dedendum_coords(m,z1, psi)
clearance = geofun.clearance(m)

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

## add 0.1 bulk of gear
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

### Steps and BC's for each step
name = ()
name = name + ('Initial',)
set = ()
bc = ()

for i in range(step):
    name = name + ('Heat' + str(i+1),)
    set = set + ('Set-' + str(i+1),)
    bc = bc + ('BC-' + str(i+1),)

########### this will change if step changes!!############################
seq = ('#100', '#80', '#40', '#20', '#10', '#8', '#4', '#2', '#1', '#80000')
##########################################################################

##BC's / steps across the whole tooth face
for i in range(1,step):
    mdb.models['Model-1'].HeatTransferStep(deltmx= maxtemp, initialInc= time[i-1], 
        maxInc= time[i-1], minInc=(time[i-1] * 10**-5), name= name[i], previous= name[i-1], 
        timePeriod=time[i-1])
    mdb.models['Model-1'].rootAssembly.Set(faces=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
        (seq[i-1], ), ), name= set[i])
    mdb.models['Model-1'].TemperatureBC(amplitude=UNSET, createStepName= name[i], 
        distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude= temp[i-1], name=
        bc[i-1], region=mdb.models['Model-1'].rootAssembly.sets[set[i]])
    if i > 1:
        mdb.models['Model-1'].boundaryConditions[bc[i-2]].deactivate(name[i])

## Break during 'rotation' before meshing again
mdb.models['Model-1'].HeatTransferStep(deltmx= maxtemp, initialInc=rest_t, 
    maxInc= rest_t, minInc=(rest_t * 10**-5), name='Rotation1', previous= name[step-1], timePeriod=
    rest_t)
mdb.models['Model-1'].boundaryConditions[bc[step-2]].deactivate('Rotation1')

## repeat BC steps
x = 1
y = 1
for i in range(step,91):

    if ((i-1)/10.0).is_integer():
        x = 1
        mdb.models['Model-1'].HeatTransferStep(deltmx= maxtemp, initialInc= time[x-1], 
        maxInc= time[x-1], minInc=(time[x-1] * 10**-5), name= ('Heat' + str(i)), previous= ('Rotation' + str(y)), 
        timePeriod=time[x-1])
        mdb.models['Model-1'].rootAssembly.Set(faces=
            mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
            (seq[x-1], ), ), name= ('Set-' + str(i+1)))
        mdb.models['Model-1'].TemperatureBC(amplitude=UNSET, createStepName= ('Heat' + str(i)), 
            distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude= temp[x-1], name=
            ('BC-' + str(i)), region=mdb.models['Model-1'].rootAssembly.sets[('Set-' + str(i+1))])        
        y = y + 1

    else:
        mdb.models['Model-1'].HeatTransferStep(deltmx= maxtemp, initialInc= time[x-1], 
            maxInc= time[x-1], minInc=(time[x-1] * 10**-5), name= ('Heat' + str(i)), previous= ('Heat' + str(i-1)), 
            timePeriod=time[x-1])
        mdb.models['Model-1'].rootAssembly.Set(faces=
            mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
            (seq[x-1], ), ), name= ('Set-' + str(i+1)))
        mdb.models['Model-1'].TemperatureBC(amplitude=UNSET, createStepName= ('Heat' + str(i)), 
            distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude= temp[x-1], name=
            ('BC-' + str(i)), region=mdb.models['Model-1'].rootAssembly.sets[('Set-' + str(i+1))])
        mdb.models['Model-1'].boundaryConditions[('BC-' + str(i-1))].deactivate(('Heat' + str(i)))

    if (i/10.0).is_integer():
        mdb.models['Model-1'].HeatTransferStep(deltmx= maxtemp, initialInc=rest_t, 
            maxInc= rest_t, minInc=(rest_t * 10**-5), name=('Rotation' + str(y)), previous= ('Heat' + str(i)), timePeriod=
            rest_t)
        mdb.models['Model-1'].boundaryConditions[('BC-' + str(i))].deactivate('Rotation' + str(y))        

    x = x + 1


### Interactions
#tip
mdb.models['Model-1'].rootAssembly.Surface(name='Int-Surf-1', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
    ('[#40000 ]', ), ))
mdb.models['Model-1'].FilmCondition(createStepName='Heat1', definition=
    EMBEDDED_COEFF, filmCoeff= h[0], filmCoeffAmplitude='', name='Int-1', 
    sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='', 
    sinkTemperature= surround_temp, surface=
    mdb.models['Model-1'].rootAssembly.surfaces['Int-Surf-1'])
#sides
mdb.models['Model-1'].rootAssembly.Surface(name='Int-Surf-2', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
    ('[#300000 ]', ), ))
mdb.models['Model-1'].FilmCondition(createStepName='Heat1', definition=
    EMBEDDED_COEFF, filmCoeff= h[1], filmCoeffAmplitude='', name='Int-2', 
    sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='', 
    sinkTemperature= surround_temp, surface=
    mdb.models['Model-1'].rootAssembly.surfaces['Int-Surf-2'])
#non-meshing tooth
mdb.models['Model-1'].rootAssembly.Surface(name='Int-Surf-3', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.getSequenceFromMask(
    ('[#23e00 ]', ), ))
mdb.models['Model-1'].FilmCondition(createStepName='Heat1', definition=
    EMBEDDED_COEFF, filmCoeff= h[2], filmCoeffAmplitude='', name='Int-3', 
    sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='', 
    sinkTemperature= surround_temp, surface=
    mdb.models['Model-1'].rootAssembly.surfaces['Int-Surf-3'])