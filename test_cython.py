# test cython
from chill_lib import PyChiller

import ovito
from ovito.io import import_file
#from ovito.data import NearestNeighborFinder
#from ovito.data import CutoffNeighborFinder
from ovito.data import ParticleProperty
from ovito.modifiers import *
import numpy as np
#from scipy.special import sph_harm


def modify(frame, input, output):
    #states = chill_plus(output.particle_properties.position.array)
    cell = output.cell.matrix.copy()
    cell = cell.transpose().flatten()
    cell = cell.astype(np.float64)
    numParticles = len(output.particle_properties.position.array)
    myChiller = PyChiller(  output.particle_properties.position.array.astype(np.float64), \
                            numParticles, \
                            cell)
    status = np.asarray(myChiller.get_status())
    output.create_user_particle_property("Chill Status", "float")
    output['Chill Status'].marray[:] = status[:];
    output['Chill Status'].changed()
    #color_property = output.create_particle_property(ParticleProperty.Type.Color)
    #color_property.marray[:] = np.asarray(((status==0).astype(np.float64), (status==3).astype(np.float64), (status==4).astype(np.float64))).transpose()
    #print (states)


if __name__ == "__main__":
    filename = "/media/henriasv/IcyBox/phd_methane_hydrates/pennyshaped_cracks_variable_methane/systematic_mw_pennycracks_VariableMethaneConcentrationInCavity-2016-05-18-154156/lmp_Nthermalize=2000.0_Nerate=10000.0_CavityMethaneLatticeSpacing=8.0_temperature=250.0_crackRadius=40.0_Nproduction=40000.0_timeStep=10.0_Nx=24_Ny=24_Nz=24_crackHeight=6.0_maxStrain=1.0545_seed=000/trajectory.lammpstrj"
    #filename = "/work/users/henriasv/molecular-simulations/systematic_mw_pennycracks_VariableMethaneConcentrationInCavity-2016-05-09-124346/lmp_Nthermalize=200.0_Nerate=1000.0_temperature=260.0_crackRadius=40.0_Nproduction=4000.0_timeStep=10.0_Nx=24_Ny=24_Nz=24_crackHeight=6.0_maxStrain=1.072_seed=000/trajectory.lammpstrj"
    node = import_file(filename, multiple_frames=True)
    node.compute()

    node.modifiers.append(SelectParticleTypeModifier(property='Particle Type', types={1}))
    node.modifiers.append(DeleteSelectedParticlesModifier())

    slices = [[5, 295], [5, 295], [144, 160]]
    offset = 6

    sliceModifiers2 = [  SliceModifier(distance=slices[0][0], inverse=True), \
                        SliceModifier(distance=slices[0][1]), \
                        SliceModifier(distance=slices[1][0],normal=[0,1,0], inverse=True), \
                        SliceModifier(distance=slices[1][1],normal=[0,1,0]), \
                        SliceModifier(distance=slices[2][0],normal=[0,0,1], inverse=True), \
                        SliceModifier(distance=slices[2][1],normal=[0,0,1]) ]

    sliceModifiers = [ SliceModifier(distance=slices[0][0]-offset, inverse=True), \
                        SliceModifier(distance=slices[0][1]+offset), \
                        SliceModifier(distance=slices[1][0]-offset,normal=[0,1,0], inverse=True), \
                        SliceModifier(distance=slices[1][1]+offset,normal=[0,1,0]), \
                        SliceModifier(distance=slices[2][0]-offset,normal=[0,0,1], inverse=True), \
                        SliceModifier(distance=slices[2][1]+offset,normal=[0,0,1]) ]

    #for sliceModifier in sliceModifiers:
    #    node.modifiers.append(sliceModifier)
    node.modifiers.append(SliceModifier(distance=slices[2][0]-offset,normal=[0,0,1], inverse=True))
    node.modifiers.append(SliceModifier(distance=slices[2][1]+offset,normal=[0,0,1]))
    node.modifiers.append(PythonScriptModifier(function=modify))
    node.modifiers.append(SliceModifier(distance=slices[2][0],normal=[0,0,1], inverse=True))
    node.modifiers.append(SliceModifier(distance=slices[2][1],normal=[0,0,1]))

    n_bins = 48
    last_frame = ovito.dataset.anim.last_frame
    modifier = BinAndReduceModifier();
    modifier.direction = BinAndReduceModifier.Direction.Vector_2
    modifier.property = "Chill Status"
    modifier.bin_count_x = n_bins
    node.modifiers.append(modifier)
    color_modifier = ColorCodingModifier(property="Chill Status", gradient=ColorCodingModifier.Rainbow(), start_value=0, end_value=4)
    node.modifiers.append(color_modifier)


    output_data = np.zeros((last_frame, n_bins))
    for frame in range(0, last_frame, 10):
        ovito.dataset.anim.current_frame = frame
        node.compute()
        output_data[frame, :] = modifier.bin_data[:]
        print(frame)
    output_data.tofile("data.bin")
