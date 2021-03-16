import meep as mp
from meep import mpb
import numpy as np
import h5py

num_bands = 6

k_points = [
            mp.Vector3(0.5, 0.5),  # M
            mp.Vector3(),          # Gamma
            mp.Vector3(0.5),       # X
            mp.Vector3(0.5, 0.5)    # M
            ]

k_points = mp.interpolate(16, k_points)

d = 0.433


crystal_material=mp.Medium(epsilon = 13.7)

geometry = [mp.Block(size=mp.Vector3(d, d, mp.inf), material=crystal_material)]

geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1))

resolution = 16

ms = mpb.ModeSolver(num_bands=num_bands,
                    default_material=mp.Medium(epsilon=2.25),
                    k_points=k_points,
                    geometry=geometry,
                    geometry_lattice=geometry_lattice,
                    resolution=resolution)



usebloch = True

def get_fields(ms, band):
    k_index = k_points.index(ms.current_k)
    hf = h5py.File('fielddata_band' + str(band) + '_kpoint' +str(k_index) + '.h5', 'w')
    efield = ms.get_efield(band, bloch_phase = usebloch)
    hf.create_dataset('efield_x_real',data = np.real(np.squeeze(efield[:,:,:,0])))
    hf.create_dataset('efield_x_imag',data = np.imag(np.squeeze(efield[:,:,:,0])))
    hf.create_dataset('efield_y_real',data = np.real(np.squeeze(efield[:,:,:,1])))
    hf.create_dataset('efield_y_imag',data = np.imag(np.squeeze(efield[:,:,:,1])))
    hf.create_dataset('efield_z_real',data = np.real(np.squeeze(efield[:,:,:,2])))
    hf.create_dataset('efield_z_imag',data = np.imag(np.squeeze(efield[:,:,:,2])))
    dfield = ms.get_dfield(band, bloch_phase = usebloch)
    hf.create_dataset('dfield_x_real',data = np.real(np.squeeze(dfield[:,:,:,0])))
    hf.create_dataset('dfield_x_imag',data = np.imag(np.squeeze(dfield[:,:,:,0])))
    hf.create_dataset('dfield_y_real',data = np.real(np.squeeze(dfield[:,:,:,1])))
    hf.create_dataset('dfield_y_imag',data = np.imag(np.squeeze(dfield[:,:,:,1])))
    hf.create_dataset('dfield_z_real',data = np.real(np.squeeze(dfield[:,:,:,2])))
    hf.create_dataset('dfield_z_imag',data = np.imag(np.squeeze(dfield[:,:,:,2])))
    hfield = ms.get_hfield(band, bloch_phase = usebloch)
    hf.create_dataset('hfield_x_real',data = np.real(np.squeeze(hfield[:,:,:,0])))
    hf.create_dataset('hfield_x_imag',data = np.imag(np.squeeze(hfield[:,:,:,0])))
    hf.create_dataset('hfield_y_real',data = np.real(np.squeeze(hfield[:,:,:,1])))
    hf.create_dataset('hfield_y_imag',data = np.imag(np.squeeze(hfield[:,:,:,1])))
    hf.create_dataset('hfield_z_real',data = np.real(np.squeeze(hfield[:,:,:,2])))
    hf.create_dataset('hfield_z_imag',data = np.imag(np.squeeze(hfield[:,:,:,2])))
    bfield = ms.get_bfield(band, bloch_phase = usebloch)
    hf.create_dataset('bfield_x_real',data = np.real(np.squeeze(bfield[:,:,:,0])))
    hf.create_dataset('bfield_x_imag',data = np.imag(np.squeeze(bfield[:,:,:,0])))
    hf.create_dataset('bfield_y_real',data = np.real(np.squeeze(bfield[:,:,:,1])))
    hf.create_dataset('bfield_y_imag',data = np.imag(np.squeeze(bfield[:,:,:,1])))
    hf.create_dataset('bfield_z_real',data = np.real(np.squeeze(bfield[:,:,:,2])))
    hf.create_dataset('bfield_z_imag',data = np.imag(np.squeeze(bfield[:,:,:,2])))

    hf.close

ms.run_tm(mpb.fix_efield_phase, mpb.fix_hfield_phase, mpb.fix_dfield_phase, mpb.fix_bfield_phase, get_fields)

banddata = ms.all_freqs

np.savetxt("bandstructure.txt", banddata)


