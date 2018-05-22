import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpp

def compute_number_of_rows(R, well_sep, section_offset):
    """
    Compute number of rows that will fit in a section that is 1/6th of a
    circle with radius R, for a given well separation.
    """
    h_section = R - (2 * section_offset + well_sep)
    row_height = well_sep * np.sqrt(3)
    return int(np.floor(h_section / row_height)), row_height

def place_wells_in_section(R, well_sep, section_offset=0.1):
    """
    Compute the locations of wells withing a section given the total radius
    of the phantom, the distance between the wells, and the offset from the
    section boundary, expressed as a fraction of the total phantom radius.
    """
    section_offset *= R
    num_rows, row_height = compute_number_of_rows(R, well_sep, section_offset)
    # If only one well will fit in a section, try reducing the section offset
    # to see if more wells can be squeezed in
    if num_rows == 1:
        num_rows, row_height = compute_number_of_rows(R, well_sep, 0)
        if num_rows == 1:
            warnings.warn(("Cannot fit multiple features in section with "
                           "feature size = %s" %(well_sep)))
    # Compute well locations given the well separation and number of rows
    num_wells = np.sum(1 + np.arange(num_rows))
    xs, ys = [], []
    for i in range(num_rows):
        rn = i + 1
        for x in np.arange(-rn, rn, 2) + 1:
            xs.append(x * well_sep)
            ys.append(-(section_offset + row_height * rn))
    return np.vstack((xs, ys)).T

def export_to_G4mac(filehandle, x, y, z, r, energy, num_evs, halfz=None):
    """
    Write Derenzo wells to a Geant4 macro file for use with 
    G4GeneralParticleSource.
    """
    srcstr =  '/gps/particle gamma\n'
    # 2D vs 3D
    if halfz is None:
        srcstr += '/gps/pos/type Plane\n'
        srcstr += '/gps/pos/shape Circle\n'
    else:
        srcstr += '/gps/pos/type Volume\n'
        srcstr += '/gps/pos/shape Cylinder\n'
    srcstr += '/gps/pos/centre %f %f %f mm\n' %(x, y, z)
    srcstr += '/gps/pos/radius %f mm\n' %(r)
    if halfz is not None:
        srcstr += '/gps/pos/halfz %f mm\n' %(halfz)
    srcstr += '/gps/ang/type iso\n'
    srcstr += '/gps/ene/type Mono\n'
    srcstr += '/gps/ene/mono %f\n' %(energy)
    srcstr += '/run/beamOn %i\n\n' %(num_evs)
    filehandle.write(srcstr)

if __name__ == "__main__":

    # Create figure with aspect ration = 1.0
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    
    # Phantom params
    R = 50.0
    cyl_height = 10.0
    gamma_en = 661.657
    num_evs = 100
    
    # Add phantom outline to plot
    derenzo_cyl = mpp.Circle((0, 0), radius=R, color="gray", alpha=0.3)
    ax.add_patch(derenzo_cyl)
    ax.set_xlim((-1.2*R, 1.2*R))
    ax.set_ylim((-1.2*R, 1.2*R))
    
    # Partition phantom into 6 sections
    hpts = np.array([(-R, 0), (R, 0)])
    plt.plot(*hpts.T, color='k')
    th = -60. * (np.pi/180)
    rot = np.array([(np.cos(th), -1*np.sin(th)),
                    (np.sin(th), np.cos(th))])
    hpts = np.dot(hpts, rot)
    plt.plot(*hpts.T, color='k')
    hpts = np.dot(hpts, rot)
    plt.plot(*hpts.T, color='k')
    
    # Generate a phantom with six sections and write out to G4 macro file
    rot_angle = -(2*np.pi) / 6.                      # 6 60-degree sections
    feature_sizes = (10.0, 8.0, 6.0, 4.0, 2.0, 1.0) # Must have 6 entries
    with open('derenzo.mac', 'w') as fh:
        for i, fs in enumerate(feature_sizes):
            # Well radius
            r = fs / 2.0
            # Determine cell locations in un-rotated frame
            locs = place_wells_in_section(R, fs)
            # Rotate into the proper cell
            th = i * rot_angle
            rot_mat = np.array([(np.cos(th), -1*np.sin(th)),
                                (np.sin(th), np.cos(th))])
            locs = np.array([np.dot(l, rot_mat) for l in locs])
            # Plot wells and write to macro file
            for xy in locs:
                cyl = mpp.Circle(xy, radius=r, color="green", alpha=0.5)
                ax.add_patch(cyl)
                export_to_G4mac(fh, xy[0], xy[1], 0, r, gamma_en, num_evs, halfz=None)
    
    plt.show()
