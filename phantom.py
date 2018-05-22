import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpp

class DerenzoPhantom(object):
    """
    Describes a cylindrical Derenzo phantom.
    """
    # Always subdivided into six sections, 60 degrees per section
    _num_sections = 6
    def __init__(self, radius, well_separations, cyl_height=0, unit="mm"):
        """
        Create a Derenzo phantom based on the input parameters.

        A circle with radius is partitioned into six equal sections.
        Each section is then populated with circular 'wells' that have a 
        diameter and spacing defined by well_separations.

        Note the kwargs 'cyl_height' and 'unit' are only relevant if phantom is
        used to generate a Geant4 general point source macro.

        Parameters
        ----------
        radius : float
            The radius of the phantom

        well_separations : array_like, 1D, with len == 6
            A sequence specifying the well diameters and separations in each
            of the six sections. 
            Sequence must contain 6 numeric elements.

        cyl_height : float, default: 0
            Depth of the wells in the phantom. All wells have the same depth.
            If 0, then the phantom is 2D.

        unit : string
            Unit of length for all relevant measurements. Default is 'mm'.
            Can be any valid unit for length from G4SystemOfUnits.hh

        Examples
        --------
        >>> radius = 50.0
        >>> well_seps = (10.0, 8.0, 6.0, 4.0, 2.0, 1.0)
        >>> my_phantom = DerenzoPhantom(radius, well_seps)
        """
        self.radius = radius
        self.well_seps = well_separations
        self.depth = cyl_height
        self.unit = unit

        # Initialize graphic (values hard-coded for now)
        self.fig = plt.figure(figsize=(6, 6))           # Aspect ratio = 1
        self.ax = self.fig.add_axes([0.1, 0.1, 0.8, 0.8])
        self.cyl_patch = mpp.Circle((0, 0), radius=self.radius, color='gray',
                                    alpha=0.3)
        self.ax.add_patch(self.cyl_patch)
        self.ax.set_xlim((-1.2*self.radius, 1.2*self.radius))
        self.ax.set_ylim((-1.2*self.radius, 1.2*self.radius))

    def show(self):
        """
        Render and display the MPL model of the phantom.
        """
        self.fig.canvas.draw()
        plt.show()

class DerenzoSection(object):
    """
    One of six sub-sections of the Derenzo phantom.
    """
    def __init__(self, phantom_radius, well_separation, rotation_deg,
                 section_offset=0.1):
        self.R = phantom_radius
        self.well_sep = well_separation
        self.r = self.well_sep / 2.0
        self.th = rotation_deg
        self.thr = self.th * (np.pi/180)
        self.rot_mat = np.array([(np.cos(self.thr), -1*np.sin(self.thr)),
                                 (np.sin(self.thr), np.cos(self.thr))])
        self.section_offset = self.R * section_offset

    @property
    def row_height(self):
        return self.well_sep * np.sqrt(3)

    @property
    def num_rows(self):
        h_section = self.R - (2 * self.section_offset + self.well_sep)
        return int(np.floor(h_section / self.row_height))

    def place_wells_in_section(self):
        """
        Method analogous to derenzo_log.place_wells_in_section
        """
        # If only one (or no) wells fit in section, try reducing section offset
        # in attempt to squeeze more in.
        #NOTE: Could add property called 'placement_policy' or something like
        # that which governs this behavior
        if self.num_rows <= 1:
            self.section_offset = 0.0
            if num_rows <= 1:
                warnings.warn(("Cannot fit multiple features in section with "
                               "feature size = %s" %(self.well_sep)))

        

if __name__ == "__main__":
    radius = 50.0
    well_seps = (10.0, 8.0, 6.0, 4.0, 2.0, 1.0)
    my_phantom = DerenzoPhantom(radius, well_seps)
    my_phantom.show()
