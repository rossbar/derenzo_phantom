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

        # Define sections
        self.sections = []
        for i, (well_sep, rot_angle) in zip(self.well_seps, 
                                            np.arange(0, 360., 360. / self._num_sections)):
            section = DerenzoSection(self.radius, well_sep)
            section.apply_rotation(rot_angle)
            self.sections.append(section)

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
    def __init__(self, phantom_radius, well_separation, section_offset=0.1):
        """
        phantom_radius = radius of containing phantom
        well_separation = distance between wells (and well diameter!)
        section_offset = fraction of phantom radius that is used as a buffer
                         zone between adjacent sections.
        """
        self.R = phantom_radius
        self.well_sep = well_separation
        self.r = self.well_sep / 2.0
        self.section_offset = self.R * section_offset
        # Determine well locations
        self.place_wells_in_section()

    @property
    def row_height(self):
        return self.well_sep * np.sqrt(3)

    @property
    def num_rows(self):
        h_section = self.R - (2 * self.section_offset + self.well_sep)
        return int(np.floor(h_section / self.row_height))

    @property
    def num_wells(self):
        return np.sum(1 + np.arange(num_rows))

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
        xs, ys = [], []
        for i in range(self.num_rows):
            rn = i + 1
            for x in range(-rn, rn, 2) + 1:
                xs.append(x * self.well_sep)
                yx.append(-(self.section_offset + self.row_height * rn))
        self.locs = np.vstack((xs, ys)).T

    def apply_rotation(self, deg):
        """
        Rotate well locations around central (z) axis by 'deg' degrees.

        deg > 0: Counter-clockwise | deg < 0: clockwise
        """
        self.rot_angle = deg
        th = -1 * deg * (np.pi / 180)
        rot_mat = np.array([(np.cos(th), -np.sin(th)),
                            (np.sin(th),  np.cos(th))])
        self.locs = np.array([np.dot(l, rot_mat) for l in self.locs])
        

if __name__ == "__main__":
    radius = 50.0
    well_seps = (10.0, 8.0, 6.0, 4.0, 2.0, 1.0)
    my_phantom = DerenzoPhantom(radius, well_seps)
    my_phantom.show()
