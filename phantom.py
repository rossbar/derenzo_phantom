import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpp

from derenzo_log import export_to_G4mac

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

        # Define sections
        self.sections = []
        for well_sep, rot_angle in zip(self.well_seps, 
                                       np.arange(0, 360., 360. / self._num_sections)):
            section = DerenzoSection(self.radius, well_sep)
            section.apply_rotation(rot_angle)
            self.sections.append(section)

        # Initialize graphic (values hard-coded for now)
        self.fig = plt.figure(figsize=(6, 6))           # Aspect ratio = 1
        self.ax = self.fig.add_axes([0.1, 0.1, 0.8, 0.8])
        self.cyl_patch = mpp.Circle((0, 0), radius=self.radius, color='gray',
                                    alpha=0.3)
        self.ax.add_patch(self.cyl_patch)
        self.ax.set_xlim((-1.2*self.radius, 1.2*self.radius))
        self.ax.set_ylim((-1.2*self.radius, 1.2*self.radius))

        # Plot well locations from all sections of the phantom
        for section in self.sections:
            section.plot_wells(self.fig, self.ax)

    @property
    def area(self):
        return np.sum([s.total_area for s in self.sections])

    @property
    def num_wells(self):
        return np.sum([s.num_wells for s in self.sections])

    def show(self):
        """
        Render and display the MPL model of the phantom.
        """
        self.fig.canvas.draw()
        plt.show()

    def export_to_G4gps_macro(self, fname, num_events, gamma_en,
                              event_mode='equal_activity'):
        """
        Generate a Geant4 general particle source macro for the phantom
        instance.

        Event mode determines how many events per well. Options are:
         - 'equal_activity': Behaves as if the wells are filled with a solution
                             of constant activity/volume, thus the events per
                             well is the total number of events divided by
                             the area ratio of the cell.
         - 'equal_counts'  : Divides the total number of events by the number
                             of wells. Usefull for accentuating smaller 
                             features for testing reconstruction.
        """
        with open(fname, 'w') as fh:
            for section in self.sections:
                # Determine the number of events to run for each cell in section
                if event_mode == 'equal_activity':
                    ne = num_events * (section.well_area / self.area)
                elif event_mode == 'equal_counts':
                    ne = num_events / self.num_wells
                elif event_mode == 'subsection_area':
                    ne = num_events * (section.total_area / self.area)
                else:
                    raise ValueError("'event_mode' not understood.")
                # TODO: add & test Poisson sampling
                # Add a GPS source for each well
                for xy in section.locs:
                    export_to_G4mac(fh, xy[0], xy[1], 0, section.r, gamma_en,
                                    ne, halfz=self.depth)
                    
    def to_image(self, img_side, num_events=1., event_mode="equal_counts", plot=False):
        """
        Generates an image source from the phantom as a numpy array.
        Event mode determines how many events per well. Options equal_activity and equal_counts are as above.
        """
        if event_mode not in ["equal_activity", "equal_counts"]: raise ValueError("'event_mode' not understood.")
        r = self.radius
        px_mm = img_side/(2*r)
        shp = (img_side, img_side)
        activity = num_events/self.area
        odd_ceil = lambda x: int( 2*(np.ceil(x)//2) + 1 )

        # Generate background image
        img = np.zeros(shp)
        img_grid = np.mgrid[:shp[0], :shp[1]]

        center = np.array(shp)/2.
        outside = np.linalg.norm(img_grid \
                                - center[:, np.newaxis, np.newaxis],
                                ord=2, axis=0) > r*px_mm
        img[outside] = np.nan

        # Generate rods
        for section in self.sections:

            # Generate a template for the section
            if event_mode == 'equal_activity': activity = num_events*section.well_area/self.area
            r_px = section.r*px_mm
            d_px = odd_ceil(2*r_px)
            sl = np.arange(d_px) - d_px//2
            mg = np.array( np.meshgrid(sl, sl) )

            sl2 = sl**2
            disk = activity*(sl2 + sl2[:, np.newaxis] <= r_px)

            # Insert the template for each rod in the section
            for loc in self.mm_to_px(section.locs, shp, r, round=True):
                index_grid = tuple(mg + loc[:, np.newaxis, np.newaxis])
                img[index_grid] += disk

        if plot:
            cmap = matplotlib.cm.get_cmap("spring").copy()
            cmap.set_bad(color='white')
            plt.imshow(img, extent=[-r, r, -r, r])

        return img

    @property
    def mm_to_px(self, xy, shp, radius, round=False):
        """
        Internal function to change from mm coords to image pixels.
        """
        alt = np.array([1., -1.])
        offset = alt*radius
        scale = alt*np.array(shp)/(2*radius)
        yx_new = (scale*(xy + offset))[:, ::-1]

        return yx_new if not round else np.rint(yx_new).astype(int)
    
    
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
        # Location for section label
        self.label_xy = np.array((0, -1.1 * self.R))

    @property
    def row_height(self):
        return self.well_sep * np.sqrt(3)

    @property
    def num_rows(self):
        h_section = self.R - (2 * self.section_offset + self.well_sep)
        return int(np.floor(h_section / self.row_height))

    @property
    def num_wells(self):
        return np.sum(1 + np.arange(self.num_rows))

    @property
    def well_area(self):
        return np.pi * self.r**2

    @property
    def total_area(self):
        return self.num_wells * self.well_area

    @property
    def label(self):
        return "%.1f mm" %(self.well_sep)

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
            if self.num_rows <= 1:
                warnings.warn(("Cannot fit multiple features in section with "
                               "feature size = %s" %(self.well_sep)))
        xs, ys = [], []
        for i in range(self.num_rows):
            rn = i + 1
            for x in np.arange(-rn, rn, 2) + 1:
                xs.append(x * self.well_sep)
                ys.append(-(self.section_offset + self.row_height * rn))
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
        # Rotate well locations
        self.locs = np.array([np.dot(l, rot_mat) for l in self.locs])
        # Rotate label location
        self.label_xy = np.dot(self.label_xy, rot_mat)

    def plot_wells(self, fig, ax):
        """
        Plot the well pattern for the given section on the input figure and
        axis handles.
        """
        # Plot wells
        for xy in self.locs:
            cyl = mpp.Circle(xy, radius=self.r, color="green", alpha=0.5)
            ax.add_patch(cyl)
        # Add label
        x, y = self.label_xy
        ax.text(x, y, self.label, horizontalalignment='center',
                verticalalignment='center', fontsize=16)

if __name__ == "__main__":
    radius = 37.0
    well_seps = (8.0, 6.0, 5.0, 4.0, 3.0, 2.0)
    my_phantom = DerenzoPhantom(radius, well_seps)
    my_phantom.show()
    my_phantom.export_to_G4gps_macro('derenzo.mac', 1e6, 661.657,
                                     event_mode='equal_activity')
