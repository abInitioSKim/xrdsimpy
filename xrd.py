"""
X-ray powder diffraction pattern simulation program

Author: Sunghyun Kim(kimsunghyun@kaist.ac.kr)

Reference:
https://chemistry.osu.edu/~woodward/ch754/lect2003/xrd_peakintensities.pdf
http://kk-sinha.blogspot.kr/2014/03/calculation-of-xrd-pattern-of-h-lufeo3.html
http://www.gsc.dicp.ac.cn/down/jx/2011/2011.5.23/2-Diffraction-A.pdf
"""
import numpy as np
from numpy.linalg import norm
from numpy import dot, cross, cos, sin, arcsin
import matplotlib.pyplot as plt
import os

class XRD(object):
    """ xrd pattern generate
    """
    def __init__(self, lattice, atoms, wavelength):
        """ constructor
            Args:
                lattice: latice matrix
                atoms: list of atoms
                wavelength: float
        """
        self.lattice = lattice
        self.rec_lattice = np.linalg.inv(lattice).T * np.pi * 2
        self.atoms = atoms
        if type(wavelength) is str:
            self.wavelength = XRD.get_wavelength(wavelength)
        else:
            self.wavelength = wavelength
        self.twotheta_list = None
        self.inten_profile = None
        self.intensity_list = None

    def get_atom_form(self, element, file_name='form.txt'):
        """ get atomic form factor function parameters
            data from
            http://lamp.tu-graz.ac.at/~hadley/ss1/crystaldiffraction/
            atomicformfactors/formfactors.php
            Args:
                element: The string indicating element (eg. 'C')
                         It should be one of the values in the first column form_file
                file_file: file_name for the atomic form factor data
        """
        with open(os.path.dirname(os.path.realpath(__file__)) + '/' +
                  file_name, 'r') as form_file:
            for line in form_file.readlines():
                line = line.split()
                if line[0] == element:
                    return [float(i) for i in line[1:]]
        print 'no data for {}'.format(element)
        return None

    @staticmethod
    def get_wavelength(anode='CuKa1'):
        """ return common x-ray wavelength
        """
        xray_wavelength = {
            'CuKa1': 1.5405981,
            'CuKa2': 1.54443,
            'CuKb1': 1.39225,
            'WLa1': 1.47642,
            'WLa2': 1.48748
        }
        return xray_wavelength[anode]

    def get_kvec_list(self, k_1, k_2, k_3, max_index=None):
        """ return possible h, k, l
        """
        from itertools import product
        wavelength = self.wavelength
        kvec_list = []

        if not max_index:
            h = int(np.ceil(1. / wavelength / norm(k_1) * 4. * np.pi))
            k = int(np.ceil(1. / wavelength / norm(k_2) * 4. * np.pi))
            l = int(np.ceil(1. / wavelength / norm(k_3) * 4. * np.pi))
            max_index = max(h, k, l)

        for (h, k, l) in product(xrange(-max_index, max_index + 1), repeat=3):
            q = h * k_1 + k * k_2 + l * k_3
            sin_twotheta = wavelength * (norm(q) / 4. / np.pi)

            if -1 <= sin_twotheta <= 1:
                kvec_list.append((h, k, l))

        return np.array(kvec_list)

    def get_xrd(self):
        """ calculate xrd pattern of structure
        """
        # lattice = structure.lattice.get_matrix()
        elements = set([atom[0] for atom in self.atoms])
        a_vec, b_vec, c_vec = self.lattice

        # volume = dot(a_vec, cross(b_vec, c_vec))

        # k_1 = 2 * np.pi * cross(b_vec, c_vec) / volume
        # k_2 = 2 * np.pi * cross(c_vec, a_vec) / volume
        # k_3 = 2 * np.pi * cross(a_vec, b_vec) / volume
        k_1, k_2, k_3 = self.rec_lattice 
        inten_profile = []
        intensity_list = []
        twotheta_list = []

        for (h, k, l) in self.get_kvec_list(k_1, k_2, k_3):
            # q = h * k_1 + k * k_2 + l * k_3
            q = np.dot([h, k, l], self.rec_lattice)
            y = norm(q) / 4. / np.pi
            s_factor_list = []
            for element in elements:
            # for element in structure.get_elements():
                a_form_parm = self.get_atom_form(element)
                atom_form = a_form_parm[0] * np.exp(-a_form_parm[1] * y ** 2) +\
                            a_form_parm[2] * np.exp(-a_form_parm[3] * y ** 2) +\
                            a_form_parm[4] * np.exp(-a_form_parm[5] * y ** 2) +\
                            a_form_parm[6] * np.exp(-a_form_parm[7] * y ** 2) +\
                            a_form_parm[8]

                for atom in [atom for atom in self.atoms
                             if atom[0] == element]:
                    pos = atom[1]
                    plane = np.array([h, k, l])
                    phase = np.dot(np.dot(pos, self.lattice),
                                   np.dot(self.rec_lattice, plane))
                    # print phase, dot(pos, plane)
                    s_factor = atom_form *\
                               np.exp(-2 * np.pi * phase * 1j)
                    s_factor_list.append(s_factor)

            if -1 <= self.wavelength * y <= 1:
                theta = arcsin(self.wavelength * y)
                if theta == np.pi / 2. or theta == 0:
                    continue
                lp_factor = (1. + cos(theta * 2) ** 2) /\
                            (8 * sin(theta) ** 2 * cos(theta))
                inten = norm(sum(s_factor_list)) ** 2 * lp_factor
                intensity_list.append(inten)

                twotheta = 2 * 180 / np.pi * theta
                twotheta_list.append(twotheta)
                inten_profile.append((h, k, l, twotheta, inten))
            else:
                print self.wavelength * y

        self.inten_profile = inten_profile
        self.twotheta_list = twotheta_list
        self.intensity_list = intensity_list

    def plot(self):
        """ Plot intensity vs two theta
        """
        def sum_intensity(x_list, y_list, tolerance=1E-5):
            """ Multiplicity treatment
                similar to histogram
                if same x in x_list, merge them
            """
            new_x_list = []
            new_y_list = []
            for i, x_value in enumerate(x_list):
                index = np.where(abs(new_x_list - x_value) <= tolerance)[0]
                if len(index) > 0:
                    new_y_list[index] += y_list[i]
                else:
                    new_x_list.append(x_list[i])
                    new_y_list.append(y_list[i])
            return new_x_list, new_y_list

        def normalize(y_list, max_value=100.):
            """ maximum y set to be 100
            """
            y_list = np.array(y_list)
            y_list = y_list / max(y_list) * max_value
            return y_list

        if not self.twotheta_list:
            self.get_xrd()
        x_list, y_list = self.twotheta_list, self.intensity_list
        x_list, y_list = sum_intensity(x_list, y_list)
        y_list = normalize(y_list)
        # print X
        # print Y
        for i in range(len(x_list)):
            plt.plot([x_list[i], x_list[i]], [0, y_list[i]], '-', c='#5ab4ac')
        plt.xlim((0, 180))
        plt.xlabel('2$\\theta$')
        plt.ylabel("Intensity")
        plt.show()


if __name__ == '__main__':
    """ example for fcc Si
    """
    from vasp_io import readCONTCAR
    # lattice = [5.43, 5.43, 5.43,
    #            90 * np.pi/180, 90 * np.pi/180, 90 * np.pi/180]
    # lat_const = 5.465
    # lattice_vecs = np.array([[0.5, 0.5, 0],
    #                          [0, 0.5, 0.5],
    #                          [0.5, 0, 0.5]])
    # atoms = [['Siv', np.array([0, 0, 0])],
    #          ['Siv', np.array([0.25, 0.25, 0.25])]]

    # lat_const, lattice_vecs, atoms = readCONTCAR('POSCAR_Si')
    # lat_const, lattice_vecs, atoms = readCONTCAR('POSCAR_MoTe2')
    lat_const, lattice_vecs, atoms = readCONTCAR('../CONTCAR_03_1')

    lattice = np.array(lattice_vecs) * lat_const
    # wavelength = XRD.get_wavelength('CuKa1')
    xrd_si = XRD(lattice, atoms, 'CuKa1')
    xrd_si.get_xrd()
    xrd_si.plot()

