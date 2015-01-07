#################################################################################
# pypegase v0.1                                                                 #
# By Colin Jacobs, colin@coljac.net - github.io/coljac/pypegase                 #
#                                                                               #
# Licensed under the MIT License:                                               #
#                                                                               #
# The MIT License (MIT)                                                         #
#                                                                               #
# Copyright (c) 2015 Colin Jacobs                                               #
#                                                                               #
# Permission is hereby granted, free of charge, to any person obtaining a copy  #
# of this software and associated documentation files (the "Software"), to deal #
# in the Software without restriction, including without limitation the rights  #
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     #
# copies of the Software, and to permit persons to whom the Software is         #
# furnished to do so, subject to the following conditions:                      #
#                                                                               #
# The above copyright notice and this permission notice shall be included in all#
# copies or substantial portions of the Software.                               #
#                                                                               #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, #
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE #
# SOFTWARE.                                                                     #
#                                                                               #
#################################################################################

from subprocess import Popen, PIPE, STDOUT
import os.path
import pickle
import glob
import numpy as np
from astropy.io import ascii
from astropy.table import Table, hstack
import math
import warnings as warn
from collections import OrderedDict
from os.path import expanduser

# Check: Binaries fraction/scenarios

if not os.environ.has_key("PEGASE_HOME"):
    warn.warn('Warning: PEGASE_HOME not set, using default')
 
pypeg_yorn = {True: 'y', False: 'n'}

def write_file(filename, s):
    with open(PEGASE.pegase_dir + filename, "w") as myfile:
        myfile.write(s)

def read_file(filename):
    with open (PEGASE.pegase_dir + filename, "r") as myfile:
        contents = myfile.read()
    return contents

pegdefaults = OrderedDict(
    imf = "IMF_Salpeter",
    imf_lower_mass = 0.1,
    imf_upper_mass = 120,
    sfr_type = "EXPONENTIAL_DECREASE",
    sfr_p1 = 1000,
    sfr_p2 = 1,
    )
mydefaults = pegdefaults.copy()

user_home = expanduser("~")
try:
    with open(user_home + "/.pypegase", "r") as f:
        for line in f:
            if line.strip() == "":
                continue
                values = [l.strip() for l in line.split('=')]
                mydefaults[values[0]] = values[1]
except IOError:
    pass

def get_default(key, defaultVal):
    if defaultVal is None:
        defaultVal = ""
    if key in mydefaults.keys():
        return type(defaultVal)(mydefaults[key])
    return defaultVal


class IMF(object):
    """
    Wrapper around an Initial Mass Function definition for use by PEGASE.
    When PEGASE.generate is called, this class will generate the correct i
    inputs to the PEGASE binaries, either by using pre-packaged IMFs
    such as Salpeter or Scalo86, or generating a custom file at runtime.
    """
    IMF_Kennicutt =         1
    IMF_Kroupa =            2
    IMF_MillerScalo =       3
    IMF_Salpeter =          4
    IMF_Scalo86 =           5
    IMF_Scalo98 =           6
    LOGNORMAL =             7
    RANA_AND_BASU =         8
    FERRINI =               9
    CUSTOM =                10

    def __init__(self, number=None, lower_mass=None, upper_mass=None, powers=None, gamma=None):
        """
        Models an IMF for PEGASE to use.
        :param number: The type of IMF - either a predefined one (e.g. IMF.IMF_Salpeter) or IMF.CUSTOM type.
        :param lower_mass: The lower bound for the IMF in M_sun, default 0.1
        :param upper_mass: The upper bound for the IMF in M_sun, default 120
        :param powers: For a custom IMF only. A list of tuples representing the lower mass for the power law and the power value.
        For instance, Scalo 86 would look like [(0.1, 3.2), (0.11, 2.455), (0.14, 2.0), ..., (41.69, -1.778)]
        :param gamma: A shortcut to creating a custom IMF, a power law with only one value for the entire mass range. A gamma of -1.35
        is equivalent to the default Salpeter IMF.
        :return:
        """
        self.number = number
        self.lower_mass = lower_mass
        self.upper_mass = upper_mass
        if number is None:
            self.number = eval('IMF.' + get_default('imf', 'Salpeter'))
        if lower_mass is None:
            self.lower_mass = get_default('imf_lower_mass', 0.1)
        if upper_mass is None:
            self.upper_mass = get_default('imf_upper_mass', 120.0)
            
        self.powers = powers
        if powers is None and gamma is not None:
            self.powers = [(self.lower_mass, gamma)]

    def formal_string(self):
        if self.number != IMF.CUSTOM:
            return None
        result = "%d\n" % len(self.powers)
        for power in self.powers:
            result += "%.2f     %.2f\n" % power
        result += "%.2f" % self.upper_mass
        return result

    def __str__(self):
        if self.number == IMF.CUSTOM:
            return "Custom IMF\n" + self.formal_string()
        return "IMF: %d, %.2f, %.2f" % (self.number,
                                        self.lower_mass, self.upper_mass)


class SNII_ejecta(object):
    """ Constants for use in SSPs.
    """
    MODEL_A = 'A'
    MODEL_B = 'B'
    MODEL_C = 'C'


class Extinction(object):
    """ Constants for use in Scenarios.
    """
    NO_EXTINCTION = 0
    EXTINCTION_SPHEROID = 1
    EXTINCTION_DISK_INC = 2
    EXTINCTION_DISK_SPECIFIC = 3


class Scenario(object):
    def __init__(self, binaries_fraction=0.05, metallicity_ism_0=0, infall=False, sfr=None,
                 metallicity_evolution=True, substellar_fraction=0, winds=False, neb_emission=True,
                 extinction=Extinction.NO_EXTINCTION, inclination=0):
        # self.name = name
        self.binaries_fraction = binaries_fraction
        self.metallicity = metallicity_ism_0
        self.infall = infall
        self.sfr = sfr
        self.metallicity_evolution = metallicity_evolution
        self.substellar_fraction = substellar_fraction
        self.winds = winds
        self.neb_emission = neb_emission
        self.extinction = extinction
        self.inclination = inclination

        if sfr is None:
            self.sfr = SFR(eval("SFR." + get_default("sfr_type", "EXPONENTIAL_DECREASE")),
                                p1=get_default("sfr_p1", 1e3),
                                p2=get_default("sfr_p2", 1.0)
                                )

                                
class SFR(object):
    """ Star formation rate to use in scenarios.
    Defaults is EXPONENTIAL_DECREASE, p1=1000, p2=1
    Note this isn't the PEGASE default, which is instant burst.
    """
    INSTANT_BURST = 0
    CONSTANT = 1
    EXPONENTIAL_DECREASE = 2
    GAS_POWER_LAW = 3

    def __init__(self, sfrtype=EXPONENTIAL_DECREASE, p1=None, p2=None):
        self.sfrtype = sfrtype
        self.p1 = p1
        self.p2 = p2
        if p1 is None and sfrtype == SFR.EXPONENTIAL_DECREASE:
            p1 = 1000
        if p2 is None and sfrtype == SFR.EXPONENTIAL_DECREASE:
            p2 = 1

class SSP(object):
    """ Parameters to be passed to SSPs.f
    """
    def __init__(self, imf=None, ejecta=SNII_ejecta.MODEL_B, winds=True):
        """

        :param imf: An IMF instance. If not provided, defaults to Salpeter 55
        :param ejecta: Defaults to SNII_ejecta.MODEL_B
        :param winds: On or off, defaults to True
        :return:
        """
        self.winds = winds
        self.ejecta = ejecta
        if winds is None:
            winds = get_default("ssp_winds", True)
        if ejecta is None:
            self.ejecta = eval('SNII_ejecta.' +
                               get_default('ssp_ejecta', 'MODEL_B'))
        self.imf = imf
        if imf is None:
            self.imf = IMF()


class PEGASE(object):    
    """ A wrapper around a set of parameters to PEGASE2 fortran code, and provides
    convenient access to the outputs.
    """

    pegase_dir = os.getenv('PEGASE_HOME', './PEGASE.2/')
    if pegase_dir[-1] != '/':
        pegase_dir = pegase_dir + '/'

    _datafiles = [
        '_SSPs.dat', '_tracksZ0.004.dat', '_colors.dat', '_tracksZ0.008.dat',
        '_scenarios.dat', '_tracksZ0.02.dat', '_tracksZ0.0001.dat',
        '_tracksZ0.05.dat', '_spectra.dat', '_tracksZ0.0004.dat',
        '_tracksZ0.1.dat'
    ]

    @classmethod
    def list_defaults(cls):
        for key, val in mydefaults.iteritems():
            print "%s = %s" % (key, val)
            
    @classmethod
    def from_file(cls, filename):
        """
        Unpickles a PEGASE object.
        :param filename: The pickled file.
        :return: A PEGASE instance.
        """
        return pickle.load(open(filename, "rb"))

    @classmethod
    def from_dir(cls, dirname):
        """
        Unpickles all files with a .peg extension and returns them as a list.
        :param dirname: The directory to search for .peg files
        :return: A list of resulting PEGASE objects
        """
        pegs = []
        files = glob.glob(dirname + "/*.peg")
        for file in files:
            pegs.append(PEGASE.from_file(file))
        return pegs

    def __init__(self, name, ssps=None, scenarios=None):
        self.name = name
        self.ssps = ssps
        self.generated = False
        self.scenarios = scenarios
        self._ssps_file = None
        self._scenarios_file = None
        self._spectra_file = []
        self._colors_file = []
        #self._colors_data = None
        if ssps is None:
            self.ssps = SSP()  # All defaults
        if scenarios is None:
            self.scenarios = [Scenario()]  # One scenario, all defaults

    def is_generated(self):
        """
        Returns the generated flag, i.e. the PEGASE binaries have been executed and the output exists ready to be
        queried.
        :return: True if the files exist.
        """
        return self.generated

    def generateSSPs(self):
        self.cleanup_files() # Remove existing files if present

        custom_IMF = False
        original_IMFs = None
        custom_IMF_file = None

        if self.ssps.imf.number == IMF.CUSTOM:
            # set a flag for cleanup
            custom_IMF = True
            original_IMFs = read_file("list_IMFs.dat")

            # write the IMF to a file
            custom_IMF_string = self.ssps.imf.formal_string()
            custom_IMF_file = "imfcustom.dat"
            write_file(custom_IMF_file, custom_IMF_string)

            # append that filename to the list_IMFs.dat
            write_file("list_IMFs.dat", original_IMFs.strip() + "\n" + custom_IMF_file + "\n")

        print "Generating SSPs..."
        result = self._call_binary('SSPs', self._get_ssps_string())
        # check it worked
        if result[1].rfind("error") > -1 or result[0] != 0:
            raise Exception("SSPs: " + result[1] + " -- " + str(result[0]))
        ssp_file = self.name + "_SSPs.dat"
        if os.path.isfile(PEGASE.pegase_dir + ssp_file):
            self._ssps_file = ssp_file
        else:
            raise Exception("SSPs file, " + ssp_file + " not created.")
        if not os.path.isfile(PEGASE.pegase_dir + self.name + "_tracksZ0.0004.dat"):
            raise Exception("Track file(s) missing.")

        print "Done."
        if custom_IMF:
            write_file("list_IMFs.dat", original_IMFs)
            try:
                os.remove(custom_IMF_file)
            except OSError as err:
                print "Can't remove " + custom_IMF_file
                print err
        return True

    def generateScenarios(self):
        print "\nGenerating scenarios..."
        result = self._call_binary('scenarios', self._get_scenario_string())
        if result[1].rfind("error") > -1 or result[0] != 0:
            raise Exception("Scenarios: " + result[1] + " -- " + str(result[0]))
        scenarios_file = self.name + "_scenarios.dat"
        if os.path.isfile(PEGASE.pegase_dir + scenarios_file):
            self._scenarios_file = scenarios_file
        else:
            raise Exception("Scenarios file, " + scenarios_file + " not created.")

        print "Done."

    def generateSpectra(self):
        print "Generating spectra..."
        result = self._call_binary('spectra', self._get_spectra_string())
        if result[1].rfind("error") > -1 or result[0] != 0:
            raise Exception("Scenarios: " + result[1] + " -- " +
                            str(result[0]))
        for i, scenario in enumerate(self.scenarios):
            spectra_file = self.name + "_spectra" + str(i+1) + ".dat"
            if os.path.isfile(PEGASE.pegase_dir + spectra_file):
                self._spectra_file.append(spectra_file)
            else:
                raise Exception("Spectra file, " + spectra_file +
                                " not created.")

        print "Done."

    def generateColors(self):
        print "Generating colors..."
        for i in range(len(self.scenarios)):
            result = self._call_binary('colors', self._get_colors_string(i+1))
            if result[1].rfind("error") > -1 or result[0] != 0:
                raise Exception("Colors: " + result[1] +
                                " -- " + str(result[0]))

            colors_file = self.name + "_colors" + str(i+1) + ".dat"
            if os.path.isfile(PEGASE.pegase_dir + colors_file):
                self._colors_file.append(colors_file)
            else:
                raise Exception("colors file, " + colors_file +
                                " not created.")
        print "Done."

    def generate(self):
        """
        Runs the respective PEGASE binaries (SSPs, scenarios, spectra and
        colors), checks the outputs exist and then returns. Sets the
        generated flag to True.
        :return: True if successful
        """
        self.generateSSPs()
        self.generateScenarios()
        self.generateSpectra()
        self.generateColors()
        self.generated = True
        return True

    def _get_ssps_string(self):
        """
        Generates the inputs to SSPs.f

        :return: String to be sent to stdin for that binary.
        """
        if self.ssps.imf.number < IMF.RANA_AND_BASU:
            imf = self.ssps.imf
            return "%d\n%.2f\n%4.2f\n%s\n%s\n%s\n\n" % (
                imf.number if imf.number != IMF.CUSTOM else 7, imf.lower_mass,
                imf.upper_mass, self.ssps.ejecta, pypeg_yorn[self.ssps.winds],
                self.name)
        else:
            return "todo"

    def _get_scenario_string(self):
        """
        Generates the inputs to scenarios.f

        :return: String to be sent to stdin for that binary.
        """

        result = "%s\n%s\n%.4f\n" % \
                 (self.name + "_scenarios.dat", self.name + "_SSPs.dat",
                  self.scenarios[0].binaries_fraction)
        
        for i, sc in enumerate(self.scenarios):
            result += "%s\n%.4f\n%s\n%d\n" % (self.name + "_spectra" +
                                              str(i + 1) + ".dat",
                                              sc.metallicity,
                                              pypeg_yorn[sc.infall],
                                              sc.sfr.sfrtype)

            if sc.sfr.sfrtype is not SFR.INSTANT_BURST:
                result += "%.4f\n%.4f\n" % (sc.sfr.p1, sc.sfr.p2)

                result += "%s\n%.4f\n%s\n%s\n%d\n" % (
                    pypeg_yorn[sc.metallicity_evolution],
                    sc.substellar_fraction,
                    pypeg_yorn[sc.winds], pypeg_yorn[sc.neb_emission],
                    sc.extinction)

                if sc.extinction == Extinction.EXTINCTION_DISK_SPECIFIC:
                    result += "%.4f\n" % sc.inclination

        return result + "end\n"

    def _get_spectra_string(self):
        """
        Generates the inputs to spectra.f

        :return: String to be sent to stdin for that binary.
        """
        return "%s\n\n" % (self.name + "_scenarios.dat")

    def _get_colors_string(self, scenario=1):
        """
        Generates the inputs to colors.f

        :return: String to be sent to stdin for that binary.
        """
        return "%s\n%s\n\n" % (self.name + "_spectra" + str(scenario)
                               + ".dat", self.name + "_colors" +
                               str(scenario) + ".dat")

    def _call_binary(self, name, inputs):
        # Spawns fortran binary and sends inputs via a pipe
        p = Popen([PEGASE.pegase_dir + name], stdout=PIPE, stdin=PIPE,
                  stderr=STDOUT, cwd=PEGASE.pegase_dir)
        stdout = p.communicate(input=inputs)
        # p.stdin.close()
        # p.communicate() p.wait()
        return (p.returncode, str(stdout))

    def colors(self, cols=None, scenario=1):
        """
        Returns an astropy table containing color information as output by colors.f.
        First column is the time (in Myr). Other columns are requested columns,
        or if none are provided, a table containing all columns in the file.
        :param cols: A list of column names to be returned, e.g. ['B-V', "g'-r'"]
        :return: An astropy.table.table.Table representing the output in xxx_colors.dat
        """
        if not self.generated:
            return None
        file_index = scenario - 1
        if not os.path.isfile(PEGASE.pegase_dir + self._colors_file[file_index]):
            raise Exception("Can't find the colors file (" +
                            self._colors_file[file_index]
                            + "). Regenerate?")

        # get the colors in a big table here
        with open (PEGASE.pegase_dir + self._colors_file[file_index], "r") as myfile:
            colordata_string = myfile.readlines()

        for i, line in enumerate(colordata_string):
            if line.rfind('*****')!=-1:
                break
        datastart = i + 2
        timesteps = int(colordata_string[i+1].strip())
        current = datastart
        allresults = None
        while current < len(colordata_string):
            ttable = ascii.read(colordata_string[current:(current+timesteps+1)], quotechar='^')
            current += (timesteps+1) # + header row

            if allresults is not None:
                del ttable['time']
                allresults = hstack([allresults, ttable])
            else:
                allresults = ttable

        return allresults if cols is None else allresults[['time'] + cols]

    def spectra(self, cols=None, scenario=1):
        """
        Returns an astropy table containing spectra information as output by spectra.f.
        First column is the time (in Myr). Other columns are requested columns, or if none are provided,
        a table containing all columns in the file.
        :param cols: A list of column names to be returned, e.g. ['B-V', "g'-r'"]
        :return: An astropy.table.table.Table representing the output in xxx_spectra.dat
        """
        if not self.generated:
            return None
        file_index = scenario - 1
        if not os.path.isfile(PEGASE.pegase_dir +
                              self._spectra_file[file_index]):
            raise Exception("Can't find the spectra file (" +
                            self._spectra_file[file_index] + "). Regenerate?")

        with open(PEGASE.pegase_dir + self._spectra_file[file_index], "r") as myfile:
            spectradata_string = myfile.readlines()

        for i, line in enumerate(spectradata_string):
            if line.rfind('*****') != -1:
                break
        datastart = i + 2
        timesteps, wavelengths_count, lines_count = [
            int(s) for s in spectradata_string[i+1].split()
            ]
        current = datastart

        wavelengths_end = int(datastart + math.ceil(wavelengths_count/5.0) - 1)
        lines_start = wavelengths_end + 1
        lines_end = int(wavelengths_end + math.ceil(lines_count/5.0))

        block_length = lines_end - datastart + 2
        
        columns = "time m_gal m_star m_wd m_nsbh m_substellar m_gas z_ism z_stars_mass \
        z_stars_bl l_bol od_v l_dust_l_bol sfr phot_lyman rate_snii rate_snia age_star_mass age_star_lbol".split()

        columns.extend(" ".join(spectradata_string[datastart:lines_end + 1]).split())
        indices = range(len(columns))

        if cols is not None:
            if not 'time' in cols:
                cols= ['time'] + cols

            indices = [columns.index(col) for col in cols]
            # But check there are no missing
            columns = cols
            
        allresults = Table(names=columns)
        
        current = lines_end + 1
        for i in range(timesteps):
            block = spectradata_string[current:current+block_length + 1]
            row = np.array(" ".join(block).split())[indices]
            allresults.add_row(row)
            current += block_length + 1

        return allresults

    def save_to_file(self, filename=None):
        """
        Saves the PEGASE instance to disk, so it can be unpickled later and provide convenient wrapper around
        a set of parameters and output files.
        :param filename: The file to pickle to.
        :return: None
        """
        if filename is None:
            filename = self.name + ".peg"
        pickle.dump(self, open(filename, "wb"))

    def cleanup_files(self):
        """
        Removes all the file associated with this run of PEGASE. Sets generated flag to False.
        :return:
        """
        # Let's not use a wildcard so there's no surprises
        for file in PEGASE._datafiles:
            try:
                os.remove(PEGASE.pegase_dir + self.name + file)
                print "Removed %s" % (self.name + file)
            except OSError: #  as err:
                # OK no matter
                pass
        self.generated = False

if __name__ == "__main__":
    pass

