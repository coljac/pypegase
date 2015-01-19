#######################################################################
# pypegase v0.1                                                       #
# By Colin Jacobs, colin@coljac.net - github.io/coljac/pypegase       #
#                                                                     #
# Copyright (c) 2015 Colin Jacobs                                     #
#                                                                     #
# Licensed under the MIT License: see license.txt                     #
#                                                                     #
#######################################################################

from subprocess import Popen, PIPE, STDOUT
import os.path
import pickle
import glob
import sys
import numpy as np
from astropy.io import ascii
from astropy.table import Table, hstack, Column
from time import sleep
import threading
import math
import warnings as warn
from collections import OrderedDict
from os.path import expanduser

# Check: Binaries fraction/scenarios

if not ("PEGASE_HOME") in os.environ.keys():
    warn.warn('Warning: PEGASE_HOME not set, using default')

pypeg_yorn = {True: 'y', False: 'n'}

pegdefaults = OrderedDict(
    imf="IMF_Salpeter",
    imf_lower_mass=0.1,
    imf_upper_mass=120,
    sfr_type="EXPONENTIAL_DECREASE",
    sfr_p1=1000,
    sfr_p2=1,
    scenario_binaries_fraction=0.05,
    scenario_metallicity_ism=0,
    scenario_infall=False,
    scenario_sfr=None,
    scenario_metallicity_evolution=True,
    scenario_substellar_fraction=0,
    scenario_winds=False,
    scenario_neb_emission=True,
    scenario_extinction="NO_EXTINCTION",
    scenario_inclination=0
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
    When PEGASE.generate is called, this class will generate the corect i
    inputs to the PEGASE binaries, either by using pre-packaged IMFs
    such as Salpeter or Scalo86, or generating a custom file at runtime.
    """
    IMF_Kennicutt = 1
    IMF_Kroupa = 2
    IMF_MillerScalo = 3
    IMF_Salpeter = 4
    IMF_Scalo86 = 5
    IMF_Scalo98 = 6
    LOGNORMAL = 7
    RANA_AND_BASU = 8
    FERRINI = 9
    CUSTOM = 10

    def __init__(self, number=None, lower_mass=None, upper_mass=None,
                 powers=None, gamma=None):
        """
        Models an IMF for PEGASE to use.
        :param number: The type of IMF - either a predefined one (e.g.
        IMF.IMF_Salpeter) or IMF.CUSTOM type.
        :param lower_mass: The lower bound for the IMF in M_sun, default 0.1
        :param upper_mass: The upper bound for the IMF in M_sun, default 120
        :param powers: For a custom IMF only. A list of tuples representing the
        lower mass for the power law and the power value.
        For instance, Scalo 86 would look like [(0.1, 3.2), (0.11, 2.455),
        (0.14, 2.0), ..., (41.69, -1.778)]
        :param gamma: A shortcut to creating a custom IMF, a power law with
        only one value for the entire mass range. A gamma of -1.35 is
        equivalent to the default Salpeter IMF.
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
        return self.__repr__()
        # if self.number == IMF.CUSTOM:
        #     return "Custom IMF\n" + self.formal_string()
        # return "IMF: %d, %.2f, %.2f" % (self.number,
        #                                 self.lower_mass, self.upper_mass)

    def __repr__(self):
        return ("%s [number=%d%s, lower_mass=%.4f, upper_mass=%.4f, " +
                "powers=%s]") % (
                    ("Custom IMF" if self.number == IMF.CUSTOM else "IMF"),
                    self.number,
                    [k for k in IMF.__dict__.keys() if IMF.__dict__[k]
                     == self.number],
                    self.lower_mass, self.upper_mass, self.powers)


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
    def __init__(self, binaries_fraction=0.05, metallicity_ism_0=0,
                 infall=False, infall_timescale=1e3, infall_gas_metallicity=0,
                 sfr=None,
                 metallicity_evolution=True, stellar_metallicity=0.02,
                 substellar_fraction=0,
                 galactic_winds=False, age_of_winds=2.001e4,
                 neb_emission=True, 
                 extinction=Extinction.NO_EXTINCTION, inclination=0.):

        self.binaries_fraction = binaries_fraction
        self.metallicity_ism_0 = metallicity_ism_0
        self.infall = infall
        self.infall_timescale = infall_timescale
        self.infall_gas_metallicity = infall_gas_metallicity
        self.sfr = sfr
        self.metallicity_evolution = metallicity_evolution
        self.stellar_metallicity = stellar_metallicity
        self.substellar_fraction = substellar_fraction
        self.galactic_winds = galactic_winds
        self.age_of_winds = age_of_winds
        self.neb_emission = neb_emission
        self.extinction = extinction
        self.inclination = inclination

        if sfr is None:
            self.sfr = SFR(eval("SFR." + get_default("sfr_type",
                                                     "EXPONENTIAL_DECREASE")),
                           p1=get_default("sfr_p1", 1e3),
                           p2=get_default("sfr_p2", 1.0))

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return (
            "Scenario [\n    binaries_fraction=%.2f, metallicity_ism_0=%.4f, infall=%s," +
            " infall_timescale=%.4f,\n    infall_gas_metallicity=%.4f,\n" +
            "    %s,\n    metallicity_evolution=%s, stellar_metallicity=%.4f,\n" +
            "    substellar_fraction=%.2f, galactic_winds=%s, age_of_winds=%.4f,\n" +
            "    neb_emission=%s, extinction=%d, inclination=%.2f\n]") % (
                self.binaries_fraction, self.metallicity_ism_0,
                pypeg_yorn[self.infall], self.infall_timescale,
                self.infall_gas_metallicity, self.sfr,
                pypeg_yorn[self.metallicity_evolution],
                self.stellar_metallicity,
                self.substellar_fraction, pypeg_yorn[self.galactic_winds],
                self.age_of_winds,  pypeg_yorn[self.neb_emission],
                self.extinction, self.inclination)

class SFR(object):
    """Star formation rate/history to use in scenarios.
    Defaults is EXPONENTIAL_DECREASE, p1=1000, p2=1
    Note this isn't the PEGASE default, which is instant burst.
    """
    FILE_SFR_AND_Z = -2
    FILE_SFR = -1
    INSTANT_BURST = 0
    CONSTANT = 1
    EXPONENTIAL_DECREASE = 2
    GAS_POWER_LAW = 3

    def __init__(self, sfrtype=EXPONENTIAL_DECREASE, p1=None, p2=None,
                 filename=None):
        self.sfrtype = sfrtype
        self.p1 = p1
        self.p2 = p2
        self.filename = filename
        if p1 is None and sfrtype == SFR.EXPONENTIAL_DECREASE:
            p1 = 1000
        if p2 is None and sfrtype == SFR.EXPONENTIAL_DECREASE:
            p2 = 1
        if p1 is None and sfrtype == SFR.CONSTANT:
            p1 = 5e-5
        if p2 is None and sfrtype == SFR.CONSTANT:
            p2 = 0.20001E05
        if p1 is None and sfrtype == SFR.GAS_POWER_LAW:
            p1 = 1
        if p2 is None and sfrtype == SFR.GAS_POWER_LAW:
            p2 = 3e3

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "SFR: [type=%d%s, p1=%.4f, p2=%.4f, filename=%s]" \
            % (self.sfrtype, [k for k in SFR.__dict__.keys() if
                              SFR.__dict__[k] == self.sfrtype],
               self.p1, self.p2, self.filename)


class SSP(object):
    """ Parameters to be passed to SSPs.f
    """
    def __init__(self, imf=None, ejecta=SNII_ejecta.MODEL_B,
                 galactic_winds=True):
        """

        :param imf: An IMF instance. If not provided, defaults to Salpeter 55
        :param ejecta: Defaults to SNII_ejecta.MODEL_B
        :param galactic_winds: On or off, defaults to True
        :return:
        """
        self.galactic_winds = galactic_winds
        self.ejecta = ejecta
        if galactic_winds is None:
            galactic_winds = get_default("ssp_winds", True)
        if ejecta is None:
            self.ejecta = eval('SNII_ejecta.' +
                               get_default('ssp_ejecta', 'MODEL_B'))
        self.imf = imf
        if imf is None:
            self.imf = IMF()

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "SSP: [%s, ejecta=%s, winds=%s]" % (
            self.imf, self.ejecta, pypeg_yorn[self.galactic_winds])


class Filter(object):
    """Represents a filter as used by PEGASE and this code to
generate colors."""
    TRANSMISSION_ENERGY = 0
    TRANSMISSION_NUMBER_PHOTONS = 1
    CALIBRATION_VEGA = 1
    CALIBRATION_AB = 2
    CALIBRATION_THUAN_GUNN = 3

    @classmethod
    def get_by_name(cls, name):
        filters = Filter.get_all()
        for filter in filters:
            if filter.name == name:
                return filter
        return None

    @classmethod
    def get_all(cls):
        filters_string = PEGASE._read_file("filters.dat")
        lines = filters_string.split("\n")
        num_filters = int(lines[0].strip())
        filters = []
        line_num = 1
        for i in range(num_filters):
            fdata = lines[line_num].split()
            txs = []
            for j in range(int(fdata[0])):
                a, b = lines[line_num + j + 1].split()
                txs.append((float(a), float(b)))
            line_num += (int(fdata[0]) + 1)
            f = Filter(fdata[3].replace("'", ""), txs,
                       transmission_type=int(fdata[1]),
                       calibration_type=int(fdata[2]))
            filters.append(f)
        return filters

    def __init__(self, name, transmissions,
                 transmission_type=TRANSMISSION_NUMBER_PHOTONS,
                 calibration_type=CALIBRATION_AB):
        """
        :param name: A name for the filter to be used in outputs
        :param responses: A list of tuples, or a 2D numpy aray, containing
        wavelength-response
        pairs.
        """
        self.name = name
        self.transmissions = transmissions
        self.transmission_type = transmission_type
        self.calibration_type = calibration_type

    def get_transmission_at_wavelength(self, wavelength):
        """Returns the transmission coefficient at wavelength, using linear
        interpolation.
        :param wavelength: Wavelength in Angstroms
        :return: Transmission coefficient (0 < t < 1)"""

        for i, trans in enumerate(self.transmissions):
            if trans[0] <= wavelength and i < len(self.transmissions) - 1 \
               and self.transmissions[i+1][0] > wavelength:
                # Linear interpolation
                t = self.transmissions[i]
                t1 = self.transmissions[i+1]
                t_coeff = ((wavelength - t[0])/(t1[0] - t[0]))*(t1[1]-t[1]) + t[1]
                return t_coeff
        return 0.

    def get_lower_bound(self):
        return self.transmissions[0][0]

    def get_upper_bound(self):
        return self.transmissions[-1][0]

    def get_midpoint(self):
        tot = 0
        for tx in self.transmissions:
            tot += tx[0]
        return tot/float(len(self.transmissions))

    def get_calibration(self):
        """Get filter calibration info. Example: filter.get_calibration()['mAV(Vega)']
        :return: A table row for the filter, with columns matching those in
        calib.dat
        """
        calib = PEGASE.get_calib()
        for f in calib:
            if f['filter'] == self.name:
                return f
        return None

    def to_string(self):
        """Returns a string representation that matches PEGASE's filter.dat
        :return: String to add to filter.dat"""
        result = "%d    %d      %d        '%s' (pypegase)\n" % (
            len(self.transmissions),
            self.transmission_type,
            self.calibration_type,
            self.name)
        for trans in self.transmissions:
            result += "%.0f    %.3f\n" % (trans[0], trans[1])
        return result

    def save_filter(self):
        """Saves the filter to filters.dat and runs calib.f. Doesn't mess with 
        colors.f and Fortran formats, though."""
        filters = PEGASE._read_file("filters.dat")
        filters_lines = filters.split("\n")
        PEGASE._write_file("filters.dat.bak", filters)
        filter_number = (int(filters_lines[0]) + 1)
        filters_lines[0] = "  %d" % filter_number
        filters = "\n".join(filters_lines)
        filters += self.to_string()
        PEGASE._write_file("filters.dat", filters)
        PEGASE._call_binary("calib", "")
        print "Filter saved with number %d." % filter_number

    def __str__(self):
        return self.__repr__()
        
    def __repr__(self):
        return "Filter: [name=%s, transmission_type=%d, calibration_type=%d, " +\
            "transmission=%s...%s]" % (
                self.name, self.transmission_type, self.calibration_type,
                self.transmissions[0], self.transmissions[-1])


class PEGASE(object):
    """ A wrapper around a set of parameters to PEGASE2 fortran code, and provides
    convenient access to the outputs.
    """

    pegase_dir = os.getenv('PEGASE_HOME', './PEGASE.2/')
    if pegase_dir[-1] != '/':
        pegase_dir = pegase_dir + '/'

    _datafiles = [
        '_SSPs.dat', '_tracksZ0.004.dat', '_tracksZ0.008.dat',
        '_scenarios.dat', '_tracksZ0.02.dat', '_tracksZ0.0001.dat',
        '_tracksZ0.05.dat', '_tracksZ0.0004.dat',
        '_tracksZ0.1.dat'  # plus colors_n_ and spectra_n_
    ]

    @classmethod
    def _write_file(cls, filename, s):
        with open(PEGASE.pegase_dir + filename, "w") as myfile:
            myfile.write(s)

    @classmethod
    def _read_file(cls, filename):
        with open(PEGASE.pegase_dir + filename, "r") as myfile:
            contents = myfile.read()
        return contents

    @classmethod
    def list_defaults(cls):
        for key, val in pegdefaults.iteritems():
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

    @classmethod
    def get_calib(cls):
        calib_data = PEGASE._read_file("calib.dat")
        calib_data = calib_data.replace("mean lambda", "mean_lambda")
        calib_data = calib_data.replace("lambda eff(Vega)", "lambda_eff(Vega)")
        table = ascii.read(calib_data, names=[
            "filter", "index", "Flambda(Vega)", "area", "mean lambda",
            "lambda eff(Vega)", "mAB(Vega)", "mTG(Vega)", "Flambda(Sun)"
        ])
        return table

    def __init__(self, name, ssps=None, scenarios=None):
        """Represents PEGASE model.
        :param name: The name which will be used as a prefix for all PEGASE 
        files.
        :param ssps: An SSP instance, wrapping the IMF and SFR. Can be omitted
        to use defaults.
        :param scenarios: A list of Scenario objects. Defaults to one scenario 
        with default parameters.
        """
        self.name = name
        self.ssps = ssps
        self.generated = False
        self.scenarios = scenarios
        self._ssps_file = None
        self._scenarios_file = None
        self._spectra_file = []
        self._colors_file = []
        # self._colors_data = None
        if ssps is None:
            self.ssps = SSP()  # All defaults
        if scenarios is None:
            self.scenarios = [Scenario()]  # One scenario, all defaults

    def is_generated(self):
        """
        Returns the generated flag, i.e. the PEGASE binaries have been 
        executed and the output exists ready to be queried.
        :return: True if the files exist.
        """
        return self.generated

    def generateSSPs(self):
        """Generates the SSP tracks, and checks to makes sure everything worked.
        There is not much reason to prefer this to generate() which will do the 
        lot.
        """
        self.cleanup_files(silent=True)  # Remove existing files if present

        custom_IMF = False
        original_IMFs = None
        custom_IMF_file = None

        if self.ssps.imf.number == IMF.CUSTOM:
            # set a flag for cleanup
            custom_IMF = True
            original_IMFs = PEGASE._read_file("list_IMFs.dat")

            # write the IMF to a file
            custom_IMF_string = self.ssps.imf.formal_string()
            custom_IMF_file = "imfcustom.dat"
            PEGASE._write_file(custom_IMF_file, custom_IMF_string)

            # append that filename to the list_IMFs.dat
            PEGASE._write_file("list_IMFs.dat", original_IMFs.strip() +
                               "\n" + custom_IMF_file + "\n")

        self._progress(0, "Generating SSPs")

        try:
            result = PEGASE._call_binary('SSPs', self._get_ssps_string())
        except Exception as e:
            print e
            if custom_IMF:
                PEGASE._write_file("list_IMFs.dat", original_IMFs)
                try:
                    os.remove(custom_IMF_file)
                except OSError as er:
                    print "Error cleanup: can't remove " + custom_IMF_file
                    print er
            raise e

        # check it worked
        if result[1].rfind("eror") > -1 or result[0] != 0:
            raise Exception("SSPs: " + result[1] + " -- " + str(result[0]))
        ssp_file = self.name + "_SSPs.dat"
        if os.path.isfile(PEGASE.pegase_dir + ssp_file):
            self._ssps_file = ssp_file
        else:
            raise Exception("SSPs file, " + ssp_file + " not created.")
        if not os.path.isfile(PEGASE.pegase_dir + self.name +
                              "_tracksZ0.0004.dat"):
            raise Exception("Track file(s) missing.")

        if custom_IMF:
            PEGASE._write_file("list_IMFs.dat", original_IMFs)
            try:
                os.remove(custom_IMF_file)
            except OSError as er:
                print "Cleanup: can't remove " + custom_IMF_file
                print er
        return True

    def generateScenarios(self):
        """Generates the ..._scenarios.dat file passed to spectra and colors.
        Call this method if you have changed any scenarios since calling
        generate() initially.
        """
        self._delete_file("_scenarios.dat")
        result = PEGASE._call_binary('scenarios', self._get_scenario_string())
        if result[1].rfind("eror") > -1 or result[0] != 0:
            raise Exception("Scenarios: " + result[1] + " -- " +
                            str(result[0]))
        scenarios_file = self.name + "_scenarios.dat"
        if os.path.isfile(PEGASE.pegase_dir + scenarios_file):
            self._scenarios_file = scenarios_file
        else:
            raise Exception("Scenarios file, " + scenarios_file +
                            " not created.")

    def generateSpectra(self):
        """Generates spectra[n].dat files, one per scenario and in the same
        order as the scenarios list. Call this manually if scenarios were
        changed since the model was generated initially."""
        for i in range(len(self._spectra_file)):
            self._delete_file("_spectra%d.dat") % i
        self._spectra_file = []

        result = PEGASE._call_binary('spectra', self._get_spectra_string())
        if result[1].rfind("eror") > -1 or result[0] != 0:
            raise Exception("Scenarios: " + result[1] + " -- " +
                            str(result[0]))
        for i, scenario in enumerate(self.scenarios):
            spectra_file = self.name + "_spectra" + str(i+1) + ".dat"
            if os.path.isfile(PEGASE.pegase_dir + spectra_file):
                self._spectra_file.append(spectra_file)
            else:
                raise Exception("Spectra file, " + spectra_file +
                                " not created.")

    def generateColors(self):
        """Generates colors[n].dat files, one per scenario and in the same
        order as the scenarios list. Call this manually if scenarios were
        changed since the model was generated initially."""
        for i in range(len(self._colors_file)):
            self._delete_file("_colors%d.dat") % i
        self._colors_file = []

        for i in range(len(self.scenarios)):
            result = PEGASE._call_binary('colors',
                                         self._get_colors_string(i+1))
            if result[1].rfind("eror") > -1 or result[0] != 0:
                raise Exception("Colors: " + result[1] +
                                " -- " + str(result[0]))

            colors_file = self.name + "_colors" + str(i+1) + ".dat"
            if os.path.isfile(PEGASE.pegase_dir + colors_file):
                self._colors_file.append(colors_file)
            else:
                raise Exception("colors file, " + colors_file +
                                " not created.")

    def _progress(self, progress, message=""):
        sys.stdout.write('\r{0} [{1}] {2}%   {3}'.format(
            "Generating: ", ('#'*int(progress/5)).ljust(20, " "),
            int(progress), message))
        sys.stdout.flush()

    def _generate(self, flagobj, *args):
        try:
            flagobj['message'] = "(Creating SSP tracks)"
            self.generateSSPs()
            flagobj['message'] = "(Creating scenarios)"
            self.generateScenarios()
            flagobj['message'] = "(Creating spectra)"
            self.generateSpectra()
            flagobj['message'] = "(Creating colors)"
            self.generateColors()
            self.generated = True
        except Exception as ex:
            flagobj['finished'] = True
            flagobj['message'] = 'Error'
            raise ex
        
        flagobj['finished'] = True

    def generate(self):
        """
        Runs the respective PEGASE binaries (SSPs, scenarios, spectra and
        colors), checks the outputs exist and then returns. Sets the
        generated flag to True. 
        :return: True if successful
        """

        flagobj = dict()
        flagobj['finished'] = False
        flagobj['message'] = "Starting"
        
        generating_thread = threading.Thread(target=self._generate,
                                             args=(flagobj,))
        generating_thread.start()
        
        total_files = 9 + len(self.scenarios)*2
        while not flagobj['finished']:
            filecount = len([name for name in os.listdir(PEGASE.pegase_dir) if os.path.isfile(
                os.path.join(PEGASE.pegase_dir, name)) and name.startswith(self.name + "_")])
            self._progress(100.0 * filecount/total_files, flagobj['message'])
            if flagobj['message'] == 'Error':
                break
            sleep(0.1)
        print ""
        return True

    def _get_ssps_string(self):
        """
        Generates the inputs to SSPs.f

        :return: String to be sent to stdin for that binary.
        """
        imf = self.ssps.imf
        return "%d\n%.2f\n%4.2f\n%s\n%s\n%s\n\n" % (
            imf.number if imf.number != IMF.CUSTOM else 7, imf.lower_mass,
            imf.upper_mass, self.ssps.ejecta, pypeg_yorn[self.ssps.galactic_winds],
            self.name)

    def _get_scenario_string(self):
        """
        Generates the inputs to scenarios.f

        :return: String to be sent to stdin for that binary.
        """

        result = "%s\n%s\n%.4f\n" % \
                 (self.name + "_scenarios.dat", self.name + "_SSPs.dat",
                  self.scenarios[0].binaries_fraction)

        for i, sc in enumerate(self.scenarios):
            result += "%s\n%.4f\n%s\n" % (self.name + "_spectra" +
                                          str(i + 1) + ".dat",
                                          sc.metallicity_ism_0,
                                          pypeg_yorn[sc.infall])
            if sc.infall:
                result += "%.4f\n%.4f\n" % (
                    sc.infall_timescale, sc.infall_gas_metallicity)

            result += "%d\n" % sc.sfr.sfrtype

            if sc.sfr.sfrtype in [
                    SFR.CONSTANT, SFR.EXPONENTIAL_DECREASE, SFR.GAS_POWER_LAW]:
                result += "%.4f\n%.4f\n" % (sc.sfr.p1, sc.sfr.p2)

            if sc.sfr.sfrtype in [SFR.FILE_SFR, SFR.FILE_SFR_AND_Z]:
                result += "%s\n" % (sc.filename)

            if sc.sfr.sfrtype != SFR.FILE_SFR_AND_Z:
                result += "%s\n" % pypeg_yorn[sc.metallicity_evolution]
                if not sc.metallicity_evolution:
                    result += "%.4f\n" % sc.stellar_metallicity

            result += "%.4f\n%s\n" % (sc.substellar_fraction,
                                      pypeg_yorn[sc.galactic_winds])
            if sc.galactic_winds:
                result += "%.4f\n" % sc.age_of_winds

            result += "%s\n%d\n" % (pypeg_yorn[sc.neb_emission], sc.extinction)

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
    @classmethod
    def _call_binary(cls, name, inputs):
        # Spawns fortran binary and sends inputs via a pipe
        p = Popen([PEGASE.pegase_dir + name], stdout=PIPE, stdin=PIPE,
                  stderr=STDOUT, cwd=PEGASE.pegase_dir)
        stdout = p.communicate(input=inputs)
        # p.stdin.close()
        # p.communicate() p.wait()
        return (p.returncode, str(stdout))

    def colors(self, cols=None, scenario=1, time_lower=None, time_upper=None):
        """
        Returns an astropy table containing color information as output by
        colors.f.
        First column is the time (in Myr). Other columns are requested columns,
        or if none are provided, a table containing all columns in the file.
        :param cols: A list of column names to be returned, e.g. ['B-V',
        "g'-r'"]
        :param time_lower: If specified, filters out rows with time below
        this value in Myr
        :param time_upper: If specified, filters out rows with the time above
        this value
        :return: An astropy.table.table.Table representing the output in
        xxx_colors.dat
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
        curent = datastart
        allresults = None
        while curent < len(colordata_string):
            ttable = ascii.read(colordata_string[curent:(curent+timesteps+1)],
                                quotechar='^')
            curent += (timesteps+1)  # + header row

            if allresults is not None:
                del ttable['time']
                allresults = hstack([allresults, ttable])
            else:
                allresults = ttable

        if time_lower is not None or time_upper is not None:
            time_lower = time_lower if time_lower is not None else 0
            time_upper = time_upper if time_upper is not None else 1e11

            for rownum in range(len(allresults)-1, -1, -1):
                if allresults['time'][rownum] < time_lower or \
                   allresults['time'][rownum] > time_upper:
                    allresults.remove_row(rownum)

        return allresults if cols is None else allresults[['time'] + cols]

    def _spectra(self, scenario=1):
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
        columns = "time m_gal m_star m_wd m_nsbh m_substellar m_gas z_ism z_stars_mass \
        z_stars_bl l_bol od_v l_dust_l_bol sfr phot_lyman rate_snii rate_snia \
        age_star_mass age_star_lbol".split()
        return timesteps, wavelengths_count, lines_count, spectradata_string, \
            datastart, columns

    def spectra(self, cols=None, scenario=1, time_lower=0, time_upper=1e10):
        """
        Returns an astropy table containing spectra information as output by
        spectra.f.
        First column is the time (in Myr). Other columns are requested columns,
        or if none are provided, a table containing all columns in the file.
        :param cols: A list of column names to be returned, e.g. ['B-V',
        "g'-r'"]
        :param time_lower: If specified, filters out rows with time below this
        value in Myr
        :param time_upper: If specified, filters out rows with the time above
        this value
        :param scenario: The number of the applicable scenario, defaults to
        the first.
        :return: An astropy.table.table.Table representing the output in
        xxx_spectra.dat
        """
        if not self.generated:
            return None
        timesteps, wavelengths_count, lines_count, spectradata_string, datastart, columns \
            = self._spectra(scenario=scenario)
        
        curent = datastart

        wavelengths_end = int(datastart + math.ceil(wavelengths_count/5.0) - 1)
        lines_start = wavelengths_end + 1
        lines_end = int(wavelengths_end + math.ceil(lines_count/5.0))

        block_length = lines_end - datastart + 2

        columns.extend(" ".join(spectradata_string[datastart:lines_end + 1]).split())
        indices = range(len(columns))

        if cols is not None:
            if not 'time' in cols:
                cols= ['time'] + cols

            indices = [columns.index(col) for col in cols]
            # But check there are no missing
            columns = cols

        allresults = Table(names=columns)

        curent = lines_end + 1
        for i in range(timesteps):
            block = spectradata_string[curent:curent+block_length + 1]
            row = np.array(" ".join(block).split())[indices]

            if (time_lower is None or int(row[0]) >= time_lower) and \
               (time_upper is None or int(row[0]) <= time_upper):
                allresults.add_row(row)
            curent += block_length + 1

        return allresults

    # TODO: Speedup
    def get_spectrum(self, time=None, wl_lower=0, wl_upper=1600000, lines=False,
                     continuum=True, scenario=1, redshift=0.0):  # , units='flux density'):
        """At timestep t=time, returns spectral data between the specified wavelengths.
        :param time: the timestep at which to return the spectum; if exact timestep does 
        not exist, raises ValueError
        :param wl_lower: Lower wavelength to include
        :param wl_upper: Upper wavelength to include
        :param lines: Whether to include spectral lines (default False)
        :param continuum: Whether to include continuum values (default True)
        :param units: 'flux density' gives erg/s/Angstrom; 'flux' gives ergs/s, 
        i.e. flux density * lambda 
        :return: An astropy table containing columns 'wavelength', 'flux'
        """
        if time is None:
            raise ValueError("No time specified")
        if not (lines or continuum):
            raise ValueError("Must specify either lines, continuum or both.")
        if lines:
            raise ValueError("Not properly implemented yet. Whoopsy.")

        timesteps, wavelengths_count, lines_count, spectradata_string, datastart, columns \
            = self._spectra(scenario=scenario)
        
        spectra = self.spectra(time_lower=time, time_upper=time)
        lambdas = []
        vals = []
        filters = (wl_lower, wl_upper)
        line_wavelengths = spectra.columns[-lines_count:]
        for col in spectra.colnames:
            try:
                if col in columns:
                    continue
                if not lines and col in line_wavelengths:
                    continue
                if not continuum and col not in line_wavelengths:
                    continue
                l = float(col)
                l *= (1 + redshift)
                if l > filters[0] and l < filters[1]:
                    lambdas.append(l)
                    # if units == 'flux':
                    #     vals.append(spectra[col][0] * l)
                    # else:
                    vals.append(spectra[col][0])
            except ValueError:
                pass  # Not a wavelength, so ignore

        result_table = Table()
        result_table['wavelength'] = Column(lambdas, unit='Angstroms',
                                            description='wavelength')
        result_table['luminosity'] = Column(vals, unit='erg.s-1.A-1')

        return result_table

    def get_integrated_color(self, timestep, filter_name, redshift=0.0,
                             normalisation='AB'):
        """Calculates a color magnitude from spectral data with a given
        filter at a particular time.
        :param timestep: The timestep from spectra.dat
        :param filter_name: The filter to use
        :param redshift: Redshift to apply to the spectrum
        :param normalisation: Default AB magnitudes, can also use 'vega'"""
        if normalisation not in ['AB', 'vega']:
            raise ValueError("Normalisation must be 'AB' or 'vega'")

        fil = Filter.get_by_name(filter_name)
        wl_lower = fil.get_lower_bound()
        wl_upper = fil.get_upper_bound()

        spectrum = self.get_spectrum(time=timestep, redshift=redshift,
                                     wl_lower=wl_lower, wl_upper=wl_upper)
        lumo_ergs_s_integrated = 0.
        integrated_t_lambda = 0.
        c = 2.99792458e10

        wavelengths = np.array(spectrum['wavelength'])
        luminosities = np.array(spectrum['luminosity'])

        transmissions = np.array([fil.get_transmission_at_wavelength(w)
                                  for w in wavelengths])

#        lumo_ergs_s_integrated = np.sum(transmissions * luminosities)
        lumo_ergs_s_integrated = np.sum(transmissions * luminosities *
                                        wavelengths)

        integrated_t_lambda = np.sum((transmissions * c)/(wavelengths**2))
#        integrated_t_lambda = fil.get_calib['area']
        # lumo_ergs_s_integrated *= fil.get_midpoint()
        return -2.5 * (np.log10(lumo_ergs_s_integrated) - np.log10(integrated_t_lambda)) \
            - 48.6

    def save_to_file(self, filename=None):
        """
        Saves the PEGASE instance to disk, so it can be unpickled later and
        provide convenient wrapper around a set of parameters and output files.
        :param filename: The file to pickle to.
        :return: None
        """
        if filename is None:
            filename = self.name + ".peg"
        pickle.dump(self, open(filename, "wb"))

    def cleanup_files(self, silent=False):
        """
        Removes all the files associated with this run of PEGASE. Sets
        generated flag to False.
        :return: None
        """
        # Let's not use a wildcard so there's no surprises
        files_to_delete = []
        files_to_delete.extend(PEGASE._datafiles)
        files_to_delete.extend(["_colors" + str(i+1) + ".dat"
                                for i in range(len(self.scenarios))])
        files_to_delete.extend(["_spectra" + str(i+1) + ".dat"
                                for i in range(len(self.scenarios))])
        for file in files_to_delete:
            self._delete_file(file)
        self.generated = False

    def _delete_file(self, file_suffix, silent=True):
        try:
            os.remove(PEGASE.pegase_dir + self.name + file_suffix)
            if not silent:
                print "Removed %s" % (self.name + file_suffix)
        except OSError:  # as err:
            pass

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        result = "PEGASE \"%s\" [\n%s\nScenarios:\n" % (self.name, self.ssps)
        for i, scen in enumerate(self.scenarios):
            result += "Scenario %d - %s\n" % (i, scen)
        return result

if __name__ == "__main__":
    pass
