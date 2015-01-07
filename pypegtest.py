import unittest
import os.path
from pypegase import PEGASE, SFR, IMF, Extinction, SNII_ejecta, read_file

FULLTEST = FALSE

class TestPyPegase(unittest.TestCase):
    files = [
        '_SSPs.dat', '_tracksZ0.004.dat', '_colors1.dat', '_tracksZ0.008.dat',
        '_scenarios.dat', '_tracksZ0.02.dat', '_tracksZ0.0001.dat', '_tracksZ0.05.dat',
        '_spectra1.dat', '_tracksZ0.0004.dat', '_tracksZ0.1.dat'
    ]

    def assertExists(self, filename, msg=None):
        if not os.path.isfile(filename):
            msg = self._formatMessage(msg, "%s does not exist" % filename)
            raise self.failureException(msg)

    def assertNotExists(self, filename, msg=None):
        if os.path.isfile(filename):
            msg = self._formatMessage(msg, "%s exists but should not." % filename)
            raise self.failureException(msg)
           
    def setUp(self):
        self.name = "utest001"
        name = self.name
        try:
            self.peg = PEGASE.from_file(name + ".peg")
        except IOError:
            self.peg = PEGASE(name)

    def test_001_PEGASEDefaults(self):
        peg = PEGASE(self.name+"new")
        self.assertEquals(peg.ssps.imf.number, IMF.IMF_Salpeter)
        self.assertTrue(peg.ssps.winds)
        self.assertEquals(peg.ssps.ejecta, SNII_ejecta.MODEL_B)
        self.assertEquals(peg.scenarios[0].binaries_fraction, 0.05)
        self.assertEquals(peg.scenarios[0].sfr.sfrtype,
                          SFR.EXPONENTIAL_DECREASE)
        self.assertEquals(peg.scenarios[0].sfr.p1, 1000)
        self.assertEquals(peg.scenarios[0].sfr.p2, 1)
        self.assertEquals(peg.scenarios[0].extinction,
                          Extinction.NO_EXTINCTION)

    def test_010_Pegase(self):
        peg = self.peg
        if FULLTEST:
            peg = PEGASE(self.name)
            
        if not peg.is_generated():
            peg.generate()
            self.assertTrue(peg.is_generated())
            peg.save_to_file(peg.name + ".peg")

        for file in TestPyPegase.files:
            self.assertExists(PEGASE.pegase_dir + peg.name + file)

        ssps_file = read_file(peg.name + "_SSPs.dat")
        self.assertEquals("utest001_tracksZ0.0001.dat \nutest001_tracksZ0.0004.dat \nutest001_tracksZ0.004.dat \n" +
        "utest001_tracksZ0.008.dat \nutest001_tracksZ0.02.dat \nutest001_tracksZ0.05.dat \nutest001_tracksZ0.1.dat \n",
                          ssps_file)

        scenarios_file = read_file(peg.name + "_scenarios.dat")
        self.assertEquals("SSPs file: utest001_SSPs.dat\nFraction of close binary systems: 0.50000E-01\n" +
        "************************************************************\n  1: utest001_spectra1.dat\nInitial metallicity: 0.00000E+00\n" +
        "No infall\nType of star formation:   2\np1: 0.10000E+04\np2: 0.10000E+01\nConsistent evolution of the stellar metallicity\n" +
        "Mass fraction of substellar objects: 0.00000E+00\nNo galactic winds\nNebular emission\nNo extinction\n",
                          scenarios_file)


        colors = peg.colors()
        self.assertIsNotNone(colors)
        self.assertAlmostEqual(colors['MBHNS'][10], .416E-04)
        self.assertAlmostEqual(colors['W(Ha)'][65], .966E-04)
        self.assertAlmostEqual(colors['3150-B'][68], 0.663)


        spectra = peg.spectra()  # 
        self.assertIsNotNone(spectra)
        self.assertAlmostEqual(spectra['9069.00'][45], .2426E32)
        self.assertAlmostEqual(spectra['9069.00'][65], .4354E25)
        self.assertAlmostEqual(spectra['95000.'][45], .2413E27)

        
    def test_015_Pickle(self):
        pass

    def test_090_CleanUp(self):
        if FULLTEST:
            peg = PEGASE("SPOD")
            peg.generate()
            peg.cleanup_files()
            for file in TestPyPegase.files:
                self.assertNotExists("SPOD_" + file)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
