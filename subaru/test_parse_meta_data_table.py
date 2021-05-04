from unittest import TestCase
import os
from subaru import suprimecam2caom2
EXAMPLE_FILE = os.path.join(suprimecam2caom2.DATADIR, "SUP_Test_data.txt")


class TestParse_meta_data_table(TestCase):

    def setUp(self):
        self.test_fobj = open(EXAMPLE_FILE)
        self.test_values="""SUPA0008102X 2002-01-09             3 --- W-C-RC     22:49:42.418 +71:06:35.81 113.32975    10.54379     42.91397     65.08745         0.308644    -0.097749     0.946142 05:51:18.296     70.0   28.683 FOCUSING         IMAG             FOCUSING                               ----       ----       ---- """.split()
        self.test_header="""FRAME_ID    DATE_OBS   FITS_SIZE[MB] WCS FILTER     RA2000       DEC2000      GALLONG      GALLAT       ECLLONG      ECLLAT             X_2000       Y_2000       Z_2000 UT_STR        EXPTIME ALTITUDE DATA_TYP         OBS_MOD          OBJECT2                           PSF_SIGMA  PSF_ELLIP  SKY_LEVEL""".split()
        self.test_dict = {}
        for idx, key in enumerate(self.test_header):
            self.test_dict[key] = self.test_values[idx]

    def test_parse_meta_data_table(self):
        self.test_fobj.seek(0)
        table = suprimecam2caom2.parse_meta_data_table(self.test_fobj.name)
        test_line = table[table['FRAME_ID']==self.test_dict['FRAME_ID']]

        for key in self.test_dict.keys():
            try:
                value = float(test_line[key][0])
                self.assertAlmostEqual(test_line[key][0], value ,places=7)
            except:
                self.assertEqual(test_line[key][0], self.test_dict[key])
