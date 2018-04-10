from __future__ import unicode_literals
import re
from caom2 import *
from astropy import time
from astropy import units
import svo
import logging
import os
import wcs2bounds
import filters
from astropy.coordinates import Angle

VERSION_DATE = time.Time(os.stat(__file__).st_mtime, format='unix')

DATADIR = __PATH__ = os.path.join(os.path.dirname(__file__), 'data')

SMOKA_OBSLOG_ENDPOINT = 'http://smoka.nao.ac.jp/status/obslog'

TELESCOPE_NAME = {"SUP": 'SUBARU',
                  "FCS": 'SUBARU',
                  "CIA": 'SUBARU',
                  'HSC': 'SUBARU',
                  'MCS': 'SUBARU'}

PROVENANCE_URLS = {"SUP": 'https://www.subarutelescope.org/Observing/Instruments/SCam/index.html',
                   "FCS": 'https://www.subarutelescope.org/Observing/Instruments/FOCAS/index.html',
                   "CIA": 'https://www.subarutelescope.org/Observing/Instruments/CIAO/index.html',
                   'HSC': 'https://www.subarutelescope.org/Observing/Instruments/HSC/index.html',
                   'MCS': 'https://www.subarutelescope.org/Observing/Instruments/MOIRCS/index.html'
                   }

INSTRUMENT_NAMES = {"MIR": ("MOIRCS", "Multi-Object Infrared Camera and Spectrograph"),
                    'HSC': ("HSC", "Hyper Suprime-Cam"),
                    "CIA": ("CIAO", "Coronagraphic Imager with Adaptive Optics"),
                    "FCS": ("FOCAS", "Faint Object Camera And Spectrograph"),
                    "HDS": ("HDS", "High Dispersion Spectrograph"),
                    "FMS": ("FMOS", "Fiber Multi Object Spectrograph"),
                    "K3D": ("Kyote3DII", "Kyoto tridimensional spectrograph II"),
                    "HIC": ("HiCIAO", "High Contrast Instrument for the Subaru Next Generation Adaptive Optics"),
                    "SUP": ("Suprime", 'Subaru Prime Focus Camera'),
                    "CAC": "CAC",
                    "COM": "Cooled Mid-Infrared Camera and Spectrograph",
                    "IRC": "Infrared Camera and Spectrograph",
                    "OHS": "Cooled Infrared Spectrograph and Camera for OHS",
                    "MCS": ("MOIRCS", "Multi-Object InfraRed Camera and Spectrograph")
                    }

INSTRUMENT_NICK_NAMES = INSTRUMENT_NAMES.keys()

PIXEL_SCALE = {"SUP": 0.2 * units.arcsecond,
               "FCS": 0.1038 * units.arcsecond,
               "CIA": 0.0213 * units.arcsecond,
               "HSC": 0.17 * units.arcsecond,
               "MCS": 0.116 * units.arcsecond}


class SMOKA(object):
    """
    A SMOKA data table based CAOM2 SimpleObservation
    """

    def __init__(self, row):
        """

        :param row: A row from the SMOKA obsLog
        :type row: Row
        """
        self.row = row
        self.bandpass_database = svo.BandpassFilterDatabase(
            telescope=TELESCOPE_NAME.get(self.instrument_name, 'SUBARU'),
            instrument=INSTRUMENT_NAMES.get(self.instrument_name, None)[0])
        filters.augment_filter_lookup_table(self.bandpass_database)

    @property
    def collection(self):
        return "SUBARU"

    @property
    def instrument_name(self):
        return u"{}".format(self.row['FRAME_ID'][0:3])

    @property
    def product_id(self):
        product_id = u"{}".format(self.row['FRAME_ID'])
        return product_id

    @property
    def observation_id(self):
        observation_id = u"{}".format(self.row['FRAME_ID'])
        observation_id = observation_id.replace("X", "0")
        observation_id = observation_id.replace("Y", "1")
        if self.instrument_name in ["SUP", "CIA", "MCS", "FCS"]:
            observation_id = observation_id.replace(self.instrument_name + "A",
                                                    self.instrument_name + "E")
        return observation_id

    @property
    def obstype(self):
        return u"{}".format(self.row['DATA_TYP'].upper())

    @property
    def target_name(self):
        return u'{}'.format(self.row['OBJECT2']).upper()

    @property
    def intent(self):
        if self.obstype == u'OBJECT':
            return ObservationIntentType.SCIENCE
        return ObservationIntentType.CALIBRATION

    @property
    def fov(self):
        """
        The Field of View of the instrument assocaited with this SMOKA row.
        :return: fov
        """
        fov = dict(FCS=(6 / 60. / 2.0, 6 / 60. / 2.0),
                   CIA=(21.8 / 60.0 / 2.0, 21.8 / 60.0 / 2.0),
                   SUPX=(34 / 60.0, 27 / 60.0),
                   SUPY=(4 * 34 / (5 * 60.0), 27 / 60.0),
                   HSC=(90/60.0, 90/60.0),
                   MCS=(4/60.0, 7/60.0))

        key = self.instrument_name
        if key == "SUP":
            "{}{}".format(key, self.product_id[-1])
        return fov.get(key, (0.1, 0.1))

    @property
    def ra(self):
        try:
            return Angle(str(self.row['RA2000']), unit='hour')
        except Exception as ex:
            logging.error(str(ex))
            logging.error("\n->"+str(self.row['RA2000'])+"<-")
            return 0.0*units.degree

    @property
    def dec(self):
        try:
            return Angle(str(self.row['DEC2000']), unit='degree')
        except Exception as ex:
            logging.error(str(ex))
            logging.error("\n->"+str(self.row['DEC2000'])+"<-")
            return 0.0*units.degree

    @property
    def position_bounds(self):
        """
        Given the RA and DEC of a SuprimeCam image return a Polygon bounding box.

        :return: bounding box
        :rtype: CoordPolygon2D
        """
        width = self.fov[0]
        height = self.fov[1]
        naxis1 = width / PIXEL_SCALE[self.instrument_name].to('degree').value
        naxis2 = height / PIXEL_SCALE[self.instrument_name].to('degree').value
        wcs_header = dict(crval1=self.ra.deg,
                          crval2=self.dec.deg,
                          crpix1=naxis1 / 2.0,
                          crpix2=naxis2 / 2.0,
                          cd1_1=PIXEL_SCALE["SUP"].to('degree').value,
                          cd2_2=PIXEL_SCALE["SUP"].to('degree').value,
                          cd1_2=0, cd2_1=0, cunit1='deg', cunit2='deg', ctype1='RA---TAN',
                          ctype2='DEC--TAN', naxis1=naxis1, naxis2=naxis2)

        return wcs2bounds.wcs2bounds(wcs_header, axes=(naxis1, naxis2))

    @property
    def position(self):
        return Position(bounds=self.position_bounds,
                        sample_size=PIXEL_SCALE[self.instrument_name].to('arcsecond').value,
                        time_dependent=False)

    @property
    def filters(self):
        """Get all filter values from the SMOKA obslog row"""
        _filters = []
        for filter_no in ["", "01", "02", "03", "04"]:
            key = "FILTER{}".format(filter_no)
            if key in self.row.colnames:
                filter_name = self.row[key]
                logging.info("Attempting to disentangle filter name: {}".format(filter_name))
                if "NONE" in filter_name or "open" in filter_name:
                    continue
                filter_name = filter_name.split('-')[-1]
                logging.info("After splitting: {}".format(filter_name))
                if len(filter_name) > 0:
                    if filter_name[-1] == "+":
                        filter_name = "SDSS_{}".format(filter_name[0])
                    if filter_name[0] == 'L':
                        filter_name = "NB{}".format(filter_name[1:])
                    logging.info("After mapping to SVO style: {}".format(filter_name))

                # The headers keywords have extra characters when compared to the SVO database.
                group = re.match('^SCFCFL[BS]?(.*?)(01)?$', filter_name)
                if group is not None:
                    filter_name = group.groups()[0]
                else:
                    # No matching filter bits.
                    continue
                _filters.append(filter_name)
                logging.info("adding name to list: {}".format(filter_name))
        return _filters

    @property
    def disperser(self):
        if "DISPERSR" not in self.row.colnames:
            return None
        return self.row['DISPERSR']

    @property
    def grism(self):
        grism_name = filters.GRISM_NAMES.get(self.disperser, None)
        return filters.DISPERSERS.get(grism_name, None)

    @property
    def energy(self):

        energy = Energy()
        energy.em_band = EnergyBand.OPTICAL
        energy.dimension = 1

        # Build the energy object
        max_wavelength = 5.0 * units.um
        min_wavelength = 0.3 * units.um

        for filter_name in self.filters:
            filter_info = self.bandpass_database[filter_name]
            if filter_info is None:
                continue
            energy.bandpass_name = u'{}'.format(filter_name)
            max_wavelength = min(max_wavelength, filter_info['wavelength_max'])
            min_wavelength = max(min_wavelength, filter_info['wavelength_min'])
            logging.debug("{}: {} {}".format(filter_name, filter_info['wavelength_min'], filter_info['wavelength_max']))
            logging.debug("min/max: {} {}".format(min_wavelength, max_wavelength))

        cwl = (min_wavelength + max_wavelength) / 2.0
        bandwidth = (max_wavelength - min_wavelength)
        resolving_power = float(cwl / bandwidth)

        if self.grism is not None:
            grism_params = self.grism.get("DEFAULT")
            for filter_name in self.filters:
                if filter_name in self.grism.keys():
                    grism_params = self.grism.get(filter_name)
            min_wavelength = max(min_wavelength, grism_params[0] * units.AA)
            max_wavelength = min(max_wavelength, grism_params[1] * units.AA)
            resolving_power = grism_params[2]
            bandwidth = max_wavelength - min_wavelength

        energy.bounds = Interval(min_wavelength.to('m').value,
                                 max_wavelength.to('m').value,
                                 samples=[shape.SubInterval(min_wavelength.to('m').value,
                                                            max_wavelength.to('m').value)])
        energy.resolving_power = float(resolving_power)
        energy.sample_size = bandwidth.to('m').value

        return energy

    @property
    def obs_mod(self):
        obs_mod = self.row['OBS_MOD']
        if self.disperser is None and obs_mod == 'SPEC':
            obs_mod = 'IMAG'
        return obs_mod

    @property
    def data_product_type(self):
        """
        compute the DataProductType based on the SMOKA obs_mod keyword.

        :return: the data product type
        :rtype: DataProductType
        """
        data_type_map = dict(IMAG=DataProductType.IMAGE,
                             SPEC=DataProductType.SPECTRUM,
                             IMAG_SINGLE=DataProductType.IMAGE,
                             MOS=DataProductType.SPECTRUM)
        return data_type_map.get(self.obs_mod, None)

    @property
    def telescope(self):
        this_telescope = Telescope(name=TELESCOPE_NAME.get(self.instrument_name, 'Subaru'))
        # These are actually the JCMT values.
        this_telescope.geo_location_x = -5461060.909
        this_telescope.geo_location_y = -2491393.621
        this_telescope.geo_location_z = 2149257.916
        return this_telescope

    @property
    def instrument(self):
        # First lets define the Observation that is this record.
        this_instrument = Instrument(name=INSTRUMENT_NAMES[self.instrument_name][0])
        this_instrument.keywords.add(str(self.row['OBS_MOD']))
        this_instrument.keywords.add(str(self.disperser))
        return this_instrument

    @property
    def target(self):
        return Target(self.target_name)

    @property
    def meta_release(self):
        return time.Time('2017-01-01 00:00:00').to_datetime()

    @property
    def data_release(self):
        return self.meta_release

    @property
    def providence(self):
        return Provenance(name='SMOKA',
                          producer='SMOKA',
                          project='SMOKA',
                          reference=PROVENANCE_URLS[self.instrument_name])

    @property
    def start_time(self):
        return time.Time("{} {}".format(self.row['DATE_OBS'], self.row['UT_STR']))

    @property
    def exptime(self):
        return self.row['EXPTIME'] * units.second

    @property
    def end_time(self):
        return time.Time(self.start_time + self.exptime)

    @property
    def time(self):
        # Build the time object.
        time_bounds = Interval(self.start_time.mjd,
                               self.end_time.mjd,
                               samples=[shape.SubInterval(self.start_time.mjd, self.end_time.mjd)])
        return Time(bounds=time_bounds,
                    dimension=1,
                    resolution=self.exptime.to('second').value,
                    sample_size=self.exptime.to('day').value,
                    exposure=self.exptime.to('second').value
                    )

    def __call__(self):

        logging.info("Creating CAOM2 observation: {}".format(self.observation_id))

        this_observation = SimpleObservation(collection=self.telescope.name,
                                             observation_id=self.observation_id,
                                             sequence_number=None,
                                             intent=self.intent,
                                             type=self.obstype,
                                             proposal=None,
                                             telescope=self.telescope,
                                             instrument=self.instrument,
                                             target=self.target,
                                             meta_release=self.meta_release
                                             )

        # Create a plane that will hold the raw data
        this_plane = Plane(self.product_id,
                           meta_release=self.meta_release,
                           data_release=self.data_release,
                           )

        this_plane.calibration_level = CalibrationLevel.RAW_STANDARD
        this_plane.data_product_type = self.data_product_type
        this_plane.provenance = self.providence
        this_plane.time = self.time
        # Build the position object.
        this_plane.position = self.position
        this_plane.energy = self.energy

        # create a reference to the data file stored at SMOKA.
        obsdate = self.start_time.iso[0:10]
        this_plane.artifacts.add(Artifact(uri='subaru:raw/{}/{}'.format(obsdate,
                                                                        this_plane.product_id),
                                          product_type=ProductType.SCIENCE,
                                          content_type='text/html',
                                          release_type=ReleaseType.DATA))

        # Add the PREVIEW artifact, stored in SMOKA, These are referenced using the FrameID value of one of the members.
        # but the frameID has an 'X' or 'Y' at end for SUP, so we need to adjust here.
        this_plane.artifacts.add(Artifact(uri='subaru:preview/{}/{}'.format(obsdate,
            this_plane.product_id.replace("X", "1").replace("Y", "1")),
            product_type=ProductType.PREVIEW,
            release_type=ReleaseType.META,
            content_type='image/png'))

        # And the plane to the observation.
        this_observation.planes.add(this_plane)

        return this_observation
