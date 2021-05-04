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
CAOM2Collection = "SUBARU"

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

# Map SMOKA instrument identifiers to instrument short names.
INSTRUMENT_NAMES = {"MIR": "MOIRCS",
                    'HSC': "HSC",
                    "CIA": "CIAO",
                    "FCS": "FOCAS",
                    "HDS": "HDS",
                    "FMS": "FMOS",
                    "K3D": "Kyote3DII",
                    "HIC": "HiCIAO",
                    "SUP": "Suprime-Cam",
                    "COM": "COMICS",
                    "IRC": "IRCS",
                    "MCS": "MOIRCS",
                    }

INSTRUMENT_NICK_NAMES = INSTRUMENT_NAMES.keys()

PIXEL_SCALE = {"SUP": 0.2 * units.arcsecond,
               "FCS": 0.1038 * units.arcsecond,
               "CIA": 0.0213 * units.arcsecond,
               "HSC": 0.17 * units.arcsecond,
               "MCS": 0.116 * units.arcsecond}

SPECTRAL_RESPONSE = {"FCS": { "min_wavelength": 400 *units.nm, "max_wavelength": 1100*units.nm},
                     "DEFAULT":{ "min_wavelength": 300 *units.nm, "max_wavelength": 5000*units.nm},
                     "HSC": {"min_wavelength": 300 * units.nm, "max_wavelength": 1050 * units.nm},
                     "CIA": {"min_wavelength": 900 * units.nm, "max_wavelength": 5500 * units.nm},
                     "SUP": { "min_wavelength": 300 *units.nm, "max_wavelength": 1050*units.nm}}


class SMOKA(object):
    """
    A SMOKA data table based CAOM2 SimpleObservation
    """

    def __init__(self, row, previous_observation=None):
        """

        :param row: A row from the SMOKA obsLog
        :param previous_observation: a previous record, which we might be adding to.
        :type row: Row
        """
        self.row = row
        self.bandpass_database = svo.BandpassFilterDatabase(
            telescope=TELESCOPE_NAME.get(self.instrument_name, 'SUBARU'),
            instrument=self.instrument_name)
        filters.augment_filter_lookup_table(self.bandpass_database)
        # the previous record, which we might be adding to.
        self.previous_observation = previous_observation

    @property
    def collection(self):
        return CAOM2Collection

    @property
    def instrument_name(self):
        return u"{}".format(self.row['FRAME_ID'][0:3])

    @property
    def product_id(self):
        return self._exposure_id.replace(self.instrument_name + "E", self.instrument_name + "A")

    @property
    def _exposure_id(self):
        """
        determine the SMOKA exposure_id from the frame_id
        :return:
        """
        exposure_id = self.row['FRAME_ID']
        if self.instrument_name in ["MSC", "FCS"]:
            # MCS and FCS are TWO CCD devices. The odd numbered FRAME_ID is the EXPOSURE_ID
            dig = re.search('(\D+)([0-9]+)', exposure_id).groups()
            base = dig[0]
            dig = dig[1]
            chip = int(dig)
            if not chip % 2 > 0 and self.instrument_name in ['MCS', 'FCS']:
                exposure_id = base+str(chip-1).zfill(len(dig))
        if self.instrument_name in ["SUP", "CIA", "MCS", "FCS"]:
            # For these instruments the 'A' in the frAme_id becomes an 'E' for the Exposure_id
            exposure_id = exposure_id.replace(self.instrument_name + "A",
                                              self.instrument_name + "E")
        return exposure_id

    @property
    def observation_id(self):
        return self._exposure_id.replace("XX", "00").replace("X", "0").replace("Y", "1")

    @property
    def type(self):
        return u"{}".format(self.row['DATA_TYP'].upper())

    @property
    def target_name(self):
        return u'{}'.format(self.row['OBJECT2']).upper()

    @property
    def intent(self):
        if self.type == u'OBJECT':
            if self.target_name in ['FLAT', 'BIAS', 'ZERO', 'COMPARISON', "DARK", 'FRINGE']:
                return ObservationIntentType.CALIBRATION
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
            key = "{}{}".format(key, self.product_id[-1])
        logging.debug("FOV KEY: {}".format(key))
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
        """
        create and return the caom2.Position object for this observation
        :return:
        """
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
                logging.debug("Attempting to disentangle filter name: {}".format(filter_name))
                if "NONE" in filter_name or "open" in filter_name:
                    continue
                filter_name = filter_name.split('-')[-1]
                logging.debug("After splitting: {}".format(filter_name))
                if len(filter_name) > 0:
                    if filter_name[-1] == "+":
                        filter_name = "SDSS_{}".format(filter_name[0])
                    if filter_name[0] == 'L':
                        filter_name = "NB{}".format(filter_name[1:])
                    logging.debug("After mapping to SVO style: {}".format(filter_name))

                # The headers keywords have extra characters when compared to the SVO database.
                group = re.match('^SCFCFL[BS]?(.*?)(01)?$', filter_name)
                if group is not None:
                    logging.debug("Group Matching: {}".format(filter_name))
                    filter_name = group.groups()[0]
                logging.debug("adding name to list: {}".format(filter_name))
                _filters.append(filter_name)
        return _filters

    @property
    def disperser(self):
        if "DISPERSR" not in self.row.colnames:
            return None
        if self.row['DISPERSR'] not in filters.GRISM_NAMES:
            return None
        return self.row['DISPERSR']

    @property
    def grism(self):
        grism_name = filters.GRISM_NAMES.get(self.disperser, None)
        if grism_name is None:
            logging.warning("Unkown Grism Name: {}".format(self.disperser))
            return None
        if grism_name not in filters.DISPERSERS:
            logging.warning("No filter details for grism: {}".format(grism_name))
            return None
        return filters.DISPERSERS[grism_name]

    @property
    def energy(self):

        energy = Energy()
        energy.em_band = EnergyBand.OPTICAL
        energy.dimension = 1

        # Build the energy object
        max_wavelength = SPECTRAL_RESPONSE.get(self.instrument_name, SPECTRAL_RESPONSE["DEFAULT"])["max_wavelength"]
        min_wavelength = SPECTRAL_RESPONSE.get(self.instrument_name, SPECTRAL_RESPONSE["DEFAULT"])["min_wavelength"]

        for filter_name in self.filters:
            logging.debug("Looking up filter information for : {}".format(filter_name))
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
            logging.debug("Grism: {}, Params: {}".format(self.grism, grism_params))
            min_wavelength = max(min_wavelength, grism_params[0] * units.AA)
            max_wavelength = min(max_wavelength, grism_params[1] * units.AA)
            logging.debug("min/max: {} {}".format(min_wavelength, max_wavelength))
            resolving_power = grism_params[2]
            bandwidth = max_wavelength - min_wavelength

        min_wavelength = min(min_wavelength, max_wavelength)
        max_wavelength = max(min_wavelength, max_wavelength)
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
        if (self.disperser is None or self.disperser == 'null') and obs_mod == 'SPEC':
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
                             IMAG_VGW=DataProductType.IMAGE,
                             SPEC=DataProductType.SPECTRUM,
                             IMAG_SINGLE=DataProductType.IMAGE,
                             COMPARISON=DataProductType.SPECTRUM,
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
        this_instrument = Instrument(name=INSTRUMENT_NAMES.get(self.instrument_name, self.instrument_name))
        this_instrument.keywords.add(str(self.row['OBS_MOD']))
        if self.disperser is not None:
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

    @property
    def obsdate(self):
        """
        Get the observation date in a format that CAOM2 wants.
        :return: obsdate
        """
        return self.start_time.iso[0:10]

    @property
    def science_artifact_uri(self):
        return 'subaru:raw/{}/{}'.format(self.obsdate,
                                         self.row['FRAME_ID'].replace("Y", "X"))

    @property
    def preview_artifact_uri(self):
        preview_id = self.product_id.replace("XX", "00").replace("X", "1").replace("Y", "1")
        return 'subaru:preview/{}/{}'.format(self.obsdate,
                                             preview_id)

    def __call__(self):

        logging.info("Creating CAOM2 observation: {}".format(self.observation_id))

        if isinstance(self.previous_observation, Observation):
            this_observation = self.previous_observation
        else:
            this_observation = SimpleObservation(collection=self.collection,
                                                 observation_id=self.observation_id)
        this_observation.meta_release = self.meta_release
        this_observation.target = self.target
        this_observation.instrument = self.instrument
        this_observation.telescope = self.telescope
        this_observation.sequence_number = None
        this_observation.intent = self.intent
        this_observation.type = self.type
        this_observation.proposal = None

        # Get or create as needed the plane that will hold the raw data
        if self.product_id not in this_observation.planes:
            this_observation.planes.add(Plane(self.product_id))
        this_plane = this_observation.planes.get(self.product_id)
        this_plane.meta_release = self.meta_release
        this_plane.data_release = self.data_release
        this_plane.calibration_level = CalibrationLevel.RAW_STANDARD
        this_plane.data_product_type = self.data_product_type
        this_plane.provenance = self.providence
        this_plane.time = self.time
        # Build the position object.
        this_plane.position = self.position
        this_plane.energy = self.energy

        # create a reference to the data file stored at SMOKA.
        # Get or create as needed the artifact that will hold the raw data
        if self.science_artifact_uri not in this_plane.artifacts:
            this_plane.artifacts.add(Artifact(self.science_artifact_uri,
                                              ProductType.SCIENCE, ReleaseType.DATA,
                                              content_type='text/html'))
        else:
            science_artifact = this_plane.artifacts.get(self.science_artifact_uri)
            science_artifact.product_type = ProductType.SCIENCE
            science_artifact.release_type = ReleaseType.DATA
            science_artifact.content_type = 'text/html'

        # Modify/Create as needed the PREVIEW artifact
        if self.preview_artifact_uri not in this_plane.artifacts:
            this_plane.artifacts.add(Artifact(uri=self.preview_artifact_uri,
                                              product_type=ProductType.PREVIEW,
                                              release_type=ReleaseType.META,
                                              content_type='image/png'))
        else:
            preview_artifact = this_plane.artifacts.get(self.preview_artifact_uri)
            preview_artifact.product_type = ProductType.PREVIEW
            preview_artifact.release_type = ReleaseType.META
            preview_artifact.content_type = 'image/png'

        return this_observation
