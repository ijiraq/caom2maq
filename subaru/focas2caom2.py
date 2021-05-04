from __future__ import unicode_literals
import argparse

from cadcutils import exceptions as cadcutils_execptions
import re
from astropy.coordinates import SkyCoord
from cadcutils import net
from caom2 import *
from astropy import time
from astropy import units
from astropy.io import ascii
import svo
import requests
import logging
import os
from tempfile import NamedTemporaryFile
from caom2repo import CAOM2RepoClient
import wcs2bounds
import tap

VERSION_DATE = time.Time('2018-04-08T07:10:00')

DATADIR = __PATH__ = os.path.join(os.path.dirname(__file__), 'data')

INSTRUMENT_NICK_NAMES = ('MIR', 'HSC', 'CIA', 'FCS', 'HDS', 'MCS', 'FMS', 'K3D',
                         'HIC', 'SUP', 'CAC', 'IRC', 'COM', 'OHS')
SMOKA_OBSLOG_ENDPOINT = 'http://smoka.nao.ac.jp/status/obslog'
TELESCOPE_NAME = {"SUP": 'SUBARU',
                  "FCS": 'SUBARU',
                  "CIA": 'SUBARU'}
PROVENANCE_URLS = {"SUP": 'https://www.subarutelescope.org/Observing/Instruments/SCam/index.html',
                   "FCS": 'https://www.subarutelescope.org/Observing/Instruments/FOCAS/index.html',
                   "CIA": 'https://www.subarutelescope.org/Observing/Instruments/CIAO/index.html'
                   }

INSTRUMENT_NAMES = {"MIR": ("MOICS", "Multi-Object Infrared Camera and Spectrograph"),
                    'HSC': ("HSC", "Hyper Suprime-Cam"),
                    "CIA": ("CIAO", "Coronagraphic Imager with Adaptive Optics"),
                    "FCS": ("FOCAS", "Faint Object Camera And Spectrograph"),
                    "HDS": ("HDS", "High Dispersion Spectrograph"),
                    "MCS": ("MOIRCS", "Multi-Object Infrared Camera and Spectrograph"),
                    "FMS": ("FMOS", "Fiber Multi Object Spectrograph"),
                    "K3D": ("Kyote3DII", "Kyoto tridimensional spectrograph II"),
                    "HIC": ("HiCIAO", "High Contrast Instrument for the Subaru Next Generation Adaptive Optics"),
                    "SUP": ("Suprime", 'Subaru Prime Focus Camera'),
                    "CAC": "CAC",
                    "COM": "Cooled Mid-Infrared Camera and Spectrograph",
                    "IRC": "Infrared Camera and Spectrograph",
                    "OHS": "Cooled Infrared Spectrograph and Camera for OHS",
                    }
DISPERSERS = {"75": {"SY47": (4700, 9100, 250),
                     "SO58": (5800, 10000, 250),
                     "DEFAULT": (3000, 11000, 250)},
              "150": {"L550": (3400, 5500, 700),
                      "SO58": (5800, 10000, 500),
                      "SY47": (4700, 9100, 500),
                      "DEFAULT": (3000, 11000, 600)},
              "300B": {"SY47": (4700, 9100, 1000),
                       "L600": (3700, 6000, 1000),
                       "DEFAULT": (3000, 11000, 1000)},
              "300R": {"SY47": (4900, 9100, 1000),
                       "SO58": (5800, 10000, 1000),
                       "L600": (3700, 5950, 2000),
                       "L550": (3400, 5250, 2000),
                       "DEFAULT": (3000, 11000, 1500)},
              "Echelle": {"I": (7300, 8700, 2500),
                          "z": (8300, 10000, 2500),
                          "DEFAULT": (3000, 11000, 2500)},
              "VPH450": {"NONE": (3800, 5250, 3000)},
              "VPH520": {"NONE": (4450, 6050, 3000)},
              "VPH650": {"SY47": (5300, 7700, 2500),
                         "DEFAULT": (3000, 11000, 2500)},
              "VPH850": {"SO58": (5800, 10350, 1500),
                         "DEFAULT": (3000, 11000, 2500)},
              "VPH900": {"SO58": (7500, 10450, 3000),
                         "DEFAULT": (3000, 11000, 2500)},
              "VPH680": {"SY47": (6450, 7350, 7500),
                         "DEFAULT": (3000, 11000, 2500)},
              "VPH800": {"SY47": (7500, 8600, 7000),
                         "DEFAULT": (3000, 11000, 2500)},
              "VPH950": {"O58": (8850, 10000, 5500),
                         "DEFAULT": (3000, 11000, 2500)}
              }
_GRISM_NAMES = (
    ("SCFCGREL01", "SCFCGRLD01", "SCFCGRMB01", "SCFCGRMR01", "SCFCGRHDEC", "SCFCGRHD45", "SCFCGRHD52", "SCFCGRHD65",
     "SCFCGRHD68", "SCFCGRHD80", "SCFCGRHD95"),
    ("75", "150", "300B", "300R", "Echelle", "VPH450", "VPH520", "VPH650", "VPH680", "VPH800", "VPH950")
)
GRISM_NAMES = {}
for idx in range(len(_GRISM_NAMES[0])):
    GRISM_NAMES[_GRISM_NAMES[0][idx]] = _GRISM_NAMES[1][idx]
PIXEL_SCALE = {"SUP": 0.2 * units.arcsecond,
               "FCS": 0.1038 * units.arcsecond,
               "CIA": 0.0213 * units.arcsecond}

FOV = dict(FCS=(6 / 60. / 2.0, 6 / 60. / 2.0),
           CIA=(21.8 / 60.0 / 2.0, 21.8 / 60.0 / 2.0),
           SUPX=(34 / 60.0, 27 / 60.0),
           SUPY=(4 * 34 / (5 * 60.0), 27 / 60.0))

FILTER_MAP = dict(SCFCFLBU01='U',
                  SCFCFLBB01='B',
                  SCFCFLBV01='V',
                  SCFCFLBR01='R',
                  SCFCFLBI01='I',
                  SCFCFLN373='N373',
                  SCFCFLN386='N386',
                  SCFCFLN487='N487',
                  SCFCFLN502='N502',
                  SCFCFLN512='N512',
                  SCFCFLN642='N642',
                  SCFCFLN670='N670',
                  SCFCFLBSZ1='SDSS_z',
                  SCFCFLSO58='O58',
                  SCFCFLSY47='Y47',
                  SCFCFLL600='L600',
                  SCFCFLL550='L550',
                  SCFCFLLC50='C50'
                  )


def get_metadata_table(instrument_name="SUP", year=2002):
    """
    retrieve the SMOKA metadata table from VOSpace.

    :param year: Year of observation records to retrieve from SMOKA
    :type year: basestring
    :param instrument_name: Instrument (short form) to retrieve observations records from SMOKA
    :type instrument_name: basestring
    :rtype: Table
    :return: table of metadata values
    """
    logging.info("Building metadata table for : {} {}".format(instrument_name, year))
    local_filename = "{}_{}.txt".format(instrument_name, year)
    if not os.access(local_filename, os.R_OK):
        smoka_obslog = "{}/{}".format(SMOKA_OBSLOG_ENDPOINT, local_filename)
        name_temp_file = NamedTemporaryFile()
        local_filename = str(name_temp_file.name)
        with open(local_filename, str('w')) as fobj:
            logging.info("Connecting to: {}".format(smoka_obslog))
            resp = requests.get(smoka_obslog)
            resp.raise_for_status()
            fobj.write(resp.content)
        name_temp_file.seek(0)
    return parse_meta_data_table(local_filename, instrument_name)


def parse_meta_data_table(local_filename, instrument_name):
    """
    Use the astropy package to return a Table representation of a SMOKA observation record.

    :param instrument_name: Name of the instrument to to get observing log for.
    :param local_filename: name of text file with SMOKA observations in a text file.
    :type local_filename: str
    :type instrument_name: str
    :return: Observation Record
    :rtype: Table
    """
    header = os.path.join(DATADIR, "{}_Header.txt".format(instrument_name))
    position_line = open(header).readlines()[1].split(' ')
    col_ends = []
    col_starts = [0]
    for col in position_line:
        if col:
            col_ends.append(col_starts[-1] + len(col))
            col_starts.append(col_ends[-1] + 1)
        else:
            col_starts[-1] += 1
    col_starts = col_starts[:-1]

    column_names = open(local_filename).readline()[1:].split()
    return ascii.read(local_filename,
                      format="fixed_width_no_header",
                      col_ends=col_ends,
                      col_starts=col_starts,
                      names=column_names,
                      fill_values=('', '----', '---'))


def position_bounds(ra, dec, instrument_name):
    """
    Given the RA and DEC of a SuprimeCam image return a Polygon bounding box.

    :param ra: RA of the centre of the FOV (degrees)
    :type ra: float
    :param dec: DEC of the centre of the FOV (degrees)
    :type dec: float
    :param instrument_name: Short form name of instrumnet, used to select w/h parameters.
    :type instrument_name: basestring
    :return: bounding box
    :rtype: CoordPolygon2D
    """
    width = FOV[instrument_name][0]
    height = FOV[instrument_name][1]
    naxis1 = width / PIXEL_SCALE[instrument_name].to('degree').value
    naxis2 = height / PIXEL_SCALE[instrument_name].to('degree').value
    wcs_header = dict(crval1=ra, crval2=dec, crpix1=naxis1 / 2.0, crpix2=naxis2 / 2.0,
                      cd1_1=PIXEL_SCALE["SUP"].to('degree').value,
                      cd2_2=PIXEL_SCALE["SUP"].to('degree').value,
                      cd1_2=0, cd2_1=0, cunit1='deg', cunit2='deg', ctype1='RA---TAN',
                      ctype2='DEC--TAN', naxis1=naxis1, naxis2=naxis2)

    return wcs2bounds.wcs2bounds(wcs_header, axes=(naxis1, naxis2))


def smoka_datarequest(frame_ids):
    """
    Make a datarequest given a SMOKA frame id.

    This is method is in here as an exmple for datalink service to follow.

    :param frame_ids: list of frames to retrieve
    :type frame_ids: list
    :return: The content of the resulting POST to the smoka datarequest service.
    """

    endpoint = """http://smoka.nao.ac.jp/datarequest"""

    data = [('action', 'Datarequest'), ('search_type', '')]
    for frame_id in frame_ids:
        data.append(('frameinfo', frame_id))

    return requests.post(endpoint, data=data).content


def build_energy(row, bandpass_database):
    """
    given a row of data from SMOKA build a CAOM2 Energy object.

    :param bandpass_database: database that contains SVO filter bandpass information.
    :param row:  SMOKA text file Row
    :type row: Table.Row
    :type bandpass_database: svo.BandPassDatabase
    :return: CAOM2.Energy
    :rtype: Energy
    """

    # Build the energy object
    filter_name = None
    filter_names = []
    max_wavelength = 5.0 * units.um
    min_wavelength = 0.3 * units.um

    if "FILTER" in row.colnames:
        try:
            filter_name = row['FILTER'].split('-')[-1]
            if filter_name[-1] == "+":
                filter_name = "SDSS_{}".format(filter_name[0])
            if filter_name[0] == 'L':
                filter_name = "NB{}".format(filter_name[1:])

            filter_info = bandpass_database[filter_name]
            max_wavelength = min(max_wavelength, filter_info['wavelength_max'])
            min_wavelength = max(min_wavelength, filter_info['wavelength_min'])
            logging.debug("{}: {} {}".format(filter_name, filter_info['wavelength_min'], filter_info['wavelength_max']))
            logging.debug("min/max: {} {}".format(min_wavelength, max_wavelength))
        except:
            pass
    else:
        for filter_number in range(1, 4):
            if 'FILTER0{}'.format(filter_number) not in row.colnames:
                continue
            raw_filter_name = row['FILTER0{}'.format(filter_number)]
            if "NONE" in raw_filter_name or "open" in raw_filter_name:
                continue

            # The headers keywords have extra characters when compared to the SVO database.
            group = re.match('^SCFCFL[BS]?(.*?)(01)?$', raw_filter_name)
            if group is not None:
                filter_name = group.groups()[0]
            else:
                # No matching filter bits.
                continue
            filter_names.append(filter_name)
            filter_info = bandpass_database[filter_name]
            if filter_info is None:
                continue

            # There can be multiple filters, they are cumulative.
            max_wavelength = min(max_wavelength, filter_info['wavelength_max'])
            min_wavelength = max(min_wavelength, filter_info['wavelength_min'])
            logging.debug("{}: {} {}".format(filter_name, filter_info['wavelength_min'], filter_info['wavelength_max']))
            logging.debug("min/max: {} {}".format(min_wavelength, max_wavelength))

    cwl = (min_wavelength + max_wavelength) / 2.0
    bandwidth = (max_wavelength - min_wavelength)
    resolving_power = float(cwl / bandwidth)

    disperse = None
    if 'DISPERSR' in row.colnames:
        disperse = row['DISPERSR']
    if disperse is not None:
        grism_name = GRISM_NAMES.get(disperse, None)
        if grism_name is not None:
            grism = DISPERSERS.get(grism_name, None)
            if grism is not None:
                grism_params = None
                for filter_name in filter_names:
                    grism_params = grism.get(filter_name, None)
                    if grism_params is not None:
                        break
                if grism_params is None:
                    grism_params = grism.get("DEFAULT")
                if grism_params is not None:
                    min_wavelength = max(min_wavelength, grism_params[0] * units.AA)
                    max_wavelength = min(max_wavelength, grism_params[1] * units.AA)
                    resolving_power = grism_params[2]
                    bandwidth = max_wavelength - min_wavelength

    energy = Energy()

    try:
        energy.bounds = Interval(min_wavelength.to('m').value,
                                 max_wavelength.to('m').value,
                                 samples=[shape.SubInterval(min_wavelength.to('m').value,
                                                            max_wavelength.to('m').value)])
        energy.resolving_power = float(resolving_power)
        energy.sample_size = bandwidth.to('m').value
    except Exception as energy_exception:
        logging.warning("Failed to build Energy for observation {}".format(energy_exception))
        pass

    energy.em_band = EnergyBand.OPTICAL
    energy.bandpass_name = u'{}'.format(filter_name)
    energy.dimension = 1

    return energy


def data_product_type(obs_mod):
    """
    compute the DataProductType based on the SMOKA obs_mod keyword.

    :param obs_mod: value of obs_mod in the SMOKA obslog file
    :return: the data product type
    :rtype: DataProductType
    """
    data_type_map = dict(IMAG=DataProductType.IMAGE,
                         SPEC=DataProductType.SPECTRUM,
                         IMAG_SINGLE=DataProductType.IMAGE,
                         MOS=DataProductType.SPECTRUM)
    return data_type_map.get(obs_mod, None)


def build_observation(smoka_meta_data_row, instrument_name='SUP'):
    """
    Build a CAOM2 observation record based on the meta data from SMOKA.

    :param smoka_meta_data_row: the SMOKA datatable to build the CAOM2 record from.
    :type smoka_meta_data_row: Table.Row
    :param instrument_name: Instrument (short form) to retrieve observations records from SMOKA
    :type instrument_name: basestring
    :return: The CAOM2 Observation
    :rtype: Observation
    """

    row = smoka_meta_data_row
    this_telescope = Telescope(name=TELESCOPE_NAME.get(instrument_name, 'Subaru'))
    # These are actually the JCMT values.
    this_telescope.geo_location_x = -5461060.909
    this_telescope.geo_location_y = -2491393.621
    this_telescope.geo_location_z = 2149257.916
    # SMOKA provides a file with 'FRAME_ID' values but the ObservationID should be the Exposure ID so we change A to E
    observation_id = u'{}'.format(row['FRAME_ID'].replace("A", "E"))
    # If the trailing letter is a 'Y' or and 'X' this tells us if there are 8 or 10 CCDs with this Exposure.
    # Those are in the obslog file but in SMOKA they want '1' and '0'.  Change that here, to be consistent.
    case = observation_id[-1]
    if case == "Y":
        observation_id = observation_id[0:-1] + "1"
    if case == "X":
        observation_id = observation_id[0:-1] + "0"

    bandpass_database = svo.BandpassFilterDatabase(telescope=TELESCOPE_NAME.get(instrument_name, 'Subaru'),
                                                   instrument=INSTRUMENT_NAMES.get(instrument_name, None)[0])
    augment_filter_lookup_table(bandpass_database)
    # First lets define the Observation that is this record.
    this_instrument = Instrument(name=INSTRUMENT_NAMES[instrument_name][0])
    this_instrument.keywords.add(str(row['OBS_MOD']))

    obstype = u"{}".format(row['DATA_TYP'].upper())
    target_name = u'{}'.format(row['OBJECT2']).upper()
    if obstype == u'OBJECT':
        intent = ObservationIntentType.SCIENCE
    else:
        intent = ObservationIntentType.CALIBRATION

    disperse = None
    if "DISPERSR" in row.colnames:
        disperse = GRISM_NAMES.get(row["DISPERSR"], None)
        if disperse is not None:
            this_instrument.keywords.add(str(disperse))

    previous_observation = tap.tap_query("SUBARU", observation_id)['lastModified']
    print(previous_observation)
    if len(previous_observation) > 0 and previous_observation['lastModified'] > VERSION_DATE:
        logging.debug("Existing entry for {} is newer than caom2 script so don't replace.".format(observation_id))
        raise ValueError("Stored Record Newer Than Script.")
    print(observation_id)

    logging.info("Creating CAOM2 observation: {}".format(observation_id))

    this_observation = SimpleObservation(collection=this_telescope.name,
                                         observation_id=observation_id, sequence_number=None,
                                         intent=intent,
                                         type=obstype,
                                         proposal=None,
                                         telescope=this_telescope,
                                         instrument=this_instrument,
                                         target=Target(name=target_name),
                                         meta_release=time.Time('2017-01-01 00:00:00').to_datetime()
                                         )

    # Create a plane that will hold the raw data
    product_id = row['FRAME_ID']
    this_plane = Plane(u"{}".format(row['FRAME_ID']),
                       meta_release=time.Time('2017-01-01 00:00:00').to_datetime(),
                       data_release=time.Time('2017-01-01 00:00:00').to_datetime()
                       )
    this_plane.calibration_level = CalibrationLevel.RAW_STANDARD
    this_plane.data_product_type = data_product_type(row['OBS_MOD'])
    logging.debug("Disperser: {}".format(disperse))
    if disperse is None:
        this_plane.data_product_type = data_product_type("IMAG")

    this_plane.provenance = Provenance(name='SMOKA',
                                       producer='SMOKA',
                                       project='SMOKA',
                                       reference=PROVENANCE_URLS[instrument_name])

    # Build the time object.
    this_time = None
    start_time = None
    try:
        start_time = time.Time("{} {}".format(row['DATE_OBS'], row['UT_STR']))
        exptime = row['EXPTIME'] * units.second
        end_time = time.Time(start_time + exptime)
        time_bounds = Interval(start_time.mjd,
                               end_time.mjd,
                               samples=[shape.SubInterval(start_time.mjd, end_time.mjd)])
        this_time = Time(bounds=time_bounds,
                         dimension=1,
                         resolution=exptime.to('second').value,
                         sample_size=exptime.to('day').value,
                         exposure=exptime.to('second').value
                         )
    except Exception as time_exception:
        logging.error("Error building Time for {}".format(observation_id))
        logging.error("{}".format(time_exception))

    this_plane.time = this_time

    # Build the position object.
    position = None
    try:
        coord = SkyCoord(row['RA2000'], row['DEC2000'], unit=('hour', 'degree'))
        if "X" in product_id:
            pbounds = position_bounds(coord.ra.degree, coord.dec.degree, "{}X".format(instrument_name))
        elif "Y" in product_id:
            pbounds = position_bounds(coord.ra.degree, coord.dec.degree, "{}Y".format(instrument_name))
        else:
            pbounds = position_bounds(coord.ra.degree, coord.dec.degree, instrument_name)

        position = Position(bounds=pbounds,
                            sample_size=PIXEL_SCALE[instrument_name].to('arcsecond').value,
                            time_dependent=False)

    except Exception as time_exception:
        logging.error("Error building Position for {}".format(observation_id))
        logging.error("{}".format(time_exception))

    this_plane.position = position

    energy = None
    try:
        energy = build_energy(row, bandpass_database)
    except Exception as time_exception:
        logging.error("Error building Energy for {}".format(observation_id))
        logging.error("{}".format(time_exception))

    this_plane.energy = energy

    # create a reference to the data file stored at SMOKA.
    obsdate = start_time.iso[0:10]
    this_plane.artifacts.add(Artifact(uri='subaru:raw/{}/{}'.format(obsdate,
                                                                    this_plane.product_id),
                                      product_type=ProductType.SCIENCE,
                                      content_type='text/html',
                                      release_type=ReleaseType.DATA))

    # Add the PREVIEW artifact, stored in SMOKA, These are referenced using the FrameID value of one of the members.
    # but the frameID has an 'X' or 'Y' at end for SUP, so we need to adjust here.
    this_plane.artifacts.add(Artifact(uri='subaru:preview/{}'.format(
        this_plane.product_id.replace("X", "1'").replace("Y", "1")),
        product_type=ProductType.PREVIEW,
        release_type=ReleaseType.META,
        content_type='image/png'))

    # And the plane to the observation.
    this_observation.planes.add(this_plane)

    return this_observation


def caom2repo(this_observation, caom2repo_client):
    """
    Put an observation into the CAOM repo service

    :param this_observation: the CAOM2 Python object to store to caom2repo service
    :param caom2repo_client:  A CAOM@ Repo client instance that will PUT this_observation to a repository.
    :return:
    """

    try:
        logging.info('Inserting observation {}'.format(this_observation.observation_id))
        if caom2repo_client is not None:
            caom2repo_client.put_observation(this_observation)
        else:
            print(this_observation)
    except cadcutils_execptions.AlreadyExistsException as repo_execption:
        logging.debug(type(repo_execption))
        logging.info('Deleting observation {}'.format(this_observation.observation_id))
        caom2repo_client.delete_observation(this_observation.collection, this_observation.observation_id)
        logging.info('Inserting observation {}'.format(this_observation.observation_id))
        caom2repo_client.put_observation(this_observation)


def main(instrument_name='SUP', year='2002', caom2repo_client=None, frame_id=None):
    observation_table = get_metadata_table(instrument_name=instrument_name,
                                           year=year)
    if frame_id is not None:
        logging.debug("Only doing: {}".format(frame_id))
        observation_table = observation_table[observation_table['FRAME_ID'] == frame_id]
        logging.debug(str(observation_table))

    for row in observation_table:
        try:
            caom2repo(build_observation(row, instrument_name=instrument_name), caom2repo_client=caom2repo_client)
        except Exception as ex:
            logging.error("Failed to build repo record for row: {}".format(row))
            logging.error(str(ex))


def augment_filter_lookup_table(bandpass_database):
    """
    Add some filters to the internal database, based on information from Subaru website.

    :param bandpass_database:
    :return:
    """

    energy_bouds = dict(O58=(580 * units.nm, 1000 * units.nm),
                        Y47=(470 * units.nm, 910 * units.nm),
                        SDSS_z=(813.850 * units.nm, 1026.852 * units.nm),
                        L600=(370 * units.nm, 600 * units.nm),
                        L550=(340 * units.nm, 525 * units.nm),
                        C50=(500 * units.nm, 11000 * units.nm),
                        Hcont=((1.573 - 0.020 / 2) * units.um, (1.573 + 0.020 / 2.0) * units.um),
                        Brgamma=((2.166 - 0.032 / 2.) * units.um, (2.166 + 0.032 / 2) * units.um))

    for key in energy_bouds:
        bandpass_database.add_static_filter(key, energy_bouds[key])


def get_repo_client(x509_certificate, repo_resource_id):
    """
    Instantiate a caom2 repository client to put records to caom2 repository with.

    :param x509_certificate:
    :param repo_resource_id:
    :return: a CAOM2RepoClient instance
    :rtype: CAOM2RepoClient
    """
    try:
        # Create a CAOM2RepoClient object.
        return CAOM2RepoClient(net.Subject(certificate=x509_certificate), resource_id=repo_resource_id)
    except Exception as repo_exception:
        logging.error("Failed to create a repo client:{}".format(repo_exception))
        return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Retrieves metadata table from SMOKA and creates CAOM2 entries")
    parser.add_argument('instrument', choices=INSTRUMENT_NAMES.keys(),
                        help="Short name of SMOKA instrument archive to process.")
    parser.add_argument('year', choices=range(1998, 2018), type=int, help="Year of observation set to load.")
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--repo-resource-id', default='ivo://cadc.nrc.ca/sc2repo',
                        help="The CAOM2REPO service to post records to")
    parser.add_argument('--certfile', default=os.path.join(os.getenv('HOME'), '.ssl/cadcproxy.pem'),
                        help="CADC X509 proxy certificat to use for authentication")
    parser.add_argument('--dry-run', action="store_true", help="Do a Dry Run and print the caom2 record to stdout.")
    parser.add_argument('--frame-id', help="Just do this frame ID from this year.", default=None)

    args = parser.parse_args()

    log_level = logging.ERROR
    if args.verbose:
        log_level = logging.INFO
    elif args.debug:
        log_level = logging.DEBUG

    logging.basicConfig(level=log_level)

    repo_client = None
    if not args.dry_run:
        # Create a CAOM2RepoClient object.
        certificate = args.certfile
        resource_id = args.repo_resource_id
        repo_client = CAOM2RepoClient(net.Subject(certificate=certificate), resource_id=resource_id)

    main(args.instrument, args.year, caom2repo_client=repo_client, frame_id=args.frame_id)
