from __future__ import unicode_literals
import argparse

import cadcutils
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

DATADIR = __PATH__ = os.path.join(os.path.dirname(__file__), 'data')
SUP_HEADER = os.path.join(DATADIR, "SUP_Header.txt")

INSTRUMENT_NICK_NAMES = ('MIR', 'HSC', 'CIA', 'FCS', 'HDS', 'MCS', 'FMS', 'K3D',
                         'HIC', 'SUP', 'CAC', 'IRC', 'COM', 'OHS')
SMOKA_OBSLOG_ENDPOINT = 'http://smoka.nao.ac.jp/status/obslog'
TELESCOPE_NAME = {"SUP": 'SUBARU'}
INSTRUMENT_NAMES = {"MIR": ("MOICS", "Multi-Object Infrared Camera and Spectrograph"),
                    'HSC': ("HSC", "Hyper Suprime-Cam"),
                    "CIA": ("CIAO", "Coronagraphic Imager"),
                    "FCS": ("FOCUS", "Faint Object Camera And Spectrograph"),
                    "HDS": ("HDS", "High Dispersion Spectrograph"),
                    "MCS": ("MOIRCS", "Multi-Object Infrared Camera and Spectrograph"),
                    "FMS": ("FMOS", "Fiber Multi Object Spectrograph"),
                    "K3D": ("Kyote3DII", "Kyoto tridimensional spectrograph II"),
                    "HIC": ("HiCIAO", "High Contrast Instrument for the Subaru Next Generation Adaptive Optics"),
                    "SUP": ("Suprime", 'Suprime-Cam'),
                    "CAC": "CAC",
                    "COM": "Cooled Mid-Infrared Camera and Spectrograph",
                    "IRC": "Infrared Camera and Spectrograph",
                    "OHS": "Cooled Infrared Spectrograph and Camera for OHS",
                    }

PIXEL_SCALE = {"SUP": 0.2 * units.arcsecond}


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
    local_filename = "{}_{}.txt".format(instrument_name, year)
    if not os.access(local_filename, os.R_OK ):
        smoka_obslog = "{}/{}".format(SMOKA_OBSLOG_ENDPOINT, local_filename)
        name_temp_file = NamedTemporaryFile()
        local_filename = str(name_temp_file.name)
        with open(local_filename, str('w')) as fobj:
            try:
                resp = requests.get(smoka_obslog)
                resp.raise_for_status()
            except Exception as ex:
                logging.error("Failed to get metadata table at: {}".format(smoka_obslog))
                logging.error(str(ex))
                return None
            fobj.write(resp.content)
        name_temp_file.seek(0)
    return parse_meta_data_table(local_filename)


def parse_meta_data_table(local_filename):
    """
    Use the astropy package to return a Table representation of a SMOKA observation record.

    :param local_filename: name of text file with SMOKA observations in a text file.
    :type local_filename: str
    :return: Observation Record
    :rtype: Table
    """
    position_line = open(SUP_HEADER).readlines()[1].split(' ')
    col_ends = []
    col_starts = [0]
    for col in position_line:
        if col:
            col_ends.append(col_starts[-1] + len(col))
            col_starts.append(col_ends[-1] + 1)
        else:
            col_starts[-1] += 1
    col_starts = col_starts[:-1]

    colnames = open(local_filename).readline()[1:].split()
    return ascii.read(local_filename,
                      format="fixed_width_no_header",
                      col_ends=col_ends,
                      col_starts=col_starts,
                      names=colnames,
                      fill_values=('', '----', '---'))


def compute_fov(ra, dec, width=5*2048, height=2*4096):
    """

    :param ra:
    :param dec:
    :return:
    """



def position_bounds(ra, dec, width=34 / 60.0, height=27 / 60.0):
    """
    Given the RA and DEC of a SuprimeCam image return a Polygon bounding box.

    :param ra: RA of the centre of the FOV (degrees)
    :type ra: float
    :param dec: DEC of the centre of the FOV (degrees)
    :type dec: float
    :param width: width of FOV (degrees)
    :param height: height of FOV (degrees)
    :return: bounding box
    :rtype: CoordPolygon2D
    """
    naxis1 = width / PIXEL_SCALE["SUP"].to('degree').value
    naxis2 = height/ PIXEL_SCALE["SUP"].to('degree').value
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

    :param row:  SMOKA text file Row
    :return: CAOM2.Energy
    :rtype: Energy
    """

    raw_filter_name = row['FILTER']

    # Build the energy object
    energy = Energy()
    try:
        filter_name = row['FILTER'].split('-')[-1]
        if filter_name[-1] == "+":
            filter_name = "SDSS_{}".format(filter_name[0])
        if filter_name[0] == 'L':
            filter_name = "NB{}".format(filter_name[1:])

        filter_info = bandpass_database[filter_name]

        cwl = (filter_info['wavelength_min'] + filter_info['wavelength_max']) / 2.0
        bandwidth = (filter_info['wavelength_max'] - filter_info['wavelength_min'])
        energy.bounds = Interval(filter_info['wavelength_min'].to('m').value,
                                 filter_info['wavelength_max'].to('m').value,
                                 samples=[shape.SubInterval(filter_info['wavelength_min'].to('m').value,
                                                            filter_info['wavelength_max'].to('m').value)])
        energy.resolving_power = float(cwl / bandwidth)
        energy.sample_size = bandwidth.to('m').value
    except:
        pass

    energy.em_band = EnergyBand.OPTICAL
    energy.bandpass_name = u'{}'.format(raw_filter_name)
    energy.dimension = 1

    return energy


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

    bandpass_database = svo.BandpassFilterDatabase(telescope=TELESCOPE_NAME.get(instrument_name, 'Subaru'),
                                                   instrument=INSTRUMENT_NAMES.get(instrument_name, None)[0])

    # First lets define the Observation that is this record.
    this_instrument = Instrument(name=INSTRUMENT_NAMES[instrument_name][1])
    this_instrument.keywords.add(str(row['OBS_MOD']))

    obstype = u'{}'.format(row['DATA_TYP'].upper())
    target_name = u'{}'.format(row['OBJECT2']).upper()
    if obstype == u'OBJECT':
        intent = ObservationIntentType.SCIENCE
    else:
        intent = ObservationIntentType.CALIBRATION

    # SMOKA provides a file with 'FRAME_ID' values but the ObservationID should be the Exposure ID so we change A to E
    observation_id = u'{}'.format(row['FRAME_ID'].replace("A","E"))
    # If the trailing letter is a 'Y' or and 'X' this tells us if there are 8 or 10 CCDs with this Exposure.
    # Those are in the obslog file but in SMOKA they want '1' and '0'.  Change that here, to be consistent.
    case = observation_id[-1]
    if case == "Y":
        observation_id = observation_id[0:-1]+"1"
    if case == "X":
        observation_id = observation_id[0:-1]+"0"

    logging.info("Creating CAOM2 observation: {}".format(observation_id))

    this_observation = SimpleObservation(collection=this_telescope.name,
                                         observation_id=observation_id,
                                         algorithm=Algorithm(),
                                         sequence_number=None,
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
    # Here we again have a trailing 'Y' or an 'X' depening on the number of frames, but now we use the 'A' identifier.
    # Leave the 'Y' and 'X' so that they can be passed to the download request page.
    this_plane = Plane(u"{}".format(row['FRAME_ID']),
                       meta_release=time.Time('2017-01-01 00:00:00').to_datetime(),
                       data_release=time.Time('2017-01-01 00:00:00').to_datetime()
                       )
    this_plane.calibration_level = CalibrationLevel.RAW_STANDARD
    this_plane.data_product_type = DataProductType.IMAGE
    this_plane.provenance = Provenance(name='SMOKA',
                                       producer='SMOKA',
                                       project='SMOKA',
                                       reference=(
                                           'https://www.subarutelescope.org/Observing/Instruments/SCam/index.html'
                                       ))

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
    except Exception as ex:
        logging.error("Error building Time for {}".format(observation_id))
        logging.error("{}".format(ex))

    this_plane.time = this_time

    # Build the position object.
    position = None
    try:
        coord = SkyCoord(row['RA2000'], row['DEC2000'], unit=('hour', 'degree'))
        if "X" in product_id:
            width = 34 / 60.0
        else:
            width = 4 * 34 / (5 * 60.0)

        position = Position(bounds=position_bounds(coord.ra.degree, coord.dec.degree, width=width),
                            sample_size=PIXEL_SCALE[instrument_name].to('arcsecond').value,
                            time_dependent=False)
    except Exception as ex:
        logging.error("Error building Position for {}".format(observation_id))
        logging.error("{}".format(ex))

    this_plane.position = position

    energy = None
    try:
        energy = build_energy(row, bandpass_database)
    except Exception as ex:
        logging.error("Error building Energy for {}".format(observation_id))
        logging.error("{}".format(ex))

    this_plane.energy = energy

    # create a reference to the data file stored at SMOKA.
    obsdate = start_time.iso[0:10]
    this_plane.artifacts.add(Artifact(uri='subaru:raw/{}/{}'.format(obsdate,
                                                                    product_id),
                                      product_type=ProductType.SCIENCE,
                                      content_type='text/html',
                                      release_type=ReleaseType.DATA))

    # Add the PREVIEW artifact, stored in SMOKA, These are refenced using the FrameID value of one of the members.
    this_plane.artifacts.add(Artifact(uri='subaru:preview/{}1'.format(product_id[0:-1]),
                                      product_type=ProductType.PREVIEW,
                                      release_type=ReleaseType.META,
                                      content_type='image/png'))

    # And the plane to the observation.
    this_observation.planes.add(this_plane)

    return this_observation


def caom2repo(this_observation, repo_client):
    """
    Put an observation into the CAOM repo service

    :param this_observation: the CAOM2 Python object to store to caom2repo service
    :return:
    """

    try:
        logging.info('Inserting observation {}'.format(this_observation.observation_id))
        if repo_client is not None:
            repo_client.put_observation(this_observation)
        else:
            print(this_observation)
    except cadcutils.exceptions.AlreadyExistsException as ex:
        logging.debug(type(ex))
        logging.info('Deleting observation {}'.format(this_observation.observation_id))
        repo_client.delete_observation(this_observation.collection, this_observation.observation_id)
        logging.info('Inserting observation {}'.format(this_observation.observation_id))
        repo_client.put_observation(this_observation)


def main(instrument_name='SUP', year='2002', repo_client=None, frame_id=None):
    observation_table = get_metadata_table(instrument_name=instrument_name,
                                           year=year)
    if frame_id is not None:
        logging.debug("Only doing: {}".format(frame_id))
        observation_table  = observation_table[observation_table['FRAME_ID']==frame_id]
        logging.debug(str(observation_table))

    for row in observation_table:
        try:
            caom2repo(build_observation(row, instrument_name=instrument_name), repo_client=repo_client)
        except Exception as ex:
            logging.error("Failed to build repo record for row: {}".format(row))
            logging.error(str(ex))



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
        try:
            # Create a CAOM2RepoClient object.
            certificate = args.certfile
            resource_id = args.repo_resource_id
            repo_client = CAOM2RepoClient(net.Subject(certificate=certificate), resource_id=resource_id)
        except Exception as ex:
            logging.error("Failed to create a repo client:{}".format(ex))


    main(args.instrument, args.year, repo_client=repo_client, frame_id=args.frame_id)
