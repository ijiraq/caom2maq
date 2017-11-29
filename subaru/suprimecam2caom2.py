from __future__ import unicode_literals
import argparse
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
                    "SUP": ("Suprime", 'Subaru Prime Focus Camera'),
                    "CAC": "CAC",
                    "COM": "Cooled Mid-Infrared Camera and Spectrograph",
                    "IRC": "Infrared Camera and Spectrograph",
                    "OHS": "Cooled Infrared Spectrograph and Camera for OHS",
                    }
PIXEL_SCALE = {"SUP": 0.2*units.arcsecond}


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

    smoka_obslog = "{}/{}_{}.txt".format(SMOKA_OBSLOG_ENDPOINT, instrument_name, year)
    local_filename = NamedTemporaryFile()
    resp = requests.get(smoka_obslog)
    resp.raise_for_status()
    local_filename.write(resp.content)
    local_filename.seek(0)
    return parse_meta_data_table(local_filename.name)

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
    points = []
    vertices = []

    # A polygon needs a set of points that define the corners and a set of vectors that define
    # the area inside those points that is the covered area.
    segment_type = SegmentType['MOVE']
    for x, y in ([-0.5, -0.5], [-0.5, 0.5], [0.5, 0.5], [0.5, -0.5]):
        xx = ra + x * width
        yy = dec + y * height
        points.append(Point(xx, yy))
        vertices.append(Vertex(xx, yy, segment_type))
        segment_type = SegmentType['LINE']

    # Close up the sample area
    vertices.append(Vertex(ra - 0.5 * width,
                           dec - 0.5 * height,
                           SegmentType['CLOSE']))

    return Polygon(points=points, samples=shape.MultiPolygon(vertices))


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

    observation_id = u'{}'.format(row['FRAME_ID'].replace("A", "E").replace("X", "0"))
    logging.info("Creating CAOM2 observation: {}".format(observation_id))

    this_observation = SimpleObservation(collection=this_telescope.name,
                                         observation_id=observation_id,
                                         algorithm=Algorithm(u'simple'),
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
    this_plane = Plane(u"{}".format(row['FRAME_ID'].replace("X", "0")), meta_release=time.Time('2017-01-01 00:00:00').to_datetime())
    this_plane.calibration_level = CalibrationLevel.RAW_STANDARD
    this_plane.data_product_type = DataProductType.IMAGE
    this_plane.provenance = Provenance(name='SMOKA',
                                       producer='SMOKA',
                                       project='SMOKA',
                                       reference=(
                                           'https://www.subarutelescope.org/Observing/Instruments/SCam/index.html'
                                       ))

    # Build the time object.
    start_time = time.Time("{} {}".format(row['DATE_OBS'], row['UT_STR']))
    exptime = row['EXPTIME'] * units.second
    end_time = time.Time(start_time + exptime)
    time_bounds = Interval(start_time.mjd,
                           end_time.mjd,
                           samples=[shape.SubInterval(start_time.mjd, end_time.mjd)])
    this_plane.time = Time(bounds=time_bounds,
                           dimension=1,
                           resolution=exptime.to('second').value,
                           sample_size=exptime.to('day').value,
                           exposure=exptime.to('second').value
                           )

    # Build the position object.
    try:
        coord = SkyCoord(row['RA2000'], row['DEC2000'], unit=('hour', 'degree'))
        position = Position(bounds=position_bounds(coord.ra.degree, coord.dec.degree),
                            sample_size=PIXEL_SCALE[instrument_name].to('arcsecond').value,
                            time_dependent=False)
    except Exception as ex:
        logging.error("{}".format(ex))
        position = None
    this_plane.position = position

    # Build the energy object
    filter_name = row['FILTER'].split('-')[-1]
    if filter_name[-1] == "+":
        filter_name = "SDSS_{}".format(filter_name[0])
    if filter_name[0] == 'L':
        filter_name = "NB{}".format(filter_name[1:])

    filter_info = bandpass_database[filter_name]

    cwl = (filter_info['wavelength_min']+filter_info['wavelength_max'])/2.0
    bandwidth = (filter_info['wavelength_max']-filter_info['wavelength_min'])
    energy = Energy()
    energy.bounds = Interval(filter_info['wavelength_min'].to('m').value,
                             filter_info['wavelength_max'].to('m').value,
                             samples=[shape.SubInterval(filter_info['wavelength_min'].to('m').value,
                                                        filter_info['wavelength_max'].to('m').value)])
    energy.bandpass_name = u'{}'.format(filter_info['bandpass_name'])
    energy.dimension = 1
    energy.em_band = EnergyBand.OPTICAL
    energy.resolving_power = float(cwl/bandwidth)
    energy.sample_size = bandwidth.to('m').value

    this_plane.energy = energy

    # create a reference to the data file stored at SMOKA.
    this_plane.artifacts.add(Artifact(uri='smoka:file/{}'.format(this_plane.product_id),
                                      product_type=ProductType.SCIENCE,
                                      release_type=ReleaseType.DATA))

    # Add the PREVIEW artifact, stored in SMOKA
    this_plane.artifacts.add(Artifact(uri='smoka:preview/{}'.format(this_plane.product_id),
                                      product_type=ProductType.PREVIEW,
                                      release_type=ReleaseType.META,
                                      content_type='image/png'))

    # And the plane to the observation.
    this_observation.planes.add(this_plane)

    return this_observation


def caom2repo(this_observation, resource_id='ivo://cadc.nrc.ca/sc2repo'):
    """
    Put an observation into the CAOM repo service

    :param this_observation: the CAOM2 Python object to store to caom2repo service
    :param resource_id: the caom2 service to put records to
    :return:
    """

    certificate = os.path.join(os.getenv('HOME'), '.ssl/cadcproxy.pem')

    repo_client = CAOM2RepoClient(net.Subject(certificate=certificate), resource_id=resource_id)

    try:
        logging.info('Inserting observation {}'.format(this_observation.observation_id))
        repo_client.put_observation(this_observation)
    except Exception as ex:
        logging.warning(str(ex))
        logging.info('Deleting observation {}'.format(this_observation.observation_id))
        repo_client.delete_observation(this_observation.collection, this_observation.observation_id)
        logging.info('Inserting observation {}'.format(this_observation.observation_id))
        repo_client.put_observation(this_observation)


def main(instrument_name='SUP', year='2002'):

    observation_table = get_metadata_table(instrument_name=instrument_name,
                                           year=year)

    for row in observation_table:
        caom2repo(build_observation(row, instrument_name=instrument_name))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Retrieves metadata table from SMOKA and creates CAOM2 entries")
    parser.add_argument('instrument', choices=INSTRUMENT_NAMES.keys(),
                        help="Short name of SMOKA instrument archive to process.")
    parser.add_argument('year', choices=range(2002, 2018), type=int, help="Year of observation set to load.")
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    log_level = logging.ERROR
    if args.verbose:
        log_level = logging.INFO
    elif args.debug:
        log_level = logging.DEBUG

    logging.basicConfig(level=log_level)

    main(args.instrument, args.year)
