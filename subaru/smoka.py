import requests
import logging
import os
from tempfile import NamedTemporaryFile
from astropy.io import ascii
import observation

DATADIR = __PATH__ = os.path.join(os.path.dirname(__file__), 'data')
SMOKA_OBSLOG_ENDPOINT = 'http://smoka.nao.ac.jp/status/obslog'


def get_metadata_table(instrument_name, year):
    """
    retrieve the SMOKA metadata table from VOSpace.

    :param year: Year of observation records to retrieve from SMOKA
    :type year: int
    :param instrument_name: Instrument (short form) to retrieve observations records from SMOKA
    :type instrument_name: basestring
    :rtype: Table
    :return: table of metadata values
    """
    if instrument_name not in observation.INSTRUMENT_NICK_NAMES:
        raise ValueError("Unkown instrument_name", instrument_name)

    logging.info("Building metadata table for : {} {}".format(instrument_name, year))
    local_filename = os.path.join(DATADIR, "{}_{}.txt".format(instrument_name, year))

    if not os.access(local_filename, os.R_OK):
        smoka_obslog = "{}/{}_{}.txt".format(SMOKA_OBSLOG_ENDPOINT, instrument_name, year)
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
