from __future__ import unicode_literals
import argparse

from cadcutils import exceptions as cadcutils_execptions
from caom2 import get_acc_meta_checksum
from caom2 import get_differences
from cadcutils import net
from caom2 import *
from caom2.obs_reader_writer import ObservationWriter
import logging
import os
from caom2repo import CAOM2RepoClient
import observation
import tap
import smoka


def caom2ingest(instrument_name, year, caom2repo_client=None, frame_id=None, force=False, dry_run=False):
    """
    Given a SUBARU instrument name and year of operation, create a caom2.Observation and
    store that into the caom2repo via the provide client.

    If a caom2repo_client is not provided then the caom2 object is written to disk in XML format.

    :param instrument_name: Short form name of Subaru instrument (e.g. SUP)
    :type instrument_name: basestring
    :param year: Year of observations to ingest
    :type year: int
    :param caom2repo_client: client instance to use when putting record to caom2repository
    :type caom2repo_client: CAOM2RepoClient
    :param frame_id: only process this frame_id from the SMOKA obslog listing.
    :type frame_id: basestring
    :param force: Reprocess this observation even if its already in the CAOM database and newer than script.
    :type force: bool
    :return: None
    """

    observation_table = smoka.get_metadata_table(instrument_name=instrument_name, year=year)

    if frame_id is not None:
        logging.debug("Only doing: {}".format(frame_id))
        observation_table = observation_table[observation_table['FRAME_ID'] == frame_id]
        logging.debug("\n"+str(observation_table))

    for row in observation_table:
        this_observation = observation.SMOKA(row)
        try:
            previous_observation = caom2repo_client.read(this_observation.collection,
                                                             this_observation.observation_id)
            f = get_acc_meta_checksum(previous_observation)
            this_observation.previous_observation = previous_observation
            this_observation_caom2_Observation = this_observation()
            g = get_acc_meta_checksum(this_observation_caom2_Observation)
            if  f == g:
                logging.info("Matches previous observation.")
                if not force:
                    logging.info("Skipping")
                    continue
                logging.info("Overwriting previous observation.")
            else:
                logging.info("Updating previous record.")
        except cadcutils_execptions.NotFoundException:
            pass
        if not dry_run:
            try:
                logging.info('Inserting observation {}'.format(this_observation.observation_id))
                caom2repo_client.put_observation(this_observation())
            except cadcutils_execptions.AlreadyExistsException as repo_execption:
                logging.debug(type(repo_execption))
                logging.info('Deleting existing observation {}'.format(this_observation.observation_id))
                caom2repo_client.delete_observation(this_observation.collection, this_observation.observation_id)
                logging.info('Inserting new version of observation {}'.format(this_observation.observation_id))
                caom2repo_client.put_observation(this_observation())
        else:
            # If we were not given a caom2repo object then write the caom2 objects to disk as an xml document.
            with open('{}.xml'.format(this_observation.observation_id), 'w') as fobj:
                ObservationWriter().write(this_observation(), fobj)
            logging.info(str(this_observation))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Retrieves metadata table from SMOKA and creates CAOM2 entries")
    parser.add_argument('instrument', choices=observation.INSTRUMENT_NICK_NAMES,
                        help="Short name of SMOKA instrument archive to process.")
    parser.add_argument('year', choices=range(1998, 2018), type=int, help="Year of observation set to load.")
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--repo', default='sc2repo',
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

    if args.repo == "tardis":
        repo_resource_id = "ivo://cadc.nrc.ca/ams"
        tap.TAP_SERVER = "http://beta.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/ams/maq/sync"
    else:
        repo_resource_id = "ivo://cadc.nrc.ca/sc2repo"
        tap.TAP_SERVER = "http://sc2.canfar.net/sc2tap/sync"

    # Create a CAOM2RepoClient object.
    certificate = args.certfile
    repo_client = CAOM2RepoClient(net.Subject(certificate=certificate), resource_id=repo_resource_id)

    caom2ingest(args.instrument, args.year, caom2repo_client=repo_client,
                frame_id=args.frame_id, force=args.force, dry_run=args.dry_run)
