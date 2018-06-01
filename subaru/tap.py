"""Query the CADC CAOM2 TAP service to determine the list of caom2 entries.

Logic:  Return all exposures associated with a collection.
"""
import requests
try:
    import cStringIO
except ImportError:
    import io as cStringIO
from astropy.io import votable
from astropy.time import Time
TAP_SERVER = "http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/tap/sync"

def tap_query(collection, observation_id, product_id=None, artifact_uri=None):
    """The __main__ part of the script"""

    tap_url = TAP_SERVER

    adql = ("SELECT observationID, artifact.lastModified as lastModified FROM caom2.Observation as observation "
            "JOIN caom2.Plane as plane on plane.obsID = observation.obsID "
            "JOIN caom2.Artifact as artifact on artifact.planeID = plane.planeID "
            "WHERE collection LIKE '{}' AND observationID LIKE '{}' ")

    adql = adql.format(collection, observation_id)

    if product_id is not None:
        adql += " AND plane.productID LIKE '{}' ".format(product_id)
    if artifact_uri is not None:
        adql += " AND artifact.uri LIKE '{}' ".format(artifact_uri)

    adql += "ORDER BY  artifact.lastModified DESC"

    # Some default parameters for that TAP service queries.
    tap_params = {'REQUEST': 'doQuery',
                  'LANG': 'ADQL',
                  'FORMAT': 'votable',
                  'QUERY': adql}

    result = requests.get(tap_url, params=tap_params)
    result.raise_for_status()
    temporary_file = cStringIO.StringIO(result.content)
    temporary_file.seek(0)
    result_table = votable.parse_single_table(temporary_file).to_table()
    if len(result_table) > 0:
        result_table['lastModified'] = Time(result_table['lastModified'])
    return result_table
