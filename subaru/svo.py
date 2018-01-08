"""
Lookup filter information using the SVO filter service.
"""
from astropy.io import votable
from cStringIO import StringIO
import fcntl
import pickle
import logging
import os
import requests

DATADIR = __PATH__ = os.path.join(os.path.dirname(__file__), 'data')

class BandpassFilterDatabase(dict):
    """
    A database of filter information.

    This class is used to look up filter information relevant to CAOM2 from a cache file (a json object). If the filter
    is not known in that database then the SVO Filter Profile Service is queried for the required information and
    the result is used to augment the filter_database.
    """

    FILTER_SERVICE_URL = "http://svo2.cab.inta-csic.es/theory/fps/fps.php"

    def __init__(self, telescope, instrument, cache_filename=None, *args, **kwargs):

        self.telescope = telescope
        self.instrument = instrument
        self._cache_filename = cache_filename
        self._filter_cache = {}
        self.init_filter_cache()
        super(BandpassFilterDatabase, self).__init__(*args, **kwargs)

    @property
    def cache_filename(self):
        """
        Name of the filter cache file.
        :return:
        """
        if self._cache_filename is None:
            self._cache_filename = "filter_cache_{}_{}.pkl".format(self.telescope, self.instrument)
            self._cache_filename = os.path.join(DATADIR, self._cache_filename)
        return self._cache_filename

    def init_filter_cache(self):
        """
        Load the filter_cache from the json file.
        :return:
        """
        if not os.access(str(self.cache_filename), os.R_OK):
            file(str(self.cache_filename), str('wb')).close()
        cfd = file(str(self.cache_filename), str('rb'))
        try:
            fcntl.lockf(cfd.fileno(), fcntl.LOCK_SH)
            p = pickle.Unpickler(cfd)
            self._filter_cache = p.load()
        except EOFError:
            self._filter_cache = {}
            pass
        finally:
            fcntl.lockf(cfd.fileno(), fcntl.LOCK_UN)
            cfd.close()

    def update_filter_cache(self):
        """
        Write the current filter_cache database to a json file.
        :return:
        """
        cfd = file(str(self.cache_filename), str('wb+'))
        p = pickle.Pickler(cfd)
        try:
            fcntl.lockf(cfd.fileno(), fcntl.LOCK_EX)
            p.dump(self._filter_cache)
        finally:
            fcntl.lockf(cfd.fileno(), fcntl.LOCK_UN)
            cfd.close()

    def __getitem__(self, k):
        """

        :param k:
        :return: filter info
        :rtype: dict
        """
        if k not in self._filter_cache.keys():
            v = self.look_up_filter(k)
            self._filter_cache[k] = v
            self.update_filter_cache()
        return self._filter_cache[k]

    def look_up_filter(self, filter_name):

        params = {'ID': "{}/{}.{}".format(self.telescope, self.instrument, filter_name)}

        try:
            votable_file = StringIO(requests.get(BandpassFilterDatabase.FILTER_SERVICE_URL,
                                                 params=params).content)
            votable_file.seek(0)
            table = votable.parse_single_table(votable_file)
        except Exception as ex:

            logging.warning(str(ex))
            return None
        wavelength_min = table.get_field_by_id('WavelengthMin')
        wavelength_max = table.get_field_by_id('WavelengthMax')
        values = {'wavelength_max': wavelength_max.value * wavelength_max.unit,
                  'wavelength_min': wavelength_min.value * wavelength_min.unit,
                  'bandpass_name': table.get_field_by_id('filterID').value}
        return values

    def add_static_filter(self, filter_name, bounds):
        """
        Add a description of a filter to the database.

        :param filter_name: Name of the filter to be used
        :param bounds: A trouble with the lowerer and upper energy bouds
        :type bounds: (Quantity, Quantity)
        :return: None
        """
        self._filter_cache[k] = {'wavelength_min': bounds[0],
                                 'wavelength_max': bounds[1],
                                 'bandpasss_name': filter_name}

