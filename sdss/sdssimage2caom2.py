from __future__ import unicode_literals
from builtins import str
from caom2 import *
from astropy.table import Table
from astropy.time import Time
from astropy import units

sdss_imaging_table = Table.read('test.fits', format='fits')

# define the wavelenght bounderies for SDSS filters (in nm):  (these are just made up, look up the correct ones.)
filter_bounds = {'u':  [348,375], 'g': [395,538], 'r':[550, 683], 'i':[673, 823], 'z': [850, 885]}

for row in sdss_imaging_table:

    # First lets define the Observation that is this record.
    instrument = Instrument(name=u'SDSS Camera')
    instrument.keywords.add(str(row['camcol']))

    observation = Observation(collection=u'SDSS',
                              observation_id=str(row['fieldID']),
                              algorithm=Algorithm(u'stack'),
                              sequence_number=None,
                              intent=ObservationIntentType.SCIENCE,
                              obs_type=u'OBJECT',
                              proposal=Proposal(proposal_id=u'SDSS',
                                                pi_name=u'SDSS',
                                                project=u'SDSS',
                                                title=u'SDSS'),
                              telescope=Telescope(name=u"Apache Point 2m",
                                                  geo_location_x=0.0,
                                                  geo_location_y=0.0,
                                                  geo_location_z=0.0,
                                                  keywords=None),
                              instrument=instrument,
                              target=Target(name=u'{}'.format(row['fieldID']),
                                            ),
                              meta_release=Time('2017-01-01 00:00:00').to_datetime()
                              )

    # Is Each filter from SDSS its own plane?  Likely.  So, planes are defined on the bandpass from CJW table.

    # Now start at the other end and define the Chunk.
    # There are not chunks by for now we must mock up chunks.  caom 2.3 will allow chunkless ingestion, future work.

    position_bounds = CoordPolygon2D(vertices=None)
    position_bounds.vertices.append(ValueCoord2D(row['raMin'], row['decMin']))
    position_bounds.vertices.append(ValueCoord2D(row['raMin'], row['decMax']))
    position_bounds.vertices.append(ValueCoord2D(row['raMax'], row['decMax']))
    position_bounds.vertices.append(ValueCoord2D(row['raMax'], row['decMin']))

    position = SpatialWCS(
        axis=CoordAxis2D(
            axis1=Axis("RA-TAN", "deg"),
            axis2=Axis("DEC-TAN", "deg"),
            bounds=position_bounds),
        coordsys="ICRS",
        equinox=2000.0)

    for bandpass_name in ['u', 'g', 'r', 'i', 'z']:


        # A quick stab at the energy bounds.
        energy_axis = wcs.CoordAxis1D(axis=wcs.Axis('WAVE', 'nm'))
        energy_axis.range = wcs.CoordRange1D(wcs.RefCoord(0.5, float(filter_bounds[bandpass_name][0])),
                                             wcs.RefCoord(1.5, float(filter_bounds[bandpass_name][1])))

        energy = SpectralWCS(axis=energy_axis,
                             specsys='TOPCENT',
                             ssysobs='TOPCENT',
                             ssyssrc='TOPCENT',
                             bandpass_name=bandpass_name)

        # Now the time bounds.
        mjd = row['tai_{}'.format(bandpass_name)]/(24*3600.0)
        start_time = Time(mjd, scale='tai', format='mjd')
        # Just some rough estimate of the drift scan time for SDSS imaging.
        end_time = Time(start_time + 143 * units.second)

        temporal_axis = wcs.CoordAxis1D(axis=wcs.Axis('UTC', 'd'))
        temporal_axis.range = wcs.CoordRange1D(wcs.RefCoord(0.5, start_time.mjd),
                                               wcs.RefCoord(0.5, end_time.mjd))
        time = TemporalWCS(axis=temporal_axis,)

        # Build a 'fake chunk' (wont need to do this in 2.3 apparently)
        chunk = Chunk(product_type=ProductType.SCIENCE,
                      position=position,
                      energy=energy,
                      time=time)


        # Add this chunk to a part.
        part = Part('image', product_type=ProductType.SCIENCE)
        part.chunks.append(chunk)


        # Stick the part in an artifact.
        URL_TO_THE_FILE = "http://blah.blah.blah/jfjfjf{}".format(bandpass_name)
        artifcat = Artifact(uri='sdss:....',
                            product_type=ProductType.SCIENCE,
                            release_type=ReleaseType.DATA)
        artifcat.parts.add(part)

        # add the artifact to the plane.
        plane = Plane(u"{}_{}".format(row['fieldID'],bandpass_name))
        plane.artifacts.add(artifcat)

        # And the plane to the observation.
        observation.planes.add(plane)
    print observation
    filename = '{}.xml'.format(row['fieldID'])
    with open(filename, str('w')) as fout:
        ObservationWriter().write(observation, fout)
