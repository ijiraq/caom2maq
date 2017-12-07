from __future__ import unicode_literals
import argparse
from cadcutils import exceptions
import numpy
from cadcutils import net
from caom2 import *
from astropy import time
from astropy.io import fits
from astropy import wcs
import logging
import os
from caom2repo import CAOM2RepoClient

# Create a CAOM2RepoClient object.
certificate = os.path.join(os.getenv('HOME'), '.ssl/cadcproxy.pem')
resource_id = 'ivo://cadc.nrc.ca/sc2repo'
repo_client = CAOM2RepoClient(net.Subject(certificate=certificate), resource_id=resource_id)

def nan2None(x):
    return not numpy.isnan(x) and x or None

def build_observation(header, override, continuum):
    """

    :type override: dict
    """


    w = wcs.WCS(header)


    spatial_wcs = SpatialWCS(axis=CoordAxis2D(axis1=Axis(ctype=unicode(w.wcs.ctype[0]), cunit=unicode(w.wcs.cunit[0])),
                                           axis2=Axis(ctype=unicode(w.wcs.ctype[1]), cunit=unicode(w.wcs.cunit[1])),
                                           function=CoordFunction2D(dimension=Dimension2D(naxis1=header['NAXIS1'],
                                                                                          naxis2=header['NAXIS2']),
                                                                    ref_coord=Coord2D(RefCoord(w.wcs.crpix[0],
                                                                                               w.wcs.crval[0]),
                                                                                      RefCoord(w.wcs.crpix[1],
                                                                                               w.wcs.crval[1])),
                                                                    cd11=w.wcs.cdelt[0],
                                                                    cd22=w.wcs.cdelt[1],
                                                                    cd12=0.0,
                                                                    cd21=0.0
                                                                    )
                                           ),
                          coordsys=unicode(w.wcs.radesys),
                          equinox=w.wcs.equinox,
                          resolution=1.97)


    spectral_wcs = SpectralWCS(axis=CoordAxis1D(axis=Axis(ctype=unicode(w.wcs.ctype[2]), cunit=unicode(w.wcs.cunit[2])),
                                          function=CoordFunction1D(naxis=header['NAXIS3'], delta=w.wcs.cdelt[2],
                                                                   ref_coord=RefCoord(w.wcs.crpix[2],
                                                                                      w.wcs.crval[2]),
                                                                   ),
                                          ),
                         specsys=unicode(w.wcs.specsys),
                         ssysobs=unicode(w.wcs.ssysobs),
                         ssyssrc=unicode(w.wcs.ssyssrc),
                         restfrq=nan2None(w.wcs.restfrq),
                         restwav=nan2None(w.wcs.restwav),
                         velosys=nan2None(w.wcs.velosys),
                         zsource=None,
                         bandpass_name=override['bandpass_name']
                         )

    telescope = Telescope(name="ALMA-12m",
                          geo_location_x=2225142.18,
                          geo_location_y=-5440307.37,
                          geo_location_z=-2481029.852
                          )

    instrument = Instrument(name="Band 3")

    target = Target(name=unicode(header['OBJECT']),
                    standard=False,
                    moving=False,
                    target_type=TargetType.OBJECT,
                    )

    proposal = Proposal(id=override['proposal_id'],
                        pi_name=override['proposal_pi'],
                        title=override['proposal_title'])

    proposal.keywords = override['proposal_keywords']

    observation = SimpleObservation(collection="ALMA",
                                    observation_id=override['observation_id'],
                                    algorithm=Algorithm(u'simple'),
                                    sequence_number=None,
                                    intent=ObservationIntentType.SCIENCE,
                                    type="OBJECT",
                                    proposal=proposal,
                                    telescope=telescope,
                                    instrument=instrument,
                                    target=target,
                                    meta_release=override['release_date']
                                    )

    provenance = Provenance(name="CASA",
                            version="CASA 4.2.2 (prerelease r30986)",
                            last_executed=time.Time(header["DATE"]).datetime,
                            reference="https://casa.nrao.edu/")

    start_mjd = time.Time(header['DATE-OBS']).mjd
    end_mjd = start_mjd + override['exptime']/24.0/3600.0
    time_bounds = Interval(start_mjd,
                           end_mjd,
                           samples=[shape.SubInterval(start_mjd, end_mjd)])

    this_time = Time(bounds=time_bounds,
                     dimension=1,
                     resolution=override['exptime'],
                     sample_size=override['exptime']/3600.0/24.0,
                     exposure=override['exptime']
                     )

    polarization_wcs = PolarizationWCS(axis=CoordAxis1D(axis=Axis(ctype=unicode(w.wcs.ctype[3]),
                                                              cunit=unicode(w.wcs.cunit[3]))))

    chunk = Chunk(product_type=ProductType.SCIENCE,
                  naxis=w.naxis,
                  position_axis_1=1,
                  position_axis_2=2,
                  position=spatial_wcs,
                  energy_axis=3,
                  energy=spectral_wcs,
                  polarization_axis=4,
                  polarization=polarization_wcs,
                  )


    part = Part(name='0', product_type=ProductType.SCIENCE)
    part.chunks.append(chunk)

    artifact = Artifact(uri=override['science.artifact_uri'],
                        product_type=ProductType.SCIENCE,
                        release_type=ReleaseType.DATA,
                        content_type=override['science.content_type'],
                        content_length=override['science.content_size'])


    artifact.parts.add(part)
    plane = Plane(product_id="source_calibrated_line_image",
                  data_release=override['release_date'],
                  meta_release=override['release_date'],
                  provenance=provenance)

    points = []
    vertices = []
    segment_type = SegmentType['MOVE']

    for x, y in ([0, 0], [1, 0], [1, 1], [0, 1 ]):
        V = w.all_pix2world(x*header['naxis1'], y*header['naxis2'], 1, 1, 1)
        points.append(Point(float(V[0]), float(V[1])))
        vertices.append(Vertex(float(V[0]), float(V[1]), segment_type))
        segment_type = SegmentType['LINE']
    vertices.append(Vertex(points[0].cval1,
                           points[0].cval2,
                           SegmentType['CLOSE']))

    polygon = shape.Polygon(points=points, samples=shape.MultiPolygon(vertices))
    position = Position(time_dependent=False,
                        bounds=polygon)
    plane.position = position

    polarization = Polarization(1,[PolarizationState.I,])
    plane.polarization = polarization

    from astropy import units, constants
    energy_upper_limit = ((header['NAXIS3'] - w.wcs.crpix[2]) * w.wcs.cdelt[2] + w.wcs.crval[2]) * w.wcs.cunit[2]
    energy_upper_limit = constants.c / energy_upper_limit
    energy_lower_limit = ((0 - w.wcs.crpix[2]) * w.wcs.cdelt[2] + w.wcs.crval[2]) * w.wcs.cunit[2]
    energy_lower_limit = constants.c / energy_lower_limit

    energy = Energy(bounds=shape.Interval(energy_lower_limit.value, energy_upper_limit.value,
                                          samples=[shape.SubInterval(energy_lower_limit.value,
                                                                     energy_upper_limit.value)]))
    plane.energy = energy
    plane.time = this_time


    plane.artifacts.add(artifact)
    # Add the PREVIEW artifact, stored in SMOKA
    plane.artifacts.add(Artifact(uri=override['preview.artifact_uri'],
                                 product_type=ProductType.PREVIEW,
                                 release_type=ReleaseType.META,
                                 content_type=override['preview.content_type']))

    plane.data_product_type = DataProductType.CUBE
    plane.calibration_level = CalibrationLevel.CALIBRATED
    observation.planes.add(plane)

    w = wcs.WCS(continuum)


    spatial_wcs = SpatialWCS(axis=CoordAxis2D(axis1=Axis(ctype=unicode(w.wcs.ctype[0]), cunit=unicode(w.wcs.cunit[0])),
                                           axis2=Axis(ctype=unicode(w.wcs.ctype[1]), cunit=unicode(w.wcs.cunit[1])),
                                           function=CoordFunction2D(dimension=Dimension2D(naxis1=header['NAXIS1'],
                                                                                          naxis2=header['NAXIS2']),
                                                                    ref_coord=Coord2D(RefCoord(w.wcs.crpix[0],
                                                                                               w.wcs.crval[0]),
                                                                                      RefCoord(w.wcs.crpix[1],
                                                                                               w.wcs.crval[1])),
                                                                    cd11=w.wcs.cdelt[0],
                                                                    cd22=w.wcs.cdelt[1],
                                                                    cd12=0.0,
                                                                    cd21=0.0
                                                                    )
                                           ),
                          coordsys=unicode(w.wcs.radesys),
                          equinox=w.wcs.equinox,
                          resolution=1.97)

    chunk = Chunk(product_type=ProductType.SCIENCE,
                  naxis=w.naxis,
                  position_axis_1=1,
                  position_axis_2=2,
                  position=spatial_wcs
                  )

    part = Part(name='0', product_type=ProductType.SCIENCE)
    part.chunks.append(chunk)

    artifact = Artifact(uri=override['science.continuum_artifact_uri'],
                        product_type=ProductType.SCIENCE,
                        release_type=ReleaseType.DATA,
                        content_type=override['science.content_type'],
                        content_length=override['science.continuum_content_size'])


    artifact.parts.add(part)
    plane = Plane(product_id="calibrated_final_cont_image",
                  data_release=override['release_date'],
                  meta_release=override['release_date'],
                  provenance=provenance)

    plane.energy = energy
    plane.position = position
    plane.polarization = polarization
    plane.time = this_time
    plane.artifacts.add(artifact)

    plane.artifacts.add(Artifact(uri=override['preview.continuum_artifact_uri'],
                                 product_type=ProductType.PREVIEW,
                                 release_type=ReleaseType.META,
                                 content_type=override['preview.content_type']))

    plane.data_product_type = DataProductType.IMAGE
    plane.calibration_level = CalibrationLevel.PRODUCT

    observation.planes.add(plane)


    # And the plane to the observation.
    observation.planes.add(plane)
    return observation


def caom2repo(this_observation):
    """
    Put an observation into the CAOM repo service

    :param this_observation: the CAOM2 Python object to store to caom2repo service
    :return:
    """

    try:
        logging.info('Inserting observation {}'.format(this_observation.observation_id))
        repo_client.put_observation(this_observation)
    except exceptions.AlreadyExistsException as ex:
        logging.info('Deleting observation {}'.format(this_observation.observation_id))
        repo_client.delete_observation(this_observation.collection, this_observation.observation_id)
        logging.info('Inserting observation {}'.format(this_observation.observation_id))
        repo_client.put_observation(this_observation)


def main(filename):
    header = fits.open(filename[0])[0].header
    continuum = fits.open(filename[1])[0].header


    override = dict(proposal_id="2013.1.00187.S",
                    proposal_title="A Survey of Ophiuchus and Orion B North",
                    bandpass_name="ALMA band 3",
                    proposal_keywords={'Low-mass star formation', 'pre-stellar cores', 'infra-red dark clouds (IRDC)'},
                    proposal_pi="Scott Schnee",
                    observation_id="A002_X9ab617_X1efb",
                    content_size=os.stat("source_calibrated_line_image_162608-24202.image.fits").st_size,
                    exptime=60.48
                    )

    override['science.artifact_uri'] = "vos://cadc.nrc.ca!vospace/helenkirk/ALMA_fits_files/{}/{}".format(override['proposal_id'], filename[0])
    override['science.continuum_artifact_uri'] = "vos://cadc.nrc.ca!vospace/helenkirk/ALMA_fits_files/{}/{}".format(override['proposal_id'], filename[1])
    override['science.continuum_content_size'] = os.stat(filename[1]).st_size
    override['science.content_type'] = 'application/fits'
    override['science.content_size'] = os.stat(filename[0]).st_size
    override['preview.artifact_uri'] = "vos://cadc.nrc.ca!vospace/helenkirk/ALMA_fits_files/{}/{}".format(override['proposal_id'],
                                                                                    filename[0].replace(".fits", ".gif"))
    override['preview.continuum_artifact_uri'] = "vos://cadc.nrc.ca!vospace/helenkirk/ALMA_fits_files/{}/{}".format(override['proposal_id'],
                                                                                    filename[1].replace(".fits", ".gif"))
    override['preview.content_type'] = 'image/png'
    override['release_date'] = time.Time(header['DATE-OBS'], scale='utc').datetime


    observation = build_observation(header, override, continuum)

    with open(str('junk.xml'), str('w')) as fobj:
        ObservationWriter().write(observation, fobj)
    caom2repo(observation)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Retrieves metadata table from SMOKA and creates CAOM2 entries")
    parser.add_argument('filename', nargs='+',
                        help="Alma fits file to process")
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    log_level = logging.ERROR
    if args.verbose:
        log_level = logging.INFO
    elif args.debug:
        log_level = logging.DEBUG

    logging.basicConfig(level=log_level)

    main(args.filename)
