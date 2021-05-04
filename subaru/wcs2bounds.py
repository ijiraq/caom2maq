from caom2 import shape, SegmentType, Point, Polygon, Vertex
from astropy import wcs
import logging
__author__ = 'jjk'


def wcs2bounds(wcs_header, axes):
    """
    Given an astropy wcs return a CAOM2 bounding box.

    :param wcs_header: a FITS header that a WCS object can be made from.
    :param axes: a tuple giving the size of the image in x/y pixels

    :rtype Polygon
    :return a polygon object that describes the boundary of the FOV.
    """
    w = wcs.WCS(wcs_header)
    try:
        stc = w.calc_footprint(axes=axes)
    except Exception as ex:
        logging.error(str(ex))
        return None
    points = []
    vertices = []

    segment_type = SegmentType['MOVE']
    for xy in stc:
        xx = xy[0]
        yy = xy[1]
        points.append(Point(xx, yy))
        vertices.append(Vertex(xx, yy, segment_type))
        segment_type = SegmentType['LINE']

    # Close up the sample area
    vertices.append(Vertex(stc[0][0],
                           stc[0][1],
                           SegmentType['CLOSE']))

    return Polygon(points=points, samples=shape.MultiPolygon(vertices))
