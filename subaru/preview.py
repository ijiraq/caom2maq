import requests
import numpy as np
from PIL import Image
from tempfile import NamedTemporaryFile
import logging
from astropy.time import Time
import cadcdata
from cadcutils import net

endpoint = "http://smoka.nao.ac.jp/thumbnail"


def get_suprimecam_mosaic_preview(date="2002-05-07", frame_id="SUPA0010598X"):
    """
    Retrieve a PNG mosiac preview from the SMOKA site
    :param date:
    :param frame_id:
    :return:
    """
    endpoint = "http://smoka.nao.ac.jp/qlis/ImagePNG"

    frame_id = frame_id.replace("X", "0")
    params = {'frameid': frame_id,
              'dateobs': date,
              'grayscale': 'linear',
              'mosaic': 'true'}

    fobj = NamedTemporaryFile()
    fobj.write(requests.get(endpoint, params=params).content)
    fobj.seek(0)
    preview = Image.open(fobj.name)
    preview = preview.rotate(90, expand=True)
    fobj.close()

    expid = frame_id.replace("SUPA", "SUPE")

    preview_filename = "{}_preview.png".format(expid)
    preview.save(preview_filename)

    thumbnail_filename = "{}_thumbnail.png".format(expid)
    preview.thumbnail((256, 256))
    preview.save(thumbnail_filename)

    return {'preview': preview_filename, 'thumbnail': thumbnail_filename}


def build_suprimecam_peviews(date="2002-05-07", frame_id="SUPA0010598X", instrument="SUP"):
    """
    Build a Suprime Cam preview given the date and frame_id number.

    :param frame_id:
    :param date:
    :param instrument:
    :return:
    """

    ham_transition_date = Time('2008-08-01')
    frame = frame_id.rstrip("X")

    mit_layout = [[0, 1, 5, 4, 9],
                  [6, 7, 2, 3, 8]]

    ham_layout = [[6, 7, 2, 1, 0],
                  [8, 9, 5, 4, 3]]

    if Time(date) < ham_transition_date:
        ccd_layout = mit_layout
    else:
        ccd_layout = ham_layout

    column = []
    for i in range(2):
        row = []
        for k in range(len(ccd_layout[i])):
            this_url = "{}/{}/{}/{}{}.jpg".format(endpoint, instrument, date, frame, ccd_layout[i][k])
            logging.debug("Getting preview from {}".format(this_url))
            print this_url
            fobj = NamedTemporaryFile()
            fobj.write(requests.get(this_url).content)
            fobj.seek(0)
            row.append(Image.open(fobj.name))
        column.append(np.hstack(row))
    expid = frame.replace("SUPA", "SUPE")
    preview = Image.fromarray(np.vstack(column)).rotate(90, expand=True)
    preview.save("{}0_preview.jpg".format(expid))
    preview.thumbnail((256, 256))
    preview.save("{}0_thumbnail.jpg".format(expid))

if __name__ == '__main__':
    import sys
    data_client = cadcdata.CadcDataClient(net.Subject())

    # Don't call our custom preview builder any more :-(
    # build_suprimecam_peviews(sys.argv[1], sys.argv[2])

    for filename in get_suprimecam_mosaic_preview(sys.argv[1], sys.argv[2]).keys():
        data_client.put_file('subaru', filename)
