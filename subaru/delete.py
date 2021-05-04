import tap
from caom2repo import CAOM2RepoClient
from cadcutils import net
import numpy 
import sys

certfile="/Users/kavelaarsj/.ssl/cadcproxy.pem"
resource_id="ivo://cadc.nrc.ca/ams"

client = CAOM2RepoClient(net.Subject(certificate=certfile), resource_id=resource_id)

tap.TAP_SERVER = "http://beta.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/ams/maq/sync"

t = tap.tap_query('SUBARU', '{}%'.foramt(sys.argv[1]), product_id=sys.argv[2]) 

for oid in numpy.unique(t['observationID']):
   client.delete_observation("SUBARU", oid)

