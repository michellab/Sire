
import os
import pycurl
import sys
from io import BytesIO

par_url = os.environ["PAR_URL"]

if par_url is None:
    raise KeyError("You must supply the PAR URL!")

try:
    upload_file = sys.argv[1]
    key = upload_file.split("/")[-1]
except:
    upload_file = "/home/sireuser/sire_devel_latest_linux.run"
    key = "sire_devel_latest_linux.run"

if not os.path.exists(upload_file):
    raise ValueError("Cannot find %s" % upload_file)

data = open(upload_file, "rb").read()

if par_url.endswith("/"):
    destination = "%s%s" % (par_url, key)
else:
    destination = "%s/%s" % (par_url, key)

buffer = BytesIO()
c = pycurl.Curl()
c.setopt(c.URL, destination)
c.setopt(c.WRITEDATA, buffer)
c.setopt(c.CUSTOMREQUEST, "PUT")
c.setopt(c.POST, True)
c.setopt(c.POSTFIELDS, data)

c.perform()
c.close()
