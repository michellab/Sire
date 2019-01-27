
import os
import Sire.Base

par_url = os.environ["PAR_URL"]

if par_url is None:
    raise KeyError("You must supply the PAR URL!")

upload_cmd = "curl"

upload_file = "/home/sireuser/sire_devel_latest_linux.run"

if not os.path.exists(upload_file):
    raise ValueError("Cannot find %s" % upload_file)

args = ("-v", "-X", "PUT" ,"-F", "'file=@%s'" % upload_file,
        "%s/sire_devel_latest_linux.run" % par_url)

process = Sire.Base.Process.run(upload_cmd, args)
process.wait()

if process.isError():
    raise ValueError("Something went wrong!")
