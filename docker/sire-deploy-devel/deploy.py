
import os
import Sire.Base

par_url = os.environ["PAR_URL"]

if par_url is None:
    raise KeyError("You must supply the PAR URL!")

upload_cmd = "curl"
args = ("-v", "-X", "PUT" ,"-F", "'file=@./sire_devel_latest_linux.run'",
        "%s/sire_devel_latest_linux.run" % par_url)

process = Sire.Base.Process.run(upload_cmd, args)
process.wait()

if process.isError():
    raise ValueError("Something went wrong!")

