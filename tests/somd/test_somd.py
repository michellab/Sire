import shlex
import subprocess
import tempfile

import sire as sr

def test_parameters():
    """A test to catch invalid SOMD parameters. Updates to add PME
       functionality have broken paramter resolution, meaning that OpenMM
       is configured to use the CUDA platform, even when the CPU platform
       is selected. This only happens for certain input parameters. Those
       used below can reproduce the issue.
    """

    # Create a temporary working directrory.
    tmp_dir = tempfile.TemporaryDirectory()

    # Load the test system.
    system = sr.load(
        ["https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/amber/ala/ala.crd",
         "https://raw.githubusercontent.com/michellab/BioSimSpace/devel/demo/amber/ala/ala.top"],
        directory=tmp_dir.name
    )

    # Write the SOMD configuration file.
    with open(f"{tmp_dir.name}/test.cfg", "w") as f:
        f.write("save coordinates = True\n"
                "minimise = True\n"
                "minimise maximum iterations = 100\n"
                "minimise tolerance = 1\n"
                "ncycles = 1\n"
                "nmoves = 1\n"
                "reaction field dielectric = 78.3\n"
                "cutoff type = cutoffperiodic\n"
                "cutoff distance = 8 angstrom\n"
                "barostat = False\n"
        )

    # Create the SOMD command.
    cmd = "somd -c ala.crd -t ala.top -C test.cfg -p CPU"

    # Run the subprocess.
    proc = subprocess.run(shlex.split(cmd), cwd=tmp_dir.name)

    # Make sure the process completed without error
    assert proc.returncode == 0
