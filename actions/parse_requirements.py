

def parse_requirements(filename):
    """Parse the requirements and return (in conda notation)
       the requirements for this system
    """
    try:
        from pip_requirements_parser import RequirementsFile
    except ImportError as e:
        print("\n\n[ERROR] ** You need to install pip-requirements-parser")
        print("Run `conda install pip-requirements-parser\n\n")
        raise e

    from pkg_resources import evaluate_marker

    reqs = RequirementsFile.from_file(filename).to_dict()["requirements"]

    deps = {}

    for req in reqs:
        name = req["name"]
        specifier = req["specifier"]
        marker = req["marker"]

        if len(specifier) == 0:
            specifier = ""
        else:
            specifier = specifier[0]

        if marker is not None:
            # check to see if this line fits this platform
            include = evaluate_marker(marker)
        else:
            include = True

        if include:
            deps[name] = specifier

    reqs = list(deps.keys())
    reqs.sort()

    result = []

    for req in reqs:
        result.append(f"{req}{deps[req]}")

    return result
