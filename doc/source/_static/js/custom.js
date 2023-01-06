function choose_version(value) {
    console.log("CHOOSE VERSION");
    console.log(value);

    console.log(window.location.href);
    var parts = window.location.href.split("sire.openbiosim.org/");
    console.log(parts);

    // should be two parts - before sire.openbiosim.org/ and after
    s = parts[0];
    s += "sire.openbiosim.org";
    s += value;

    if (parts[1].startsWith("versions/")) {
        var parts2 = parts[1].split("/");
        parts2.shift();
        parts2.shift();
        parts[1] = parts2.join("/");
    }

    s += parts[1];

    console.log(`Redirecting to ${s}`);
    window.location = s;
}

function fill_versions() {
    var xmlhttp = new XMLHttpRequest();

    href = window.location.href;

    console.log(href);

    xmlhttp.onreadystatechange = function () {
        if (this.readyState == 4 && this.status == 200) {
            var versions = JSON.parse(this.responseText);
            console.log(versions)

            var s = "<div class=\"version_box\">" +
                "<select id=\"versions\" onchange=\"choose_version(this.value)\">";

            for (var i in versions) {
                var version = versions[i][0];
                var path = versions[i][1];
                if (href.includes(path)) {
                    // we are looking at this version now
                    s += `<option value="${path}" selected>${version}</option>`;
                } else {
                    s += `<option value="${path}">${version}</option>`;
                }
            }

            s += "</select></div>";

            document.getElementById("version_box").innerHTML = s;
        }
        else
        {
            console.log("INCORRECT RESPONSE!");
            console.log(this.status);

            var s = "<div class=\"version_box\">" +
                    "<select id=\"versions\" onchange=\"choose_version(this.value)\">";

            s += `<option value="${path}" selected>main</option>`;
        }

        s += "</select></div>";

        document.getElementById("version_box").innerHTML = s;
    };

    try
    {
        xmlhttp.open("GET", "https://sire.openbiosim.org/versions.json",
                     true);
        xmlhttp.send();
    }
    catch(error)
    {
        console.log("FAILED TO GET VERSIONS!");
        console.log(error);
    }
}

window.onload = fill_versions;
