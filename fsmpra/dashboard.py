"""Bundle the local viewer and measured pilot data into a portable HTML file."""

import json
from pathlib import Path

from . import SCHEMA_VERSION
from .io import new_output, read_json


def export_dashboard(pilot_path, output):
    pilot = read_json(pilot_path)
    from .compare import COLLECTION_SCHEMA, check_compatible
    if pilot.get("schema_version") == COLLECTION_SCHEMA:
        if not pilot.get("datasets"):
            raise ValueError("Empty comparison collection")
        check_compatible(pilot["datasets"])
    elif pilot.get("schema_version") != SCHEMA_VERSION or not pilot.get("variants"):
        raise ValueError("Expected a nonempty pilot measurement bundle")
    assets = Path(__file__).resolve().parent.parent / "dashboard"
    template = (assets / "index.html").read_text()
    text = json.dumps(pilot, separators=(",", ":"), allow_nan=False)
    # JSON is placed in a script element; escape HTML delimiters even though
    # the script is inert JSON. Names/provenance can contain arbitrary text.
    text = text.replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
    template = template.replace("/* STYLE_PLACEHOLDER */", (assets / "style.css").read_text())
    template = template.replace("/* APP_PLACEHOLDER */", (assets / "app.js").read_text())
    template = template.replace('<script id="pilot-data" type="application/json">null</script>',
                                '<script id="pilot-data" type="application/json">' + text + '</script>')
    out = new_output(output)
    path = out / "index.html"
    path.write_text(template)
    return {"html": str(path), "bytes": path.stat().st_size, "mode": "local_self_contained"}
