from branca.element import MacroElement
from folium.elements import JSCSSMixin
from branca.element import Template

class EasyPrint(JSCSSMixin, MacroElement):
    """
    Adds an EasyPrint button to your map for downloading printable versions.

    Parameters
    ----------
    position : str
        Change the position of the button, can be:
        'topleft', 'topright', 'bottomright' or 'bottomleft'.
        Default: 'topleft'.
    title : str
        The title of the print button.
        Default: 'Print map'.
    filename : str
        The name of the output file.
        Default: 'map'.
    exportOnly : bool
        If true, hides the button on the map and triggers printing programmatically.
        Default: False.
    sizeModes : list
        List of size modes to enable. Options are 'Current', 'A4Landscape', 'A4Portrait'.
        Default: ['Current'].
    """

    _template = Template(
        """
        {% macro script(this, kwargs) %}
            L.easyPrint(
                {{ this.options|tojson }}
            ).addTo({{ this._parent.get_name() }});
        {% endmacro %}
        """
    )

    default_js = [
        (
            "leaflet.easyPrint.js",
            "https://cdn.jsdelivr.net/npm/leaflet-easyprint/dist/bundle.min.js",
        )
    ]
    default_css = []  # No CSS required for this plugin

    def __init__(
        self,
        position="topleft",
        title="Print map",
        filename="STARGEN_snapshot",
        exportOnly=True,
        sizeModes=["A4Landscape"],
        **kwargs
    ):
        super().__init__()
        self._name = "EasyPrint"
        self.options = {
            "position": position,
            "title": title,
            "filename": filename,
            "exportOnly": exportOnly,
            "sizeModes": sizeModes,
        }

        self.options.update(kwargs)
