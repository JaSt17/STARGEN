from branca.element import MacroElement
from folium.elements import JSCSSMixin
from branca.element import Template, MacroElement


class BigImage(JSCSSMixin, MacroElement):
    """
    Adds a BigImage button to your map to enable image download functionality.

    Parameters
    ----------
    position : str
        Change the position of the button, can be:
        'topleft', 'topright', 'bottomright' or 'bottomleft'
        Default: 'topleft'.
    filename : str
        The name of the downloaded image file.
        Default: 'map'.
    download_link : bool
        Whether to show the download link.
        Default: True.
    hide_control : bool
        Whether to hide the control button after download.
        Default: False.
    custom_options : dict
        Additional options for the BigImage control.
    """

    _template = Template(
        """
        {% macro script(this, kwargs) %}
            L.control.bigImage(
                {{ this.options|tojson }}
            ).addTo({{ this._parent.get_name() }});
        {% endmacro %}
        """
    )

    default_js = [
        (
            "Control.BigImage.js",
            "https://cdn.jsdelivr.net/gh/pasichnykvasyl/Leaflet.BigImage/dist/Leaflet.BigImage.min.js",
        )
    ]
    default_css = [
        (
            "Control.BigImage.css",
            "https://cdn.jsdelivr.net/gh/pasichnykvasyl/Leaflet.BigImage/dist/Leaflet.BigImage.min.css",
        )
    ]
    
    def __init__(
        self,
        position="topleft",
        inputTitle="Choose a scale",
        title="Download image",
        minScale=1,
        **kwargs
    ):
        super().__init__()
        self._name = "BigImage"
        self.options = {
            "position": position,
            "inputTitle": inputTitle,
            "title": title,
            "minScale": minScale,
            "maxScale": minScale+5,
        }

        self.options.update(kwargs)

