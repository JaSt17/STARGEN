import h3
import folium
from folium import Map, Element
from branca.element import Template, MacroElement
import matplotlib.colors as mcolors
import base64
from folium.plugins import AntPath
import matplotlib.colors as mcolors


def split_hexagon_if_needed(hexagon):
    """
    Splits a hexagon if it crosses the antimeridian.

    A hexagon crosses the antimeridian if the difference between its maximum 
    and minimum longitudes is greater than 180 degrees. This function checks 
    for such a condition and splits the hexagon into two parts if necessary.

    Parameters:
        hexagon (str): The H3 index of the hexagon to be checked and potentially split.

    Returns:
        list: A list containing one or two tuples of coordinates. If the hexagon 
              does not cross the antimeridian, the list contains one tuple of 
              coordinates. If it does cross the antimeridian, the list contains 
              two tuples of coordinates representing the split hexagon.
    """
    # Get the boundary as a list of latitude-longitude pairs
    boundary = h3.cell_to_boundary(hexagon)
    longitudes = [lon for lat, lon in boundary]

    # Check if the hexagon crosses the antimeridian
    if max(longitudes) - min(longitudes) > 180:
        first_hex = []
        second_hex = []

        # Split the hexagon into two parts
        for lat, lon in boundary:
            if lon <= 0:  # Western hemisphere
                first_hex.append((lat, lon + 360))  # Adjust longitude for continuity
                second_hex.append((lat, lon))
            else:  # Eastern hemisphere
                first_hex.append((lat, lon))
                second_hex.append((lat, lon - 360))  # Adjust longitude for continuity

        return [tuple(first_hex), tuple(second_hex)]
    else:
        return [tuple(boundary)]


def draw_sample_hexagons(hex_dict, samples_per_hexagon, m=None, color='grey', zoom_start=1, show_samples_per_hexagon=True):
    """
    Draws hexagons on a map, displaying only the borders for hexagons that contain samples.

    Parameters:
        hex_dict (dict): A dictionary where keys are hexagon H3 indices and 
                         values are internal average sample distances.
        samples_per_hexagon (dict): A dictionary where keys are hexagon H3 indices and 
                                    values are the number of samples within the hexagon.
        m (folium.Map, optional): An existing Folium map object to plot on. 
                                  If None, a new map is created. Defaults to None.
        color (str, optional): The color of the hexagon borders. Defaults to 'grey'.
        zoom_start (int, optional): The initial zoom level of the map. Defaults to 1.

    Returns:
        folium.Map: The map object with the plotted hexagons.
    """
    if m is None:
        m = folium.Map(location=(0.0, 0.0), tiles="Esri worldstreetmap", zoom_start=zoom_start)

    for hexagon, sample_distance in hex_dict.items():
        parts = split_hexagon_if_needed(hexagon)
        for part in parts:
            polygon = folium.Polygon(
                locations=part,
                weight=1,
                color=color,
                fill_opacity=0.0,
                fill=True
            )
            polygon.add_child(folium.Tooltip(f"Internal scaled genetic distance: {sample_distance}"))
            polygon.add_to(m)

            if show_samples_per_hexagon and hexagon in samples_per_hexagon:
                # Calculate the center of the polygon
                latitudes = [point[0] for point in part]
                longitudes = [point[1] for point in part]
                center_lat = sum(latitudes) / len(latitudes)
                center_lon = sum(longitudes) / len(longitudes)

                # Add a marker at the center with the number of samples
                folium.Marker(
                    location=(center_lat, center_lon),
                    icon=folium.DivIcon(html=f'<div style="font-size: 12px; color: grey;">{samples_per_hexagon[hexagon]}</div>')
                ).add_to(m)

    return m


def draw_hexagons(hexagons, m=None, color='white', zoom_start=1, value=None, opacity=0.3, imputed=False):
    """
    Draws hexagons on a map.

    Parameters:
        hexagons (list): A list of hexagon H3 indices to be plotted.
        m (folium.Map, optional): An existing Folium map object to plot on. 
                                  If None, a new map is created. Defaults to None.
        color (str, optional): The fill color of the hexagons. Defaults to 'white'.
        zoom_start (int, optional): The initial zoom level of the map. Defaults to 1.
        value (str, optional): The value to display in the tooltip. Defaults to None.
        opacity (float, optional): The fill opacity of the hexagons. Defaults to 0.5.
        imputed (bool, optional): Whether the value is imputed. Adds "(Imputed)" 
                                  to the tooltip if True. Defaults to False.

    Returns:
        folium.Map: The map object with the plotted hexagons.
    """
    if m is None:
        m = folium.Map(location=(0.0, 0.0), tiles="Esri worldstreetmap", zoom_start=zoom_start)

    for hexagon in hexagons:
        parts = split_hexagon_if_needed(hexagon)
        for part in parts:
            polygon = folium.Polygon(
                locations=part,
                weight=0,
                color=None,
                fill_color=color,
                fill_opacity=opacity,
                fill=True
            )
            # Add imputed to tooltip if `imputed` is True
            tooltip_text = f"{value} (Imputed)" if imputed else str(value)
            polygon.add_child(folium.Tooltip(tooltip_text))
            polygon.add_to(m)

    return m


from matplotlib import cm, colors as mcolors

def get_color_gradient():
    """
    Returns a colormap for gradient coloring.
    """
    return cm.get_cmap('coolwarm')  # Example colormap

def draw_hexagons_with_values(hex_dict, m=None, zoom_start=1, threshold=0.0, imputed=False, opacity=0.5):
    """
    Draws hexagons on a map with values determining their fill color.

    Parameters:
        hex_dict (dict): A dictionary where keys are hexagon H3 indices and 
                         values are the distance values determining color and tooltip.
        m (folium.Map, optional): An existing Folium map object to plot on. 
                                  If None, a new map is created. Defaults to None.
        zoom_start (int, optional): The initial zoom level of the map. Defaults to 1.
        threshold (float, optional): The minimum value required to plot a hexagon. Defaults to 0.0.
        imputed (bool, optional): Whether the values are imputed. Adds "(Imputed)" 
                                  to the tooltip if True. Defaults to False.
        opacity (float, optional): The fill opacity of the hexagons. Defaults to 0.5.

    Returns:
        folium.Map: The map object with the plotted hexagons.
    """
    hexagons = hex_dict.keys()
    values = hex_dict.values()

    cmap = get_color_gradient()

    for hexagon, value in zip(hexagons, values):
        if value < threshold:
            continue

        # Normalize the value for the colormap (assuming values are between -1 and 1)
        normalized_value = (value + 1) / 2
        color = mcolors.to_hex(cmap(normalized_value))  # Convert to a hex color

        m = draw_hexagons(
            [hexagon],
            m,
            color=color,
            zoom_start=zoom_start,
            value=value,
            opacity=opacity,
            imputed=imputed
        )

    return m


def draw_barriers(barriers_dict, m=None, zoom_start=1, threshold=0.0):
    """
    Draws barriers on a map with colors based on their values.

    Parameters:
        barriers_dict (dict): A dictionary where keys are barrier coordinates 
                              (list of tuples) and values are the distance values 
                              for determining color and tooltip.
        m (folium.Map, optional): An existing Folium map object to plot on. 
                                  If None, a new map is created. Defaults to None.
        zoom_start (int, optional): The initial zoom level of the map. Defaults to 1.
        threshold (float, optional): The minimum value required to plot a barrier. Defaults to 0.0.

    Returns:
        folium.Map: The map object with the plotted barriers.
    """
    if m is None:
        m = folium.Map(location=(0.0, 0.0), tiles="Esri worldstreetmap", zoom_start=zoom_start)

    cmap = get_color_gradient()

    for barrier, value in barriers_dict.items():
        if value < threshold:
            continue

        # Normalize value for colormap
        normalized_value = (value + 1) / 2  # Adjust based on expected value range
        color = mcolors.to_hex(cmap(normalized_value))

        try:
            # Create a PolyLine for the barrier
            polyline = folium.PolyLine(
                locations=barrier,
                color=color,
                weight=3,  # Line thickness
                opacity=0.7  # Line transparency
            )
            # Add a tooltip showing the value
            tooltip_text = f"Value: {value:.2f}"
            polyline.add_child(folium.Tooltip(tooltip_text))
            polyline.add_to(m)
        except Exception as e:
            print(f"Error drawing barrier {barrier}: {e}")

    return m


def draw_migration_for_time_bin(time_bin, m, color="green"):
    """
    Draw migration paths for hexagon pairs within a specified time bin on a given map.

    Parameters:
        time_bin (dict): A dictionary where keys are tuples of hexagon H3 indices (hex1, hex2)
                         and values are the migration distances between them.
        m (folium.Map): An existing Folium map object to add the migration paths to.
        color (str, optional): The color of the migration paths. Defaults to "green".

    Returns:
        folium.Map: The map object with the migration paths added.
    """

    def adjust_for_antimeridian(midpoint1, midpoint2):
        """Adjusts midpoints for the antimeridian crossing."""
        if midpoint1[1] < midpoint2[1]:
            midpoint1_adj = (midpoint1[0], midpoint1[1] + 360)
            midpoint2_adj = (midpoint2[0], midpoint2[1] - 360)
        else:
            midpoint1_adj = (midpoint1[0], midpoint1[1] - 360)
            midpoint2_adj = (midpoint2[0], midpoint2[1] + 360)
        return [[midpoint1_adj, midpoint2], [midpoint1, midpoint2_adj]]

    for pair, distance in time_bin.items():
        hex1, hex2 = pair
        midpoint1 = h3.cell_to_latlng(hex1)  
        midpoint2 = h3.cell_to_latlng(hex2)  

        # Handle antimeridian crossing
        if abs(midpoint1[1] - midpoint2[1]) > 180:
            lines = adjust_for_antimeridian(midpoint1, midpoint2)
        else:
            lines = [[midpoint1, midpoint2]]
        
        # Add paths to the map
        for line in lines:
            ant_path = AntPath(
                locations=line,
                color=color,
                reverse=True,
                dash_array=[10, 20],  # Dashed path
                delay=800  # Animation delay
            )
            ant_path.add_child(folium.Tooltip(f'{distance} (Migration Distance)'))
            ant_path.add_to(m)

    return m


def get_color_gradient():
    """
    Create a custom colormap with a gradient of colors ranging from sand yellow to orange to dark red.

    Returns:
    cmap : LinearSegmentedColormap
        A matplotlib colormap object with the specified color gradient.
    """
    
    # Define the colors for the colormap
    colors = [
        (0.0, 0.93, 0.79, 0.69),  # Sand yellow
        (0.5, 1.0, 0.65, 0.0),    # Orange
        (1.0, 0.55, 0.0, 0.0)     # Dark red
    ]

    # Normalize colors to be between 0 and 1
    normalized_colors = [(value, (r, g, b)) for value, r, g, b in colors]
    
    # Create the colormap
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_color_gradient", normalized_colors, N=256)
    
    return cmap


def add_legend(m):
    """
    Adds a draggable legend to the provided folium map.

    The legend includes:
    - Symbols representing different types of areas and routes.
    - A color gradient representing scaled genetic distances.

    Parameters:
    m (folium.Map, optional): An existing Folium map object to plot on. 

    Returns:
    m : The map object with the legend added.
    """
    template = """
    {% macro html(this, kwargs) %}
    <div id='maplegend' class='maplegend' 
        style='position: absolute; z-index: 9999; background-color: rgba(255, 255, 255, 0.5);
        border-radius: 6px; padding: 10px; font-size: 10.5px; width: 180px; height: 110px; right: 20px; top: 20px; cursor: move;'>     
    <div class='legend-scale'>
    <ul class='legend-labels'>
        <li><svg height="12" width="12">
            <polygon points="5,0 10,3.33 10,8.67 5,12 0,8.67 0,3.33" style="fill:none;opacity: 0.5;stroke:black" />
            </svg>Area with Samples</li>
        <li><svg height="12" width="12">
            <polygon points="5,0 10,3.33 10,8.67 5,12 0,8.67 0,3.33" style="fill:black;opacity: 0.6;stroke:none" />
            </svg>Isolated Population</li>
        <li><svg height="12" width="10"><line x1="0" y1="2" x2="10" y2="10" style="stroke:green;stroke-width:2" /></svg>Possible Migration Route</li>
    </ul>
    </div>
    <div class='legend-gradient'>
        <span style="font-weight: bold;">Scaled Genetic Distances (log2)</span>
        <span style='background: linear-gradient(to right, 
            rgb(237, 201, 175) 0%,     /* Sand yellow */
            rgb(255, 165, 0) 50%,      /* Orange */
            rgb(139, 0, 0) 100%        /* Dark red */
        );
        width: 100%; height: 10px; display: block;'></span>
        <div style='display: flex; justify-content: space-between;'>
            <span>-1</span>
            <span>0</span>
            <span>1</span>
        </div>
    </div>
    </div> 
    <style type='text/css'>
    .maplegend .legend-scale ul {margin: 0; padding: 0; color: #0f0f0f;}
    .maplegend .legend-scale ul li {list-style: none; line-height: 18px; margin-bottom: 1.5px;}
    .maplegend ul.legend-labels li span {float: left; height: 12px; width: 12px; margin-right: 4.5px;}
    .maplegend ul.legend-labels li svg {margin-right: 4.5px;}
    </style>
    <script type='text/javascript'>
        dragElement(document.getElementById('maplegend'));

        function dragElement(element) {
            var pos1 = 0, pos2 = 0, pos3 = 0, pos4 = 0;
            if (document.getElementById(element.id + "header")) {
        
                document.getElementById(element.id + "header").onmousedown = dragMouseDown;
            } else {
            
                element.onmousedown = dragMouseDown;
            }

            function dragMouseDown(e) {
                e = e || window.event;
                e.preventDefault();
                pos3 = e.clientX;
                pos4 = e.clientY;
                document.onmouseup = closeDragElement;
                document.onmousemove = elementDrag;
            }

            function elementDrag(e) {
                e = e || window.event;
                e.preventDefault();
                pos1 = pos3 - e.clientX;
                pos2 = pos4 - e.clientY;
                pos3 = e.clientX;
                pos4 = e.clientY;
                element.style.top = (element.offsetTop - pos2) + "px";
                element.style.left = (element.offsetLeft - pos1) + "px";
            }

            function closeDragElement() {
                // stop moving when mouse button is released:
                document.onmouseup = null;
                document.onmousemove = null;
            }
        }
    </script>
    {% endmacro %}
    """
    macro = MacroElement()
    macro._template = Template(template)

    macro.add_to(m)
    return m
    
    
    