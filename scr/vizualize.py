from h3 import h3
import folium
from folium import Map, Element
from func import normalize_distances
from branca.element import Template, MacroElement
import matplotlib.colors as mcolors
import base64
from folium.plugins import AntPath


# function that splits a hexagon if it crosses the antimeridian
def split_hexagon_if_needed(hexagon):
    boundary = h3.h3_to_geo_boundary(hexagon, geo_json=False)
    longitudes = [lon for lat, lon in boundary]

    # Check if the hexagon crosses the antimeridian
    if max(longitudes) - min(longitudes) > 180:
        boundary = list(boundary)
        first_hex = list()
        second_hex = list()
        # make get two hexagons from the original one
        for i in range(len(boundary)):
            if boundary[i][1] <= 0:
                first_hex.append((boundary[i][0], boundary[i][1] + 360))
                second_hex.append((boundary[i][0], boundary[i][1]))
            if boundary[i][1] > 0:
                first_hex.append((boundary[i][0], boundary[i][1]))
                second_hex.append((boundary[i][0], boundary[i][1] - 360))
        # return two tuples of coordinates of the two hexagons
        return [tuple(first_hex), tuple(second_hex)]
    else:
        # return the original hexagon
        return [boundary]
    
# function that empty hexagons on a map with only the boarders for the hexagons which hold sampels
def draw_sample_hexagons(hex_dict, m=None, color='grey', zoom_start=1):
    # Create a map if it is not provided
    if m is None:
        m = folium.Map(location=(0.0, 0.0), tiles="Esri worldstreetmap", zoom_start=zoom_start)

    # Plot hexagons
    for hexagon in hex_dict.keys():
        # split the hexagon if it crosses the antimeridian
        parts = split_hexagon_if_needed(hexagon)
        for part in parts:
            polygon = folium.Polygon(
                locations=part,
                weight=1,
                color=color,
                fill_opacity=0.0,
                fill=True
            )
            # add tooltip with internal hexagon distance
            polygon.add_child(folium.Tooltip(f"Internal average sample distance: {hex_dict[hexagon]}"))
            polygon.add_to(m)
    return m


# function that draws hexagons on a map
def draw_hexagons(hexagons, m=None, color='darkgreen', zoom_start=1, value=None, opacity=0.5, imputed=False):
    # Create a map if it is not provided
    if m is None:
        m = folium.Map(location=(0.0, 0.0), tiles="Esri worldstreetmap", zoom_start=zoom_start)

    # Plot hexagons
    for hexagon in hexagons:
        # split the hexagon if it crosses the antimeridian
        parts = split_hexagon_if_needed(hexagon)
        for part in parts:
            polygon = folium.Polygon(
                locations=part,
                weight=1,
                color=None,
                fill_color=color,
                fill_opacity=opacity,
                fill=True
            )
            if value >= 0:
                if imputed:
                    polygon.add_child(folium.Tooltip(f"{value} (Imputed)"))
                else:
                    polygon.add_child(folium.Tooltip(value))
            
            polygon.add_to(m)
    return m

def draw_hexagons_with_values(hex_dict, m=None, zoom_start=1, threshold=0.0, imputed=False):
    hexagons = hex_dict.keys()
    values = hex_dict.values()
    
    # get color gradient
    cmap = get_color_gradient()

    # write the values to the center of each hexagon
    for hexagon, value in zip(hexagons, values):
        if value < threshold:
            continue
        col = mcolors.to_hex(cmap(value))
        m = draw_hexagons([hexagon], m, color=col, zoom_start=zoom_start, value=value, imputed=imputed)
    return m


def draw_barriers(barriers_dict, m=None, zoom_start=1, threshold=0.0):
    barriers = barriers_dict.keys()
    values = barriers_dict.values()
    
    # get color gradient
    cmap = get_color_gradient()
    
    for barrier, value in zip(barriers, values):
        if value < threshold:
            continue
        col = mcolors.to_hex(cmap(value))
        barrier = list(barrier)
        # try to draw the barrier
        try:
            polyline = folium.PolyLine(barrier, color = col)
            polyline.add_child(folium.Tooltip(value))
            polyline.add_to(m)
        except:
            pass
    return m

# this function takes a time bin and a map and draws all neighboring lines for the hexagons in that time bin
def draw_migration_for_time_bin(time_bin, m, color="green"):
    # Loop through all pairs of hexagons in the time bin
    for pair, distance in time_bin.items():
        hex1, hex2 = pair
        
        # Get the midpoints of both hexagons
        midpoint1 = h3.h3_to_geo(hex1)
        midpoint2 = h3.h3_to_geo(hex2)
        
        # Check if the points are on opposite sides of the antimeridian
        if abs(midpoint1[1] - midpoint2[1]) > 180:
            if midpoint1[1] < midpoint2[1]:
                midpoint1_adj = (midpoint1[0], midpoint1[1] + 360)
                midpoint2_adj = (midpoint2[0], midpoint2[1] - 360)
            else:
                midpoint1_adj = (midpoint1[0], midpoint1[1] - 360)
                midpoint2_adj = (midpoint2[0], midpoint2[1] + 360)
            lines = [[midpoint1_adj, midpoint2], [midpoint1, midpoint2_adj]]
        else:
            lines = [[midpoint1, midpoint2]]
        
        # Loop over all lines and draw them on the map
        for line in lines:
            # draw ant paths which show the migration routes
            ant_path = AntPath(locations=line, color='green', reverse='true', dash_array=[10, 20], delay=800)
            ant_path.add_child(folium.Tooltip(f'{distance} (Migration Distance)'))
            ant_path.add_to(m)
    
    return m

def get_color_gradient():
    # Define the colors for the colormap
    colors = [
        (0.95, 0.87, 0.55),  # Sand color yellow
        (0.92, 0.78, 0.48),  # Intermediate yellow-brown
        (0.87, 0.69, 0.38),  # Yellow-brown tone
        (0.78, 0.62, 0.34),  # Intermediate brown
        (0.66, 0.72, 0.73),  # Light blue-gray
        (0.57, 0.76, 0.85),  # Intermediate light blue
        (0.48, 0.82, 0.95),  # Light blue
        (0.42, 0.76, 0.91),  # Intermediate turquoise
        (0.36, 0.70, 0.95)   # Turquoise blue
    ]

    # Create the colormap
    cmap = mcolors.LinearSegmentedColormap.from_list("costum_color_gradient", colors)
    
    return cmap

def add_legend(m):
    # Create the legend template as an HTML element
    template = """
    {% macro html(this, kwargs) %}
    <div id='maplegend' class='maplegend' 
        style='position: absolute; z-index: 9999; background-color: rgba(255, 255, 255, 0.5);
        border-radius: 6px; padding: 10px; font-size: 10.5px; width: 170px; height: 110px; right: 20px; top: 20px; cursor: move;'>     
    <div class='legend-scale'>
    <ul class='legend-labels'>
        <li><svg height="12" width="12">
            <polygon points="5,0 10,3.33 10,8.67 5,12 0,8.67 0,3.33" style="fill:none;opacity: 0.5;stroke:black" />
            </svg>Area with Samples</li>
        <li><svg height="12" width="12">
            <polygon points="5,0 10,3.33 10,8.67 5,12 0,8.67 0,3.33" style="fill:red;opacity: 0.7;stroke:none" />
            </svg>Isolated Population</li>
        <li><svg height="12" width="10"><line x1="0" y1="2" x2="10" y2="10" style="stroke:green;stroke-width:2" /></svg>Possible Migration Route</li>
    </ul>
    </div>
    <div class='legend-gradient'>
        <span style="font-weight: bold;">Genetic Distances</span>
        <span style='background: linear-gradient(to right, 
            rgb(242, 222, 140),  /* Sand color yellow */
            rgb(235, 199, 123),  /* Intermediate yellow-brown */
            rgb(222, 176, 97),   /* Yellow-brown tone */
            rgb(199, 159, 86),   /* Intermediate brown */
            rgb(169, 184, 186),  /* Light blue-gray */
            rgb(145, 194, 218),  /* Intermediate light blue */
            rgb(123, 209, 242),  /* Light blue */
            rgb(107, 194, 232),  /* Intermediate turquoise */
            rgb(92, 178, 242)    /* Turquoise blue */
        );
 width: 100%; height: 10px; display: block;'></span>
        <div style='display: flex; justify-content: space-between;'>
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


    macro = MacroElement()
    macro._template = Template(template)

    macro.add_to(m)
    return m
    