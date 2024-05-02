from h3 import h3
import folium
from folium import Map, Element
from func import normalize_distances
import matplotlib.colors as mcolors
import base64


# function that draws hexagons on a map
def draw_hexagons(hexagons, m=None, color='darkgreen', zoom_start=1, value=None, opacity=0.5, imputed=False):
    # Create a map if it is not provided
    if m is None:
        m = folium.Map(location=(0.0, 0.0), tiles="Esri worldstreetmap", zoom_start=zoom_start)

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
            if value:
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
        polyline = folium.PolyLine(barrier, color = col)
        polyline.add_child(folium.Tooltip(value))
        polyline.add_to(m)
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
            polyline = folium.PolyLine(locations=line, color=color)
            polyline.add_child(folium.Tooltip(distance))
            polyline.add_to(m)
    
    return m

def get_color_gradient():
    # Define the colors for the colormap
    colors =[(1, 0.9, 0.7), (0.8, 0.5, 0.2), (0.4, 0.5, 0.9), (0, 0.1, 0.4)]
    # Create the colormap
    cmap = mcolors.LinearSegmentedColormap.from_list("costum_color_gradient", colors)
    
    return cmap
    