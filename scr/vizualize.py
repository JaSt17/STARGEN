from h3 import h3
import folium
from folium import Map, Element
from func import normalize_distances
import matplotlib.colors as mcolors
import base64


# function that draws hexagons on a map
def draw_hexagons(hexagons, m=None, color='orange', zoom_start=1, value=None):
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
                fill_opacity=0.4,
                fill=True
            )
            if value:
                polygon.add_child(folium.Tooltip(value))
            
            polygon.add_to(m)
    return m

def draw_hexagons_with_values(hex_dict, m=None, zoom_start=1, threshold=0.0):
    hexagons = hex_dict.keys()
    values = hex_dict.values()
    
    # create a color gradient to color the lines based on the normalized distance
    colors = [(1, 0.2, 0), (0, 0, 0.5)]  # Dark blue to orange
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_darkblue_to_orange", colors)

    # write the values to the center of each hexagon
    for hexagon, value in zip(hexagons, values):
        col = mcolors.to_hex(cmap(value))
        m = draw_hexagons([hexagon], m, color=col, zoom_start=zoom_start, value=value)
    return m

def draw_barriers(barriers_dict, m=None, zoom_start=1, threshold=0.0):
    barriers = barriers_dict.keys()
    values = barriers_dict.values()
    
    colors = [(1, 0.2, 0, 0.2), (0, 0, 0.5, 0.2)]  # Dark blue to orange with reduced opacity
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_darkblue_to_orange", colors)
    
    for barrier, value in zip(barriers, values):
        col = mcolors.to_hex(cmap(value))
        barrier = list(barrier)
        polyline = folium.PolyLine(barrier, color = col)
        polyline.add_child(folium.Tooltip(value))
        polyline.add_to(m)
    return m

# this function takes two hexagons and a map and draws a line between the two midpoints of the hexagons
def draw_borders(hexagon1, hexagon2, m, color, distance=None):
    # get the midpoint of both hexagons
    mid1 = h3.h3_to_geo(hexagon1)
    mid2 = h3.h3_to_geo(hexagon2)
    
    # check if the points are on the opposite sides of the antimeridian
    if abs(mid1[1] - mid2[1]) > 180:
        # if they are, add 360 to the longitude of the point with the smaller longitude
        if mid1[1] < mid2[1]:
            mid1 = (mid1[0], mid1[1] + 360)
        else:
            mid2 = (mid2[0], mid2[1] + 360)
            
    # Draw line between the two midpoints
    line = [mid1, mid2]
    folium.PolyLine(locations=line,
                    color=color,
                    weight=5).add_to(m)
    
        # Calculate the midpoint of the line
    midpoint = [sum(x)/len(x) for x in zip(*line)]

    # Create a custom icon for displaying text
    icon = folium.DivIcon(html=f'<div style="font-size: 12pt">{distance:.2f}</div>')

    # Add a marker at the midpoint with the custom text icon
    folium.Marker(location=midpoint, icon=icon).add_to(m)
    
    return m

# this function takes a time bin and a map and draws all neighboring lines for the hexagons in that time bin
def draw_all_boarders_for_time_bin(time_bin, m, color="red", threshold=0.0):
    # get the normalized distance for the color gradient
    normalized_time_bin = normalize_distances(time_bin)
    
    # create a color gradient to color the lines based on the normalized distance
    colors = [(1,1,0), (255,0,0)]
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_darkred_to_yellow", colors)
    
    # loop through all hexagons in the time bin
    for pair in time_bin:
        # only draw the line if the distance is below the threshold
        if time_bin[pair] >= threshold:
            # get the color for the line based on the normalized distance
            col = mcolors.to_hex(cmap(normalized_time_bin[pair]))
            hex1, hex2 = list(pair)
            m = draw_borders(hex1, hex2, m, color= col, distance=time_bin[pair])
    return m