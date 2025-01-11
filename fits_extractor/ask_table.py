# fits_extractor library
# Author : Nicolas Obrier

# ask_table module

from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
import numpy as np
from shapely.geometry import Polygon, Point
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from mocpy import MOC, WCS
import matplotlib.pyplot as plt
from astropy.wcs import WCS as wcsastropy


def is_point_in_polygon(point, polygon):
    """
    Checks if a celestial point is inside a polygon defined by celestial corners.

    Parameters:
    point (SkyCoord): Celestial coordinates (RA, DEC) of the point to check.
    polygon (SkyCoord): Celestial coordinates (RA, DEC) of the vertices of the polygon.

    Returns:
    bool: True if the point is inside the polygon, otherwise False.
    """
    # Convert the polygon vertices into a list of (RA, DEC) tuples
    polygon_coords = [(ra, dec) for ra, dec in zip(polygon.ra.deg, polygon.dec.deg)]

    # Convert the point into a tuple (RA, DEC)
    point_coords = (point.ra.deg, point.dec.deg)

    # Utility function to subtract vectors
    def vector_subtract(p1, p2):
        return np.array([p1[0] - p2[0], p1[1] - p2[1]])

    # Initialize the sign of the first cross product
    n = len(polygon_coords)
    P1 = np.array(polygon_coords[0])
    P2 = np.array(polygon_coords[1])
    edge_vector = vector_subtract(P2, P1)
    point_vector = vector_subtract(point_coords, P1)
    initial_scalar = np.cross(edge_vector, point_vector)

    # Iterate through all edges of the polygon
    for i in range(n):
        # Current and next vertices (with modulo to loop back to the last vertex)
        P1 = np.array(polygon_coords[i])
        P2 = np.array(polygon_coords[(i + 1) % n])

        # Vectors for the edge and the point
        edge_vector = vector_subtract(P2, P1)
        point_vector = vector_subtract(point_coords, P1)

        # Calculate the cross product
        scalar = np.cross(edge_vector, point_vector)

        # If the sign differs from the initial one, the point is outside the polygon
        if initial_scalar * scalar < 0:
            return False

    # If all cross products have the same sign, the point is inside the polygon
    return True

def parse_polygon_stcs(polygon_str):
    """
    Parses an STCS string to extract the corner coordinates in RA and DEC.

    Example input: 
    "Polygon ICRS 83.980060 -69.309782 83.979642 -69.229658 83.753698 -69.229658 83.753280 -69.309782"

    Parameters:
    polygon_str (str): The STCS string representing the polygon with RA and DEC coordinates.

    Returns:
    SkyCoord: A SkyCoord object containing the RA and DEC coordinates of the polygon corners.
    """
    # Split the STCS string into parts and extract the reference frame
    parts = polygon_str.split()
    ref_frame = parts[1]  # For example: "ICRS"
    
    # Convert the RA and DEC coordinates from the string into a NumPy array, reshaped to pairs of (RA, DEC)
    coords = np.array(parts[2:], dtype=float).reshape(-1, 2)
    
    # Create SkyCoord objects for RA and DEC with units of degrees
    ra = coords[:, 0] * u.deg
    dec = coords[:, 1] * u.deg
    
    # Return the SkyCoord object for the polygon corners, with the appropriate frame
    return SkyCoord(ra, dec, frame=ref_frame.lower())

def fov_center_from_corners(corners):
    """
    Calculates the center and field of view (FoV) from the corner coordinates.

    Parameters:
    corners (SkyCoord): A SkyCoord object containing the RA and DEC coordinates of the corners.

    Returns:
    tuple: A tuple containing the FoV in degrees and the SkyCoord object representing the center of the corners.
    """
    # Calculate the mean RA and DEC to find the center of the corners
    center_ra = corners.ra.mean()
    center_dec = corners.dec.mean()
    center = SkyCoord(center_ra, center_dec, unit="deg", frame=corners.frame)

    # Calculate the distance from the center to each corner
    distances = center.separation(corners)
    
    # Find the maximum distance to define the FoV
    max_distance = distances.max()

    # Define the FoV as four times the maximum distance from the center
    fov_im = 4 * max_distance
    
    return fov_im, center


def is_polygon_intersecting_circle(corners, circle_center, radius):
    """
    Check if a polygon intersects a circle.

    Parameters:
        corners (list of tuples): List of (RA, DEC) coordinates defining the polygon vertices.
        circle_center (SkyCoord): Center of the circle.
        radius (Angle): Radius of the circle.

    Returns:
        bool: True if the polygon intersects the circle, False otherwise.
    """
    # Create a Shapely Polygon
    polygon = Polygon(corners)

    # Create a Shapely Point and Buffer (circle)
    circle = Point(circle_center.ra.deg, circle_center.dec.deg).buffer(radius.deg)

    # Check for intersection
    return polygon.intersects(circle)

def load_moc(row):
    """Load the MOC from a row."""
    moc_string = row.get("MOC")
    if not moc_string:
        raise ValueError("MOC string is missing.")
    return MOC.from_string(moc_string)


def parse_polygon(row):
    """Extract the corners and compute the FOV and center."""
    corners = parse_polygon_stcs(row["Polygon"])
    fov_im, center = fov_center_from_corners(corners)
    return corners, fov_im, center


def extract_wcs_parameters(row):
    """Extract and validate WCS parameters from the row."""
    required_params = ["NAXIS1", "NAXIS2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2"]
    wcs_params = {param: row.get(param, None) for param in required_params}
    if None in wcs_params.values():
        raise ValueError("Missing or invalid WCS parameters.")

    wcs_params["CDELT1"] = row.get("CDELT1", 1)
    wcs_params["CDELT2"] = row.get("CDELT2", 1)
    wcs_params["CTYPE"] = row.get("CTYPE1", "SIN")
    if "---" in wcs_params["CTYPE"]:
        wcs_params["CTYPE"] = wcs_params["CTYPE"].split("---")[1]

    cd1_2 = row.get("CD1_2", 0)
    cd1_1 = row.get("CD1_1", 1)
    angle = np.arctan(float(cd1_2) / float(cd1_1)) if cd1_2 else 0
    wcs_params["ANGLE"] = Angle(angle * 180 / np.pi, u.degree)
    wcs_params["FRAME"] = "icrs"  # Default frame

    return wcs_params


def create_wcs_object(fig, fov, center_coord, wcs_params):
    """Create a WCS object using mocpy."""
    return WCS(
        fig,
        fov=fov,
        center=center_coord,
        coordsys=wcs_params["FRAME"],
        rotation=wcs_params["ANGLE"],
        projection=wcs_params["CTYPE"],
    )


def plot_moc_and_points(ax, moc, wcs, coord_in, coord_out, title):
    """Plot the MOC and add the points to the WCS plot."""
    moc.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="red", linewidth=1)
    moc.border(ax=ax, wcs=wcs, alpha=1, color="red")

    # Add points
    ax.scatter(
        coord_in.ra.deg, coord_in.dec.deg, transform=ax.get_transform("icrs"),
        s=100, c="green", marker="x", label="Point IN"
    )
    ax.scatter(
        coord_out.ra.deg, coord_out.dec.deg, transform=ax.get_transform("icrs"),
        s=100, c="red", marker="x", label="Point OUT"
    )

    # Add labels, grid, and legend
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")
    ax.set_title(title)
    ax.grid(color="black", linestyle="dotted")
    ax.legend(loc="upper right")


def plot_row_moc_points(row, coord_in, coord_out):
    """Main function to plot the MOC and points for a single row."""
    try:
        moc = load_moc(row)
        corners, fov_im, center = parse_polygon(row)
        wcs_params = extract_wcs_parameters(row)

        center_coord = SkyCoord(
            wcs_params["CRVAL1"], wcs_params["CRVAL2"], unit="deg", frame=wcs_params["FRAME"]
        )

        fig = plt.figure(figsize=(15, 10))
        with create_wcs_object(fig, fov_im, center_coord, wcs_params) as wcs_new:
            ax = fig.add_subplot(111, projection=wcs_new)
            plot_moc_and_points(
                ax, moc, wcs_new, coord_in, coord_out,
                f"Coverage of {row['OBJECT']} with Points IN/OUT"
            )
        plt.show()

    except Exception as e:
        print(f"An error occurred: {e}")

def plot_coverage(row):
    """Plot the MOC coverage of a row."""
    try:
        # Extract data
        moc = load_moc(row)
        corners, fov_im, center = parse_polygon(row)
        wcs_params = extract_wcs_parameters(row)

        center_coord = SkyCoord(
            wcs_params["CRVAL1"], wcs_params["CRVAL2"], unit="deg", frame=wcs_params["FRAME"]
        )

        # Create plot
        fig = plt.figure(figsize=(15, 10))
        with create_wcs_object(fig, fov_im, center_coord, wcs_params) as wcs_new:
            ax = fig.add_subplot(111, projection=wcs_new)
            moc.fill(ax=ax, wcs=wcs_new, alpha=0.5, fill=True, color="red", linewidth=1)
            moc.border(ax=ax, wcs=wcs_new, alpha=1, color="red")
            # Add labels and title
            ax.set_xlabel("RA")
            ax.set_ylabel("DEC")
            ax.set_title(f"Coverage of {row['OBJECT']}")
            ax.grid(color="black", linestyle="dotted")

        plt.show()

    except Exception as e:
        print(f"An error occurred while processing {row.get('OBJECT', 'Unknown')}: {e}")


def rows_containing_point(fits_table, coord):
    """Find rows whose polygons contain the given point."""
    rows = []
    for index, row in enumerate(fits_table):
        try:
            corners = parse_polygon_stcs(row["Polygon"])
            if is_point_in_polygon(coord, corners):
                rows.append(row)
        except Exception as e:
            print(f"An error occurred while processing row {index}: {e}")
    return rows


def rows_intersecting_circle(fits_table, circle_center, radius):
    """Find rows whose polygons intersect a circle defined by center and radius."""
    rows = []
    for index, row in enumerate(fits_table):
        try:
            corners = parse_polygon_stcs(row["Polygon"])
            # Convert corners to a list of (RA, DEC) tuples
            corners_list = [(corner.ra.deg, corner.dec.deg) for corner in corners]
            if is_polygon_intersecting_circle(corners_list, circle_center, radius):
                rows.append(row)
        except Exception as e:
            print(f"An error occurred while processing row {index}: {e}")
    return rows

# Function to plot the union of MOCs with a search circle
def plot_union_of_mocs_with_circle(fits_table, circle_center, circle_radius):
    """
    Plot the union of all MOCs in the FITS table and overlay a search circle.

    Parameters:
        fits_table (Table): FITS table containing MOCs.
        circle_center (SkyCoord): Center of the search circle.
        circle_radius (Quantity): Radius of the search circle in angular units.
    """
    # Initialize moc_union with the MOC from the first row
    try:
        moc_union = load_moc(fits_table[0])
    except Exception as e:
        raise ValueError(f"Error extracting MOC from the first row: {e}")

    # Compute the union of the remaining MOCs
    for index in range(1, len(fits_table)):
        try:
            moc = load_moc(fits_table[index])
            moc_union = moc_union + moc
        except Exception as e:
            print(f"Error processing row {index}: {e}")

    # Define a WCS for plotting
    wcs = wcsastropy(
        {
            "naxis": 2,
            "naxis1": 3240,
            "naxis2": 1620,
            "crpix1": 1620.5,
            "crpix2": 810.5,
            "cdelt1": -0.1,
            "cdelt2": 0.1,
            "ctype1": "RA---AIT",
            "ctype2": "DEC--AIT",
        },
    )

    # Create the figure and axes with the Elliptical frame
    fig = plt.figure(figsize=(21, 14))
    ax = fig.add_subplot(1, 1, 1, projection=wcs, frame_class=EllipticalFrame)

    # Plot the union of all MOCs
    moc_union.fill(
        ax=ax,
        wcs=wcs,
        alpha=0.5,
        fill=True,
        color="blue",
        linewidth=0,
        label="Union of MOCs",
    )
    moc_union.border(ax=ax, wcs=wcs, alpha=1, color="black")

    # Generate circle coordinates in RA/Dec
    theta = np.linspace(0, 2 * np.pi, 100)
    ra_circle = circle_center.ra.deg + circle_radius.to_value(u.deg) * np.cos(theta) / np.cos(circle_center.dec.radian)
    dec_circle = circle_center.dec.deg + circle_radius.to_value(u.deg) * np.sin(theta)

    # Convert RA/Dec to the world coordinate system for plotting
    circle_coords = SkyCoord(ra=ra_circle * u.deg, dec=dec_circle * u.deg, frame="icrs")
    ra_proj = circle_coords.ra.wrap_at(180 * u.deg).deg
    dec_proj = circle_coords.dec.deg

    # Plot the circle in the Mollweide projection
    ax.plot(
        ra_proj,
        dec_proj,
        transform=ax.get_transform("world"),
        color="red",
        linewidth=1.5,
        label="Search Circle"
    )

    # Add labels, legend, and grid
    ax.legend()
    ax.set_aspect(1.0)
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.title("All MOCs with Search Circle")
    plt.grid(color="black", linestyle="dotted")

    # Show the plot
    plt.show()