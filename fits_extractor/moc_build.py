# fits_extractor library
# Author : Nicolas Obrier

# moc_build module

from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from astropy.wcs import WCS


def image_corners_in_deg(header):
    """
    Calculates the sky coordinates of the four corners of the image using the WCS information from the header.

    Parameters:
    header (Header): The FITS header containing WCS information.

    Returns:
    SkyCoord: A SkyCoord object containing the RA and DEC of the four corners.
    """
    # Create a WCS object to handle the image projection
    wcs = WCS(header, naxis=2)
    
    # Get the image dimensions
    nx, ny = header['NAXIS1'], header['NAXIS2']
    
    # Define the pixel coordinates of the four corners
    pixel_corners = np.array([[0, 0], [0, ny-1], [nx-1, ny-1], [nx-1, 0]])
    
    # Convert the pixel coordinates to sky coordinates (RA, DEC)
    sky_corners = wcs.pixel_to_world(pixel_corners[:, 0], pixel_corners[:, 1])
    
    return sky_corners

def stcs_construct(corners, ref="ICRS"):
    """
    Constructs an STC (Space-Time Coordinate) string representation of a polygon given its corner coordinates.

    Parameters:
    corners (SkyCoord): A SkyCoord object containing the RA and DEC of the polygon corners.
    ref (str): The reference frame for the coordinates, default is "ICRS".

    Returns:
    str: The STCS string representation of the polygon.
    """
    # Create the STC representation by formatting the RA and DEC of each corner
    stcs_representation = f"Polygon {ref} " + " ".join(f"{ra:.6f} {dec:.6f}" for ra, dec in zip(corners.ra.deg, corners.dec.deg))
    
    return stcs_representation

def fov_center_moc(corners):
    """
    Calculates the center of the field of view (FOV) and the maximum FOV distance.

    Parameters:
    corners (SkyCoord): A SkyCoord object containing the RA and DEC of the four corners of the FOV.

    Returns:
    tuple: The field of view (FOV) in degrees and the SkyCoord object representing the center of the FOV.
    """
    # Calculate the center of the image by averaging the RA and DEC of the corners
    center_ra = corners.ra.mean()
    center_dec = corners.dec.mean()
    center = SkyCoord(center_ra, center_dec, unit="deg", frame="icrs")

    # Calculate the separation (distance) between the center and each corner
    distances = center.separation(corners)

    # Find the maximum distance to define the FOV
    max_distance = distances.max()

    # Set the FOV as four times the maximum distance
    fov_im = 4 * max_distance

    return fov_im, center

def calculate_moc_depth(fov_deg, k=10):
    """
    Calculates the optimal MOC depth based on the field of view (FOV) in degrees.

    Parameters:
    fov_deg (float): The field of view in degrees.
    k (int): A constant factor to adjust the resolution, default is 10.

    Returns:
    int: The optimal MOC depth, constrained to be between 0 and 29.
    """
    # Calculate the target resolution in degrees
    target_resolution = fov_deg / k
    
    # Calculate the MOC depth based on the resolution
    moc_depth = int(np.ceil(np.log2(180 / (target_resolution * np.pi))))
    
    # Return the MOC depth constrained to a valid range (0 to 29)
    return max(0, min(moc_depth, 29))