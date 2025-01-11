# fits_extractor library
# Author : Nicolas Obrier

# table_build module

from astropy.io.fits import Header
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import re
import requests
from astroquery.simbad import Simbad


def extract_header(fits_file: str):
    """
    Extracts the header of the primary HDU from a FITS file.

    Parameters:
    fits_file (str): Path to the FITS file.

    Returns:
    header (Header): Header object of the primary HDU.
    """
    with fits.open(fits_file) as hdul:
        header = hdul[0].header
    return header

def basic_metadata_image(header: Header, keywords_selected: list) -> dict:
    """Get selected metadata from a fits image
    :param path_file: path to the fits file
    :param keywords_selected: list of keywords to extract
    :return: dictionary of the metadata selected for the fits file
    """
    metadata = {}

    for keyword in keywords_selected:
        if keyword in header:
            metadata[keyword] = header[keyword]
        else :
            metadata[keyword] = None
    return metadata

def searching_metadata_image(header: Header, keywords_selected: list, metadata: dict) -> dict:
    """
    Searches the FITS header for keywords or patterns based on selected keywords
    and updates the metadata dictionary with found values.

    Parameters:
    - header: Header object from astropy.io.fits
    - keywords_selected: List of keywords to search for
    - metadata: Dictionary to store found metadata

    Returns:
    - Updated metadata dictionary
    """
    for keyword in keywords_selected:
        if metadata.get(keyword) is None:  # Check if the keyword is not already in metadata

            # Split the keyword by '-' to separate the parts
            parts = keyword.split('-')
            if len(parts) == 2:
                # Create the regex pattern to match the two parts with any characters in between
                pattern = re.compile(f"{re.escape(parts[0])}.*{re.escape(parts[1])}")

                # Search in the header
                for key in header.keys():
                    if pattern.match(key):  # If a match is found
                        metadata[keyword] = header[key]  # Update the metadata with the value from the header
                        break  # Stop searching after the first match

    return metadata

def transform_pc_to_cd(header: Header, metadata: dict) -> dict:
    """
    Converts PC matrix elements to CD matrix elements in the metadata.
    
    If neither PC nor CD are present, sets CD to an identity matrix.

    Parameters:
    header (Header): FITS header object containing the transformation matrix.
    metadata (dict): Dictionary containing metadata such as NAXIS and CDELT.

    Returns:
    dict: Updated metadata with CD matrix elements.
    """
    # Number of axes in the FITS data
    naxis = metadata["NAXIS"]

    # Check if CD1_1 is missing in metadata, which may indicate absence of CD matrix
    if metadata.get("CD1_1") is None:
        # Check if any PC matrix keywords exist in the header
        has_pc_matrix = any(key.startswith("PC") for key in header.keys())

        # If PC matrix exists, convert PC to CD
        if has_pc_matrix:
            for key in header.keys():
                if key.startswith("PC"):
                    # Extract the indices i and j from the PC key (e.g., PC1_1 -> i=1, j=1)
                    i, j = map(int, key[2:].split('_'))

                    # Construct the corresponding CD key (e.g., PC1_1 -> CD1_1)
                    cd_key = f"CD{i}_{j}"

                    # Check if CDELT{i} is present in metadata and ensure i, j are within valid range
                    if f"CDELT{i}" in metadata and 1 <= i <= naxis and 1 <= j <= naxis:
                        # Convert PC matrix element to CD matrix element
                        metadata[cd_key] = metadata[f"CDELT{i}"] * header[key]

        else:
            # If neither PC nor CD matrices are found, set CD to identity matrix
            for i in range(1, naxis + 1):
                for j in range(1, naxis + 1):
                    cd_key = f"CD{i}_{j}"

                    # Set the diagonal elements to 1 (identity matrix), off-diagonal to 0
                    metadata[cd_key] = 1.0 if i == j else 0.0

    return metadata

def wcs_to_coord(header: Header) -> tuple:
    """
    Extracts the RA (Right Ascension) and DEC (Declination) coordinates of an object
    from a FITS header using WCS (World Coordinate System) information. Supports 
    multi-dimensional WCS data and uses the first two axes for 2D projections.

    Parameters:
    header (Header): The FITS header containing WCS information.

    Returns:
    tuple: A tuple containing the RA and DEC coordinates.
    """
    # Number of axes in the FITS image (NAXIS)
    naxis = header['NAXIS']
    
    # Check if the image has more than 2 axes (multi-dimensional WCS)
    if naxis > 2:
        # If the image has more than 2 axes, we will use the first two axes for RA/DEC
        # print(f"Warning: Image has {naxis} dimensions. Using the first two axes for RA/DEC.")
        # Initialize WCS with just the first two axes
        wcs = WCS(header, naxis=2)
    else:
        # For 2D WCS, initialize WCS normally
        wcs = WCS(header)
    
    # Get the default pixel positions for CRPIX1 and CRPIX2, usually the center of the image
    # If these values are not available in the header, default to 0
    x_pixel, y_pixel = header.get('CRPIX1', 0), header.get('CRPIX2', 0)
    
    # Convert pixel coordinates to world coordinates (RA, DEC)
    # The '1' at the end of the function call indicates 1-based indexing for the pixel positions
    x_coord, y_coord = wcs.wcs_pix2world(x_pixel, y_pixel, 1)
    
    return x_coord, y_coord

def homogenized_metadata(header: Header, metadata: dict, keywords_selected: list, descriptions: dict) -> dict:
    """
    Homogenizes the metadata of a FITS file by performing transformations and adding descriptions.
    This includes transforming PC matrix elements to CD matrix elements, setting default values 
    for missing CDELT1 and CDELT2, and adding descriptions for each metadata field.

    Parameters:
    header (Header): The FITS header containing the information to homogenize.
    metadata (dict): The metadata dictionary to be homogenized.
    keywords_selected (list): List of selected keywords to extract from the FITS header.
    descriptions (dict): A dictionary mapping keywords to their descriptions.

    Returns:
    dict: A dictionary containing the homogenized metadata with values and descriptions.
    """
    # Step 1: Extract basic metadata from the FITS header based on selected keywords
    metadata = basic_metadata_image(header, keywords_selected)

    # Step 2: Search for additional metadata in the FITS header and update the metadata dictionary
    metadata = searching_metadata_image(header, keywords_selected, metadata)

    # Step 3: Transform PC matrix elements to CD matrix elements if they are present in the header
    metadata = transform_pc_to_cd(header, metadata)

    # Step 4: Set default values for CDELT1 and CDELT2 if they are not present in the metadata
    # These represent the pixel scale along each axis (e.g., degrees/pixel for RA/DEC)
    if metadata.get("CDELT1") is None:
        metadata["CDELT1"] = 1  # Default value for the first axis scale
    if metadata.get("CDELT2") is None:
        metadata["CDELT2"] = 1  # Default value for the second axis scale

    # Step 5: Add descriptions for each metadata field using the provided descriptions dictionary
    homo_metadata = {}
    for key, value in metadata.items():
        # If a description is available for the keyword, use it; otherwise, provide a default message
        description = descriptions.get(key, "No description available")
        
        # Store both the value and description in the output dictionary
        homo_metadata[key] = {
            'value': value,
            'description': description
        }

    # Return the homogenized metadata
    return homo_metadata

def object_check(metadata):
    """
    Checks if the "OBJECT" value in metadata is recognized by the SESAME service (via NED, Simbad, Vizier).
    If found, the "OBJECT" value is standardized with the primary identifier.

    Parameters:
    metadata (dict): Dictionary containing an "OBJECT" key with a "value" sub-key, representing the object name.

    Returns:
    dict: Updated metadata dictionary with the "OBJECT" value standardized if a match is found in SESAME.
    """
    #Flag to know if the function returned a new corrected name
    find = False
    
    # URL to query SESAME service for object identifiers
    url_sesame = "https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/NSV"  # Query identifiers from NED, Simbad, Vizier (NSV)

    # Make a copy of the metadata to avoid modifying the original
    metadata_checked = metadata.copy()

    # Ensure that the "OBJECT" key exists and has a non-empty value
    if "OBJECT" in metadata_checked and metadata_checked["OBJECT"]["value"]:
        # Retrieve the object name and standardize it (remove spaces, convert to lowercase)
        name = metadata_checked["OBJECT"]["value"]
        name_lower = name.lower().replace(" ", "")
        
        try:
            # Request object information from SESAME
            response = requests.get(f"{url_sesame}?{name_lower}")
            response.raise_for_status()  # Raise an exception if the HTTP request failed

            identifiers = []  # List to store possible identifiers from the response
            primary_id = None  # Variable to store the primary identifier

            # Parse the SESAME response
            for line in response.text.splitlines():
                line = line.strip()

                # If the line starts with "%I.0", it's the primary identifier
                if line.startswith("%I.0"):
                    primary_id = re.sub(r"\s+", "", line.split(" ", 1)[1])  # Clean up the identifier
                # If the line starts with "%I", it's another identifier
                elif line.startswith("%I"):
                    identifier = re.sub(r"\s+", "", line.split(" ", 1)[1])  # Clean up the identifier
                    identifiers.append(identifier)

            # Check if the object name matches any of the identifiers
            if any(name_lower == identifier.lower() for identifier in identifiers):
                # If a match is found, update the "OBJECT" value with the primary identifier
                metadata_checked["OBJECT"]["value"] = primary_id
                find == True
            elif primary_id and name_lower == primary_id.lower():
                # If the name matches the primary identifier, print a message
                print(f"{name_lower} is already the main identifier.")
            else:
                # If no match is found, print a message
                print(f"{name} not found in SESAME identifiers.")

        except requests.RequestException as e:
            # Handle errors in the HTTP request (e.g., connection issues)
            print(f"Error accessing SESAME service: {e}")
        except Exception as e:
            # Catch any other unexpected errors
            print(f"An unexpected error occurred: {e}")
    else:
        # If the "OBJECT" value is not found or is empty, print a message
        print('No "OBJECT" value in the metadata provided.')

    # Return the updated metadata (whether modified or not)
    return (metadata_checked, find)

def object_name_coords(metadata):
    """
    Searches for the name of an astronomical object in SIMBAD using RA and DEC coordinates
    provided in the metadata dictionary. If an object is found, it updates the metadata with
    the object name.

    Parameters:
    metadata (dict): Dictionary containing "RA" and "DEC" keys with their respective values 
                     (Right Ascension and Declination in degrees).

    Returns:
    dict: Updated metadata dictionary with the "OBJECT" value, if the object is found in SIMBAD.
    """
    # Initialize the Simbad query interface
    simbad = Simbad()

    try:
        # Step 1: Create a SkyCoord object using RA and DEC values from metadata
        # The units are specified in degrees (u.deg)
        coord = SkyCoord(metadata["RA"]["value"], metadata["DEC"]["value"], unit=(u.deg, u.deg))

        # Step 2: Query SIMBAD for the object at the given coordinates
        result_table = simbad.query_region(coord)

        # Step 3: Check if any object was found in the query result
        if result_table is not None and len(result_table) > 0:
            # Step 4: Update the "OBJECT" value in metadata with the name of the first object found
            # print(result_table[0])
            metadata["OBJECT"]["value"] = result_table[0]["main_id"]  # Extract object name
        else:
            # If no object was found, print a message
            print("No object found for the given coordinates.")

    except Exception as e:
        # Catch any errors (e.g., connection issues or incorrect data) and print an error message
        print(f"An error occurred: {e}")

    # Return the updated metadata, which may or may not include the object name
    return metadata