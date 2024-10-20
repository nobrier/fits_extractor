# To do : 
# - Create a function that homogenized the result (units, Object name are conform (Sesame ?))

# - Create a dataset (have the relevant metadata in it)
# - Append the MOC and the polygon that define image coverage 
#        - Use astropy to compute the Polygon (ie: the image region which is rectangular) and its MOC from the WCS header

# - Build a function which tests if a point is in a polygon (See Annexes in the subject)

# - Make a Python library (pip).

from astropy.io import fits

NGC_4486_path = "../data/NGC_4486_MIPS_M1.fits"
UGC_09618_2MASS_H_path = "../data/UGC_09618_2MASS_H.fits"
SN1987A_87_smHB_path = "../data/SN1987A_87_smHB.fits"



all_keywords = ["NAXIS", "NAXIS1", "NAXIS2", "OBJECT", "RA", "DEC", "RADESYS", "DATE_OBS", "MJD", "MJD_OBS", "EXPTIME", "INSTRUME", "TELESCOP", "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CDELT1", "CDELT2", "CD1_1", "CD1_2", "CD2_1", "CD2_2", "PC1_1", "PC1_2", "PC2_1", "PC2_2", "CROTA1", "CROTA2"]


def metadata_image(path_file: str, keyword_extracted: list) -> dict:
    """Get selected metadata from a fits image
    :param path_file: path to the fits file
    :param keyword_extracted: keywords/cards extracted from the whole metadata of the fits file
    :return: dictionnary of the metadata selected for the fits file
    """
    fits_file = fits.open(path_file)
    metadata ={}

    for i, hdu in enumerate(fits_file):
        header = fits_file[0].header 
        for keyword, value in header.items():
            for j in range (len(keyword_extracted)):
                if keyword == keyword_extracted[j]:
                    metadata[keyword] = header[keyword]                     
                   
    fits_file.close()
    return metadata

descriptions_dict = {
    "NAXIS": "Integer specifying the number of axes in data",
    "NAXIS1": "Number of elements along the first axis of data",
    "NAXIS2": "Number of elements along the second axis of data",
    "OBJECT": "Name of the observed object",
    "RA": "Right Ascension of the object",
    "DEC": "Declination of the object",
    "RADESYS": "Reference system for RA and DEC",
    "DATE_OBS": "Date of the observation",
    "MJD": "Modified Julian Date of the observation",
    "MJD_OBS": "[days] MJD for the observation",
    "EXPTIME": "[sec] Exposure time of the observation",
    "INSTRUME": "Name of the instrument used",
    "TELESCOP": "Name of the telescope used",
    "CTYPE1": "Type for the first coordinate axis",
    "CTYPE2": "Type for the second coordinate axis",
    "CUNIT1": "Units for CRVAL and CDELT for Axis 1",
    "CUNIT2": "Units for CRVAL and CDELT for Axis 2",
    "CRVAL1": "World coordinate value at the reference point",
    "CRVAL2": "World coordinate value at the reference point",
    "CRPIX1": "Location of reference point in pixels (Axis 1)",
    "CRPIX2": "Location of reference point in pixels (Axis 2)",
    "CDELT1": "Increment of world coordinate for Axis 1",
    "CDELT2": "Increment of world coordinate for Axis 2",
    "CD1_1": "Matrix element for linear transformation (1,1)",
    "CD1_2": "Matrix element for linear transformation (1,2)",
    "CD2_1": "Matrix element for linear transformation (2,1)",
    "CD2_2": "Matrix element for linear transformation (2,2)",
    "PC1_1": "Transformation matrix element (1,1)",
    "PC1_2": "Transformation matrix element (1,2)",
    "PC2_1": "Transformation matrix element (2,1)",
    "PC2_2": "Transformation matrix element (2,2)",
    "CROTA1": "Rotation from the standard coordinate system",
    "CROTA2": "Rotation from the standard coordinate system"
}

homo_keywords = ["NAXIS", "NAXIS1", "NAXIS2", "OBJECT", "RA", "DEC", "RADESYS", "DATE_OBS", "MJD_OBS", "EXPTIME", "INSTRUME", "TELESCOP", "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CDELT1", "CDELT2", "CD1_1", "CD1_2", "CD2_1", "CD2_2", "CROTA1", "CROTA2"]

def homogenized_metadata(metadata: dict, homogenized_keywords: list, descriptions: dict) -> dict:
    """Homogenize metadata of a fits file (PC -> CD, units,...)
    :param metadata: metadata to homogenize
    :param homogenized_keywords: keywords/cards that are homogeneous for every fits file
    :return: dictionnary of the metadata homogenized
    """
    homo_metadata = {}
    if "PC1_1" in metadata: #PCi_i to CDi_i transformation
        metadata["CD1_1"] = metadata["PC1_1"] * metadata["CDELT1"]
        metadata["CD1_2"] = metadata["PC1_2"] * metadata["CDELT1"]
        metadata["CD2_1"] = metadata["PC2_1"] * metadata["CDELT2"]
        metadata["CD2_2"] = metadata["PC2_2"] * metadata["CDELT2"]
    #############################
    #Add units homogenization
    #############################    
    for keyword in homo_keywords:
        if keyword in metadata:  
            value = metadata[keyword]
            description = descriptions[keyword]
            homo_metadata[keyword] = {
                'value': value,
                'description': description
            }
    return homo_metadata

# def object_check(metadata: dict, ...) -> dict:
    

            

print(fits.getheader(SN1987A_87_smHB_path))   

print("----------------------------------------")

print(homogenized_metadata(metadata_image(SN1987A_87_smHB_path, all_keywords), homo_keywords, descriptions_dict))
