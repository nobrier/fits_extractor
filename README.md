# fits_extractor Library

`fits_extractor` is a Python library designed for astronomers and astrophysicists to easily manipulate and analyze FITS (Flexible Image Transport System) files. The library provides tools for extracting metadata, handling MOCs (Multi-Order Coverage maps), and performing operations related to celestial geometry and field-of-view visualization.

## Features

- **Metadata Extraction and Homogenization:**
  - Extract relevant keywords from FITS headers.
  - Convert and homogenize WCS-related transformations (e.g., PC matrix to CD matrix).
  - Add detailed descriptions to metadata for clarity.

- **MOC and FOV Management:**
  - Calculate the corners of an image in sky coordinates.
  - Compute the optimal MOC depth based on field-of-view dimensions.
  - Generate STC representations of celestial regions.

- **Object Identification:**
  - Use the SIMBAD service to identify celestial objects from coordinates.
  - Verify and standardize object names using the SESAME service.

- **WCS Utilities:**
  - Transform pixel coordinates to celestial coordinates.
  - Extract RA/DEC from FITS headers with multi-dimensional WCS handling.

## Project Structure

The project follows a modular design with the following structure:

```
fits_extractor
├── README.md
├── pyproject.toml
├── requirements.txt
├── setup.py
├── data      # FITS files for testing (read-only)
│   ├── 2013.1.00034.S_SB_X5_GB_X6_MB_X7_midz_cell10_25342_sci.spw0_1_2_3.cont.I.image.fits
│   ├── 2013.1.01292.S_SB_X4eb_GB_X4ec_MB_X4ed_2-38011_sci.spw0_1_2_3.cont.I.image.fits
│   ├── 5GHz_n_f.fits
│   ├── G327.617-0.364_I4.fits
│   ├── G351.632-0.459_atlasgal.fits
│   ├── G351.702+0.672_atlasgal.fits
│   ├── G9_POLIN.fit
│   ├── N2.20100426.52760.fits
│   ├── NGC_4486_MIPS_M1.fits
│   ├── SN1987A_87_smHB.fits
│   ├── SN1987A_cut_35_smHB.fits
│   ├── SNaverage5-0.fits
│   ├── UGC_09618_2MASS_H.fits
│   ├── UGC_09618_S_2MASS_H.fits
│   ├── av_galcen_2mass.fits
│   ├── av_galcen_spitzer.fits
│   ├── id12_GAL-Survey-GC.fpsf.fits
│   ├── imageih.fit
│   ├── imagerf.fit
│   ├── u.fit
│   └── whsky072.fit
├── fits_extractor
│   ├── __init__.py
│   ├── ask_table.py                            # Modules for Asking the table created
│   ├── moc_build.py                            # Modules for MOC and FOV calculations
│   └── table_build.py                          # Metadata extraction and transformation
```

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/nobrier/fits_extractor.git
   cd fits_extractor
   ```

2. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Example 1: Extracting Metadata

```python
from fits_extractor.table_build import extract_header, homogenized_metadata

header = extract_header("example.fits")
keywords = ["NAXIS", "CDELT1", "CDELT2", "OBJECT"]
descriptions = {
    "NAXIS": "Number of axes in the FITS file",
    "CDELT1": "Pixel scale along the RA axis (degrees/pixel)",
    "CDELT2": "Pixel scale along the DEC axis (degrees/pixel)",
    "OBJECT": "Name of the observed object"
}
metadata = homogenized_metadata(header, {}, keywords, descriptions)
print(metadata)
```

### Example 2: Generating a MOC

```python
from fits_extractor.moc_build import image_corners_in_deg, fov_center_moc, calculate_moc_depth
from astropy.io import fits

# Load FITS header
with fits.open("example.fits") as hdul:
    header = hdul[0].header

# Calculate image corners
corners = image_corners_in_deg(header)

# Determine FOV center and size
fov_size, center = fov_center_moc(corners)

# Compute MOC depth
moc_depth = calculate_moc_depth(fov_size.deg)
print(f"MOC depth: {moc_depth}")
```

### Example 3: Verifying Object Names

```python
from fits_extractor.table_build import object_check

metadata = {"OBJECT": {"value": "M31"}}
updated_metadata, found = object_check(metadata)
if found:
    print("Object name corrected:", updated_metadata["OBJECT"]["value"])
else:
    print("Object name could not be verified.")
```

## Dependencies

- `astropy`
- `numpy`
- `astroquery`
- `requests`
- `matplotlib`
- `mocpy`
- `ipyaladin`
- `shapely`

## Author

Nicolas Obrier

