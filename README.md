# SoundingDataPlotter (v0.1)
Group of Python scripts to plot University of Wyoming (UWYO) and processed OWLETS-2 sounding data.

WARNING: A lot of these scripts were made for personal use and are not optimized. Use at your own risk...

These scripts are meant to be run in the same directory. The hierarchy as shown here is what you should emulate on your device for them to work.

## Folder/File Descriptions
The sounding data in this repository are from the [OWLETS-2](https://www-air.larc.nasa.gov/cgi-bin/ArcView/owlets.2018?SONDE=1) campaign conducted during the summer of 2018. These soundings were performed at Hart Miller Island (HMI), Howard University Beltsville (HUBV), and the University of Maryland, Baltimore County (UMBC).

**Full Sounding Data:**

`Raw Soundings` - All of the sounding data converted from ICARRT to csv, separated by each launch with columns holding a different measurement variable. The data between each sounding may vary due to a different type of model used, so be sure to check the "Header" column for more details.

**Processed Data:**
`Filenames.csv` - Contains the locations for all of the processed sounding files.

`Derived MR` - These are the soundings where I derived the mixing ratio using the Arden Buck Eqaution

`Measured MR` - Soundings where the mixing ratio was already included

## File Descriptions

In my processed data folders, there are 7 columns:

O3 - Ozone; Displayed in parts per billion by volume (ppbv)

Theta_K - Potential Temperature; Displayed in Kelvin (K)

RH - Relative Humidity; Displayed in %

q - Specific Humidity; Displayed in g/kg

N - Atmospheric Refractivity; (unitless)

Alt - Altitude; Displayed in meters (m)

malr - Moist Adiabatic Lapse Rate; Displayed in K/km

Ri - Richardson Index (*Note: Treat as preliminary data, may be inaccurate*); (unitless)

## Scripts
The scripts described below require the following 3rd party packages: matplotlib, metpy, numpy, pandas, siphon, and scipy. They are compatible with >= Python 3.8.5.

`OzoneDataProcessor.py`: Plots the vertical profiles of ozone, potential temperature, relative humidity, specific humidity, atmospheric refractivity, and stability using the processed files described above. In addition, this script will also detect the planetary boundary layer using a methodology similar to the one described in [Wang et al. (2014)](https://amt.copernicus.org/articles/7/1701/2014/amt-7-1701-2014.pdf). 
An example of one of the created plots is shown below.

![HUBV29-Jun-2018 172732](https://user-images.githubusercontent.com/94017926/141014186-5fe346dd-ec26-4ad8-a1fe-d188bd487da8.jpg)

`PlotSoundings.py`: Similar to `OzoneDataProcessor.py`, but contains functions rather than being a single script. `Example.py` shows how to use some of these functions to process a UWYO Sounding.

Feel free to submit a pull request for new features and/or bug fixes.
