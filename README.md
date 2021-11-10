# SoundingDataPlotter (v0.1)
Group of Python scripts to plot UWYO and processed OWLETS-2 sounding data.

WARNING: These scripts were made by an amatuer wannabe data science poser who only slightly knows what he is doing. Use at your own risk...

These scripts are meant to be run in the same directory. The hierarchy as shown here is what you should emulate on your device for them to work.

## Folder Descriptions
`Helpers` - Contains the helper files used to create the visualizations and process the data.

The sounding files in this repository are from the [OWLETS-2](https://www-air.larc.nasa.gov/cgi-bin/ArcView/owlets.2018?SONDE=1) campaign conducted during the summer of 2018. These soundings were performed at Hart Miller Island (HMI), Howard University Beltsville (HUBV), and the University of Maryland, Baltimore County (UMBC).

**Full Sounding Data:**

`Raw Soundings` - All of the sounding data converted from ICARRT to csv, separated by each launch with columns holding a different measurement variable. The data between each sounding may vary due to a different type of model used, so be sure to check the "Header" column for more details.

**Processed Data:**

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
The scripts described below require the following 3rd party packages: matplotlib, metpy, numpy, pandas, and scipy. They are compatible with >= Python 3.8.5.

`OzoneDataProcessor.py`: Plots the vertical profiles of ozone, potential temperature, relative humidity, specific humidity, atmospheric refractivity, and stability using the processed files described above. An example of one of the created plots is shown below.

![HUBV29-Jun-2018 172732](https://user-images.githubusercontent.com/94017926/141014186-5fe346dd-ec26-4ad8-a1fe-d188bd487da8.jpg)


`SoundingPlotter.py`: Similar to `OzoneDataProcessor.py`, but contains functions rather than being a single script. 





Feel free to submit a pull request for new features and/or bug fixes.
