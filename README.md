# SoundingDataPlotter (v0.1)
Group of Python scripts to plot UWYO and processed OWLETS-2 sounding data.

WARNING: These scripts were made by an amatuer wannabe data science poser who only slightly knows what he is doing. Use at your own risk...

These scripts are meant to be run in the same directory. The hierarchy as shown here is what you should emulate on your device for them to work.

## File Descriptions
These files are from the OWLETS-2 campaign run in the summer of 2018.

'Helpers' - Contains the helper files used to create the visualizations and process the data.

**Full Sounding Data**

`Raw Soundings` - All of the sounding data, separated by each launch with columns holding a different measurement variable. The data between each sounding may vary due to a different type of model used, so be sure to check the "Header" column for more details.

**Processed Data**

`Derived MR` - These are the soundings where I derived the mixing ratio using the Arden Buck Eqaution

'Measured MR' - Soundings where the mixing ratio was already included

In my processed data folders, there are 7 columns:

O3 - Ozone; Displayed in parts per billion by volume (ppbv)

Theta_K - Potential Temperature; Displayed in Kelvin (K)

RH - Relative Humidity; Displayed in %

q - Specific Humidity; Displayed in g/kg

N - Atmospheric Refractivity; (unitless)

Alt - Altitude; Displayed in meters (m)

Note that 

