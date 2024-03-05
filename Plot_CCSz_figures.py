import h5py
import hdf5plugin
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Install the packages using the command "pip install -r Python_requirements.txt" in your terminal

# These MZA files were converted from Agilent .d IM-MS data files:
#   For IM spectra, the m/z dimension is stored as indexes (mzbins) of a single full m/z array common for all spectra in the file.
#   IonMobilityBin = 0, represents the total frame spectrum (TFS) or summed spectra in the frame (i.e., ignores IM dimension).

basePath = "E:/mza/CCSz-demo" # <--- Update this path according to the location in your computer
outputPath = "Python_OutputFigures"
mzafiles = ["DT.mza", "DT-linear_CCSz.mza", "SLIM.mza", "SLIM-2nd-order_CCSz.mza", "SLIM-3rd-order_CCSz.mza"]
myYtitles = ["DTIMS AT (ms)", "DTIMS CCS/z", "SLIM AT (ms)", "SLIM 2nd order CCS/z", "SLIM 3rd order CCS/z"]

tics = []
for f in range(0, len(mzafiles)):
    fullPath = os.path.join(basePath, "RawDataMza", mzafiles[f])
    mza = h5py.File(fullPath, 'r')
    print(fullPath)
    #print(list(mza.keys())) # h5py.File acts like a Python dictionary

    # Reading Metadata table:
    metadata = mza["Metadata"]
    # Remove TFS: IonMobilityBin = 0, represents the total frame spectrum (TFS) or summed spectra in the frame (i.e., ignores IM dimension)
    metadata = metadata[metadata["IonMobilityBin"] > 0]

    # Get array of m/z values (common for all spectra in the file) and map mzbins to m/z:
    full_mz = mza["Full_mz_array"][:]

    mzs = []
    ats = []
    points = []

    for k in range(0, metadata.size - 1):

        scan = metadata["Scan"][k]
        # Get the path to the data
        mzapath = str(metadata["MzaPath"][k], 'utf-8')
        # Get the bins of m/z values and their intensities
        mzbins = mza["Arrays_mzbin" + mzapath + "/" + str(scan)][:]
        mz_array = np.array([full_mz[i] for i in mzbins])
        intensity_array = mza["Arrays_intensity" + mzapath + "/" + str(scan)][:]

        for j in range(0, mz_array.size - 1):
            points.append((mz_array[j],
                        intensity_array[j],
                    metadata["IonMobilityTime"][k]))

    # After everything, remember to close the data file
    mza.close()

    df = pd.DataFrame(points, columns=["mz", "intensity", "at"])
    df = df[df["intensity"] > max(df["intensity"])*0.008]
    
    # Scatter plot: --------------------

    fig, ax = plt.subplots()

    ax.scatter(x=df["mz"],
                y=df["at"],
                c=np.log10(df["intensity"]),
                cmap='turbo',
                marker='.',
                s=5,
                linewidths=0.9,
                edgecolors='face')

    # Add label names
    ax.set_xlabel("m/z", fontsize=17)
    ax.set_ylabel(myYtitles[f], fontsize=17)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    plt.savefig(os.path.join(basePath, outputPath, mzafiles[f].replace('.mza', '.png')))
    #plt.show()
    #plt.close()
    if 'CCS' in myYtitles[f]:
        dftic = df.groupby('at')['intensity'].sum().reset_index()
        dftic['intensity'] = dftic['intensity']/max(dftic['intensity'])
        dftic['method'] = myYtitles[f]
        tics.append(dftic)

# Plot with overlaid total ion intensities: --------------

tics = pd.concat(tics)
tics = tics[tics['at'] > 180]

# Define unique colors for each method (can be adjusted to your preference)
colors = {method: color for method, color in zip(tics['method'].unique(), plt.cm.tab10.colors)}

plt.figure()

# Plotting each method's data with its respective color
for method, group in tics.groupby('method'):
    plt.plot(group['at'], group['intensity'], label=method, color=colors[method])

plt.ylabel("Normalized total ion intensity")
plt.xlabel("CCS/z")
plt.legend(title="Method")
plt.savefig(os.path.join(basePath, outputPath, "TICs.png"))
