# CalibSCE

This code package produces a 3D SCE spatial calibration for MicroBooNE using a set of near-crossing t0-tagged cosmic muon tracks from data or MC (CalibSCE.cpp).  Addtional applications are included to smooth the resultant 3D spatial distortion maps (MakeSmoothHists.cpp), produce associated smooth E field maps (MakeSmoothHists.cpp), and extract systematic biases on these measurements using a set of laser tracks from the MicroBooNE UV laser system (MakeSystVar.cpp).

First one must obtain input cosmic muon track files (data and MC) from the FNAL machines and place them in the "data" folder (or change the file path in the application code if running on the FNAL machines).  These files are currently located here:

```/uboone/data/users/mrmooney/MC_Cosmics.root```

```/uboone/data/users/mrmooney/Data_Run1_EXTBNB.root```

Another preliminary step is building the code to produce the calibration applications:

```$ make```

To build the initial 3D SCE spatial calibration maps, run the following command:

```$ ./CalibSCE```

Aside from file names/paths at the top of the code, you may need to change "isMC" to true/false if you want to run on MC/data.  Also "numCalibTracks" should be changed if you want to run on a different number of tracks (nominal calibration used 200k tracks, but this could take many hours to run).  Finally, "VelRatio" should be changed from 0.992 to 1.0 for data in the "main" function if more modern MC is used including 1.098 m/ms drift velocity in the simulation/reconstruction chain, as older MicroBooNE reconstruction used 1.114 m/ms leading to the need to apply an additional correction factor to data first in the original SCE calibration.  Also "ShiftedCathode" should have a factor that changes from 2.5818 to 2.5524 for data in that case.

Note the "output_siminterp_MicroBooNE_4p5_gap.root" file that is included in the "data" folder, which includes the cathode correction to be used prior to the crossing-track calibration.  This can be reproduced by following the logic of the MicroBooNE cosmics-based SCE paper and interpolating the 2D face distortion maps across the cathode.

To produce smooth spatial distortion maps and associated E field maps from the output of the previous step, run the following command:

```$ ./MakeSmoothHists```

Again, note that "isMC" may need to be changed to true/false if you want to run on MC/data, and file names/paths may require changing as well.

Finally, to produce the systematic bias maps, run the following command:

```$ ./MakeSystVar```

Once again, note that "isMC" may need to be changed to true/false if you want to run on MC/data, and file names/paths may require changing.

This code makes use of interpolation tools included in the Fade3D library; note that the makefile (GNUmakefile) will likely need to be modified to point to the location of your local Fade3D installation files.  More information here:  http://www.geom.at, bkorn@geom.at

For any questions about this code, please contact Mike Mooney:  mrmooney@colostate.edu