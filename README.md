
This R package contains code supporting the analysis of [facial shape and motion](http://people.bath.ac.uk/jjf23/face/index.html)
by [Julian Faraway](http://people.bath.ac.uk/jjf23/)

You can install it (provided you have loaded the `devtools` package) with:

```
install_github("julianfaraway/brinla")
```

The functions available are:

```
arraymean                           Array mean
asymmetry                           Asymmetry scores
averecface                          Average faces and rotate so head is upright
bridgefix                           Location shift onto template
cleansm                             Process a facial motion to fill-in missing values and report
                                    diagnostics
constructfaceseq                    Array of faces constructed from geodesic and profile information
constructscorefaceseq               Array of faces constructed from geodesic and score informationTitle
extmfconfig                         Extract information about the motion capture file
faceGPA                             Generalized Procrustes Analysis
faceOPA                             Ordinary Procrustes Analysis
faceframe                           Construct a single frame suitable for movie construction
facemovie                           Make a facial motion movie
facenorm                            Compute size of face
first_PCA_traj                      extract the first PCA, select start and max, return initial and max
                                    frames
misfixlin                           Fill in missing values using linear interpolation
multifaceplot                       Plot an array of shapes with a different color for each marker
parproc                             Partial Procrustes rotation
plotface                            Plot a face or a pair of faces
read_face_file                      Read in a face file
readlmk                             Read in an LMK file
reflectface                         Reflect face
rotface                             Rotate a face matrix in SD
```

The package is under development and lacks full documentation or testing. If you are interested in using it,
please contact me.

