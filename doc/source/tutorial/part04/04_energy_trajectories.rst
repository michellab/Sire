======================================
Measuring Energies Across Trajectories
======================================

Just as you can make measurements across frames of a trajectory,
you can also calculate energies across frames of a trajectory.

For example, let's now load a trajectory for the ``aladip`` system...

>>> mols = sr.load(sr.expand(sr.tutorial_url, ["ala.top", "ala.traj"]))

If you are in a Jupyter notebook or similar you can view the trajectory
using

>>> mols.view()

.. image:: images/04_04_01.jpg
   :alt: Viewing the aladip trajectory in NGLView

You can play the trajectory by clicking the play button.

There are multiple ways that you can calculate the energy of each
frame of the trajectory. You could do it yourself in a loop, e.g.

>>> for frame in mols.trajectory()[0:5]:
...     print(frame.energy())
-5919.77 kcal mol-1
-5974.75 kcal mol-1
-5928.36 kcal mol-1
-5917.37 kcal mol-1
-5968.4 kcal mol-1

calculates the total energy of the molecules for the first five frames of
the trajectory.

You could also calculate this by using the ``apply`` function, as you
did for measurements.

>>> print(mols.trajectory()[0:5].apply("energy"))
[-5919.77 kcal mol-1, -5974.75 kcal mol-1, -5928.36 kcal mol-1, -5917.37 kcal mol-1, -5968.4 kcal mol-1]

You can use the name of the function, or can supply you own custom function, e.g.

>>> print(mols.trajectory()[0:5].apply(
...                lambda frame: frame[0].energy(frame["water"])) )
[-52.8993 kcal mol-1, -46.7754 kcal mol-1, -45.9996 kcal mol-1, -53.2822 kcal mol-1, -57.4076 kcal mol-1]

calculates the energy between the first molecule and the water molecules
for the first five frames.

The simplest way though is to the use the :func:`~sire.mol.Trajectory.energy`
function of :class:`~sire.mol.Trajectory`.

>>> print(mols.trajectory()[0:5].energy())
   frame  time    1-4_LJ  1-4_coulomb          LJ  ...   dihedral  improper  intra_LJ  intra_coulomb        total
0      0   0.2  3.509838    44.810452  912.765546  ...   9.800343  0.485545 -1.311255     -45.398213 -5919.765132
1      1   0.4  2.700506    47.698455  949.719749  ...  11.776295  1.131481 -1.617496     -48.137251 -5972.683109
2      2   0.6  2.801076    43.486411  845.766625  ...  11.614774  0.124729 -1.103966     -44.458050 -5926.976442
3      3   0.8  3.365638    47.483966  949.699492  ...  11.383852  0.339333 -0.983872     -48.191510 -5917.369422
4      4   1.0  3.534830    48.596027  899.781622  ...  10.214994  0.255331 -1.699613     -48.393884 -5968.440737
[5 rows x 13 columns]

This performs the calculation for each frame, returning the result as a pandas
dataframe. If you are in a Jupyter notebook or similar, you can create
a graph of this using the ``pretty_plot`` function (or the standard pandas
plotting functions), e.g. plotting the energy for every frame gives

>>> mols.trajectory().energy().pretty_plot()

.. image:: images/04_04_02.jpg
   :alt: Graph showing energy of each frame of the aladip trajectory.

Calculating energies of selections
==================================

Just as for single-point energies, you can also calculate energies of
selections across a trajectory.

For example, here is the energy for the first molecule for every frame;

>>> print(mols[0].trajectory().energy())
     frame        time    1-4_LJ  1-4_coulomb      angle      bond   dihedral  improper  intra_LJ  intra_coulomb      total
0        0    0.200000  3.509838    44.810452   7.570059  4.224970   9.800343  0.485545 -1.311255     -45.398213  23.691739
1        1    0.400000  2.700506    47.698455  12.470519  2.785874  11.776295  1.131481 -1.617496     -48.137251  28.808384
2        2    0.600000  2.801076    43.486411  11.607753  2.023439  11.614774  0.124729 -1.103966     -44.458050  26.096167
3        3    0.800000  3.365638    47.483966   6.524609  0.663454  11.383852  0.339333 -0.983872     -48.191511  20.585470
4        4    1.000000  3.534830    48.596027   6.517530  2.190370  10.214994  0.255331 -1.699613     -48.393885  21.215584
..     ...         ...       ...          ...        ...       ...        ...       ...       ...            ...        ...
495    495   99.199997  2.665994    42.866319  11.339087  4.172684   9.875872  0.356887 -1.584092     -44.220004  25.472748
496    496   99.400002  3.062467    44.852774   9.268408  1.878366  10.548897  0.327064 -1.814718     -44.419098  23.704160
497    497   99.599998  3.530233    44.908117  10.487378  4.454670  10.223964  1.006034 -0.692972     -44.902059  29.015366
498    498   99.800003  3.511116    42.976288   9.017446  0.809064  10.841436  0.518190 -1.862433     -43.205031  22.606076
499    499  100.000000  3.768998    41.625135  13.629923  1.089916  11.889372  0.846805 -1.897328     -44.306436  26.646386
[500 rows x 11 columns]

while here is the interaction energy between the first molecule and
all of the water molecules

>>> print(mols[0].trajectory().energy(mols["water"]))
     frame        time         LJ    coulomb      total
0        0    0.200000 -16.737847 -36.161409 -52.899257
1        1    0.400000  -9.298626 -37.461305 -46.759931
2        2    0.600000 -18.198259 -27.847419 -46.045678
3        3    0.800000 -17.793692 -35.488619 -53.282310
4        4    1.000000 -12.810445 -44.586756 -57.397201
..     ...         ...        ...        ...        ...
495    495   99.199997 -11.317037 -36.419693 -47.736730
496    496   99.400002 -14.172688 -29.708265 -43.880953
497    497   99.599998  -7.688854 -36.490591 -44.179445
498    498   99.800003 -12.603816 -35.767549 -48.371365
499    499  100.000000  -6.053596 -42.322420 -48.376016
[500 rows x 5 columns]

and here is the energy between the carbon atoms in the first molecule
and all of the non-carbon atoms in the first molecule

>>> print(mols[0]["element C"].trajectory().energy(mols[0]["not element C"]))
     frame        time    1-4_LJ  1-4_coulomb      angle      bond   dihedral  improper  intra_LJ  intra_coulomb      total
0        0    0.200000  1.424470    54.597537   7.512649  4.194598   9.800343  0.485545 -0.844606     -15.571905  61.598631
1        1    0.400000  0.791583    52.969358  12.008645  1.821351  11.776295  1.131481 -0.937170     -17.901480  61.660062
2        2    0.600000  1.533581    53.594240  11.531951  1.623215  11.614774  0.124729 -0.885362     -16.394078  62.743051
3        3    0.800000  1.644253    53.721723   6.478095  0.614500  11.383852  0.339333 -0.950609     -19.388795  53.842352
4        4    1.000000  1.609999    53.472406   6.515561  1.965713  10.214994  0.255331 -1.010526     -21.909819  51.113659
..     ...         ...       ...          ...        ...       ...        ...       ...       ...            ...        ...
495    495   99.199997  1.324959    51.009025  11.173632  3.007010   9.875872  0.356887 -0.907750     -17.684035  58.155600
496    496   99.400002  1.310871    51.685638   8.855937  1.300661  10.548897  0.327064 -0.896125     -22.184101  50.948842
497    497   99.599998  2.073824    51.499993  10.461966  2.345364  10.223964  1.006034 -0.867677     -20.452924  56.290544
498    498   99.800003  1.357558    50.525535   8.838606  0.306832  10.841436  0.518190 -0.940949     -18.327720  53.119488
499    499  100.000000  1.783374    50.771940  13.388107  0.649792  11.889372  0.846805 -0.883830     -19.271354  59.174206
[500 rows x 11 columns]

Calculating energies for each view in a selection
=================================================

Similarly, you can use the :func:`~sire.mol.Trajectory.energies` function
to calculate the energies of each view in a selection. For example,
to calculate the energy between each residue in the first molecule and all
the water molecules use

>>> print(mols[0].residues().trajectory().energies(mols["water"]))
     frame        time  LJ(ACE:1)  LJ(ALA:2)  ...  coulomb(NME:3)  total(ACE:1)  total(ALA:2)  total(NME:3)
0        0    0.200000  -6.223427  -5.550336  ...       -5.258177    -16.052065    -26.624932    -10.222261
1        1    0.400000  -3.207325  -3.958487  ...       -7.924708    -12.427712    -24.274698    -10.057524
2        2    0.600000  -5.751530  -6.979127  ...       -3.016906    -15.462394    -22.098776     -8.484508
3        3    0.800000  -5.575469  -6.994782  ...       -4.398916    -17.252672    -26.407282     -9.622357
4        4    1.000000  -4.911498  -3.529181  ...       -6.049042    -19.537946    -27.440446    -10.418809
..     ...         ...        ...        ...  ...             ...           ...           ...           ...
495    495   99.199997  -2.556662  -7.445725  ...       -6.133000    -20.918159    -19.370920     -7.447650
496    496   99.400002  -5.112405  -5.134994  ...       -2.569671    -18.692191    -18.693802     -6.494959
497    497   99.599998  -2.955196  -3.097448  ...       -4.703180    -18.708102    -19.131953     -6.339391
498    498   99.800003  -3.724783  -6.256671  ...       -4.199227    -18.245600    -23.304175     -6.821591
499    499  100.000000  -4.956662  -1.798866  ...       -6.793895    -22.883995    -19.400058     -6.091963
[500 rows x 11 columns]

or to calculate the energy between each individual carbon atom in the
first molecule, and all non-carbon atoms, use

>>> print(mols[0].atoms("element C").trajectory().energies(mols[0]["not element C"]))
     frame        time  1-4_LJ(C:15)  1-4_LJ(C:5)  ...  total(CA:9)  total(CB:11)  total(CH3:19)  total(CH3:2)
0        0    0.200000      0.172479     0.145310  ...     0.734309     10.422320      12.830321     -7.916957
1        1    0.400000      0.068533     0.293878  ...     0.833780     10.329925      12.680069     -5.459349
2        2    0.600000      0.277433     0.119361  ...     2.105670      9.903695      12.342441     -8.597333
3        3    0.800000      0.218850     0.293046  ...     0.648526      8.335720      11.774262     -5.820661
4        4    1.000000      0.272544     0.131482  ...     1.613490      8.462497      11.444502     -5.993019
..     ...         ...           ...          ...  ...          ...           ...            ...           ...
495    495   99.199997      0.084804     0.092928  ...     2.501087      8.260509      13.692395     -8.903780
496    496   99.400002      0.301822     0.055070  ...     0.860013      7.782723      12.373809     -9.595075
497    497   99.599998      0.287670     0.167273  ...     0.810351      8.201226      12.855330    -10.161720
498    498   99.800003      0.247033     0.071368  ...     0.636738      8.162127      12.441049    -11.020026
499    499  100.000000      0.762106     0.001875  ...     0.772173      8.378594      12.675099     -9.996188
[500 rows x 46 columns]

or you can calculate the energy between the first two residues of the
first molecule, and the last residue of that molecule,

>>> print(mols[0].residues()[0:2].trajectory().energies(mols[0].residues()[-1]))
     frame        time  1-4_LJ(ALA:2)  ...  intra_coulomb(ALA:2)  total(ACE:1)  total(ALA:2)
0        0    0.200000       0.804882  ...            -21.180433     -0.623493     11.417398
1        1    0.400000       0.345863  ...            -19.744563     -0.887233     11.519219
2        2    0.600000       0.570195  ...            -20.630107     -0.486087     10.305274
3        3    0.800000       0.583037  ...            -20.152031     -0.744566     11.144469
4        4    1.000000       0.865922  ...            -20.661027     -0.736126     11.418229
..     ...         ...            ...  ...                   ...           ...           ...
495    495   99.199997       1.541341  ...            -25.465191      0.631658      8.578779
496    496   99.400002       0.881839  ...            -24.581280      0.339605     10.761385
497    497   99.599998       1.974929  ...            -25.758041      0.272158     12.628546
498    498   99.800003       1.860815  ...            -25.768444      0.550181      9.597360
499    499  100.000000       1.335519  ...            -25.910358      0.396877      8.686331
[500 rows x 17 columns]

As before, all of the above are returned as pandas dataframes, and so can be
plotted using the standard pandas functions, or via the ``pretty_plot`` function
that is added by sire when the energies are calculated, e.g.

>>> mols[0].residues()[0:2].trajectory().energies(mols[0].residues()[-1]).pretty_plot()

.. image:: images/04_04_03.jpg
   :alt: Graph showing the energy between the first two and last residues of the first molecule from each frame of the aladip trajectory.

Note how the column names for these dataframes include both the energy component
and an identifier for the view. You can programmatically generate these
column names using the :func:`sire.colname` function e.g.

>>> print(sr.colname(mols[0].residues()[0], "total"))
total(ACE:1)

gives the column name for the total energy component for the first
residue.

Calculating energies where the selection changes per frame
==========================================================

The :func:`~sire.mol.Trajectory.energy` and :func:`~sire.mol.Trajectory.energies`
functions calculate energies using static selections for each frame. This
means that the selection is made before calculating the energy, and then
this selection is used for each energy calculation.

This is a problem if the selection is based on something that would change
in each frame, e.g. for example the coordinates of the atoms.
For example, here we try to calculate the energy between the first molecule
and all water molecules that are within 5 Å.

>>> print(mols[0].trajectory().energy(mols["water and molecule within 5 of molidx 0"]))
     frame        time         LJ    coulomb      total
0        0    0.200000 -13.634133 -33.715126 -47.349259
1        1    0.400000  -6.016929 -36.199230 -42.216159
2        2    0.600000 -14.898939 -27.797556 -42.696495
3        3    0.800000 -14.223504 -34.152557 -48.376061
4        4    1.000000  -9.138862 -43.100886 -52.239748
..     ...         ...        ...        ...        ...
495    495   99.199997  -0.078865  -3.048407  -3.127272
496    496   99.400002  -1.594843  -0.130453  -1.725296
497    497   99.599998  -1.258401  -3.576412  -4.834813
498    498   99.800003  -1.418459  -4.838913  -6.257371
499    499  100.000000  -2.376892  -1.286287  -3.663179
[500 rows x 5 columns]

This clearly hasn't done what we want, because the energy trends towards
zero over time. This is because the selection of water molecules is fixed,
and was not updated for every frame. The selected water molecules diffused
away, and so the energy trended towards zero.

You need to manually update the selection if you want it to update
for every frame. You could do this via a loop, e.g.

>>> for frame in mols.trajectory():
...     print(frame[0].energy(frame["water and molecule within 5 of molidx 0"]))
-47.3493 kcal mol-1
-43.1736 kcal mol-1
-44.3147 kcal mol-1
-49.283 kcal mol-1
-53.6223 kcal mol-1
...
-40.8416 kcal mol-1
-38.6555 kcal mol-1
-38.1635 kcal mol-1
-42.4407 kcal mol-1
-42.367 kcal mol-1

or by using ``apply``

>>> print(mols.trajectory().apply(
...          lambda frame: frame[0].energy(frame["water and molecule within 5 of molidx 0"])) )
[-47.3493 kcal mol-1, -43.1736 kcal mol-1, -44.3147 kcal mol-1, -49.283 kcal mol-1,
...
 -38.6555 kcal mol-1, -38.1635 kcal mol-1, -42.4407 kcal mol-1, -42.367 kcal mol-1]


