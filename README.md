
<center>
<img src="man/figures/Lipic_Logo.png" width=350 />
</center>

# Getting Started

In lipidomics shotgun and HILIC-MS methods, monoisotopologue (M+0) peak
can be used for quantitative purposes. However, intrinsic M+0 intensity
differences should be taken into account when ions with different
isotopic pattern distribution are considered. This is commonly referred
to as the “Type I” isotopic effect. Particularly, it becomes noticeable
as the carbon atom difference increases. Furthermore, isotopic pattern
overlap is another major issue in lipid class-related mass spectra.
Indeed, lipids differing by only the number of C=C bonds frequently
occur within the same class. Here, the monoisotopologue M+0 peak
intensity is artificially enhanced by M+2 and M+4 isotopologues of
lipids with one or two additional double bonds, respectively. This kind
of interference is commonly known as either the “double bond overlapping
effect” or the “Type II” effect. Moreover, Type II correction should
account for complications in High Resolution mass spectra. These are due
to the misalignment between the M+0 and the interfering isotopologue
peaks. As a consequence, Type II correction can no longer be performed
by simply subtracting the expected peak intensities of the M+2 and M+4
interfering isotopologues to that of the convoluted M+0 peak. The
**`LIPIC`** package has been specifically designed for performing both
Type I and Type II corrections in lipid class high resolution mass
spectra. Nevertheless, the correction can be applied to lower resolution
spectra as well.

### Dependencies on other R Packages

LIPIC is based on functions from other R packages, namely
[enviPat](https://cran.rstudio.com/web/packages/enviPat/index.html) for
isotope fine structure calculation and
[plotly](https://cran.r-project.org/web/packages/plotly/index.html) for
the interactive graphical output. Hence, the user needs to install the
aforementioned R packages using the following code lines:

``` r
install.packages("enviPat")
install.packages("plotly")
```

Moreover, the `pick.peaks` function (from the ChemometricsWithR package)
has been added to the LIPIC package.

### Install LIPIC

The user can install the development version of `LIPIC` from GitHub. For
this purpose, the `devtools` package must have been formerly installed.

``` r
install.packages("devtools")
devtools::install_github("AnCaste/LIPIC")
```

### Uninstall LIPIC

If no longer needed, the LIPIC package can be removed in any moment
using the following code line:

``` r
remove.packages("LIPIC")
```

### Load LIPIC

Prior to further operations, the content of the Lipic Package can be
loaded into the R environment using the `library()` function, as shown
below.

``` r
library("LIPIC")
```

## LIPIC Package main Functions: Lipic and Lipic\_It

Both Type I and Type II corrections can be performed using either the
`Lipic_It` or the `Lipic` functions. The choice of a function rather
than the other depends upon the user’s preliminary information. Indeed,
the reliability of the `LIPIC` correction workflow depends on **`dmz`**
values. Here, the “dmz” symbol refers to the absolute error in the
measurement of each M+0 peak *m/z*. Specifically, `dmz` is obtained by
the difference between the measured M+0 peak *m/z* and its theoretical
value.

Typically, `dmz` values are initially unknown. In this case, the user
should refer to the `Lipic_It` function. `Lipic_It` can calculate an
estimate of `dmz` for each M+0 peak using a recursive approach. When
convergence is reached, the `dmz` values referred to the last iteration
are used to build the output.

Conversely, if the user already retrieved information about `dmz`, he
could use the `Lipic` function specifying the corresponding values in
the input data frame, as described in the next section.

### Lipic and Lipic\_It: Input Arguments

`Lipic` and `Lipic_It` share the same input arguments, namely the
**input data frame** (*data*), the lipid ions **charge state** (*z*, the
sign must be specified as well) and the **operating resolving power of
the mass spectrometer** (*Res*). The values of the remaining parameters
(i.e. *span*, *ppm* and *xseq*) are set as default. Check `Lipic` and
`Lipic_It` help pages for further details.

``` r
Lipic_It(data,z,Res,span=25,ppm=10,xseq=0.001)
## A typical code line to run the Lipic_It function. The input data frame,z and Res need to be specified by the user.

Lipic(data,z,Res,span=25,ppm=10,xseq=0.001)
## A typical code line to run the Lipic function. The input data frame,z and Res need to be specified by the user.
```

#### The Input Data Frame

This section is focused on the **input data frame** structure and how to
build it properly. The `LIPIC` package offers an example of input data
frame, namely `CLmix`. It contains specific information for all the
identified the sum compositions attributed to the corresponding
\[M-2H\]<sup>2−</sup> ions in the direct infusion ESI(-)-MS spectrum of
a Cardiolipin mix purchased from Avanti Polar Lipids. Such spectrum is
embedded in the LIPIC package as a data frame named `CLmixSpectrum`.
Both `CLmix` and `CLmixSpectrum` can be loaded into the R environment
using the `data()` function as shown below.

``` r
data("CLmix")
## Load the CLmix data frame.

data("CLmixSpectrum")
## Load the spectrum of the Cardiolipin mix.
```

The `CLmixSpectrum` can be exported as a Tab Separated Value .txt file
using the `write.table()` function.

``` r
write.table(CLMixSpectrum,"C:/Documents/Spectra/CLMixSpectrum.txt",sep="\t",row.names=FALSE)
## The file path ("C:/Documents/Spectra/CLMixSpectrum.txt") is illustrative only. The you need to replace it with the file path leading to the position in which you want the .txt to be saved.
```

Let us now focus on the data frame content. The first two columns must
be filled with the sum composition and the chemical formula for all the
identified **lipids ions** in the mass spectrum. The remaining three
columns refer to the corresponding M+0 peak parameters, i.e the m/z, dmz
and intensity values, respectively.

**N.B.** The elements in the input data frame must be listed following
the ascending order of the “measured *m/z*” values. This is a mandatory
requirement for LIPIC functions to perform a reliable correction.

As shown below, when dmz values are not known, the dmz column should be
filled with zeros. Than, the `Lipic_It` function has to be employed for
Type I and Type II corrections. As will be further pointed out,
`Lipic_It` will return an estimate for each dmz.

| Sum Composition |   Formula    | Measured m/z | dmz | Intensity |
|:---------------:|:------------:|:------------:|:---:|:---------:|
|     CL 70:8     | C79H136O17P2 |   709.4623   |  0  |  1016080  |
|     CL 70:7     | C79H138O17P2 |   710.4693   |  0  | 12176987  |
|     CL 70:6     | C79H140O17P2 |   711.4784   |  0  | 23126606  |
|     CL 70:5     | C79H142O17P2 |   712.4829   |  0  | 10350815  |
|     CL 70:4     | C79H144O17P2 |   713.4917   |  0  |  2470070  |
|     CL 71:7     | C80H140O17P2 |   717.4781   |  0  |  2258730  |
|     CL 71:6     | C80H142O17P2 |   718.4862   |  0  |  2275952  |
|     CL 71:5     | C80H144O17P2 |   719.4941   |  0  |  1831880  |
|     CL 71:4     | C80H146O17P2 |   720.4996   |  0  |  1241680  |
|    CL 72:10     | C81H136O17P2 |   721.4612   |  0  |  1009175  |
|     CL 72:9     | C81H138O17P2 |   722.4708   |  0  | 21706585  |
|     CL 72:8     | C81H140O17P2 |   723.4779   |  0  | 526420104 |
|     CL 72:7     | C81H142O17P2 |   724.4802   |  0  | 291071934 |
|     CL 72:6     | C81H144O17P2 |   725.4915   |  0  | 79097022  |
|     CL 72:5     | C81H146O17P2 |   726.4956   |  0  | 28824140  |
|     CL 72:4     | C81H148O17P2 |   727.5041   |  0  |  4486245  |
|    CL 74:10     | C83H140O17P2 |   735.4760   |  0  | 11987946  |
|     CL 74:9     | C83H142O17P2 |   736.4820   |  0  | 20176370  |
|     CL 74:8     | C83H144O17P2 |   737.4901   |  0  | 15311198  |
|     CL 74:7     | C83H146O17P2 |   738.4978   |  0  |  8180636  |
|     CL 74:6     | C83H148O17P2 |   739.5029   |  0  |  3781301  |

Conversely, if the user had preliminary information about dmz values, he
could use them to compile the fourth column of the input data frame. In
such case, the `Lipic` function can be used for both the isotopic
corrections.

#### The VisualSpec Function as a Tool for Building the Input Data Frame

The user can easily retrieve the *m/z* and intensity values of each M+0
peak using the `VisualSpec` function included in the LIPIC package.
Specifically, `VisualSpec` performs a peak peaking operation powered by
the `pick.peaks`.

`pick.peaks` belongs to the ChemometricsWithR package, but has been
included in LIPIC.

`VisualSpec` returns a data frame in which the *m/z* and intensity
values are listed for all the identified peaks. Furthermore, it
generates a `plotly` interactive version of the mass spectrum.

As shown below, `VisualSpec` requires in input the file path to the Tab
Separated Value .txt file of the mass spectrum. The remaining parameters
(i.e. mzmin and mzmax) are set to 0. These default values should not be
changed if the full *m/z* range needs to be considered. Otherwise, if
the user is interested in a part of the mass spectrum, he should
indicate the corresponding boundary values in order to focus the
visualization and the peak picking operation to a specific spectral
region. Check the `VisualSpec` help page for further information.

``` r
VisualSpec=function(path,mzmin=0,mzmax=0,span=25)
```

The application of `VisualSpec` to the Cardiolipin Mix is shown below.

``` r
VisualSpec("C:/Documents/Spectra/CLMixSpectrum.txt")
## The file path ("C:/Documents/Spectra/CLMixSpectrum.txt") is illustrative only. You need to replace it with the file path in which your .txt file is stored.
```

Click [here](https://chart-studio.plotly.com/~AnCaste/3/#/) to see the
interactive output.

#### The Output Data Frame

The **output data frame** (ODF) returned by both `Lipic` and `Lipic_It`
can be considered as an extension of the input data frame. Particularly,
when `Lipic_It` is employed, the `dmz` column of the ODF is
automatically compiled by the `dmz` values calculated by the function
itself. Conversely, when `Lipic` is used, the `dmz` column shows the
same `dmz` values indicated by the user in the input data frame.

In particular, when the `CLmix` dataset is used as input of the
`Lipic_It` function, the following output data frame is obtained:

``` r
ODF=Lipic_It(CLmix,-2,75000)
## Here, the output data frame generated by the Lipic_It function is assigned to the "ODF" variable.
```

| Sum Composition |   Formula    | Measured m/z |     dmz      | Measured Intensity |     Type II     | Type II & Type I | Lipid Class Adj-RA (%) | Lipid Class UnAdj-RA (%) | Most Abundant Lipid Adj-RA (%) | Most Abundant Lipid UnAdj-RA (%) |
|:---------------:|:------------:|:------------:|:------------:|:------------------:|:---------------:|:----------------:|:----------------------:|:------------------------:|:------------------------------:|:--------------------------------:|
|     CL 70:8     | C79H136O17P2 |  709.462338  | -0.000848516 |    1016079.965     |  1016079.9650   |   2516235.426    |      0.1252130837      |       0.0929741276       |          0.1897999332          |           0.1888221388           |
|     CL 70:7     | C79H138O17P2 |  710.469288  | -0.001898516 |    12176987.420    |  11974702.3341  |   29661138.020   |      1.4759996297      |       1.1144838519       |          2.2373431188          |           2.2634170382           |
|     CL 70:6     | C79H140O17P2 |  711.478384  | 0.000197484  |    23126605.970    |  20255145.5254  |   50183176.551   |      2.4972187499      |       2.1171202257       |          3.7853228915          |           4.2996818509           |
|     CL 70:5     | C79H142O17P2 |  712.482868  | -0.000318516 |    10350815.020    |  3061781.2415   |   7587463.842    |      0.3775679076      |       0.9477807060       |          0.5723232872          |           1.9248578567           |
|     CL 70:4     | C79H144O17P2 |  713.491739  | 0.000552484  |    2470070.396     |  1413423.4277   |   3503438.164    |      0.1743383355      |       0.2262259087       |          0.2642647515          |           0.4594445898           |
|     CL 71:7     | C80H140O17P2 |  717.478096  | -0.001090516 |    2258730.451     |  2267125.7973   |   5677667.645    |      0.2825324954      |       0.2090112171       |          0.4282671362          |           0.4244830907           |
|     CL 71:6     | C80H142O17P2 |  718.486184  | 0.000997484  |    2275952.222     |  1655627.1830   |   4147215.370    |      0.2063740220      |       0.2106531810       |          0.3128249416          |           0.4278177725           |
|     CL 71:5     | C80H144O17P2 |  719.494110  | 0.000923484  |    1831880.435     |  1373308.2461   |   3440818.470    |      0.1712222499      |       0.1695905689       |          0.2595413406          |           0.3444232794           |
|     CL 71:4     | C80H146O17P2 |  720.499573  | -0.000613516 |    1241680.246     |   845610.9472   |   2119161.328    |      0.1054538546      |       0.1149778107       |          0.1598485875          |           0.2335096515           |
|    CL 72:10     | C81H136O17P2 |  721.461196  | -0.001990516 |    1009174.578     |  1009174.5780   |   2553482.050    |      0.1270665528      |       0.0943503790       |          0.1926094505          |           0.1916171823           |
|     CL 72:9     | C81H138O17P2 |  722.470787  | -0.000399516 |    21706585.310    |  21499157.4820  |   54411117.049   |      2.7076098214      |       2.0298715193       |          4.1042369390          |           4.1224875305           |
|     CL 72:8     | C81H140O17P2 |  723.477884  | -0.001302516 |   526420104\.100   | 523708139\.9260 |  1325730406.349  |     65.9710894284      |      49.2389971904       |         100.0000000000         |          100.0000000000          |
|     CL 72:7     | C81H142O17P2 |  724.480154  | -0.003032516 |   291071933\.900   | 102823948\.0437 |  260351364\.829  |     12.9556228700      |      27.2318254486       |         19.6383339767          |          55.3054022268           |
|     CL 72:6     | C81H144O17P2 |  725.491457  | 0.000270484  |    79097022.300    |  40463772.8341  |  102478237\.296  |      5.0995292291      |       7.4017812630       |          7.7299454554          |          15.0323558263           |
|     CL 72:5     | C81H146O17P2 |  726.495593  | -0.003593516 |    28824140.490    |  14617892.8089  |   37029660.761   |      1.8426725750      |       2.6979392478       |          2.7931516532          |           5.4792733438           |
|     CL 72:4     | C81H148O17P2 |  727.504142  | -0.001044516 |    4486244.889     |  -5596529.7140  |      0.000       |      0.0000000000      |       0.4200088707       |          0.0000000000          |           0.8530004563           |
|    CL 74:10     | C83H140O17P2 |  735.476037  | -0.003149516 |    11987945.850    |  11987945.8500  |   31006567.124   |      1.5429509671      |       1.1456831502       |          2.3388289938          |           2.3267800231           |
|     CL 74:9     | C83H142O17P2 |  736.482025  | -0.004161516 |    20176370.440    |  16883633.3425  |   43679183.405   |      2.1735665870      |       1.9286902323       |          3.2947259259          |           3.9169973850           |
|     CL 74:8     | C83H144O17P2 |  737.490078  | -0.002108516 |    15311198.400    |  9947342.6919   |   25740503.697   |      1.2809007497      |       1.4639626100       |          1.9416092121          |           2.9731771432           |
|     CL 74:7     | C83H146O17P2 |  738.497823  | -0.002363516 |    8180636.426     |  4761241.0078   |   12323379.741   |      0.6132368867      |       0.7823617279       |          0.9295539789          |           1.5889067052           |
|     CL 74:6     | C83H148O17P2 |  739.502867  | -0.004319516 |    3781300.769     |  2094541.0227   |   5422483.508    |      0.2698340045      |       0.3617107631       |          0.4090185668          |           0.7346022132           |

If compared to the input data frame, the ODF shows 6 additional columns.
Among these, the content of the 8th column (i.e. **“lipid class
Adj-RA(%)”**) is particularly important for quantitative purposes.
Specifically, the term “lipid class (%)RA” refers to the percentage
ratio between the M+0 peak intensity and the sum of all the M+0
intensities detected for all the sum composition identified for a
specific lipid class. Here, the *“Adj”* notation indicates that the
considered intensities have been corrected for both Type I and Type II
isotopic effect. If equal ionization yields are assumed, the “lipid
class Adj-RA(%)” values can be considered an estimate of the percentage
molar fraction of the corresponding sum compositions. On the other hand,
the *“UnAdj-RA(%)”* notation in **“lipid class UnAdj-RA(%)”** (9th
column) indicates that such relative abundances are calculated when only
Type I correction is performed on measured intensities. Therefore, the
comparison between the values in the 8th and 9th column allows the user
to perceive the extent of Type II correction. Moreover, such comparison
is graphically enhanced by the **interactive bar chart** output powered
by `plotly`. Click [here](https://chart-studio.plotly.com/~AnCaste/7)
for a preview on the bar chart created by `Lipic_It` when using the
`CLmix` as input data frame.

Let us focus on the 10th and 11th columns. Here, the term “most abundant
lipid (%)RA” refers to the percentage ratio between the M+0 peak
intensity and highest M+0 peak intensities detected among all the sum
composition. Again, when the *“Adj”* notation is used, these intensities
are corrected for both Type I and Type II isotope effect. Conversely,
the *“UnAdj”* term indicates that only Type I correction has been
performed on measured intensities.

The user can check the effects of isotope correction by comparing the
values in the 5th, 6th and 7th columns. Specifically, the 5th column
(“Measured Intensity”) shows the M+0 peak measured intensities
(i.e. those indicated in the input data frame). The 6th shows the
corresponding values after Type II correction, whereas the 7th column
reports both Type I and Type II correction results. In particular, the
values in the 7th column correspond to a reliable estimate of the
overall intensity distributed over the whole isotopic pattern. The user
can refer to them for quantitative purposes.

Both `Lipic` and `Lipic_It` offer another validation tool. If required
by the user, an interactive spectra comparison is generated in the R
Studio viewer pane.

``` r
ODF=Lipic_It(CLmix,-2,75000)
Do you wish to compare the real spectrum with the simulated ones? (yes/no): yes
## Type yes to ask for spectra comparison.

Please, paste .txt mass spectrum file path: "C:/Documents/Spectra/CLMixSpectrum.txt" 
## The file path ("C:/Documents/Spectra/CLMixSpectrum.txt") is illustrative only. You need to replace it with the file path lin which your spectrum .txt file is stored.
```

The interactive plot displays the overlap between the real mass spectrum
(green line) and two simulated mass spectra (blue and red lines). The
**corrected simulated spectrum** (blue line) is calculated considering
the effect of both Type I and Type II corrections. Conversely, in the
**uncorrected simulated spectrum** (red line) the consequences of Type
II isotopic effect are neglected. The user can choose which spectra to
be shown in the overlaid visualization by selecting them using the
interactive legend. You can do it by clicking on the coloured line and
making the corresponding spectrum disappear.
[Here](https://chart-studio.plotly.com/~AnCaste/11/#/) you can practice
with the interactive spectra comparison generated by `Lipic_It` for the
main ion clusters detected in the Cardiolipin mix mass spectrum.

**N.B.** The simulated spectra are calculated starting from theoretical
*m/z* values. As a consequence, they appear shifted in respect to the
corresponding peaks in the real spectrum. Thus, the user should focus on
peak height and shape for assessing the correction reliability.

## References

For LIPIC theoretical background see:

-   A. Castellaneta, I. Losito, D. Coniglio, B. Leoni, P. Santamaria, M.
    A. Di Noia, L. Palmieri, C. D. Calvano, T.R.I. Cataldi: “*LIPIC: an
    automated workflow to account for isotopologues-related
    interferences in electrospray ionization high resolution mass
    spectra of phospholipids*”, Journal of the American Society for Mass
    Spectrometry, 2021. DOI: <https://doi.org/10.1021/jasms.1c00008>.

For packages on which LIPIC is based see:

-   M. Loos, C. Gerber, F. Corona, J. Hollender, H. Singer: *Accelerated
    Isotope Fine Structure Calculation Using Pruned Transition Trees*.
    Anal. Chem. 87, 5738–5744 (2015). DOI:
    <https://doi.org/10.1021/acs.analchem.5b00941>.

-   Plotly Technologies Inc.: Collaborative data science,
    <https://plotly.com/> (2015).

-   R. Wehrens: *Chemometrics withR*. Springer Verlag, Berlin Heidelberg
    (2020). DOI: <https://doi.org/10.1007/978-3-642-17841-2>.

## Acknowledgments

Thanks to Prof. Ilario Losito, Dr. Davide Coniglio, Prof. Cosima Damiana
Calvano and Prof. Tommaso Cataldi for their contribution to the
development of the LIPIC theoretical background. Prof. Pietro
Santamaria, Dott. Beniamino Leoni, Prof. Luigi Palmieri and Dr. Maria
Antonietta Di Noia are acknowledged for providing samples employed to
obtain data for LIPIC testing.
