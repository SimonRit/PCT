# PCT list-mode data format

PCT uses its own data format to store list-mode proton CT data as [MetaImage](https://itk.org/Wiki/ITK/MetaIO/Documentation) files (`.mhd` or `.mha` extension). The data corresponds to so-called proton pairs information (position, direction, energy, ...), i.e. the information recorded by the entrance and exit detectors which are paired, e.g.  by a coincidence electronics or by an algorithm.

Some PCT executables produce such data, e.g. `pctpairprotons` from [ROOT](https://root.cern/) data produced by a [GATE](https://github.com/OpenGATE/opengate) simulation, whereas other ones take such data as input, e.g. `pctbinning`. PCT uses an [ITK](https://itk.org/) image internally (`itk::Image`) to read from / write to these [MetaImage](https://itk.org/Wiki/ITK/MetaIO/Documentation) files.

Proton pairs are stored in 2D images of 3D float vectors, in which the first dimension is a series of 5 or 6 vectors of 3 floats each, one per proton pair, and the second dimension corresponds to the number of proton pairs. Each vector stores the following data:

1. $(u_{\text{in}},v_{\text{in}},w_{\text{in}})$ position of the proton at the entrance detector;
2. $(u_{\text{out}},v_{\text{out}},w_{\text{out}})$ position of the proton at the exit detector;
3. $(\dot u_{\text{in}}, \dot v_{\text{in}}, \dot w_{\text{in}})$ direction of the proton at the entrance detector (as unit vector);
4. $(\dot u_{\text{out}}, \dot v_{\text{out}}, \dot w_{\text{out}})$ direction of the proton at the exit detector (as unit vector);
5. $(e_{\text{in}},e_{\text{out}},t)$ where
    - $e_{\text{in}}$ and $e_{\text{out}}$ are the proton energy at the entrance and exit detectors, respectively. **Important note: if $e_{\text{in}}=0$, then $e_{\text{out}}$ will directly be interpreted as the water equivalent path length (WEPL) for the proton pair.**
    - $t$ is not used in PCT but can be used to store some useful scalar, e.g., the time-of-flight or the [GATE](https://github.com/OpenGATE/opengate) `TrackID` in `pctpairprotons`.
6. Optionally, mainly for the work in http://doi.org/10.1088/0031-9155/61/9/3258 $(\text{creatorProcess}, \text{nuclearProcess}, \text{order})$ are used in simulations to indicate
    - the process which created the exit particle,
    - whether the particle did encounter a nuclear interaction,
    - the number of particle interactions.

Note that PCT assumes that the proton beam goes along the $w$ axis in the positive direction ($w_{\text{in}} < w_{\text{out}}$).
