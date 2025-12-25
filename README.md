# FastFourierBeamformers

Fourier beamformers for pulse-echo ultrasound based on a k-space representation of the image

The goal of the Fourier beamformer is to reconstruct the ultrasound image more quickly than a delay-and-sum beamformer by leveraging underlying Stolt mappings from the frequency-domain representation of the recorded channel data to the k-space (i.e., wavenumber space) representation of the image.  Early work applies Fourier beamforming to plane-wave imaging [1][2] and multistatic synthetic aperture (i.e., full-matrix capture—FMC) [3], each with its own set of Stolt mappings.  Recent work extended the Stolt’s mapping used for the multistatic setup [3] to a focused-transmit scenario [4] with the caveat that the focal depth must be the same for all transmit beams, and each transmit focus must be directly underneath a corresponding transducer element.  Many conventional focused-transmit sequences fall outside this limited use case.  Because different imaging sequences require different Stolt mappings [4], and in many cases, a closed-form Stolt mapping may not be known, Fourier beamforming is often restricted to a limited set of use cases.  

Fourier beamforming currently lacks a unifying framework that applies to all transmit sequences.  [Retrospective encoding for conventional ultrasound sequences (REFoCUS)](https://github.com/nbottenus/REFoCUS) [5] could be applied to any transmit sequence to recover multistatic data, which has a known Stolt mapping [3].  Alternatively, [Reverse-Time Migration (RTM)](https://github.com/rehmanali1994/FourierDomainBeamformer) generalizes to all transmit sequences without requiring an intermediate multistatic representation [6]; therefore, this repository also presents a generalized Stolt mapping based on a k-space representation of RTM.  Following the spirit of [Garcia's work on f-k migration [2]](https://github.com/rehmanali1994/Plane_Wave_Ultrasound_Stolt_F-K_Migration.github.io), my goal was to make both the previous Fourier beamformers and the RTM-based approach more broadly accessible.  Rather than assert that any particular Fourier beamforming method is optimal in all situations, this repository aims to provide an accessible implementation of each Fourier beamformer to clarify the connection between the mathematical theory and its practical implementation for purely educational purposes.  

# Citing this Work
If you use the code/algorithm for research, please cite the SPIE conference proceedings paper: 

> Ali, R. "Fast Fourier Beamforming for Arbitrary Ultrasound Imaging Sequences Based on a K-Space Implementation of Reverse-Time Migration". Medical Imaging 2026: Ultrasonic Imaging and Tomography. SPIE, 2026.

You can reference a static version of this code by its DOI number: [![DOI](https://zenodo.org/badge/1122896318.svg)](https://doi.org/10.5281/zenodo.18056898)

See the SPIE poster below ![](FastFourierBeamforming_Poster.pdf)

# Code and Sample Datasets
Functions for conventional delay-and-sum (DAS) beamformers are contained in [Functions_BeamformDAS](https://github.com/rehmanali1994/FastFourierBeamformers/tree/main/Functions_BeamformDAS), and functions for Fourier beamforming are contained in [Functions_BeamformFastFourier](https://github.com/rehmanali1994/FastFourierBeamformers/tree/main/Functions_BeamformFastFourier).  Example scripts for multistatic synthetic aperture ([Examples_FMC](https://github.com/rehmanali1994/FastFourierBeamformers/tree/main/Examples_FMC)), plane-wave ([Examples_PlaneWave](https://github.com/rehmanali1994/FastFourierBeamformers/tree/main/Examples_PlaneWave)), and focused-transmit ([Examples_FocTx](https://github.com/rehmanali1994/FastFourierBeamformers/tree/main/Examples_FocTx)) imaging are provided.  

**Please download the sample data (FieldII.mat) under the [releases](https://github.com/rehmanali1994/FastFourierBeamformers/releases) tab for this repository, and place that data in the ([Examples_FMC](https://github.com/rehmanali1994/FastFourierBeamformers/blob/main/Examples_FMC)) directory. Alternatively you can run [Examples_FMC/GenFullSynthDataPhantom.m](https://github.com/rehmanali1994/FastFourierBeamformers/blob/main/Examples_FMC/GenFullSynthDataPhantom.m) to recreate FieldII.mat—this will require the [Field II [7][8]](https://field-ii.dk/) simulation software.**

Every script with a name formatted `recon_*.m` demonstrates a specific reconstruction technique. The multistatic FMC dataset `FieldII.mat` must be converted to a plane-wave dataset (`FieldII_PW.mat` will be created by [Examples_PlaneWave/ConvertFMC2PW.m](https://github.com/rehmanali1994/FastFourierBeamformers/blob/main/Examples_PlaneWave/ConvertFMC2PW.m)) or a focused-transmit dataset (`FieldII_FocTx.mat` will be created by [Examples_FocTx/ConvertFMC2FocTx.m](https://github.com/rehmanali1994/FastFourierBeamformers/blob/main/Examples_FocTx/ConvertFMC2FocTx.m)) before the `recon_*.m` scripts in the [Examples_PlaneWave](https://github.com/rehmanali1994/FastFourierBeamformers/tree/main/Examples_PlaneWave) or [Examples_FocTx](https://github.com/rehmanali1994/FastFourierBeamformers/tree/main/Examples_FocTx) directories can be run.  Each conversion script applies transmit beamforming ([Functions_BeamformDAS/focus_fs_to_TxBeam.m](https://github.com/rehmanali1994/FastFourierBeamformers/blob/main/Functions_BeamformDAS/focus_fs_to_TxBeam.m)) to the multistatic FMC dataset to simulate the specific transmit sequence (plane wave or focused transmit).

# Sample Results

[Field II [7][8]](https://field-ii.dk/) simulations of radio-frequency channel data for multistatic synthetic aperture, plane-wave, and walking-aperture focused transmit sequences with the same acquisition parameter described in [6] were used to compare delay-and-sum (DAS) beamforming [9], the previous Fourier beamformers [1]-[4], a multistatic FMC Stolt migration [3] after REFoCUS [5], and the RTM-based Fourier beamformer proposed in this work.  For plane-wave imaging we found that [1] generally provided better resolution and contrast than [2] and became the Fourier beamformer used for comparison.  For focused transmit imaging, [4] became the basis for comparison.  

![](GitHubResults.png)

# References

[1] Cheng, J., & Lu, J. Y. (2006). Extended high-frame rate imaging method with limited-diffraction beams. IEEE transactions on ultrasonics, ferroelectrics, and frequency control, 53(5), 880-899.

[2] Garcia, D., Le Tarnec, L., Muth, S., Montagnon, E., Porée, J., & Cloutier, G. (2013). Stolt's fk migration for plane wave ultrasound imaging. IEEE transactions on ultrasonics, ferroelectrics, and frequency control, 60(9), 1853-1867.

[3] Hunter, A. J., Drinkwater, B. W., & Wilcox, P. D. (2008). The wavenumber algorithm for full-matrix imaging using an ultrasonic array. IEEE transactions on ultrasonics, ferroelectrics, and frequency control, 55(11), 2450-2462.

[4] Jiang, C., Liu, C., Zhan, Y., & Ta, D. (2022). The spectrum-beamformer for conventional B-mode ultrasound imaging system: principle, validation, and robustness. Ultrasonic Imaging, 44(2-3), 59-76.

[5] Ali, R., Herickhoff, C. D., Hyun, D., Dahl, J. J., & Bottenus, N. (2019). Extending retrospective encoding for robust recovery of the multistatic data set. IEEE transactions on ultrasonics, ferroelectrics, and frequency control, 67(5), 943-956.

[6] Ali, R. (2021). Fourier-based synthetic-aperture imaging for arbitrary transmissions by cross-correlation of transmitted and received wave-fields. Ultrasonic imaging, 43(5), 282-294.

[7] Jensen, J. A. (1996). Field: A program for simulating ultrasound systems. In 10TH NORDICBALTIC CONFERENCE ON BIOMEDICAL IMAGING, VOL. 4, SUPPLEMENT 1, PART 1: 351--353.

[8] Jensen, J. A., & Svendsen, N. B. (1992). Calculation of pressure fields from arbitrarily shaped, apodized, and excited ultrasound transducers. IEEE transactions on ultrasonics, ferroelectrics, and frequency control, 39(2), 262-267.

[9] Bae, M. H., & Jeong, M. K. (2000). A study of synthetic-aperture imaging with virtual source elements in B-mode ultrasound imaging systems. IEEE transactions on ultrasonics, ferroelectrics, and frequency control, 47(6), 1510-1519.
