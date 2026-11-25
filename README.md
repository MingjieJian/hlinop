# HLINOP

- New version with bug fix thanks to Richard Hoppe, 11 July 2023; see history.txt
- hbop_auto.f added 17 December 2025; small changes to hbop.f
- November 2026: neutral-perturber term switchable (off by default); HLINPROF uses the full Stehle & Hutcheon tables and HBOP uses HLINPROF; HLINOP broadening mechanisms convolved (unit area); test suite; see history.txt

A set of codes for computing hydrogen line profiles and opacities developed by Paul Barklem and Nik Piskunov.  This code was first presented at IAU 210, Modelling Stellar Atmospheres, in Uppsala June 2002 (Barklem & Piskunov 2003). There a need for codes to compute these opacities for atellar atmosphere applications. Codes such as HLINOP and BALMER9 have been provided by Kurucz and Peterson (http://kurucz.harvard.edu/) in the 1970's based on the broadening theories of the time.  Since then, understanding of H line broadening has improved and there is a need to update these codes. Three codes have been produced, namely:

HLINPROF: For detailed, accurate calculation of hydrogen line profiles from the Stark-broadening tables of Stehle & Hutcheon (1999), suitable for detailed analysis of Balmer and other lines. The table `stehle_long.dat` covers Lyman (n_up 3-30; Lyman alpha is left to HLINOP), Balmer (3-30) and Paschen (4-30) lines, temperatures 2500 K - 1.26e6 K and electron densities 1e10 - 3.2e19 cm^-3 (fewer densities for higher n_up, up to each line's Inglis-Teller limit). It is built from the CDS catalogue J/A+AS/140/93 with `tools/make_stehle_table.py` (little-endian; see the script for the format). Above a line's last tabulated density the profile of that density is scaled quasistatically (width proportional to n_e^(2/3)); below 1e10 cm^-3 the first profile is used; other lines and temperatures outside the table go to HLINOP. The tabulated profiles include the thermal Doppler broadening of H; since November 2026 HLINPROF convolves them with the remaining Doppler broadening (the Doppler width passed in beyond the thermal one, e.g. microturbulence) and with the Lorentz profile of the radiative, self and helium widths, which were previously added beyond 0.2 A and so increased the area by up to 4 per cent (`HLINOP_SET_CONVOLUTION(0)` restores the earlier behaviour). The table file names are compiled in from an include file `DATA.FILES` (see `SYNTH/Makefile` for how to generate it). The older six-line tables (`stehle.dat.*`) and `hlinprof_nofiles.f` are not used by this version.

HLINOP: For approximate, quick calculation of any line of neutral H, suitable for model atmosphere calculations.  Based on the original Kurucz and Peterson version. Since November 2026 the broadening mechanisms are convolved instead of added: the Stark profile convolves the ion quasistatic profile with the electron profile (impact core and quasistatic electron wing), and the line profile convolves the Doppler profile (with fine structure) with the Lorentz (resonance, van der Waals, natural) and Stark profiles. The earlier code had profile areas of 1.7-2.1; against the Stark-only profiles of Stehle & Hutcheon (1999, Balmer and Paschen) and Lemke (1997, VCS, Brackett) the median half widths and peaks now agree to 2-5 per cent. `CALL HLINOP_SET_CONVOLUTION(0)` restores the earlier profiles; Lyman alpha keeps them for now.

HBOP: A wrapper around hlinprof.f and hlinop.f to implement the occupation probability formalism of Daeppen, Anderson and Milhalas (1987) and thus account for the merging of lines into the continuum. The neutral-perturber (excluded-volume) part of the occupation probability is off by default; `CALL HBOP_SET_EXCLUDED_VOLUME(1.0)` before calling HBOP restores the original hard-sphere treatment of Hummer & Mihalas (1988), and values between 0 and 1 scale the excluded volume. Line profiles come from HLINPROF by default (with the HLINOP cutoffs at the series limits), so HBOP needs hlinprof.f and its table; `CALL HBOP_SET_LINE_PROFILE(0)` uses HLINOP for all lines as before.  As bound-free and bound-bound opacities are merged, the distinction between line and continuous opacity becomes blurred, so one should be careful not to include things twice (including Rayleigh scattering off Lyman alpha). For *efficient* implementation you may have to adapt this, in particular so that populations are not recomputed.

HBOP_AUTO: This code is a wrapper to HBOP, which will give the total bound-x (bound-bound and bound-free) neutral H opacity at a given wavelength, from the occupation probability formalism, WITHOUT the need for a line list.  This wrapper will simply produce a line list including all lines up to a highest principal quantum number NMAX, in the relevant spectral region (i.e. at the interrogation wavelength).  

The main codes are in the root directory.  An example of their application to the SYNTH spectral synthesis program by Nik Piskunov is also provided.

Views, comments, bug reports and suggested improvements all welcome: paul.barklem@physics.uu.se

Please cite the DOI if you make use of the code in research:
[![DOI](https://zenodo.org/badge/21607/barklem/hlinop.svg)](https://zenodo.org/badge/latestdoi/21607/barklem/hlinop)

References:

- Barklem, P. S., and N. Piskunov. “Hydrogen Balmer Lines as Probes of Stellar Atmospheres,” Vol. 210, 2003. http://adsabs.harvard.edu/abs/2003IAUS..210P.E28B.

- Dappen, Werner, Lawrence Anderson, and Dimitri Mihalas. “Statistical Mechanics of Partially Ionized Stellar Plasma - The Planck-Larkin Partition Function, Polarization Shifts, and Simulations of Optical Spectra.” The Astrophysical Journal 319 (August 1, 1987): 195–206. doi:10.1086/165446.

- Stehle, C., and R. Hutcheon. “Extensive tabulations of Stark broadened hydrogen line profiles.” Astronomy and Astrophysics Supplement Series 140 (1999): 93. doi:10.1051/aas:1999118. CDS catalogue J/A+AS/140/93.

- Lemke, M. “Extended VCS Stark broadening tables for hydrogen - Lyman to Brackett series.” Astronomy and Astrophysics Supplement Series 122 (1997): 285. doi:10.1051/aas:1997134.
