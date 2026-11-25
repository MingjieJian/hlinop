C***********************************************************************
C  Hydrogen bound opacity.  Designed to give the total neutral H opacity 
C  (bound-bound and bound free) at a given wavelength for a given list of 
C  lines at specified plasma conditions.  Occupation probability formalism
C  is used. Appropriate for applications like model stellar atmosphere 
C  calculations.
C
C  Needs the code for normalised line profiles HLINOP which can 
C  be used separately.
C
C  Paul Barklem, Nik Piskunov, Uppsala August 2003. 
C  Minor updates December 2024
C
C  A number of parts from or based on code by Deane Peterson and Bob 
C  Kurucz. 
C
C  Some contributions and significant input from Kjell Eriksson.
C
C  We would very much appreciate bug reports to: 
C  paul.barklem@physics.uu.se 
C
C  Table of Contents:
C  ------------------
C 
C  REAL HBOP(WAVE,N,NLO,NUP,WAVEH,NH,NHE,NE,T,DOP,NPOP,NL,TOTAL,CONTIN)   
C  REAL FUNCTION WCALC(NH, NE, NHW, NS, TEMP)	
C  SUBROUTINE HBOP_SET_EXCLUDED_VOLUME(SCALE)
C  SUBROUTINE HBOP_SET_LINE_PROFILE(MODE)
C  BLOCK DATA HBOPEVD
C  REAL FUNCTION HBF(WAVE,NLO)
C  REAL*8 FUNCTION HBNSTAR(WAVEV,NLO)
C  SUBROUTINE LOGINT(X,Y,N,XI,YI)
C  SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
C
C  NEEDS HLINPROF (hlinprof.f, with its table files) for the line 
C  profiles by default; HBOP_SET_LINE_PROFILE(0) uses HLINOP only.
C
C  NEEDS HLINOP which contains
C  REAL FUNCTION HLINOP(WAVE,NBLO,NBUP,WAVEH,T,XNE,H1FRC,HE1FRC,DOPPLE)
C  REAL FUNCTION HLCUT(WAVE,NBLO,NBUP,WAVEH)
C  REAL FUNCTION HPROFL(N,M,WAVE,WAVEH,T,XNE,H1FRC,HE1FRC,DOPPH)
C  REAL FUNCTION STARK1(N,M,WAVE,WAVEH,T,XNE)
C  REAL FUNCTION HFNM(N,M)
C  REAL FUNCTION VCSE1F(X)
C  REAL FUNCTION SOFBET(B,P,N,M)
C  REAL*8 FUNCTION AIRVAC(W)
C  REAL*8 FUNCTION VACAIR(W)  
C
C  Variables to be set / aware of before use:
C  ------------------------------------------
C
C  HBOP  : NLEVELS, NBF
C  WCALC : neutral-perturber excluded volume, set at run time with
C          CALL HBOP_SET_EXCLUDED_VOLUME(SCALE); default 0 (off)
C  HBOP  : line profiles, set at run time with
C          CALL HBOP_SET_LINE_PROFILE(MODE); default 1 (HLINPROF)
C  HLINOP: CUTNEXT, REFLECT
C  HPROFL: RAYLCUT
C
C***********************************************************************

      SUBROUTINE HBOP(WAVE, N, NLO, NUP, WAVEH, NH, NHE, NE, T, DOP,
     *                  NPOP, NL, TOTAL, CONTIN)
C
C  Returns the total absorption coefficient due to H bound levels at WAVE, 
C  for a given line list employing the occupation probability formalism 
C  of Daeppen, Anderson & Milhalas ApJ 1987 (DAM).
C
C  The extension of the occupation probability formalism to non-LTE by 
C  Hubeny, Hummer and Lanz A&A 1994 (HHL) is implemented (though see 
C  comments below).
C
C  To use the code in LTE:
C         set NL = 0 
C         NPOP is then not used, population computed internally
C  To use the code in non-LTE:
C         set NL = number of states to treat in non-LTE
C         set NPOP = array of number densities for each state
C         states above NL will be treated in LTE if needed
C                            
C  Note that Hubeny, Hummer and Lanz 1994 have extended the occupation 
C  probability formalism to non-LTE, though in a somewhat unclear manner 
C  in my opinion.  The main point to be aware of is that w(HHL) does not 
C  equal w(DAM).  The first is a conditional probability and the second 
C  is unconditional.  The problem is that the definition in sect 2.1 of HHL 
C  and derivation for w in HHL Appendix A is for the unconditional 
C  probability w(DAM).  In the rest of the paper it's conditional. The paper 
C  doesn't really make this difference clear IMO. 
C
C  Thus w(j)(HHL) = P(j|i) = P(i & j)/P(i) = w(j)/w(i) for j>i
C                                          = 1         for j<i
C  where w(m) is as computed in App A which is as defined in DAM
C  and sect 2.1 of HHL.  Since HHL assume w(i)=1 this leads to no problems.
C  
C 
C  wave   = interrogation wavelength in A               real*8
C  n      = number of lines                             integer
C  nlo    = array of lower principal quantum numbers    int array
C  nup    = array of upper principal quantum numbers    int array
C  waveh  = array of central wavelengths                real*8 array 
C  t      = kinetic temperature in K
C  ne     = electron number density in cm-3
C  nh     = number density of H I in cm-3
C  nhe    = number density of He I in cm-3
C  dop    = reduced Doppler width delta_lambda / lambda_0
C         = reduced Doppler width delta_nu / nu_0
C  npop   = number density of each level in cm-3
C  nl     = number of levels for which populations are pre-specified 
C           -- higher levels are assumed in LTE
C  total  = returns the total (line + continuous) absorption coefficient
C  contin = returns the continuous absorption coefficient
C
C  Important parameters
C  
C  nlevels = number of levels summed in partition functions and computed
C            populations and occupation probabilities. Must be larger or 
C            equal to the highest level in the line list and NL.
C  nbf     = highest level for which the bound-free component is 
C            accounted for (limited to 15 at present). 
C
C  The neutral-perturber (excluded-volume) part of the occupation 
C  probability is off by default; see HBOP_SET_EXCLUDED_VOLUME.
C
C  Paul Barklem, Uppsala, August 2003
C  Dec 2024, error catch for NUP > NLEVELS added
C  Nov 2026, excluded-volume switch (HBOP_SET_EXCLUDED_VOLUME); the 
C            populations are recomputed when the setting changes
C  Nov 2026, line broadening by H with the total H I density NH
C            (was the ground-state population NP(1))
C  Nov 2026, dissolved fraction from the code's own level energies
C            in vacuum wavelengths (HBNSTAR)
C  Nov 2026, levels whose red-wing taper has underflowed are skipped
C            (Inf/NaN beyond about 5 micron)
C  Nov 2026, line profiles from HLINPROF (Stehle & Hutcheon tables for
C            Lyman beta and up, Balmer and Paschen to n_up = 30; HLINOP
C            otherwise) with the HLINOP cutoffs; HBOP_SET_LINE_PROFILE
C
      IMPLICIT NONE
      INTEGER N, NLO(*), NUP(*), NLEVELS, I, J, NBF, NL, FF,FN
      PARAMETER (NLEVELS = 100)  
      PARAMETER (NBF = 6)  
      REAL NH, NHE, NE, T, CONTIN, DOP, TOTAL
      REAL NHL, NHEL, NEL, TL
      REAL*8 W(NLEVELS), G(NLEVELS), WGE(NLEVELS)
      REAL NLTE(NLEVELS), NPOP(*), NP(NLEVELS)
      REAL Z, H, C, HC, K, KT, LINE, SIGMA, CHI, SF
      REAL IONH, X, HFNM, FNM(NLEVELS,NLEVELS), HLINOP, HBF
      REAL TS, TF
      REAL*8 WAVE, WAVEH(*), REDCUT, WAVEV, HBNSTAR, AIRVAC
      REAL*8 EHYD(NLEVELS), CONTH(NLEVELS), WCALC, D, WSTAR, TAPER
      LOGICAL FIRST
      REAL EXCLUDED_VOLUME_SCALE, EVSL
      COMMON /HBOPEV/ EXCLUDED_VOLUME_SCALE
      INTEGER LINE_PROFILE_MODE
      COMMON /HBOPLP/ LINE_PROFILE_MODE
      EXTERNAL HBOPEVD
      REAL PROF, HLCUT, WCUTF
      REAL*8 CAA
      PARAMETER (CAA=2.99792458D18)   ! speed of light in A/s
      SAVE 
      DATA FIRST/.TRUE./, EVSL/-1./
      PARAMETER (H=6.62607E-27, C=2.9979E10, K=1.38065E-16)
      PARAMETER (IONH=2.17991e-11) 
      PARAMETER (HC=H*C)
      PARAMETER (SF=0.0265) ! pi*e*e/m/c [cgs]  sigma = sf * f * norm-prof
C
      IF (FIRST) THEN
C
C  Compute energy levels in the H atom, and store f-values
C
        EHYD(1) =      0.000D0
        EHYD(2) =  82259.105D0
        EHYD(3) =  97492.302D0
        EHYD(4) = 102823.893D0
        EHYD(5) = 105291.651D0
        EHYD(6) = 106632.160D0
        EHYD(7) = 107440.444D0
        EHYD(8) = 107965.051D0
        DO 1 I = 9, NLEVELS
 1      EHYD(I) = 109678.764D0 - 109677.576D0/(I*I)
        DO 2 I = 1, NLEVELS
 2      CONTH(I) = 109678.764D0 - EHYD(I)
        FIRST = .FALSE.
        DO 4 I = 1, NLEVELS-1
        DO 3 J = I+1, NLEVELS
 3      FNM(I,J) = HFNM(I,J)
 4      CONTINUE
      ENDIF 
C
C  Compute partition function and populations in LTE.
C  Need LTE populations even in non-LTE case
C  for calculating stimulated free-bound emission
C
C  Create population array NP for all levels 
C
      IF ((ABS(NH-NHL).LT.1E-3*NH).AND.(ABS(NHE-NHEL).LT.1E-3*NHE).AND.
     +   (ABS(NE-NEL).LT.1E-3*NE).AND.(ABS(T-TL).LT.1E-3*T).AND.
     +   (EXCLUDED_VOLUME_SCALE.EQ.EVSL))
     *   GOTO 16
      Z = 0.
      KT = K*T
      DO 10 I = 1, NLEVELS
      W(I) = WCALC(NH, NE, NHE, FLOAT(I), T)
      G(I) = 2.*I*I
      WGE(I) = W(I)*G(I)*DEXP(-HC*EHYD(I)/KT)
 10   Z = Z + WGE(I)
      DO 15 I = 1, NLEVELS
 15   NLTE(I) = NH*WGE(I)/Z  ! LTE populations
      NHL = NH
      NHEL = NHE
      NEL = NE
      TL = T
      EVSL = EXCLUDED_VOLUME_SCALE
 16   CONTINUE
C
      DO 20 I = 1, NLEVELS
      IF (I .LE. NL) THEN 
        NP(I) = NPOP(I)        ! non-LTE populations for states below NL
      ELSE 
        NP(I) = NLTE(I)        ! LTE populations otherwise
      ENDIF
 20   CONTINUE
C
C  Compute line opacity components
C
      LINE = 0.
      DO 30 I = 1, N
      IF (NUP(I).GT.NLEVELS) STOP 'HBOP : NUP > NLEVELS'           
      IF (LINE_PROFILE_MODE.EQ.1) THEN
C  HLINPROF (per A) converted to per Hz; HLINOP cutoffs applied to the
C  tabulated profiles (HLINOP, used for the other lines, has them)
        WCUTF = HLCUT(WAVE,NLO(I),NUP(I),WAVEH(I))
        PROF = 0.
        IF (WCUTF.GT.0.) THEN
          CALL HLINPROF(WAVE,WAVEH(I),T,NE,NLO(I),NUP(I),NH,NHE,
     *                  DOP,PROF)
          PROF = PROF * WAVE * WAVE / CAA
          IF ((NLO(I).LE.3).AND.(NUP(I).LE.30).AND.
     *        .NOT.((NLO(I).EQ.1).AND.(NUP(I).EQ.2))) PROF = PROF*WCUTF
        ENDIF
      ELSE
        PROF = HLINOP(WAVE,NLO(I),NUP(I),WAVEH(I),T,NE,NH,NHE,DOP)
      ENDIF
      SIGMA = SF * FNM(NLO(I),NUP(I)) * PROF
      CHI = 0.
      IF (W(NLO(I)).GT.0.)  ! populations would be zero also
     +   CHI = SIGMA * ( W(NUP(I))/W(NLO(I))*NP(NLO(I)) - 
     -                 NP(NUP(I))*G(NLO(I))/G(NUP(I)) )  !stim em term
 30   LINE = LINE + CHI
C
C  Compute the continuous components - sum over NBF lowest states
C     
      CONTIN = 0.
C  The continuum is computed in vacuum wavelengths; WAVE is in air 
C  above 2000 A (as for the lines, see HLINOP)
      WAVEV = WAVE
      IF (WAVE.GT.2000.D0) WAVEV = AIRVAC(WAVE)
      DO 40 I = 1, NBF
C
C  Compute dissolved fraction for the wavelength for this lower level
C
C  D does not decline fast enough to overcome the bound-free extrapolation
C  and this can lead to large dissolved contibutions very far from 
C  the bound-free jump, particularly from the Lyman continuum in cool stars.
C  This is certainly unphysical as the extrapolation will be invalid 
C  for such cases (which corresponds basically to H2+ if the perturber 
C  is a proton).  I apply an arbitrary tapering redward of the alpha line.
C 
CC      CHI = 0.     
CC      IF (WAVE.GT.(1.d8/CONTH(I))) THEN 
      D = 1.
      TAPER = 1.
C  Effective quantum number n* of the dissolved upper state reached 
C  from level I at this (vacuum) wavelength, from the same level 
C  energies as the lines (Nov 2026; before, IONH = R_infinity was used,
C  which put the n = 2 edge 2 A blueward of the series limit, and the 
C  air wavelength was used as if in vacuum)
      X = HBNSTAR(WAVEV, I)
      IF(X .GT. 0.) THEN
        WSTAR = WCALC(NH, NE, NHE, X, T)
        D = 1.0-(WSTAR/W(I))
        REDCUT = 1.D8/(EHYD(I+1)-EHYD(I))
        IF (WAVEV .GT. REDCUT) TAPER = DEXP((REDCUT-WAVEV)*1.D-2)
      ENDIF	
C  Far redward of the alpha line the taper underflows to 0 while the 
C  extrapolated cross section HBF may overflow (0 * Inf = NaN beyond 
C  about 5 micron before Nov 2026): skip such levels
      CHI = 0.
      IF (TAPER.GT.1.D-30) CHI = (NP(I) - NLTE(I)*DEXP(-HC/(WAVEV*1E-8)
     *       /KT)) * D * HBF(WAVEV, I) * TAPER
CC      ENDIF
 40   CONTIN = CONTIN + CHI
C
C  Sum up and we're done
C
      TOTAL = CONTIN + LINE
C
      RETURN
      END

    
C***********************************************************************
 
      REAL*8 FUNCTION WCALC(NH, NE, NHE, NS, TEMP)	
C  	
C  Computes the occupation probability for a H atom in 
C  state with effective principal quantum number NS in a plasma
C  enviroment with NH, NE, NHE number densities of H, ions,
C  and He atoms respectively.  This code assumes the perturbing
C  neutral H and He are in the ground state, (noting the hard
C  sphere approximation used is quite crude anyway) and that ions 
C  are predominantly singly ionized (ie. N(Z=1 ions) = Ne).
C
C  See eqn 4.71 Hummer & Milhalas (1988) 
C  Ions are now treated via Appendix A Hubeny, Hummer & Lanz (1994)
C  which is a fit to detailed calculation including correlations,
C  thus the temperature dependence
C
C  Sizes of neutral atoms adopted are sqrt(<r^2>)
C
C  The excluded volume of the neutral perturbers is multiplied by
C  EXCLUDED_VOLUME_SCALE (set with HBOP_SET_EXCLUDED_VOLUME):
C  0 = neutral term off (default), 1 = Hummer & Mihalas eq 4.71.
C  
C  Coded by Paul Barklem and Kjell Eriksson, Aug 2003
C 	
      IMPLICIT NONE
      REAL NH, NE, NHE, NS, TEMP
      REAL*8 IONH, A0, E, PI, K, CHI, RIH, NEUTR, ION, NS2, NS4 
      REAL*8 F, A, BETAC, X, WNEUTR, WION, X1, X2      
      REAL EXCLUDED_VOLUME_SCALE
      COMMON /HBOPEV/ EXCLUDED_VOLUME_SCALE
      PARAMETER (IONH=2.17991E-11) 
      PARAMETER (A0=5.29177E-9)
      PARAMETER (E=4.803207E-10) 
      PARAMETER (PI=3.1415927)
C
      NS2 = NS*NS
      NS4 = NS2*NS2
C  
C  Neutral perturbers
C
      CHI=IONH/(NS2)
      RIH=SQRT(2.5*NS4 + 0.5*NS2)*A0
      X1=RIH + 1.73*A0
      X2=RIH + 1.02*A0
      NEUTR=NH*X1*X1*X1 + NHE*X2*X2*X2
      WNEUTR = DEXP(-4.18879 * EXCLUDED_VOLUME_SCALE * NEUTR)  ! 4*pi/3
C
C  Charged perturbers
C      
      K=1.
      IF(NS .gt. 3.) THEN
C       K=5.33333333*(NS/(NS+1.))**2 * (NS + 7./6.)/(NS2+NS+0.5)
        K=5.33333333*NS/(NS+1.)/(NS+1.)
      ENDIF
c      ION  = NE*16.*(E*E/CHI/DSQRT(K))**3
c      WION = DEXP(-4.*PI/3. * (ION))
      IF ((NE.GT.10.).AND.(TEMP.GT.10.)) THEN  ! just in case!!! 
        A = 0.09 * DEXP(0.16667d0*LOG(NE)) / SQRT(TEMP)
        X = DEXP(3.15*DLOG(1.+A))
        BETAC = 8.3E14 * EXP(-0.66667*LOG(NE)) * K / NS4
        F = 0.1402d0*X*BETAC*BETAC*BETAC /
     /             (1.+0.1285*X*BETAC*DSQRT(BETAC))
        WION = F/(1.d0+F)
      ELSE 
        WION = 1.0d0     
      ENDIF 
C      
      WCALC = WION * WNEUTR 
      RETURN
      END


C***********************************************************************

      SUBROUTINE HBOP_SET_EXCLUDED_VOLUME(SCALE)
C
C  Sets the neutral-perturber part of the occupation probability used 
C  by WCALC and hence by HBOP (lines, level populations and the 
C  dissolved-level pseudo-continuum).
C
C  SCALE multiplies the excluded volume of neutral H and He perturbers
C  (hard spheres, Hummer & Mihalas 1988 eq 4.71):
C     SCALE = 0.  neutral term off (default if this is never called)
C     SCALE = 1.  original hard-sphere formula of Hummer & Mihalas
C  Values in between scale the volume (the exponent of w_neutral).
C  The charged-particle part is not affected.
C
C  Call it before HBOP; the setting stays until it is called again.
C
      REAL SCALE
      REAL EXCLUDED_VOLUME_SCALE
      COMMON /HBOPEV/ EXCLUDED_VOLUME_SCALE
      EXTERNAL HBOPEVD
      IF (SCALE.LT.0.) STOP 'HBOP_SET_EXCLUDED_VOLUME: SCALE < 0'
      EXCLUDED_VOLUME_SCALE = SCALE
      RETURN
      END

C***********************************************************************

      SUBROUTINE HBOP_SET_LINE_PROFILE(MODE)
C
C  Selects the line profiles used by HBOP:
C     MODE = 1  HLINPROF (default): Stehle & Hutcheon (1999) tables for
C               Lyman (n_up 3-30), Balmer and Paschen (n_up to 30), 
C               HLINOP for the other lines; needs hlinprof.f and its 
C               table files.  The HLINOP cutoffs are applied.
C     MODE = 0  HLINOP for all lines (the behaviour before Nov 2026)
C
      INTEGER MODE
      INTEGER LINE_PROFILE_MODE
      COMMON /HBOPLP/ LINE_PROFILE_MODE
      EXTERNAL HBOPEVD
      IF ((MODE.NE.0).AND.(MODE.NE.1)) 
     *   STOP 'HBOP_SET_LINE_PROFILE: MODE must be 0 or 1'
      LINE_PROFILE_MODE = MODE
      RETURN
      END

C***********************************************************************

      BLOCK DATA HBOPEVD
C
C  Defaults: neutral-perturber excluded volume off; line profiles from
C  HLINPROF.
C
      REAL EXCLUDED_VOLUME_SCALE
      COMMON /HBOPEV/ EXCLUDED_VOLUME_SCALE
      INTEGER LINE_PROFILE_MODE
      COMMON /HBOPLP/ LINE_PROFILE_MODE
      DATA EXCLUDED_VOLUME_SCALE /0./
      DATA LINE_PROFILE_MODE /1/
      END
     

C***********************************************************************

      REAL FUNCTION HBF(WAVE, NLO)
C
C  Returns the bound-free cross section due to H in level NLO
C  at wavelength WAVE (in Angstroms).  
C
C  It gives an extrapolation redward of usual bound free jumps
C  so that we can use the occupation probability formalism
C  for "broadening" the jump due to nearby particles
C
C  Data (Gaunt factors - Karzas & Latter) from MARCS code
C
      INTEGER N(15), NLO
      REAL L1(12), G1(12), L2(12), G2(12), L3(14), G3(14)
      REAL L4(14), G4(14), L5(14), G5(14), L6(14), G6(14)
      REAL L7(15), G7(15), L8(15), G8(15), L9(6), G9(6)
      REAL L10(16), G10(16), L11(6), G11(6), L12(6), G12(6)
      REAL L13(6), G13(6), L14(6), G14(6), L15(6), G15(6)
      REAL G, W
      REAL*8 WAVE
C
      DATA N/ 12, 12, 14, 14, 14, 14, 15, 15, 6, 16, 6, 6, 6, 6, 6/
C
      DATA L1/   182.,   241.,   300.,   356.,   408.,   456.,
     *           538.,   631.,   729.,   821.,   877.,   912./
      DATA G1/  .9939,  .9948,  .9856,  .9721,  .9571,  .9423,
     *          .9157,  .8850,  .8531,  .8246,  .8076,  .7973/
      DATA L2/   215.,   301.,   398.,   503.,   614.,   729.,
     *           965.,  1313.,  1824.,  2525.,  3144.,  3647./
      DATA G2/  1.043,  1.058,  1.063,  1.061,  1.056,  1.049,
     *          1.033,  1.009,  .9755,  .9342,  .9011,  .8761/
      DATA L3/   222.,   316.,   424.,   545.,   677.,   821.,
     *          1132.,  1641.,  2525.,  4103.,  6034.,  6933.,
     *          7529.,  8207./
      DATA G3/  1.053,  1.072,  1.081,  1.083,  1.081,  1.078,
     *          1.068,  1.051,  1.024,  .9854,  .9458,  .9293,
     *          .9189,  .9075/
      DATA L4/   224.,   321.,   434.,   561.,   703.,   859.,
     *          1205.,  1799.,  2918.,  5252.,  8895., 10998.,
     *         12576., 14588./
      DATA G4/  1.056,  1.077,  1.087,  1.091,  1.091,  1.089,
     *          1.082,  1.069,  1.047,  1.014,  .9736,  .9542,
     *          .9408,  .9247/
      DATA L5/   226.,   324.,   438.,   569.,   715.,   877.,
     *          1241.,  1882.,  3144.,  6034., 11397., 15095.,
     *         18235., 22794./
      DATA G5/  1.058,  1.080,  1.090,  1.095,  1.095,  1.094,
     *          1.089,  1.078,  1.060,  1.030,  .9925,  .9719,
     *          .9563,  .9358/
      DATA L6/   226.,   325.,   441.,   573.,   722.,   887.,
     *          1262.,  1931.,  3282.,  6564., 13448., 18916.,
     *         24120., 32915./
      DATA G6/  1.059,  1.081,  1.092,  1.097,  1.098,  1.097,
     *          1.092,  1.083,  1.067,  1.041,  1.006,  .9851,
     *          .9682,  .9436/
      DATA L7/   227.,   326.,   442.,   576.,   726.,   894.,
     *          1275.,  1961.,  3372.,  6933., 15095., 22347.,
     *         29992., 36616., 44694./
      DATA G7/  1.060,  1.082,  1.093,  1.098,  1.099,  1.099,
     *          1.095,  1.086,  1.071,  1.047,  1.015,  .9952,
     *          .9776,  .9641,  .9494/
      DATA L8/   227.,   326.,   443.,   578.,   729.,   897.,
     *          1284.,  1982.,  3433.,  7196., 16398., 25326.,
     *         35615., 45361., 58446./
      DATA G8/  1.060,  1.082,  1.094,  1.099,  1.100,  1.100,
     *          1.096,  1.088,  1.074,  1.052,  1.022,  1.003,
     *          .9853,  .9708,  .9540/
      DATA L9/   227.,   901.,  7383., 40703., 67537., 74126./
      DATA G9/  1.060,  1.101,  1.056,  .9917,  .9632,  .9576/
      DATA L10/  227.,   327.,   445.,   580.,   732.,   903., 
     *          1294.,  2006.,  3507.,  7529., 18235., 29992.,
     *         45588., 63316., 72940., 91175./
      DATA G10/ 1.060,  1.083,  1.095,  1.100,  1.102,  1.101,
     *          1.098,  1.091,  1.078,  1.058,  1.032,  1.014,
     *          .9970,  .9814,  .9737,  .9606/
      DATA L11/   82.,   905.,  7636., 49822., 96995.,111189./
      DATA G11/ .9379,  1.102,  1.060,  1.002,  .9703,  .9632/
      DATA L12/   82.,   905.,  7720., 53950.,112562.,132138./
      DATA G12/ .9379,  1.102,  1.062,  1.005,  .9732,  .9653/ 
      DATA L13/   82.,   906.,  7793., 57343.,130250.,154534./
      DATA G13/ .9380,  1.102,  1.063,  1.008,  .9758,  .9672/
      DATA L14/   82.,   907.,  7846., 60381.,147056.,178775./
      DATA G14/ .9380,  1.103,  1.064,  1.011,  .9780,  .9689/
      DATA L15/   82.,   908.,  7887., 63316.,162813.,207216./
      DATA G15/ .9380,  1.103,  1.064,  1.014,  .9803,  .9703/
C
      W = WAVE
      G = 0.
      IF (NLO.EQ.1) CALL LOGINT(L1,G1,N(1),W,G)
      IF (NLO.EQ.2) CALL LOGINT(L2,G2,N(2),W,G)
      IF (NLO.EQ.3) CALL LOGINT(L3,G3,N(3),W,G)
      IF (NLO.EQ.4) CALL LOGINT(L4,G4,N(4),W,G)
      IF (NLO.EQ.5) CALL LOGINT(L5,G5,N(5),W,G)
      IF (NLO.EQ.6) CALL LOGINT(L6,G6,N(6),W,G)
      IF (NLO.EQ.7) CALL LOGINT(L7,G7,N(7),W,G)
      IF (NLO.EQ.8) CALL LOGINT(L8,G8,N(8),W,G) 
      IF (NLO.EQ.9) CALL LOGINT(L9,G9,N(9),W,G)
      IF (NLO.EQ.10) CALL LOGINT(L10,G10,N(10),W,G)
      IF (NLO.EQ.11) CALL LOGINT(L11,G11,N(11),W,G)
      IF (NLO.EQ.12) CALL LOGINT(L12,G12,N(12),W,G)
      IF (NLO.EQ.13) CALL LOGINT(L13,G13,N(13),W,G)
      IF (NLO.EQ.14) CALL LOGINT(L14,G14,N(14),W,G)
      IF (NLO.EQ.15) CALL LOGINT(L15,G15,N(15),W,G)
C
      HBF = 1.044E-26 * G * WAVE**3.d0 / (NLO**5.d0)
      RETURN
      END


C***********************************************************************
      REAL*8 FUNCTION HBNSTAR(WAVEV, NLO)
C
C  Effective principal quantum number n* of the upper state reached 
C  from level NLO by a photon of vacuum wavelength WAVEV (A), with the 
C  level energies used for the lines (HLINOP, HBOP):
C     1/n*^2 = (CONTH(NLO) - 1e8/WAVEV) / R_H,
C  R_H = 109677.576 cm^-1.  Returns 0 blueward of the series limit 
C  (photoionisation, no dissolution).
C
      IMPLICIT NONE
      INTEGER NLO, I
      REAL*8 WAVEV, Q, EHYD(100), CONTH(100)
      LOGICAL FIRST
      SAVE FIRST, CONTH
      DATA FIRST/.TRUE./
      IF (FIRST) THEN
        EHYD(1) =      0.000D0
        EHYD(2) =  82259.105D0
        EHYD(3) =  97492.302D0
        EHYD(4) = 102823.893D0
        EHYD(5) = 105291.651D0
        EHYD(6) = 106632.160D0
        EHYD(7) = 107440.444D0
        EHYD(8) = 107965.051D0
        DO 1 I = 9, 100
 1      EHYD(I) = 109678.764D0 - 109677.576D0/(I*I)
        DO 2 I = 1, 100
 2      CONTH(I) = 109678.764D0 - EHYD(I)
        FIRST = .FALSE.
      ENDIF
      Q = CONTH(NLO) - 1.D8/WAVEV
      IF (Q.LE.0.D0) THEN
        HBNSTAR = 0.D0
      ELSE
        HBNSTAR = SQRT(109677.576D0/Q)
      ENDIF
      RETURN
      END

C***********************************************************************

      SUBROUTINE LOGINT(X, Y, N, XI, YI)
C
C  Does log interpolation and extrapolation of array Y
C  of length N.  Returns value YI at XI.
C  Assumes monotonically increasing X.
C  Bug fix 20230711 thanks to Richard Hoppe
C  Changed to log-linear 20230712
C
      IMPLICIT NONE
      INTEGER N, I, J
      REAL X(*), Y(*), XI, YI
      REAL XL(N), YL(N), XIL, YIL, XS(3), YLS(3), DY        
C
      J = 0
      DO 10 I = 1, N
      YL(I) = LOG10(Y(I))
      IF (X(I).LT.XI) J = I
 10   CONTINUE
      IF (J.EQ.0) J = 1
      IF (J.GT.N-3) J = N-2
C
C  polynomial 3 point interpolation or extrapolation
C     
      XS(1) = X(J)
      XS(2) = X(J+1)
      XS(3) = X(J+2)
      YLS(1) = YL(J)
      YLS(2) = YL(J+1)
      YLS(3) = YL(J+2)
      CALL POLINT(XS, YLS, 3, XI, YIL, DY)
      YI = 10.**YIL
C
      RETURN
      END
   

C***********************************************************************
C  Numerical recipes polynomial interpolation routine

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      INTEGER I,N,NMAX,M,NS
      REAL XA(*),YA(*),X,Y,DY,DIF,DIFT,HO,HP,W,DEN
      PARAMETER (NMAX=10) 
      DIMENSION C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.d0)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END


C***********************************************************************
