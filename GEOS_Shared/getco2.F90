!-------------------------------------------------------------------------
!NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS!
!-------------------------------------------------------------------------
!BOP
!
! !INTERFACE:

 REAL FUNCTION GetCO2(Year, DayOfYear)

! !DESCRIPTION:
!
!  Given the year and day-of-year, this function returns the RCP60 CO2 concentration in 
!  mole fraction (volume mixing ratio).  If Year is less than 1765, the value for 1765
!  is returned.  If Year is greater than 2150, the value for 2150 is returned.  In the
!  original dataset, the value for 2150 is used for all years through 2500.  We choose
!  to truncate the list at 2151 for this application.
!
!  DayOfYear is expected to have a value of 1.00 at 0:00 UTC Jan 1.
!
!  In-line documentation from the source dataset is reproduced below:
! RCP6.0_MIDYEAR__CONCENTRATIONS____________________________
!
!CONTENT:           CMIP5 RECOMMENDATIONS FOR ANNUAL AVERAGE, GLOBAL MEAN CONCENTRATIONS.
!RUN:               RCP6, FINAL RELEASE, 30 May 2010
!RCP6.0 CONTACT:    AIM group, Toshihiko Masui (masui@nies.go.jp) & Kenichi Matsumoto (matsumoto.kenichi@nies.go.jp)
!DATE:              30/5/2010 01:02:33 (updated description) 
!MAGICC-VERSION:    6.3.09, 25 November 2009
!FILE PRODUCED BY:  RCP Concentration Calculation & Data Group, M. Meinshausen, S. Smith, K. Riahi, D. van Vuuren et al.
!DOCUMENTATION:     M. Meinshausen, S. Smith et al. "The RCP GHG concentrations and their extension from 1765 to 2500", in prep., Climatic Change.
!CMIP5 INFO:        http://cmip-pcmdi.llnl.gov/cmip5/
!RCP DATABASE:      http://www.iiasa.ac.at/web-apps/tnt/RcpDb
!FURTHER INFO:      For data sources, aknowledgements and further information, see http://www.pik-potsdam.de/~mmalte/rcps
!NOTE:              RCP6 starts 2005; 20th century data and earlier is provided for convenience; 
!
! !REVISION HISTORY:
!  29 Oct 2010  Nielsen, adapted for GEOS-5.
! 
!EOP
! ---------------------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Year
  INTEGER, INTENT(IN) :: DayOfYear

  REAL :: f,i,n
  INTEGER :: previous,current,next
  
  INTEGER, PARAMETER :: firstYear = 1764
  INTEGER, PARAMETER :: finalYear = 2151
  INTEGER, PARAMETER :: tableLength = finalYear-firstYear+1

  REAL, SAVE :: CO2ppmv(tableLength) = (/                                278.052, &
   278.052,  278.106,  278.220,  278.343,  278.471,  278.600,  278.733,  278.869, &
   279.009,  279.153,  279.302,  279.457,  279.618,  279.782,  279.943,  280.097, &
   280.243,  280.382,  280.518,  280.657,  280.803,  280.957,  281.118,  281.282, &
   281.443,  281.598,  281.747,  281.891,  282.031,  282.167,  282.299,  282.427, &
   282.551,  282.671,  282.787,  282.899,  283.007,  283.111,  283.211,  283.307, &
   283.400,  283.490,  283.578,  283.661,  283.735,  283.797,  283.847,  283.889, &
   283.926,  283.963,  284.001,  284.043,  284.086,  284.129,  284.167,  284.198, &
   284.223,  284.244,  284.263,  284.281,  284.300,  284.320,  284.340,  284.360, &
   284.380,  284.400,  284.385,  284.280,  284.125,  283.975,  283.825,  283.675, &
   283.525,  283.425,  283.400,  283.400,  283.425,  283.500,  283.600,  283.725, &
   283.900,  284.075,  284.225,  284.400,  284.575,  284.725,  284.875,  285.000, &
   285.125,  285.275,  285.425,  285.575,  285.725,  285.900,  286.075,  286.225, &
   286.375,  286.500,  286.625,  286.775,  286.900,  287.000,  287.100,  287.225, &
   287.375,  287.525,  287.700,  287.900,  288.125,  288.400,  288.700,  289.025, &
   289.400,  289.800,  290.225,  290.700,  291.200,  291.675,  292.125,  292.575, &
   292.975,  293.300,  293.575,  293.800,  294.000,  294.175,  294.325,  294.475, &
   294.600,  294.700,  294.800,  294.900,  295.025,  295.225,  295.500,  295.800, &
   296.125,  296.475,  296.825,  297.200,  297.625,  298.075,  298.500,  298.900, &
   299.300,  299.700,  300.075,  300.425,  300.775,  301.100,  301.400,  301.725, &
   302.075,  302.400,  302.700,  303.025,  303.400,  303.775,  304.125,  304.525, &
   304.975,  305.400,  305.825,  306.300,  306.775,  307.225,  307.700,  308.175, &
   308.600,  309.000,  309.400,  309.750,  310.000,  310.175,  310.300,  310.375, &
   310.375,  310.300,  310.200,  310.125,  310.100,  310.125,  310.200,  310.325, &
   310.500,  310.750,  311.100,  311.500,  311.925,  312.425,  313.000,  313.600, &
   314.225,  314.848,  315.500,  316.272,  317.075,  317.795,  318.397,  318.925, &
   319.647,  320.647,  321.605,  322.635,  323.902,  324.985,  325.855,  327.140, &
   328.677,  329.742,  330.585,  331.747,  333.272,  334.848,  336.525,  338.360, &
   339.728,  340.793,  342.198,  343.783,  345.283,  346.797,  348.645,  350.737, &
   352.487,  353.855,  355.017,  355.885,  356.777,  358.128,  359.837,  361.462, &
   363.155,  365.323,  367.348,  368.865,  370.467,  372.522,  374.760,  376.812, & !1997-2004
   378.812,  380.828,  382.777,  384.800,  386.935,  389.072,  391.167,  393.240, &
   395.297,  397.346,  399.387,  401.418,  403.431,  405.425,  407.400,  409.360, & !2013-2020
   411.297,  413.219,  415.144,  417.083,  419.036,  421.004,  422.978,  424.950, &
   426.916,  428.876,  430.832,  432.807,  434.831,  436.916,  439.068,  441.286, & !2029-2036
   443.567,  445.903,  448.282,  450.698,  453.150,  455.645,  458.182,  460.762, &
   463.405,  466.120,  468.907,  471.768,  474.692,  477.670,  480.697,  483.777, & !2045-2052
   486.916,  490.103,  493.338,  496.642,  500.022,  503.483,  507.023,  510.634, &
   514.305,  518.027,  521.797,  525.619,  529.486,  533.400,  537.381,  541.443, & !2061-2068
   545.589,  549.820,  554.129,  558.486,  562.867,  567.272,  571.701,  576.146, &
   580.606,  585.105,  589.653,  594.257,  598.918,  603.538,  608.020,  612.363, & !2077-2084
   616.572,  620.648,  624.583,  628.381,  632.065,  635.649,  639.141,  642.597, &
   646.061,  649.515,  652.951,  656.364,  659.754,  663.107,  666.423,  669.723, & !2093-2100
   673.022,  676.291,  679.501,  682.645,  685.711,  688.693,  691.589,  694.400, &
   697.114,  699.732,  702.276,  704.764,  707.202,  709.595,  711.934,  714.205, & !2109-2116
   716.401,  718.519,  720.560,  722.510,  724.368,  726.158,  727.896,  729.591, &
   731.244,  732.848,  734.389,  735.862,  737.263,  738.592,  739.834,  740.989, & !2125-2132
   742.080,  743.124,  744.128,  745.096,  746.019,  746.882,  747.677,  748.400, &
   749.052,  749.617,  750.093,  750.506,  750.871,  751.198,  751.489,  751.735, & !2141-2148
   751.922,  751.999,  751.999 /)

! Establish location in table for current, previous and next year.
! ----------------------------------------------------------------
  current  = Year-firstYear+1
  current  = MAX(current,1)
  current  = MIN(current,tableLength)

  previous = current-1
  previous = MAX(previous,1)
  previous = MIN(previous,tableLength)
  IF(Year > finalYear) previous = tableLength

  next     = current+1
  next     = MAX(next,1)
  next     = MIN(next,tableLength)
  IF(Year < firstYear) next = 1

! Divide the year into halves.
! ----------------------------
  IF(dayOfYear <= 182) THEN
   i = CO2ppmv(previous)
   f = CO2ppmv(current)
   n = dayOfYear+183
  ELSE
   i = CO2ppmv(current)
   f = CO2ppmv(next)
   n = dayOfYear-183
  END IF

! Linear interpolation to the given day-of-year.
! ----------------------------------------------
  GetCO2 = i + (f-i)*n/365.00

! Convert to mole fraction (volume mixing ratio).
! -----------------------------------------------
  GetCO2 = GetCO2*1.00E-06

 END FUNCTION GetCO2
