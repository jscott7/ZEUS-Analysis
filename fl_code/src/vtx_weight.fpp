C ===================================================
      Subroutine VTX_WEIGHT(ZVTX_TRUE,weight,
     &                      genver,rewver,input_init)
C ===================================================
C
C VTX reweighting routine to reweight the 94, 95, 96 or 97 MC 
C vertex distribution to the "true" measured one.
C
C INPUT:  ZVTX_TRU   -> generated z-vertex position in cm
C         INIT       -> 1 redo initialisation (only needed if 
C                         genver or rewver is changed while running)
C                       0 keep weights as chosen before
C         GENVER     -> version ID with which MC was generated
C                       941 = default for 1994 (MOZART 1994.1)
C                       950 = default for 1995 (MOZART 1995.0)
C                       961 = preliminary version for 1996(MOZART 1996.1)
C                       963 = default for 1996 (MOZART 1996.3)
C                       971 = default for 1997 (MOZART 1997.1)
C         REWVER     -> "true" vertex distribution to which you
C                       want to reweight (94, 95, 96 or 97)
C                 
C OUTPUT: weight     -> weight for this particular event
C
C  Version 1.2     minor change input to ensure that code
C                  compiles on every machine                03/04/1998
C
C  Version 1.1     added GENZ961                            02/04/1998
C  Version 1.0                                              20/02/1998
C
C A.Quadt & O.C.Ruske, 03.04.1998
C -------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER I, BIN_NR, genver,rewver, init, input_init
      Save    genver_old, rewver_old
      INTEGER genver_old, rewver_old
      DATA    genver_old, rewver_old /-99,-99/
      INTEGER bin_from, bin_to
      Save    Lfirst
      Logical Lfirst
      Data    Lfirst /.true./
      Save    norm_from, norm_to, nrb_from,nrb_to
      Save    zmax_from, zmax_to
      REAL    ZVTX_TRUE, weight, norm_from, norm_to, nrb_from,nrb_to
      REAL    zmax_from, zmax_to
      Save    meas94,       meas95,       meas96,       meas97
      Real    meas94(100),  meas95(100),  meas96(100),  meas97(100)
      Save    genz941,      genz950,      genz961,      genz963
      Save    genz971
      Real    genz941(200), genz950(150), genz961(100), genz963(100)
      Real    genz971(100)
      Save    WeightBin,      FromVtx,      ToVtx
      Real    WeightBin(200), FromVtx(200), ToVtx(200)

C  ===========================================================
C  Measured Z-distribution for all NVTX 1994 e+ EVTAKE runs ==
C  ===========================================================
      DATA MEAS94/
     &2253.616,    2772.894,    3351.101,     3977.799,     4637.664,    
     &5310.759,    5973.310,    5375.280,     9178.290,     10587.50,    
     &8407.360,    6265.280,    7356.890,     8099.640,     7256.050,    
     &11186.80,    7658.040,    7534.090,     5697.790,     7291.150,    
     &5927.467,    7455.276,    4021.941,     6965.573,     2280.805,    
     &7005.310,    7857.960,    7047.220,     8212.410,     7160.200,    
     &11856.47,    9811.567,    6179.087,     8398.749,     13113.03,    
     &10284.34,    12016.84,    16228.54,     24638.52,     31936.35,    
     &60287.22,    69228.40,    93075.77,     139897.9,     160183.1,    
     &211703.0,    262892.7,    299200.0,     352136.9,     367695.0,    
     &386798.0,    392920.0,    401509.9,     380594.2,     321290.1,    
     &298059.4,    240984.1,    236068.2,     152960.5,     113835.6,    
     &110088.9,    64152.11,    48067.30,     33249.91,     28504.94,    
     &16845.44,    13199.49,    13215.79,     8062.620,     6444.747,    
     &6709.360,    7331.500,    6718.550,     6538.150,     6184.320,    
     &5887.853,    5643.243,    8096.084,     8503.365,     7214.989,    
     &7620.790,    10583.97,    11541.80,     11732.49,     12976.55,    
     &15050.75,    13402.05,    10694.12,     11455.63,     8322.610,    
     &6003.127,    5549.040,    4563.369,     4644.078,     3873.537,    
     &3285.053,    2233.732,    3102.934,     2073.627,     1183.998/    



C  ========================================================
C  Measured Z-distribution for all NVTX 1995 EVTAKE runs ==
C  ========================================================
      DATA MEAS95/
     &1400.092,    1998.451,     2757.402,     3677.697,    4741.559,    
     &5909.297,    7119.014,     10055.60,     8909.770,    11335.00,   
     &10506.20,    11172.80,     9813.810,     11973.30,    8413.410,    
     &8199.690,    6832.660,     6312.740,     4517.930,    3705.390,    
     &3617.010,    5121.481,     7328.377,     5846.690,    6119.808,    
     &8864.630,    10158.10,     8636.650,     11525.10,    12161.10,    
     &14369.82,    15126.54,     18976.22,     20694.10,    27991.00,    
     &32853.72,    47845.44,     69966.48,     92986.34,    118343.6,    
     &171052.9,    194775.9,     251306.4,     296319.8,    324992.9,    
     &386444.0,    453032.1,     418799.0,     461921.8,    470572.8,    
     &459054.0,    394037.1,     367554.9,     325387.5,    284211.8,    
     &225379.2,    196587.1,     141269.8,     101947.5,    78004.27,    
     &59126.25,    35091.99,     31231.27,     20331.42,    17660.28,    
     &13143.76,    9780.028,     9748.893,     9678.068,    9051.649,    
     &7211.760,    6471.380,     6746.660,     6881.380,    8710.100,    
     &7391.978,    7581.057,     9019.000,     8299.194,    9536.899,    
     &11168.12,    10901.65,     10241.77,     10166.15,    9635.363,    
     &9140.095,    8688.360,     7302.819,     5781.897,    6476.838,    
     &4684.256,    5145.475,     4051.022,     3958.733,    3173.436,    
     &2504.409,    1922.743,     2460.962,     1703.078,    949.9960/    


C  ========================================================
C  Measured Z-distribution for all NVTX 1996 EVTAKE runs ==
C  ========================================================
      DATA MEAS96/
     &3077.634,    4328.571,     5918.110,     7865.613,    10162.33,    
     &12763.36,    15582.88,     19155.90,     20326.50,    25175.60,    
     &27249.80,    28402.40,     27896.80,     28070.20,    28912.60,    
     &28952.40,    20304.00,     21441.20,     21053.10,    15449.50,    
     &17632.92,    15364.44,     11864.99,     12083.16,    16559.48,    
     &16399.60,    19151.10,     18099.00,     17770.80,    20901.20,    
     &18231.70,    25113.94,     30812.88,     36933.86,    48212.26,    
     &53991.52,    72219.53,     112622.3,     144722.4,    187638.6,    
     &280786.8,    312197.0,     407826.3,     495234.1,    549975.0,    
     &651727.8,    733226.3,     734187.0,     770264.3,    765316.1,    
     &715932.0,    636808.4,     595553.5,     524728.8,    457364.0,    
     &350559.8,    296207.7,     222493.0,     156201.4,    108320.0,    
     &96181.39,    60206.02,     41144.57,     33285.95,    25926.79,    
     &21580.63,    18678.97,     17750.34,     17382.94,    16021.42,    
     &15224.60,    14125.50,     13308.90,     15697.20,    18307.40,    
     &19532.45,    21020.20,     25253.20,     29047.18,    27029.20,    
     &27108.95,    36082.92,     29820.23,     29712.89,    27066.61,    
     &28029.63,    22992.95,     20589.24,     19837.89,    17927.89,    
     &12104.86,    11577.32,     8800.497,     6940.002,    5320.172,    
     &4960.657,    4165.944,     2846.996,     2255.428,    2327.490/    

C  ========================================================
C  Measured Z-distribution for all NVTX 1997 EVTAKE runs ==
C  ========================================================
      DATA MEAS97/
     &454.4782,    625.5173,     851.9411,     1148.215,    1531.372,
     &2021.073,    2639.529,     4498.350,     6408.490,    4871.900,
     &7341.990,    7532.750,     12279.00,     12551.00,    16602.80,
     &18452.90,    21526.60,     25559.10,     28554.30,    25837.40,
     &21634.93,    29865.21,     24230.90,     29487.62,    38087.41,
     &29224.70,    30275.30,     32065.10,     32065.10,    32893.00,
     &29329.45,    37007.73,     36282.48,     49529.93,    60172.87,
     &68668.46,    88470.56,     110471.3,     172039.8,    212925.7,
     &303320.0,    336487.9,     381415.8,     516386.8,    539362.0,
     &637074.3,    750307.4,     736348.9,     776420.3,    759701.8,
     &774530.4,    659178.0,     595411.8,     553047.3,    456119.0,
     &356898.1,    308234.0,     249850.3,     180879.2,    135685.2,
     &108409.0,    73638.17,     50762.87,     40345.84,    35374.24,
     &31253.31,    23566.51,     23458.65,     19092.42,    21737.77,
     &22591.30,    18362.50,     17297.10,     13552.20,    14742.90,
     &14738.06,    14410.27,     14237.54,     16725.56,    15631.31,
     &15409.08,    16394.14,     16506.64,     16576.87,    16505.72,
     &14831.28,    13036.32,     12677.44,     12246.52,    12150.92,
     &10586.56,    11162.15,     9814.348,     9887.035,    5220.086,
     &5597.416,    3342.341,     3227.303,     1712.133,    1551.577/



C  ====================================
C  === GENZ94   MOZART distribution ===
C  ====================================
      DATA genz941/
     +       254.0,     246.0,     272.0,       263.0,       268.0,
     +       277.0,     263.0,     269.0,       301.0,       289.0,
     +       283.0,     311.0,     320.0,       330.0,       353.0,
     +       362.0,     347.0,     419.0,       378.0,       357.0,
     +       419.0,     391.0,     375.0,       342.0,       395.0,
     +       375.0,     411.0,     395.0,       395.0,       426.0,
     +       399.0,     423.0,     421.0,       409.0,       445.0,
     +       428.0,     425.0,     443.0,       440.0,       428.0,
     +       460.0,     436.0,     440.0,       461.0,       432.0,
     +       392.0,     416.0,     461.0,       453.0,       459.0,
     +       423.0,     431.0,     433.0,       452.0,       439.0,
     +       484.0,     478.0,     485.0,       497.0,       449.0,
     +       494.0,     504.0,     486.0,       541.0,       483.0,
     +       542.0,     549.0,     549.0,       544.0,       588.0,
     +       629.0,     642.0,     733.0,       831.0,       870.0,
     +       967.0,    1131.0,    1303.0,      1530.0,      1793.0,
     +      1965.0,    2393.0,    2601.0,      3120.0,      3527.0,
     +      4074.0,    4752.0,    5313.0,      6122.0,      6939.0,
     +      7700.0,    8371.0,    9488.0,     10360.0,     11070.0,
     +     12019.0,   13080.0,   13691.0,     14578.0,     15007.0,
     +     15433.0,   15829.0,   16127.0,     15872.0,     15785.0,
     +     15803.0,   15249.0,   14874.0,     14446.0,     13680.0,
     +     12755.0,   11911.0,   10949.0,     10122.0,      9149.0,
     +      8359.0,    7443.0,    6634.0,      5920.0,      5174.0,
     +      4425.0,    3906.0,    3239.0,      2767.0,      2378.0,
     +      2021.0,    1773.0,    1509.0,      1302.0,      1134.0,
     +      1007.0,     841.0,     751.0,       659.0,       644.0,
     +       549.0,     504.0,     442.0,       472.0,       449.0,
     +       437.0,     462.0,     400.0,       411.0,       431.0,
     +       380.0,     409.0,     366.0,       389.0,       416.0,
     +       393.0,     413.0,     395.0,       395.0,       440.0,
     +       390.0,     444.0,     438.0,       444.0,       456.0,
     +       478.0,     473.0,     520.0,       556.0,       603.0,
     +       572.0,     558.0,     578.0,       616.0,       614.0,
     +       656.0,     549.0,     569.0,       560.0,       569.0,
     +       485.0,     500.0,     448.0,       422.0,       406.0,
     +       351.0,     355.0,     285.0,       266.0,       227.0,
     +       205.0,     201.0,     177.0,       205.0,       160.0,
     +       149.0,     157.0,     134.0,       110.0,       124.0,
     +        83.0,      79.0,      70.0,        75.0,        51.0 /


C  ====================================
C  === GENZ95.0 MOZART distribution ===
C  ====================================
      DATA genz950/
     .    0,     2,     0,     1,     3,
     .    0,     2,     1,     1,     3,
     .    5,     1,     3,     1,     3,
     .    5,     3,     4,     2,     3,
     .    3,     3,     4,     4,     4,
     .    2,     6,     5,     8,     4,
     .    4,     4,     4,    12,     6,
     .    8,    17,    18,    23,    35,
     .   49,    57,    67,    45,    74,
     .  185,   235,   314,   360,   387,
     .  442,   425,   565,   729,   783,
     .  938,  1128,  1296,  1427,  1685,
     . 2030,  2557,  3620,  5427,  8238,
     .12874, 19344, 28008, 39029, 51466,
     .65284, 77281, 86538, 92192, 92252,
     .88513, 79067, 66675, 52949, 39920,
     .28547, 19574, 12840,  8319,  5474,
     . 3555,  2389,  1836,  1437,  1334,
     . 1261,  1055,  1003,   966,   942,
     . 1044,  1046,   944,   784,   730,
     .  604,   527,   408,   301,   199,
     .   87,    61,    59,    42,    27,
     .   36,    28,    17,    14,    11,
     .    5,    11,     3,     9,     6,
     .    8,     4,     3,     6,     0,
     .    4,     3,     4,     3,     4,
     .    3,     2,     3,     1,     2,
     .    1,     2,     1,     1,     2,
     .    1,     0,     1,     1,     2,
     .    2,     0,     1,     2,     1/

C  ====================================
C  === GENZ96.1 MOZART distribution ===
C  ====================================
      DATA genz961/
     +   1.0,      2.4,      2.6,      3.5,      4.0,
     +   4.9,      6.7,      7.4,      10.5,     16.3,
     +  19.1,     21.8,     22.2,     23.1,     23.6,
     +  24.0,     26.7,     26.0,     25.2,     25.1,
     +  25.6,     24.0,     23.5,     22.6,     20.5,
     +  18.9,     17.7,     18.7,     21.9,     22.9,
     +  30.4,     37.5,     41.1,     52.1,     66.9,
     +  85.1,    107.9,    180.8,    242.6,    327.0,
     + 489.2,    560.5,    683.7,    799.8,    903.4,
     +1051.0,   1181.4,   1185.5,   1208.7,   1144.7,
     +1111.6,   1029.2,    906.7,    779.8,    625.6,
     + 477.5,    361.9,    311.5,    223.3,    155.3,
     + 110.2,     66.9,     54.8,     46.7,     38.2,
     +  26.4,     22.2,     19.0,     18.6,     18.2,
     +  17.7,     17.8,     18.5,     20.0,     22.4,
     +  23.3,     24.2,     28.8,     33.2,     34.8,
     +  34.7,     33.6,     33.1,     31.5,     28.4,
     +  25.0,     22.8,     18.6,     13.0,      9.3,
     +   8.8,      7.6,      6.5,      4.9,      4.8,
     +   4.0,      3.9,      2.8,      2.2,      1.2/

C  ====================================
C  === GENZ96.3 MOZART distribution ===
C  ====================================
      DATA genz963/
     &  0.0     , 0.0     , 8.1     , 0.0     , 
     &  12.0    , 30.6    , 0.0     , 21.2    , 
     &  5.1     , 22.8    , 48.2    , 34.8    , 
     &  17.2    , 13.4    , 38.0    , 42.3    , 
     &  48.0    , 20.3    , 23.0    , 31.3    , 
     &  15.3    , 13.6    , 11.9    , 7.9     , 
     &  28.7    , 25.6    , 23.4    , 24.3    , 
     &  28.8    , 29.2    , 26.9    , 37.0    , 
     &  44.7    , 53.4    , 57.0    , 89.8    , 
     &  125.4   , 160.0   , 241.9   , 299.5   , 
     &  404.4   , 524.5   , 640.3   , 782.5   , 
     &  879.7   , 1001.9  , 1062.2  , 1143.4  , 
     &  1173.4  , 1119.6  , 1081.8  , 952.3   , 
     &  935.7   , 809.5   , 644.3   , 532.4   , 
     &  391.3   , 297.9   , 235.4   , 148.9   , 
     &  113.8   , 82.0    , 55.7    , 40.7    , 
     &  39.4    , 26.7    , 21.0    , 18.1    , 
     &  23.7    , 19.0    , 17.6    , 25.2    , 
     &  15.9    , 23.9    , 22.4    , 29.8    , 
     &  29.4    , 36.0    , 38.7    , 37.0    , 
     &  34.4    , 40.9    , 46.3    , 51.4    , 
     &  43.7    , 42.8    , 28.6    , 29.6    , 
     &  21.7    , 18.7    , 15.1    , 13.1    , 
     &  14.3    , 6.7     , 8.4     , 4.8     , 
     &  1.6     , 1.5     , 3.0     , 0.8     /

C  ====================================
C  === GENZ97.1 MOZART distribution ===
C  ====================================
      DATA genz971/
     &  0.0     , 0.0     , 0.0     , 4.2     , 
     &  9.9     , 10.3    , 5.5     , 11.6    , 
     &  5.6     , 12.5    , 34.0    , 27.6    , 
     &  13.5    , 10.3    , 14.2    , 28.2    , 
     &  29.9    , 16.5    , 28.5    , 15.3    , 
     &  24.3    , 21.9    , 26.3    , 14.4    , 
     &  47.5    , 41.0    , 38.2    , 50.6    , 
     &  48.0    , 45.7    , 46.3    , 50.2    , 
     &  52.3    , 70.3    , 79.6    , 106.8   , 
     &  135.3   , 153.8   , 243.4   , 305.1   , 
     &  382.5   , 499.5   , 598.6   , 728.0   , 
     &  831.8   , 984.5   , 1044.0  , 1121.6  , 
     &  1134.9  , 1108.1  , 1095.9  , 983.4   , 
     &  958.8   , 847.4   , 668.3   , 556.3   , 
     &  427.1   , 325.6   , 276.5   , 191.0   , 
     &  136.3   , 100.4   , 79.7    , 51.0    , 
     &  45.8    , 34.5    , 30.6    , 23.9    , 
     &  27.9    , 22.7    , 23.8    , 25.6    , 
     &  18.9    , 18.5    , 18.4    , 20.3    , 
     &  13.3    , 14.2    , 19.2    , 15.5    , 
     &  16.9    , 15.7    , 19.2    , 21.8    , 
     &  19.8    , 13.4    , 11.8    , 15.4    , 
     &  11.1    , 10.3    , 10.5    , 7.2     , 
     &  11.2    , 8.3     , 5.9     , 3.1     , 
     &  1.4     , 1.1     , 0.9     , 0.5     /


C ----------------------
C --- Initialisation ---
C ----------------------
      init = input_init

      IF (genver_old.NE.genver .OR.
     &    rewver_old.NE.rewver) THEN
        genver_old = genver
        rewver_old = rewver
        init       = 1
      ENDIF
      If(Lfirst.OR.Init.EQ.1) Then
        write(6,*) ' '
        Lfirst    = .False.
        Init      = 0
        norm_from = 0.0
        norm_to   = 0.0

        IF (genver.EQ.941) THEN
          nrb_from  = 200.0
          zmax_from = 100.0
          DO I=1,INT(nrb_from)
            FromVtx(I)    = GENZ941(I)
            norm_from     = norm_from + 
     &                      FromVtx(I)*2.0*zmax_from/nrb_from
          ENDDO
        ENDIF
        IF (genver.EQ.950) THEN
          nrb_from  = 150.0
          zmax_from = 200.0
          DO I=1,INT(nrb_from)
            FromVtx(I)    = GENZ950(I)
            norm_from     = norm_from + 
     &                      FromVtx(I)*2.0*zmax_from/nrb_from
          ENDDO
        ENDIF
        IF (genver.EQ.961) THEN
          nrb_from  = 100.0
          zmax_from = 100.0
          DO I=1,INT(nrb_from)
            FromVtx(I)    = GENZ961(I)
            norm_from     = norm_from + 
     &                      FromVtx(I)*2.0*zmax_from/nrb_from
          ENDDO
        ENDIF
        IF (genver.EQ.963) THEN
          nrb_from  = 100.0
          zmax_from = 100.0
          DO I=1,INT(nrb_from)
            FromVtx(I)    = GENZ963(I)
            norm_from     = norm_from + 
     &                      FromVtx(I)*2.0*zmax_from/nrb_from
          ENDDO
        ENDIF
        IF (genver.EQ.971) THEN
          nrb_from  = 100.0
          zmax_from = 100.0
          DO I=1,INT(nrb_from)
            FromVtx(I)    = GENZ971(I)
            norm_from     = norm_from + 
     &                      FromVtx(I)*2.0*zmax_from/nrb_from
          ENDDO
        ENDIF

        IF (rewver.eq.94) THEN
          nrb_to    = 100.0
          zmax_to   = 100.0
          DO I=1,INT(nrb_to)
            ToVtx(I)     = MEAS94(I)
            norm_to      = norm_to   + 
     &                     ToVtx(I)*2.0*zmax_to/nrb_to
          ENDDO
        ENDIF
        IF (rewver.eq.95) THEN
          nrb_to    = 100.0
          zmax_to   = 100.0
          DO I=1,INT(nrb_to)
            ToVtx(I)     = MEAS95(I)
            norm_to      = norm_to   + 
     &                     ToVtx(I)*2.0*zmax_to/nrb_to
          ENDDO
        ENDIF
        IF (rewver.eq.96) THEN
          nrb_to    = 100.0
          zmax_to   = 100.0
          DO I=1,INT(nrb_to)
            ToVtx(I)     = MEAS96(I)
            norm_to      = norm_to   + 
     &                     ToVtx(I)*2.0*zmax_to/nrb_to
          ENDDO
        ENDIF
        IF (rewver.eq.97) THEN
          nrb_to    = 100.0
          zmax_to   = 100.0
          DO I=1,INT(nrb_to)
            ToVtx(I)     = MEAS97(I)
            norm_to      = norm_to   + 
     &                     ToVtx(I)*2.0*zmax_to/nrb_to
          ENDDO
        ENDIF

        IF ((norm_to.LE.0.OR.norm_from.LE.0.0).AND.rewver.NE.97) THEN
          WRITE(*,*) 'WARNING: ERROR in VTX weighting routine',
     &               norm_to, norm_from
          STOP
        ENDIF
        WRITE(*,*)
        WRITE(*,91) 
        WRITE(*,93)
        WRITE(*,94)
        WRITE(*,95)
        WRITE(*,96)
        WRITE(*,93)
        WRITE(*,92) GENVER, rewver
        WRITE(*,93)
        WRITE(*,91)
        WRITE(999,91) 
        WRITE(999,93)
        WRITE(999,94)
        WRITE(999,95)
        WRITE(999,96)
        WRITE(999,93)
        WRITE(999,92) GENVER, rewver
        WRITE(999,93)
        WRITE(999,91)

        WRITE(*,*)
 91     FORMAT('  +--------------------------------+')
 92     FORMAT('  | reweighted from ',I3,'  to  ',I2,'    |')
 93     FORMAT('  |                                |')
 94     FORMAT('  | VTX reweighting routine 94-97  |')
 95     FORMAT('  | by A.Quadt & O.C.Ruske         |')
 96     FORMAT('  | (quadt@desy.de, ruske@desy.de) |')
      Endif

C ---------------------------------------
C --- Determine weight for this event ---
C ---------------------------------------
      IF (ABS(ZVTX_TRUE).GT.100.0) THEN
        weight = 1.0
        RETURN

      ELSE
        bin_from = 1+INT((ZVTX_TRUE+zmax_from)/
     &                   (2.0*zmax_from/nrb_from))
        bin_to   = 1+INT((ZVTX_TRUE+zmax_to)/
     &                   (2.0*zmax_to/nrb_to))
        IF (BIN_from.LE. 0)            BIN_from =           1
        IF (BIN_to  .LE. 0)            BIN_to   =           1
        IF (BIN_from.GT.INT(nrb_from)) BIN_from = INT(nrb_from)
        IF (BIN_to  .GT.INT(nrb_to))   BIN_to   = INT(nrb_to)
        IF (FromVtx(bin_from).NE.0.0) THEN
          Weight   = (ToVtx  (bin_to)  /norm_to)    / 
     &               (FromVtx(bin_from)/norm_from)
        ELSE
          Weight   = 1.0
        ENDIF

      ENDIF

      Return
      End

