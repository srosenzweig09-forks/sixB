COM_ENERGY = 13000.0 # GeV
CROSS_SECTION = 1 # pb
PROCESS = 'HiggsBSM:gg2H2 = on'
SLHA_TABLE = """BLOCK SPINFO
     1   FeynHiggs
     2   2.12.0
     2   built on ott 13, 2016
BLOCK MODSEL
         1                  0   # Model
         2                  1   # GridPts
         3                  0   # Content
         4                  0   # RPV
         5                  0   # CPV
         6                  0   # FV
BLOCK SMINPUTS
         1     1.28952828E+02   # invAlfaMZ
         2     1.16637000E-05   # GF
         3     1.19000000E-01   # AlfasMZ
         4     9.11876000E+01   # MZ
         5     4.16000000E+00   # Mb
         6     1.73200000E+02   # Mt
         7     1.77703000E+00   # Mtau
        11     5.10998902E-04   # Me
        13     1.05658357E-01   # Mmu
        21     6.00000000E-03   # Md
        22     3.00000000E-03   # Mu
        23     9.50000000E-02   # Ms
        24     1.28600000E+00   # Mc
BLOCK MINPAR
         3     5.00000000E+00   # TB
BLOCK EXTPAR
         0     0.00000000E+00   # Q
         1     9.54716519E+01   # M1
         2     2.00000000E+02   # M2
         3     1.50000000E+03   # M3
        11     1.54000000E+03   # At
        12     1.54000000E+03   # Ab
        13     1.54000000E+03   # Atau
        23     2.00000000E+02   # MUE
        25     5.00000000E+00   # TB
        26     1.40000000E+02   # MA0
        27     1.61442988E+02   # MHp
        31     5.00000000E+02   # MSL(1)
        32     5.00000000E+02   # MSL(2)
        33     1.00000000E+03   # MSL(3)
        34     5.00000000E+02   # MSE(1)
        35     5.00000000E+02   # MSE(2)
        36     1.00000000E+03   # MSE(3)
        41     1.50000000E+03   # MSQ(1)
        42     1.50000000E+03   # MSQ(2)
        43     1.00000000E+03   # MSQ(3)
        44     1.50000000E+03   # MSU(1)
        45     1.50000000E+03   # MSU(2)
        46     1.00000000E+03   # MSU(3)
        47     1.50000000E+03   # MSD(1)
        48     1.50000000E+03   # MSD(2)
        49     1.00000000E+03   # MSD(3)
BLOCK MASS
   1000012     4.96147383E+02   # MSf(1,1,1)
   1000011     5.02124333E+02   # MSf(1,2,1)
   2000011     5.01706018E+02   # MSf(2,2,1)
   1000002     1.49910023E+03   # MSf(1,3,1)
   2000002     1.49962019E+03   # MSf(2,3,1)
   1000001     1.50108899E+03   # MSf(1,4,1)
   2000001     1.50018986E+03   # MSf(2,4,1)
   1000014     4.96147383E+02   # MSf(1,1,2)
   1000013     5.02149324E+02   # MSf(1,2,2)
   2000013     5.01681028E+02   # MSf(2,2,2)
   1000004     1.49910022E+03   # MSf(1,3,2)
   2000004     1.49962131E+03   # MSf(2,3,2)
   1000003     1.50109009E+03   # MSf(1,4,2)
   2000003     1.50018877E+03   # MSf(2,4,2)
   1000016     9.98079268E+02   # MSf(1,1,3)
   1000015     1.00046976E+03   # MSf(1,2,3)
   2000015     1.00145112E+03   # MSf(2,2,3)
   1000006     8.76514793E+02   # MSf(1,3,3)
   2000006     1.13484824E+03   # MSf(2,3,3)
   1000005     9.99672858E+02   # MSf(1,4,3)
   2000005     1.00226018E+03   # MSf(2,4,3)
        25     1.10184245E+02   # Mh0
        35     1.53285865E+02   # MHH
        36     1.40000000E+02   # MA0
        37     1.61075400E+02   # MHp
   1000022     8.32246583E+01   # MNeu(1)
   1000023     1.46162492E+02   # MNeu(2)
   1000025    -2.07032524E+02   # MNeu(3)
   1000035     2.73117025E+02   # MNeu(4)
   1000024     1.38016688E+02   # MCha(1)
   1000037     2.71807047E+02   # MCha(2)
   1000021     1.50000000E+03   # MGl
BLOCK DMASS
         0     1.73200000E+02   # Q
        25     8.39496840E-01   # Delta Mh0
        35     2.27280933E-01   # Delta MHH
        36     0.00000000E+00   # Delta MA0
        37     1.63478690E-01   # Delta MHp
BLOCK NMIX
     1   1     9.02741035E-01   # ZNeu(1,1)
     1   2    -1.76906456E-01   # ZNeu(1,2)
     1   3     3.40847919E-01   # ZNeu(1,3)
     1   4    -1.93869609E-01   # ZNeu(1,4)
     2   1    -3.94247799E-01   # ZNeu(2,1)
     2   2    -6.93610991E-01   # ZNeu(2,2)
     2   3     4.67997695E-01   # ZNeu(2,3)
     2   4    -3.80066603E-01   # ZNeu(2,4)
     3   1     8.11445475E-02   # ZNeu(3,1)
     3   2    -1.12684390E-01   # ZNeu(3,2)
     3   3    -6.81879284E-01   # ZNeu(3,3)
     3   4    -7.18163235E-01   # ZNeu(3,4)
     4   1    -1.51798743E-01   # ZNeu(4,1)
     4   2     6.89137234E-01   # ZNeu(4,2)
     4   3     4.47036346E-01   # ZNeu(4,3)
     4   4    -5.49732225E-01   # ZNeu(4,4)
BLOCK UMIX
     1   1    -6.25444296E-01   # UCha(1,1)
     1   2     7.80268821E-01   # UCha(1,2)
     2   1     7.80268821E-01   # UCha(2,1)
     2   2     6.25444296E-01   # UCha(2,2)
BLOCK VMIX
     1   1    -7.80268821E-01   # VCha(1,1)
     1   2     6.25444296E-01   # VCha(1,2)
     2   1     6.25444296E-01   # VCha(2,1)
     2   2     7.80268821E-01   # VCha(2,2)
BLOCK STAUMIX
     1   1    -6.27000005E-01   # USf(1,1)
     1   2     7.79019251E-01   # USf(1,2)
     2   1     7.79019251E-01   # USf(2,1)
     2   2     6.27000005E-01   # USf(2,2)
BLOCK STOPMIX
     1   1     7.08166916E-01   # USf(1,1)
     1   2    -7.06045055E-01   # USf(1,2)
     2   1     7.06045055E-01   # USf(2,1)
     2   2     7.08166916E-01   # USf(2,2)
BLOCK SBOTMIX
     1   1    -4.89397966E-01   # USf(1,1)
     1   2     8.72060566E-01   # USf(1,2)
     2   1     8.72060566E-01   # USf(2,1)
     2   2     4.89397966E-01   # USf(2,2)
BLOCK ALPHA
              -6.63896461E-01   # Alpha
BLOCK DALPHA
               2.25792252E-02   # Delta Alpha
BLOCK HMIX Q= -0.99900000E+03
         1     2.00000000E+02   # MUE
         2     5.00000000E+00   # TB
BLOCK MSOFT Q=  0.00000000E+00
         1     9.54716519E+01   # M1
         2     2.00000000E+02   # M2
         3     1.50000000E+03   # M3
        31     5.00000000E+02   # MSL(1)
        32     5.00000000E+02   # MSL(2)
        33     1.00000000E+03   # MSL(3)
        34     5.00000000E+02   # MSE(1)
        35     5.00000000E+02   # MSE(2)
        36     1.00000000E+03   # MSE(3)
        41     1.50000000E+03   # MSQ(1)
        42     1.50000000E+03   # MSQ(2)
        43     1.00000000E+03   # MSQ(3)
        44     1.50000000E+03   # MSU(1)
        45     1.50000000E+03   # MSU(2)
        46     1.00000000E+03   # MSU(3)
        47     1.50000000E+03   # MSD(1)
        48     1.50000000E+03   # MSD(2)
        49     1.00000000E+03   # MSD(3)
BLOCK AE Q=  0.00000000E+00
     1   1     0.00000000E+00   # Af(1,1)
     2   2     0.00000000E+00   # Af(2,2)
     3   3     1.54000000E+03   # Af(3,3)
BLOCK AU Q=  0.00000000E+00
     1   1     0.00000000E+00   # Af(1,1)
     2   2     0.00000000E+00   # Af(2,2)
     3   3     1.54000000E+03   # Af(3,3)
BLOCK AD Q=  0.00000000E+00
     1   1     0.00000000E+00   # Af(1,1)
     2   2     0.00000000E+00   # Af(2,2)
     3   3     1.54000000E+03   # Af(3,3)
BLOCK YE Q=  0.00000000E+00
     1   1     1.49657094E-05   # Yf(1,1)
     2   2     3.09443379E-03   # Yf(2,2)
     3   3     5.20441717E-02   # Yf(3,3)
BLOCK YU Q=  0.00000000E+00
     1   1     1.75722993E-05   # Yf(1,1)
     2   2     7.53265896E-03   # Yf(2,2)
     3   3     1.01450741E+00   # Yf(3,3)
BLOCK YD Q=  0.00000000E+00
     1   1     1.74829084E-04   # Yf(1,1)
     2   2     2.76809607E-03   # Yf(2,2)
     3   3     1.19892640E-01   # Yf(3,3)
BLOCK VCKMIN
         1     2.25300000E-01   # lambda
         2     8.08000000E-01   # A
         3     1.32000000E-01   # rhobar
         4     3.41000000E-01   # etabar
BLOCK MSL2 Q=  0.00000000E+00
     1   1     2.50000000E+05   # MSL2(1,1)
     2   2     2.50000000E+05   # MSL2(2,2)
     3   3     1.00000000E+06   # MSL2(3,3)
BLOCK MSE2 Q=  0.00000000E+00
     1   1     2.50000000E+05   # MSE2(1,1)
     2   2     2.50000000E+05   # MSE2(2,2)
     3   3     1.00000000E+06   # MSE2(3,3)
BLOCK MSQ2 Q=  0.00000000E+00
     1   1     2.25000000E+06   # MSQ2(1,1)
     2   2     2.25000000E+06   # MSQ2(2,2)
     3   3     1.00000000E+06   # MSQ2(3,3)
BLOCK MSU2 Q=  0.00000000E+00
     1   1     2.25000000E+06   # MSU2(1,1)
     2   2     2.25000000E+06   # MSU2(2,2)
     3   3     1.00000000E+06   # MSU2(3,3)
BLOCK MSD2 Q=  0.00000000E+00
     1   1     2.25000000E+06   # MSD2(1,1)
     2   2     2.25000000E+06   # MSD2(2,2)
     3   3     1.00000000E+06   # MSD2(3,3)
BLOCK TE Q=  0.00000000E+00
     1   1     0.00000000E+00   # Tf(1,1)
     2   2     0.00000000E+00   # Tf(2,2)
     3   3     8.01480244E+01   # Tf(3,3)
BLOCK TU Q=  0.00000000E+00
     1   1     0.00000000E+00   # Tf(1,1)
     2   2     0.00000000E+00   # Tf(2,2)
     3   3     1.56234141E+03   # Tf(3,3)
BLOCK TD Q=  0.00000000E+00
     1   1     0.00000000E+00   # Tf(1,1)
     2   2     0.00000000E+00   # Tf(2,2)
     3   3     1.84634666E+02   # Tf(3,3)
BLOCK SELMIX
     1   1     9.99999260E-01   # UASf(1,1)
     1   4    -1.21690415E-03   # UASf(1,4)
     2   2     9.72951221E-01   # UASf(2,2)
     2   5    -2.31010653E-01   # UASf(2,5)
     3   3    -6.27000005E-01   # UASf(3,3)
     3   6     7.79019251E-01   # UASf(3,6)
     4   1     1.21690415E-03   # UASf(4,1)
     4   4     9.99999260E-01   # UASf(4,4)
     5   2     2.31010653E-01   # UASf(5,2)
     5   5     9.72951221E-01   # UASf(5,5)
     6   3     7.79019251E-01   # UASf(6,3)
     6   6     6.27000005E-01   # UASf(6,6)
BLOCK USQMIX
     1   1     9.99999997E-01   # UASf(1,1)
     1   4     7.69624283E-05   # UASf(1,4)
     2   2     9.99457413E-01   # UASf(2,2)
     2   5     3.29375163E-02   # UASf(2,5)
     3   3     7.08166916E-01   # UASf(3,3)
     3   6    -7.06045055E-01   # UASf(3,6)
     4   1    -7.69624283E-05   # UASf(4,1)
     4   4     9.99999997E-01   # UASf(4,4)
     5   2    -3.29375163E-02   # UASf(5,2)
     5   5     9.99457413E-01   # UASf(5,5)
     6   3     7.06045055E-01   # UASf(6,3)
     6   6     7.08166916E-01   # UASf(6,6)
BLOCK DSQMIX
     1   1     9.99997553E-01   # UASf(1,1)
     1   4    -2.21213987E-03   # UASf(1,4)
     2   2     9.99388672E-01   # UASf(2,2)
     2   5    -3.49611623E-02   # UASf(2,5)
     3   3    -4.89397966E-01   # UASf(3,3)
     3   6     8.72060566E-01   # UASf(3,6)
     4   1     2.21213987E-03   # UASf(4,1)
     4   4     9.99997553E-01   # UASf(4,4)
     5   2     3.49611623E-02   # UASf(5,2)
     5   5     9.99388672E-01   # UASf(5,5)
     6   3     8.72060566E-01   # UASf(6,3)
     6   6     4.89397966E-01   # UASf(6,6)
BLOCK CVHMIX
     1   1     9.65448398E-01   # UH(1,1)
     1   2     2.60594304E-01   # UH(1,2)
     1   3     0.00000000E+00   # UH(1,3)
     2   1    -2.60594304E-01   # UH(2,1)
     2   2     9.65448398E-01   # UH(2,2)
     2   3     0.00000000E+00   # UH(2,3)
     3   1     0.00000000E+00   # UH(3,1)
     3   2     0.00000000E+00   # UH(3,2)
     3   3     1.00000000E+00   # UH(3,3)
DECAY        25     2.15784596E-02   # Gamma(h0)
     2.48701500E-04   2        22        22   # BR(h0 -> photon photon)
     3.09762530E-05   2        22        23   # BR(h0 -> photon Z)
     4.66178562E-04   2        23        23   # BR(h0 -> Z Z)
     4.92108626E-03   2       -24        24   # BR(h0 -> W W)
     4.62172352E-03   2        21        21   # BR(h0 -> gluon gluon)
     8.02598301E-09   2       -11        11   # BR(h0 -> Electron electron)
     3.56996356E-04   2       -13        13   # BR(h0 -> Muon muon)
     1.02488815E-01   2       -15        15   # BR(h0 -> Tau tau)
     2.17631612E-08   2        -2         2   # BR(h0 -> Up up)
     3.01402414E-03   2        -4         4   # BR(h0 -> Charm charm)
     1.33486419E-06   2        -1         1   # BR(h0 -> Down down)
     3.35230005E-04   2        -3         3   # BR(h0 -> Strange strange)
     8.83514903E-01   2        -5         5   # BR(h0 -> Bottom bottom)
DECAY        35     4.95432152E-02   # Gamma(HH)
    -6.49316549E-05   2        22        22   # BR(HH -> photon photon)
    -2.49868289E-04   2        22        23   # BR(HH -> photon Z)
    -7.81308519E-03   2        23        23   # BR(HH -> Z Z)
    -8.21428340E-02   2       -24        24   # BR(HH -> W W)
    -7.38973430E-03   2        21        21   # BR(HH -> gluon gluon)
    -7.83692054E-09   2       -11        11   # BR(HH -> Electron electron)
     3.48621678E-04   2       -13        13   # BR(HH -> Muon muon)
    -9.99685826E-02   2       -15        15   # BR(HH -> Tau tau)
    -7.51795949E-09   2        -2         2   # BR(HH -> Up up)
    -1.04156648E-03   2        -4         4   # BR(HH -> Charm charm)
    -1.22014736E-06   2        -1         1   # BR(HH -> Down down)
    -3.06422266E-04   2        -3         3   # BR(HH -> Strange strange)
    -8.00666631E-01   2        -5         5   # BR(HH -> Bottom bottom)
    -6.48731778E-06   2        23        36   # BR(HH -> Z A0)
DECAY        36     6.45284315E-02   # Gamma(A0)
     5.89232343E-06   2        22        22   # BR(A0 -> photon photon)
     2.49514053E-06   2        22        23   # BR(A0 -> photon Z)
     5.43852250E-04   2        21        21   # BR(A0 -> gluon gluon)
     8.51765916E-09   2       -11        11   # BR(A0 -> Electron electron)
     3.78894526E-04   2       -13        13   # BR(A0 -> Muon muon)
     1.08749042E-01   2       -15        15   # BR(A0 -> Tau tau)
     5.15940424E-10   2        -2         2   # BR(A0 -> Up up)
     7.21724477E-05   2        -4         4   # BR(A0 -> Charm charm)
     1.34985248E-06   2        -1         1   # BR(A0 -> Down down)
     3.38995870E-04   2        -3         3   # BR(A0 -> Strange strange)
     8.89845960E-01   2        -5         5   # BR(A0 -> Bottom bottom)
     6.13371698E-05   2        23        25   # BR(A0 -> Z h0)
DECAY        37     8.38900517E-03   # Gamma(Hp)
     7.94454390E-08   2       -11        12   # BR(Hp -> Electron nu_e)
     3.39653714E-03   2       -13        14   # BR(Hp -> Muon nu_mu)
     9.60534196E-01   2       -15        16   # BR(Hp -> Tau nu_tau)
     1.08807956E-05   2        -1         2   # BR(Hp -> Down up)
     1.21751076E-04   2        -3         2   # BR(Hp -> Strange up)
     7.58339279E-05   2        -5         2   # BR(Hp -> Bottom up)
     2.71691715E-05   2        -1         4   # BR(Hp -> Down charm)
     3.31890820E-03   2        -3         4   # BR(Hp -> Strange charm)
     1.06189156E-02   2        -5         4   # BR(Hp -> Bottom charm)
     1.26717620E-02   2        -5         6   # BR(Hp -> Bottom top)
     8.64828509E-03   2        24        25   # BR(Hp -> W h0)
     3.46732569E-06   2        24        35   # BR(Hp -> W HH)
     5.72214217E-04   2        24        36   # BR(Hp -> W A0)
DECAY         6     1.37288978E+00   # Gamma(top)
     9.98824052E-01   2         5        24   # BR(top -> bottom W)
     1.17594765E-03   2         5        37   # BR(top -> bottom Hp)
"""

import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
                         pythiaPylistVerbosity = cms.untracked.int32(1),                        
                         filterEfficiency = cms.untracked.double(1),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         SLHATableForPythia8 = cms.string('%s' % SLHA_TABLE),
                         comEnergy = cms.double(COM_ENERGY),
                         crossSection = cms.untracked.double(CROSS_SECTION),                         
                         maxEventsToPrint = cms.untracked.int32(1),
                         PythiaParameters = cms.PSet(
                             pythia8CommonSettingsBlock,
                             pythia8CUEP8M1SettingsBlock,
                             processParameters = cms.vstring(
                                 'Higgs:useBSM = on', 
                                 PROCESS, 
                                 'SLHA:allowUserOverride = off', 
                                 'SLHA:minMassSM = 100.', 
                                 'PhaseSpace:mHatMin = 56.0'
                             ),
                             parameterSets = cms.vstring(
                                 'pythia8CommonSettings',
                                 'pythia8CUEP8M1Settings',
                                 'processParameters'
                                 )
                             )
                         )

ProductionFilterSequence = cms.Sequence(generator)
