チュートリアル
=================

Heisenberg模型の計算
---------------------

ここではHeisenberg模型のGreen関数の計算を実施します．
入力ファイルとしてapp/sample/denovoにあるnamelist.defを使用します．
このファイルは以下の通りです．

.. code-block:: none

  &filename
    inham = ""
    invec = ""
  /
  &ham
  Jx = 1d0
  Jy = 0d0
  Jz = 1d0
  Dz = 1d0
  /
  &cg
    maxloops = 100
    convfactor = 6
  /
  &dyn
    calctype = "normal"
    nomega = 100
  !  omegamin = (-2d0, 0.1d0)
  !  omegamax = ( 1d0, 0.1d0)
  !  omegamin = (2d0, 0d0)
  !  omegamax = (3d0, 0d0)
    outrestart = .TRUE.
  /

ここで計算している系のハミルトニアンは以下の通りです．

.. math::

   \begin{align}
     {\hat H} = \sum_{i}
     \left(
     \begin{matrix}
       {\hat S}_{i x} & {\hat S}_{i y} & {\hat S}_{i z}
     \end{matrix}
     \right)
     \left(
     \begin{matrix}
       1 & 1 & 0 \\
       -1 & 0 & 0 \\
       0 & 0 & 1
     \end{matrix}
     \right)
     \left(
     \begin{matrix}
       {\hat S}_{i+1 x} \ {\hat S}_{i+1 y} \ {\hat S}_{i+1 z}
     \end{matrix}
     \right)\end{align}

Shiftk.outを以下の通り実行します．

.. code-block:: sh

  cd app/sample/denovo
  ShiftK.out namelist.def

上記コマンドにより，ShiftK.outの計算が実施され，以下のように出力されます．

.. code-block:: none

  fname namelist.def
     Open input file namelist.def

     Number of processes :            1


   ##########  Input FileName  ##########

     Hamiltonian file : 
     Excited state file : 

   ##########  Input Parameter for Hamiltonian ##########

       TOTAL Number of sites :            4
       LOCAL Number of sites :            4
                          Jx :    1.0000000000000000     
                          Jy :    0.0000000000000000     
                          Jz :    1.0000000000000000     
                          Dz :    1.0000000000000000     
     Dim. of Hamiltonian :           16

   ##########  Input Parameter for CG Iteration ##########

                     Method : bicg
     Maximum number of loop :          100
      Convergence Threshold :    9.9999999999999995E-007

   ##########  Calculation of the Ground State ##########


     Compute Maximum energy

      Iter      Residual      Energy
         1    0.75200E+00    0.77981E+00
         2    0.42765E+00    0.93541E+00
         3    0.28818E+00    0.11544E+01
         4    0.31071E+00    0.12711E+01
         5    0.74436E-01    0.13687E+01
         6    0.15358E-01    0.14134E+01
         7    0.52412E-02    0.14141E+01
         8    0.35121E-02    0.14142E+01
         9    0.14002E-02    0.14142E+01
        10    0.32273E-03    0.14142E+01
        11    0.11296E-03    0.14142E+01
        12    0.30276E-04    0.14142E+01
        13    0.18218E-04    0.14142E+01
        14    0.65172E-05    0.14142E+01
        15    0.14131E-05    0.14142E+01
        16    0.31249E-06    0.14142E+01
       E_max =    1.4142135623724179     

     Compute Minimum energy

      Iter      Residual      Energy
         1    0.99446E+00    0.71706E+00
         2    0.13481E+01   -0.34893E+00
         3    0.46768E+00   -0.16316E+01
         4    0.23050E+00   -0.19035E+01
         5    0.60898E-01   -0.19769E+01
         6    0.31219E-01   -0.19975E+01
         7    0.20137E-02   -0.19998E+01
         8    0.15690E-03   -0.20000E+01
         9    0.82760E-04   -0.20000E+01
        10    0.33788E-04   -0.20000E+01
        11    0.68419E-05   -0.20000E+01
        12    0.24782E-05   -0.20000E+01
        13    0.24909E-06   -0.20000E+01
       E_min =   -1.9999999999980642     

   ##########  Generate Right Hand Side Vector ##########


   ##########  Input Parameter for Spectrum  ##########

              Max. of Omega :          (1.4142135623724179,3.41421356237048210E-002)
              Min. of Omega :         (-1.9999999999980642,3.41421356237048210E-002)
              Num. of Omega :          100
           Calculation type : normal                                                                                                                                                                                                                                                          

   ##########  CG Initialization  ##########


   #####  BiCG Iteration  #####

      iter status1 status2 status3      Residual       Proj. Res.
         1       1       0      40    0.10330E+02    0.25000E+00
         2       2       0      20    0.43613E+01    0.13281E-14
         3       3       0      59    0.19410E-05    0.20025E-15
         4      -4       0      59    0.88591E-06    0.29334E-14
     Converged in iteration            4
  
   ##########  Output Restart Parameter  ##########

     Num. of Iteration (Current Run) :            4
     Current Omega_Seed :   0.24633E-03  0.34142E-01
  
   ##########  Output Restart Vector  ##########

     Dim. of Residual vector :           16

   #####  Done  #####

上記計算により，以下のファイルが出力されます．

- residual.dat
- output/ResVec.dat0
- output/TriDiagComp.dat
- output/dynamicalG.dat


QMR法によるHeisenberg模型の計算
----------------------------------

ここでは前節のHeisenberg模型でのGreen関数の計算をQMR法にて実行します．
下記のファイルをnamelist.defとして準備します．
CG法との違いとしては、&cg以下にmethod = "shifted_qmr_sym"を入れます．

.. code-block:: none

  &filename
    inham = ""
    invec = ""
  /
  &ham
  Jx = 1d0
  Jy = 0d0
  Jz = 1d0
  Dz = 1d0
  /
  &cg
    method = "shifted_qmr_sym"
    maxloops = 100
    convfactor = 6
  /
  &dyn
    calctype = "normal"
    nomega = 100
  !  omegamin = (-2d0, 0.1d0)
  !  omegamax = ( 1d0, 0.1d0)
  !  omegamin = (2d0, 0d0)
  !  omegamax = (3d0, 0d0)
    outrestart = .TRUE.
  /


ShiftK.outを以下の通り実行します．

.. code-block:: sh

    ShiftK.out namelist.def

上記コマンドにより，ShiftK.outの計算が実施され，以下のように出力されます．

.. code-block:: none

  fname namelist.def
     Open input file namelist.def

     Number of processes :            1


   ##########  Input FileName  ##########

     Hamiltonian file : 
     Excited state file : 

   ##########  Input Parameter for Hamiltonian ##########

       TOTAL Number of sites :            4
       LOCAL Number of sites :            4
                          Jx :    1.0000000000000000     
                          Jy :    0.0000000000000000     
                          Jz :    1.0000000000000000     
                          Dz :    1.0000000000000000     
     Dim. of Hamiltonian :           16

   ##########  Input Parameter for CG Iteration ##########

                     Method : shifted_qmr_sym
     Maximum number of loop :          100
      Convergence Threshold :    9.9999999999999995E-007

   ##########  Calculation of the Ground State ##########


     Compute Maximum energy

      Iter      Residual      Energy
         1    0.79370E+00    0.80884E+00
         2    0.45477E+00    0.99598E+00
         3    0.25991E+00    0.12425E+01
         4    0.22845E+00    0.13313E+01
         5    0.78465E-01    0.13864E+01
         6    0.72022E-02    0.14131E+01
         7    0.33248E-02    0.14142E+01
         8    0.13621E-02    0.14142E+01
         9    0.48926E-03    0.14142E+01
        10    0.10202E-03    0.14142E+01
        11    0.29717E-04    0.14142E+01
        12    0.89912E-05    0.14142E+01
        13    0.93757E-05    0.14142E+01
        14    0.15761E-05    0.14142E+01
        15    0.36535E-06    0.14142E+01
       E_max =    1.4142135623725789     
  
     Compute Minimum energy

      Iter      Residual      Energy
         1    0.81095E+00    0.70090E+00
         2    0.10797E+01   -0.43226E-01
         3    0.60811E+00   -0.10935E+01
         4    0.47092E+00   -0.15142E+01
         5    0.27814E+00   -0.17787E+01
         6    0.19555E+00   -0.19561E+01
         7    0.28018E-01   -0.19915E+01
         8    0.48870E-02   -0.19998E+01
         9    0.32955E-02   -0.20000E+01
        10    0.10375E-02   -0.20000E+01
        11    0.72808E-03   -0.20000E+01
        12    0.10413E-03   -0.20000E+01
        13    0.16623E-04   -0.20000E+01
        14    0.13642E-04   -0.20000E+01
        15    0.35541E-05   -0.20000E+01
        16    0.21036E-05   -0.20000E+01
        17    0.40609E-06   -0.20000E+01
       E_min =   -1.9999999999988347     

   ##########  Generate Right Hand Side Vector ##########


   ##########  Input Parameter for Spectrum  ##########

              Max. of Omega :          (1.4142135623725789,3.41421356237141399E-002)
              Min. of Omega :         (-1.9999999999988347,3.41421356237141399E-002)
              Num. of Omega :          100
           Calculation type : normal                                                                                                                                                                                                                                                          

   ##########  CG Initialization  ##########


   #####  BiCG Iteration  #####

      iter status1 status2 status3      Residual       Proj. Res.
         1       1       0       0    0.77439E-07    0.29318E+00
         2       2       0       0    0.61400E+00    0.44281E-01
         3       3       0       0   -0.11051E+00    0.25475E-01
         4       4       0       0    0.29000E+00    0.17172E-07
         5       5       0       0   -0.16824E+00    0.17172E-07
         6       6       0       0   -0.34423E+00    0.17172E-07
         7      -7       0       0    0.51036E+00    0.17172E-07
     Converged in iteration            7

   ##########  Output Restart Parameter  ##########

     Num. of Iteration (Current Run) :            7
     Current Omega_Seed :   0.00000E+00  0.00000E+00

   ##########  Output Restart Vector  ##########

     Dim. of Residual vector :           16

   #####  Output Restart File  #####

     Restart File is written.

   #####  Done  #####

上記計算により，以下のファイルが出力されます．

- residual.dat
- output/ResVec.dat0
- output/TriDiagComp.dat
- output/dynamicalG.dat
- output/qmr_restart.dat


ハミルトニアンを読み込んでのQMR法の実行
-------------------------------------------

ここではハミルトニアンと状態ベクトルを受け取りQMR法で
Green関数の計算を実施します。
以下のファイルを準備します。

- namelist.def
- zvo_Ham.dat
- zvo_Excited.dat

各ファイルの内容は下記の通りです。

**namelist.def**

.. code-block:: none

  &filename
    inham = "zvo_Ham.dat"
    invec = "zvo_Excited.dat"
    hamtype = "real"
  /
  &ham
  /
  &cg
    method = "shifted_qmr_sym"
    maxloops = 100
    convfactor = 6
  /
  &dyn
    calctype = "normal"
    nomega = 100
  !  omegamin = (-2d0, 0.1d0)
  !  omegamax = ( 1d0, 0.1d0)
  !  omegamin = (2d0, 0d0)
  !  omegamax = (3d0, 0d0)
    outrestart = .TRUE.
  /

**zvo_Ham.dat**

.. code-block:: none

  %%MatrixMarket matrix coordinate complex hermitian
  16 16 20 
  1 1 1.000000 0.000000
  3 2 0.500000 0.000000
  5 3 0.500000 0.000000
  6 4 0.500000 0.000000
  6 6 -1.000000 0.000000
  7 6 0.500000 0.000000
  9 2 0.500000 0.000000
  9 5 0.500000 0.000000
  10 6 0.500000 0.000000
  11 4 0.500000 0.000000
  11 7 0.500000 0.000000
  11 10 0.500000 0.000000
  11 11 -1.000000 0.000000
  12 8 0.500000 0.000000
  13 6 0.500000 0.000000
  13 11 0.500000 0.000000
  14 12 0.500000 0.000000
  15 8 0.500000 0.000000
  15 14 0.500000 0.000000
  16 16 1.000000 0.000000


**zvo_Excited.dat**

.. code-block:: none

  16
        -0.0000000000000000        0.0000000000000000
         0.0000000000000000        0.0000000000000000
        -0.0000000000000000        0.0000000000000000
         0.0000000000000000        0.0000000000000000
         0.0000000000000000       -0.0000000000000000
         0.0000000000000000        0.0000000000000000
         0.2473723858348939        0.1487959544466372
         0.0000000000000000        0.0000000000000000
        -0.0000000000000000        0.0000000000000000
         0.0000000000000000        0.0000000000000000
        -0.4947447716697879       -0.2975919088932743
         0.0000000000000000        0.0000000000000000
         0.2473723858348939        0.1487959544466372
         0.0000000000000000        0.0000000000000000
        -0.0000000000000000       -0.0000000000000000
         0.0000000000000000        0.0000000000000000


Shiftk.outを以下の通り実行します．

.. code-block:: sh

    ShiftK.out namelist.def

上記コマンドにより，ShiftK.outの計算が実施され，以下のように出力されます．

.. code-block:: none

  fname namelist.def
     Open input file namelist.def

     Number of processes :            1


   ##########  Input FileName  ##########

     Hamiltonian file : zvo_Ham.dat
     Excited state file : zvo_Excited.dat

   ##########  Input Parameter for CG Iteration ##########

                     Method : shifted_qmr_sym
     Maximum number of loop :          100
      Convergence Threshold :    9.9999999999999995E-007

   ##########  Input Hamiltonian  ##########

             Dim. of Hamiltonian :           16          16
     Num. of Non-Zero Components :           20

   ##########  Input Right Hand Side Vector ##########

     Dim. of RHS vector :           16

   ##########  Input Parameter for Spectrum  ##########

              Max. of Omega :          (1.0000000000000000,1.00000000000000002E-002)
              Min. of Omega :          (0.0000000000000000,1.00000000000000002E-002)
              Num. of Omega :          100
           Calculation type : normal                                                                                                                                                                                                                                                          

   ##########  CG Initialization  ##########


   #####  BiCG Iteration  #####

      iter status1 status2 status3      Residual       Proj. Res.
         1       1       0       0   -0.00000E+00    0.60594E+00
         2       2       0       0    0.00000E+00    0.29566E+00
         3       3       0       0    0.00000E+00    0.24733E+00
         4      -4       0       0   -0.00000E+00    0.10295E-13
     Converged in iteration            4

   ##########  Output Restart Parameter  ##########

     Num. of Iteration (Current Run) :            4
     Current Omega_Seed :   0.00000E+00  0.00000E+00

   ##########  Output Restart Vector  ##########

     Dim. of Residual vector :           16

   #####  Output Restart File  #####

     Restart File is written.

   #####  Done  #####

上記計算により，以下のファイルが出力されます．

- residual.dat
- timer_inham.dat
- output/ResVec.dat0
- output/TriDiagComp.dat
- output/dynamicalG.dat
- output/qmr_restart.dat

