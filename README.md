# Instruction to run ProtoDUNE-SP Calibration Codes
## Login to dunegpvm
I will login to dunegpvm08 for this instruction.
```
$ ssh -Y <user_name>@dunegpvm08.fnal.gov
```

## Setup an area
Then, move to your app directory.
```
$ cd /dune/app/users/<user_name>
```
Make a directory for protodune calibration codes.
```
$ mkdir ProtoDUNE
$ cd ProtoDUNE
$ mkdir calib
$ cd calib
```
We will use dunesw version `v09_49_00d00 -q e20:prof`.
First,
```
$ source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
```
Then, make an dunesw area
```
$ mrb newDev -v v09_49_00d00 -q e20:prof
```
You will see messages like bellow.
```
building development area for larsoft v09_49_00d00 -q e20:prof


The following configuration is defined:
  The top level directory is .
  The source code directory will be under .
  The build directory will be under .
  The local product directory will be under .

MRB_BUILDDIR is /dune/app/users/sungbino/ProtoDUNE/test_for_Abbey/build_slf7.x86_64
MRB_SOURCE is /dune/app/users/sungbino/ProtoDUNE/test_for_Abbey/srcs 
INFO: cannot find larsoft/v09_49_00d00/releaseDB/base_dependency_database
      or larsoftcode/v09_49_00d00/releaseDB/base_dependency_database
      mrb checkDeps and pullDeps will not have complete information

IMPORTANT: You must type
    source /dune/app/users/sungbino/ProtoDUNE/test_for_Abbey/localProducts_larsoft_v09_49_00d00_e20_prof/setup
NOW and whenever you log in
```
Then,
```
$ source localProducts_larsoft_v09_49_00d00_e20_prof/setup 
```
You will see messages like
```
MRB_PROJECT=larsoft
MRB_PROJECT_VERSION=v09_49_00d00
MRB_QUALS=e20:prof
MRB_TOP=/dune/app/users/sungbino/ProtoDUNE/test_for_Abbey
MRB_SOURCE=/dune/app/users/sungbino/ProtoDUNE/test_for_Abbey/srcs
MRB_BUILDDIR=/dune/app/users/sungbino/ProtoDUNE/test_for_Abbey/build_slf7.x86_64
MRB_INSTALL=/dune/app/users/sungbino/ProtoDUNE/test_for_Abbey/localProducts_larsoft_v09_49_00d00_e20_prof

PRODUCTS=/dune/app/users/sungbino/ProtoDUNE/test_for_Abbey/localProducts_larsoft_v09_49_00d00_e20_prof:/cvmfs/dune.opensciencegrid.org/products/dune:/cvmfs/larsoft.opensciencegrid.org/products:/cvmfs/larsoft.opensciencegrid.org/packages:/cvmfs/fermilab.opensciencegrid.org/products/common/db/
CETPKG_INSTALL=/dune/app/users/sungbino/ProtoDUNE/test_for_Abbey/localProducts_larsoft_v09_49_00d00_e20_prof
```
Let's copy "protoduneana"
```
$ mrb g protoduneana
```
You will see messages like
```
Cloning into 'protoduneana'...
X11 forwarding request failed on channel 0
remote: Enumerating objects: 8516, done.
remote: Counting objects: 100% (597/597), done.
remote: Compressing objects: 100% (208/208), done.
remote: Total 8516 (delta 434), reused 508 (delta 389), pack-reused 7919
Receiving objects: 100% (8516/8516), 2.99 MiB | 16.45 MiB/s, done.
Resolving deltas: 100% (6403/6403), done.
Updating files: 100% (487/487), done.
NOTICE: Adding protoduneana to CMakeLists.txt file
```
Let's synchronize protoduneana's version with the dunssw
```
$ cd src/protoduneana
$ git tag
```
You will see a list of tags.
```
$ git checkout v09_49_00d00
```
You will see message like
```
Note: switching to 'v09_49_00d00'.

You are in 'detached HEAD' state. You can look around, make experimental
changes and commit them, and you can discard any commits you make in this
state without impacting any branches by switching back to a branch.

If you want to create a new branch to retain commits you create, you may
do so (now or later) by using -c with the switch command. Example:

  git switch -c <new-branch-name>

Or undo this operation with:

  git switch -

Turn off this advice by setting config variable advice.detachedHead to false

HEAD is now at 831b48f Merge branch 'release/v09_49_00d00'
```
Then,
```
$ git checkout -b v09_49_00d00
```
If you check the branc, you should see
```
$ git branch -vv
  develop      d97577d [origin/develop] Merge tag 'v09_75_00d00' into develop
* v09_49_00d00 831b48f Merge branch 'release/v09_49_00d00'
```
Let's build the area
```
$ cd ../
$ mrbsetenv
```
You will see message like
```
The working build directory is /dune/app/users/sungbino/ProtoDUNE/test_for_Abbey/build_slf7.x86_64
The source code directory is /dune/app/users/sungbino/ProtoDUNE/test_for_Abbey/srcs
----------- check this block for errors -----------------------
INFO: mrb v6_08_00 requires cetmodules >= 2.31.00 to run: attempting to configure...v3_21_01 OK
----------------------------------------------------------------
To inspect build variable settings, execute /dune/app/users/sungbino/ProtoDUNE/test_for_Abbey/build_slf7.x86_64/cetpkg_info.sh

Please use "buildtool" (or "mrb b") to configure and build MRB project "larsoft", e.g.:

  buildtool -vTl [-jN]

See "buildtool --usage" (short usage help) or "buildtool -h|--help"
(full help) for more details.
```
Then, you can build the area,
```
$ mrb i -j4
```
You should see the following message after building the area.
```
------------------------------------
INFO: stage install SUCCESS for MRB project larsoft v09_49_00d00
------------------------------------
```

## Run
Please note that you can use these commands to run codes in the area for next login's
```
$ cd /dune/app/users/<user_name>/ProtoDUNE/calib
$ source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
$ setup dunesw v09_49_00d00 -q e20:prof
$ source /dune/app/users/sungbino/ProtoDUNE/larsoft/localProducts_larsoft_v09_49_00d00_e20_prof/setup
$ mrbslp
$ mrbsetenv
```

Let's assume that you are on `/dune/app/users/<user_name>/ProtoDUNE/calib`. Go to the directory where we have codes for calibration
```
$ cd srcs/protoduneana/protoduneana/singlephase/michelremoving/
```
Make a txt file named as `input_run5387.txt` with contents.
```
root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2020/detector/physics/PDSPProd4/00/00/53/87/np04_run005387_PDSPProd4_michelremoving_merged_0_200.root
root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2020/detector/physics/PDSPProd4/00/00/53/87/np04_run005387_PDSPProd4_michelremoving_merged_200_400.root
root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2020/detector/physics/PDSPProd4/00/00/53/87/np04_run005387_PDSPProd4_michelremoving_merged_400_600.root
root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2020/detector/physics/PDSPProd4/00/00/53/87/np04_run005387_PDSPProd4_michelremoving_merged_600_800.root
root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/protodune-sp/root-tuple/2020/detector/physics/PDSPProd4/00/00/53/87/np04_run005387_PDSPProd4_michelremoving_merged_800_950.root
```
Note that the ProtoDUNE-SP's Run5387 is used as general reference for Run1 calibration.
Get a valid access to those root files using commands bellow (I recommand make a .sh file in your home direcoty with those commands for your convenience).
```
$ kx509
$ export EXPERIMENT=dune
$ export ROLE=Analysis
$ voms-proxy-init -rfc -noregen -voms dune:/dune/Role=$ROLE -valid 24:00
```

### Uniformity correction : YZ correction
Run the code at `/dune/app/users/<user_name>/ProtoDUNE/calib/srcs/protoduneana/protoduneana/singlephase/michelremoving/`.
```
$ make_yz_correction input_run5387.txt 2
```
It will take several minutes with messages like
```
2
michelremoving2/Event
Plugin version SecClnt v5.1.0 is incompatible with secztn v5.5.3 (must be <= 5.1.x) in sec.protocol libXrdSecztn-5.so
Process Run 5387
0/122143
10000/122143
20000/122143
30000/122143
40000/122143
50000/122143
60000/122143
70000/122143
80000/122143
90000/122143
100000/122143
110000/122143
120000/122143
*************** Calculating the local median dQ/dx values for each Y-Z cell ******************
********************** Calculating fractional dQ/dx corrections for each Y-Z cell ********************
******************** Calculating corrected dQ/dx value for each Y-Z cell **********************
********************** Calculating fractional dQ/dx corrections for each Y-Z cell ********************
******************** Calculating corrected dQ/dx value for each Y-Z cell **********************
********************** Calculating fractional dQ/dx corrections for each Y-Z cell ********************
******************** Calculating corrected dQ/dx value for each Y-Z cell **********************
Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
crossing tracks 147931
*************** Y_Z_Correction_make_class.C macro has ended ******************
```
Then, you will have a new file (`YZcalo_mich2_r5387.root`)in the directory.
```
-bash-4.2$ ls -lht
total 384K
-rw-r--r-- 1 sungbino dune 1009K May 30 10:39 YZcalo_mich2_r5387.root
-rw-r--r-- 1 sungbino dune   993 May 30 10:31 input_run5387.txt
.
.
.
```

### Uniformity correction : X correction
Before running, we should modify/check codes for several values.
Open `./Xcalo/protoDUNE_X_calib.C`, and modify
- Line 29 : modify/checkthe liquid argon density value
- Line 30 - 31 : modify/check the modified box model parameters
Then, build the area again,
```
$ mrb i -j4
```
Let's run at `/dune/app/users/<user_name>/ProtoDUNE/calib/srcs/protoduneana/protoduneana/singlephase/michelremoving/`,
```
$ make_x_correction input_run5387.txt 2 1
```
It will take several minute with messages bellow.
```
michelremoving2/Event
SCE on
efield at the anode neg0.432389
efield at the anode pos0.449398
Plugin version SecClnt v5.1.0 is incompatible with secztn v5.5.3 (must be <= 5.1.x) in sec.protocol libXrdSecztn-5.so
0/122143
10000/122143
20000/122143
30000/122143
40000/122143
50000/122143
60000/122143
70000/122143
80000/122143
90000/122143
100000/122143
110000/122143
120000/122143
*************** Calculating the local median dQ/dx values for each Y-Z cell ******************
**************** Calculating fractional correction for each x cell *********************
**************** Calculating XYZ corrected dQ/dx values ********************
**************** Calculating fractional correction for each x cell *********************
**************** Calculating XYZ corrected dQ/dx values ********************
**************** Calculating fractional correction for each x cell *********************
**************** Calculating XYZ corrected dQ/dx values ********************
Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
*************** X_Correction_make_class.C macro has ended ******************
```
with several new files,
```
-bash-4.2$ ls -lht
total 1.5M
-rw-r--r-- 1 sungbino dune    13 May 30 10:56 global_median_0_r5387.txt
-rw-r--r-- 1 sungbino dune    13 May 30 10:56 global_median_1_r5387.txt
-rw-r--r-- 1 sungbino dune    13 May 30 10:56 global_median_2_r5387.txt
-rw-r--r-- 1 sungbino dune  5.9K May 30 10:56 globalmedians_cathanode_r5387.root
-rw-r--r-- 1 sungbino dune   17K May 30 10:56 Xcalo_mich2_r5387.root
.
.
.
```
