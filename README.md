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
We will use dunesw version v09_49_00d00 -q e20:prof.
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
