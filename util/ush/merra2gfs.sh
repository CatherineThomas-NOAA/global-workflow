#!/bin/sh


TMPDIR="/scratch4/NCEPDEV/stmp3/Catherine.Thomas/tmpmerra" #temporary run directory
SRCDIR="/scratch4/NCEPDEV/da/save/Catherine.Thomas/fv3gfs/branches/fv3gfs.master.20181120/util/sorc/merra2gfs.fd" #utility source code directory
CASE="C192"  #output resolution, "C192"
LEVS="128"   #number of levels, "128"
CDATE="2017011500" #date, "YYYYMMDDHH"
INPUT="/scratch4/NCEPDEV/da/noscrub/Catherine.Thomas/GEOS5_2018/MERRA/MERRA2_400.inst3_3d_asm_Nv.20170115.nc4" #input merra file, full path

mkdir -p $TMPDIR
cd $TMPDIR
cp $SRCDIR/* .

ifort -o g5.out -check bounds -traceback Sjl_remap.f90 fv3wam_grids_mod.f90  Read_nc4merra8_qiqlqv.f90 G5_2_FV3.f90 -132 -L/apps/netcdf/4.3.0-intel/lib -lnetcdf -lnetcdff -I/apps/netcdf/4.3.0-intel/include -mcmodel medium -shared-intel

 cat > file_res << EOF
$CASE
$LEVS
$CDATE
$INPUT
EOF

./g5.out

