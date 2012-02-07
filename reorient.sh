#!/bin/sh -x
#From Mike Hodge via Mark Jenkinson
export FSLOUTPUTTYPE=NIFTI
input=$1
newim=$2
im=${input}_std
fslreorient2std $input $im 
dx=`$FSLDIR/bin/fslval $im pixdim1`;
dy=`$FSLDIR/bin/fslval $im pixdim2`;
dz=`$FSLDIR/bin/fslval $im pixdim3`;
ori=`$FSLDIR/bin/fslorient $im`;

if [ $ori = RADIOLOGICAL ] ; then
  dx=`echo "$dx * -1" | bc -l`;
fi

$FSLDIR/bin/fslhd $im | grep sto_xyz | sed 's/sto_xyz:.[    ]*//' > ${newim}_oldsform.mat

cp -p ${im}.nii ${newim}.nii

$FSLDIR/bin/convert_xfm -omat ${newim}_oldsforminv.mat -inverse ${newim}_oldsform.mat
origx=`sed -n 1p ${newim}_oldsforminv.mat | awk '{ print $4 }'`
origy=`sed -n 2p ${newim}_oldsforminv.mat | awk '{ print $4 }'`
origz=`sed -n 3p ${newim}_oldsforminv.mat | awk '{ print $4 }'`
origx=`echo "-1 * $origx * $dx" | bc -l`
origy=`echo "-1 * $origy * $dy" | bc -l`
origz=`echo "-1 * $origz * $dz" | bc -l`
$FSLDIR/bin/fslorient -setsform $dx 0 0 $origx 0 $dy 0 $origy 0 0 $dz $origz 0 0 0 1 $newim 
$FSLDIR/bin/fslorient -copysform2qform $newim

$FSLDIR/bin/fslhd $newim | grep sto_xyz | sed 's/sto_xyz:.[         ]*//' > ${newim}_newsform.mat

$FSLDIR/bin/fslhd -x ${newim}.nii > ${newim}.xml