#! /bin/bash

RPATH=/nfs-6/userdata/stop/output_V00-02-32_2012

declare -a Samples=(ttdl_powheg ttsl_powheg ttfake_powheg w1to4jets data_diel data_dimu ttV diboson triboson data_mueg DY1to4Jtot tW_lepsl tW_lepdl tW_lepfake)
#declare -a Samples=(ttdl_powheg ttsl_powheg w1to4jets data_muo data_ele data_diel data_dimu ttV diboson triboson tWall_lep data_mueg DY1to4Jtot tW_lepsl tW_lepdl)
#declare -a Samples=(T2tt_250_0 T2tt_350_0 T2tt_450_0 T2tt_300_5 T2tt_300_100 ttdl_powheg ttsl_powheg w1to4jets data_muo data_ele data_diel data_dimu ttV diboson triboson tWall data_mueg DYStitchtot)
#declare -a Samples=(T2tt_250_0 T2tt_350_0 T2tt_450_0 T2tt_300_5 T2tt_300_100 ttdl_powheg ttsl_powheg)

for SAMPLE in ${Samples[@]}
  do root -b -q -l do.C\(\"$RPATH\",\"$SAMPLE\"\) &
done
wait
