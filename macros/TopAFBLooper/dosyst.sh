#! /bin/bash
RPATHttbar=/nfs-6/userdata/linacre/babies

declare -a ttbarSamples=(ttdl_FullLept_mass169_5_8TeV_mcatnlo ttdl_FullLept_scaleup_8TeV_mcatnlo ttdl_FullLept_scaledown_8TeV_mcatnlo ttdl_FullLept_mass175_5_8TeV_mcatnlo ttdl_noCorr_8TeV_mcatnlo ttsl_noCorr_8TeV_mcatnlo ttfake_noCorr_8TeV_mcatnlo ttdl_powheg_SCbug ttsl_powheg_SCbug ttfake_powheg_SCbug ttdl_powheg_SCfix ttsl_powheg_SCfix ttfake_powheg_SCfix)

for SAMPLE in ${ttbarSamples[@]}
  do root -b -q -l do.C\(\"$RPATHttbar\",\"$SAMPLE\"\) &
done

wait
