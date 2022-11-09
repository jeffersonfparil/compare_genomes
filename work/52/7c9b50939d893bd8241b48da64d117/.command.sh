#!/usr/bin/env bash
cd /data-weedomics-3
wget 'http://data.pantherdb.org/ftp/panther_library/current_release/PANTHER17.0_hmmscoring.tgz'
tar -xvzf PANTHER17.0_hmmscoring.tgz
rm PANTHER17.0_hmmscoring.tgz
mv target/ PantherHMM_17.0/
cd PantherHMM_17.0/
wget 'http://data.pantherdb.org/ftp/hmm_classifications/current_release/PANTHER17.0_HMM_classifications'
grep -v ':SF' PANTHER17.0_HMM_classifications > Panther17.0_HMM_familyIDs.txt
