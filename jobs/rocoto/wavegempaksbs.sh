#!/bin/ksh -x

###############################################################
echo
echo "=============== START TO SOURCE FV3GFS WORKFLOW MODULES ==============="
. $HOMEgfs/ush/load_fv3gfs_modules.sh
status=$?
[[ $status -ne 0 ]] && exit $status


###############################################################
echo
echo "=============== START TO SOURCE RELEVANT CONFIGS ==============="
configs="base wave wavegempaksbs"
for config in $configs; do
    . $EXPDIR/config.${config}
    status=$?
    [[ $status -ne 0 ]] && exit $status
done


###############################################################
echo
echo "=============== START TO RUN WAVE GEMPAK_SBS ==============="
# Execute the JJOB
$HOMEgfs/jobs/JWAVE_GEMPAK_SBS
status=$?
exit $status

###############################################################
# Force Exit out cleanly
if [ ${KEEPDATA:-"NO"} = "NO" ] ; then rm -rf $DATAROOT ; fi
exit 0
