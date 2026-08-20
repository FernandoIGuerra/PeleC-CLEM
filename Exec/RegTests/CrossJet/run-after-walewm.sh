#!/bin/bash
# Wait for the running WALE+WM job (mpiexec pid 4176900) to finish, then run the
# non-reacting Smagorinsky case on 64 ranks. Detached with setsid/nohup so it
# survives the terminal and the agent session.
#
# Output layout (plt/chk paths come from the .inp files themselves):
#   ResultsNonReact/CLEM-WALE-WM/      plt_nrWM*, chk_nrWM*, run log, datlog
#   ResultsNonReact/CLEM-Smagorinsky/  plt_nrS*,  chk_nrS*,  run log, datlog
#
# NOT chained: the WALE-without-wall-model case. Run 04 already provides it
# (WALE, no WM, chemically inert for all 0.94 ms). To add it anyway, append a
# third block below using inputs-case1-nonreacting-CLEM-WALE.inp.
cd /home/pc08/Fernando/CLEMNew/PeleC-CLEM/Exec/RegTests/CrossJet || exit 1

WM=ResultsNonReact/CLEM-WALE-WM
NS=ResultsNonReact/CLEM-Smagorinsky
mkdir -p "$WM" "$NS"

echo "[chain] $(date '+%F %T')  waiting on WALE+WM (pid 4176900)"
while kill -0 4176900 2>/dev/null; do sleep 60; done
echo "[chain] $(date '+%F %T')  WALE+WM exited"
sleep 30

# Collect the finished WALE+WM run's log and datlog next to its plotfiles.
mv -f nonReactWaleWM.log "$WM"/ 2>/dev/null
mv -f datlog_nrWM        "$WM"/ 2>/dev/null
echo "[chain] $(date '+%F %T')  WALE+WM log and datlog moved into $WM"

echo "[chain] $(date '+%F %T')  starting Smagorinsky on 64 ranks"
mpiexec -n 64 PeleC3d.gnu.MPI.ex \
    inputs-case1-nonreacting-CLEM-smagorinsky.inp \
    > "$NS"/nonReactSmagorinsky.log 2>&1
st=$?
mv -f datlog_nrS "$NS"/ 2>/dev/null
echo "[chain] $(date '+%F %T')  Smagorinsky exited with status $st"
