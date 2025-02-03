#!/bin/bash 
start=$(date +%s.%N)
end=$(date +%s.%N)
function computeMinutes()
{
dt=$(echo "$end - $start" | bc) 
dd=$(echo "$dt/86400" | bc) 
dt2=$(echo "$dt-86400*$dd" | bc) 
dh=$(echo "$dt2/3600" | bc) 
dt3=$(echo "$dt2-3600*$dh" | bc) 
dm=$(echo "$dt3/60" | bc) 
ds=$(echo "$dt3-60*$dm" | bc) 
}
source /opt/OpenFOAM/OpenFOAM-v2012/etc/bashrc; 
cd /scratch/gopalanh/simulations/SCOUT/snap/ 
cd Case2_NJC_ThermalComfort 
echo "command,status,execution time (in minutes)" | cat > log.caseMonitor
start=$(date +%s.%N)
blockMesh > log.blockMesh
end=$(date +%s.%N)
computeMinutes 
start=$(date +%s.%N)
if tail -n 4 log.blockMesh | grep "End"
then 
 echo "blockMesh,1,$dm" >> log.caseMonitor
 decomposePar -force -copyZero > log.decomposePar 2>&1 
else 
 echo "blockMesh,0,-1" >> log.caseMonitor
exit 1
fi 
end=$(date +%s.%N) 
computeMinutes
start=$(date +%s.%N)
if tail -n 4 log.blockMesh | grep "End"
then 
 echo "blockMesh,1,$dm" >> log.caseMonitor
 surfaceFeatureExtract  > log.surfaceFeatureExtract 2>&1 
else 
 echo "blockMesh,0,-1" >> log.caseMonitor
exit 1
fi 
end=$(date +%s.%N) 
computeMinutes
start=$(date +%s.%N)
if tail -n 4 log.decomposePar | grep "End"
then 
 echo "decomposePar,1,$dm" >> log.caseMonitor
mpirun -np 24 snappyHexMesh -parallel -overwrite > log.snappyHexMesh 2>&1 
else 
 echo "decomposePar,0,-1" >> log.caseMonitor
exit 1
fi 
end=$(date +%s.%N)
computeMinutes
start=$(date +%s.%N)
if tail -n 4 log.snappyHexMesh | grep "Finalising parallel run" 
then 
 echo "snappyHexMesh,1,$dm" >> log.caseMonitor
mpirun -np 24 topoSet -parallel   > log.topoSet 2>&1 
else 
 echo "snappyHexMesh,0,-1" >> log.caseMonitor
exit 1
fi 
end=$(date +%s.%N)
computeMinutes 
start=$(date +%s.%N)
if tail -n 4 log.topoSet | grep "Finalising parallel run" 
then 
 echo "topoSet,1,$dm" >> log.caseMonitor
mpirun -np 24 createPatch -parallel -overwrite > log.createPatch 2>&1 
else 
 echo "topoSet,0,-1" >> log.caseMonitor
exit 1
fi 
end=$(date +%s.%N)
computeMinutes 
mv log.topoSet log.topoSet.1 
start=$(date +%s.%N)
if tail -n 4 log.createPatch | grep "Finalising parallel run" 
then 
 echo "createPatch,1,$dm" >> log.caseMonitor
mpirun -np 24 topoSet -parallel -dict system/topoSetDict1  > log.topoSet 2>&1 
else 
 echo "createPatch,0,-1" >> log.caseMonitor
exit 1
fi 
end=$(date +%s.%N)
computeMinutes 
mv log.createPatch log.createPatch.1 
start=$(date +%s.%N)
if tail -n 4 log.topoSet | grep "Finalising parallel run" 
then 
 echo "topoSet,1,$dm" >> log.caseMonitor
mpirun -np 24 createPatch -parallel -overwrite -dict system/createPatchDict1  > log.createPatch 2>&1 
else 
 echo "topoSet,0,-1" >> log.caseMonitor
exit 1
fi 
end=$(date +%s.%N)
computeMinutes 
mv log.topoSet log.topoSet.2 
start=$(date +%s.%N)
if tail -n 4 log.createPatch | grep "Finalising parallel run" 
then 
 echo "createPatch,1,$dm" >> log.caseMonitor
mpirun -np 24 unsteadyThermalFoam -parallel > log.unsteadyThermalFoam 2>&1 
else 
 echo "createPatch,0,-1" >> log.caseMonitor
exit 1
fi 
mv log.createPatch log.createPatch.2 
end=$(date +%s.%N)
computeMinutes
start=$(date +%s.%N)
if tail -n 4 log.unsteadyThermalFoam | grep "Finalising parallel run" 
then 
 echo "unsteadyThermalFoam,1,$dm" >> log.caseMonitor
mpirun -np 24 postProcess -parallel -func predefinedSampleDict > log.postProcess 2>&1 
else 
 echo "unsteadyThermalFoam,0,-1" >> log.caseMonitor
exit 1
fi 
end=$(date +%s.%N) 
computeMinutes 
if tail -n 4 log.postProcess | grep "Finalising parallel run" 
then 
 echo "postProcess,1,$dm" >> log.caseMonitor
else 
 echo "postProcess,0,-1" >> log.caseMonitor
exit 1
fi 
