/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////             StoR            ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



OFv1806
cd $WM_PROJECT_DIR/src/combustionModels
tree $WM_PROJECT_DIR/src/combustionModels
grep -r EDC



foam
cp -r --parents src/combustionModels $WM_PROJECT_USER_DIR
cd $WM_PROJECT_USER_DIR/src/combustionModels

sed -i s/FOAM_LIBBIN/FOAM_USER_LIBBIN/g Make/files
tail Make/files
wmake

ldd `which sprayFoam` | grep libcombustionModels.so

cp -r EDC StoR
mv StoR/EDC.C StoR/StoR.C
mv StoR/EDC.H StoR/StoR.H
mv StoR/EDCs.C StoR/StoRs.C
sed -Ei 's/EDC/StoR/g' StoR/StoR*.*

vi Make/files

wmakeLnInclude -u .
wmake




///////////////////////////  STOR CASE CREATION  ////////////////////////////////////////////

pwd

cp -r $FOAM_TUTORIALS/lagrangian/sprayFoam/aachenBomb $FOAM_RUN/aachenBombStoR
cp -r $FOAM_TUTORIALS/combustion/reactingFoam/RAS/SandiaD_LTS/constant/combustionProperties $FOAM_RUN/aachenBombStoR/constant


sed -Ei 's/EDC/StoR/g' /home/shijiex/OpenFOAM/shijiex-v1806/run/aachenBombStoR/constant/combustionProperties
vi $FOAM_RUN/aachenBombStoR/constant/combustionProperties


blockMesh -case $FOAM_RUN/aachenBombStoR
sprayFoam -case $FOAM_RUN/aachenBombStoR




/////////////////////////      STOR MODIFICATION    ///////////////////////////////////////////

cp $FOAM_RUN/StoR/StoR* $WM_PROJECT_DIR/src/combustionModels/StoR

vi StoR/StoRs.C

vi StoR/StoR.H


wmake
blockMesh -case $FOAM_RUN/aachenBombStoR
sprayFoam -case $FOAM_RUN/aachenBombStoR
more log.sprayFoam





















/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////              CCM            ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

OFv1806
foam
cp -r --parents src/thermophysicalModels/chemistryModel $WM_PROJECT_USER_DIR
cd $WM_PROJECT_USER_DIR/src/thermophysicalModels/chemistryModel

sed -i s/FOAM_LIBBIN/FOAM_USER_LIBBIN/g Make/files
tail Make/files
wmake

ldd `which sprayFoam` | grep libchemistryModel.so

cp -r chemistryModel/StandardChemistryModel chemistryModel/CCMChemistryModel
mv chemistryModel/CCMChemistryModel/StandardChemistryModel.C chemistryModel/CCMChemistryModel/CCMChemistryModel.C
mv chemistryModel/CCMChemistryModel/StandardChemistryModel.H chemistryModel/CCMChemistryModel/CCMChemistryModel.H
mv chemistryModel/CCMChemistryModel/StandardChemistryModelI.H chemistryModel/CCMChemistryModel/CCMChemistryModelI.H
sed -Ei 's/StandardChemistryModel/CCMChemistryModel/g' chemistryModel/CCMChemistryModel/CCMChemistryModel*.*
sed -Ei 's/standard/CCM/g' chemistryModel/CCMChemistryModel/CCMChemistryModel*.*


grep -r TDAC

cp $FOAM_RUN/CCM/BasicChemistryModels.C chemistryModel/BasicChemistryModel
cp $FOAM_RUN/CCM/basicChemistryModelTemplates.C chemistryModel/basicChemistryModel
cp $FOAM_RUN/CCM/makeChemistrySolverTypes.H chemistrySolver/chemistrySolver

cp $FOAM_RUN/CCM/CCMChemistryModel.C chemistryModel/CCMChemistryModel
cp $FOAM_RUN/CCM/CCMChemistryModel.H chemistryModel/CCMChemistryModel
cp $FOAM_RUN/CCM/CCMChemistryModelI.H chemistryModel/CCMChemistryModel

touch chemistryModel/makeChemistryModel.H
wmakeLnInclude -u .
wmake



///////////////////////////  CCM CASE CREATION  ////////////////////////////////////////////

cp -r $FOAM_TUTORIALS/lagrangian/sprayFoam/aachenBomb $FOAM_RUN/aachenBombCCM
cp -r $FOAM_TUTORIALS/combustion/reactingFoam/laminar/counterFlowFlame2D/constant/combustionProperties $FOAM_RUN/aachenBombCCM/constant

vi $FOAM_RUN/aachenBombCCM/constant/chemistryProperties
cp $FOAM_RUN/CCM/CCMProperties $FOAM_RUN/aachenBombCCM/constant

blockMesh -case $FOAM_RUN/aachenBombCCM
sprayFoam -case $FOAM_RUN/aachenBombCCM



















/////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////             StoR with CCM CASE            ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

foam
cp -r $FOAM_TUTORIALS/lagrangian/sprayFoam/aachenBomb $FOAM_RUN/aachenBombCCMStoR
cd $FOAM_RUN/aachenBombCCMStoR

vi constant/combustionProperties

vi constant/chemistryProperties

cp $FOAM_RUN/CCM/CCMProperties $FOAM_RUN/aachenBombCCMStoR/constant
vi constant/CCMProperties

touch Allrun
vi Allrun
/*******************************************************************************************/
#!/bin/sh
cd ${0%/*} || exit 1				# Run from this directory
. $WM_PROJECT_DIR/bin/tools/RunFunctions	# Tutorial run functions

runApplication blockMesh
runApplication $(getApplication)

#---------------------------------------------------------------------
/*******************************************************************************************/


touch Allclean
vi Allclean
/*******************************************************************************************/
#!/bin/sh
cd ${0%/*} || exit 1				# Run from this directory
. $WM_PROJECT_DIR/bin/tools/CleanFunctions	# Tutorial clean functions

cleanCase

#---------------------------------------------------------------------
/*******************************************************************************************/


chmod 755 All*
ll All*

./Allrun
more $FOAM_RUN/aachenBombCCMStoR/log.sprayFoam
tail $FOAM_RUN/aachenBombCCMStoR/log.sprayFoam


cp -r $FOAM_TUTORIALS/lagrangian/sprayFoam/aachenBomb $FOAM_RUN/aachenBomb
cd $FOAM_RUN/aachenBomb
cp $FOAM_RUN/aachenBombCCMStoR/All* $FOAM_RUN/aachenBomb
./Allclean
./Allrun
