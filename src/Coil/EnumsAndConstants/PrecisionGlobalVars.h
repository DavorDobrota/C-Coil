#ifndef GENERAL_COIL_PROGRAM_PRECISIONGLOBALVARS_H
#define GENERAL_COIL_PROGRAM_PRECISIONGLOBALVARS_H

const int g_minPrimLengthIncrements = 2;
const int g_minPrimThicknessIncrements = 2;
const int g_minPrimAngularIncrements = 2;

const int g_minSecLengthIncrements = 2;
const int g_minSecThicknessIncrements = 2;
const int g_minSecAngularIncrements = 2;

const int g_baseLayerIncrementsCPU = 10;
const int g_baseLayerIncrementsGPU = 10;

const double g_primAngularWeightModifier = 1.414213562;
const double g_primLinearWeightModifier = 1.414213562;

const double g_secAngularWeightModifier = 1;
const double g_secLinearWeightModifier = 1;

const double g_thinCoilApproximationRatio = 1e-7;
const double g_zAxisApproximationRatio = 1e-12;
const double g_commonPlaneApproximationRatio = 1e-14;

#endif //GENERAL_COIL_PROGRAM_PRECISIONGLOBALVARS_H
