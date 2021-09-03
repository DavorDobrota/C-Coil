#ifndef GENERAL_COIL_PROGRAM_PRECISIONGLOBALVARS_H
#define GENERAL_COIL_PROGRAM_PRECISIONGLOBALVARS_H

const int g_minPrimLengthIncrements = 6;
const int g_minPrimThicknessIncrements = 6;
const int g_minPrimAngularIncrements = 6;

const int g_minSecLengthIncrements = 4;
const int g_minSecThicknessIncrements = 4;
const int g_minSecAngularIncrements = 4;

const int g_baseLayerIncrements = 10;
const double g_angularWeightModifier = 1;
const double g_primLinearWeightModifier = 1.41;

const double g_thinCoilApproximationRatio = 1e-9;
const double g_zAxisApproximationRatio = 1e-12;

const int g_defaultChunkSize = 2;

#endif //GENERAL_COIL_PROGRAM_PRECISIONGLOBALVARS_H
