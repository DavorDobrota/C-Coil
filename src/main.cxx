#include <cstdio>
#include <cmath>

#include "Coil.h"
#include "CoilGroup.h"
#include "Test.h"


int main()
{
//    testCoilGroupMTvsMTD(32, 50'000);

//    double temperatureK = 1000.0;
//    double voltageU = 5.0;
//    double initialZ = 0.4;
//
//    double mass = 9.10938e-31;
//    double charge = 1.60218e-19;
//    double boltzmann = 1.38065e-23;
//
//    double timeStep = 1e-12;
//
//    double radialVelocity = std::sqrt(2 * temperatureK * boltzmann / mass);
//    double zAxisVelocity = std::sqrt((boltzmann * temperatureK + 2 * charge * voltageU) / mass);
//    vec3::CoordVector3 position = vec3::CoordVector3(vec3::CYLINDRICAL, initialZ, 0, 0);
//
//    Coil fieldCoil = Coil(0.05, 0.05, 0.05, 2500);
//
//    vec3::FieldVector3 magneticField = fieldCoil.computeBFieldVector(position);
//
//    double magneticForceZ = 0.0;
//
//    int cnt = 0;
//
//    while (position.component1 > 0.0)
//    {
//        position.component2 = mass * radialVelocity / (charge * magneticField.zComponent);
//
//        magneticField = fieldCoil.computeBFieldVector(position);
//
//        magneticForceZ = charge * radialVelocity * magneticField.xComponent;
//        zAxisVelocity -= (magneticForceZ * timeStep) / mass;
//
//        position.component1 -= (zAxisVelocity * timeStep + 0.5 * magneticForceZ * timeStep * timeStep);
//        cnt++;
//        if (cnt % 100 == 0)
//            printf("%d : %.7g %.7g  %.7g\n", cnt, position.component1, position.component2, zAxisVelocity);
//    }

//    double temperatureK = 1000.0;
//    double voltageU = 0.0;
//    double initialZ = 0.000000000001;
//
//    double mass = 9.10938e-31;
//    double charge = 1.60218e-19;
//    double boltzmann = 1.38065e-23;
//
//    double timeStep = 1e-16;
//
//    double xAxisVelocity = std::sqrt(temperatureK * boltzmann / mass);
//    double yAxisVelocity = std::sqrt(temperatureK * boltzmann / mass);
//    double zAxisVelocity = 1.0 * std::sqrt((temperatureK * boltzmann + 2 * charge * voltageU) / mass);
//
//    vec3::FieldVector3 positionVec = vec3::FieldVector3(0.0, 0.0, initialZ);
//    vec3::FieldVector3 velocityVec = vec3::FieldVector3(xAxisVelocity, yAxisVelocity, zAxisVelocity);
//
//    Coil fieldCoil = Coil(0.05, 0.05, 0.05, 2500,4, PrecisionFactor(1.0));
//    vec3::FieldVector3 magneticField;
//    vec3::FieldVector3 magneticForce;
//
//    int cnt = 0;
//    double simSpaceLimit = 1.0;
//
//    while (positionVec.zComponent > 0.0 && positionVec.zComponent <= simSpaceLimit)
//    {
//        if (cnt % 1000000 == 0)
//            printf("%d : %.7g %.7g %.7g |  %.7g %.7g %.7g | %.7g\n", cnt,
//                   positionVec.xComponent, positionVec.yComponent, positionVec.zComponent,
//                   velocityVec.xComponent, velocityVec.yComponent, velocityVec.zComponent,
//                   std::sqrt(velocityVec.xComponent * velocityVec.xComponent +
//                             velocityVec.yComponent * velocityVec.yComponent +
//                             velocityVec.zComponent * velocityVec.zComponent));
//
//        magneticField = fieldCoil.computeBFieldVector(vec3::CoordVector3::convertToCoordVector(positionVec));
//
//        magneticForce = vec3::FieldVector3::crossProduct(velocityVec, magneticField) * charge;
//        velocityVec += magneticForce * (timeStep / mass);
//        positionVec += velocityVec * timeStep + magneticForce * (0.5 * timeStep * timeStep / mass);
//
//        cnt++;
//    }

//    double wireThickness = 0.0003;
//
//    for (int i = 1; i < 200; ++i)
//    {
//        Coil coil = Coil(0.0125,  wireThickness * 4, i * wireThickness, i * 4);
//        printf("%.15g\n", coil.computeAndSetSelfInductance(PrecisionFactor(8.0), CPU_MT));
//    }
//
//    Coil coil = Coil(0.01, 0.009, 0.05, 450);
//    printf("%.10g\n", coil.computeAndSetSelfInductance(PrecisionFactor(10.0), CPU_MT));
//
//    coil = Coil(0.01, 0.009, 0.05, 480);
//    printf("%.10g\n", coil.computeAndSetSelfInductance(PrecisionFactor(10.0), CPU_MT));


//    Coil coil = Coil(0.1, 0.1, 0.1, 10000);
//    printf("%.15g\n", coil.computeAndSetSelfInductance(PrecisionFactor(10.0), CPU_MT));

    testCoilGroupMTDInductanceAndForce(12);

    return 0;
}
