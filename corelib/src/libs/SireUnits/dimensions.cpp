
#include "dimensions.h"
#include "temperature.h"

namespace SireUnits
{

namespace Dimension
{

template class PhysUnit<0,0,0,0,0,0,0>;
template class PhysUnit<1,0,0,0,0,0,0>;
template class PhysUnit<1,0,0,0,0,-1,0>;
template class PhysUnit<0,1,0,0,0,0,0>;
template class PhysUnit<0,0,1,0,0,0,0>;// Time;
template class PhysUnit<0,0,0,1,0,0,0>;// Charge;
template class PhysUnit<0,0,0,1,0,-1,0>;// MolarCharge;
template class PhysUnit<0,0,0,0,1,0,0>;// Temperature;
template class PhysUnit<0,0,0,0,0,1,0>;// Quantity;
template class PhysUnit<0,0,0,0,0,0,1>;// Angle;
template class PhysUnit<0,2,0,0,0,0,0>;// Area;
template class PhysUnit<0,3,0,0,0,0,0>;// Volume;
template class PhysUnit<0,3,0,0,0,-1,0>;// MolarVolume;
template class PhysUnit<0,1,-1,0,0,0,0>;// Velocity;
template class PhysUnit<0,0,-1,0,0,0,1>;// AngularVelocity;
template class PhysUnit<0,1,-2,0,0,0,0>;// Acceleration;
template class PhysUnit<0,0,-2,0,0,0,1>;// AngularAcceleration;
template class PhysUnit<1,2,-2,0,0,0,0>;// Energy;
template class PhysUnit<1,2,-2,0,0,-1,0>;// MolarEnergy;
template class PhysUnit<1,2,-3,0,0,0,0>;// Power;
template class PhysUnit<1,2,-3,0,0,-1,0>;// MolarPower;
template class PhysUnit<1,-3,0,0,0,0,0>;// Density;
template class PhysUnit<1,-3,0,0,0,-1,0>;// MolarDensity;
template class PhysUnit<1,1,-2,0,0,0,0>;// Force;
template class PhysUnit<1,-1,-2,0,0,0,0>;// Pressure;
template class PhysUnit<0,0,-1,1,0,0,0>;// Current;
template class PhysUnit<-1,-2,2,2,0,0,0>;// Capacitance;
template class PhysUnit<1,2,-2,-1,0,0,0>;// Potential;

}

}
