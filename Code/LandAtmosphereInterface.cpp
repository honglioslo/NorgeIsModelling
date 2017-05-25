#include "LandAtmosphereInterface.h"

LandAtmosphereInterface::LandAtmosphereInterface() :
    evaporationModelling('E')
{
}

LandAtmosphereInterface::~LandAtmosphereInterface()
{
}

void  LandAtmosphereInterface::SetEvaporationModelling(char value)
{
    evaporationModelling = value;
}
char  LandAtmosphereInterface::GetEvaporationModelling() const
{
    return evaporationModelling;
}


void LandAtmosphereInterface::SetLandAtmosphereInterfaceType(char value)
{
    SetEvaporationModelling(value);
}
