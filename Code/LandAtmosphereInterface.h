#pragma once
class LandAtmosphereInterface
{
public:
    void SetLandAtmosphereInterfaceType(char value);
    void SetEvaporationModelling(char value);
    char GetEvaporationModelling() const;
    LandAtmosphereInterface();
    ~LandAtmosphereInterface();

private:
    char evaporationModelling;
};
