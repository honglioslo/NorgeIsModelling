#pragma once

class HbvAquifer;
class KwaAquifer;
class Glacier;
class Lake;

class Albedo
{
public:
    Albedo();
    ~Albedo();
    void SurfaceAlbedo(Lake *ptrLake, HbvAquifer *ptrHbvAquifer, KwaAquifer *ptrKwaAquifer, Glacier *ptrGlacier);
};
