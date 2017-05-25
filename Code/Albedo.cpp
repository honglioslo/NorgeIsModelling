#include "Albedo.h"
#include "Glacier.h"
#include "GlacierIce.h"
#include "HbvAquifer.h"
#include "HBV.h"
#include "InputElement.h"
#include "Lake.h"
#include "HbvAquifer.h"
#include "KwaAquifer.h"
#include "KinematicWave.h"
#include "LakeWaterBalance.h"
#include "ParametersLandSurface.h"
#include "Vegetation.h"
#include "Snow.h"

#include "stdafx.h"

using namespace std;

Albedo::Albedo()
{
}


Albedo::~Albedo()
{
}


void Albedo::SurfaceAlbedo(Lake *ptrLake, HbvAquifer *ptrHbvAquifer, KwaAquifer *ptrKwaAquifer, Glacier *ptrGlacier)
{
    int i;
    double precipitation, temperature, accumulationTemperature, leafAreaIndex, interceptionStore, snowCoverFraction;
    double albedoSnow, albedoVegetation, albedoSoil, albedoSurface, albedoReduction;
    double minSnowAlbedo = 0.2, maxSnowAlbedo = 0.8;
    LANDSURFACE landSurfaceType;

    // Lake albedo
    if (ptrLake)
    {
        albedoSurface = ptrLake->GetLakeWaterBalance()->GetLandSurfacePar()->GetALBEDO();
        ptrLake->SetSurfaceAlbedo(albedoSurface);
    }
    else
    {

        // Read values from previous time step
        if (ptrHbvAquifer)
        {
            landSurfaceType = ptrHbvAquifer->GetLandSurfaceType();
            albedoVegetation = ptrHbvAquifer->GetHBV()->GetLandSurfacePar()->GetALBEDO();
            precipitation = ptrHbvAquifer->GetHBV()->GetInputElement()->GetInput(0);
            temperature = ptrHbvAquifer->GetHBV()->GetInputElement()->GetInput(1);
            accumulationTemperature = ptrHbvAquifer->GetHBV()->GetLandSurfacePar()->GetACC_TEMP();
            leafAreaIndex = ptrHbvAquifer->GetVegetation()->GetLeafAreaIndex();
            interceptionStore = ptrHbvAquifer->GetVegetation()->GetInterceptionStore();
            albedoSnow = ptrHbvAquifer->GetSnow()->GetSnowAlbedo();
            snowCoverFraction = ptrHbvAquifer->GetSnow()->GetSnowCoverFraction();
        }
        else if (ptrKwaAquifer)
        {
            landSurfaceType = ptrKwaAquifer->GetLandSurfaceType();
            albedoVegetation = ptrKwaAquifer->GetKinematicWave()->GetLandSurfacePar()->GetALBEDO();
            precipitation = ptrKwaAquifer->GetKinematicWave()->GetInputElement()->GetInput(0);
            temperature = ptrKwaAquifer->GetKinematicWave()->GetInputElement()->GetInput(1);
            accumulationTemperature = ptrKwaAquifer->GetKinematicWave()->GetLandSurfacePar()->GetACC_TEMP();
            leafAreaIndex = ptrKwaAquifer->GetVegetation()->GetLeafAreaIndex();
            interceptionStore = ptrKwaAquifer->GetVegetation()->GetInterceptionStore();
            albedoSnow = ptrKwaAquifer->GetSnow()->GetSnowAlbedo();
            snowCoverFraction = ptrKwaAquifer->GetSnow()->GetSnowCoverFraction();
        }
        else if (ptrGlacier)
        {
            landSurfaceType = ptrGlacier->GetLandSurfaceType();
            albedoVegetation = ptrGlacier->GetGlacierIce()->GetLandSurfacePar()->GetALBEDO();
            precipitation = ptrGlacier->GetGlacierIce()->GetPrecipitation();
            temperature = ptrGlacier->GetGlacierIce()->GetTemperature();
            accumulationTemperature = ptrGlacier->GetGlacierIce()->GetLandSurfacePar()->GetACC_TEMP();
            leafAreaIndex = ptrGlacier->GetVegetation()->GetLeafAreaIndex();
            interceptionStore = ptrGlacier->GetVegetation()->GetInterceptionStore();
            albedoSnow = ptrGlacier->GetSnow()->GetSnowAlbedo();
            snowCoverFraction = ptrGlacier->GetSnow()->GetSnowCoverFraction();
        }
        else
        {
            cout << endl << " Not legal surfaceElement type in SurfaceAlbedo \n";
            //exit(1);
        }

        //  Snow albedo
        if (temperature < accumulationTemperature && precipitation > 1.0e-4)
        {
            albedoSnow = maxSnowAlbedo;
            snowCoverFraction = 1.0;
        }
        else if (snowCoverFraction > 0.0)
        {
            if (temperature < 0.05)
            {
                albedoReduction = 0.02 * 0.05;
            }
            else
            {
                albedoReduction = 0.02 * temperature;
            }
            albedoSnow = albedoSnow - (albedoSnow - 0.5) * albedoReduction;
            if (albedoSnow < 0.2)
            {
                albedoSnow = minSnowAlbedo;
            }
        }
        else
        {
            albedoSnow = 0.0;
        }
        // Vegetation albedo
	/*        if (landSurfaceType == M_BARE_SOIL)
        {
            if (precipitation < 1.0e-4)
            {
                albedoVegetation = albedoVegetation * 2.0;
            }
        }
        else if (landSurfaceType != ROCK && landSurfaceType != GLACIER && landSurfaceType != M_ROCK && landSurfaceType != M_URBAN &&
                 landSurfaceType != M_WATER && leafAreaIndex <= 4.0)
        {
            if (precipitation < 1.0e-4)
            {
                albedoSoil = albedoVegetation * 2.0;
            }
            else
            {
                albedoSoil = albedoVegetation;
            }
            albedoVegetation = albedoSoil + 0.25 * leafAreaIndex * (albedoVegetation - albedoSoil);
        }
        // Snow and vegetation albedo for conifers
        if (landSurfaceType == M_CONIFERS && snowCoverFraction > 0.0)
        {
            if (interceptionStore > maximumInterceptedWater * leafAreaIndex)
            {
                albedoVegetation = albedoSnow * (1.0 - vegetationDensity) + albedoVegetation * vegetationDensity;
            }
	    }*/
        // Weighted mean of snow and vegetation albedo
        albedoSurface = albedoSnow * snowCoverFraction + albedoVegetation * (1 - snowCoverFraction);

        // Store values from this time step
        if (ptrHbvAquifer)
        {
            ptrHbvAquifer->GetSnow()->SetSnowAlbedo(albedoSnow);
            ptrHbvAquifer->SetSurfaceAlbedo(albedoSurface);
        }
        else if (ptrKwaAquifer)
        {
            ptrKwaAquifer->GetSnow()->SetSnowAlbedo(albedoSnow);
            ptrKwaAquifer->SetSurfaceAlbedo(albedoSurface);
        }
        else if (ptrGlacier)
        {
            ptrGlacier->GetSnow()->SetSnowAlbedo(albedoSnow);
            ptrGlacier->SetSurfaceAlbedo(albedoSurface);
        }
    }
}
