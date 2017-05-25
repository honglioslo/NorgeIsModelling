#include "Snow.h"
#include "ParametersLandSurface.h"

Snow::Snow() :
//  temp(0.0),
    snowAlbedo(0.8),
    snowStore(0.0),
    meltWater(0.0),
    snowWaterEquivalentChange(0.0),
    waterOutput(0.0),
    snowCoverFraction(0.0)
{
    int i;
    for (i = 0; i < numberSnowClasses; i++)
    {
        distSnowStore[i] = 0.0;
        distMeltWater[i] = 0.0;
        distWaterOutput[i] = 0.0;
    }
    SetGeneralPar(0);
    SetLandSurfacePar(0);
    //  SetInputTimeSeries(0);
    SetInputElement(0);
    SetLandScapeElement(0);
}

Snow::~Snow()
{
}

void  Snow::SetSnowAlbedo(double value)
{
    snowAlbedo = value;
}
double  Snow::GetSnowAlbedo() const
{
    return snowAlbedo;
}
void  Snow::SetGeneralPar(ParametersGeneral *parObj)
{
    commonPar = parObj;
}
ParametersGeneral * Snow::GetGeneralPar() const
{
    return commonPar;
}
void  Snow::SetLandSurfacePar(ParametersLandSurface *parObj)
{
    landSurfacePar = parObj;
}
ParametersLandSurface * Snow::GetLandSurfacePar()
{
    return landSurfacePar;
}
//  void SetInputTimeSeries(InputTimeSeries *inTimeSeriesObj) { inTimeSeries = inTimeSeriesObj; }
//  InputTimeSeries *GetInputTimeSeries() const { return inTimeSeries; }
void  Snow::SetInputElement(InputElement *inElementObj)
{
    inElement = inElementObj;
}
InputElement * Snow::GetInputElement() const
{
    return inElement;
}
void  Snow::SetLandScapeElement(DistributedElement *theElement)
{
    landScapeElement = theElement;
}
DistributedElement * Snow::GetLandScapeElement() const
{
    return landScapeElement;
}

void Snow::SetSnowStore(double value)
{
    //  cout << " SetSnowStore 1 " << value << "  " << snowStore << "  " << meltWater << endl;
    meltWater = value * landSurfacePar->GetMAX_REL();
    snowStore = value - meltWater;
    if (snowStore < 0.0) snowStore = 0.0;
    for (int i = 0; i < numberSnowClasses; i++)
    {
        distMeltWater[i] = value * landSurfacePar->GetMAX_REL();
	distSnowStore[i] = value - distMeltWater[i];
	if (distSnowStore[i] < 0.0) distSnowStore[i] = 0.0;
        //    cout << " SetSnowStore " << i << "  distSnowStore " << distSnowStore[i] << "  distMeltWater " << distMeltWater[i] << endl;
    }
    //  cout << " SetSnowStore 2 " << value << "  " << snowStore << "  " << meltWater << endl;
}

double Snow::GetSnowCoverFraction() const
{
    return snowCoverFraction;
}

double Snow::GetSnowStore() const
{
    return snowStore;
}

double Snow::GetMeltWater() const
{
    return meltWater;
}

double Snow::GetSnowWaterEquivalentChange() const
{
    return snowWaterEquivalentChange;
}

double Snow::GetWaterOutput() const
{
    return waterOutput;
}


void Snow::WaterBalance(int timeStep, double waterInput, double temp)
{
    int i;
    double waterAdded;                                /*  Water added to snow store (m)  */
    double snowMelt = 0.0;                              /*  Snowmelt produced in current time step (m)  */
    double reFreeze = 0.0;                              /*  Refreeze of meltwater in snow (m)  */
    double initialSwe;                                /*  Snow water equivalent at start of time step (m)  */

    //  temp = GetInputElement()->GetInput(1);
    //  cout << "timeStep " << timeStep << "             Snow waterInput " << waterInput << "    Temperature " << temp << endl;
    initialSwe = snowStore + meltWater;

    // Uniform snow distribution
    if (landSurfacePar->GetCV_SNOW() == 0.0)
    {

        /*  Snow accumulation  */
        if (temp < landSurfacePar->GetACC_TEMP() && waterInput > 0.0)
        {
            snowStore = snowStore + waterInput;
            waterAdded = 0.0;
        }
        else
        {
            waterAdded = waterInput;
        }

        if (snowStore > 0.0)
        {
            /*  Water added to meltwater from water input  */
            meltWater = meltWater + waterAdded;

            /*  Refreeze of meltwater in snow  */
            if (temp < landSurfacePar->GetMELT_TEMP() && meltWater > 0.0)
            {
                reFreeze = landSurfacePar->GetFREEZE_EFF() * landSurfacePar->GetSNOW_MELT_RATE() * (landSurfacePar->GetMELT_TEMP() - temp);
                if (reFreeze > meltWater)
                {
                    reFreeze = meltWater;
                }
                meltWater = meltWater - reFreeze;
                snowStore = snowStore + reFreeze;
            }

            /*  Snow store and snow melt  */
            if (temp > landSurfacePar->GetMELT_TEMP() && snowStore > 0.0)
            {
                snowMelt = landSurfacePar->GetSNOW_MELT_RATE() * (temp - landSurfacePar->GetMELT_TEMP());
                if (snowMelt > snowStore)
                {
                    snowMelt = snowStore;
                }
                snowStore = snowStore - snowMelt;
                meltWater = meltWater + snowMelt;
            }

            /* Water output from snow store */
            if (meltWater > landSurfacePar->GetMAX_REL() * snowStore)
            {
                waterOutput = meltWater - landSurfacePar->GetMAX_REL() * snowStore;
                meltWater = meltWater - waterOutput;
            }
            else
            {
                waterOutput = 0.0;
            }
        }
        else
        {
            waterOutput = waterInput;
        }

        /*  Fraction of area covered by snow  */
        if (snowStore > 0.0)
        {
            snowCoverFraction = 1.0;
        }
        else
        {
            snowCoverFraction = 0.0;
        }
    }

    // Lognormal snow distribution
    else
    {

        for (i = 0; i < numberSnowClasses; i++)
        {

            /*  Snow accumulation  */
            if (temp < landSurfacePar->GetACC_TEMP() && waterInput > 0.0)
            {
                distSnowStore[i] = distSnowStore[i] + waterInput * landSurfacePar->GetSNOW_WEIGHT(i);
                waterAdded = 0.0;
            }
            else
            {
                waterAdded = waterInput;
            }

            if (distSnowStore[i] > 0.0)
            {
                /*  Water added to meltwater from water input  */
                distMeltWater[i] = distMeltWater[i] + waterAdded;

                /*  Refreeze of meltwater in snow  */
                if (temp < landSurfacePar->GetMELT_TEMP() && distMeltWater[i] > 0.0)
                {
                    reFreeze = landSurfacePar->GetFREEZE_EFF() * landSurfacePar->GetSNOW_MELT_RATE() * (landSurfacePar->GetMELT_TEMP() - temp);
                    if (reFreeze > distMeltWater[i])
                    {
                        reFreeze = distMeltWater[i];
                    }
                    distMeltWater[i] = distMeltWater[i] - reFreeze;
                    distSnowStore[i] = distSnowStore[i] + reFreeze;
                }

                /*  Snow store and snow melt  */
                if (temp > landSurfacePar->GetMELT_TEMP() && distSnowStore[i] > 0.0)
                {
                    snowMelt = landSurfacePar->GetSNOW_MELT_RATE() * (temp - landSurfacePar->GetMELT_TEMP());
                    if (snowMelt > distSnowStore[i])
                    {
                        snowMelt = distSnowStore[i];
                    }
                    distSnowStore[i] = distSnowStore[i] - snowMelt;
                    distMeltWater[i] = distMeltWater[i] + snowMelt;
                }

                /* Water output from snow store */
                if (distMeltWater[i] > landSurfacePar->GetMAX_REL() * distSnowStore[i])
                {
                    distWaterOutput[i] = distMeltWater[i] - landSurfacePar->GetMAX_REL() * distSnowStore[i];
                    distMeltWater[i] = distMeltWater[i] - distWaterOutput[i];
                }
                else
                {
                    distWaterOutput[i] = 0.0;
                }
            }
            else
            {
                distWaterOutput[i] = waterInput;
            }
        }

        /*  Accumulate values from distributed snow classes  */
        snowStore = 0.0;
        waterOutput = 0.0;
        meltWater = 0.0;
        for (i = 0; i < numberSnowClasses; i++)
        {
            //      cout << " snow.h " << i << "  distSnowStore " << distSnowStore[i] << "  distMeltWater " << distMeltWater[i] << endl;
            snowStore = snowStore + distSnowStore[i] * probNorm[i];
            waterOutput = waterOutput + distWaterOutput[i] * probNorm[i];
            meltWater = meltWater + distMeltWater[i] * probNorm[i];
        }

        /*  Fraction of area covered by snow  */
        snowCoverFraction = 0.0;
        for (i = numberSnowClasses - 1; i >= 0; i--)
        {
            if (distSnowStore[i] > 0.0)
            {
                snowCoverFraction = snowCoverFraction + probNorm[i];
            }
        }

    }

    snowWaterEquivalentChange = snowStore + meltWater - initialSwe;

}

