#include "Routing.h"
#include "DistributedElement.h"
#include "SubCatchment.h"
#include "Dew.h"
#include "DateTime.h"

void TraverseWaterCourseTravelTime(SubCatchment * const thisSubCatchment, double travelTime, double * waterVelocity, ofstream &fout)
{
    int i;
    double velocity, spaceDist;
    DistributedElement * thisElement;

    // Travel time for landscape element
    if (thisSubCatchment->GetLandScapeElement()->GetFlowDirection() == 1 || thisSubCatchment->GetLandScapeElement()->GetFlowDirection() == 4 || thisSubCatchment->GetLandScapeElement()->GetFlowDirection() == 16 || thisSubCatchment->GetLandScapeElement()->GetFlowDirection() == 64)
    {
        spaceDist = sqrt(thisSubCatchment->GetLandScapeElement()->GetArea());
    }
    else
    {
        spaceDist = sqrt(thisSubCatchment->GetLandScapeElement()->GetArea()) * sqrt(2.0);
    }
    velocity = waterVelocity[0];           // m/s
    travelTime = travelTime + spaceDist / (velocity * sqrt(tan(thisSubCatchment->GetLandScapeElement()->GetSlopeAngle() * acos(-1.0) / 180.0)));
    thisSubCatchment->GetLandScapeElement()->SetTravelTime(travelTime);

    // All upstream sub-catchment/watercourse elements
    for (i = 0; i < thisSubCatchment->GetNumUpStream(); i++)
    {
        TraverseWaterCourseTravelTime(thisSubCatchment->GetUpStream(i), travelTime, waterVelocity, fout);
    }

    // All upslope landscape elements
    thisElement = thisSubCatchment->GetLandScapeElement();
    while (thisElement)
    {
        TraverseLandScapeTravelTime(thisElement, travelTime, waterVelocity, fout);
        thisElement = thisElement->GetNextElement();
    }
}


void TraverseLandScapeTravelTime(DistributedElement * const thisElement, double travelTime, double * waterVelocity, ofstream &fout)
{
    int i;
    double velocity, spaceDist;

    // Travel time for landscape element
    if (thisElement->GetFlowDirection() == 1 || thisElement->GetFlowDirection() == 4 || thisElement->GetFlowDirection() == 16 || thisElement->GetFlowDirection() == 64)
    {
        spaceDist = sqrt(thisElement->GetArea());
    }
    else
    {
        spaceDist = sqrt(thisElement->GetArea()) * sqrt(2.0);
    }
    velocity = waterVelocity[0] * 60.0;           // m/s -> m/minute
    if (thisElement->GetTravelTime() < 0)
    {
        travelTime = travelTime + spaceDist / (velocity * sqrt(tan(thisElement->GetSlopeAngle() * acos(-1.0) / 180.0)));
        thisElement->SetTravelTime(travelTime);
    }

    for (i = 0; i < thisElement->GetNumUpLand(); i++)
    {
        TraverseLandScapeTravelTime(thisElement->GetUpLandFlow(i), travelTime, waterVelocity, fout);
    }
}


void RouteSourceToSinkDischarge(DistributedElement * Dew, int numLand, double * sourceToSinkDischarge, int initialTimeSteps, int numberTimeSteps, int timeStep, int secondsPerTimeStep)
{
    int i, travelTimeSteps;

    for (i = 0; i < numLand; i++)
    {
        if (Dew[i].GetTravelTime() > 0)
        {
            travelTimeSteps = (int)Dew[i].GetTravelTime() / secondsPerTimeStep;
            if (timeStep + travelTimeSteps < initialTimeSteps + numberTimeSteps)
            {
                sourceToSinkDischarge[timeStep + travelTimeSteps] = sourceToSinkDischarge[timeStep + travelTimeSteps] + Dew[i].GetDischarge();
            }
        }
    }
}


void WriteSourceToSinkDischarge(SubCatchment * const thisSubCatchment, double * sourceToSinkDischarge, DateTime startSimulationTime,
                                DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps, int secondsPerTimeStep)
{
    FILE *fpOut;
    char fileName[100];
    int i, k, timeStep;
    DateTime datetime;

    sprintf(fileName, "sts_%08d.var", thisSubCatchment->GetIdentifier());
    if ((fpOut = fopen(fileName, "w")) == NULL)
    {
        printf("\n Filen %s ikke funnet!\n\n", fileName);
        exit(1);
    }
    timeStep = initialTimeSteps;
    for (datetime = startSimulationTime; datetime <= endSimulationTime; datetime += secondsPerTimeStep)
    {
        if (sourceToSinkDischarge[timeStep] != missingData)
        {
            fprintf(fpOut, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), sourceToSinkDischarge[timeStep]);
            timeStep++;
        }
    }
}


void TraverseRouteWaterCourse(SubCatchment * const thisSubCatchment, int timeStep, int secondsPerTimeStep, ofstream &fout)
{
    int i;
    double accumulatedRouteInput = 0.0;
    double routeInput1, routeInput2, routeOutput1, routeOutput2;
    double qRef, yRef, cRef, spaceDist, courantNumber, epsilon;
    DistributedElement * thisElement;

    for (i = 0; i < thisSubCatchment->GetNumUpStream(); i++)
    {
        TraverseRouteWaterCourse(thisSubCatchment->GetUpStream(i), timeStep, secondsPerTimeStep, fout);
        accumulatedRouteInput = accumulatedRouteInput + thisSubCatchment->GetUpStream(i)->GetAccumulatedDischarge(timeStep);
    }
    accumulatedRouteInput = accumulatedRouteInput + thisSubCatchment->GetAccumulatedInFlow(timeStep);

    // Routing initial values
    routeInput2 = accumulatedRouteInput;
    if (timeStep > 0)
    {
        routeInput1 = thisSubCatchment->GetAccumulatedRouteInput();
        routeOutput1 = thisSubCatchment->GetAccumulatedDischarge(timeStep - 1);
    }
    else
    {
        routeInput1 = (thisSubCatchment->GetAccumulatedDischarge(0)) / 2.0;
        routeOutput1 = routeInput1;
    }

    // Routing for river elements
    if (thisSubCatchment->GetManning() > 0.0)
    {
        // Routing distance
        if (thisSubCatchment->GetLandScapeElement()->GetFlowDirection() == 1 || thisSubCatchment->GetLandScapeElement()->GetFlowDirection() == 4 || thisSubCatchment->GetLandScapeElement()->GetFlowDirection() == 16 || thisSubCatchment->GetLandScapeElement()->GetFlowDirection() == 64)
        {
            spaceDist = sqrt(thisSubCatchment->GetLandScapeElement()->GetArea());
        }
        else
        {
            spaceDist = sqrt(thisSubCatchment->GetLandScapeElement()->GetArea()) * sqrt(2.0);
        }

        // Routing with reference values
        qRef = (routeInput1 + routeOutput1 + routeInput2) / 3.0;
        if (qRef != 0.0)
        {
            cRef = RoutReferenceValues(thisSubCatchment->GetManning(), thisSubCatchment->GetRiverSlope(), thisSubCatchment->GetWidth(), qRef);
            courantNumber = spaceDist / (cRef * secondsPerTimeStep);
            if (courantNumber < 1.0)
            {
                epsilon = 0.5 * (1 - qRef / (cRef * thisSubCatchment->GetRiverSlope() * thisSubCatchment->GetWidth() * spaceDist));
                routeOutput2 = RoutOutputValues(cRef, epsilon, secondsPerTimeStep, spaceDist, routeInput1, routeInput2, routeOutput1);
            }
            else
            {
                cout << " Courant condition not satisfied when In1 is " << routeInput1 << " In2 is " << routeInput2 << " Ou1 is " << routeOutput1 << endl;
                routeOutput2 = routeInput2;
            }
        }
        else
        {
            routeOutput2 = 0.0;
        }
    }
    // Lake element, water is routed directly to outlet without change in storage
    else
    {
        routeOutput2 = routeInput2;
    }

    if (routeOutput2 < 0.0)
    {
        routeOutput2 = 0.0;
    }

    // Store fluxes
    thisSubCatchment->SetAccumulatedDischarge(timeStep, routeOutput2);
    thisSubCatchment->SetAccumulatedRouteInput(accumulatedRouteInput);
}


double RoutReferenceValues(double manning, double riverSlope, double width, double Q)
{
    int i;
    double eps = 1.0e-6;
    double y0 = 1.0, y1 = 2.0, difference;
    double areaY, perimeterY, qrefY, crefY, diffY;

    difference = y1 - y0;
    while (difference > eps)
    {
        y0 = y1;
        areaY = width * y0;
        perimeterY = width + 2 * y0;
        qrefY = RoutReferenceDischarge(manning, riverSlope, areaY, perimeterY);
        crefY = RoutReferenceWaterLevel(manning, riverSlope, width, areaY, perimeterY);
        diffY = width * crefY;
        y1 = y0 - (qrefY - Q) / diffY;
        difference = y1 - y0;
    }
    areaY = width * y1;
    perimeterY = width + 2 * y1;
    crefY = RoutReferenceWaterLevel(manning, riverSlope, width, areaY, perimeterY);
    return crefY;
}


double RoutReferenceDischarge(double manning, double riverSlope, double area, double perimeter)
{
    return sqrt(riverSlope) * pow(area, 5.0 / 3.0) / (manning * pow(perimeter, 2.0 / 3.0));
}


double RoutReferenceWaterLevel(double manning, double riverSlope, double width, double area, double perimeter)
{
    return (5.0 / 3.0) * (sqrt(riverSlope)) * pow(area, 2.0 / 3.0) * (1 - (4.0 / 5.0) * area / (width * perimeter)) / (manning * pow(perimeter, 2.0 / 3.0));
}


double RoutOutputValues(double cRef, double epsilon, double secondsPerTimeStep, double spaceDist, double routeInput1, double routeInput2, double routeOutput1)
{
    double c1p, c2p, c3p, out2;

    c1p = (cRef * secondsPerTimeStep - 2 * spaceDist * epsilon) / (cRef * secondsPerTimeStep + 2 * spaceDist * (1 - epsilon));
    c2p = (cRef * secondsPerTimeStep + 2 * spaceDist * epsilon) / (cRef * secondsPerTimeStep + 2 * spaceDist * (1 - epsilon));
    c3p = 1 - c1p - c2p;
    out2 = c1p * routeInput2 + c2p * routeInput1 + c3p * routeOutput1;
    return out2;
}

