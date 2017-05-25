#include "GlacierElements.h"
#include "DistributedElement.h"
#include "Glacier.h"
#include "GlacierIce.h"
#include "ParametersGeneral.h"
#include "ParametersGlacierRetreat.h"
#include "SubCatchment.h"
#include "stdafx.h"
#include "DateTime.h"

GlacierElements::GlacierElements(SubCatchment ** OutletList, int numOut) :
    numWatcOut(numOut)
{
    Outlet = OutletList;
}

GlacierElements::~GlacierElements()
{
}

void GlacierElements::BuildGlacierLists()
{
    int i;
    for (i = 0; i < numWatcOut; i++)
    {
        TraverseGlacierSubCatchment(Outlet[i], &glacierElementsList);
    }
}

void GlacierElements::SortGlacierLists()
{
    SortAllGlacierLists(&glacierElementsList);
}

void GlacierElements::RemoveElementsGlacierLists()
{
    RemoveElementsAllGlacierLists(&glacierElementsList);
}


void GlacierElements::TraverseGlacierSubCatchment(SubCatchment * const thisSubCatchment, list <DistributedElement *> *glacierElementsList)
{
    int i;
    bool first = true;
    DistributedElement *thisElement, *lastDew = 0, *distributedElementList = 0;
    for (i = 0; i < thisSubCatchment->GetNumUpStream(); i++)
    {
        TraverseGlacierSubCatchment(thisSubCatchment->GetUpStream(i), glacierElementsList);
    }
    thisElement = thisSubCatchment->GetLandScapeElement();
    while (thisElement)
    {
        TraverseGlacierLandScape(thisElement, &distributedElementList, &lastDew);
        thisElement = thisElement->GetNextElement();
    }
    // Sort and store list of glacier elements
    if (distributedElementList)
    {
        SortOneGlacierList(&distributedElementList);
        (*glacierElementsList).push_back(distributedElementList);
    }
}


void GlacierElements::TraverseGlacierLandScape(DistributedElement * thisElement, DistributedElement ** distributedElementList, DistributedElement ** lastDew)
{
    int i;
    for (i = 0; i < thisElement->GetNumUpLand(); i++)
    {
        TraverseGlacierLandScape(thisElement->GetUpLandFlow(i), distributedElementList, lastDew);
    }
    // Build list of glacier elements
    if (thisElement->GetGlacier())
    {
        if (!*distributedElementList)
        {
            *distributedElementList = thisElement;
            *lastDew = thisElement;
        }
        else
        {
            (*lastDew)->SetNextGlacierElement(thisElement);
            *lastDew = thisElement;
        }
    }
}


void GlacierElements::SortAllGlacierLists(list <DistributedElement *> *glacierElementsList)
{
    list <DistributedElement *>::iterator listIterator;
    for (listIterator = (*glacierElementsList).begin(); listIterator != (*glacierElementsList).end(); listIterator++)
    {
        if (*listIterator)
        {
            SortOneGlacierList(&(*listIterator));
        }
    }
}


// Sort list of landscape elements with glaciers
void GlacierElements::SortOneGlacierList(DistributedElement ** distributedElementList)
{
    int i, iHong;
    bool first = true;
    DistributedElement *thisDew, *endDew, *tempDew, *thisNext, *endNext, *thisPrevious, *endPrevious;
    thisDew = *distributedElementList;
    thisPrevious = 0;
    while (thisDew->GetNextGlacierElement())
    {
        endPrevious = thisDew;
        endDew = thisDew->GetNextGlacierElement();
        while (endDew) 
        { 
            if (endDew->GetGlacier()->GetGlacierSurfaceElevation() < thisDew->GetGlacier()->GetGlacierSurfaceElevation())
            {
                thisNext = thisDew->GetNextGlacierElement();
                endNext = endDew->GetNextGlacierElement();
                tempDew = thisDew;
                thisDew = endDew;
                endDew = tempDew;
                if (thisDew == thisNext)
                {
                    thisDew->SetNextGlacierElement(endDew);
                }
                else
                {
                    thisDew->SetNextGlacierElement(thisNext);
                    endPrevious->SetNextGlacierElement(endDew);
                }
                endDew->SetNextGlacierElement(endNext);
                if (thisPrevious)
                {
                    thisPrevious->SetNextGlacierElement(thisDew);
                }
            }
            endPrevious = endDew;
            endDew = endDew->GetNextGlacierElement();

        }
        if (first)
        {
            *distributedElementList = thisDew;
        }
        thisPrevious = thisDew;
        thisDew = thisDew->GetNextGlacierElement();
        first = false;
    }
}


void GlacierElements::RemoveElementsAllGlacierLists(list <DistributedElement *> *glacierElementsList)
{
    list <DistributedElement *>::iterator listIterator;
    for (listIterator = (*glacierElementsList).begin(); listIterator != (*glacierElementsList).end(); listIterator++)
    {
        //    if (*listIterator) cout << "\n1 RemoveElementsAllGlacierLists (*listIterator) " << (*listIterator) << endl;
        if (*listIterator)
        {
            RemoveElementsOneGlacierList(&(*listIterator));
        }
    }
    //  Test
    /*  DistributedElement *thisDew;
    list <DistributedElement *>::iterator testIterator;
    for (testIterator=(*glacierElementsList).begin(); testIterator!=(*glacierElementsList).end(); testIterator++) {
    thisDew = (*testIterator);
    if (*testIterator) cout << "2 RemoveElementsAllGlacierLists (*testIterator) " << (*testIterator) << endl;
    while (thisDew) {
    cout << "3 RemoveElementsAllGlacierLists thisDew " << thisDew << "  " << thisDew->GetGeoIndex() << endl;
    thisDew=thisDew->GetNextGlacierElement();
    }
    }*/
    //  Test end
}


// Remove ice free elements from list of landscape elements with glaciers original not used by extended Huss algorithm
/*void GlacierElements::RemoveElementsOneGlacierList(DistributedElement ** distributedElementList)
{
    int i;
    DistributedElement *thisDew, *thisPrevious;
    thisDew = *distributedElementList;
    *distributedElementList = 0;
    thisPrevious = 0;
    while (thisDew)
    {
        //    cout << "1 RemoveElementsOneGlacierList thisDew " << thisDew << "  " << thisDew->GetGeoIndex() << "  " << thisDew->GetGlacier()->GetGlacierIceThickness() << "  " << thisDew->GetGlacier()->GetGlacierSurfaceElevation() << endl;
        if (thisDew->GetGlacier()->GetGlacierIceThickness() > 0.0)
        {
            //      cout << "2 RemoveElementsOneGlacierList thisDew " << thisDew << "  " << thisDew->GetGeoIndex() << "  " << thisDew->GetGlacier()->GetGlacierIceThickness() << "  " << thisDew->GetGlacier()->GetGlacierSurfaceElevation() << endl;
            if (thisPrevious)
            {
                thisPrevious->SetNextGlacierElement(thisDew);
            }
            else
            {
                *distributedElementList = thisDew;
            }
            thisPrevious = thisDew;
        }
        else
        {
            thisDew->GetGlacier()->SetGlacierIceAreaFraction(0.0);
        }
        thisDew = thisDew->GetNextGlacierElement();
    }
    //  Test
    /*  thisDew=*distributedElementList;
    while (thisDew) {
    cout << "3 RemoveElementsOneGlacierList thisDew " << thisDew << "  " << thisDew->GetGeoIndex() << "  " << thisDew->GetGlacier()->GetGlacierIceThickness() << "  " << thisDew->GetGlacier()->GetGlacierSurfaceElevation() << endl;
    thisDew=thisDew->GetNextGlacierElement();
    }*/
    //  Test end
    /*}*/


// Set glacier ice are fraction of ice free elements equal to zero in extended Huss algorithm
void GlacierElements::RemoveElementsOneGlacierList(DistributedElement ** distributedElementList)
{
    int i;
    DistributedElement *thisDew;
    thisDew = *distributedElementList;
    while (thisDew)
    {
      if (thisDew->GetGlacier()->GetGlacierIceThickness() <= 0.0)
      {
	thisDew->GetGlacier()->SetGlacierIceAreaFraction(0.0);
      }
      thisDew = thisDew->GetNextGlacierElement();
    }
}


void GlacierElements::SetThisYearAnnualGlacierValues(DateTime datetime)
{
    DistributedElement *thisDew;
    list <DistributedElement *>::iterator listIterator;
    for (listIterator = glacierElementsList.begin(); listIterator != glacierElementsList.end(); listIterator++)
    {
        thisDew = (*listIterator);
        while (thisDew)
        {
            //      cout << "1 SetThisYearAnnualGlacierValues thisDew " << thisDew << "  " << thisDew->GetGeoIndex() << "  " << thisDew->GetGlacier()->GetGlacierIceAreaFraction() << endl;
            thisDew->SetThisYearAnnualGlacierValues(datetime);
            thisDew = thisDew->GetNextGlacierElement();
        }
    }
    //cout << "\n";
}


/*  Algorithm for redistribution of glacier surface based on Huss et al, Hydrology and Earth System Sciences, 2010 & Huss & Hock, Frontiers in Earth Science, 2015 */
void GlacierElements::GlacierSurfaceElevationReDistribution()
{
    int numberGlacierElements;
    int i, indexGlacierIceAdvance, numberGlacierIceAdvance;
    double a, b, c, gamma, scaleFactor;
    double elevationMin, elevationMax, elevationNormalized, elevationChangeNormalized;
    double elevationNew, iceThicknessNew;
    double densityIce, massBalanceSum, massBalanceSum2, meanMassBalance, areaSum, areaElevationSum;
    double iceThicknessMin, areaFractionNew, iceVolumeNew;
    double positiveMassBalanceExcess, increaseThreshold, iceVolumeExcess;
    DistributedElement *thisDew;
    list <DistributedElement *>::iterator listIterator;
    for (listIterator = glacierElementsList.begin(); listIterator != glacierElementsList.end(); listIterator++)
    {
        // Test output
        /*        cout << "11111111111111111111\n";
	thisDew = (*listIterator);
        while (thisDew)
        {
          cout << "elevation thickness " << thisDew->GetGeoIndex() << "  " << thisDew->GetElevation() << "  " << thisDew->GetGlacier()->GetGlacierSurfaceElevation() << "  " << thisDew->GetGlacier()->GetGlacierIceThickness() << endl;  
          thisDew = thisDew->GetNextGlacierElement();
	  }*/
        // Test output end
        thisDew = (*listIterator);
        if (thisDew)
        {
            densityIce = thisDew->GetGeneralPar()->GetDENSITY_ICE();
            a = thisDew->GetGlacier()->GetGlacierIce()->GetGlacierRetreatPar()->GetA();
            b = thisDew->GetGlacier()->GetGlacierIce()->GetGlacierRetreatPar()->GetB();
            c = thisDew->GetGlacier()->GetGlacierIce()->GetGlacierRetreatPar()->GetC();
            gamma = thisDew->GetGlacier()->GetGlacierIce()->GetGlacierRetreatPar()->GetGAMMA();
            increaseThreshold = thisDew->GetGlacier()->GetGlacierIce()->GetGlacierRetreatPar()->GetINCREASE_THRESH();
            numberGlacierIceAdvance = thisDew->GetGlacier()->GetGlacierIce()->GetGlacierRetreatPar()->GetNUMBER_ADVANCE();
            indexGlacierIceAdvance = -1;
            DistributedElement * iceAdvanceElement = new DistributedElement[numberGlacierIceAdvance];
            while (thisDew)
            {
              //                if (thisDew->GetGlacier()->GetGlacierIceThickness() <= 0.0) 
                if (thisDew->GetGlacier()->GetGlacierIceAreaFraction() == 0.0) 
                {
                    if (indexGlacierIceAdvance < numberGlacierIceAdvance-1) 
                    {
                        indexGlacierIceAdvance++;
                    }
                    if (indexGlacierIceAdvance == 0) 
                    {
                        iceAdvanceElement[0] = *thisDew;
                    }
                    else
                    {
                        for (i=indexGlacierIceAdvance; i>=1; i--)
                        {      
                            iceAdvanceElement[i] = iceAdvanceElement[i-1];
                        }
                        iceAdvanceElement[0] = *thisDew;
                    }
                }
            thisDew = thisDew->GetNextGlacierElement();
            }
            // Huss et al, eq. 1
            elevationMin = (-1) * largeMissingData;
            elevationMax = largeMissingData;
            iceThicknessMin = (-1) * largeMissingData;
            thisDew = (*listIterator);
            while (thisDew)
            {
              //              if (thisDew->GetGlacier()->GetGlacierIceThickness() > 0.0) 
              if (thisDew->GetGlacier()->GetGlacierIceAreaFraction() > 0.0) 
              {
                if (thisDew->GetGlacier()->GetGlacierSurfaceElevation() < elevationMin)
                {
                    elevationMin = thisDew->GetGlacier()->GetGlacierSurfaceElevation();
                }
                if (thisDew->GetGlacier()->GetGlacierSurfaceElevation() > elevationMax)
                {
                    elevationMax = thisDew->GetGlacier()->GetGlacierSurfaceElevation();
                }
                if (thisDew->GetGlacier()->GetGlacierIceThickness() < iceThicknessMin)
                {
                    iceThicknessMin = thisDew->GetGlacier()->GetGlacierIceThickness();
                }
              }
              thisDew = thisDew->GetNextGlacierElement();
            }
            massBalanceSum = 0.0;
            areaSum = 0.0;
            areaElevationSum = 0.0;
            numberGlacierElements = 0;
            thisDew = (*listIterator);
            while (thisDew)
            {
              //              if (thisDew->GetGlacier()->GetGlacierIceThickness() > 0.0) 
              if (thisDew->GetGlacier()->GetGlacierIceAreaFraction() > 0.0) 
              {
                numberGlacierElements++;
                if (elevationMax != elevationMin)
                {
                    elevationNormalized = (elevationMax - thisDew->GetGlacier()->GetGlacierSurfaceElevation()) / (elevationMax - elevationMin);
                    elevationChangeNormalized = pow((elevationNormalized + a), gamma) + b * (elevationNormalized + a) + c;
                }
                else
                {
                    elevationChangeNormalized = thisDew->GetGlacier()->GetGlacierSurfaceElevationChangeNormalized();
                }
                thisDew->GetGlacier()->SetGlacierSurfaceElevationChangeNormalized(elevationChangeNormalized);
                massBalanceSum = massBalanceSum + thisDew->GetAnnualMassBalance() * thisDew->GetArea() * thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
                areaSum = areaSum + thisDew->GetArea() * thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
                areaElevationSum = areaElevationSum + elevationChangeNormalized * thisDew->GetArea() * thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
              }
              thisDew = thisDew->GetNextGlacierElement();
            }
            if (numberGlacierElements < 0)
            {
                cout << " Number of glacier elements error:  numberGlacierElements = " << numberGlacierElements << endl;
                exit(1);
            }
            // Control scale factor Huss et al
            if (areaElevationSum != 0.0)
            {
                scaleFactor = massBalanceSum / (areaElevationSum * densityIce);
            }
            else
            {
                scaleFactor = 0.0;
            }
	    //	    cout << " Huss et al eq 1 ferdig \n";
            // Huss et al, eq. 2 and eq. 3
            massBalanceSum2 = 0.0;
            iceVolumeExcess = 0.0;
            thisDew = (*listIterator);
            while (thisDew)
            {
              //              if (thisDew->GetGlacier()->GetGlacierIceThickness() > 0.0) 
              if (thisDew->GetGlacier()->GetGlacierIceAreaFraction() > 0.0) 
              {
                if (scaleFactor * thisDew->GetGlacier()->GetGlacierSurfaceElevationChangeNormalized() <= increaseThreshold)
                {
                  elevationNew = thisDew->GetGlacier()->GetGlacierSurfaceElevation() + scaleFactor * thisDew->GetGlacier()->GetGlacierSurfaceElevationChangeNormalized();
                  iceThicknessNew = thisDew->GetGlacier()->GetGlacierIceThickness() + scaleFactor * thisDew->GetGlacier()->GetGlacierSurfaceElevationChangeNormalized();
                }
                else
                {
                  iceVolumeExcess = iceVolumeExcess + (scaleFactor * thisDew->GetGlacier()->GetGlacierSurfaceElevationChangeNormalized() - increaseThreshold) * thisDew->GetArea() * thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
                  elevationNew = thisDew->GetGlacier()->GetGlacierSurfaceElevation() + increaseThreshold;
                  iceThicknessNew = thisDew->GetGlacier()->GetGlacierIceThickness() + increaseThreshold;
		  //		  cout << " iceVolumeExcess " << iceVolumeExcess << endl;
                }
                if (iceThicknessNew < 0.0)
                {
                    elevationNew = elevationNew - iceThicknessNew;
                    iceThicknessNew = 0.0;
                }
                thisDew->GetGlacier()->SetGlacierSurfaceElevation(elevationNew);
                thisDew->GetGlacier()->SetGlacierIceThickness(iceThicknessNew);
                thisDew->GetGlacier()->SetGlacierIceVolume(iceThicknessNew * thisDew->GetArea()*thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0);
                massBalanceSum2 = massBalanceSum2 + thisDew->GetGlacier()->GetGlacierSurfaceElevationChangeNormalized() * thisDew->GetArea() * thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
              }
              thisDew = thisDew->GetNextGlacierElement();
            }
            // Control mass balance sum
            if (elevationMax != elevationMin)
            {
                massBalanceSum2 = scaleFactor * densityIce * massBalanceSum2;
                if (massBalanceSum < massBalanceSum2 - epsilon || massBalanceSum > massBalanceSum2 + epsilon)
                {
                    cout << " Mass balance error:  mass balance sum initial = " << massBalanceSum << " mass balance sum final = " << massBalanceSum2 << endl;
                    exit(1);
                }
            }
	    //	    cout << " Huss et al eq 2 og 3 ferdig \n";
	    //	    cout << " iceVolumeExcess " << iceVolumeExcess << endl;
            // Distribute excess ice on ice free elements
            if (iceVolumeExcess > 0.0)
            {
            for (i=0; i<=indexGlacierIceAdvance; i++)
                {
		  //		  cout << i << "  " << iceVolumeExcess << endl;
                    if (iceVolumeExcess > 0.0)
                    {
                        areaFractionNew = iceAdvanceElement[i].GetGlacier()->GetAreaFraction();
                        iceVolumeNew = iceThicknessMin * iceAdvanceElement[i].GetArea() * areaFractionNew / 100.0;
			//			cout << " areaFractionNew " << areaFractionNew << endl;
			//			cout << " iceVolumeNew " << iceVolumeNew << endl;
			//			cout << " iceVolumeExcess " << iceVolumeExcess << endl;
                        if (iceVolumeNew > iceVolumeExcess)
                        {
                            iceVolumeNew = iceVolumeExcess;
                            iceVolumeExcess = 0.0;
                        }
                        else
                        {
                            iceVolumeExcess = iceVolumeExcess - iceVolumeNew;
                        }
                        iceThicknessNew = iceVolumeNew / (iceAdvanceElement[i].GetArea() * areaFractionNew / 100.0);
                        elevationNew = iceAdvanceElement[i].GetGlacier()->GetGlacierSurfaceElevation() + iceThicknessNew;
                        iceAdvanceElement[i].GetGlacier()->SetGlacierIceAreaFraction(areaFractionNew);
                        iceAdvanceElement[i].GetGlacier()->SetGlacierSurfaceElevation(elevationNew);
                        iceAdvanceElement[i].GetGlacier()->SetGlacierIceThickness(iceThicknessNew);
                        iceAdvanceElement[i].GetGlacier()->SetGlacierIceVolume(iceVolumeNew);
			//			cout << " iceVolumeExcess " << iceVolumeExcess << endl;
                    }
                }
            }
            // Distribute excess ice on all elements with ice
            if (iceVolumeExcess > 0.0)
            {
                areaSum = 0.0;
                thisDew = (*listIterator);
                while (thisDew)
                {
                  //              if (thisDew->GetGlacier()->GetGlacierIceThickness() > 0.0) 
                  if (thisDew->GetGlacier()->GetGlacierIceAreaFraction() > 0.0) 
                  {
                      areaSum = areaSum + thisDew->GetArea() * thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0;
                  }
		  thisDew = thisDew->GetNextGlacierElement();
                }
                thisDew = (*listIterator);
                while (thisDew)
                {
                  //              if (thisDew->GetGlacier()->GetGlacierIceThickness() > 0.0) 
                  if (thisDew->GetGlacier()->GetGlacierIceAreaFraction() > 0.0) 
                  {
		    //		      elevationNew = thisDew->GetGlacier()->GetGlacierSurfaceElevation() + (iceVolumeExcess*thisDew->GetArea()*thisDew->GetGlacier()->GetGlacierIceAreaFraction()/100.0/areaSum) / (thisDew->GetArea()*thisDew->GetGlacier()->GetGlacierIceAreaFraction()/100.0);
		    //		      iceThicknessNew = thisDew->GetGlacier()->GetGlacierIceThickness() + (iceVolumeExcess*thisDew->GetArea()*thisDew->GetGlacier()->GetGlacierIceAreaFraction()/100.0/areaSum) / (thisDew->GetArea()*thisDew->GetGlacier()->GetGlacierIceAreaFraction()/100.0);
                      elevationNew = thisDew->GetGlacier()->GetGlacierSurfaceElevation() + iceVolumeExcess / areaSum;
                      iceThicknessNew = thisDew->GetGlacier()->GetGlacierIceThickness() + iceVolumeExcess / areaSum;
                      if (iceThicknessNew < 0.0)
                      {
                          elevationNew = elevationNew - iceThicknessNew;
                          iceThicknessNew = 0.0;
                      }
                      thisDew->GetGlacier()->SetGlacierSurfaceElevation(elevationNew);
                      thisDew->GetGlacier()->SetGlacierIceThickness(iceThicknessNew);
                      thisDew->GetGlacier()->SetGlacierIceVolume(iceThicknessNew * thisDew->GetArea()*thisDew->GetGlacier()->GetGlacierIceAreaFraction() / 100.0);
                  }
                  thisDew = thisDew->GetNextGlacierElement();
                }
            }
	    /*            cout << "22222222222222222222\n";
            cout << "indexGlacierIceAdvance " << indexGlacierIceAdvance << endl;
            for (i=0; i<=indexGlacierIceAdvance; i++)
            {
              cout << "elevation thickness " << i << "  " << iceAdvanceElement[i].GetGeoIndex() << "  " << iceAdvanceElement[i].GetElevation() << "  " << iceAdvanceElement[i].GetGlacier()->GetGlacierSurfaceElevation() << "  " << iceAdvanceElement[i].GetGlacier()->GetGlacierIceThickness() << endl;  
	      }*/
            delete[] iceAdvanceElement;
        }
    }
    // Test output
    /*    for (listIterator = glacierElementsList.begin(); listIterator != glacierElementsList.end(); listIterator++)
    {
        cout << "33333333333333333333\n";
	thisDew = (*listIterator);
        while (thisDew)
        {
            cout << "elevation thickness " << thisDew->GetGeoIndex() << "  " << thisDew->GetElevation() << "  " << thisDew->GetGlacier()->GetGlacierSurfaceElevation() << "  " << thisDew->GetGlacier()->GetGlacierIceThickness() << endl;  
            thisDew = thisDew->GetNextGlacierElement();
	}
	}*/
    // Test output end
}

