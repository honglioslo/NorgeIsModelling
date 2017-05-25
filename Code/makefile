OBJS = Dew.o Albedo.o Atmosphere.o CTime.o Date.o DateTime.o DistributedElement.o Glacier.o GlacierElements.o GlacierIce.o HbvAquifer.o HBV.o InputElement.o InputTimeSeries.o KinematicWave.o KwaAquifer.o Lake.o LakeWaterBalance.o LandAtmosphereInterface.o MeteorologicalStations.o ModelControl.o Parameters.o ParametersGeneral.o ParametersGlacierRetreat.o ParametersKiWa.o ParametersLandSurface.o ParametersSubSurfaceHbv.o PenmanMonteith.o Routing.o SelectedHbvTimeSeriesElements.o SelectedKiWaTimeSeriesElements.o SelectedSubCatchmentTimeSeriesElements.o Snow.o stdafx.o SubCatchment.o TotalReservoirStorage.o Util.o Vegetation.o

C++C=g++

dew: 	$(OBJS)
	$(C++C) $(C++FLAGS) $(OBJS) -o dew 

clean:
	rm -f $(OBJS) *.o core
