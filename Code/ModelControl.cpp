#include "ModelControl.h"

ModelControl::ModelControl() :
    modelRun('M'),
    hierarchy('H'),
    inputFormat('F'),
    evaporationModelling('E'),
    glacierModelling('G'),
    routingType('R')
{
}

ModelControl::~ModelControl()
{
}

void  ModelControl::SetModelRun(char value)
{
    modelRun = value;
}
char  ModelControl::GetModelRun() const
{
    return modelRun;
}
void  ModelControl::SetHierarchy(char value)
{
    hierarchy = value;
}
char  ModelControl::GetHierarchy() const
{
    return hierarchy;
}
void  ModelControl::SetInputFormat(char value)
{
    inputFormat = value;
}
char  ModelControl::GetInputFormat() const
{
    return inputFormat;
}
void  ModelControl::SetEvaporationModelling(char value)
{
    evaporationModelling = value;
}
char  ModelControl::GetEvaporationModelling() const
{
    return evaporationModelling;
}
void  ModelControl::SetGlacierModelling(char value)
{
    glacierModelling = value;
}
char  ModelControl::GetGlacierModelling() const
{
    return glacierModelling;
}
void  ModelControl::SetRoutingType(char value)
{
    routingType = value;
}
char  ModelControl::GetRoutingType() const
{
    return routingType;
}


void ModelControl::SetModelControl(char mrun, char hier, char form, char evap, char glac, char rout)
{
    SetModelRun(mrun);
    SetHierarchy(hier);
    SetInputFormat(form);
    SetEvaporationModelling(evap);
    SetGlacierModelling(glac);
    SetRoutingType(rout);
}
