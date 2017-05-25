#pragma once
class ModelControl
{
public:
    void SetModelControl(char mrun, char hier, char form, char evap, char glac, char rout);
    void SetModelRun(char value);
    char GetModelRun() const;
    void SetHierarchy(char value);
    char GetHierarchy() const;
    void SetInputFormat(char value);
    char GetInputFormat() const;
    void SetEvaporationModelling(char value);
    char GetEvaporationModelling() const;
    void SetGlacierModelling(char value);
    char GetGlacierModelling() const;
    void SetRoutingType(char value);
    char GetRoutingType() const;
    ModelControl();
    ~ModelControl();

private:
    char modelRun;
    char hierarchy;
    char inputFormat;
    char evaporationModelling;
    char glacierModelling;
    char routingType;
};
