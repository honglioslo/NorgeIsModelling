#include "Util.h"

#include "stdafx.h"

#include "Dew.h"
using namespace std;

// Power function
double power(double base, double exponent)
{
    if (base > 0.0)
    {
        return pow(base, exponent);
    }
    else if (base == 0.0)
    {
        return 0.0;
    }
    else if (base >= -epsilon)
    {
        return 0.0;
    }
    else
    {
        printf("\n\n    *****     function power     base = %f     exponent = %f     *****\n\n", base, exponent);
        return 0.0;
    }
}


// Leap year
int leapYear(int year)
{
    if (year % 400 == 0)
    {
        return 1;
    }
    else if (year % 100 == 0)
    {
        return 0;
    }
    else if (year % 4 == 0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


// Find day number
int dayNumber(int year, int month, int day)
{
    int accDays[] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
    if (month < 3)
    {
        return accDays[month - 1] + day;
    }
    else
    {
        return accDays[month - 1] + day + leapYear(year);
    }
}


// Day number to date
void dayNo2Date(int dayNo, int year, int * month, int * day)
{
    int i;
    int accDays[] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
    if (dayNo > 365 + leapYear(year))
    {
        dayNo = dayNo - 365 - leapYear(year);
    }
    if (leapYear(year))
        for (i = 2; i < 12; i++)
        {
            accDays[i]++;
        }
    i = 12;
    while (dayNo <= accDays[i - 1])
    {
        i--;
    }
    *month = i;
    *day = dayNo - accDays[i - 1];
}


// Potential evapotranspiration
double potentialEvap(double temp, double epotPar)
{
    if (temp > 0.0)
    {
        return epotPar * temp;
    }
    else
    {
        return 0.0;
    }
}


// Actual evapotranspiration HBV
double HBVTranspSoilEvap(double soilMoist, double temp, double epotPar, double fieldCapacity, double fcDel)
{
    double epot;
    epot = potentialEvap(temp, epotPar);
    if (soilMoist > fcDel * fieldCapacity)
    {
        return epot;
    }
    else if (soilMoist > 0.0)
    {
        return epot * soilMoist / (fieldCapacity * fcDel);
    }
    else
    {
        return 0.0;
    }
}


// Actual evapotranspiration KinematicWave
double KiWaTranspSoilEvap(double teta, double temp, double epotPar, double eactPar, double tSat0, double wiltPoint)
{
    /* Condition : eactPar * (tSat0-wiltPoint) - wiltPoint > 0 */
    double epot;
    epot = potentialEvap(temp, epotPar);
    if (teta > eactPar * (tSat0 - wiltPoint))
    {
        return epot;
    }
    else if (teta > wiltPoint)
    {
        return epot * (teta - wiltPoint) / (eactPar * (tSat0 - wiltPoint) - wiltPoint);
    }
    else
    {
        return 0.0;
    }
}
