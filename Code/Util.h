#pragma once
double power(double base, double exponent);
int leapYear(int year);
int dayNumber(int year, int month, int day);
void dayNo2Date(int dayNo, int year, int * month, int * day);
double potentialEvap(double temp, double epotPar);
double HBVTranspSoilEvap(double soilMoist, double temp, double epotPar, double fieldCapacity, double fcDel);
double KiWaTranspSoilEvap(double teta, double temp, double epotPar, double eactPar, double tSat0, double wiltPoint);


