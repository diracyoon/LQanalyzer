#ifndef HISTUTILS_H__
#define HISTUTILS_H__

class TH2;
class TH1;
class THStack;

void FixOverUnderFlows(TH1* hist, float max_x=10000.);
void FixOverFlows(TH1* hist, float max_X);
void FixUnderFlows(TH1* hist);

double IntegralStack(THStack* stack);

#endif // HISTUTILS_H__
