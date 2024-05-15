
#include "distributions.h"


RandomSampler::RandomSampler()
    : gen(rd()), 
      dist_consumption(CONSUMPTION_MIN, CONSUMPTION_MAX),
      dist_avg_return(ROI_MEAN, ROI_STDEV),
      dist_wage(WAGE_MU, WAGE_SIGMA),
      dist_wealth(WEALTH_MU, WEALTH_SIGMA) {}

double RandomSampler::generateConsumption()  { return dist_consumption(gen); }
double RandomSampler::generateAvgReturn()    { return dist_avg_return(gen);  }
double RandomSampler::generateWage()         { return dist_wage(gen);        }
double RandomSampler::generateWealth()       { return dist_wealth(gen);      }
