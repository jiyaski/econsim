
#include "distributions.h"
#include <cmath> 


RandomSampler::RandomSampler()
    : gen(rd()), 
      dist_age(AGE_MIN, AGE_MAX), 
      dist_min_consumption(MIN_CONSUMPTION_MU, MIN_CONSUMPTION_SIGMA), 
      dist_consumption(CONSUMPTION_MIN, CONSUMPTION_MAX),
      dist_annual_ROI(ROI_MEAN, ROI_STDEV),
      dist_wage(WAGE_MU, WAGE_SIGMA),
      dist_init_wealth(INIT_WEALTH_MU, INIT_WEALTH_SIGMA) {}

double RandomSampler::generateAge()             { return dist_age(gen);             } 
double RandomSampler::generateWage()            { return dist_wage(gen);            }
double RandomSampler::generateAnnualROI()       { return dist_annual_ROI(gen);      }
double RandomSampler::generateConsumption()     { return dist_consumption(gen);     }
double RandomSampler::generateMinConsumption()  { return dist_min_consumption(gen); } 
double RandomSampler::generateInitWealth() {
    double point = dist_init_wealth(gen); 
    return std::exp(point) - 2 * INIT_WEALTH_10TH_PERCENTILE; 
}
