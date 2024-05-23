
#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <random>


class RandomSampler {
public:
    RandomSampler();  // Constructor 

    double generateAge(); 
    double generateMinConsumption(); 
    double generateConsumption(); 
    double generateAnnualROI();
    double generateWage();
    double generateInitWealth();

private:
    // const values are generally empirically determined by separate data analysis 
    static constexpr int AGE_MIN = 0; 
    static constexpr int AGE_MAX = 49; 
    static constexpr double MIN_CONSUMPTION_MU = 9.2103; 
    static constexpr double MIN_CONSUMPTION_SIGMA = 0.69078; 
    static constexpr double CONSUMPTION_MIN = 0.0000005;
    static constexpr double CONSUMPTION_MAX = 0.000005;
    static constexpr double ROI_MEAN = 0.06;
    static constexpr double ROI_STDEV = 0.03;
    static constexpr double WAGE_MU = 10.854;
    static constexpr double WAGE_SIGMA = 0.7564;
    static constexpr double INIT_WEALTH_MU = 10.412; 
    static constexpr double INIT_WEALTH_SIGMA = 1.0273; 
    static constexpr double INIT_WEALTH_10TH_PERCENTILE = -11270; 

    std::random_device rd;
    std::mt19937 gen;

    std::uniform_int_distribution<> dist_age; 
    std::lognormal_distribution<> dist_min_consumption; 
    std::uniform_real_distribution<> dist_consumption;
    std::normal_distribution<> dist_avg_return;
    std::lognormal_distribution<> dist_wage;
    std::normal_distribution<> dist_init_wealth;  // actual dist is lognormal but requires manual processing 
};


#endif // DISTRIBUTIONS_H
