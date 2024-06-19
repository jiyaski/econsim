
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
    // const values are generally empirically determined by separate data analysis, then tweaked 
    static constexpr double AGE_MIN = 0; 
    static constexpr double AGE_MAX = 50; 
    static constexpr double MIN_CONSUMPTION_MU = 9.7860;          // 9.7860 
    static constexpr double MIN_CONSUMPTION_SIGMA = 0.46052;      // 0.46052
    static constexpr double CONSUMPTION_MIN = 0.0000005;         // 0.0000005
    static constexpr double CONSUMPTION_MAX = 0.000005;           // 0.000005
    static constexpr double ROI_MEAN = 0.06;                     // 0.06
    static constexpr double ROI_STDEV = 0.03;                     // 0.03
    static constexpr double WAGE_MU = 10.954;                      // 10.954
    static constexpr double WAGE_SIGMA = 0.7564;                  // 0.7564
    static constexpr double INIT_WEALTH_MU = 10.412;              // 10.412 
    static constexpr double INIT_WEALTH_SIGMA = 1.0273;             // 1.0273 
    static constexpr double INIT_WEALTH_10TH_PERCENTILE = -11270; 

    std::random_device rd;
    std::mt19937 gen;

    std::uniform_real_distribution<> dist_age; 
    std::lognormal_distribution<> dist_min_consumption; 
    std::uniform_real_distribution<> dist_consumption;
    std::normal_distribution<> dist_annual_ROI;
    std::lognormal_distribution<> dist_wage;
    std::normal_distribution<> dist_init_wealth;  // actual dist is lognormal but requires manual processing 
};


#endif // DISTRIBUTIONS_H
