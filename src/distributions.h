
#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <random>


class RandomSampler {
public:
    RandomSampler();  // Constructor 

    double generateConsumption(); 
    double generateAvgReturn();
    double generateWage();
    double generateWealth();

private:
    // const values are generally empirically determined by separate data analysis 
    static constexpr double CONSUMPTION_MIN = 0.0000006;
    static constexpr double CONSUMPTION_MAX = 0.000006;
    static constexpr double ROI_MEAN = 0.06;
    static constexpr double ROI_STDEV = 0.03;
    static constexpr double WAGE_MU = 10.854;
    static constexpr double WAGE_SIGMA = 0.7564;
    static constexpr double INIT_WEALTH_MU = 10.412; 
    static constexpr double INIT_WEALTH_SIGMA = 1.0273; 

    std::random_device rd;
    std::mt19937 gen;

    std::uniform_real_distribution<> dist_consumption;
    std::normal_distribution<> dist_avg_return;
    std::lognormal_distribution<> dist_wage;
    std::lognormal_distribution<> dist_wealth;
};


#endif // DISTRIBUTIONS_H
