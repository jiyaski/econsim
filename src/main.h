
#include "distributions.h" 

#include <iostream>
#include <string>


class Agent {
public: 

    Agent(); 
    Agent(double age, double wealth, double wage, double annual_ROI, 
            double consumption_param, double min_consumption); 

    void update(); 
    std::string toString(); 

    double age; 
    double wealth; 
    double wage_param; 
    double annual_ROI; 
    double consumption_param; 
    double min_consumption; 
}; 


enum Timestep {
    DAY = 365, 
    WEEK = 52, 
    MONTH = 12, 
    YEAR = 1
}; 



void write_csv(
        const std::vector<std::vector<double>>& data, 
        const std::string& filename
); 