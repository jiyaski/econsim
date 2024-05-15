
#include "distributions.h" 

#include <iostream>


class Agent {
public: 
    Agent(double init_wealth, double wage, double avg_return, double consumption); 

    void update(); 
    void print(); 

    double wealth; 
    double wage; 
    double avg_return; 
    double c_param;    // consumption param
}; 



void write_csv(
        const std::vector<std::vector<double>>& data, 
        const std::string& filename
); 