
#include "main.h"

#include <string>
#include <cstring>
#include <fstream> 
#include <iostream> 
#include <vector>
#include <unistd.h>


int main(int argc, char* argv[]) {

    RandomSampler randomSampler; 
    std::vector<Agent> agents; 
    int num_agents = 100;   // defaults if no CLI arg provided 
    int timestep = YEAR;    // 
    int lifetime = 50;      // 

    // parse CLI args 
    int opt; 
    while ((opt = getopt(argc, argv, "n:t:l:")) != -1) {
        switch (opt) {
        case 'n': 
            num_agents = std::atoi(optarg); 
            break; 
        case 't': 
            if (strcmp(optarg, "day") == 0) {
                timestep = DAY;
            } else if (strcmp(optarg, "week") == 0) {
                timestep = WEEK;
            } else if (strcmp(optarg, "month") == 0) {
                timestep = MONTH;
            } else if (strcmp(optarg, "year") == 0) {
                timestep = YEAR;
            } else {
                std::cerr << "Invalid timestep argument. Use 'day', 'week', 'month', or 'year'.\n";
                return 1;
            }
        case 'l': 
            lifetime = std::atoi(optarg); 
            break; 
        default:
            std::cerr << "Usage: " << argv[0] << " [-n num_agents] [-t timestep] [-l lifetime]\n";
            return 1;
        }
    }

    // initialize agents 
    for (int i = 0; i < num_agents; i++) {
        Agent agent(
            randomSampler.generateAge(), 
            randomSampler.generateInitWealth(), 
            randomSampler.generateWage(), 
            randomSampler.generateAnnualROI(), 
            randomSampler.generateConsumption(), 
            randomSampler.generateMinConsumption()
        ); 
        agents.push_back(agent); 
    }

    // run simulation and collect results 
    std::vector<std::vector<double>> wealth_data(lifetime, std::vector<double>(num_agents)); 
    for (int t = 0; t < lifetime; t++) {
        for (int i = 0; i < num_agents; i++) {
            wealth_data.at(t).at(i) = agents.at(i).wealth; 
            agents.at(i).update(); 
        }
    }
    
    // write out data 
    write_csv(wealth_data, "temp.csv"); 
    


    return 0; 
}




Agent::Agent(int age, double wealth, double wage_param, double annual_ROI, 
                double consumption_param, double min_consumption) : 
        age(age), wealth(wealth), wage_param(wage_param), annual_ROI(annual_ROI), 
        consumption_param(consumption_param), min_consumption(min_consumption) {}; 



void Agent::update() {
    double temp = wealth * std::pow(1 + annual_ROI, 1 / )
}

void Agent::print() {
    std::cout << "Agent:\t wealth: " << wealth << ",\twage: " << wage << ",\tavg_return: " 
    << avg_return << ",\tc_param: " << c_param << std::endl; 
}



void write_csv(const std::vector<std::vector<double>>& data, const std::string& filename) {

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) // Check if delimiter is needed
                file << "\t";
        }
        file << "\n"; // End of line for each row
    }

    file.close();
    if (!file) { std::cerr << "An error occurred during file write to: " << filename << std::endl; } 
    else { std::cout << "File written successfully: " << filename << std::endl; }
}