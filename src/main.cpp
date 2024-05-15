
#include "main.h"

#include <string>
#include <fstream> 
#include <iostream> 
#include <vector>
#include <unistd.h>


int main(int argc, char* argv[]) {

    RandomSampler randomSampler; 
    std::vector<Agent> agents; 
    int num_agents = 100;   // defaults if no CLI arg provided 
    int max_timestep = 50;  // 

    // parse CLI args 
    int opt; 
    while ((opt = getopt(argc, argv, "n:t:")) != -1) {
        switch (opt) {
        case 'n': 
            num_agents = std::atoi(optarg); 
            break; 
        case 't': 
            max_timestep = std::atoi(optarg); 
            break; 
        }
    }

    // initialize agents 
    for (int i = 0; i < num_agents; i++) {
        Agent agent(
            randomSampler.generateWealth(), 
            randomSampler.generateWage(), 
            randomSampler.generateAvgReturn(), 
            randomSampler.generateConsumption()
        ); 
        agents.push_back(agent); 
    }

    // run simulation and collect results 
    std::vector<std::vector<double>> wealth_data(max_timestep, std::vector<double>(num_agents)); 
    for (int t = 0; t < max_timestep; t++) {
        for (int i = 0; i < num_agents; i++) {
            wealth_data.at(t).at(i) = agents.at(i).wealth; 
            agents.at(i).update(); 
        }
    }
    
    // write out data 
    write_csv(wealth_data, "temp.csv"); 
    


    return 0; 
}





Agent::Agent(double init_wealth, double wage, double avg_return, double c_param) : 
        wealth(init_wealth), wage(wage), avg_return(avg_return), c_param(c_param) {}; 

void Agent::update() {
    // todo: make documentation for why this is the update rule 
    wealth = wealth * (1 - (1 / (1 + c_param * wealth))) * (1 + avg_return) + wage; 
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