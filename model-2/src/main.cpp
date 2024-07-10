
#include "main.h"

#include <string>
#include <cstring>
#include <algorithm> 
#include <fstream> 
#include <iostream> 
#include <vector>
#include <unistd.h>

RandomSampler randomSampler; 

static constexpr double MAX_DEBT_WAGE_MULTIPLIER = 1.5; 
int num_agents; 
int timestep; 
int lifetime; 
bool debug; 


int main(int argc, char* argv[]) {

    std::vector<Agent> agents; 
    num_agents = 1000;      // defaults if no CLI arg provided 
    timestep = YEAR;        // 
    lifetime = 50;          // 
    debug = false;          // 

    // parse CLI args 
    int opt; 
    while ((opt = getopt(argc, argv, "n:t:l:d")) != -1) {
        switch (opt) {
        case 'n': 
            num_agents = std::atoi(optarg); 
            break; 
        case 't': 
            if       (strcmp(optarg, "day") == 0)    { timestep = DAY;   } 
            else if  (strcmp(optarg, "week") == 0)   { timestep = WEEK;  } 
            else if  (strcmp(optarg, "month") == 0)  { timestep = MONTH; } 
            else if  (strcmp(optarg, "year") == 0)   { timestep = YEAR;  } 
            else {
                std::cerr << "Invalid timestep argument. Use 'day', 'week', 'month', or 'year'.\n";
                return 1;
            }
            break; 
        case 'l': 
            lifetime = std::atoi(optarg); 
            break; 
        case 'd': 
            debug = true; 
            break; 
        default:
            std::cerr << "Usage: " << argv[0] << " [-n num_agents] [-t timestep] [-l lifetime]\n";
            return 1;
        }
    }

    // initialize agents 
    for (int i = 0; i < num_agents; i++) {
        Agent agent;
        agents.push_back(agent); 
    }

    // run simulation and collect results 
    std::vector<double> curr(num_agents); 
    std::vector<std::vector<double>> wealth_data(lifetime, std::vector<double>(num_agents)); 
    for (int t = 0; t < lifetime * timestep; t++) {

        for (int i = 0; i < num_agents; i++) {
            Agent& agent = agents.at(i); 

            // retire & respawn agents who get too old 
            if (agent.age > lifetime) {
                Agent new_agent; 
                new_agent.age = 0; 
                agent = new_agent; 
            }

            curr.at(i) = agent.wealth; 
            agent.update(); 
        }

        // record results 
        if (t % timestep == 0) { wealth_data.at(t / timestep) = curr; } 
        if (debug) { std::cout << "t=" << t << ":\t\t" << agents.at(0).toString() << std::endl; } 
    }
    
    // write out data 
    write_csv(wealth_data, "data/results.csv"); 
    return 0; 
}



Agent::Agent() : Agent(
        randomSampler.generateAge(), 
        randomSampler.generateInitWealth(), 
        randomSampler.generateWage(), 
        randomSampler.generateAnnualROI(), 
        randomSampler.generateConsumption(), 
        randomSampler.generateMinConsumption()
) {}; 


Agent::Agent(double age, double wealth, double wage_param, double annual_ROI, 
                double consumption_param, double min_consumption) : 
        age(age), wealth(wealth), wage_param(wage_param), annual_ROI(annual_ROI), 
        consumption_param(consumption_param), min_consumption(min_consumption) {}; 


void Agent::update() {
    double timedelta = 1 / static_cast<float>(timestep); 
    double annual_wage = wage_param * (0.6667 + 0.01333 * age); 
    wealth = wealth * std::pow(1 + annual_ROI, timedelta) + annual_wage / timestep; 
    double consumption = std::max( 
            min_consumption / timestep, 
            wealth / (timestep * (1 + consumption_param * wealth)) ); 
    double min_allowed_wealth = -MAX_DEBT_WAGE_MULTIPLIER * annual_wage; 

    wealth = std::max(min_allowed_wealth, wealth - consumption); 
    age += timedelta; 
}


std::string Agent::toString() {
    double annual_wage = wage_param * (0.667 + 0.0133 * age); 
    return "age: "          + std::to_string(age) +
        ",\twealth: "       + std::to_string(wealth) +
        ",\twage: "         + std::to_string(annual_wage) +
        ",\tannual ROI: "   + std::to_string(annual_ROI) +
        ",\tc_param: "      + std::to_string(consumption_param);
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