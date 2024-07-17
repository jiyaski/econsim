
#include <string>
#include <random>



// bounds for the distribution of Maximum Possible Surplus (MPS) and `r` constants for each firm 
static constexpr float LOG_MPS_MIN = 6; 
static constexpr float LOG_MPS_MAX = 8; 
static constexpr float LOG_R_MIN = 2; 
static constexpr float LOG_R_MAX = 6; 


enum Timestep {
    DAY = 365, 
    WEEK = 52, 
    MONTH = 12, 
    YEAR = 1
}; 


struct SimConfig {
    int num_firms = 10; 
    float ease_of_entry = 1; 
    std::string compl_sign = "any"; 
    float compl_min = 0.0f; 
    float compl_max = 0.1f; 
    bool flat_demand = false; 
    bool uniform_demand = false; 
    bool uniform_costs = false; 

    Timestep timestep = Timestep::MONTH; 
    int runtime = 50; 
    bool record = 0; 
    int recording_step = 1; 
    bool log = 0;
    int warmup = 0; 

    std::string to_string() const; 
}; 


struct Firm { 
    double profit = 0; 
    double deadweight_loss = 0; 
    double consumer_surplus = 0; 
    int exit_counter = 0; 
}; 


void step(); 

void add_good(std::mt19937& gen); 
void remove_good(int k); 
std::string profits_as_str(); 
void parse_options(int argc, char* argv[]); 

double get_price(size_t k, double quantity); 
double get_marginal_revenue(size_t k, double quantity); 
double get_cost(size_t k, double quantity); 
double get_marginal_cost(size_t k, double quantity); 
double get_avg_total_cost(size_t k, double quantity); 
double get_avg_var_cost(size_t k, double quantity); 

double get_opt_quantity(size_t k); 
double get_social_opt_quantity(size_t k); 

