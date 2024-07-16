
#include "pch.h"
#include "main.h" 

#include <cstring>
#include <iostream> 
#include <fstream>
#include <vector> 
#include <unistd.h>
#include <memory>


namespace po = boost::program_options;
SimConfig sim;

Eigen::VectorXd quants(0);   // quantities 
Eigen::VectorXd prices(0); 
Eigen::VectorXd quants_0(0); 
Eigen::VectorXd prices_0(0); 
Eigen::MatrixXd Elas(0, 0);  // elasticity matrix 
Eigen::MatrixXd Elas_Inv(0, 0); 
Eigen::Matrix<double,Eigen::Dynamic,4> Prod_Costs(0, 4); 
int N;  // number of firms 

std::vector<std::unique_ptr<Firm>> firms; 



int main(int argc, char* argv[]) {

    parse_options(argc, argv);
    std::cout << "\n" << sim.to_string() << "\n\n"; 

    std::random_device rd;
    std::mt19937 gen(rd()); 
    N = 0; 

    for (int i = 0; i < sim.num_firms; i++) {
        add_good(gen); 
    }

    std::cout << "Elasticities: --------------------------------\n" << Elas 
            << "\n\nquants_0: " << quants_0.transpose() 
            << "\n\nprices_0: " << prices_0.transpose() 
            << "\n\nProd_Costs: \n" << Prod_Costs << "\n\n"; 
    
    for (int i = 0; i < 10; i++) {
        step(); 
        std::cout << "profits (t=" << i << "): \t" << profits_as_str() << "\n"; 
    }
    std::cout << std::endl; 


    return 0;
}


void step() {

    for (int i = 0; i < N; i++) {

        // remove any firms that have been unprofitable long enough to exit 
        if (firms.at(i)->exit_counter > 3 * sim.timestep) {
            remove_good(i); 
            continue; 
        }

        // find the optimal production quantity. a0, a1, a2 are params of the quadratic equation MC = MR 
        double opt_quant; 
        double marginal_revenue = Elas_Inv.row(i).dot(quants_0 - quants) - quants(i) / Elas(i, i); 
        double a2 = 3 * Prod_Costs(i, 3); 
        double a1 = 2 * Prod_Costs(i, 2); 
        double a0 = Prod_Costs(i, 1) - marginal_revenue; 
        double discriminant = a1 * a1 - 4 * a2 * a0; 

        // no real roots -> MC > MR always -> shutdown production; otherwise use rightmost root 
        if (discriminant <= 0) { opt_quant = 0; } 
        else { opt_quant = (-a1 + sqrt(discriminant)) / (2 * a2); } 

        // calculate stuff based on optimal quantity, to make shutdown/exit decisions 
        double opt_cost = calculate_cost(i, opt_quant); 
        double opt_avg_total_cost = (opt_quant != 0) ? (opt_cost / opt_quant) : 0; 
        double opt_avg_var_cost = (opt_quant != 0) ? ((opt_cost - Prod_Costs(i, 0)) / opt_quant) : 0; 
        Eigen::VectorXd temp_quants = quants; 
        temp_quants(i) = opt_quant; 
        double opt_price = Elas_Inv.row(i).dot(quants_0 - temp_quants); 

        // shutdown/exit logic 
        if (opt_price < opt_avg_var_cost) { opt_quant = 0; } 
        if (opt_price < opt_avg_total_cost) { firms.at(i)->exit_counter++; } 

        // calculate actual next quantity; recalculate other stuff based on it 
        quants(i) = std::min(opt_quant, quants(i) + quants_0(i) / 365 / 2); 
        prices(i) = Elas_Inv.row(i).dot(quants_0 - quants); 
        double cost = calculate_cost(i, quants(i)); 
        double avg_total_cost = (quants(i) != 0) ? (cost / quants(i)) : 0; 
        double avg_var_cost = (quants(i) != 0) ? ((cost - Prod_Costs(i, 0)) / quants(i)) : 0; 
        
        // calculate relevant data to store 
        double profit = (prices(i) - avg_total_cost) * quants(i); 
        // todo: calculate consumer surplus & deadweight loss (need numerical integration) 

        // update firm 
        firms.at(i)->profit += profit; 
    }
}


/**
 * cost of producing a given quantity of the `i`th good 
 * 
 * @param ind - index of the good's row in the Prod_Costs matrix 
 * @param quantity - quantity of the good that would be produced 
 */
double calculate_cost(int ind, double quantity) {
    double q_term = 1; 
    double result = 0; 
    for (int i = 0; i <= 3; i++) {
        result += Prod_Costs(ind, i) * q_term; 
        q_term *= quantity; 
    }
    return result / sim.timestep; 
}


std::string profits_as_str() {
    std::string result = ""; 
    for (const auto& firm : firms) {
        result += std::to_string(firm->profit) + " \t"; 
    }
    return result; 
}


/**
 * adds a good to the market. 
 * 
 * @param gen - pseudorandom generator to use for initializing the good's properties 
 */
void add_good(std::mt19937& gen) {

    // a + rand*(b-a) transforms the uniform distribution to produce rand outputs between a and b 
    std::uniform_real_distribution<> dist(0, 1); 

    // generate new entry for quants_0 and prices_0 vectors 
    double quant_0, price_0, MPS; 
    if (sim.uniform_demand && N >= 1) {
        quant_0 = quants_0(0); 
        price_0 = prices_0(0); 
        MPS = quant_0 * price_0 / 2; 
    } else {
        MPS = std::pow(10, LOG_MPS_MIN + dist(gen) * (LOG_MPS_MAX - LOG_MPS_MIN)) / sim.timestep; 
        double r = std::pow(10, LOG_R_MIN + dist(gen) * (LOG_R_MAX - LOG_R_MIN)) / sim.timestep; 
        quant_0 = 2 * r; 
        price_0 = MPS / r; 
    }
    quants_0.conservativeResize(N + 1); 
    prices_0.conservativeResize(N + 1); 
    quants_0(N) = quant_0; 
    prices_0(N) = price_0; 


    //// add new row/col to Elasticity matrix & recalculate its inverse 
    double ela_ii = sim.flat_demand ? 0 : quant_0 / price_0; 
    Elas.conservativeResize(N+1, N+1); 
    Elas(N, N) = ela_ii; 

    for (int i = 0; i < N; i++) { 
        if (sim.flat_demand || sim.compl_sign == "zero") { 
            Elas(N, i) = 0; 
            Elas(i, N) = 0;
            continue; 
        }

        // if sim.compl_sign = "pos", we'll just leave these values as they are
        Elas(N, i) =  (sim.compl_min + dist(gen) * (sim.compl_max - sim.compl_min)) * ela_ii; 
        Elas(i, N) = ((sim.compl_min + dist(gen) * (sim.compl_max - sim.compl_min)) * Elas(i, i)); 

        if (sim.compl_sign == "neg") {
            Elas(N, i) = -Elas(N, i); 
            Elas(i, N) = -Elas(i, N); 
        } 
        else if (sim.compl_sign == "any") {
            double rand = dist(gen); 
            if (rand < 0.25)       { /* do nothing */ } 
            else if (rand < 0.5)   { Elas(N, i) = -Elas(N, i); } 
            else if (rand < 0.75)  { Elas(i, N) = -Elas(i, N); } 
            else                   { Elas(N, i) = -Elas(N, i);    Elas(i, N) = -Elas(i, N); }
        }
    }
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(Elas); 
    Elas_Inv = cod.pseudoInverse(); 

    // update quantity & price vectors
    quants.conservativeResize(N + 1); 
    prices.conservativeResize(N + 1); 
    quants(N) = 0; 
    prices(N) = Elas_Inv.row(N).dot(quants_0 - quants);

    // generate new production cost parameters (c0->fixed costs; c1->linear costs; c2->quadratic; c3->cubic) 
    if (sim.uniform_costs && N >= 1) {
        Prod_Costs.conservativeResize(N + 1, Eigen::NoChange); 
        Prod_Costs.row(N) = Prod_Costs.row(0); 
    } else {
        double c0 = dist(gen) * MPS / 5; 
        double c3 = dist(gen) * 16 / 3 * price_0 / (quant_0*quant_0); 
        std::vector<double> intervals = {0.0, quant_0}; 
        std::vector<double> weights = {1.0, 0.0};  // max probability at 0, decreasing to zero at quant_0 
        std::piecewise_linear_distribution<> dist2(intervals.begin(), intervals.end(), weights.begin()); 

        double vertex_x = dist2(gen); 
        Eigen::VectorXd quants_temp = quants; 
        quants_temp(N) = vertex_x; 
        double vertex_y_max = Elas_Inv.row(N).dot(quants_0 - quants_temp); 
        double vertex_y = dist(gen) * vertex_y_max; 

        double c2 = -3 * c3 * vertex_x; 
        double c1 = vertex_y + 3 * c3 * vertex_x  * vertex_x; 
        Prod_Costs.conservativeResize(N + 1, Eigen::NoChange); 
        Prod_Costs.row(N) << c0 / sim.timestep, c1 / sim.timestep, c2 / sim.timestep, c3 / sim.timestep; 
    }

    firms.push_back(std::make_unique<Firm>(Firm{})); 
    N++; 
}


/**
 * removes a good from the market. 
 * 
 * @param k - index of the good to remove
 */
void remove_good(int k) {

    Eigen::VectorXd quants_temp(N - 1); 
    Eigen::VectorXd prices_temp(N - 1); 
    Eigen::VectorXd quants_0_temp(N - 1); 
    Eigen::VectorXd prices_0_temp(N - 1); 
    Eigen::MatrixXd Elas_temp_1(N, N - 1); 
    Eigen::MatrixXd Elas_temp_2(N - 1, N - 1); 
    Eigen::Matrix<double, Eigen::Dynamic, 4> Prod_Costs_temp(N - 1, 4); 

    quants_temp << quants.head(k), quants.tail(N - k - 1); 
    prices_temp << prices.head(k), prices.tail(N - k - 1); 
    quants_0_temp << quants_0.head(k), quants_0.tail(N - k - 1); 
    prices_0_temp << prices_0.head(k), prices_0.tail(N - k - 1); 
    Elas_temp_1 << Elas.leftCols(k), Elas.rightCols(N - k - 1); 
    Elas_temp_2 << Elas_temp_1.topRows(k), Elas_temp_1.bottomRows(N - k - 1); 
    Prod_Costs_temp << Prod_Costs.topRows(k), Prod_Costs.bottomRows(N - k - 1); 

    quants = quants_temp; 
    prices = prices_temp; 
    quants_0 = quants_0_temp; 
    prices_0 = prices_0_temp; 
    Elas = Elas_temp_2; 
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(Elas); 
    Elas_Inv = cod.pseudoInverse(); 
    Prod_Costs = Prod_Costs_temp; 

    firms.erase(firms.begin() + k); 
    N--; 
}


void parse_options(int argc, char* argv[]) {
    
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message") 
        ("config,c", po::value<std::string>(), "(str) config file to use") 
        ("timestep,t", po::value<std::string>(), "(str) timestep (day, week, month, or year)") 
        ("runtime,T", po::value<int>(), "(int) how long the sim should run, in years") 
        ("record,r", "write sim data to output files") 
        ("rec_step,R", po::value<std::string>(), "(int) output data every ___ timesteps") 
        ("log,l", "output debug and other log data") 
        ("warmup,w", po::value<int>(), "(int) number of warmup timesteps (no data is collected for these)")

        ("num_firms,f", po::value<int>(), "(int) number of firms in the simulation")
        ("ease_of_entry,e", po::value<float>(), "(float) ease of entry for new firms")
        ("compl_sign,s", po::value<std::string>(), "(string) possible signs for complementarity values (pos, neg, zero, or any)")
        ("compl_min,m", po::value<float>(), "(float) min complementarity magnitude")
        ("compl_max,M", po::value<float>(), "(float) max complementarity magnitude")
        ("flat_demand,d", "forces elasticity matrix to be zero, resulting in constant demand curve")
        ("uniform_demand,u", "q_0, p_0 entries for all firms will be same")
        ("uniform_costs,U", "all firms have same production cost parameters")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        exit(1);
    }

    // Read options from the config file if provided.
    if (vm.count("config")) {
        std::string config_file = "config/" + vm["config"].as<std::string>() + ".cfg"; 
        std::ifstream ifs(config_file.c_str()); 
        if (ifs) {
            po::store(po::parse_config_file(ifs, desc), vm);
        } else {
            std::cerr << "Cannot open config file: " << config_file << "\n";
            exit(1);
        }
    }

    po::notify(vm);

    if (vm.count("timestep")) {
        std::string timestep_str = vm["timestep"].as<std::string>(); 
        if       (timestep_str == "day"  )  { sim.timestep = Timestep::DAY;   } 
        else if  (timestep_str == "week" )  { sim.timestep = Timestep::WEEK;  } 
        else if  (timestep_str == "month")  { sim.timestep = Timestep::MONTH; } 
        else if  (timestep_str == "year" )  { sim.timestep = Timestep::YEAR;  } 
        else  { 
            std::cerr << "Invalid timestep: " << timestep_str << "\n"; 
            exit(1); 
        }
    }

    if  (vm.count("runtime"))  { sim.runtime = vm["runtime"].as<int>();         }
    if  (vm.count("record"))   { sim.record = true;                             }
    if  (vm.count("rec_step")) { sim.recording_step = vm["rec_step"].as<int>(); }
    if  (vm.count("log"))      { sim.log = true;                                }
    if  (vm.count("warmup"))   { sim.warmup = vm["warmup"].as<int>();           }

    if  (vm.count("num_firms"))      { sim.num_firms = vm["num_firms"].as<int>();           }
    if  (vm.count("ease_of_entry"))  { sim.ease_of_entry = vm["ease_of_entry"].as<float>(); }
    if  (vm.count("compl_sign"))     { sim.compl_sign = vm["compl_sign"].as<std::string>(); } 
    if  (vm.count("compl_min"))      { sim.compl_min = vm["compl_min"].as<float>();         }
    if  (vm.count("compl_max"))      { sim.compl_max = vm["compl_max"].as<float>();         }
    if  (vm.count("flat_demand"))    { sim.flat_demand = true;                              }
    if  (vm.count("uniform_demand")) { sim.uniform_demand = true;                           }
    if  (vm.count("uniform_costs"))  { sim.uniform_costs = true;                            }
}


std::string SimConfig::to_string() const {
    return "Simulation Config Parameters: \n=============================" 
            + std::string("\nnum_firms: ") + std::to_string(num_firms) 
            + std::string("\nease_of_entry: ") + std::to_string(ease_of_entry) 
            + std::string("\ncompl_sign: ") + compl_sign
            + std::string("\ncompl_min: ") + std::to_string(compl_min) 
            + std::string("\ncompl_max: ") + std::to_string(compl_max) 
            + std::string("\nflat_demand: ") + std::to_string(flat_demand) 
            + std::string("\nuniform_demand: ") + std::to_string(uniform_demand) 
            + std::string("\nuniform_costs: ") + std::to_string(uniform_costs) 
            + std::string("\ntimestep: ") + std::to_string(timestep) 
            + std::string("\nruntime: ") + std::to_string(runtime) 
            + std::string("\nrecord: ") + std::to_string(record) 
            + std::string("\nrecording_step: ") + std::to_string(recording_step) 
            + std::string("\nlog: ") + std::to_string(log) 
            + std::string("\nwarmup: ") + std::to_string(warmup); 
}



