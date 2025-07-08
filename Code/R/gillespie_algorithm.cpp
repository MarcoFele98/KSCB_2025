// Learn more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

#include <Rcpp.h>

using namespace Rcpp; // very bad using namespace Bjarne forgive :(

// [[Rcpp::export]]
NumericMatrix gillespie_simulation(
    IntegerVector &states,
    const double &max_simulation_time,
    const int &number_reagents,
    const int &number_reactions,
    const double &omega,
    const NumericVector &reaction_rates,
    const IntegerMatrix &stoichiometry_reagents,
    const IntegerMatrix &stoichiometry_products
) 
{
  // Checks
  if(number_reagents != stoichiometry_reagents.ncol() | 
     number_reagents != stoichiometry_products.ncol()) {
    stop("Error: incorrect number of reagents");
  }
  if(number_reactions != stoichiometry_reagents.nrow() | 
     number_reactions != stoichiometry_products.nrow()) {
    stop("Error: incorrect number of reactions");
  }
  
  // Initialize simulation 
  states = clone(states); // want to pass by value rather than by reference (which is default behaviour)
  double time = 0;
  NumericMatrix output_long(1000000, number_reagents + 3); // return the number of individuals in each state, the time, the waiting time, the event occurred. I first create a big datatset and then crop based on need... I do not know if this is the best approach rather than dynamically allocating memory (?). I trade-off memory usage (only while simulation is running though) with rapidity
  for(int reagent = 0; reagent < number_reagents; reagent++) {
    output_long.at(0, reagent) = states.at(reagent);
  }
  output_long.at(0, number_reagents + 2) = NA_REAL; // set event occurred to NA. Remember that the first entry is 0.
  
  // Run Gillespie algorithm
  int iteration = 0;
  while(time < max_simulation_time) {
    // Checks
    if(iteration == 1000000) {
      warning("Warning: simulation stopped because reached more than 1'000'000 events");
      break;
    }
    if(sum(states) > 10000) {
      warning("Warning: simulation stopped because species reached more than 10'000 entities");
      break;
    }
      
    // Find propensity of reaction 
    NumericVector propensities(number_reactions);
    for(int reaction_index = 0; reaction_index < number_reactions; reaction_index++) {
      IntegerVector reagents_needed(number_reagents);
      reagents_needed = stoichiometry_reagents(reaction_index, _);
      
      // Check if there are sufficient reagents to make the reaction occur
      IntegerVector remaining_reagents(number_reagents); 
      remaining_reagents = states - reagents_needed;
      if(is_false(all(remaining_reagents >= 0))) {
        //Rprintf("Impossible reaction number %i \n", reaction_index);
        continue;
      } // If the reaction can happen and propensity will be changed to a value bigger than zero in the next section
      
      // Calculate the product of factorial simplifications for every reagent participating in the reaction
      double reaction_product = 1; // just a dummy variable
      for(int reagent_index = 0; reagent_index < number_reagents; reagent_index++) {
        int stoichiometry_reagent = reagents_needed.at(reagent_index);
        if(stoichiometry_reagent != 0) { // if the reagent is present in the reaction
          // Calculate factorial simplification of: states[reagent_index]! / (states[reagent_index] - reagents_needed[reagent_index])!
          double factorial_simplification = 1;
          for(int number_of_reactants = states.at(reagent_index); 
              number_of_reactants > states.at(reagent_index) - stoichiometry_reagent; 
              number_of_reactants--) {
            factorial_simplification *= number_of_reactants;
          }
          reaction_product *= factorial_simplification;
        } // If reagent stoichiometry is zero do not multiply anything for that reagent
      }
      
      propensities.at(reaction_index) = reaction_rates.at(reaction_index) * pow(omega, 1 - sum(reagents_needed)) * reaction_product; 
    }
    if(is_true(all(propensities < 0.00001))) { // stop the simulation if no reaction can occur any more
      warning("Warning: simulation stopped because no reaction could occur anymore");
      break;
    }
    
    // Draw waiting time
    double total_propensities = sum(propensities);
    double waiting_time;
    waiting_time = R::rexp(1 / total_propensities);
    
    // Draw event
    int event;
    NumericVector normalized_propensities(number_reactions);
    normalized_propensities = propensities / total_propensities; 
    event = as<int>(sample(number_reactions, 1, false, normalized_propensities, false));
      
    // Update states
    states = states - stoichiometry_reagents(event, _) + stoichiometry_products(event, _);

    // Update time and save results
    iteration += 1;
    time += waiting_time;
    for(int reagent_index = 0; reagent_index < number_reagents; reagent_index++) {
      output_long.at(iteration, reagent_index) = states.at(reagent_index);
    }
    output_long.at(iteration, number_reagents) = time; 
    output_long.at(iteration, number_reagents + 1) = waiting_time; 
    output_long.at(iteration, number_reagents + 2) = event; 
    }
  
  // Crop the results
  NumericMatrix output(iteration + 1, // iteration starts from zero
                       number_reagents + 3); // return the number of individuals in each state, the time, the waiting time, the event occurred
  for(int column = 0; column < number_reagents + 3; column++) {
    for(int row = 0; row <= iteration; row++) {
      output.at(row, column) = output_long.at(row, column);
    }
  }
  
  // Give names to output
  StringVector names(number_reagents + 3);
  for(int reagent_index = 0; reagent_index < number_reagents; reagent_index++) {
    names.at(reagent_index) = "species_" + std::to_string(reagent_index);
  }
  names.at(number_reagents) = "time";
  names.at(number_reagents + 1) = "waiting_time";
  names.at(number_reagents + 2) = "event";
  
  colnames(output) = names;
  
  return(output);
} 

