/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 50;
  std::default_random_engine gen;
  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

	for(int i = 0; i < num_particles; i++){
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);		
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;

		particles.push_back(particle);
		weights.push_back(1.0);
	}

	is_initialized = true;

  }

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
	normal_distribution<double> std_pos_x(0 , std_pos[0]);
	normal_distribution<double> std_pos_y(0 , std_pos[1]);
	normal_distribution<double> std_pos_theta(0 , std_pos[2]);
  
    std::default_random_engine gen;
  
	// perform a movement step for all particles	
	for(int i = 0; i < num_particles; i++){

		if(fabs(yaw_rate) == 0){
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		} 
		else{
			particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + (yaw_rate*delta_t)) - sin(particles[i].theta)); 
			particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate*delta_t)));
			particles[i].theta += yaw_rate * delta_t; 
		}

		// add gaussian noise to all the particles
		particles[i].x += std_pos_x(gen);
		particles[i].y += std_pos_y(gen);
		particles[i].theta += std_pos_theta(gen);
	}	
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations, Particle &particle) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
    std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;

	// For each observation...
	for(unsigned j = 0; j < observations.size(); ++j) {
		int index, landmark_id;
		double min_dist = 1.0e99;

		// ...find the minimum distance between current observation and
		// landmarks in sensor range only
		double x1 = observations[j].x;
		double y1 = observations[j].y;

		for(unsigned k = 0; k < predicted.size(); ++k) {
			double x2 = predicted[k].x;
			double y2 = predicted[k].y;

			double current_dist = dist(x1, y1, x2, y2);

			if(current_dist < min_dist) {
				min_dist = current_dist;
				index = k;
				landmark_id = predicted[k].id;
			}
		}

		// Assign the nearest landmark to current observation
		observations[j].id = index;

		associations.push_back(landmark_id);
		sense_x.push_back(x1);
		sense_y.push_back(y1);
	}

	// Set assocations for visualization purpose only
//	particle = SetAssociations(particle, associations, sense_x, sense_y);
	
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
// For each particle do following
	for(int i = 0; i < num_particles ; i++){ 
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;

		std::vector<LandmarkObs> map_observations; 
		for(unsigned int j = 0; j < observations.size(); ++j){
		  //  cout<<"x: "<<observations[j].x<<"\ty: "<<observations[j].y<<"Before id"<<observations[j].id<<endl;
			double t_x = p_x + (cos(p_theta)*observations[j].x) - (sin(p_theta)*observations[j].y);
			double t_y = p_y + (sin(p_theta)*observations[j].x) + (cos(p_theta)*observations[j].y);
			int id = observations[j].id;			
			map_observations.push_back(LandmarkObs{ id, t_x, t_y });
		}
		
		std::vector<LandmarkObs> predicted;
		// Get list of landmarks in the sensor range from the particle
		for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j){
			double l_x = map_landmarks.landmark_list[j].x_f;
			double l_y = map_landmarks.landmark_list[j].y_f;
			int id = map_landmarks.landmark_list[j].id_i;
			if (dist(p_x,p_y,l_x,l_y) <= sensor_range){
				predicted.push_back(LandmarkObs{id, l_x, l_y});				
			}
		}
		
		// Associate observations with nearest landmarks
		dataAssociation(predicted, map_observations, particles[i]);
		// Reinitialize weight
		particles[i].weight = 1.0;	

		// Calculate weight of the particle using multi-variate gaussian
		for(int j = 0; j < map_observations.size(); ++j){
			double map_x = map_observations[j].x;
			double map_y = map_observations[j].y;
			//cout<<"x: "<<observations[j].x<<"\ty: "<<observations[j].y<<"After id"<<observations[j].id<<endl;
			int lm_id = map_observations[j].id;
			/*double true_lm_x , true_lm_y;
			for(int k = 0; k < predicted.size(); k++){			
				if (predicted[k].id == lm_id){
					true_lm_x = predicted[k].x;
					true_lm_y = predicted[k].y;					
					break;
					
				}
			}	*/
			double ldm_x = predicted[lm_id].x;
			double ldm_y = predicted[lm_id].y;
			double std_x = std_landmark[0];
			double std_y = std_landmark[1];

			double weight = ( 1/(2*M_PI*std_x*std_y)) * exp( -( (pow(map_x-ldm_x,2)/(2*pow(std_x, 2))) + ((pow(map_y-ldm_y,2)/(2*pow(std_y, 2))))));
			if(weight>0){
				particles[i].weight *= weight;	
			}			
		}
		weights[i]=(particles[i].weight);
	}
	

  	// Calculate sum of all weights for normalization
  /*	double sum_of_weights = 0.0;
  	for(int j = 0; j < num_particles; j++){
		sum_of_weights += particles[j].weight;
  	}

  	// Normalize weights
  	for(int j = 0; j < num_particles; j++){
  		particles[j].weight /= sum_of_weights ;
		weights[j] = particles[j].weight;
  	}*/
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
	// Resample particles with replacement with probability proportional to their weight
	default_random_engine gen;
	std::discrete_distribution<int> d(weights.begin(), weights.end());
	std::vector<Particle> new_particles;

	for(int i = 0; i < num_particles; i++) {
		new_particles.push_back(particles[d(gen)]);
	}

	particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}