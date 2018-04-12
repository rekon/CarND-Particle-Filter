/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	default_random_engine gen;
	double std_x = std[0], std_y = std[1], std_theta = std[2];

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	num_particles = 100;
	for( int i=0; i<num_particles; i++ ){
		Particle p;

		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1;

		particles.push_back(p);
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	double std_x = std[0], std_y = std[1], std_theta = std[2];

	normal_distribution<double> dist_x(0, std_x);
	normal_distribution<double> dist_y(0, std_y);
	normal_distribution<double> dist_theta(0, std_theta);

	for( int i=0; i<num_particles; i++ ){
		double noise_x = dist_x(gen);
		double noise_y = dist_y(gen);
		double noise_theta = dist_theta(gen);
		Particle p = particles[i];
		if( yaw_rate > 0.0001){
			p.x += (velocity / yaw_rate) * (sin(p.theta + yaw_rate*delta_t) - sin(p.theta)) + dist_x(gen);
			p.y += (velocity / yaw_rate) * (-cos(p.theta + yaw_rate*delta_t) + cos(p.theta)) + dist_y(gen);
		} else {
			p.x += (velocity  * (cos(p.theta) * delta_t) + dist_x(gen);
			p.y += (velocity  * (sin(p.theta) * delta_t) + dist_y(gen);
		
		}
		p.theta += yaw_rate*delta_t + dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	int n_obs = observations.size(), n_pdt = predicted.size(), index;
	double min_dist, current_dist;

	for( int i=0; i<n_obs; i++ ){
		min_dist = std::numeric_limits<double>::infinity();
		index = 0;
		for( int j=0; j<n_pdt; j++){
			current_dist = dist( observations[i].x, observations[i].y,predicted[j].x, predicted[j].y);
			if( min_dist > current_dist ){
				min_dist = current_dist;
				index = j;
			}
		}
		observations[i] = predicted[ index ];
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	std::vector<LandmarkObs> observations_in_world_coordinates, owc_copy;
	Particle p;

	for (int k=0; k<num_particles; k++ ){
		observations_in_world_coordinates.clear();
		owc_copy.clear();
		p=particles[k];

		for( int i=0; i<observations.size(); i++){
			LandmarkObs obs;
			obs.id = i
			obs.x = p.x + cos(p.theta) * observations[i].x - sin(p.theta) * observations[i].y;
			obs.y = p.y + sin(p.theta) * observations[i].x + cos(p.theta) * observations[i].y;

			observations_in_world_coordinates.push_back( obs );
			owc_copy.push_back( obs );
		}

		dataAssociation( map_landmarks, observations_in_world_coordinates );

		p.weight = 1;
		for( int i=0; i<owc_copy.size(); i++ ){
			p.weight *= multivariate_gauss( owc_copy[i].x, observations_in_world_coordinates[i].x, std_landmark[0],
											owc_copy[i].y, observations_in_world_coordinates[i].y, std_landmark[1] );
		}

	}
	

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
