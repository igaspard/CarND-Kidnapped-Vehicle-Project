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
	num_particles = 1000;

	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i = 0; i < num_particles; ++i) {
		Particle particle;
		particle.id 		= i;
		particle.x 			= dist_x(gen);
		particle.y	 		= dist_y(gen);
		particle.theta 	= dist_theta(gen);
		particle.weight	= 1.0;

		particles.push_back(particle);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	normal_distribution<double> n_x(0, std_pos[0]);
	normal_distribution<double> n_y(0, std_pos[1]);
	normal_distribution<double> n_theta(0, std_pos[2]);

	for (auto& p: particles) {
		if (fabs(yaw_rate) > 0.001) {
				p.x += velocity / yaw_rate * ( sin(p.theta + yaw_rate*delta_t) - sin(p.theta) );
				p.y += velocity / yaw_rate * ( cos(p.theta) - cos(p.theta + yaw_rate*delta_t) );

				p.theta += yaw_rate * delta_t;
		}
		else {
				p.x += velocity * delta_t * cos(p.theta);
				p.y += velocity * delta_t * sin(p.theta);
		}
		p.x += n_x(gen);
		p.y += n_y(gen);

		p.theta += n_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	for (auto& obs: observations) {
		double min_dist = numeric_limits<double>::max();
		for (auto& pred: predicted) {
			double cur_dist = dist(obs.x, obs.y, pred.x, pred.y);
			if (cur_dist < min_dist) {
				min_dist = cur_dist;
				obs.id = pred.id;
			}
		}
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
	const double gauss_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
	const double std_landmark_x2 = std_landmark[0] * std_landmark[0];
	const double std_landmark_y2 = std_landmark[1] * std_landmark[1];
	for (auto& p: particles) {
		// Homogenous Transformation from vehice coordinates to map coordinates
		vector<LandmarkObs> trans_obs;
		for (auto& obs: observations) {
			LandmarkObs tobs;
			tobs.x 	= p.x + cos(p.theta)*obs.x - sin(p.theta)*obs.y;
			tobs.y 	= p.y + sin(p.theta)*obs.x + cos(p.theta)*obs.y;
			trans_obs.push_back(tobs);
		}
		// Collect the map landmark predictions
		vector<LandmarkObs> predictions;
		for (auto& lm: map_landmarks.landmark_list) {
			LandmarkObs tobs;
			tobs.x 	= lm.x_f;
			tobs.y 	= lm.y_f;
			tobs.id	= lm.id_i;
			double distance = dist(p.x, p.y, tobs.x, tobs.y);
			if (distance < sensor_range) {
				predictions.push_back(tobs);
			}
		}
    dataAssociation(predictions, trans_obs);
		// Multivariate-Gaussian probability density
    p.weight = 1.0;
		for (auto& tobs: trans_obs) {
			double x_obs = tobs.x;
			double y_obs = tobs.y;

			int asso_prediction = tobs.id;
			double mu_x, mu_y;
			for (auto& pred: predictions) {
				if (pred.id == asso_prediction) {
					mu_x = pred.x;
					mu_y = pred.y;
					break;
				}
			}
			double x_diff = x_obs - mu_x;
			double y_diff = y_obs - mu_y;
			double x_diff_2 = x_diff*x_diff;
			double y_diff_2 = y_diff*y_diff;

			double multi = gauss_norm * exp( -(x_diff_2/2*std_landmark_x2 + y_diff_2/2*std_landmark_y2) );
			p.weight *= multi;
		}
    weights.push_back(p.weight);
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	discrete_distribution<int> distribution(weights.begin(), weights.end());

	vector<Particle> resample_particles;
	for (int i = 0; i < num_particles; ++i) {
		resample_particles.push_back(particles[distribution(gen)]);
	}
	particles = resample_particles;

	weights.clear();
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
