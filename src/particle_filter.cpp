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

/*
 * Convert particle coordinate to map coordinate.
 * @param (xp, yp, heading) x, y coordinates and heading angle of the particle.
 * @param (xc, yc) obserbations in the particle coordinates.
 * @param (c_x, c_y) observations in the map coordinates.
 */
inline void convert_coordinate(double xp,
							   double yp,
							   double heading,
							   double xc,
							   double yc,
							   double &xm,
							   double &ym)
{
	std::cout << "convert_coordinate" << std::endl;
	xm = xp + cos(heading) * xc - sin(heading) * yc;
	ym = yp + sin(heading) * xc + cos(heading) * yc;
}

inline double calculateWeight(double xp,
							  double yp,
							  double xl,
							  double yl,
							  double std_x,
							  double std_y)
{
	std::cout << "calculateWeight" << std::endl;
	double denom = 2.0 * M_PI * std_x * std_y;
	double dx = xp - xl;
	double dy = yp - yl;
	double dxSqrt = dx * dx;
	double dySqrt = dy * dy;
	double exponent = ((dxSqrt / (2.0 * std_x * std_x)) + (dySqrt / (2.0 * std_y * std_y))) * -1.0;
	double numera = exp(exponent);
	return numera/denom;
}

void ParticleFilter::init(double x, double y, double theta, double std[])
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	std::cout << "init" << std::endl;
	num_particles = 100;
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	for (int i = 0; i < num_particles; ++i)
	{
		double sample_x, sample_y, sample_theta;
		sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);
		Particle p;
		p.id = i;
		p.x = sample_x;
		p.y = sample_y;
		p.theta = sample_theta;
		p.weight = 1.0;
		weights.push_back(p.weight);
		particles.push_back(p);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	std::cout << "prediction" << std::endl;
	default_random_engine gen;
	for (auto &p : particles)
	{
		double new_theta = p.theta + yaw_rate * delta_t;
		double cal_x = 0.0;
		double cal_y = 0.0;
		if (yaw_rate == 0) {
			cal_x = p.x + (velocity * delta_t * sin(p.theta));
			cal_y = p.y + (velocity * delta_t * cos(p.theta));
		} else {
			cal_x = p.x + (velocity * (sin(new_theta) - sin(p.theta))) / yaw_rate;
			cal_y = p.y + (velocity * (cos(p.theta) - cos(new_theta))) / yaw_rate;
		}
		normal_distribution<double> dist_x(cal_x, std_pos[0]);
		normal_distribution<double> dist_y(cal_y, std_pos[1]);
		normal_distribution<double> dist_theta(new_theta, std_pos[2]);
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations)
{
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	std::cout << "dataAssociation" << std::endl;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
								   const std::vector<LandmarkObs> &observations, const Map &map_landmarks)
{
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
	std::cout << "updateWeights" << std::endl;
	for (auto& p : particles)
	{
		double x_nearest;
		double y_nearest;
		double minDist = sensor_range * 10000000.0;

		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;

		for (auto obs : observations)
		{
			double xm = 0;
			double ym = 0;
			convert_coordinate(p.x, p.y, p.theta, obs.x, obs.y, xm, ym);

			for (auto lmark : map_landmarks.landmark_list)
			{
				double norm = dist(xm, ym, lmark.x_f, lmark.y_f);
				if (norm < minDist)
				{
					minDist = norm;
					x_nearest = lmark.x_f;
					y_nearest = lmark.y_f;
					associations.clear();
					sense_x.clear();
					sense_y.clear();
					associations.push_back(lmark.id_i);
					sense_x.push_back(lmark.x_f);
					sense_y.push_back(lmark.y_f);
				}
				if (norm == minDist) {
					associations.push_back(lmark.id_i);
					sense_x.push_back(lmark.x_f);
					sense_y.push_back(lmark.y_f);
				}
			}
		}
		if (minDist > sensor_range)
		{
			std::cout << "Too large distance to nearest landmark. " << minDist << endl;
		}
		double weight = calculateWeight(p.x, p.y, x_nearest, y_nearest, std_landmark[0], std_landmark[1]);
		p.weight = weight;
		weights.push_back(weight);
		p.associations = associations;
		p.sense_x = sense_x;
		p.sense_y = sense_y;
		SetAssociations(p, associations, sense_x, sense_y);
	}
}

void ParticleFilter::resample()
{
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::cout << "resample" << std::endl;
	default_random_engine gen;
	discrete_distribution<> dd(weights.begin(), weights.end());
	for (auto& p : particles) {
		double w = weights.at(dd(gen));
		p.weight = w;
	}
}

Particle ParticleFilter::SetAssociations(Particle &particle, const std::vector<int> &associations,
										 const std::vector<double> &sense_x, const std::vector<double> &sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates
	std::cout << "SetAssociations" << std::endl;
	particle.associations = associations;
	particle.sense_x = sense_x;
	particle.sense_y = sense_y;
	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	std::cout << "getAssociations" << std::endl;
	vector<int> v = best.associations;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1); // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1); // get rid of the trailing space
	return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	std::cout << "getSenseY" << std::endl;
	vector<double> v = best.sense_y;
	stringstream ss;
	copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
	string s = ss.str();
	s = s.substr(0, s.length() - 1); // get rid of the trailing space
	return s;
}
