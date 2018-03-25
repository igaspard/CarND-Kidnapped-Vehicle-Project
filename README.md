# Kidnapped Vehicle Project - Gaspard Shen
In this Kidnapped Vehicle project, I implement the particle filter to demonstrate how the localization problem can be handled by this concept. By using the GPS data to initialize the car location and predicted new location. In the end, got result pass the criteria and the errors are **x.116, y.108 yaw.004** by using **1000** particles and runtime is **55.82 sec**.

## Implementation

![](/output/FlowChart.png)

1. Initialization
  - At first, set the number of particles. Initialize all particles to first position based on estimates of x, y, theta and their uncertainties from GPS and all weights to 1. Add random Gaussian noise to each particle. In `ParticleFilter::init()`.
2. Prediction
  - Add measurements to each particle and add random Gaussian noise in `ParticleFilter::prediction()`.
3. Update Step (Update Weights)
  - Here we update the weights of each particle using a **mult-variate Gaussian distribution**.
  - The observations are given in the VEHICLE'S coordinate system and our particles are located according to the MAP'S coordinate system. We need to transform between the two systems based on the **Homogenous Transformation**.
  - Based on the Map Landmark positions in the sensor range as given, associate the observation with implement the helper function `ParticleFilter::dataAssociation()` to find out each the closest predicted measurement to particular landmark.
  - Finish the calculation of Multivariate-Gaussian probability density in the `ParticleFilter::updateWeights()`.
4. Resample
  - `ParticleFilter::resample()` resample particles with replacement with probability proportional to their weight.

## Results

In the beginning, I use the 100 particles and pass the criteria as error **x.164, y.148, yaw.005**. Below I find the out increase the number of particles from 100, 500, 1000 and up to 2000. Look like when the number of particles over 1000, the error can't lower than **x.116, y.108 yaw.004** and runtime is **55.82 sec**. When 2000 particles, the runtime become 76.16 sec but the error is similar to 1000. So in the end i choose 1000 as final numbers.
**Number of Particles 1000**
![](/output/Result_NumOfParticles1000.png)
**Number of Particles 100**
![](/output/Result_NumOfParticles100.png)
**Number of Particles 500**
![](/output/Result_NumOfParticles500.png)
**Number of Particles 2000**
![](/output/Result_NumOfParticles2000.png)
