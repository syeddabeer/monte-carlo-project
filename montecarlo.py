# Monte Carlo Simulation
# Project code Electrical Power systems Reliability
# Dabeeruddin Syed

# Importing required packages
import numpy
import matplotlib.pyplot as plt
import math
import scipy.stats
SEED = 4

# Class for components that have two possible objects
class TwoState:
    def __init__(self, failure_rate, repair_rate, random_state, name=''):
        self._failure_rate = failure_rate/8760
        self.repair_rate = repair_rate/8760
        self.random_state = random_state
        self.name = name
        # Gives a name to the state (i.e. 1,2,3 ....)

        # state True indicates the component is UP
        self.state = True
        self.ttc = self.time_to_change()

    def __repr__(self):
        return self.name

    @property
    def failure_rate(self):
        return self._failure_rate

    @property
    def ttc(self):
        return self._ttc

    @ttc.setter
    def ttc(self, value):
        self._ttc = value

    def time_to_change(self):
        #  Calculation of the time to change based on the state of component
        if self.state:
            rate = self.failure_rate
        else:
            rate = self.repair_rate

        ttc = time_to_change(randomnum=self.random_state.rand(), rho=rate)

        return ttc

    def change_state_update_ttc(self):
        #  Update Time to Change and invert the state of a component
        self.state = not self.state
        self.ttc = self.time_to_change()

    def update_ttc(self, delta):
        #  Updating the time to change
        #  Basically changing the time to change with delta value
        #  Delta is the time that has passed since the actual ttc value is set
        ttc = self.ttc - delta
        if ttc < 0:
            ttc = 0

        self.ttc = ttc

#  The function that calculates if the load is served by G1, G2, G3 and TL1, TL2
def load_served(generators, transmission_lines, load):

    #  Calculation of the Generation by the Generators 1 and 2
    #  Generators 1 and 2 serving to the load is effected by the state of transmission line
    #  Generator 3 is directly serving the load without transmission line in middle
    Generation_by_G1_G2 = 0
    for g in generators[0:2]:
        if g.state:
            Generation_by_G1_G2 = Generation_by_G1_G2 + 75


    #  Calculation of the Transmission capacity in transmission system
    transmission_capacity = 0
    for t in transmission_lines:
        if t.state:
            transmission_capacity = transmission_capacity + 100


    #  One of the components of power that contributes to load is the
    #  minimum of the Generation by G1 and G2 and the total capacity of
    #  the transmission system.
    #  The other component of generation comes directly from G3 generator
    capacity_by_G1_G2_TL1_TL2 = min(Generation_by_G1_G2, transmission_capacity)

    #  Calculation of amount of load served by Generator 3
	#  This amount is directly served to Load. Doesnot depend on transmission capacity
    G3 = 0
    if generators[2].state:
        G3 = G3 + 75

    # Total capacity is capacity_by_G1_G2_TL1_TL2 and capacity_by_G3
    total_capacity = capacity_by_G1_G2_TL1_TL2 + G3
    return total_capacity >= load

#  Calculates the time till the failure of the next component from the random number generated
#  Formula used is below
#  x = -ln(z) / rho
def time_to_change(*, rho, randomnum):
    return -numpy.log(randomnum) / rho

# Class to generate the Transmission Line
# This is being done separately for the transmission line as
# transmission line is effected by weather conditions and
# have different failure rates in different weather conditions
class Line(TwoState):
    def __init__(self, weather, failure_rate_normal_weather, failure_rate_adverse_weather, repair_rate, random_state, name=''):
        # Creating an object of the class.
        self.weather = weather

        # Track the weather's previous state.
        self.prev_weather_state = weather.state

        # Convert the yearly transition rates to hourly transition rates
        self.failure_rate_normal_weather = failure_rate_normal_weather/8760
        self.failure_rate_adverse_weather = failure_rate_adverse_weather/8760

        super().__init__(failure_rate=self.failure_rate, repair_rate=repair_rate,
                         random_state=random_state, name=name)

    @TwoState.failure_rate.getter
    def failure_rate(self):
        if self.weather.state:
            return self.failure_rate_normal_weather
        else:
            return self.failure_rate_adverse_weather

    @TwoState.ttc.getter
    def ttc(self):
        if self.prev_weather_state != self.weather.state:
            self.ttc = self.time_to_change()
            # Reset the weather's previous state.
            self.prev_weather_state = self.weather.state
        return self._ttc

#  Class for the Daily Load conditions given in the problem
#  Daily load has five possible values
class Load:
    # Daily load curve start times of different load values
    Starting_Time = numpy.array([0, 4, 8, 12, 16])

    # The corresponding load values versus the starting time
    LOAD = numpy.array([60, 105, 205, 105, 60])

    def __init__(self):
        #  Time starts at 0
        self.index = 0
        self.state = self.LOAD[self.index]
        self.ttc = self.time_to_change()
        self.name = 'Load'

    def __repr__(self):
        return self.name

    def time_to_change(self):
        # Calculation of the time difference between now and the next state
        try:
            ttc = self.Starting_Time[self.index + 1] - self.Starting_Time[self.index]
        except IndexError:
            ttc = 24 - self.Starting_Time[self.index]
        return ttc

    def change_state_update_ttc(self):
        # Updating the index
        index = self.index + 1

        # When the index has reached the end of array, we reset the index to 0
        if index == self.Starting_Time.shape[0]:
            index = 0

        self.index = index
        self.state = self.LOAD[self.index]
        self.ttc = self.time_to_change()

    def update_ttc(self, delta):
        self.ttc = self.ttc - delta


#  Main Function of the Monte Carlo Simulation
#  All the transition rates are provided in per year
def main():
    # Creation of the random state
    random_state = numpy.random.RandomState(SEED)

	# Generator Parameters in per year
    generator_failure_rate = 36.5
    generator_repair_rate = 1095

    # Create list of three identical generators
    # There are three Generators
    # Random numbers are to be generated for three of them
    generators = [
        TwoState(failure_rate=generator_failure_rate, repair_rate=generator_repair_rate,
                 random_state=random_state, name='G{}'.format(i + 1))
                 for i in range(3)]

    #  Weather Parameters Given
    weather_failure_rate = 43.8
    weather_repair_rate = 438
    weather = TwoState(failure_rate=weather_failure_rate, repair_rate=weather_repair_rate,
                        random_state=random_state, name='Weather')


    #  Transmission Line Parameters Given
    failure_rate_normal_weather = 10
    failure_rate_adverse_weather = 100
    repair_rate = 1095

    #  Creation of transmission lines states
    #  This uses the Transmission lines class defined
    #  There are two transmission lines
    #  Random numbers are to be generated for both of them
    transmission_lines = [
        Line(weather=weather, failure_rate_normal_weather=failure_rate_normal_weather, failure_rate_adverse_weather=failure_rate_adverse_weather,
             repair_rate=repair_rate, random_state=random_state, name='T{}'.format(i+1))
             for i in range(2)]

    # Load initialize
	# Calling the Load Class
    load = Load()

    # Monte Carlo Simulation Steps
    # Algorithm:
    # 1. Generate the Random numbers
    # 2. Use the random numbers to calculate the time to change
    # 3. Invert the state of component with least time to change.
    # 4. Generate another random number for the component with inverted state
    # 5. Calculate the time to change from random number and update the TTC to this time
    # 6. Update the time to Change for the remaining components (decrease by delta)
    # 7. If there was loss of load, then update the running time of unserved load.
    # 8. Determine whether or not load is served.
    # 9. Repeat the steps until the coefficient of variance is satisfied

    # Determine the type of components that exist
    components = [*generators, weather, *transmission_lines, load]

    # Time is set to 0 initially
    time = 0
    total_time = 0
    # Initialize a variable to count the time of loss of load
    time_of_loss_of_load = 0
    # Initial load_served_status is true at time = 0
    load_served_status = True
    # Tracking load_served_status for the frequency of loss of load
    previous_load_served_status = True
    loss_of_load_count = 0
    time_of_loss_of_load_array = numpy.array([])
    time_of_loss_of_load_average = numpy.array([])
    year = 0
    coeff_variation_array = numpy.array([])

    for iteration in range(1, 1000000, 1):
        components.sort(key=lambda x: x.ttc)
        this_component = components[0]
        delta = this_component.ttc
        if delta < 0:
            raise ValueError('delta is < 0 ?')
        if not load_served_status:
            pass

            time_of_loss_of_load = time_of_loss_of_load + delta

            if previous_load_served_status:
                loss_of_load_count = loss_of_load_count + 1

        previous_load_served_status = load_served_status

        pre_event_state = this_component.state

        this_component.change_state_update_ttc()

        if time > 8760:
            year = year + 1
            # adds up the time to the previous total time
            # total time gives the total simulation time elapsed at the end
            total_time = total_time + time
            # reset time for the current year
            # time gives the current time in the current year
            time = 0
            time_of_loss_of_load_array = numpy.append(time_of_loss_of_load_array, time_of_loss_of_load)
            time_of_loss_of_load = 0

            tus_avg = time_of_loss_of_load_array.mean()
            time_of_loss_of_load_average = numpy.append(time_of_loss_of_load_average, tus_avg)

            coeff_variation = numpy.std(time_of_loss_of_load_average) / numpy.mean(time_of_loss_of_load_array)
            coeff_variation_array = numpy.append(coeff_variation_array, coeff_variation)

            #  Convergence Criteria COV is less than or equal to 5 %
            #  Additional convergence criteria used is run simulations for 100 years
            if year > 50 and coeff_variation <= 0.05:
                print('Convergence Criteria of COV 5% is met.')
                break

        else:
            time = time + delta

        for comp in components[1:]:
            comp.update_ttc(delta)

        load_served_status = load_served(generators, transmission_lines, load.state)

    #  Final Results
    print('#####################################################')
    print('Loss of Load Probability:', time_of_loss_of_load_array.sum()/total_time)
    print('Frequency of Load Loss:', loss_of_load_count / (total_time/8760), '/year')
    print('Coefficient of covariance', coeff_variation)
    print('total_time elapsed:', total_time/8760, 'years')
    print('#####################################################')

if __name__ == '__main__':
    main()