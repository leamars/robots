robots
======

Simulates an evolutionary algorithm, shown to produce communication and lying between "robots". Written in CUDA, to take advantage of mass-parallelization. If requested, I could reformat it into regular C. It would be very slow.

This is an implementation of the simulator described in this paper:
http://stephane.magnenat.net/data/Evolutionary%20Conditions%20for%20the%20Emergence%20of%20Communication%20in%20Robots%20-%20Dario%20Floreano,%20Sara%20Mitri,%20St%C3%A9phane%20Magnenat,%20Laurent%20Keller%20-%20Current%20Biology%20-%202007.pdf

================

Summary of paper:
A robot has 3 outputs: left motor, right motor, light brightness
The two independent motors control the two wheels on a robot, allowing movement
A robot has 10 inputs: 8 visual, 1 food sensor, 1 poisin sensor
The 8 visual sensors capture data from all 360 degrees around the robot, looking for light produced by other robots and edibles.
A robot has 30 bytes of DNA: each is a weight for the connection between any given output and input.
When a robot takes in it's 10 inputs, each output sums the weighted connections, and acts upon those outputs.

DNA is initially produced randomaly.
10 robots are placed in an arena.
There are two edibles, food and poisin.
Every second that a robot is near food, it gains a point.
Every second that a robot is near poisin, it loses a point.
Only 5 robots can be around an edible at a time.
The wheels move the robots around the arena.
Multiple "arena" scenarios play out, using different random DNA.

The top 20% of scorers of a generation are selected,
their DNA is randomly mated,
and is used to produce another generaton of robots.

The paper shows that after 20 generations, robots develop communication through light to help each other find food.
By the 50th generation, robots lie about the location of the food, and trick each other into eating poisin.

================

More details soon.