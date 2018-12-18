# Introduction to the Cyclic train timetabling program
Cyclic train timetabling based on extended time-discretized time-space network modeling framework

### Input datasets include:

(1) Input_node.csv and input_link.csv contain the railway network structure, including stations and sections;

(2) input_node.csv specifies the given line plan information, including frequency, origin station, destinatnion station, minimum and maximum dwell times in the station, minimum and maximum running times in the section.

### How to visualize the cyclic train diagram:

The C++ program for the cyclic train timetabling will generate two output data files as the input file of the cyclic train diagram visualization tool in Python, namely, output_agent_LR.csv and output_agent_ADMM.csv.

An example of the cyclic train diagram for the illustrative example 2 is shown bellow:

![1545116041](https://user-images.githubusercontent.com/45885665/50137006-c56ba800-02d4-11e9-8c04-919a73cf90ba.jpg)
